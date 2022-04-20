#  The MIT License
#  
#  Copyright (c) 2019 Georg Schnabel 
#  
#  Permission is hereby granted, free of charge, 
#  to any person obtaining a copy of this software and 
#  associated documentation files (the "Software"), to 
#  deal in the Software without restriction, including 
#  without limitation the rights to use, copy, modify, 
#  merge, publish, distribute, sublicense, and/or sell 
#  copies of the Software, and to permit persons to whom 
#  the Software is furnished to do so, 
#  subject to the following conditions:
#  
#  The above copyright notice and this permission notice 
#  shall be included in all copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
#  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
#  ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
#  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#  

#' Setup Cluster for TALYS
#'
#' Setup a computing cluster for TALYS calculations
#'
#' @param clust object as created by the package \code{clusterSSH}
#'
#' @return
#' list with functions to perform TALYS calculations
#' @export
#'
#' @import clusterSSH data.table
initClusterTALYS <- function(clust, talysExe,
                             runOpts=NULL, calcsPerJob=100) {

  defaults <- list(runOpts=runOpts, calcsPerJob=calcsPerJob)

  # create function performTALYS with string in talysExe incepted
  performTALYS <- eval(parse(text=deparse(substitute(function(input) {

    stopifnot(is.list(input))
    library(digest)
    library(TALYSeval)
    library(data.table)

    # error handling
    stackvarname <- paste0("stackinfo_",paste0(sample(letters, 10), collapse=""))
    conditionFun <- function(cond) {
      stacktrace <- sys.calls()
      stacktraceChar <- sapply(stacktrace,deparse)
      startPos <- match("dummyFun()",stacktraceChar) + 1
      endPos <- length(stacktrace)-2
      stacktrace <- stacktrace[startPos:endPos]
      #stacktrace[[endPos]][[1]] <- as.name(get("funAttrs",.GlobalEnv)$fun.name)
      assign(stackvarname, list(condition=cond, stack=stacktrace), .GlobalEnv)
      return(cond)
    }

    # create the temporary directory
    globalTempdir <- Sys.getenv("TMPDIR")
    if (globalTempdir=="") globalTempdir <- getwd()
    dir.create(globalTempdir, showWarnings=FALSE)
    # create the temporary directory for the current calculation
    cnt <- 0
    succ <- FALSE
    while (!succ && (cnt<-cnt+1) < 10) {
      proposedDirname <- sprintf("tmpcalc_%s", paste0(sample(letters,10),collapse=""))
      proposedPathname <- file.path(globalTempdir, proposedDirname)
      succ <- dir.create(proposedPathname, showWarnings = FALSE)
    }
    if (!succ) stop(paste0("unable to create temporary directory in ", globalTempdir))
    basedir <- proposedPathname
    setwd(basedir)

    talysMod <- createModelTALYS()
    result <- list()
    for (idx in seq_along(input))
    {
      curSaveDir <- input[[idx]]$saveDir
      curCalcIdx <- input[[idx]]$calcIdx

      result[[idx]] <- list()
      curInp <- input[[idx]]$input
      result[[idx]]$input <- curInp

      dummyFun <- function() {
        talysMod$prepare(basedir,curInp)
      }
      tryCatch(withCallingHandlers(dummyFun(),
                                   error = conditionFun),
               error=function(e) e)

      if (!exists(stackvarname, .GlobalEnv)) {

        # execute it
        cmdstr <- paste0(talysExe," < input > output")
        system(cmdstr)

        # save talys output
        numLines <- as.integer(system("wc -l output | awk '{print $1}'", intern=TRUE))
        if (numLines <= 60)
          outputSummary <- system("cat output", intern=TRUE)
        else
        {
          cmdstr <- "head -n 15 output; echo --- LINES SKIPPED ---; tail -n 15 output"
          outputSummary <- system(cmdstr, intern=TRUE)
        }
        result[[idx]]$output <- outputSummary

        # handle errors gracefully
        dummyFun <- function() {
          talysMod$read(".", copy(input[[idx]]$outspec), packed=FALSE)
        }
        curRes <- tryCatch(withCallingHandlers(dummyFun(),
                                               error = conditionFun),
                           error=function(e) e)
      }

      if (exists(stackvarname,.GlobalEnv)) {
        result[[idx]]$error <- get(stackvarname, .GlobalEnv)$condition
        result[[idx]]$stacktrace <- get(stackvarname, .GlobalEnv)$stack
      } else {
        result[[idx]]$result <- curRes
      }

      # clean up the calculation
      if (!is.null(curSaveDir)) {
        tarfile <- sprintf("calc.%04d.tar", curCalcIdx)
        saveFilePath <- file.path(curSaveDir, tarfile)
        tarcmd <- paste0('tar -czf ', tarfile,' *')
        movecmd <- paste0('rsync --remove-source-files -av ', tarfile, ' ', saveFilePath)
        if (system(tarcmd, intern=FALSE) != 0)
          stop(paste0("Problem with: ", tarcmd))
        if (system(movecmd, intern=FALSE) != 0)
          stop(paste0("Problem with: ", movecmd))
      }

      unlink(list.files(basedir))
    }
    # delete the temporary calculation directory
    unlink(basedir, recursive=TRUE)
    result
  }, env=list(talysExe=talysExe)))))

  isRunningTALYS <- function(jobList,combine=TRUE) {
    stopifnot(length(jobList)>=1)
    runStatus <- rep(FALSE,length(jobList))
    for (idx in seq(length(jobList),1)) {
      runStatus[idx] <- clust$isRunning(jobList[[idx]])
      if (runStatus[idx] && isTRUE(combine))
        return(TRUE)
    }
    if (isTRUE(combine))
      FALSE
    else
      runStatus
  }

  getResultTALYS <- function(jobList, selection=TRUE) {

    stopifnot(isTRUE(selection == TRUE) || is.numeric(selection))

    if (isTRUE(selection))
    {
      unlist(lapply(jobList,clust$result),recursive=FALSE)
    }
    else
    {
      numTasksPerJob <- sapply(jobList, clust$numberOfTasks)
      taskStartIdx <- cumsum(c(0, numTasksPerJob))
      jobIdx <- findInterval(selection, taskStartIdx, left.open=TRUE)
      assignDt <- data.table(jobIdx = jobIdx, globalTaskIdx = selection)
      assignDt[, taskIdx := globalTaskIdx - taskStartIdx[jobIdx]]
      assignDt[, result := list(clust$result(jobList[[.BY[["jobIdx"]]]], taskIdx)), by="jobIdx"]
      setkey(assignDt, globalTaskIdx)
      results <- assignDt[J(selection), result]
      unlist(results, recursive=FALSE)
    }
  }

  runTALYS <- function(inpSpecList, outSpec, runOpts=NULL,
                       calcsPerJob=NULL, pollTime=0, saveDir=NULL) {

    if (!is.list(inpSpecList) || !all(sapply(inpSpecList, is.list)))
      stop("inpSpecList must be a list of TALYS inputs")
    if (!is.data.table(outSpec) && !(is.list(outSpec) &&
                                    all(sapply(outSpec, is.data.table)) &&
                                    length(inpSpecList) == length(outSpec)))
      stop(paste0("outSpec must be either a datatable or ",
                  "a list of datatables of the same length as inpSpecList"))
    stopifnot(is.numeric(pollTime))

    if (is.null(runOpts))
      runOpts <- defaults$runOpts
    if (is.null(calcsPerJob))
      calcsPerJob <- defaults$calcsPerJob

    # split requested calculations into several jobs
    if (is.data.table(outSpec))
      input <- lapply(seq_along(inpSpecList),function(i)
        list(input=inpSpecList[[i]], outspec=outSpec,
             saveDir=saveDir, calcIdx=i))
    else if (is.list(outSpec))
    {
      input <- mapply(function(i, x, y) {
        list(input=x, outspec=y, saveDir = saveDir, calcIdx=i)
      }, i=seq_along(inpSpecList), x=inpSpecList, y=outSpec, SIMPLIFY = FALSE)
    }
    else
      stop("should not happen")

    groupfact <- rep(seq_len(ceiling(length(inpSpecList)/calcsPerJob)),
                     each=calcsPerJob)[seq_along(input)]
    inputList <- split(input,groupfact)
    # start the jobs
    jobList <- replicate(length(inputList),NULL,simplify=FALSE)
    for (jobIdx in seq_along(inputList)) {
      jobList[[jobIdx]] <- clust$start(performTALYS,inputList[[jobIdx]],runOpts)
    }
    if (pollTime > 0) {
      while(isRunningTALYS(jobList)) {
        Sys.sleep(pollTime)
      }
    }
    jobList
  }

  list(run=runTALYS,result=getResultTALYS,isRunning=isRunningTALYS)
}
