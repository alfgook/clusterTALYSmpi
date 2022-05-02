#  The MIT License
#  
#  Copyright (c) 2019 Alf Göök 
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
#'
#' @return
#' list with functions to perform TALYS calculations
#' @export
#'
#' @import data.table
#' @useDynLib clusterTALYSmpi initalize_mpi
#' @useDynLib clusterTALYSmpi finalize_mpi
#' @useDynLib clusterTALYSmpi start_mpi_workers

initClusterTALYSmpi <- function(talysExe="talys", runOpts=NULL) {

  # will need to resolve this, should not have to specify a complete path to the .so file
  #dyn.load(paste0("/home/alf/programs/talys-mpi/runTALYSmpi/start_mpi_workers", .Platform$dynlib.ext))
  .C("initalize_mpi")

  close <- function(env) {
    .C("finalize_mpi")
  }

  ee <- environment()
  reg.finalizer(ee, close, onexit = TRUE)

  defaults <- list(runOpts=runOpts)
  theResults <- NA

  stopifnot(("maxNumCPU" %in% names(runOpts)),is.numeric(runOpts$maxNumCPU))

  runTALYS <- function(inpSpecList, outSpec, runOpts=NULL, calcsPerJob=NULL, pollTime=0, saveDir=NULL) {

    if (!is.list(inpSpecList) || !all(sapply(inpSpecList, is.list)))
      stop("inpSpecList must be a list of TALYS inputs")
    if (!is.data.table(outSpec) && !(is.list(outSpec) &&
                                    all(sapply(outSpec, is.data.table)) &&
                                    length(inpSpecList) == length(outSpec)))
      stop(paste0("outSpec must be either a datatable or ",
                  "a list of datatables of the same length as inpSpecList"))

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

    if (is.null(runOpts))
      runOpts <- defaults$runOpts

    if (is.data.table(outSpec))
      input <- lapply(seq_along(inpSpecList),function(i)
        list(input=inpSpecList[[i]], outspec=outSpec,saveDir=saveDir, calcIdx=i, calcDir=""))
    else if (is.list(outSpec))
    {
      input <- mapply(function(i, x, y) {
        list(input=x, outspec=y, saveDir = saveDir, calcIdx=i, calcDir="")
      }, i=seq_along(inpSpecList), x=inpSpecList, y=outSpec, SIMPLIFY = FALSE)
    }
    else
      stop("should not happen")

    jobList <- replicate(length(input),NULL,simplify=FALSE)
    talysMod <- createModelTALYS()
    for (jobIdx in seq_along(input)) {

      # create the temporary directory
      globalTempdir <- runOpts$TMPDIR
      if (globalTempdir=="") globalTempdir <- Sys.getenv("TMPDIR")
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
      
      curInp <- input[[jobIdx]]$input
      dummyFun <- function() {
        talysMod$prepare(basedir,curInp)
      }
      tryCatch(withCallingHandlers(dummyFun(), error = conditionFun), error=function(e) e)

      jobList[[jobIdx]] <- basedir
      input[[jobIdx]]$calcDir <- basedir
    }

    #create the command string
    #nbrCPU <- 32 # I need to figure out where to take this from
    #.C("start_mpi_workers",as.integer(nbrCPU),as.character("runTALYSmpi"),as.character(do.call(paste,jobList)));
    #directories <- do.call(paste,jobList)

    #run the jobs
    base_wd <- getwd()

    if(("bindir" %in% names(runOpts))) {
      bin_path <- runOpts$bindir
    } else {
      bin_path <- "/usr/local/bin"
    }

    if(("maxNumCPU" %in% names(runOpts))) {
      maxNumCPU <- runOpts$bindir
    } else {
      maxNumCPU <- length(jobList)
    }

    .C("start_mpi_workers",
        worker_program = as.character("runTALYSmpi"),
        job_list = as.character(jobList),
        number_of_jobs = as.integer(length(jobList)),
        number_of_workers = as.integer(maxNumCPU),
        talys_exe = as.character(talysExe),
        bin_path = as.character(bin_path));

    #get the results
    resultList <- replicate(length(input),NULL,simplify=FALSE)
    for (jobIdx in seq_along(input)) {
      #resultList[[jobIdx]] <- readResults(inputList[[jobIdx]])

      curSaveDir <- input[[jobIdx]]$saveDir
      curCalcIdx <- input[[jobIdx]]$calcIdx

      curInp <- input[[jobIdx]]$input
      resultList[[jobIdx]]$input <- curInp

      # save talys output
      setwd(input[[jobIdx]]$calcDir)
      numLines <- as.integer(system("wc -l output | awk '{print $1}'", intern=TRUE))
      if (numLines <= 60)
        outputSummary <- system("cat output", intern=TRUE)
      else {
        cmdstr <- "head -n 15 output; echo --- LINES SKIPPED ---; tail -n 15 output"
        outputSummary <- system(cmdstr, intern=TRUE)
      }
      resultList[[jobIdx]]$output <- outputSummary

      # handle errors gracefully
      dummyFun <- function() {
        #talysMod$finalize(input[[jobIdx]]$calcDir)
        talysMod$read(input[[jobIdx]]$calcDir, copy(input[[jobIdx]]$outspec), packed=FALSE)
      }
      curRes <- tryCatch(withCallingHandlers(dummyFun(),
                                             error = conditionFun),
                         error=function(e) e)

      if (exists(stackvarname,.GlobalEnv)) {
        resultList[[jobIdx]]$error <- get(stackvarname, .GlobalEnv)$condition
        resultList[[jobIdx]]$stacktrace <- get(stackvarname, .GlobalEnv)$stack
      } else {
        resultList[[jobIdx]]$result <- curRes
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

      unlink(list.files(input[[jobIdx]]$calcDir))
      unlink(input[[jobIdx]]$calcDir, recursive=TRUE)
    }
    # delete the temporary calculation directory
    unlink(basedir, recursive=TRUE)

    setwd(base_wd)
    theResults <<- resultList
    jobList

    # it would make more sense, and be simpler to just return the result here,
    # but to be compatible with Georgs clusterTalys.R the result will be returned
    # from the getResults function. This is so I can use the same controlling 
    # scripts for both the present clusterTALYSmpi and the original clusterTALYS codes.
    # All that is needed is to change the TalysHnd in the config file
  }

  getResults <- function(jobList, selection=TRUE) {
    return(theResults)
  }

  isRunningTALYS <- function(jobList,combine=TRUE) {
    # dummy function to be compatible with georgs clusterTALYS, will always return FALSE
    
    stopifnot(length(jobList)>=1)
    runStatus <- rep(FALSE,length(jobList))
    #for (idx in seq(length(jobList),1)) {
    #  runStatus[idx] <- clust$isRunning(jobList[[idx]])
    #  if (runStatus[idx] && isTRUE(combine))
    #    return(TRUE)
    #}
    if (isTRUE(combine))
      FALSE
    else
      runStatus
  }

  splitJob <- function(inpSpecList, outSpec, runOpts=NULL, calcsPerJob=NULL, pollTime=0, saveDir=NULL) {
    if (!is.list(inpSpecList) || !all(sapply(inpSpecList, is.list)))
      stop("inpSpecList must be a list of TALYS inputs")
    if (!is.data.table(outSpec) && !(is.list(outSpec) &&
                                    all(sapply(outSpec, is.data.table)) &&
                                    length(inpSpecList) == length(outSpec)))
      stop(paste0("outSpec must be either a datatable or ",
                  "a list of datatables of the same length as inpSpecList"))

      if (is.null(runOpts))
      runOpts <- defaults$runOpts

      if (is.data.table(outSpec) || (is.list(outSpec) && length(outSpec)<runOpts$maxNumCPU) ) { # single job
        jobList <- runTALYS(inpSpecList,outSpec,runOpts=runOpts)
      } else if (is.list(outSpec)) {
        # split requested calculations into several jobs
        InpChunks <- split(inpSpecList, ceiling(seq_along(lst)/runOpts$maxNumCPU))
        outChunks <- split(outSpec, ceiling(seq_along(lst)/runOpts$maxNumCPU))
      } else {
        stop("should not happen")
      }

      inputList <- split(input,groupfact)
      # start the jobs
      jobList <- replicate(length(inputList),NULL,simplify=FALSE)
      for (jobIdx in seq_along(inputList)) {
        jobList[[jobIdx]] <- runTALYS(inputList[[jobIdx]],runOpts=runOpts)
        runTALYS <- function(inpSpecList, outSpec, runOpts=NULL, calcsPerJob=NULL, pollTime=0, saveDir=NULL)
      }

      jobList
  }

  list(run=splitJob,result=getResults,isRunning=isRunningTALYS,close=close)
}
