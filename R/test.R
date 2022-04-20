testfunc <- function(...) {
    nbrCPU <- 4 # I need to figure out where to take this from
	#cmdstr <- paste0("mpirun -np ", nbrCPU, " runTALYSmpi ")
	#job_str <- do.call(paste,...)
	cmdstr <- paste0("mpirun -np ", nbrCPU, " runTALYSmpi ", do.call(paste,...))
	print(cmdstr)
	print("------")
}