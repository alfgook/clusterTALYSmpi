### clusterTALYS - R package

The package `cluserTALYS` enables launching calculations
in parallel over an SSH connection on computers with multiple processors,
clusters of work stations, or scientific computing clusters.

## Requirements

This package makes use of the functionality of the R package
[clusterSSH](https://github.com/gschnabel/clusterSSH) to 
communicate with the cluster, which must therefore be installed
together with its dependencies. Some of the required packages
currently use the commands `rsync`, `ssh`, `ps`, `mkpipe` and 
pipes. This means that `clusterSSH` due to its dependencies
can only be used under UNIX-like systems.
Further, the package [TALYSeval](https://github.com/gschnabel/TALYSeval)
must be installed on the computing nodes of the cluster and a  
shared network file system availabe used for communication.

## Installation

In order to install the package with all its dependencies,
execute the following commands in a terminal:
```
mkdir inst
cd inst

git clone https://github.com/gschnabel/interactiveSSH.git
git clone https://github.com/gschnabel/rsyncFacility.git
git clone https://github.com/gschnabel/remoteFunctionSSH.git
git clone https://github.com/gschnabel/clusterSSH.git
git clone https://github.com/gschnabel/clusterTALYS.git
git clone https://github.com/gschnabel/TALYSeval.git

R CMD INSTALL interactiveSSH
R CMD INSTALL rsyncFacility.git
R CMD INSTALL remoteFunctionSSH
R CMD INSTALL clusterSSH
R CMD INSTALL clusterTALYS
R CMD INSTALL TALYSeval
```

Further, R must be installed on the cluster and the 
package TALYSeval installed there.
Depending on the cluster configuration, the following on
the login node may do the job:
```
mkdir inst
cd inst

git clone https://github.com/gschnabel/TALYSeval.git
R CMD INSTALL TALYSeval
```

## Basic usage

First, we need to initialize the functionality of 
package `remoteFunctionSSH` to lauch R functions on 
a remote computer (e.g., the login node):
```
remObj <- initSSH(ssh_login, ssh_pw,
                  tempdir.loc = calcdir_loc,
                  tempdir.rem = calcdir_rem)
```
The variables `ssh_login` and `ssh_pw` have to be initialized
by the SSH credentials of the user. The variable `calcdir_loc`
must contain an existing path on the local machine, preferably 
empty. The variable `calcdir_rem` contains an existing path on the
remote machine, whose content mirrors that of `calcdir_loc`. 

Next, we need to set up the functionality of the package
`clusterSSH` to launch code given as an R function in parallel:
```
clustObj <- initCluster(functions_multinode, remObj)
```
The first variable `functions_multinode` is part of the clusterSSH
package and contains functions to manage parallel computations on
a cluster of worker nodes with a shared network filesystem but
without a task scheduler, such as the Sun Grid Engine (SGE).

In addition, a script must be manually started on the worker nodes
that should participate in computations. 
```
startNodeController_<- remObj$createRemoteFunction(
    clustObj$startNodeController,
    fun.name = "startNodeController"
)
startNodeController()
```
The last command will print out a command string, say `cmdstr`, 
that you should launch on the worker nodes. It may be convenient 
idea to execute `nohup cmdstr &` instead of just `cmdstr` to
install the controller scripts on the worker nodes permanently. 

Finally, the functionality of the `clusterTALYS` package needs to
be initialized. This can be done via the command
```
talysClustObj <- initClusterTALYS(clustObj, "talys")
```
The second argument is for the situation where the talys executable
is in of the directories specified in the environment variable `PATH`.
Otherwise, the full path to the executable must be specified.

A calculation can be now launched in the following way
(here with 56Fe and the neutron-induced total cross section
as an example):
```
library(data.table)

paramList1 <- list(projectile = "n",
                   element = "Fe",
                   mass = 56,
                   energy = c(1,2,3))

paramList2 <- list(projectile = "p",
                   element = "Fe",
                   mass = 56,
                   energy = c(1,2,3))

outSpec <- data.table(REAC = "CS/TOT",
                      L1 = seq(1, 3, length = 10),
                      L2 = 0, L3 = 0) 

severalParamLists <- list(paramList1, paramList2)

runObj <- talysClustObj$run(severalParamLists, outSpec)
```
In this example two calculations are launched in parallel.
The latter function `talysClustObj$run` returns immediately
without waiting for the termination of the calculations.
A function that does wait is `talysClustObj$eval`.

The termination of the calculations can be checked with
`talysClustObj$isRunning` and in case this function
returns `FALSE`, results retrieved by
```
result <- talysClustObj$result(runObj)
```
