#!/bin/bash

export CCPSIM_HOME=$(dirname $(readlink -f $0))

export PATH=$CCPSIM_HOME/bin/:$PATH
## export environment variable for local R repository
export R_LIBS=$CCPSIM_HOME/bin/R_libs
export CIRI_SIM_EXEC=$CCPSIM_HOME/bin/CIRI_simulator.pl

$CCPSIM_HOME/bin/scons -f $CCPSIM_HOME/scons/main.py $1 CIRI_SIM_EXEC=$CIRI_SIM_EXEC
