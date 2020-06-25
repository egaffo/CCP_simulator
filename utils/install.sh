#!/bin/bash

export CCPSIM_HOME=$(dirname $(readlink -f $0))/../

mkdir -p $CCPSIM_HOME/bin

export PATH=$CCPSIM_HOME/bin/:$PATH
## export environment variable for local R repository
export R_LIBS=$CCPSIM_HOME/tools/R_libs

## install Scons
cd $CCPSIM_HOME/bin;
## -N check if newer version is online; -n just check whether the file exists and does not send HTTP requests
wget  -nc http://prdownloads.sourceforge.net/scons/scons-local-3.1.2.tar.gz;
mkdir -p scons_dir
tar -xf scons-local-3.1.2.tar.gz -C scons_dir ;
ln -sf $CCPSIM_HOME/bin/scons_dir/scons.py $CCPSIM_HOME/bin/scons;
cd -

cd $CCPSIM_HOME/bin
ln -sf ../utils/*.py .
ln -sf ../utils/*.R .
cd -

## install other tools
$CCPSIM_HOME/bin/scons -i -f $CCPSIM_HOME/scons/install.py $1
