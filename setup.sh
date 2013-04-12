#!/bin/sh
acqu_dir=$(dirname $(readlink -f $0))
export acqu_sys=$acqu_dir/acqu_core
export acqu=$acqu_dir/acqu_user
export CALIB=$acqu_dir/CaLib
export OSCAR=$acqu_dir/OSCAR
# use a widespread compiler
# only make in acqu_core and acqu_user need this
export CCCOMP=g++
# acqu_dir is not needed
acqu_dir=