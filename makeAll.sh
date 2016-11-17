#!/bin/sh

#Check correct numnber of command line arguments
if (( $# > 1 )); then
   echo "usage: makeAll.sh [c compiler choice: defaults to gcc if not specified]"
  exit 1
fi

if (( $# < 1 )); then
    compiler=gfortran
    echo "Defaulting to gfortran compiler"
fi



#Check correct numnber of command line arguments
if [[ $# -eq 1 ]]; then
    compiler=$1
   echo "Using compiler $compiler."
fi

make GF=$compiler

cp GLaMM_model ~/bin

