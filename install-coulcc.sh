#!/bin/bash

set -e # exit on error

usage()
{
  echo "nQMCC Install Script"
  echo ""
  echo "USAGE:"
  echo "./install-nqmcc.sh [--help|-h] | [-p|--prefix=PREFIX] | [-m|--machine=MACHINE] | [-g|--debug=DEBUG] | [-v|--VERBOSE=VERBOSE] "
  echo ""
  echo " --help             Display this help text"
  echo " --prefix=PREFIX    Install binary in 'PREFIX/bin'"
  echo "                    Default prefix='Path/To/nQMCC/build/bin'"
  echo " --debug=DEBUG      Compile with -g"
  echo "                    Default off"
  echo " --verbose=VERBOSE  Print Compiler information"
  echo "                    Default off"
  echo " --machine=MACHINE  sets compiler and flags"
  echo "                    Default MACHINE='gnu'"
  echo "                    Available:          '"
  echo "                               gnu  '"
  echo "                               intel'"
  echo "                               cray'"
  echo ""
}

PREFIX="$(pwd)/build"
MACHINE="gnu"
DEBUG=""
VERBOSE=""
PHI_R4=0
while [ "$1" != "" ]; do
  PARAM=$(echo "$1" | awk -F= '{print $1}')
  VALUE=$(echo "$1" | awk -F= '{print $2}')
  case $PARAM in
    -h | --help)
      usage
      exit
      ;;
    -p | --prefix)
      PREFIX=$VALUE
      ;;
    -m | --machine)
      MACHINE=$VALUE
      ;;
    -g | --debug)
      DEBUG=" -g "
      ;;
    -v | --verbose)
      VERBOSE=" --verbose"
      ;;
    *)
      echo "ERROR: unknown parameter \"$PARAM\""
      usage
      exit 1
      ;;
  esac
  shift
done

set -u # error on use of undefined variable

case $MACHINE in
  gnu)
  CMK_FC=${CMK_FC:-"gfortran"}
  #CMK_FC=${CMK_FC:-"gfortran"}
  CMK_FLAG="-Wall -cpp -O3 -std=f2023 -ffree-line-length-512 -fbacktrace -fcheck=all -pedantic"
  ;;
  improv)
  CMK_FC=${CMK_FC:-"ifx"}
  CMK_FLAG="-warn -cpp -O3"
  ;;
  cray)
  CMK_FC=${CMK_FC:-"ftn"}
  CMK_FLAG="-O3 -e Z"
  ;;

  *)
  echo "ERROR: unknown machine \"$MACHINE\""
      usage
      exit 1
      ;;
esac

CI=${CI:-"false"} # GitHub Actions workflows set CI=true
CMK_FLAG=${CMK_FLAG}" ${DEBUG} ${VERBOSE}"

#buildit
cmake -S . -B ${PREFIX} -DCMAKE_Fortran_COMPILER=${CMK_FC}  \
                        -DCMAKE_Fortran_FLAGS="${CMK_FLAG}"\
                        -DCOULCC_MASTER_PROJECT="ON"\
