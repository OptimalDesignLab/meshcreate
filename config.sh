#!/bin/bash

mkdir -pv ./install
installdir="`pwd`/install"
echo $installdir

cd ./src
mkdir -vp ./build
cd ./build

cmake .. \
  -DSCOREC_PREFIX=$SCOREC_PREFIX \
  -DCMAKE_CXX_COMPILER="mpicxx" \
  -DCMAKE_C_FLAGS="-O2 -g -Wall" \
  -DCMAKE_CXX_FLAGS="-O2 -g -Wall -std=c++11" \
  -DENABLE_THREADS=OFF \
  -DCMAKE_INSTALL_PREFIX=$installdir
#  -DIS_TESTING=True
#  -DENABLE_ZOLTAN=


#  -DCMAKE_C_COMPILER="mpicc" \
