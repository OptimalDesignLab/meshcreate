#!/bin/bash

startdir=`pwd`
cd ./src/build
make -j 4
make install

cd $startdir
cp -vsn ./install/bin/* ./
