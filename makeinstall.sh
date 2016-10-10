#!/bin/bash

startdir=`pwd`
cd ./src/build
make
make install

cd $startdir
cp -vsn ./install/bin/* ./
