#!/bin/bash

# $1 = number of parts to partition the mesh into
# $2 = name of output file (no extension)
start_dir=`pwd`

if [ "$#" -le 0 ]; then
  echo "too few arguments provided"
  exit 1
fi

if [ "$#" -le 1 ]; then
  outname="abc2"
else
  outname=$2
fi

cd ./meshfiles

mpirun -np $1 zsplit abc.dmg abc.smb $outname.smb $1
cp -v abc.dmg $outname.dmg

mpirun -np $1 verify $outname.dmg $outname.smb

cd $start_dir
