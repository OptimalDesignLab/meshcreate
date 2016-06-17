#!/bin/bash

# $1 = split fname
# $2 = np
# $3 = refined fname

start_dir=`pwd`
cd ./meshfiles
mkmodel abc.smb abc.dmg
# split
mpirun -np $2 zsplit abc.dmg abc.smb $1.smb $2
echo "finished splitting original mesh"
cp -v ./abc.dmg ./$1.dmg
echo "finished making model of split mesh"
mpirun -np $2 render $1.dmg $1.smb $1
echo "finished rendering split mesh"

# refine
mpirun -np $2 uniform $1.dmg $1.smb $3.smb
echo "finished uniform refinement"
cp -v ./$1.dmg ./$3.dmg
echo "finished making model"
mpirun -np $2 render $3.dmg $3.smb $3
echo "finished rendering refined mesh"

cd $start_dir
