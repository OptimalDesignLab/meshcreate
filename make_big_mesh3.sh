#!/bin/bash

# make a set of parallel 3D meshes
# this directly calls all executables, because the BGQ can't run a shell
# script as a job

# the output file name is of the form:
# "$name_prefix"_$idx_0.smb
# where "$name_prefix", $idx are defined below

# variable to determine the number of elements

first=13  # number of elements in each direction in the first mesh
inc=1  # increment for the number of elements
last=13  # number of elements in the last mesh
nnodes=8  # number of BGQ nodes
nproc=512
name_prefix="cube_p3_gamma_"  # the file name prefix

# current mesh index
idx=1

for i in `seq $first $inc $last`;
do
    echo "numel = $i"
    outname="$name_prefix"
    # temporarily disable idx
#    outname="$name_prefix"_"$idx"_
#    echo "outname = $outname"

    # make serial mesh meshfiles/abc.smb
    srun -N 1 ./a1tet $i $i $i  # make abc.smb

    start_dir=`pwd`
    cd ./meshfiles

    # make dmg
    srun -N 1 mkmodel abc.smb abc.dmg

    # split it
    srun -N $nnodes -n $nproc --overcommit zsplit abc.dmg abc.smb $outname.smb $nproc
    cp -v ./abc.dmg $outname.dmg # copy model file to new name


    # verify it
    srun -N $nnodes -n $nproc --overcommit verify $outname.dmg $outname.smb

    # increment idx
    idx=$(expr $idx + 1)

    # move back into starting directory
    cd $start_dir
done
