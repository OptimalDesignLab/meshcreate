#!/bin/bash

# make a set of serial 3D meshes
# the output file name is of the form:
# "$name_prefix"_$idx_0.smb
# where "$name_prefix", $idx are defined below

# variable to determine the number of elements
first=5  # number of elements in each direction in the first mesh
inc=1  # increment for the number of elements
last=25  # number of elements in the last mesh
nproc=4
name_prefix="pcube3"  # the file name prefix
perturb=0  # control perturbation of mesh points
echo "perturb = $perturb"

# current mesh index
idx=1

for i in `seq $first $inc $last`;
do
	echo "numel = $i"
#        ./a1tri2_polar $i $i
        ./make_parallel_mesh3.sh $i $i $i $perturb "$name_prefix"_"$idx"_ $nproc
#        cp -v ./meshfiles/abc0.smb ./meshfiles/"$name_prefix"_"$idx"_0.smb

        idx=$(expr $idx + 1)
done
