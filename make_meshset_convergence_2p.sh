#!/bin/bash

# make a set of parallel meshes
# the output file name is of the form:
# "$name_prefix"_$idx_$nproc
# where "$name_prefix", $idx, and $nproc are defined below

# variable to determine the number of elements
first=10  # number of elements in each direction in the first mesh
inc=5  # increment for the number of elements
last=15  # number of elements in the last mesh
nproc=16  # number of processors to partition the mesh into
name_prefix="vortexcurve10_p1_p"  # the file name prefix
perturb=0  # control perturbation of mesh points
echo "perturb = $perturb"

# current mesh index
idx=1

for i in `seq $first $inc $last`;
do
	echo "numel = $i"
#        ./a1tri2_polar $i $i
        ./make_parallel_mesh2.sh $i $i $perturb "$name_prefix"_"$idx"_ $nproc
#        cp -v ./meshfiles/"$name_prefix"_"$idx"* ~/meshfiles/

        idx=$(expr $idx + 1)
done
