#!/bin/bash

# make a set of parallel meshes
# variable to determine the number of elements
first=10
inc=1
last=30
nproc=1
name_prefix="square4"
# control perturbation of mesh points
perturb=0
echo "perturb = $perturb"

# current mesh index
idx=1

for i in `seq $first $inc $last`;
do
	echo "numel = $i"
#        ./a1tri2_polar $i $i
        ./make_parallel_mesh.sh $i $i $perturb "$name_prefix"_"$idx"_ $nproc
        cp -v ./meshfiles/"$name_prefix"_"$idx"* ~/meshfiles/

        idx=$(expr $idx + 1)
done
