#!/bin/bash

# make a set of parallel meshes for a strong scaling study
# the output file name is of the form:
# "$name_prefix"_$idx_$nproc
# where "$name_prefix", $idx, and $nproc are defined below

# variable to determine the number of elements
first=1  # number of processes to start with
inc=1  # increment for the number of processes
last=32  # final number of processes
name_prefix="strong_small"  # the file name prefix
perturb=0  # control perturbation of mesh points
echo "perturb = $perturb"

numEl=200

# make the mesh
./a1tri2 $numEl $numEl $perturb
./make_dmg.sh

# now split it repeatedly
# current mesh index
idx=1
i=$first
#for i in `seq $first $inc $last`;
while [ $i -le $last ];
do
	echo "nprocs = $i"
#        ./a1tri2_polar $i $i
#        ./make_parallel_mesh2.sh $i $i $perturb "$name_prefix"_"$idx"_ $nproc
         ./split.sh $i "$name_prefix"_"$idx"_ $nproc

#        cp -v ./meshfiles/"$name_prefix"_"$idx"* ~/meshfiles/

        i=$(($i + $inc))
        idx=$(expr $idx + 1)
done



