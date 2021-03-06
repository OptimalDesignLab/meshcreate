#!/bin/bash
# variable to determine the number of elements
first=4
inc=4
last=24
name_prefix="squarevortex_bc_s"
# control perturbation of mesh points
perturb=0
echo "perturb = $perturb"

# current mesh index
idx=1

for i in `seq $first $inc $last`;
do
	echo "numel = $i"
#        ./a1tri2_polar $i $i
	./a1tri2 $i $i $perturb
	cp -v ./meshfiles/abc0.smb ./meshfiles/"$name_prefix"_"$idx"_0.smb

        idx=$(expr $idx + 1)
done
