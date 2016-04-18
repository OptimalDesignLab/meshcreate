#!/bin/bash
# variable to determine the number of elements
first=12
inc=1
last=32

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
	cp -v ./meshfiles/0.smb ~/pdesolver/mesh_files/square3_"$idx"0.smb

        idx=$(expr $idx + 1)
done
