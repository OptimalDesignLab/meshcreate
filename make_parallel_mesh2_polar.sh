#!/bin/bash

# $1 number of rectangles in x direction
# $2 number of rectangles in y direction
# $3 perturb or not
# $4 name of output mesh
# $5 number of parts to partition into

./a1tri2_polar $1 $2 $3
./make_dmg.sh
./split.sh $5 $4

