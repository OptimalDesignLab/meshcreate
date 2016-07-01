#!/bin/bash

# $1 number of cubes in x direction
# $2 number of cubes in y direction
# $3 number of cubes in the z direction
# $4 perturb or not
# $5 name of output mesh
# $6 number of parts to partition into

./a1tet $1 $2 $3 $4
./make_dmg.sh
./split.sh $6 $5
