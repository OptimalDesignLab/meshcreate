#!/bin/bash

# $1 = .su2 file
# $2 = number of elements in the x direction
# $3 = number of elements in the y direction

julia ./extract_su2.jl $1  # create coords.dat
./a1tri2_su2 $2 $3
