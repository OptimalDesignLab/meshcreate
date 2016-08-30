# createSMBMesh
This repo contains some small programs useful for creating structured meshes
for Pumi.

  `a1tri2`: creates a mesh of triangular elements on a square domain
  `a1tri2_polar`: creates a mesh of triangular elements on a quarter 
                  circle domain
  `a1tet`: creates a regular mesh of tetrahedra on a square domain

## Compiling
  The script `build.scorec.sh2` is used to build the programs.  To compile
  `a1tri2.cc` for example, it is invoked as follows:
```
  ./build.scorec.sh2 a1tri2
```

## Running
### `a1tri2` (serial)
  This executable takes 3 arguments: the number of rectangles in the x
  direction, the number of rectangles in the y direction, and a flag 
  (either 0 or 1) determining whether or not the the vertices of the mesh
  are perturbed ( 0 = no perturbation).  Each rectangle is divided into
  2 triangles to produce a mesh of triangular elements.

### `a1tri2_polar` (serial)
  This executable takes the first two argument that `a1tri2` requires.
  Vertex perturbation is not supported

### `a1tet` (serial)
  This execuable takes 3 required arguments and one optional argument.  The
  first 3 specify the number of cubes in the x, y, and z directions, 
  respectively.  Each cube is divided into 6 tetrahedra.  The fourth option 
  controls vertex perturbation.  If not specified, no perturbation is applied

### Scripts
  Several shell scripts are provided to automate common tasks.  All scripts
  have a name, followed by a number and possibly a letter.  The number 
  describes the dimensional of the creates meshes, and the letter describes
  whether the resulting mesh is serial or parallel.  For example,
  `make_meshset_convergence_3s.sh` performs the action of 
  `make_meshset_convergence.sh` (described below), resulting in a 3 dimensional
  serial mesh.

### Parallel Meshes
  The easiest way to create a parallel mesh is to use the script 
  `make_parallel_mesh.sh`.  It takes 5 arguments, the first three are the
   same as `a1tri21`, the fourth is the name of the output mesh file (without
  extensions), and the fifth is the number of processors to partition the mesh
  into.  Pumi will create a both a `$4.dmg` file and a set of `"$4""$rank".smb`
  files, where `rank` is the mpi rank that will load the mesh file.
  Note that these scripts require the `mkmodel`, `zsplit`, and `verify` 
  programs to be found via the `$PATH` variable.  They are included with both
  the Scorec installation of Pumi and the local build of PumiInterface.

  Alternatively, the script `refine_and_split.sh` takes an exising serial mesh,
  partitions it into the specified number of processes, and then performs 
  uniform mesh refinement to the partitioned mesh.  This doubles the number 
  of elements in each direction (ie. in 2d the resulting mesh has 4x the number
  of elements).  See the top of the script for arguments.

### Making Sets of Meshes
  Two additional scripts exist for making sets of meshes.  
  `make_meshset_convergence.sh` makes a set of meshes with increasing 
  numbers of elements.  It takes no arguments, see the comments at the top
  of the script for details on how to specify the sizing and naming of the 
  meshes.

  The script `make_meshset_strong.sh` produces a set of meshes with the 
  same number of elements divided among an increasing number of processes.
  See the comments at the top of the script to specify sizing and partition.

## Geometry classification
  The mesh generators correctly classify the meshes on geometric entities.
  Specifically, the corners, edges, and interior of the domain are all 
  geometric entities, and the mesh entities that lie on those geometric
  entities are classified correctly.

  Geometric entities are number from 0 to n, starting from the bottom left
  corner, going counter-clockwise

## Periodic Meshes
  Periodic mesh can be created by setting the boolean variable *periodic near
  the top of each source file.  Do not directly set the fields of the
  Periodic struct!
