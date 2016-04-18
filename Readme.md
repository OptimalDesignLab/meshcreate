# createSMBMesh
This repo contains some small programs useful for creating structured meshes
for Pumi.

  `a1tri2`: creates a mesh of triangular elements on a square domain
  `a1tri2_polar`: creates a mesh of triangular elements on a quarter 
                  circle domain
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

### Parallel Meshes
  The easiest way to create a parallel mesh is to use the script 
  `make_parallel_mesh.sh`.  It takes 5 arguments, the first three are the
   same as `a1tri21`, the fourth is the name of the output mesh file (without
  extensions), and the fifth is the number of processors to partition the mesh
  into.  Pumi will create a both a `$4.dmg` file and a set of `"$4""$rank".smb`
  files, where `rank` is the mpi rank that will load the mesh file.
  Note that these scripts require the `mkmodel`, `zsplit`, and `verify` 
  programs to be found via the `$PATH` variable.  They are not (currently) 
  part of the non-Scorec installations of Pumi.

### Making Sets of Meshes
  Two additional scripts exist for making sets of mesh with increasing numbers
  of elements.  In serial, `make_meshes.sh` can be used to make a set of meshes.
  The script takes no arguments, see the comments at the top of the script 
  for details on how to specify the sizes of the meshes.
  In parallel, `make_meshset.sh` can be used.  It calles `make_parallel_mesh.sh`
  internally to create each mesh.  As with `make_meshes.sh`, it takes no
  arguments and there is a description of how to specify the relevent parameters
  at the top of the file.

## Geometry classification
  The mesh generators correctly classify the meshes on geometric entities.
  Specifically, the corners, edges, and interior of the domain are all 
  geometric entities, and the mesh entities that lie on those geometric
  entities are classified correctly.

  Geometric entities are number from 0 to n, starting from the bottom left
  corner, going counter-clockwise
