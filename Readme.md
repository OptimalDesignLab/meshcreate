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
### `a1tri2`
  This executable takes 3 arguments: the number of rectangles in the x
  direction, the number of rectangles in the y direction, and a flag 
  (either 0 or 1) determining whether or not the the vertices of the mesh
  are perturbed ( 0 = no perturbation).  Each rectangle is divided into
  2 triangles to produce a mesh of triangular elements.

### `a1tri2_polar`
  This executable takes the first two argument that `a1tri2` requires.
  Vertex perturbation is not supported


## Geometry classification
  The mesh generators correctly classify the meshes on geometric entities.
  Specifically, the corners, edges, and interior of the domain are all 
  geometric entities, and the mesh entities that lie on those geometric
  entities are classified correctly.

  Geometric entities are number from 0 to n, starting from the bottom left
  corner, going counter-clockwise
