# createSMBMesh
This repo contains some small programs useful for creating structured meshes
for Pumi.

  `a1tri2`: creates a mesh of triangular elements on a square domain

  `a1tri2_polar`: creates a mesh of triangular elements on a quarter 
                  circle domain

  `a1tet`: creates a regular mesh of tetrahedra on a square domain


## Compiling
  CMake is used internally to build all executables and install them to 
  the `repo/install` directory.  They are also symlinked into the current
  directory for ease of use.
  The `config.sh` script configures the installation, and the `makeinstall.sh`
  script builds and installs all executables.
  The environmental variable `SCOREC_PREFIX` can be used to specify 
  a Pumi installation to link against.  If unspecified, CMakes
  `find_package` function will be used to locate one.  Note that the 
  `use_julialib.sh` script in PumiInterface will set `SCOREC_PREFIX` if needed
  to ensure this package is linking against the same Pumi installation as 
  PumiInterface.  In general, it is not required to link against the same
  Pumi, although it seems like the most convenient thing to do.

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
### 2D:
  The mesh generators correctly classify the meshes on geometric entities.
  Specifically, the corners, edges, and interior of the domain are all 
  geometric entities, and the mesh entities that lie on those geometric
  entities are classified correctly.

  Geometric entities are number from 0 to n, starting from the bottom left
  corner, going counter-clockwise

### 3D:
  The geometric faces are classified as follows:
  xy plane at zmin: 0
  xz plane at ymin: 1
  yz plane at xmax: 2
  xz plane at ymax: 3
  yz plane at ymin: 4
  xy plane at zmax: 5

  The edges are classified as:
  xmin to xmax at ymin, zmin: 0
  ymin to ymax at xmax, zmin: 1
  xmax to xmin at ymax, zmin: 2
  ymax to ymin at xmin, zmin: 3

  Edges 4 through 7 correspond to edges 0 through 3, respectively, at zmax.

  zmin to zmax at xmin, ymin: 8
               at xmax, ymin: 9
               at xmax, ymax: 10
               at xmin, ymax: 11

  The vertices are classified as follows:
  (xmin, ymin, zmin): 0
  (xmax, ymin, zmin): 1
  (xmax, ymax, zmin): 2
  (xmin, ymax, zmin): 3

  Vertices 4 through 7 correspond to vertices 0 through 3, respectively at zmax.

## Periodic Meshes
  Periodic mesh can be created by setting the boolean variable *periodic near
  the top of each source file.  Do not directly set the fields of the
  Periodic struct!

## Other Utilities
  The `renderit` program loads a specified mesh and model and writes a vtk file.
  The `print_geo_coords` program loads a mesh and model and prints the coordinates
  of all nodes on a specified geometric entity.
