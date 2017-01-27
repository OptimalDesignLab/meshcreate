#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfShape.h>
#include <cstdlib>
#include <cassert>
#include <iostream>


void checkConnectivity(apf::Mesh* m, int check_dim);

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << ": <model> <mesh>" << std::endl;
    MPI_Finalize();
    return 1;
  }

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
    std::cout << "running strict verifier" << std::endl;

  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1], argv[2]);

  // run the standard verifier
  m->verify();

  // check for dangling elements
  for (int i=0; i < m->getDimension(); ++i)
    checkConnectivity(m, i);


  // cleanup
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
} // function main


// this function checks that every element is connected to
// another element via an entity of dimension check_dim
void checkConnectivity(apf::Mesh* m, int check_dim)
{

  int mesh_dim = m->getDimension();
  apf::MeshEntity* e;
  apf::Downward down;  // downward adjacent entities of dim check_dim
  apf::DynamicArray<apf::MeshEntity*> adjacencies;  // upward adjacent element
  apf::MeshIterator* it = m->begin(mesh_dim);

  int elnum = 0;
  while ( (e = m->iterate(it)) )
  {
    int ndown = m->getDownward(e, check_dim, down);

    int num_up = 0;
    for (int i=0; i < ndown; ++i)
    {
      m->getAdjacent(down[i], mesh_dim, adjacencies);
      num_up += adjacencies.getSize();
    }

    if (num_up <= ndown)
    {
      std::cerr << "element " << elnum << "is dangling by a dimension " << check_dim << "entity" << std::endl;
      std::abort();
    }

    elnum += 1;

  }  // end while

} // function checkConnectivity
