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

#include <apfNumbering.h>

#include <queue>
#include <vector>

void checkConnectivity(apf::Mesh* m, int check_dim);
void checkCompactness(apf::Mesh* m, int d1, int d2);

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
  int dim = m->getDimension();

  // run the standard verifier
  m->verify();

  // check for dangling elements
  for (int i=0; i < dim; ++i)
    checkConnectivity(m, i);

  checkCompactness(m, dim, dim - 1);

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


// check if every entity of dimension d1 is reachable via entities of
// dimension d2
void checkCompactness(apf::Mesh* m, int d1, int d2)
{

  int num_d1 = m->count(d1);
  apf::Adjacent adj;
  std::queue<apf::MeshEntity*> que;
  std::vector<bool> seen_d1(num_d1, false);  // entities of dimension d1 already seen  

  apf::Numbering* n = apf::numberOwnedDimension(m, "compactnessnums", d1);

  apf::MeshIterator* it = m->begin(d1);
  apf::MeshEntity* e = m->iterate(it);
  int currel_num = apf::getNumber(n, e, 0, 0);
  seen_d1[currel_num] = true;

  // put first element in que
  que.push(e);

  // dont' need the iterater anymore
  m->end(it);

  while ( que.size() > 0 )
  {
    e = que.front();
    que.pop();
    currel_num = apf::getNumber(n, e, 0, 0);
    seen_d1[currel_num] = true;

    // add adjacent elements to que if they have not been seen yet
    apf::getBridgeAdjacent(m, e, d2, d1, adj);

    for (unsigned i=0; i < adj.size(); i++)
    {
      e = adj[i];
      currel_num = apf::getNumber(n, e, 0, 0);
      if ( !seen_d1[currel_num])
        que.push(e);
    }

  }  // end while


  // check results
  for (int i=0; i < num_d1; i++)
    if ( !seen_d1[i] )
    {
      std::cout << "entity " << i << " of dimension " << d1 << " is not reachable" << std::endl;
      std::abort();
    }

  // free memory
  apf::destroyNumbering(n);

}  // end checkCompactness
