#include <gmi_mesh.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfMDS.h>
#include <PCU.h>
#include <parma.h>
#ifdef HAVE_SIMMETRIX
#include <gmi_sim.h>
#include <SimUtil.h>
#endif
#include <apfZoltan.h>
#include <pcu_util.h>
#include <cstdlib>
#include <iostream>

namespace {

const char* modelFile = 0;
const char* meshFile = 0;
const char* outFile = 0;
int partitionFactor = 1;  // number of parts to split into
int orig_parts = 1;  // number of parts meshFile is partitioned into

void freeMesh(apf::Mesh* m)
{
  m->destroyNative();
  apf::destroyMesh(m);
}

apf::Migration* getPlan(apf::Mesh* m)
{
  // TODO: make this a global splitter?
  apf::Splitter* splitter = apf::makeZoltanSplitter(
      m, apf::GRAPH, apf::PARTITION, false);
  apf::MeshTag* weights = Parma_WeighByMemory(m);
  apf::Migration* plan = splitter->split(weights, 1.05, partitionFactor);
  apf::removeTagFromDimension(m, weights, m->getDimension());
  m->destroyTag(weights);
  delete splitter;
  return plan;
}

bool switchToOriginals()
{
  int self = PCU_Comm_Self();

  int groupRank;
  int group;
  bool isOriginal;
  if ( self < orig_parts)
  {
    groupRank = self;
    group = 0;
    isOriginal = true;
  } else  // assign all other processes their own communicator
  {
    groupRank = 0;
    group = self - orig_parts + 1;
    isOriginal = false;
  }

  MPI_Comm groupComm;
  MPI_Comm_split(MPI_COMM_WORLD, group, groupRank, &groupComm);
  PCU_Switch_Comm(groupComm);

  return isOriginal;
}

void switchToAll()
{
  MPI_Comm prevComm = PCU_Get_Comm();
  PCU_Switch_Comm(MPI_COMM_WORLD);
  MPI_Comm_free(&prevComm);
  PCU_Barrier();
}

void getConfig(int argc, char** argv)
{
  if ( argc != 6 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <outMesh> <factor> <orig_parts>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  modelFile = argv[1];
  meshFile = argv[2];
  outFile = argv[3];
  partitionFactor = atoi(argv[4]);
  orig_parts = atoi(argv[5]);
  PCU_ALWAYS_ASSERT(partitionFactor <= PCU_Comm_Peers());
}

}

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
#ifdef HAVE_SIMMETRIX
  SimUtil_start();
  Sim_readLicenseFile(0);
  gmi_sim_start();
  gmi_register_sim();
#endif
  gmi_register_mesh();
  getConfig(argc,argv);

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  gmi_model* g = 0;
  g = gmi_load(modelFile);
  apf::Mesh2* m = 0;
  apf::Migration* plan = 0;
  bool isOriginal = switchToOriginals();

  // load mesh on first orig_parts processes
  if (isOriginal) {
    m = apf::loadMdsMesh(g, meshFile);
    plan = getPlan(m);
  }

  // now split it into all the processes
  switchToAll();
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0)
    std::cout << "finished switching to all" << std::endl;

  m = apf::repeatMdsMesh(m, g, plan, partitionFactor);

  if (myrank == 0)
    std::cout << "finished repeating Mds Mesh" << std::endl;

  Parma_PrintPtnStats(m, "");
  m->writeNative(outFile);
  freeMesh(m);
#ifdef HAVE_SIMMETRIX
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  SimUtil_stop();
#endif
  PCU_Comm_Free();
  MPI_Finalize();
}
