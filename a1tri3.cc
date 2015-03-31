#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

int main(int argc, char** argv)
{
  // init
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, false);
/*    
  // apply field shape
  apf::FieldShape* linear2 = apf::getLagrange(1);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);
*/

  // create a single quad
/*
  apf::Vector3 coords1 (0.0, 0.0, 0.0);
  apf::Vector3 coords2 ( 1.0, 0, 0.0);
  apf::Vector3 coords3 ( 1.0, 1.0, 0.0);
  apf::Vector3 coords4 (0.0, 1.0, 0.0);
*/
  apf::Vector3 coords1 (-1.0, -1.0, 0.0);
  apf::Vector3 coords2 ( 1.0, -1.0, 0.0);
  apf::Vector3 coords3 ( 1.0, 1.0, 0.0);
  apf::Vector3 coords4 (-1.0, 1.0, 0.0);


  apf::MeshEntity* vert1 = m->createVert(0);
  m->setPoint(vert1, 0, coords1);

  apf::MeshEntity* vert2 = m->createVert(0);
  m->setPoint(vert2, 0, coords2);

  apf::MeshEntity* vert3 = m->createVert(0);
  m->setPoint(vert3, 0, coords3);

  apf::MeshEntity* vert4 = m->createVert(0);
  m->setPoint(vert4, 0, coords4);

  apf::MeshEntity* vertices1[] = {vert1, vert2, vert4};
  apf::MeshEntity* vertices2[] = {vert2, vert3, vert4};

  apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vertices1);
  apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vertices2);
/*
  apf::FieldShape* linear2 = apf::getLagrange(1);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);
*/


  // build, verify  mesh
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m);
  std::cout << "finished deriving model" << std::endl;
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  m->verify();
  std::cout << "verified" << std::endl;

  apf::FieldShape* linear2 = apf::getLagrange(2);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);


  // check that the mesh actually has a FieldShape with the right name
  apf::FieldShape* m_shape = m->getShape();
  std::cout << "mesh shape name = " << m_shape->getName() << std::endl;

  
  // write output and clean up
  apf::writeVtkFiles("outQuad", m);
  m->writeNative("/users/creanj/atest/meshfiles/");

  // build a mesh element
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e = m->iterate(it);
  apf::MeshElement* e_el = apf::createMeshElement(m, e);
  int numI = apf::countIntPoints(e_el, 2);  // segfault
  std::cout << numI << " integration points required" << std::endl;


  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

