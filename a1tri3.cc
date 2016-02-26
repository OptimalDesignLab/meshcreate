#include <iostream>
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
  apf::Downward edge_verts;
  int model_dim = 0;  // the dimension of the model entity to classify
                      // a newly created mesh entity on
  apf::ModelEntity* model_entity;  // the model entity itself


  // create the vertices
  std::cout << "creating vertices" << std::endl;
  model_entity = m->findModelEntity(model_dim, 0);
  apf::MeshEntity* vert1 = m->createVert(model_entity);
  m->setPoint(vert1, 0, coords1);

  model_entity = m->findModelEntity(model_dim, 1);
  apf::MeshEntity* vert2 = m->createVert(model_entity);
  m->setPoint(vert2, 0, coords2);

  model_entity = m->findModelEntity(model_dim, 2);
  apf::MeshEntity* vert3 = m->createVert(model_entity);
  m->setPoint(vert3, 0, coords3);

  model_entity = m->findModelEntity(model_dim, 3);
  apf::MeshEntity* vert4 = m->createVert(model_entity);
  m->setPoint(vert4, 0, coords4);

  apf::MeshEntity* vertices1[] = {vert1, vert2, vert4};
  apf::MeshEntity* vertices2[] = {vert2, vert3, vert4};

  // create the edges
  std::cout << "creating edges" << std::endl;
  model_dim = 1;

  edge_verts[0] = vertices1[0];
  edge_verts[1] = vertices1[1];
  model_entity = m->findModelEntity(model_dim, 0);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices1[1];
  edge_verts[1] = vertices1[2];
  model_entity = m->findModelEntity(2, 0);  // classified on geometric face
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices1[2];
  edge_verts[1] = vertices1[0];
  model_entity = m->findModelEntity(model_dim, 1);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices2[0];
  edge_verts[1] = vertices2[1];
  model_entity = m->findModelEntity(model_dim, 3);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices2[1];
  edge_verts[1] = vertices2[2];
  model_entity = m->findModelEntity(model_dim, 4);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

/*
  // create faces
  model_dim = 2;
  model_entity = m->findModelEntity(model_dim, 0);
  std::cout << "Creating first face" << std::endl;
  m->createEntity(apf::Mesh::TRIANGLE, model_entity, vertices1);
  std::cout << "creating second face" << std::endl;
  m->createEntity(apf::Mesh::TRIANGLE, model_entity, vertices2);
*/



  // create the elements
  model_dim = 2;
  model_entity = m->findModelEntity(model_dim, 0);

  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices1);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices2);

  // build, verify  mesh
/*
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m);
  std::cout << "finished deriving model" << std::endl;
*/
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
  m->writeNative("/users/creanj/meshcreate/meshfiles/");

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

