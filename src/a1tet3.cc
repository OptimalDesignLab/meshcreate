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
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 3, false);
  
  apf::Vector3 coords1 (0, 0, 0);
  apf::Vector3 coords2 (1, 0, 0);
  apf::Vector3 coords3 (0, 1, 0);
  apf::Vector3 coords4 (0, 1, 1);
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

  apf::MeshEntity* vertices[] = {vert1, vert2, vert3, vert4};

  // create the edges
  std::cout << "creating edges" << std::endl;
  model_dim = 1;

  edge_verts[0] = vertices[0];
  edge_verts[1] = vertices[1];
  model_entity = m->findModelEntity(model_dim, 0);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[1];
  edge_verts[1] = vertices[2];
  model_entity = m->findModelEntity(model_dim, 1);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[2];
  edge_verts[1] = vertices[0];
  model_entity = m->findModelEntity(model_dim, 2);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[0];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 3);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[1];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 4);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[2];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 5);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);


  // create faces
  std::cout << "creating faces" << std::endl;
  apf::MeshEntity* face_verts[3];
  model_dim = 3;

  face_verts[0] = vertices[0];
  face_verts[1] = vertices[1];
  face_verts[2] = vertices[2];
  model_entity = m->findModelEntity(model_dim, 0);
  std::cout << "about to create first face" << std::endl;
  std::cout << "model_dim = " << model_dim << std::endl;
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);
  std::cout << "finished creating first face" << std::endl;

  model_dim = 2;
  face_verts[0] = vertices[0];
  face_verts[1] = vertices[1];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 1);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

  face_verts[0] = vertices[1];
  face_verts[1] = vertices[2];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 2);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

  face_verts[0] = vertices[0];
  face_verts[1] = vertices[2];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 3);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);


  // create the elements
  model_dim = 3;
  model_entity = m->findModelEntity(model_dim, 0);
  apf::buildElement(m, model_entity, apf::Mesh::TET, vertices);

  // create a second element vert
  model_dim = 0;
  apf::Vector3 coords5 (0, 1, -1);
  model_entity = m->findModelEntity(model_dim, 4);
  vert4 = m->createVert(model_entity);
  m->setPoint(vert4, 0, coords5);

  vertices[3] = vert4;

  // create edges
  model_dim = 1;
  edge_verts[0] = vertices[0];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 6);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[1];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 7);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  edge_verts[0] = vertices[2];
  edge_verts[1] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 8);
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

  // create faces
  model_dim = 2;
  face_verts[0] = vertices[0];
  face_verts[1] = vertices[1];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 4);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

  face_verts[0] = vertices[1];
  face_verts[1] = vertices[2];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 5);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

  face_verts[0] = vertices[0];
  face_verts[1] = vertices[2];
  face_verts[2] = vertices[3];
  model_entity = m->findModelEntity(model_dim, 6);
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

  // create the elements
  model_dim = 3;
  model_entity = m->findModelEntity(model_dim, 0);
  apf::buildElement(m, model_entity, apf::Mesh::TET, vertices);



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
  apf::writeVtkFiles("outTet", m);
  m->writeNative("/users/creanj/meshcreate/meshfiles/");

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

