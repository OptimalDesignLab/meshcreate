#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <stdio.h>
#include <iostream>
#include <cmath>  // trig functions
#include <stdbool.h>
//#include "funcs1.h"
#include "apfSBPShape.h"

// create a mesh for a sector of a circle
// all angles measured in radians
// argv[1] = number of elements in r direction
// argv[2] = number of elements in theta direction
// See a1tri2 for description of geometric classification
int main(int argc, char** argv)
{

  std::cout << "argc = " << argc << std::endl;

  for (int i = 0; i < argc; ++i)
  {
    std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
    int val = 0;
    sscanf(argv[i], "%d", &val);
    std::cout << "val = " << val << std::endl;
  }

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, false);
/*
  apf::FieldShape* m_shape = m->getShape();
  const char shape_name[] = m_shape->getName();
  std::cout << "shape name = " << shape_name << std::endl;
*/
  // Problem 3
  std::cout << "Problem 3 output: " << std::endl;
  double pi = 3.14159265;
  int numElr = 0;  // number of elements in r direction
  sscanf(argv[1], "%d", &numElr);
  int numEltheta = 0;  // nuber of elements in theta direction
  sscanf(argv[2], "%d", &numEltheta);
  double r_range = 2.0;  // rmax - rmin
  double theta_range = M_PI/2.0;  // theta max - theta min
  double r_spacing = r_range/numElr;  // spacing of el radially
  double theta_spacing = theta_range/numEltheta;  // spacing of elements circumfrentially
  double r_0 = 1.0;  // r coordinate of lower left corner of current element
  double theta_0 = 0.0 ; // theta coordinate of lower left corner of current element
  double r_i = r_0;
  double theta_i = theta_0;

  apf::MeshEntity* vertices[numElr+1][numEltheta+1];  // hold pointers to all vertices
  apf::MeshEntity* vertices_i[3];  // hold vertices for a particular element
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point
//  apf::FieldShape* linear2 = apf::getSBPQuadratic();
  int model_dim = 0;  // the dimension of the model entity to classify
                      // a newly created mesh entity on
  int model_tag = 0; // the tag of the model entity

  apf::ModelEntity* model_entity;  // the model entity itself

/*
  // for linear meshes
  apf::FieldShape* linear2 = apf::getLagrange(1);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);
*/

  std::cout << "numElr = " << numElr << std::endl;
  std::cout << "numEltheta = " << numEltheta << " theta_spacing = " << theta_spacing << std::endl;


  // create all vertices
  for (int j = 0; j < (numEltheta+1); ++j)  // loop from inner radius to outer
  {
    for (int i = 0; i < (numElr+1); ++i) // loop from theta0 to theta max, counter clockwise
    {
 
       std::cout << "creating point at r = " << r_i << " , theta = " << theta_i << std::endl;
       double x_i = r_i*cos(theta_i);
       double y_i = r_i*sin(theta_i);
       coords_i[0] = x_i;
       coords_i[1] = y_i;
 
       // figure out what model entity to classify this vertex on 
       if ( i == 0 && j == 0 )
       {
         model_dim = 0;
         model_tag = 0;
       } else if ( i == (numElr) && j == 0 )
       {
         model_dim = 0;
         model_tag = 1;
       } else if ( i == (numElr) && j == (numEltheta) )
       {
         model_dim = 0;
         model_tag = 2;
       } else if ( i == 0 && j == (numEltheta) )
       {
         model_dim = 0;
         model_tag = 3;
       } else if ( j == 0)   // bottom row
       {
         model_dim = 1;
         model_tag = 0;
       } else if ( i == (numElr) )  // right column
       {
         model_dim = 1;
         model_tag = 1;
       } else if ( j == (numEltheta))  // top row
       {
         model_dim = 1;
         model_tag = 2;
       } else if ( i == 0 )  // left column
       {
         model_dim = 1;
         model_tag = 3;

       } else  // classify on face
       {
         model_dim = 2; 
         model_tag = 0;
       }

       model_entity = m->findModelEntity(model_dim, model_tag);
      
       vertices[i][j] = m->createVert(model_entity);
       m->setPoint(vertices[i][j], 0, coords_i);
//       theta_i = theta_i + theta_spacing;
       r_i = r_i + r_spacing;
    }
//    theta_i = theta_0;
    theta_i = theta_i + theta_spacing;  // increment y_i
    r_i = r_0;  // reset x_i to beginning of row
//    r_i = r_i + r_spacing;
  }  

  // create the boundary mesh edges, classifying them correctly
  int j = 0;  // row index (bottom to top)
  int i = 0;  // column index (left to right)
  model_dim = 1;
  apf::Downward edge_verts;

  // bottom row
  std::cout << "Creating bottom row edges" << std::endl;
  model_tag = 0;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( i = 0; i < numElr; ++i)
  {
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i + 1][j];

    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // right column
  std::cout << "Creating right column edges" << std::endl;
  i = numElr;
  model_tag = 1;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEltheta; ++j)
  {
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i][j + 1];
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // top row
  std::cout << "Creating top row edges" << std::endl;
  j = numEltheta;
  model_tag = 2;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( i = 0; i < numElr; ++i)
  {
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i + 1][j];
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // left column
  std::cout << "Creating left column edges" << std::endl;
  i = 0;
  model_tag = 3;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEltheta; ++j)
  {
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i][j + 1];
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }



  // classify all MeshEntities on geometric face
  model_dim = 2;
  model_tag = 0;
  model_entity = m->findModelEntity(model_dim, model_tag);
  
  // bools to control which edges get created
  bool do_bottom = false;
  bool do_diag = false;
  bool do_left = false;
  bool do_right = false;
  bool do_top = false;

  // build element from verticies
  for (int j = 0; j < numEltheta; ++j)
  {
    for (int i = 0; i < numElr; ++i)
    {
      int el_num = i + numEltheta*j;
//      std::cout << "creating element i = " << i << " , j = " << j << std::endl;
      // get correct vertices
      vertices_i[0] = vertices[i][j];
      vertices_i[1] = vertices[i+1][j];
      vertices_i[2] = vertices[i][j+1];
      std::cout << "Element " << el_num << " has verticies "; 
      std::cout << i + ((numElr+1)*j) << " , " << i+1 + ((numElr+1)*j) << " , " << i+1 + ((numElr+1)*(j+1));
      std::cout << " , " << i + ((numElr+1)*(j+1)) << std::endl;
//      std::cout << "assembled vertices" << std::endl;
//
      // create any edges that were not created before
      std::cout << "creating edges for first triangle" << std::endl;
       if ( i == 0 && j == 0 ) // bottom left
       {
         std::cout << "creating edges for bottom left corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElr-1) && j == 0 ) // bottom right
       {
         std::cout << "creating edges for bottom right corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElr-1) && j == (numEltheta-1) )  // top right
       {

         std::cout << "creating edges for top right corner" << std::endl;
         do_diag = true;
       } else if ( i == 0 && j == (numEltheta-1) )  // top left
       {
         std::cout << "creating edges for top left corner" << std::endl;
         do_diag = true;
       } else  // interior
       {

         std::cout << "creating interior mesh edges" << std::endl;
         do_diag = true;
       }

       // actually create the vertices
       if (do_bottom)
       {
         std::cout << "creating bottom edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[1];
         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_diag)
       {
         std::cout << "creating diagonal edge" << std::endl;
         edge_verts[0] = vertices_i[1];
         edge_verts[1] = vertices_i[2];
         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_left)
       {
         std::cout << "creating left edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[2];
         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }


      // reset flags
      do_bottom = false;
      do_diag = false;
      do_left = false;
      do_right = false;
      do_top = false;




      // counterclockwise ordering
      apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices_i);
      vertices_i[0] = vertices[i+1][j];
     vertices_i[1] = vertices[i+1][j+1];
      vertices_i[2] = vertices[i][j+1];

      std::cout << "Creating edges for second triangle" << std::endl;
      if ( i == 0 && j == 0 && (numElr > 1 && numEltheta > 1) ) // bottom left
       {
         std::cout << "creating edges for bottom left corner" << std::endl;
         do_right = true;
         do_top = true;
       // bottom left, 1 element in Y
       } else if (i == 0 && j == 0 && (numElr > 1 && numEltheta == 1) )
       {
         
         std::cout << "creating edges for bottom left corner with 1 element";
         std::cout << " in y direction" << std::endl;
         do_right = true;
       // bottom left, 1 element in X
       } else if ( i == 0 && j == 0 && (numElr == 1 && numEltheta > 1) )
       {
         std::cout << "creating edges for bottom left corner with 1 element";
         std::cout << " in x direction" << std::endl;
         do_top = true;
       // bottom left, 1 element in X and Y
       } else if ( i == 0 && j == 0 && (numElr == 1 && numEltheta == 1) )
       {

         std::cout << "creating edges for bottom left corner with 1 element";
         std::cout << " in x and y  direction" << std::endl;


       } else if ( i == (numElr-1) && j == 0 && numEltheta > 1 )  // bottom right
       {
         std::cout << "creating edges for bottom right corner" << std::endl;
         do_top = true;

       } else if ( i == (numElr-1) && j == 0 && numEltheta == 1)
       {
         std::cout << "creating edges for bottom right corner with 1 element";
         std::cout << " in the y direction" << std::endl;



       } else if ( i == (numElr-1) && j == (numEltheta-1) ) // top right
       {

         std::cout << "creating edges for top right corner" << std::endl;
         // do nothing

       } else if ( i == 0 && j == (numEltheta-1) )  // top left
       {
         std::cout << "creating edges for top left corner" << std::endl;
         do_right = true;

       } else if ( j == 0)   // bottom row
       {
         std::cout << "creating edges for bottom row" << std::endl;
         do_right = true;
         do_top = true;

       } else if ( i == (numElr-1) )  // right column
       {
         std::cout << "creating edges for right column" << std::endl;
         do_top = true;

       } else if ( j == (numEltheta-1))  // top row
       {
         std::cout << "creating edges for top row" << std::endl;
         do_right = true;

       } else if ( i == 0 )  // left column
       {
         std::cout << "creating edges for left column" << std::endl;
         do_right = true;
         do_top = true;

       } else  // do both
       {
         std::cout << "creating interior mesh edges" << std::endl;
         do_right = true;
         do_top = true;
       }

       if (do_right)
       {
         std::cout << "creating right edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[1];
         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_top)
       {
         std::cout << "creating top edge" << std::endl;
         edge_verts[0] = vertices_i[1];
         edge_verts[1] = vertices_i[2];
         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }

      // reset flags
      do_bottom = false;
      do_diag = false;
      do_left = false;
      do_right = false;
      do_top = false;


      apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices_i);

/*
      vertices_i[1] = vertices[i][[j+1];
      apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vertices_i);
*/
     }
  }
/*
//  apf::FieldShape* linear2 = apf::getSBPQuadratic();
    apf::FieldShape* linear2 = apf::getLagrange(1);
//  const char shape_name[] = linear2->getName();
  std::cout << "shape name = " << linear2->getName() << std::endl;
  m->init(linear2);
*/
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

 // for quadratic meshes
//  apf::FieldShape* linear2 = apf::getSerendipity();
    apf::FieldShape* linear2 = apf::getSBPShape(1);
//  apf::FieldShape* linear2 = apf::getLagrange(2);
  apf::changeMeshShape(m, linear2, true);  // last argument should be true for second order

  std::cout << "changed mesh shape" << std::endl;
  apf::FieldShape* m_shape = m->getShape();
//  const char shape_name[] = m_shape->getName();
  std::cout << "mesh shape name = " << m_shape->getName() << std::endl;

//  apf::EntityShape* m_entity_quad = m_shape->getEntityShape(apf::Mesh::QUAD);
/*
  // get values
  apf::Vector3 xi(-0.25, -0.25, 0);
  apf::NewArray<double> vals;
  m_entity_quad->getValues(xi, vals);
  std::cout << "values at (-0.25. -0.25, 0) = " << vals[0] << " , " << vals[1] << " , " << vals[2] << " , " << vals[3] << std::endl;

  // get gradients
  apf::NewArray<apf::Vector3> vals2;
  m_entity_quad->getLocalGradients(xi, vals2);
  std::cout << "gradients at (-0.25. -0.25, 0) = " << vals2[0] << " , " << vals2[1] << " , " << vals2[2] << " , " << vals2[3] << std::endl;
*/

/*  
  // count nodes
  int numNodes = m_entity_quad->countNodes();
  std::cout << "number of nodes = " << numNodes << std::endl;

  // check number of nodes for each type of entity
  bool nodecnt[4];
  nodecnt[0] = m_shape->countNodesOn(apf::Mesh::VERTEX);
  nodecnt[1] = m_shape->countNodesOn(apf::Mesh::EDGE);
  nodecnt[2] = m_shape->countNodesOn(apf::Mesh::TRIANGLE);
  nodecnt[3] = m_shape->countNodesOn(apf::Mesh::TET);
//  nodecnt[3] = m_shape->countNodesOn(apf::Mesh::QUAD);
  std::cout << "nodecounts: " << nodecnt[0] << " , " << nodecnt[1] << " , " << nodecnt[2] << " , " << nodecnt[3] << std::endl;
*/

  // write output and clean up
  apf::writeVtkFiles("outTri", m);
  m->writeNative("/users/creanj/meshcreate/meshfiles/");
/*
  apf::MeshIterator* it = m->begin(2);
  apf::MeshEntity* e = m->iterate(it);
  apf::MeshElement* e_el = apf::createMeshElement(m, e);
  int numI = apf::countIntPoints(e_el, 5);
  std::cout << numI << " integration points required" << std::endl;
*/

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}

