#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <iostream>
#include <stdbool.h>
//#include "funcs1.h"
#include "apfSBPShape.h"

void checkMesh(apf::Mesh* m);


// the program creates a cartesian mesh composed to triangles.
// The mesh is created by subdividing each rectangle into two triangles
// Arguments:
//   1: number of rectangles in x direction
//   2: number of rectangles in y direction
//   3: whether or not to perturb the node locations
//
// Conventions:
//   This mesh generator correctly classifies mesh entities onto geometric
//   entities.  The geometric entities are number starting from 0 in the 
//   bottom left corner, going counterclockwise.
//     
int main(int argc, char** argv)
{
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
  if (argc < 2 || argc > 4)
  {
    std::cerr << "Error: wrong number of arguments" << std::endl;
    return 1;
  }
  // Problem 3
  std::cout << "Problem 3 output: " << std::endl;
  int numElx = 0;  // number of elements in x direction
  sscanf(argv[1], "%d", &numElx);
  int numEly = 0;  // nuber of elements in y direction
  sscanf(argv[2], "%d", &numEly);
  int apply_pert = 0;
  if (argc == 4)
  {
    sscanf(argv[3], "%d", &apply_pert);
  }
  

  double xmin = 1.0;
  double ymin = 1.0;
  double xdist = 2;  // xmax - xmin
  double ydist = 2;  // ymax - ymin
  double x_spacing = xdist/numElx;  // spacing of el
  double y_spacing = ydist/numEly;
  double x_0 = xmin;  // x coordinate of lower left corner of current element
  double y_0 = ymin ; // y coordinate of lower left corner of current element
  double x_i = x_0;
  double y_i = y_0;
  double x_inner = x_i;
  double y_inner = y_i;
  double pert_fac = 10*M_PI;
  double pert_mag = 0.1;
  apf::MeshEntity* vertices[numElx+1][numEly+1];  // hold pointers to all vertices
  apf::MeshEntity* vertices_i[3];  // hold vertices for a particular element
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point
  int model_dim = 0;  // the dimension of the model entity to classify
                      // a newly created mesh entity on
  int model_tag = 0; // the tag of the model entity

  apf::ModelEntity* model_entity;  // the model entity itself
//  apf::FieldShape* linear2 = apf::getSBPQuadratic();
/*
  // for linear meshes
  apf::FieldShape* linear2 = apf::getLagrange(1);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);
*/


  std::cout << "Creating " << numElx << " by " << numEly << " mesh" << std::endl;

  // create all vertices
  for (int j = 0; j < (numEly+1); ++j)  // loop up columns (bottom to top)
  {
    for (int i = 0; i < (numElx+1); ++i) // loop over rows (left to right)
    {
       std::cout << "creating point at x = " << x_inner << " , y = " << y_inner << std::endl;
       coords_i[0] = x_inner;
       coords_i[1] = y_inner;
      
       // figure out what model entity to classify this vertex on 
       if ( i == 0 && j == 0 )
       {
         model_dim = 0;
         model_tag = 0;
       } else if ( i == (numElx) && j == 0 )
       {
         model_dim = 0;
         model_tag = 1;
       } else if ( i == (numElx) && j == (numEly) )
       {
         model_dim = 0;
         model_tag = 2;
       } else if ( i == 0 && j == (numEly) )
       {
         model_dim = 0;
         model_tag = 3;
       } else if ( j == 0)   // bottom row
       {
         model_dim = 1;
         model_tag = 0;
       } else if ( i == (numElx) )  // right column
       {
         model_dim = 1;
         model_tag = 1;
       } else if ( j == (numEly))  // top row
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
       x_inner = x_inner + x_spacing + pert_mag*sin(apply_pert*pert_fac*x_inner);
       y_inner = y_i + pert_mag*sin(apply_pert*pert_fac*x_inner);
    }

    y_i = y_i + y_spacing;  // increment y_i
    x_i = x_0;  // reset x_i to beginning of row

    y_inner = y_i;// + pert_mag*sin(apply_pert*pert_fac*y_i);
    x_inner = x_i; // + pert_mag*sin(apply_pert*pert_fac*x_i);
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
  for ( i = 0; i < numElx; ++i)
  {
    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i + 1][j];

    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // right column
  std::cout << "Creating right column edges" << std::endl;
  i = numElx;
  model_tag = 1;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEly; ++j)
  {
    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i][j + 1];
    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // top row
  std::cout << "Creating top row edges" << std::endl;
  j = numEly;
  model_tag = 2;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( i = 0; i < numElx; ++i)
  {
    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i + 1][j];
    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // left column
  std::cout << "Creating left column edges" << std::endl;
  i = 0;
  model_tag = 3;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEly; ++j)
  {
    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i][j + 1];
    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
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
  for (int j = 0; j < numEly; ++j)
  {
    for (int i = 0; i < numElx; ++i)
    {
      std::cout << "\ncreating element " << i << ", " << j << std::endl;
//      int el_num = i + numElx*j;
//      std::cout << "creating element i = " << i << " , j = " << j << std::endl;
      // get correct vertices
      vertices_i[0] = vertices[i][j];
      vertices_i[1] = vertices[i+1][j];
      vertices_i[2] = vertices[i][j+1];
//      std::cout << "Element " << el_num << " has verticies "; 
//      std::cout << i + ((numElx+1)*j) << " , " << i+1 + ((numElx+1)*j) << " , " << i+1 + ((numElx+1)*(j+1));
//      std::cout << " , " << i + ((numElx+1)*(j+1)) << std::endl;
//      std::cout << "assembled vertices" << std::endl;

      // create any edges that were not created before
      std::cout << "creating edges for first triangle" << std::endl;
       if ( i == 0 && j == 0 ) // bottom left
       {
         std::cout << "creating edges for bottom left corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElx-1) && j == 0 ) // bottom right
       {
         std::cout << "creating edges for bottom right corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElx-1) && j == (numEly-1) )  // top right
       {

         std::cout << "creating edges for top right corner" << std::endl;
         do_diag = true;
       } else if ( i == 0 && j == (numEly-1) )  // top left
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
      std::cout << "about to create first triangle" << std::endl;
      apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices_i);
//      m->createEntity(apf::Mesh::TRIANGLE, model_entity, vertices_i);

      // do other half of rectangle
      vertices_i[0] = vertices[i+1][j];
      vertices_i[1] = vertices[i+1][j+1];
      vertices_i[2] = vertices[i][j+1];

      std::cout << "Creating edges for second triangle" << std::endl;
      if ( i == 0 && j == 0 ) // bottom left
       {
         std::cout << "creating edges for bottom left corner" << std::endl;
         do_right = true;
         do_top = true;

       } else if ( i == (numElx-1) && j == 0 )  // bottom right
       {
         std::cout << "creating edges for bottom right corner" << std::endl;
         do_top = true;

       } else if ( i == (numElx-1) && j == (numEly-1) ) // top right
       {

         std::cout << "creating edges for top right corner" << std::endl;
         // do nothing

       } else if ( i == 0 && j == (numEly-1) )  // top left
       {
         std::cout << "creating edges for top left corner" << std::endl;
         do_right = true;

       } else if ( j == 0)   // bottom row
       {
         std::cout << "creating edges for bottom row" << std::endl;
         do_right = true;
         do_top = true;

       } else if ( i == (numElx-1) )  // right column
       {
         std::cout << "creating edges for right column" << std::endl;
         do_top = true;

       } else if ( j == (numEly-1))  // top row
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



       std::cout << "about to create second triangle" << std::endl;
//      m->createEntity(apf::Mesh::TRIANGLE, model_entity, vertices_i);
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
  checkMesh(m);
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
  // count nodes
//  int numNodes = m_entity_quad->countNodes();
//  std::cout << "number of nodes = " << numNodes << std::endl;
/*
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
  m->writeNative("./meshfiles/");
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

// check that the number of downward and upward adjacencies is correct
void checkMesh(apf::Mesh* m)
{
  // check verts -> edges
  apf::MeshIterator* it = m->begin(0);
  apf::MeshEntity* e;
//  apf::Up ups;  // upward adjacencies
  apf::Adjacent adj;  // multi level upward adjacencies
  apf::Downward down;
  int cnt = 0;  // current entity number
  int cnt_i = 0;  // number of adjacent entities

  while ( (e = m->iterate(it)) )
  {
    // check verts -> edges
    cnt_i = m->countUpward(e);
    if (cnt_i < 2)
    {
      std::cerr << "vertex " << cnt << " has too few edges" << std::endl;
    } else
    {
      std::cout << "vertex -> edge pass " << std::endl;
    }

    // check vert -> face
    m->getAdjacent(e, 2, adj);
    if (adj.size() > 6)
    {
      std::cerr << "vertex " << cnt << " has too many faces" << std::endl;
    } else
    {
      std::cout << "vertex -> face pass " << std::endl;
    }

    ++cnt;
  }

  it = m->begin(1);
  cnt = 0;

  while ( (e = m->iterate(it)) )
  {
    // check edge -> verts
    cnt_i = m->getDownward(e, 0, down);
    if ( cnt_i != 2)
    {
      std::cerr << "edge " << cnt << " has too few verts" << std::endl;
    } else
    {
      std::cout << "edge -> vert pass" << std::endl;
    }

    // check edge -> face
    m->getAdjacent(e, 2, adj);
    if (adj.size() > 2)
    {
      std::cerr << "edge " << cnt << " has too many faces" << std::endl;
    } else
    {
      std::cout << "edge -> face pass" << std::endl;
    }

    ++cnt;

  }

  it = m->begin(2);
  cnt = 0;

  while ( (e = m->iterate(it)) )
  {
    // check face -> edges
    cnt_i = m->getDownward(e, 1, down);
    if ( cnt_i != 3)
    {
      std::cerr << "face " << cnt << " has incorrect number of edges" << std::endl;
    } else
    {
      std::cout << "face -> edge pass" << std::endl;
    }

    // check face -> verts
    cnt_i = m->getDownward(e, 0, down);
    if (cnt_i != 3)
    {
      std::cerr << "face " << cnt << " has incorrect number of vertices" << std::endl;
    } else
    {
      std::cout << "face -> vert pass" << std::endl;
    }

    ++cnt;
  }


} // end function
