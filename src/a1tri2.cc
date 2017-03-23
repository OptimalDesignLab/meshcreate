#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <cstdlib>
#include <iostream>
#include <stdbool.h>
#include <limits.h>
//#include "funcs1.h"
//#include "apfSBPShape.h"

struct _Periodic {
  bool x;
  bool y;
};
typedef struct _Periodic Periodic;

struct _Counts {
  int numElx;
  int numEly;
};
typedef struct _Counts Counts;

struct _DomainSize {
  double xmin;
  double xmax;
  double ymin;
  double ymax;
  int numElx;
  int numEly;
};
typedef struct _DomainSize DomainSize;

void checkMesh(apf::Mesh* m);
void setMatches(apf::Mesh2*m, apf::MeshEntity*** verts, Periodic periodic, Counts counts);
apf::MeshEntity* getEdge(apf::Mesh* m, apf::MeshEntity* v1, apf::MeshEntity* v2);
void printMatches(apf::Mesh* m, Periodic periodic);

/* functions for applying coordinate mapping */
void remapCoordinates(apf::Mesh2* m, DomainSize domainsize);
void nondimensionalize(DomainSize domainsize, apf::Vector3& coords);
void redimensionalize(DomainSize domainsize, apf::Vector3& coords);
void mapFunction(DomainSize domainsize, apf::Vector3& coords);

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
  
  long numElx_l = numElx;
  long numEly_l = numEly;
  long nvert = (numElx_l+1)*(numEly_l+1);
  long nedges = numElx_l*(numEly_l+1) + numEly_l*(numElx_l+1) + numElx_l*numEly_l;
  long nel = numElx_l*numEly_l*2;

  std::cout << "expected entity counts: " << std::endl;
  std::cout << "  Vertices: " << nvert << std::endl;
  std::cout << "  Edges: " << nedges << std::endl;
  std::cout << "  Elements: " << nel << std::endl;

  if (nvert > INT_MAX)
  {
    std::cerr << "Error: number of vertices will overflow 32-bit integer" << std::endl;
    return 1;
  }

  if (nedges > INT_MAX)
  {
    std::cerr << "Error: number of edges will overflow 32-bit integer" << std::endl;
    return 1;
  }
  if (nel > INT_MAX)
  {
    std::cerr << "Error: number of elements will overflow 32-bit integer" << std::endl;
    return 1;
  }

  Counts counts = {numElx, numEly};
  bool xperiodic = true;  // make x direction periodic
  bool yperiodic = true; // make y direction periodic
  // making x direction direction periodic mean setting edges along y 
  // axis to match, hence the reversal
  Periodic periodic = {yperiodic, xperiodic};

  int coord_order = 1;  // coordinate field polynomial order

  double xmin = 0;
  double ymin = -5;
  double xdist = 20;  // xmax - xmin
  double ydist = 10;  // ymax - ymin
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

  DomainSize domainsize = {xmin, xmin + xdist, ymin, ymin + ydist, numElx, numEly};

  bool isMatched = false;
  if (periodic.x || periodic.y)
    isMatched = true;

  std::cout << "about to allocate memory" << std::endl;
  apf::MeshEntity*** vertices = (apf::MeshEntity***)std::calloc(numElx+1, sizeof(apf::MeshEntity**));
  for (int i = 0; i < (numElx+1); ++i)
  {
    vertices[i] = (apf::MeshEntity**)std::calloc(numEly+1, sizeof(apf::MeshEntity*));
  }
//  apf::MeshEntity* vertices[numElx+1][numEly+1];  // hold pointers to all vertices
  std::cout << "finished allocating memory" << std::endl;
  apf::MeshEntity* vertices_i[3];  // hold vertices for a particular element
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point
  int model_dim = 0;  // the dimension of the model entity to classify
                      // a newly created mesh entity on
  int model_tag = 0; // the tag of the model entity

  apf::ModelEntity* model_entity;  // the model entity itself

  std::cout << "Creating " << numElx << " by " << numEly << " mesh" << std::endl;

  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, isMatched);

  // create all vertices
  for (int j = 0; j < (numEly+1); ++j)  // loop up columns (bottom to top)
  {
    for (int i = 0; i < (numElx+1); ++i) // loop over rows (left to right)
    {
//       std::cout << "creating point at x = " << x_inner << " , y = " << y_inner << std::endl;
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
//  std::cout << "Creating bottom row edges" << std::endl;
  model_tag = 0;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( i = 0; i < numElx; ++i)
  {
//    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i + 1][j];

//    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // right column
//  std::cout << "Creating right column edges" << std::endl;
  i = numElx;
  model_tag = 1;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEly; ++j)
  {
//    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[0] = vertices[i][j];
    edge_verts[1] = vertices[i][j + 1];
//    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // top row
//  std::cout << "Creating top row edges" << std::endl;
  j = numEly;
  model_tag = 2;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( i = 0; i < numElx; ++i)
  {
//    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i + 1][j];
//    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
    m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
  }

  // left column
//  std::cout << "Creating left column edges" << std::endl;
  i = 0;
  model_tag = 3;
  model_entity = m->findModelEntity(model_dim, model_tag);
  for ( j = 0; j < numEly; ++j)
  {
//    std::cout << " creating edge at i = " << i << " j = " << j << std::endl;
    edge_verts[1] = vertices[i][j];
    edge_verts[0] = vertices[i][j + 1];
//    std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
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
//      std::cout << "\ncreating element " << i << ", " << j << std::endl;
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
//      std::cout << "creating edges for first triangle" << std::endl;
       if ( i == 0 && j == 0 ) // bottom left
       {
//         std::cout << "creating edges for bottom left corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElx-1) && j == 0 ) // bottom right many elements
       {
//         std::cout << "creating edges for bottom right corner" << std::endl;
         do_diag = true;
       } else if ( i == (numElx-1) && j == (numEly-1) )  // top right
       {
//         std::cout << "creating edges for top right corner" << std::endl;
         do_diag = true;
       } else if ( i == 0 && j == (numEly-1) )  // top left
       {
//         std::cout << "creating edges for top left corner" << std::endl;
         do_diag = true;
       } else  // interior
       {

//         std::cout << "creating interior mesh edges" << std::endl;
         do_diag = true;
       }

       // actually create the vertices
       if (do_bottom)
       {
//         std::cout << "creating bottom edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[1];
//         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_diag)
       {
//         std::cout << "creating diagonal edge" << std::endl;
         edge_verts[0] = vertices_i[1];
         edge_verts[1] = vertices_i[2];
//         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_left)
       {
//         std::cout << "creating left edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[2];
//         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }


      // reset flags
      do_bottom = false;
      do_diag = false;
      do_left = false;
      do_right = false;
      do_top = false;



      // counterclockwise ordering
//      std::cout << "about to create first triangle" << std::endl;
      apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices_i);
//      m->createEntity(apf::Mesh::TRIANGLE, model_entity, vertices_i);

      // do other half of rectangle
      vertices_i[0] = vertices[i+1][j];
      vertices_i[1] = vertices[i+1][j+1];
      vertices_i[2] = vertices[i][j+1];
/*
      std::cout << "numElx = " << numElx << std::endl;
      std::cout << "numEly = " << numEly << std::endl;
      std::cout << "i = " << i << std::endl;
      std::cout << "j = " << j << std::endl;
      std::cout << "Creating edges for second triangle" << std::endl;
*/
      // bottom left, many elements
      if ( i == 0 && j == 0 && (numElx > 1 && numEly > 1) ) 
       {
/*
         std::cout << " i = " << i << ", j = " << j << ", numElx = " << numElx;
         std::cout << " numEly = " << numEly << std::endl;
         std::cout << "creating edges for bottom left corner" << std::endl;
*/
         do_right = true;
         do_top = true;
       // bottom left, 1 element in Y
       } else if (i == 0 && j == 0 && (numElx > 1 && numEly == 1) )
       {         
//         std::cout << "creating edges for bottom left corner with 1 element";
//         std::cout << " in y direction" << std::endl;
         do_right = true;
       // bottom left, 1 element in X
       } else if ( i == 0 && j == 0 && (numElx == 1 && numEly > 1) )
       {
//         std::cout << "creating edges for bottom left corner with 1 element";
//         std::cout << " in x direction" << std::endl;
         do_top = true;
       // bottom left, 1 element in X and Y
       } else if ( i == 0 && j == 0 && (numElx == 1 && numEly == 1) )
       {
//         std::cout << "creating edges for bottom left corner with 1 element";
//         std::cout << " in x and y  direction" << std::endl;
       } else if ( i == (numElx-1) && j == 0 && numEly > 1 )  // bottom right
       {
//         std::cout << "creating edges for bottom right corner" << std::endl;
         do_top = true;

       } else if ( i == (numElx-1) && j == 0 && numEly == 1)
       {
//         std::cout << "creating edges for bottom right corner with 1 element";
//         std::cout << " in the y direction" << std::endl;
       } else if ( i == (numElx-1) && j == (numEly-1) ) // top right
       {
//         std::cout << "creating edges for top right corner" << std::endl;
         // do nothing

       } else if ( i == 0 && j == (numEly-1) )  // top left
       {
//         std::cout << "creating edges for top left corner" << std::endl;
         do_right = true;

       } else if ( j == 0)   // bottom row
       {
//         std::cout << "creating edges for bottom row" << std::endl;
         do_right = true;
         do_top = true;

       } else if ( i == (numElx-1) )  // right column
       {
//         std::cout << "creating edges for right column" << std::endl;
         do_top = true;
       } else if ( j == (numEly-1))  // top row
       {
//         std::cout << "creating edges for top row" << std::endl;
         do_right = true;

       } else if ( i == 0 )  // left column
       {
//         std::cout << "creating edges for left column" << std::endl;
         do_right = true;
         do_top = true;

       } else  // do both
       {
//         std::cout << "creating interior mesh edges" << std::endl;
         do_right = true;
         do_top = true;
       }

       if (do_right)
       {
//         std::cout << "creating right edge" << std::endl;
         edge_verts[0] = vertices_i[0];
         edge_verts[1] = vertices_i[1];
//         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }
       if (do_top)
       {
//         std::cout << "creating top edge" << std::endl;
         edge_verts[0] = vertices_i[1];
         edge_verts[1] = vertices_i[2];
//         std::cout << " with vertices " << edge_verts[0] << " " << edge_verts[1] << std::endl;
         m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);
       }

      // reset flags
      do_bottom = false;
      do_diag = false;
      do_left = false;
      do_right = false;
      do_top = false;

      apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices_i);

     }
  }

  // set periodic boundaries if needed
  setMatches(m, vertices, periodic, counts);
 
  // make the coordinate field the requested order
  apf::FieldShape* linear2 = apf::getLagrange(coord_order);
  apf::changeMeshShape(m, linear2, true);  // last argument should be true for second order
  std::cout << "remapping coordinates" << std::endl;
  remapCoordinates(m, domainsize);

  std::cout << "changed mesh shape" << std::endl;
  apf::FieldShape* m_shape = m->getShape();
  std::cout << "mesh shape name = " << m_shape->getName() << std::endl;

 
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  checkMesh(m);
  m->verify();
  std::cout << "verified" << std::endl;


//  printMatches(m, periodic);
  // write output and clean up
  apf::writeVtkFiles("outTri", m);
  m->writeNative("./meshfiles/abc.smb");

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
//      std::cout << "vertex -> edge pass " << std::endl;
    }

    // check vert -> face
    m->getAdjacent(e, 2, adj);
    if (adj.size() > 6)
    {
      std::cerr << "vertex " << cnt << " has too many faces" << std::endl;
    } else
    {
//      std::cout << "vertex -> face pass " << std::endl;
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
//      std::cout << "edge -> vert pass" << std::endl;
    }

    // check edge -> face
    m->getAdjacent(e, 2, adj);
    if (adj.size() > 2)
    {
      std::cerr << "edge " << cnt << " has too many faces" << std::endl;
    } else
    {
//      std::cout << "edge -> face pass" << std::endl;
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
//      std::cout << "face -> edge pass" << std::endl;
    }

    // check face -> verts
    cnt_i = m->getDownward(e, 0, down);
    if (cnt_i != 3)
    {
      std::cerr << "face " << cnt << " has incorrect number of vertices" << std::endl;
    } else
    {
//      std::cout << "face -> vert pass" << std::endl;
    }

    ++cnt;
  }


} // end function

// set matched entities for periodic boundaries
void setMatches(apf::Mesh2*m, apf::MeshEntity*** verts, Periodic periodic, Counts counts)
{
  apf::MeshEntity* e1;
  apf::MeshEntity* e2;
  apf::MeshEntity* edge1;
  apf::MeshEntity* edge2;
  if (periodic.x)
  {
    std::cout << "setting x direction matches" << std::endl;
    // set matching vertices
    for (int i = 0; i < (counts.numElx+1); ++i)
    {
      e1 = verts[i][0];
      e2 = verts[i][counts.numEly];
      m->addMatch(e1, 0, e2);
      m->addMatch(e2, 0, e1);
    }

    // now do edges
    for (int i = 0; i < counts.numElx; ++i)
    {
      e1 = verts[i][0];
      e2 = verts[i+1][0];
      edge1 = getEdge(m, e1, e2);

      e1 = verts[i][counts.numEly];
      e2 = verts[i+1][counts.numEly];
      edge2 = getEdge(m, e1, e2);
      m->addMatch(edge1, 0, edge2);
      m->addMatch(edge2, 0, edge1);
    }

  }
 
  if (periodic.y)
  {
    std::cout << "setting y direction matches" << std::endl;
    for (int i = 0; i < (counts.numEly+1); ++i)
    {
      e1 = verts[0][i];
      e2 = verts[counts.numElx][i];
      m->addMatch(e1, 0, e2);
      m->addMatch(e2, 0, e1);
    }

    for (int i = 0; i < counts.numEly; ++i)
    {
      e1 = verts[0][i];
      e2 = verts[0][i+1];
      edge1 = getEdge(m, e1, e2);

      e1 = verts[counts.numElx][i];
      e2 = verts[counts.numElx][i+1];
      edge2 = getEdge(m, e1, e2);
      m->addMatch(edge1, 0, edge2);
      m->addMatch(edge2, 0, edge1);
    }
  }

  if (periodic.x && periodic.y)
  {
    e1 = verts[0][0];
    e2 = verts[counts.numElx][counts.numEly];
    m->addMatch(e1, 0, e2);
    m->addMatch(e2, 0, e1);

    e1 = verts[counts.numElx][0];
    e2 = verts[0][counts.numEly];
    m->addMatch(e1, 0, e2);
    m->addMatch(e2, 0, e1);
  }

} // function setMatches

// get edge defined by 2 vertices
apf::MeshEntity* getEdge(apf::Mesh* m, apf::MeshEntity* v1, apf::MeshEntity* v2)
{
  static apf::Up up1; static apf::Up up2;
  // get edge common to v1 and v2
  m->getUp(v1, up1);
  m->getUp(v2, up2);

  for (int i=0; i < up1.n; ++i)
    for (int j=0; j < up2.n; ++j)
    {
      if (up1.e[i] == up2.e[j])
        return up1.e[i];
    }

  return NULL;

}  // function getEdge

void printMatches(apf::Mesh* m, Periodic periodic)
{
  std::cout << std::endl;
  apf::Sharing* shr = getSharing(m);
  if (periodic.x)
  {
    for (int dim=0; dim < 3; ++dim)
    {
      std::cout << "\nchecking entities of dimension " << dim << std::endl;
      apf::MeshIterator* it = m->begin(dim);
      apf::MeshEntity* e;
      while ( (e = m->iterate(it)) )
      {
        apf::Matches matches;
        m->getMatches(e, matches);
        for (apf::Matches::iterator it = matches.begin(); it != matches.end(); ++it)
          std::cout << "entity " << e << " is matched to " << it->entity << std::endl;
        if (dim == 0)
          std::cout << "is entity owned: " << shr->isOwned(e) << std::endl;
      }
    }
  }
} // function printMatches
     
// this function applies a mapping to the coordinate field of the form:
//    x -> x* -> f(x*) -> x'
// where x is the original coordinate, x* is a nondimensionalized coordinate in
// the range [0, 1], f() is a function that maps x* to another value in
// the range [0, 1], and x' converts the (nondimensional) output of f() back
// to a dimensional value
// An inner function is called that determines f().  This function can be
// modified by the user, but remapCoordinates itself should not be modified.
void remapCoordinates(apf::Mesh2* m, DomainSize domainsize)
{
  apf::FieldShape* fshape = m->getShape();
  apf::MeshIterator* it;
  apf::MeshEntity* e;
  int maxdim = m->getDimension();
  apf::Vector3 coords;

//  int dim = 1;  // only do edges

  for (int dim = 0; dim <= maxdim; dim++)
  {
    if ( !fshape->hasNodesIn(dim) )
      continue;

    // if we get here, there are nodes in this dimension
    it = m->begin(dim);
    while ( (e = m->iterate(it)) )
    {
      int nnodes = fshape->countNodesOn(m->getType(e));

      for (int node = 0; node < nnodes; node++)
      {
        m->getPoint(e, node, coords);
        nondimensionalize(domainsize, coords);
        mapFunction(domainsize, coords);
        redimensionalize(domainsize, coords);
        m->setPoint(e, node, coords);
      }
    }  // end while loop

  }  // end for dim


}  // function remapCoordinates

void nondimensionalize(DomainSize domainsize, apf::Vector3& coords)
{
  double x = coords[0];
  double y = coords[1];

  const double delta_x = domainsize.xmax - domainsize.xmin;
  const double delta_y = domainsize.ymax - domainsize.ymin;

  x = (x - domainsize.xmin)/delta_x;
  y = (y - domainsize.ymin)/delta_y;

  coords[0] = x;
  coords[1] = y;

}  // function nondimensionalize

void redimensionalize(DomainSize domainsize, apf::Vector3& coords)
{
  double x = coords[0];
  double y = coords[1];

  const double delta_x = domainsize.xmax - domainsize.xmin;
  const double delta_y = domainsize.ymax - domainsize.ymin;

  x = x*delta_x + domainsize.xmin;
  y = y*delta_y + domainsize.ymin;

  coords[0] = x;
  coords[1] = y;

}  // function redimensionalize

// a function f that performs mapping x* - f(x*), where x* is the nondimensional
// coordinates (x and y)
void mapFunction(DomainSize domainsize, apf::Vector3& coords)
{
/*
  // calculate the element size in the x and y directions in non-dimensional
  // coordinates
  double hx = 1/( (double)domainsize.numElx);
  double hy = 1/( (double)domainsize.numEly);
  double nwaves = 5;  // number of waves over the domain (periods of sin)
  double amplitude = 0.2; // fraction of hx used as perturbation magnitude

  // apply a small amplitude sin wave
  coords[0] = coords[0] + amplitude*hx*sin( M_PI*coords[0]*2*nwaves);
  coords[1] = coords[1] + amplitude*hy*sin( M_PI*coords[1]*2*nwaves);
*/  
  /*
  coords[0] = sin( (M_PI*0.5)*coords[0] );
  coords[1] = sin( (M_PI*0.5)*coords[1] );
  */
}  // function mapFunction
