#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <iostream>
#include <stdbool.h>
#include <limits.h>
//#include "funcs1.h"
#include "apfSBPShape.h"
#include "apfSBPShape3.h"

// 3D array
struct _Array3 {
  apf::MeshEntity** ptr;
  int m;
  int n;
  int p;
};
typedef struct _Array3 Array3;
 


struct _Sizes {
  int numElx;
  int numEly;
  int numElz;
};
typedef struct _Sizes Sizes;

struct _Geom {
  int model_tag;
  int model_dim;
};
typedef struct _Geom Geom;

struct _VertIdx {
  int i;
  int j;
  int k;
};

typedef struct _VertIdx VertIdx;

// addition on VertIdx
VertIdx add(const VertIdx v1, const VertIdx v2);
VertIdx add(const VertIdx v1, const int arr[3]);

#define NTETS 6
#define NEDGES 19
#define NFACES 18
int tets[NTETS][4][3];  // define vertices used for constructing tetrahedrons

int edges[NEDGES][2][3];  // describe vertices defining each edge of the tets
int faces[NFACES][3][3];  // describe the vertices defining each face of hte tets

#define NREAREDGES 12
int rear_edge_idx[NREAREDGES]; // indices of the rear edges in the edges list
bool create_edge[NREAREDGES];  // whether or not to create the edge 

int rear_face_idx[6];  // indices of rear faces in the faces list
bool create_face[6];  // whether or not to create the faces

// 3D array functions
//-----------------------------------------------------------------------------
int getInd3(const int i, const int j, const int k, const int imax, const int jmax, const int kmax)
{
  return i + j*imax + imax*jmax*k;
}

apf::MeshEntity* getindex3(Array3 arr, const int i, const int j, const int k)
{
 const int idx = getInd3(i, j, k, arr.m, arr.n, arr.p);
 return arr.ptr[idx];
}

apf::MeshEntity* getindex3(Array3 arr, VertIdx v)
{
  return getindex3(arr, v.i, v.j, v.k);
}

void setindex3(Array3 arr, const int i, const int j, const int k, 
              apf::MeshEntity* val)
{
 const int idx = getInd3(i, j, k, arr.m, arr.n, arr.p);
// printf("setting arr %i %i %i = linear index %i to value %10.6f\n", i, j, k, idx, val);
 arr.ptr[idx] = val;
}

void setindex3(Array3 arr, VertIdx v,
               apf::MeshEntity* val)
{
  setindex3(arr, v.i, v.j, v.k, val);
}

void zero3(Array3 arr)
{
  int len = arr.m*arr.n*arr.p;
  for (int i = 0; i < len; ++i)
  {
    arr.ptr[i] = NULL;
  }
}

//-----------------------------------------------------------------------------


// declare the verts used to break the cube into tetrahedra
void declareTets();

// get the classification of vertices at indices (i, j, k) in the array of vertices
Geom getVertClassification(int i, int j, int k, Sizes sizes);
Geom getVertClassification(VertIdx v, Sizes sizes);
// populate the global variable edges by inspecting tets
void extractEdges();

// determine if the edge defined by the two vertices is already in the edges list
bool checkEdges(int vert1[3], int vert2[3], int nedges);

// determine if the face defined by the 3 vertices is already in the faces list
bool checkFace(int vert1[3], int vert2[3], int vert3[3], int nfaces);

// check of an entry in either edges or faces contains the vertex
bool contains(int list_entry[][3], const int dim1, int vert1[3]);

// populate the faces global variable by inspecting the tets variable
void extractFaces();

// populate the rear_edge_idx and create_edge arrays
void extractRearEdges();

// identiry which edges need to be created for a given element
void identifyEdges(VertIdx start);

// populate rear_face_idx and create_face arrays
void extractRearFaces();

// populate create_face
void identifyFaces(VertIdx start);

// create the edges on a brick
void createEdges(apf::Mesh2* m, Sizes sizes, VertIdx start, apf::MeshEntity**** verts);

// create a single edge
void createEdge(apf::Mesh2* m, apf::MeshEntity* edge_verts[2], apf::ModelEntity* model_entity);

// create all the faces on a cube
void createFaces(apf::Mesh2* m, Sizes sizes, VertIdx start, apf::MeshEntity**** verts);

// create a single face
void createFace(apf::Mesh2* m, apf::MeshEntity* edge_verts[3], apf::ModelEntity* model_entity);

// get the model entity given the geometric classification of the vertices
apf::ModelEntity* getModelEntity(apf::Mesh* m, Geom g_class[], const int ng);

// get a vert from the array of all vertices
apf::MeshEntity* getVert(VertIdx pos, apf::MeshEntity**** verts);

// get the vertices of a tet
void getTetVerts(VertIdx start, int tet[4][3], apf::MeshEntity**** verts_all,  apf::MeshEntity* verts[4]);

// copy array of length 3
void copyArray(int src[3], int dest[3]);

// print array to stdout
void printArray(int vals[3]);

// print edges array
void printEdges(int nedges);

// translate binary offsets into integer
int getVertNum(int vert[3]);

// find index of entry in array
int contains(int vals[], const int len, int  val);

// print the verts array
void printVerts(apf::MeshEntity**** verts, Sizes sizes);

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
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 3, false);
/*
  apf::FieldShape* m_shape = m->getShape();
  const char shape_name[] = m_shape->getName();
  std::cout << "shape name = " << shape_name << std::endl;
*/
  if (argc < 3 || argc > 5)
  {
    std::cerr << "Error: wrong number of arguments" << std::endl;
    return 1;
  }
  // Problem 3
  std::cout << "Problem 3 output: " << std::endl;
  int numElx = 0;  // number of elements in x direction
  sscanf(argv[1], "%d", &numElx);
  int numEly = 0;  // number of elements in y direction
  sscanf(argv[2], "%d", &numEly);
  int numElz = 0; // number of elements in the z direction
  sscanf(argv[3], "%d", &numElz);
  int apply_pert = 0;
  if (argc == 5)
  {
    sscanf(argv[4], "%d", &apply_pert);
  }
  Sizes sizes = {numElx, numEly, numElz};
  declareTets();
  extractEdges();
  extractFaces();
  extractRearEdges();
  extractRearFaces();

  // TODO: fix this 
  long numElx_l = numElx;
  long numEly_l = numEly;
  long nvert = (numElx_l+1)*(numEly_l+1);
  long nedges = numElx_l*(numEly_l+1) + numEly_l*(numElx_l+1) + numElx_l*numEly_l;
  long nel = numElx_l*numEly_l + 1;

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

  double xmin = 1.0;
  double ymin = 1.0;
  double zmin = 1.0;
  double xdist = 2;  // xmax - xmin
  double ydist = 2;  // ymax - ymin
  double zdist = 2;
  double x_spacing = xdist/numElx;  // spacing of el
  double y_spacing = ydist/numEly;
  double z_spacing = zdist/numElz;
  double x_0 = xmin;  // x coordinate of lower left corner of current element
  double y_0 = ymin ; // y coordinate of lower left corner of current element
  double z_0 = zmin;
  /* nominal locations of current point */
  double x_i = x_0;
  double y_i = y_0;
  double z_i = z_0;
  /* possibly perturbed locations of current point */
  double x_inner = x_i;
  double y_inner = y_i;
  double z_inner = z_i;
  double pert_fac = 10*M_PI;
  double pert_mag = 0.1;

  std::cout << "SIZE_MAX = " << SIZE_MAX << std::endl;
  std::cout << "about to allocate memory" << std::endl;
  apf::MeshEntity**** vertices = (apf::MeshEntity****) calloc(numElx+1, sizeof(apf::MeshEntity***));
  for (int i = 0; i < (numElx+1); ++i)
  {
    vertices[i] = (apf::MeshEntity***) calloc(numEly+1, sizeof(apf::MeshEntity**));
    for (int j = 0; j < (numEly+1); ++j)
    {
      vertices[i][j] = (apf::MeshEntity**) calloc(numElz+1, sizeof(apf::MeshEntity*));
    }
    
  }
//  apf::MeshEntity* vertices[numElx+1][numEly+1];  // hold pointers to all vertices
  std::cout << "finished allocating memory" << std::endl;
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point
  Geom geo;
  
  apf::ModelEntity* model_entity;  // the model entity itself


  std::cout << "Creating " << numElx << " by " << numEly << "by " << numElz; 
  std::cout << " mesh" << std::endl;

  // create all vertices
  for (int k = 0; k < (numElz+1); ++k)
  {
    for (int j = 0; j < (numEly+1); ++j)  // loop up columns (bottom to top)
    {
      for (int i = 0; i < (numElx+1); ++i) // loop over rows (left to right)
      {
         std::cout << "creating point at x = " << x_inner << " , y = " << y_inner << " , z = " << z_inner << std::endl;
         coords_i[0] = x_inner;
         coords_i[1] = y_inner;
         coords_i[2] = z_inner;

         geo = getVertClassification(i, j, k, sizes);
         std::cout << "model_dim = " << geo.model_dim << std::endl;
         std::cout << "model_tag = " << geo.model_tag << std::endl;
         model_entity = m->findModelEntity(geo.model_dim, geo.model_tag);
         std::cout << "model entity = " << model_entity << std::endl;
           
         vertices[i][j][k] = m->createVert(model_entity);
         std::cout << "vertex = " << vertices[i][j][k] << std::endl;
         m->setPoint(vertices[i][j][k], 0, coords_i);
         x_inner = x_inner + x_spacing + pert_mag*sin(apply_pert*pert_fac*x_inner);
         y_inner = y_i + pert_mag*sin(apply_pert*pert_fac*x_inner);
         z_inner = z_i + pert_mag*sin(apply_pert*pert_fac*y_inner);
      }

      y_i = y_i + y_spacing;  // increment y_i
      x_i = x_0;  // reset x_i to beginning of row

      y_inner = y_i;// + pert_mag*sin(apply_pert*pert_fac*y_i);
      x_inner = x_i; // + pert_mag*sin(apply_pert*pert_fac*x_i);
      z_inner = z_i;
    }
    z_i = z_i + z_spacing;
    x_i = x_0;  // reset to beginning of row
    y_i = y_0; // reset to beginning of row

    y_inner = y_i;
    x_inner = x_i;
    z_inner = z_i;

  }

  printVerts(vertices, sizes);

  // classify all MeshEntities on geometric face
  int model_dim = 3;
  int model_tag = 0;
  model_entity = m->findModelEntity(model_dim, model_tag);
  VertIdx start;
  apf::MeshEntity* vertices_i[4];  // hold vertices for a particular element
  // build element from verticies
  for (int k = 0; k < numElz; ++k)
  {
    for (int j = 0; j < numEly; ++j)
    {
      for (int i = 0; i < numElx; ++i)
      {
        std::cout << "\ncreating element " << i << ", " << j << ", " << k << std::endl;
  //      int el_num = i + numElx*j;
  //      std::cout << "creating element i = " << i << " , j = " << j << std::endl;
        start.i = i; start.j = j; start.k = k;
        identifyEdges(start);
        createEdges(m, sizes, start, vertices);
        identifyFaces(start);
        createFaces(m, sizes, start, vertices);

        for (int v = 0; v < 6; ++v)
        {
          getTetVerts(start, tets[v], vertices, vertices_i);
          apf::buildElement(m, model_entity, apf::Mesh::TET, vertices_i);
        }


      }
    }
  }

  // build, verify  mesh
/*
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m);
  std::cout << "finished deriving model" << std::endl;
*/
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
//  checkMesh(m);
  m->verify();
  std::cout << "verified" << std::endl;

 // for quadratic meshes
//  apf::FieldShape* linear2 = apf::getSerendipity();
    apf::FieldShape* linear2 = apf::getSBP3Shape(1);
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
  m->writeNative("./meshfiles/abc.smb");
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


void declareTets()
{
  int vals[] = {0, 0, 0}; copyArray(vals, tets[0][0]);
  int vals2[] = {1, 0, 0}; copyArray(vals2, tets[0][1]);
  int vals3[] = {1, 1, 1}; copyArray(vals3, tets[0][2]);
  int vals4[] = {1, 0, 1}; copyArray(vals4, tets[0][3]);

  int vals5[] = {1, 1, 1}; copyArray(vals5, tets[1][0]);
  int vals6[] = {0, 0, 0}; copyArray(vals6, tets[1][1]);
  int vals7[] = {1, 0, 1}; copyArray(vals7, tets[1][2]);
  int vals8[] = {0, 0, 1}; copyArray(vals8, tets[1][3]);

  int vals9[] = {0, 0, 0}; copyArray(vals9, tets[2][0]);
  int vals10[] = {1, 0, 0}; copyArray(vals10, tets[2][1]);
  int vals11[] = {1, 1, 0}; copyArray(vals11, tets[2][2]);
  int vals12[] = {1, 1, 1}; copyArray(vals12, tets[2][3]);

  int vals13[] = {0, 1, 0}; copyArray(vals13, tets[3][0]);
  int vals14[] = {0, 0, 0}; copyArray(vals14, tets[3][1]);
  int vals15[] = {1, 1, 0}; copyArray(vals15, tets[3][2]);
  int vals16[] = {1, 1, 1}; copyArray(vals16, tets[3][3]);

  int vals17[] = {0, 1, 0}; copyArray(vals17, tets[4][0]);
  int vals18[] = {0, 0, 0}; copyArray(vals18, tets[4][1]);
  int vals19[] = {1, 1, 1}; copyArray(vals19, tets[4][2]);
  int vals20[] = {0, 1, 1}; copyArray(vals20, tets[4][3]);

  int vals21[] = {0, 0, 0}; copyArray(vals21, tets[5][0]);
  int vals22[] = {1, 1, 1}; copyArray(vals22, tets[5][1]);
  int vals23[] = {0, 1, 1}; copyArray(vals23, tets[5][2]);
  int vals24[] = {0, 0, 1}; copyArray(vals24, tets[5][3]);

}

Geom getVertClassification(VertIdx v, Sizes sizes)
{
  return getVertClassification(v.i, v.j, v.k, sizes);
}

Geom getVertClassification(int i, int j, int k, Sizes sizes)
{
 int model_dim;
 int model_tag;
 int numElx = sizes.numElx;
 int numEly = sizes.numEly;
 int numElz = sizes.numElz;
 // figure out what model entity to classify this vertex on 
/* corners in x-y plane */
if ( i == 0 && j == 0 && k == 0 )
{
  model_dim = 0;
  model_tag = 0;
} else if ( i == (numElx) && j == 0 && k == 0 )
{
  model_dim = 0;
  model_tag = 1;
} else if ( i == (numElx) && j == (numEly) && k == 0)
{
  model_dim = 0;
  model_tag = 2;
} else if ( i == 0 && j == (numEly) && k == 0 )
{
  model_dim = 0;
  model_tag = 3;
/* corners in z = numElz plane */
} else if ( i == 0 && j == 0 && k == numElz ) // corners in z=numElz plane
{
  model_dim = 0;
  model_tag = 4;
} else if ( i == (numElx) && j == 0 && k == numElz )
{
  model_dim = 0;
  model_tag = 5;
} else if ( i == (numElx) && j == (numEly) && k == numElz)
{
  model_dim = 0;
  model_tag = 6;
} else if ( i == 0 && j == (numEly) && k == numElz )
{
  model_dim = 0;
  model_tag = 7;
/*  domain edges in x-y plane */
} else if ( j == 0 && k == 0)   // bottom row
{
  model_dim = 1;
  model_tag = 0;
} else if ( i == (numElx) && k == 0)  // right column
{
  model_dim = 1;
  model_tag = 1;
} else if ( j == (numEly) && k == 0)  // top row
{
  model_dim = 1;
  model_tag = 2;
} else if ( i == 0 && k == 0 )  // left column
{
  model_dim = 1;
  model_tag = 3;
/*  domain edges in z = numElz plane */
/*  counterclockwise from the origin */
} else if ( j == 0 && k == numElz)
{
  model_dim = 1;
  model_tag = 4;
} else if ( i == (numElx) && k == numElz)
{
  model_dim = 1;
  model_tag = 5;
} else if ( j == (numEly) && k == numElz)
{
  model_dim = 1;
  model_tag = 6;
} else if ( i == 0 && k == numElz )
{
  model_dim = 1;
  model_tag = 7;
/*  domain edges in parallel to z axis */
/*  counterclockwise starting at the origin */
} else if ( i == 0 && j == 0)
{
  model_dim = 1;
  model_tag = 8;
} else if ( i == numElx && j == 0)
{
  model_dim = 1;
  model_tag = 9;
} else if ( i == numElx && j == numEly)
{
  model_dim = 1;
  model_tag = 10;
} else if ( i == 0 && j == numEly )
{
  model_dim = 1;
  model_tag = 11;
/* domain faces */
} else if ( k == 0 )  // bottom face (x-y plane)
{
  model_dim = 2;
  model_tag = 0;
} else if ( j == 0)  // x-z plane
{
  model_dim = 2;
  model_tag = 1;
} else if ( i == numElx) // parallel to y-z plane
{
  model_dim = 2;
  model_tag = 2;
} else if ( j == numEly)  // parallel to x-z plane
{
  model_dim = 2;
  model_tag = 3;
} else if ( i == 0)  // y-z plane
{
  model_dim = 2;
  model_tag = 4;
} else if ( k == numElz)  // top face (parallel to x-y plane)
{
  model_dim = 2;
  model_tag = 5;

} else  // classify on face
{
  model_dim = 3; 
  model_tag = 0;
}

  Geom geo = {model_tag, model_dim};
  return geo;
}

// determine the unique set of edges as defined by their vertices
void extractEdges()
{
  std::cout << "extracting edges" << std::endl;
// algorithm: loop over vertices in tet, form edge between current vertex and
//            all previous vertices, check if edge is already in list
//            add it if not
  for (int i=0; i < NFACES; ++i)
  {
    for (int j=0; j < 2; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        edges[i][j][k] = 0;
      }
    }
  }


  int nedges=0;
  for (int i=0; i < NTETS; ++i)
  {
    for (int j=0; j < 4; ++j)
    {
      for (int k=(j+1); k < 4; ++ k)
      {
        // check against list of existing edges
        bool already_added = checkEdges(tets[i][j], tets[i][k], nedges);

        if (!already_added)
        {
          if (nedges >= NEDGES)
            std::cerr << "  Warning: too many edges detected" << std::endl;

          copyArray(tets[i][j], edges[nedges][0]);
          copyArray(tets[i][k], edges[nedges][1]); 
          ++nedges;
        }
      }
    }
  }

  if (nedges != NEDGES)
  {
    std::cerr << "Warning: wrong number of edges detected" << std::endl;
    throw nedges;
  } else
  {
    std::cout << "Correct number of edges detected" << std::endl;
  }
}

// check if the edges list already contains the edge defined by 
// these two verts
bool checkEdges(int vert1[3], int vert2[3], int nedges)
{
  bool found1;
  bool found2;
  for (int i=0; i < nedges; ++i)
  {
    // check if edges[i] contains vert1
    found1 = contains(edges[i], 2, vert1);
    found2 = contains(edges[i], 2, vert2);

    if (found1 && found2)
      return true;
  }

  return false;

}

// get the index in the edge array containing the specified edge
int findEdge(int edge[2][3])
{
  int idx = -1;
  bool found1;
  bool found2;
  for (int i=0; i < NEDGES; ++i)
  {
    found1 = contains(edges[i], 2, edge[0]);
    found2 = contains(edges[i], 2, edge[1]);
    if (found1 && found2)
    {
      idx = i;
      break;
    }
  }

  return idx;
}


bool checkFace(int vert1[3], int vert2[3], int vert3[3], int nfaces)
{
  bool found1;
  bool found2;
  bool found3;

  for (int i=0; i < NFACES; ++i)
  {
    found1 = contains(faces[i], 3, vert1);
    found2 = contains(faces[i], 3, vert2);
    found3 = contains(faces[i], 3, vert3);

    if (found1 && found2 && found3)
      return true;
  }
  
  return false;
}

int findFace(int face[3][3])
{
  int idx = -1;
  bool found1;
  bool found2;
  bool found3;

  for (int i=0; i < NFACES; ++i)
  {
    found1 = contains(faces[i], 3, face[0]);
    found2 = contains(faces[i], 3, face[1]);
    found3 = contains(faces[i], 3, face[2]);

    if (found1 && found2 && found3)
    {
      idx = i;
      break;
    }
  }

  return idx;
}



// check if the specified entry in th list contains a given
// vertex (in either position)
bool contains(int list_entry[][3], const int dim1, int vert1[3])
{
  bool matched;
  for (int i = 0; i < dim1; ++i)
  {

    matched = true;
    for (int j=0; j < 3; ++j)
    {
      matched = matched && (vert1[j] == list_entry[i][j]);
    }
    if (matched)
      return matched;
  }

  return matched;
}

// get the list of triangles 
void extractFaces()
{
  std::cout << "extracting faces" << std::endl;
// algorithm: loop over vertices in tet, form face between current vertex and
//            all sets of 2 future vertices, check if face is already in list
//            add it if not

  for (int i=0; i < NFACES; ++i)
  {
    for (int j=0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k)
      {
        faces[i][j][k] = 0;
      }
    }
  }


  bool foundface;
  int nfaces = 0;
  for (int i=0; i < NTETS; ++i)
  {
    for (int v1=0; v1 < 4; ++v1)
    {
      for (int v2=(v1+1); v2 < 4; ++v2)
      {
        for (int v3=(v2+1); v3 < 4; ++v3)
        {
          foundface = checkFace(tets[i][v1], tets[i][v2], tets[i][v3], nfaces);
          if (!foundface)
          {
            if (nfaces >= NFACES)
              std::cerr << "Warning: too many faces detected" << std::endl;

            copyArray(tets[i][v1], faces[nfaces][0]); 
            copyArray(tets[i][v2], faces[nfaces][1]); 
            copyArray(tets[i][v3], faces[nfaces][2]); 
            ++nfaces;
          }
        }
      }
    }
  }

  if ( nfaces != NFACES)
  {
    std::cerr << "Warning: wrong number of faces detected" << std::endl;
    throw nfaces;
  } else
  { 
    std::cout << "Correct number of faces detected" << std::endl;
  }
}


void extractRearEdges()
{
  int rear_edges[NREAREDGES][2][3] ={ 
    // edge on axes of the cube
    { {0, 0, 0}, {1, 0, 0} },
    { {0, 0, 0}, {0, 1, 0} },
    { {0, 0, 0}, {0, 0, 1} },
  
    // j = 0 edges
    { {0, 0, 0}, {1, 0, 1} },
    { {1, 0, 0}, {1, 0, 1} },
    { {0, 0, 1}, {1, 0, 1} },
    // i = 0 edges
    { {0, 0, 0}, {0, 1, 1} },
    { {0, 1, 0}, {0, 1, 1} },
    { {0, 0, 1}, {0, 1, 1} },
    // k = 0 edges
    { {0, 0, 0}, {1, 1, 0} },
    { {1, 0, 0}, {1, 1, 0} },
    { {0, 1, 0}, {1, 1, 0} },
  };

  for (int i = 0; i < NREAREDGES; ++i)
  {
    rear_edge_idx[i] = findEdge(rear_edges[i]);
    create_edge[i] = false;
    if (rear_edge_idx[i] < 0)
    {
      std::cerr << "rear edge " << i << " not found" << std::endl;
      throw i;
    }
    std::cout << "rear edge " << i << " idx = " << rear_edge_idx[i] << std::endl;

  }
}

//  set the values in create_edge
//  start should be the indices of the vertex at the origin of the cube
void identifyEdges(VertIdx start)
{
  for (int i=0; i < NREAREDGES; ++i)
  {
    create_edge[i] = false;
  }

  const int i = start.i;
  const int j = start.j;
  const int k = start.k;
  // select which edges to create
  if (j == 0 && k == 0)
    create_edge[0] = true;
  if (i == 0 && k == 0)
    create_edge[1] = true;
  if ( i == 0 && j == 0)
    create_edge[2] = true;
  if ( j == 0 )
  {
    create_edge[3] = true;
    create_edge[4] = true;
    create_edge[5] = true;
  }
  if ( i == 0)
  {
    create_edge[6] = true;
    create_edge[7] = true;
    create_edge[8] = true;
  }
  if ( k == 0)
  {
    create_edge[9] = true;
    create_edge[10] = true;
    create_edge[11] = true;
  }

}

    

void extractRearFaces()
{
  int rear_faces[6][3][3] = {
    // face 2
    { {0, 0, 0}, {1, 0, 1}, {0, 0, 1} },
    { {0, 0, 0}, {1, 0, 0}, {1, 0, 1} },
     // face 5
    { {0, 0, 0}, {0, 1, 0}, {0, 1, 1} },
    { {0, 0, 0}, {0, 1, 1}, {0, 0, 1} }, 
    // face 1
    { {0, 0, 0}, {1, 1, 0}, {1, 0, 0} },
    { {0, 0, 0}, {1, 1, 0}, {0, 1, 0} },
  };

  for (int i = 0; i < 6; ++i)
  {
    rear_face_idx[i] = findFace(rear_faces[i]);
    create_face[i] = false;
    if (rear_face_idx[i] < 0)
    {
      std::cerr << "rear face " << i << " not found" << std::endl;
      throw i;
    }

    std::cout << "rear face " << i << " idx = " << rear_face_idx[i] << std::endl;
  }
}

void identifyFaces(VertIdx start)
{
  for (int i=0; i < 6; ++i)
  {
    create_face[i] = false;
  }

  // select which faces to create
  if (start.j == 0)
  {
    create_face[0] = true;
    create_face[1] = true;
  }
  if (start.i == 0)
  {
    create_face[2] = true;
    create_face[3] = true;
  }
  if (start.k == 0)
  {
    create_face[4] = true;
    create_face[5] = true;
  }
}


// create the needed edges on an element
// TODO: need Mesh*, need sizes
void createEdges(apf::Mesh2* m, Sizes sizes, VertIdx start, apf::MeshEntity**** verts)
{
  std::cout << "\ncreating edges" << std::endl;
  std::cout << "start idx = " << start.i << ", " << start.j << ", " << start.k << std::endl;
  int idx;
  VertIdx vertidx;
  apf::MeshEntity* verts_i[2];
  Geom g_class[2]; // geometric classification of the vertices
  apf::ModelEntity* model_entity;
  for (int i = 0; i < NEDGES; ++i)
  {
    std::cout << "\nconsidering edge " << i << std::endl;
    idx = contains(rear_edge_idx, NREAREDGES, i);

    if (idx > 0)  // if found
    {
      std::cout << "this is a rear edge" << std::endl;
      if (!create_edge[idx])
      {
        std::cout << "skipping this edge" << std::endl;
        continue;
      }
    }

    std::cout << "creating edge" << std::endl;
    // if we get here, we should create the edge
    for (int j=0; j < 2; ++j)
    {
      vertidx = add(start, edges[i][j]);

      std::cout << "vert idx = " << vertidx.i << ", " << vertidx.j << ", " << vertidx.k << std::endl;
      verts_i[j] = getVert(vertidx, verts);
      std::cout << "vert " << j << " = " << verts_i[j] << std::endl;
      g_class[j] = getVertClassification(vertidx, sizes);
    }

    model_entity = getModelEntity(m, g_class, 2);
    createEdge(m, verts_i, model_entity);
  }  

}


void createEdge(apf::Mesh2* m, apf::MeshEntity* edge_verts[2], apf::ModelEntity* model_entity) 
{
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

}

void createFaces(apf::Mesh2* m, Sizes sizes, VertIdx start, apf::MeshEntity**** verts)
{
  std::cout << "\ncreating faces" << std::endl;
  int idx;
  VertIdx vertidx;
  apf::MeshEntity* verts_i[3];
  Geom g_class[3]; // geometric classification of the vertices
  apf::ModelEntity* model_entity;
  for (int i = 0; i < NFACES; ++i)
  {
    std::cout << "\nconsidering face " << i << std::endl;
    idx = contains(rear_face_idx, 6, i);

    if (idx > 0)  // if found
    {
      std::cout << "face is a rear face" << std::endl;
      if (!create_face[idx])
      {
        std::cout << "not creating face" << std::endl;
        continue;
      }
    }
    std::cout << "creating face " << std::endl;
    // if we get here, we should create the edge
    for (int j=0; j < 3; ++j)
    {

      vertidx = add(start, faces[i][j]);
      std::cout << "vert idx = " << vertidx.i << ", " << vertidx.j << ", " << vertidx.k << std::endl;
      verts_i[j] = getVert(vertidx, verts);
      g_class[j] = getVertClassification(vertidx, sizes);
      std::cout << "vert " << j << " = " << verts_i[j] << std::endl;
    }
    std::cout << "finished getting geometric classification" << std::endl;

    model_entity = getModelEntity(m, g_class, 3);
    std::cout << "model_entity = " << model_entity << std::endl;
    std::cout << "creating face" << std::endl;
    createFace(m, verts_i, model_entity);
    std::cout << "finished creating face" << std::endl;
  }  

}


void createFace(apf::Mesh2* m, apf::MeshEntity* face_verts[3], apf::ModelEntity* model_entity) 
{
  std::cout << "face_verts = " << face_verts[0] << ", " << face_verts[1];
  std::cout << ", " << face_verts[2] << std::endl;
  std::cout << "model entity = " << model_entity << std::endl;
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);
}




// get the ModelEntity that the edge defind by g1 and g2 should be classified
// on
apf::ModelEntity* getModelEntity(apf::Mesh* m, Geom g_class[], const int ng)
{
  int model_dim;
  int model_tag;

  int max_dim = 0;
  int min_tag = INT_MAX;
  int n_max_dim = 0;

  // find the maximum dimension
  for (int i=0; i < ng; ++i)
  {
    model_dim = g_class[i].model_dim;
    if (model_dim > max_dim)
    {
      max_dim = model_dim;
      n_max_dim = 1;
    } else if (model_dim == max_dim)
    {
      ++n_max_dim;
    }
  }

  // now find the minimum tag on that dimension
  for (int i=0; i < ng; ++i)
  {
    if (g_class[i].model_dim == max_dim)
    {
      model_tag = g_class[i].model_tag;
      if (model_tag < min_tag)
        min_tag = model_tag;
    }
  }

  std::cout << "model_dim = " << max_dim << std::endl;
  std::cout << "model_tag = " << min_tag << std::endl;
  // max_dim and min_tag now fully describe the geometric entity
  apf::ModelEntity* me = m->findModelEntity(max_dim, min_tag);
  std::cout << "model entity = " << me << std::endl;
  return me;
}


// get a vert from the array of all vertices
apf::MeshEntity* getVert(VertIdx pos, apf::MeshEntity**** verts)
{
  apf::MeshEntity* vert = verts[pos.i][pos.j][pos.k];

  return vert;
}

void getTetVerts(VertIdx start, int tet[4][3], apf::MeshEntity**** verts_all,  apf::MeshEntity* verts[4])
{
  VertIdx pos;
  for (int i=0; i < 4; ++i)
  {
    pos = add(start, tet[i]);
    verts[i] = getVert(pos, verts_all);
  }
}


void copyArray(int src[3], int dest[3])
{
  for (int i = 0; i < 3; ++ i)
    dest[i] = src[i];
}

void printArray(int vals[3])
{
  std::cout << vals[0] << ", " << vals[1] << ", " << vals[2];
}

void printEdges(int nedges)
{
  std::cout << "edges = " << std::endl;
  for (int i=0; i < nedges; ++i)
  {
    std::cout << "  edge " << i << ": v1 = " << getVertNum(edges[i][0]);
    std::cout << ", v2 = " << getVertNum(edges[i][1]); 
    std::cout << std::endl;
  }
}

int getVertNum(int vert[3])
{
  return vert[0] + 2*vert[1] + 4*vert[2];
}

int contains(int vals[], const int len, int  val)
{
  int idx = -1;
  for (int i = 0; i < len; ++i)
  {
    if (vals[i] == val)
    {
      idx = i;
      break;
    }
  }

  return idx;
}

VertIdx add(const VertIdx v1, const VertIdx v2)
{
  VertIdx result = {v1.i + v2.i, v1.j + v2.j, v1.k + v2.k};
  return result;
}

VertIdx add(const VertIdx v1, const int arr[3])
{
  VertIdx result = {v1.i + arr[0], v1.j + arr[1], + v1.k + arr[2]};
  return result;
}

void printVerts(apf::MeshEntity**** verts, Sizes sizes)
{
  std::cout << "verts: " << std::endl;
  for (int i=0; i < (sizes.numElx+1); ++i)
    for (int j=0; j < (sizes.numEly+1); ++j)
      for(int k=0; k < (sizes.numElz+1); ++k)
      {
        std::cout << "vert[" << i << "][" << j << "][" << k << "] = " << verts[i][j][k] << std::endl;
      }
}
