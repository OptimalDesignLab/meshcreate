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
#include <sys/time.h>
#include <sys/resource.h>

//#include "funcs1.h"
//#include "apfSBPShape.h"
//#include "apfSBPShape3.h"

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

bool isequal(Geom g1, Geom g2);

struct _VertIdx {
  int i;
  int j;
  int k;
};

typedef struct _VertIdx VertIdx;

struct _Counts {
  long numVert;
  long numEdge;
  long numFace;
  long numEl;
};

typedef struct _Counts Counts;

// addition on VertIdx
VertIdx add(const VertIdx v1, const VertIdx v2);
VertIdx add(const VertIdx v1, const int arr[3]);
void add(const VertIdx v1, const int arr[3], VertIdx &vout);
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

const int face_edges[6][4] = {
                         {0, 1, 2, 3},
                         {0, 9, 4, 8},
                         {1, 10, 5, 9},
                         {2, 11, 6, 10},
                         {3, 8, 7, 11},
                         {4, 5, 6, 7},
                       };
const int face_verts[6][4] = {
                           {0, 1, 2, 3},
                           {0, 1, 4, 5},
                           {1, 2, 5, 6},
                           {2, 3, 6, 7},
                           {0, 3, 4, 7},
                           {4, 5, 6, 7},
                         };
// map [vert1][vert2] -> {model_dim, model_tag}
const int vert_only_classification[8][8][2] = 
  { 
    // vert1 to all others
    { {-1, -1}, {1, 0}, {2, 0}, {1, 3}, {1, 8}, {2, 1}, {3, 0}, {2, 4} },
    // vert2 to all others
    { {1, 0}, {-1, -1}, {1, 1}, {2, 0}, {2, 1}, {1, 9}, {2, 2}, {3, 0} },
    // vert3 to all others
    { {2, 0}, {1, 1}, {-1, -1}, {1, 2}, {3, 0}, {2, 2}, {1, 10}, {2, 3} },
    // vert4 to all others
    { {1, 3}, {2, 0}, {1, 2}, {-1, -1}, {2, 4}, {3, 0}, {2, 3}, {1, 11} },
    // vert5 to all others
    { {1, 8}, {2, 1}, {3, 0}, {2, 4}, {-1, -1}, {1, 4}, {2, 5}, {1, 7} },
    // vert6 to all others
    { {2, 1}, {1, 9}, {2, 2}, {3, 0}, {1, 4}, {-1, -1}, {1, 5}, {2, 5} },
    // vert7 to all others
    { {3, 0}, {2, 2}, {1, 10}, {2, 3}, {2, 5}, {1, 5}, {-1, -1}, {1,6} },
    // vert8 to all others
    { {2, 4}, {3, 0}, {2, 3}, {1, 11}, {1, 7}, {2, 5}, {1, 6}, {-1, -1} }
  };

// map geometric verts to geometric edges
const int vert_edges[8][3] = 
  {
    {0, 3, 8},
    {0, 1, 9},
    {1, 2, 10},
    {2, 3, 11},
    {4, 7, 8},
    {4, 5, 9},
    {5, 6, 10},
    {6, 7, 11},
  };
// declare the verts used to break the cube into tetrahedra
void declareTets();

// get the classification of vertices at indices (i, j, k) in the array of vertices
Geom getVertClassification(int i, int j, int k, Sizes sizes);
Geom getVertClassification(VertIdx v, Sizes sizes);

// get classification of a new edge
Geom getEdgeClassification(apf::Mesh* m, VertIdx _v1, VertIdx _v2, 
                           apf::MeshEntity* verts[2], Sizes sizes);


// get the classification of an edge defined by 2 verts
// the edge must already exist
Geom getEdgeClassification(apf::Mesh* m, apf::MeshEntity* verts[2]);
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

Geom getMaxGeometry(apf::Mesh* m, Geom g_class[], const int ng);

// get the model entity
apf::ModelEntity* getModelEntity(apf::Mesh* m, Geom g_class);

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
int contains(const int vals[], const int len, int  val);

// print the verts array
void printVerts(apf::MeshEntity**** verts, Sizes sizes);

// get time with microsecond precision
double get_time();

void calcCentroid(apf::Mesh* m, apf::MeshEntity** verts,  double centroid[3]);

void checkMesh(apf::Mesh* m, Sizes sizes, Counts counts);

void checkMeshFaces(apf::Mesh* m);

class FaceCallback : public apf::BuildCallback
{
  public:
    apf::Mesh2* m_closure;
    void call(apf::MeshEntity* e)
    {
      int type = m_closure->getType(e);
      const char* tname =  apf::Mesh::typeName[type];
      int dim = apf::Mesh::typeDimension[type];

      std::cout << "new entity of type " << tname << " (dimension " << dim << ") created" << std::endl;
    }
};
class RegionCallback : public apf::BuildCallback
{
  public:
    apf::Mesh2* m_closure;
    void call(apf::MeshEntity* e)
    {
      int type = m_closure->getType(e);
      const char* tname =  apf::Mesh::typeName[type];
      int dim = apf::Mesh::typeDimension[type];
      if (dim == 2)
      {
        apf::Downward down;
        m_closure -> getDownward(e, 0, down);
        double centroid[3];
        calcCentroid(m_closure, down, centroid);
        std::cout << "created unexpected triangles with centroid = " << centroid[0];
        std::cout << ", " << centroid[1] << ", " << centroid[2] << std::endl;
        std::cout << "face vertices = " << down[0] << ", " << down[1] << ", ";
        std::cout << down[2] << std::endl;
      } else
      {


        std::cout << "new entity of type " << tname << " (dimension " << dim << ") created" << std::endl;
      }
    }
};


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
  long numElz_l = numElz;
  long nvert = (numElx_l+1)*(numEly_l+1)*(numElz_l+1);
  long nedges =  // cube edges
                 numElx_l*(numEly_l+1)*(numElz_l+1)
                + numEly_l*(numElx_l+1)*(numElz_l+1)
                + numElz_l*(numElx_l+1)*(numEly_l+1) 
                // cube face edges
                + numElx_l*numElz_l*(numEly_l+1)
                + numEly_l*numElz_l*(numElx_l+1)
                + numElx_l*numEly_l*(numElz_l+1)
                // cube interior
                + numElx_l*numEly_l*numElz_l;
  long nfaces = numElx_l*numEly_l*(numElz_l+1)*2
                + numEly_l*numElz_l*(numElx_l+1)*2
                + numElx_l*numElz_l*(numEly_l+1)*2
                + numElx_l*numEly_l*numElz_l*6;

  long nel = numElx_l*numEly_l*numElz_l*6;

  Counts counts = {nvert, nedges, nfaces, nel};

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

  if (nfaces > INT_MAX)
  {
    std::cerr << "Error: number of faces will overflow 32-bit integer" << std::endl;
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

//  std::cout << "SIZE_MAX = " << SIZE_MAX << std::endl;
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


  std::cout << "Creating " << numElx << " by " << numEly << " by " << numElz; 
  std::cout << " mesh" << std::endl;

  // create all vertices
  double t_vert_start = get_time();
  for (int k = 0; k < (numElz+1); ++k)
  {
    for (int j = 0; j < (numEly+1); ++j)  // loop up columns (bottom to top)
    {
      for (int i = 0; i < (numElx+1); ++i) // loop over rows (left to right)
      {
//         std::cout << "creating point at x = " << x_inner << " , y = " << y_inner << " , z = " << z_inner << std::endl;
         coords_i[0] = x_inner;
         coords_i[1] = y_inner;
         coords_i[2] = z_inner;

         geo = getVertClassification(i, j, k, sizes);
//         std::cout << "model_dim = " << geo.model_dim << std::endl;
//         std::cout << "model_tag = " << geo.model_tag << std::endl;
         model_entity = m->findModelEntity(geo.model_dim, geo.model_tag);
//         std::cout << "model entity = " << model_entity << std::endl;
           
         vertices[i][j][k] = m->createVert(model_entity);
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

  double vert_telapsed = get_time() - t_vert_start;
  std::cout << "created vertices in " << vert_telapsed << " seconds" << std::endl;
//  printVerts(vertices, sizes);
//  RegionCallback cb;
//  cb.m_closure = m;
  // classify all MeshEntities on geometric face
  int model_dim = 3;
  int model_tag = 0;
  model_entity = m->findModelEntity(model_dim, model_tag);
  VertIdx start;
  apf::MeshEntity* vertices_i[4];  // hold vertices for a particular element
  // build element from verticies
  double t_el_start = get_time();
  for (int k = 0; k < numElz; ++k)
  {
    for (int j = 0; j < numEly; ++j)
    {
      for (int i = 0; i < numElx; ++i)
      {
//        std::cout << "\ncreating element " << i << ", " << j << ", " << k << std::endl;
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

  double el_telapsed = get_time() - t_el_start;
  std::cout << "created all elements in " << el_telapsed << " seconds " << std::endl;
  // build, verify  mesh
/*
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m);
  std::cout << "finished deriving model" << std::endl;
*/
  double t_verify_start = get_time();
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  checkMeshFaces(m);
  checkMesh(m, sizes, counts);
  m->verify();
  std::cout << "verified" << std::endl;
  double verify_telapsed = get_time() - t_verify_start;
  std::cout << "total verification time: " << verify_telapsed << " seconds" << std::endl;
 // for quadratic meshes
//  apf::FieldShape* linear2 = apf::getSerendipity();
//    apf::FieldShape* linear2 = apf::getSBP3Shape(1);
  apf::FieldShape* linear2 = apf::getLagrange(1);
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
  apf::writeVtkFiles("outTet", m);
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

void checkMesh(apf::Mesh* m, Sizes sizes, Counts counts)
{
  int numElx = sizes.numElx;
  int numEly = sizes.numEly;
  int numElz = sizes.numElz;

  int numVert = m->count(0);
  int ret_stat = 0;
  if ( numVert != counts.numVert)
  {
    std::cerr << "Number of vertices is incorrect: expected: ";
    std::cerr << counts.numVert << ", got " << numVert << std::endl;
    ret_stat += 1;
  }

  int numEdge = m->count(1);
  if (numEdge != counts.numEdge)
  {
    std::cerr << "Number of edges is incorrect: expected: ";
    std::cerr << counts.numEdge << ", got " << numEdge << std::endl;
    ret_stat += 1;
  }

  int numFace = m->count(2);
  if (numFace != counts.numFace)
  {
    std::cerr << "Number of faces is incorrect: expected: ";
    std::cerr << counts.numFace << ", got " << numFace << std::endl;
    ret_stat += 1;
  }

  int numEl = m->count(3);
  if (numEl != counts.numEl)
  {
    std::cerr << "Number of elements is incorrect: expected: ";
    std::cerr << counts.numEl << ", got " << numEl << std::endl;
    ret_stat += 1;
  }

  // check classification
  int classification[4][4];
  for (int i=0; i < 4; ++i)
    for (int j=0; j < 4; ++j)
      classification[i][j] = 0;

  apf::MeshEntity* e;
  apf::ModelEntity* me;
  int model_dim;
  for (int dim=0; dim < 4; ++dim)
  {
    apf::MeshIterator* it = m->begin(dim);
    while ( (e = m->iterate(it)) )
    {
      me = m->toModel(e);
      model_dim = m->getModelType(me);
      classification[dim][model_dim] = classification[dim][model_dim] + 1;
    }
  }

  // calculate expected classification
  int verts_on_verts = 8;
  int verts_on_edges = (numElx -1)*4 + (numEly-1)*4 + (numElz-1)*4;
  int verts_on_faces = (numElx-1)*(numElz-1)*2
                       +  (numEly-1)*(numElz-1)*2
                       + (numElx-1)*(numEly-1)*2;
  int verts_on_regions = (numElx-1)*(numEly-1)*(numElz-1);

  int edges_on_verts = 0;
  int edges_on_edges = numElx*4 + numEly*4 + numElz*4;
  int edges_on_faces = 2*(numElx*(numEly-1) + numEly*(numElx-1) + numElx*numEly)
                     + 2*(numEly*(numElz-1) + numElz*(numEly-1) + numEly*numElz)
                     + 2*(numElx*(numElz-1) + numElz*(numElx-1) + numElx*numElz)
                     ;
  int edges_on_regions = counts.numEdge - edges_on_verts - edges_on_edges
                         - edges_on_faces;

  int faces_on_verts = 0;
  int faces_on_edges = 0;
  int faces_on_faces = numElx*numEly*4 + numEly*numElz*4 + numElx*numElz*4;
  int faces_on_regions = counts.numFace - faces_on_faces;

  int regions_on_verts = 0;
  int regions_on_edges = 0;
  int regions_on_faces = 0;
  int regions_on_regions = counts.numEl;

  int expected_class[4][4] = { 
    {verts_on_verts, verts_on_edges, verts_on_faces, verts_on_regions},
    {edges_on_verts, edges_on_edges, edges_on_faces, edges_on_regions},
    {faces_on_verts, faces_on_edges, faces_on_faces, faces_on_regions},
    {regions_on_verts, regions_on_edges, regions_on_faces, regions_on_regions}
  };

  // check the total number of classified entities is correct
  for (int dim=0; dim < 4; ++dim)
    for (int class_dim=0; class_dim < 4; ++class_dim)
    {
      if (classification[dim][class_dim] !=  expected_class[dim][class_dim])
      {
        std::cerr << "Wrong number of dimension " << dim << " MeshEntities";
        std::cerr << " classified on model dimension " << class_dim;
        std::cerr << ", expected: " << expected_class[dim][class_dim];
        std::cerr << ", got: " << classification[dim][class_dim] << std::endl;
        ret_stat += 1;
      }
    }


  // now check the right number of entities is classified on each edge and face
  int edge0 = numElx-1;
  int edge1 = numEly-1;
  int edge2 = edge0;
  int edge3 = edge1;
  int edge4 = edge0;
  int edge5 = edge1;
  int edge6 = edge0;
  int edge7 = edge1;
  int edge8 = numElz-1;
  int edge9 = edge8;
  int edge10 = edge8;
  int edge11 = edge8;

  int expected_vert_edgeclass[] = {edge0, edge1, edge2, edge3, edge4, edge5, 
                                   edge6, edge7, edge8, edge9, edge10, edge11};

  int vert_edgeclass[12];
  for (int i=0; i < 12; ++i)
    vert_edgeclass[i] = 0;

  apf::MeshIterator* it = m->begin(0);
  int model_tag;
  while ( (e = m->iterate(it)) )
  {
    me = m->toModel(e);
    model_dim = m->getModelType(me);
    model_tag = m->getModelTag(me);
    if (model_dim == 1)
      vert_edgeclass[model_tag] = vert_edgeclass[model_tag] + 1;
  }

  for (int i=0; i <12; ++i)
  {
    if (vert_edgeclass[i] != expected_vert_edgeclass[i])
    {
      std::cerr << "incorrect number of vertices classified on model edge";
      std::cerr << i << ", expected: " << expected_vert_edgeclass[i];
      std::cerr << ", got: " << vert_edgeclass[i] << std::endl;
      ++ret_stat;
    }
  }



  // edges on edges
  edge0 = numElx;
  edge1 = numEly;
  edge2 = edge0;
  edge3 = edge1;
  edge4 = edge0;
  edge5 = edge1;
  edge6 = edge0;
  edge7 = edge1;
  edge8 = numElz;
  edge9 = edge8;
  edge10 = edge8;
  edge11 = edge8;

  int expected_edge_edgeclass[] = {edge0, edge1, edge2, edge3, edge4, edge5, 
                                   edge6, edge7, edge8, edge9, edge10, edge11};

  int edge_edgeclass[12];
  for (int i=0; i < 12; ++i)
    edge_edgeclass[i] = 0;

  // edges on faces
  int face0 = numElx*(numEly-1) + numEly*(numElx-1) + numElx*numEly;
  int face1 = numElx*(numElz-1) + numElz*(numElx-1) + numElx*numElz;
  int face2 = numEly*(numElz-1) + numElz*(numEly-1) + numEly*numElz;
  int face3 = face1;
  int face4 = face2;
  int face5 = face0;

  int expected_edge_faceclass[] = { face0, face1, face2, face3, face4, face5};
  int edge_faceclass[6];
  for (int i=0; i < 6; ++i)
    edge_faceclass[i] = 0;

  it = m->begin(1);
  while ( (e = m->iterate(it)) )
  {
    me = m->toModel(e);
    model_dim = m->getModelType(me);
    model_tag = m->getModelTag(me);
    if (model_dim == 1)
    {
      edge_edgeclass[model_tag] = edge_edgeclass[model_tag] + 1;
    } else if (model_dim == 2)
    {
      edge_faceclass[model_tag] = edge_faceclass[model_tag] + 1;
    }
  }

  for (int i=0; i < 12; ++i)
  {
    if (expected_edge_edgeclass[i] != edge_edgeclass[i])
    {
      std::cerr << "incorrect number of edges classified on model edge ";
      std::cerr << i << ", expected: " << expected_edge_edgeclass[i];
      std::cerr << ", got: " << edge_edgeclass[i] << std::endl;
      ++ret_stat;
    }
  }

  for (int i=0; i < 6; ++i)
  {
    if (expected_edge_faceclass[i] != edge_faceclass[i])
    {
      std::cerr << "incorrect number of edges classified on model face ";
      std::cerr << i << ", expected: " << expected_edge_faceclass[i];
      std::cerr << ", got: " << edge_faceclass[i] << std::endl;
      ++ret_stat;
    }
  }


  // faces on faces
  face0 = numElx*numEly*2;
  face1 = numElx*numElz*2;
  face2 = numEly*numElz*2;
  face3 = face1;
  face4 = face2;
  face5 = face0;


  int expected_face_faceclass[] = { face0, face1, face2, face3, face4, face5};
  int face_faceclass[6];
  for (int i=0; i < 6; ++i)
    face_faceclass[i] = 0;

  it = m->begin(2);
  while ( (e = m->iterate(it)) )
  {
    me = m->toModel(e);
    model_dim = m->getModelType(me);
    model_tag = m->getModelTag(me);
    if (model_dim == 2)
      face_faceclass[model_tag] = face_faceclass[model_tag] + 1;
  }

  for (int i=0; i < 6; ++i)
  {
    if ( expected_face_faceclass[i] != face_faceclass[i])
    {
      std::cerr << "incorrect number of faces classified on model face ";
      std::cerr << i << ", expected: " << expected_face_faceclass[i];
      std::cerr << ", got: " << face_faceclass[i] << std::endl;
      ++ret_stat;
    }
  }

  if (ret_stat > 0)
  {
    std::cerr << "Error: " << ret_stat << " checkMesh tests failed" << std::endl;
    abort();
  }

} // end function


void calcCentroid(apf::Mesh* m, apf::MeshEntity* e, double centroid[3])
{
  apf::Downward verts;
  m->getDownward(e, 0, verts);
  calcCentroid(m, verts, centroid);

}

void calcCentroid(apf::Mesh* m, apf::MeshEntity** verts,  double centroid[3])
{
  apf::Vector3 p1;
  apf::Vector3 p2;
  apf::Vector3 p3;
  m->getPoint(verts[0], 0, p1);
  m->getPoint(verts[1], 0, p2);
  m->getPoint(verts[2], 0, p3);

  for (int i =0; i < 3; ++i)
  {
    centroid[i] = (p1[i] + p2[i] + p3[i])/3;
  }
}

void checkMeshFaces(apf::Mesh* m)
{
  apf::MeshEntity* e;
  apf::Up up;
  apf::ModelEntity* me;
  int model_dim;
  int model_tag;
  double centroid[3];

  int facenum = 0;
  apf::MeshIterator* it = m->begin(2);

  while ( (e = m->iterate(it)) )
  {
    me = m->toModel(e);
    model_dim = m->getModelType(me);
    model_tag = m->getModelTag(me);
    m->getUp(e, up);
    if (model_dim == 2) // equal order classification
    {
      calcCentroid(m, e, centroid);
      if (up.n != 1)
      {
        std::cout << "face " << facenum << " has " << up.n;
        std::cout << " adjacent regions, expected 1";
        std::cout << ", centroid = " << centroid[0] << ", " << centroid[1];
        std::cout << ", " << centroid[2] << std::endl;
        std::cout << "model_dim = " << model_dim << ", model_tag = " << model_tag << std::endl;
        std::cout << "model entity = " << me << std::endl;
      }
    } else if (model_dim == 3)
    {
      if (up.n != 2)
      {
        std::cout << "face " << facenum << " has " << up.n;
        std::cout << " adjacent regions, expected 2";
        std::cout << ", centroid = " << centroid[0] << ", " << centroid[1];
        std::cout << ", " << centroid[2] << std::endl;
        std::cout << "model_dim = " << model_dim << ", model_tag = " << model_tag << std::endl;
        std::cout << "model entity = " << me << std::endl;
      }
    }

    ++facenum;

  }
}

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



Geom getEdgeClassification(apf::Mesh* m, VertIdx _v1, VertIdx _v2, 
                           apf::MeshEntity* verts[2], Sizes sizes)
{
  Geom _g1 = getVertClassification(_v1, sizes);
  Geom _g2 = getVertClassification(_v2, sizes);
  // sort by model dimension
  Geom g1, g2;
  VertIdx v1, v2;

  if (_g1.model_dim > _g2.model_dim)
  {
    g1 = _g2;
    g2 = _g1;
    v1 = _v2;
    v2 = _v1;
  } else
  {
    g1 = _g1;
    g2 = _g2;
    v1 = _v1;
    v2 = _v2;
  }
/*
  std::cout << "g1 dim, tag = " << g1.model_dim << ", " << g1.model_tag << std::endl;
  std::cout << "g2 dim, tag = " << g2.model_dim << ", " << g2.model_tag << std::endl;
  std::cout << "v1 idx = " << v1.i << ", " << v1.j << ", " << v1.k << std::endl;
  std::cout << "v2 idx = " << v2.i << ", " << v2.j << ", " << v2.k << std::endl;
*/

  int model_dim = -1;
  int model_tag = -1;

  if (g1.model_dim == g2.model_dim && g1.model_tag != g2.model_tag)
  {
//    std::cout << "same dimension different tag" << std::endl;
    if (g1.model_dim == 0) // vert-vert
    {
//      std::cout << "vert-vert" << std::endl;
      model_dim = vert_only_classification[g1.model_tag][g2.model_tag][0];
      model_tag = vert_only_classification[g1.model_tag][g2.model_tag][1];
      
    } else if (g1.model_dim == 1) // edge-edge
    {
//      std::cout << "edge-edge" << std::endl;
      // figure out which face it is
      bool foundface = false;
      int idx1, idx2;
      for (int i=0; i < 6; ++i)
      {
        idx1 = contains(face_edges[i], 4, g1.model_tag);
        idx2 = contains(face_edges[i], 4, g2.model_tag);
        if (idx1 >= 0 && idx2 >= 0)
        {
//          std::cout << "found face " << i << std::endl;
          model_dim = 2;
          model_tag = i;
          foundface = true;
          break;
        }
      }

      if (!foundface)
      {
        model_dim = 3;
        model_tag = 0;
      }

    } else if (g1.model_dim == 2) // face-face
    {
      model_dim = 3;
      model_tag = 0;

    } else   //region-region
    {
      std::cerr << "edge classification failed" << std::endl;
      abort();
    }
    // needed for single dimensions
  } else if (g1.model_dim == 0 && g2.model_dim == 2)
  {
    int idx = contains(face_verts[g2.model_tag], 4, g1.model_tag);
    if (idx >= 0)
    {
      model_dim = g2.model_dim;
      model_tag = g2.model_tag;
    } else
    {
      model_dim = 3;
      model_tag = 0;
    }

  } else if (g1.model_dim == 0 && g2.model_dim == 1)
  {
    // check if the vert is collinear with the edge
    int idx = contains(vert_edges[g1.model_tag], 3, g2.model_tag);
    if (idx >= 0)
    {
      model_dim = g2.model_dim;
      model_tag = g2.model_tag;
    } else  // figure out the face that is bounded by the edge and the vert
    {  // this is really slow, but this code path only gets exercised for 
       // small meshes
      bool found = false;
      for (int vert_edge=0; vert_edge < 3; ++vert_edge)
      {
        for (int face=0; face < 6; ++face)
        {
          // check if this face contains both the edge an an edge that
          // the vert could belong to
          int idx1 = contains(face_edges[face], 4, g2.model_tag);
          int idx2 = contains(face_edges[face], 4, vert_edges[g1.model_tag][vert_edge]);
          if (idx1 >= 0 && idx2 >= 0)
          {
            model_dim = 2;
            model_tag = face;
            found = true;
          }
        }
      }

      if (!found)
      {
        model_dim = 3;
        model_tag = 0;
      }

    }

  } else if ( g1.model_dim == 1 && g2.model_dim == 2)
  {
    // check if face contains the edge
    // if yes, getMaxGeometry
    // if no, classify on region
    int idx;
    idx = contains(face_edges[g2.model_tag], 4, g1.model_tag);
    if (idx >= 0)
    {
      model_dim = g2.model_dim;
      model_tag = g2.model_tag;
    } else
    {
      model_dim = 3;
      model_tag = 0;
    }

  } else  // general case
  {
//    std::cout << "general case " << std::endl;
    Geom g_class[] = {g1, g2};
    Geom g = getMaxGeometry(m, g_class, 2);
    model_dim = g.model_dim;
    model_tag = g.model_tag;
  }

//  std::cout << "model_tag = " << model_tag << ", model_dim = " << model_dim << std::endl;
  Geom g = {model_tag, model_dim};
  return g;

}

// get the classification of an already existing edge
Geom getEdgeClassification(apf::Mesh* m, apf::MeshEntity* verts[2])
{

  Geom geo;
  apf::MeshEntity* e = apf::findUpward(m, apf::Mesh::EDGE, verts);
  apf::ModelEntity* me = m->toModel(e);
  geo.model_dim = m->getModelType(me);
  geo.model_tag = m->getModelTag(me);
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
  bool matched = true;
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
    { {0, 0, 0}, {0, 0, 1} },
    { {0, 0, 0}, {1, 0, 0} },
    { {0, 0, 0}, {0, 1, 0} },
  
    // i = 0 edges
    { {0, 1, 0}, {0, 1, 1} },
    { {0, 0, 1}, {0, 1, 1} },
    { {0, 0, 0}, {0, 1, 1} },
    // j = 0 edges
    { {1, 0, 0}, {1, 0, 1} },
    { {0, 0, 1}, {1, 0, 1} },
    { {0, 0, 0}, {1, 0, 1} },
    // k = 0 edges
    { {1, 0, 0}, {1, 1, 0} },
    { {0, 1, 0}, {1, 1, 0} },
    { {0, 0, 0}, {1, 1, 0} },
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
  if (i == 0 && j == 0)
    create_edge[0] = true;
  if (j == 0 && k == 0)
    create_edge[1] = true;
  if ( i == 0 && k == 0)
    create_edge[2] = true;
  if ( i == 0 )
  {
    create_edge[3] = true;
    create_edge[4] = true;
    create_edge[5] = true;
  }
  if ( j == 0)
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
//  std::cout << "\ncreating edges" << std::endl;
//  std::cout << "start idx = " << start.i << ", " << start.j << ", " << start.k << std::endl;
  int idx;
  VertIdx v1, v2, vertidx;
  apf::MeshEntity* verts_i[2];
  Geom g; // geometric classification of the edge
  apf::ModelEntity* model_entity;
  for (int i = 0; i < NEDGES; ++i)
  {
//    std::cout << "\nconsidering edge " << i << std::endl;
    idx = contains(rear_edge_idx, NREAREDGES, i);

    if (idx >= 0)  // if found
    {
//      std::cout << "this is a rear edge" << std::endl;
      if (!create_edge[idx])
      {
//        std::cout << "skipping this edge" << std::endl;
        continue;
      }
    }

//    std::cout << "creating edge" << std::endl;

    // if we get here, we should create the edge
    add(start, edges[i][0], v1);
    add(start, edges[i][1], v2);
//    v1 = add(start, edges[i][0]);
//    v2 = add(start, edges[i][1]);

    for (int j=0; j < 2; ++j)
    {
//      vertidx = add(start, edges[i][j]);
      add(start, edges[i][j], vertidx);
      verts_i[j] = getVert(vertidx, verts);
    }

    g = getEdgeClassification(m, v1, v2, verts_i, sizes);

    model_entity = getModelEntity(m, g);
    createEdge(m, verts_i, model_entity);
  }  

}


void createEdge(apf::Mesh2* m, apf::MeshEntity* edge_verts[2], apf::ModelEntity* model_entity) 
{
  /*
  apf::MeshEntity* existing = apf::findUpward(m, apf::Mesh::EDGE, edge_verts);
  if (existing) {
    std::cerr << "createEdge error: edge already exists!\n";
    abort();
  }
  */
  m->createEntity(apf::Mesh::EDGE, model_entity, edge_verts);

}

int tri_edges[3][2] = { {0, 1}, {1, 2}, {0, 2} };

void createFaces(apf::Mesh2* m, Sizes sizes, VertIdx start, apf::MeshEntity**** verts)
{
//  std::cout << "\ncreating faces" << std::endl;
  int idx;
  VertIdx vertidx;
  apf::MeshEntity* verts_i[3];
  apf::MeshEntity* edge_verts[2];
  Geom g_class[3]; // geometric classification of the vertices
  Geom g;  // geometric classification of the face
//  double centroid[3];
  apf::ModelEntity* model_entity;
  for (int i = 0; i < NFACES; ++i)
  {
//    std::cout << "\nconsidering face " << i << std::endl;
    idx = contains(rear_face_idx, 6, i);

    if (idx >= 0)  // if found
    {
//      std::cout << "face is a rear face" << std::endl;
      if (!create_face[idx])
      {
//        std::cout << "not creating face" << std::endl;
        continue;
      }
    }
//    std::cout << "creating face " << std::endl;
    // if we get here, we should create the face
    for (int j=0; j < 3; ++j)
    {
      add(start, faces[i][j], vertidx);
//      vertidx = add(start, faces[i][j]);
      verts_i[j] = getVert(vertidx, verts);
//      std::cout << "vert " << j << " idx " << vertidx.i << ", " << vertidx.j;
//      std::cout << ", " << vertidx.k << std::endl;
    }

    for (int j=0; j < 3; ++ j)
    {
      edge_verts[0] = verts_i[tri_edges[j][0]];
      edge_verts[1] = verts_i[tri_edges[j][1]];
      g_class[j] = getEdgeClassification(m, edge_verts);
//      std::cout << "edge " << j << " dim, tag = " << g_class[j].model_dim;
//      std::cout << ", " << g_class[j].model_tag << std::endl;
    }

//    calcCentroid(m, verts_i, centroid);
//    std::cout << "creating face with centroid = " << centroid[0];
//    std::cout << ", " << centroid[1] << ", " << centroid[2] << std::endl;

    g = getMaxGeometry(m, g_class, 3);
/*
    if (g.model_dim < 2)
    {
      std::cout << "Warning: classifying face on model dim " << g.model_dim << std::endl;
    }
*/
    model_entity = getModelEntity(m, g);
    createFace(m, verts_i, model_entity);
  }  

}

void createFace(apf::Mesh2* m, apf::MeshEntity* face_verts[3], apf::ModelEntity* model_entity) 
{
/*
  double centroid[3];
  calcCentroid(m, face_verts, centroid);

  std::cout << "face_verts = " << face_verts[0] << ", " << face_verts[1];
  std::cout << ", " << face_verts[2] << std::endl;
  std::cout << "model entity = " << model_entity << std::endl;
  std::cout << "centroid = " << centroid[0] << ", " << centroid[1] << ", ";
  std::cout << centroid[2] << std::endl;
*/
/*
  apf::MeshEntity* e;
  FaceCallback cb;
  cb.m_closure = m;
*/
  apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, face_verts);

//  model_entity = m->toModel(e);
}




// get the geometry that the edge defind by g1 and g2 should be classified
// on
Geom getMaxGeometry(apf::Mesh* m, Geom g_class[], const int ng)
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

  Geom geo = {min_tag, max_dim};
  return geo;
}

apf::ModelEntity* getModelEntity(apf::Mesh* m, Geom g_class)
{
  int max_dim = g_class.model_dim;
  int min_tag = g_class.model_tag;
  // max_dim and min_tag now fully describe the geometric entity
  apf::ModelEntity* me = m->findModelEntity(max_dim, min_tag);
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
    add(start, tet[i], pos);
//    pos = add(start, tet[i]);
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

int contains(const int vals[], const int len, int  val)
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

void add(const VertIdx v1, const int arr[3], VertIdx &vout)
{
  vout.i = v1.i + arr[0];
  vout.j = v1.j + arr[1];
  vout.k = v1.k + arr[2];
}

bool isequal(Geom g1, Geom g2)
{
  return g1.model_dim == g2.model_dim && g1.model_tag == g2.model_tag;
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

double get_time()
{
  struct timeval t;
  gettimeofday(&t, NULL);
  return t.tv_sec + t.tv_usec*1e-6;
}
