#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdbool.h>
#include <limits.h>
#include <cassert>

// function finds the line starting with sentinel, which must be of the form
// sentinel= 123  and returns the integer value to the right of the equals
// sign
// The input stream is left on that line, so the next call to std::getline
// will get the first data line of the section
int getSectionCount(std::ifstream& ifs, std::string sentinel, bool reset_stream = true)
{

  // rewind the stream back to the beginning
  // this is potentially inefficient, but ensures the sentinel will be
  // found regardless of the order the sections appear in the mesh file
  assert (ifs.tellg() != -1);
  if ( reset_stream )
  {
    std::cout << "currently at position " << ifs.tellg() << std::endl;
//    ifs.clear();
    ifs.seekg(0);
  }


  // find the NPOIN line
  std::string currline;
//  std::string sentinel ("NPOIN");
  char numel_str[256];
  int numel = 0;
  while ( std::getline(ifs, currline) )
  {
//    std::cout << "line = " << currline << std::endl;
    // if line starts with sentinel
    if ( currline.compare(0, sentinel.length(), sentinel) == 0 )
    {
      // find delimiter
      std::string::size_type idx = currline.find("=");
//      std::cout << "equals sign at index " << idx << std::endl;

      // extract the number after the equals sign
      for (size_t i=idx+1; i < currline.length(); ++i)
        numel_str[i - idx - 1] = currline[i];

      sscanf(numel_str, "%d", &numel);
      break;
    }  // end if

  }  // end while

  return numel;
}  // end function

std::vector<apf::MeshEntity*> createVerts(std::ifstream& ifs, apf::Mesh2* m)
{

  std::string sent ("NPOIN");
  int numVert = getSectionCount(ifs, sent);
  std::vector<apf::MeshEntity*> verts(numVert);
  std::string currline;

  double x, y;
  int vertnum;
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point

  int model_dim = 2;  // classify everything on the interior of the face
  int model_tag = 0;
  apf::ModelEntity* model_entity = m->findModelEntity(model_dim, model_tag);

  for (int i=0; i < numVert; ++i)
  {
    std::getline(ifs, currline);

    sscanf(currline.c_str(), "%le %le %d", &x, &y, &vertnum);
    coords_i[0] = x;
    coords_i[1] = y;

    verts[vertnum] = m->createVert(model_entity);
    m->setPoint(verts[vertnum], 0, coords_i);
  }

  return verts;

}

void createElements(std::ifstream& ifs, apf::Mesh2* m, std::vector<apf::MeshEntity*> verts)
{

  // TODO: assert mesh is 2 dimensional (NDIME)
  //
  std::string sent ("NELEM");
  std::string currline;
  int numel = getSectionCount(ifs, sent);
  int model_dim = 2;  // classify everything on the interior of the face
  int model_tag = 0;
  apf::ModelEntity* model_entity = m->findModelEntity(model_dim, model_tag);

  int eltype, vnum1, vnum2, vnum3, elnum;
  apf::MeshEntity* vertices[3];

  for (int i=0; i < numel; ++i)
  {
    std::getline(ifs, currline);
    sscanf(currline.c_str(), "%d %d %d %d %d", &eltype, &vnum1, &vnum2, &vnum3, &elnum);
    assert( eltype == 5);  // triangular element

    vertices[0] = verts[vnum1];
    vertices[1] = verts[vnum2];
    vertices[2] = verts[vnum3];

    apf::buildElement(m, model_entity, apf::Mesh::TRIANGLE, vertices);

  }

}  // end function

// reclassifies the boundary edges to be dimension 1 geometric entities
// Specifically, each of NMARK boundary conditions gets its own geometry,
// numbered from 0 to NMARK-1, in the order they appear in the file
void createBoundaries(std::ifstream& ifs, apf::Mesh2* m, std::vector<apf::MeshEntity*> verts)
{


  std::string sent ("NMARK");
  int nmark = getSectionCount(ifs, sent);
  std::cout << "nmark = " << nmark << std::endl;
  std::string sent2 ("MARKER_ELEMS");
  std::string currline;
  apf::Up edges1;
  apf::Up edges2;

  for (int geo_num=0; geo_num < nmark; ++geo_num)
  {
    std::cout << "geo_num = " << geo_num << std::endl;

    // get the model entity for this geometry
    apf::ModelEntity* model_entity = m->findModelEntity(1, geo_num);
    // seek without resetting the stream
    int nedges = getSectionCount(ifs, sent2, false);
    std::cout << "nedges = " << nedges << std::endl;
    int eltype, vnum1, vnum2;
    apf::MeshEntity *v1, *v2;
    for (int i=0; i < nedges; ++i)
    {
      std::cout << "i = " << i << std::endl;
      std::getline(ifs, currline);
      sscanf(currline.c_str(), "%d %d %d", &eltype, &vnum1, &vnum2);
      assert(eltype == 3);  // edge

      v1 = verts[vnum1];
      v2 = verts[vnum2];

      // get the edge connecting the two vertices
      m->getUp(v1, edges1);
      m->getUp(v2, edges2);
      apf::MeshEntity* edge;
      bool foundEdge = false;
      // N^2 search, but N is small
      for (int edge1=0; edge1 < edges1.n; ++edge1)
        for (int edge2=0; edge2 < edges2.n; ++edge2)
          if ( edges1.e[edge1] == edges2.e[edge2])
          {
            edge = edges1.e[edge1];
            foundEdge = true;
          }


      assert(foundEdge);
      // change the edge classification
      m->setModelEntity(edge, model_entity);
      // change the vertex classification too
      m->setModelEntity(v1, model_entity);
      m->setModelEntity(v2, model_entity);
    }  // end loop i
  }  // end loop geo_num

}  // end function



int main(int argc, char** argv)
{

  if (argc != 2)
  {
    std::cerr << "Usage: " << argv[0] << " fname.su2" << std::endl;
    return 1;
  }

  // open the vertex coordinate file
  std::ifstream ifs(argv[1], std::ifstream::in);
  if (!ifs.good())
  {
    std::cerr << "Error: coords.dat does not exist" << std::endl;
    return 1;
  }

  MPI_Init(&argc,&argv);
  PCU_Comm_Init();

  std::string sent ("NPOIN");
  int numel = getSectionCount(ifs, sent);
  std::cout << "numel = " << numel << std::endl;

  int numel2 = getSectionCount(ifs, sent);
  std::cout << "numel2 = " << numel2 << std::endl;

  // create the mesh
  gmi_register_null();
  gmi_model* g = gmi_load(".null");
  bool isMatched = false;
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 2, isMatched);

  auto verts = createVerts(ifs, m);
  createElements(ifs, m, verts);
  std::cout << "creating boundaries" << std::endl;
  createBoundaries(ifs, m, verts);
  std::cout << "finished creating boundaries" << std::endl;

  // finish the mesh
  std::cout << "changed mesh shape" << std::endl;
  apf::FieldShape* m_shape = m->getShape();
  std::cout << "mesh shape name = " << m_shape->getName() << std::endl;

 
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  m->verify();
  std::cout << "verified" << std::endl;

  apf::writeVtkFiles("outTri", m);
  m->writeNative("./meshfiles/abc.smb");

  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();


  return 0;
}
