#include <apf.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfNumbering.h>
#include <apfShape.h>

//#include "funcs1.h"
#include "apfSBPShape.h"

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
  // Problem 3
  std::cout << "Problem 3 output: " << std::endl;
  // specify number of cubes in each direction
  // each cube is divided into 6 tets
  // the subdivions is the sameon used for Marching Tetrahedra
  // projected onto the xyz planes, you get two triangles in x and y
  // and one in z
  int numElx = 1;
  int numEly = 1; 
  int numElz = 1;
  double xrange = 2.0;
  double yrange = 2.0;
  double zrange = 2.0;
  double x_spacing = xrange/numElx;  // spacing of el
  double y_spacing = yrange/numEly;
  double z_spacing = zrange/numElz;
  double x_0 = -1.0;  // x coordinate of lower left corner of starting element
  double y_0 = -1.0 ; // y coordinate of lower left corner of starting element
  double z_0 = -1.0;  // z coordinate of lower left corning of starting elemnt
  double x_i = x_0;  // current vertex coordinates
  double y_i = y_0;
  double z_i = z_0;
  apf::MeshEntity* vertices[numElx+1][numEly+1][numElz+1];  // hold pointers to all vertices
  apf::MeshEntity* tet1[4];  // hold vertices for a particular element
  apf::MeshEntity* tet2[4];
  apf::MeshEntity* tet3[4];
  apf::MeshEntity* tet4[4];
  apf::MeshEntity* tet5[4];
  apf::MeshEntity* tet6[4];
  apf::MeshEntity** all_tets[6];  // array of tet1 - tet5
  apf::Vector3 coords_i(0.0,0.0,0.0);  // hold coordinates of each point
//  apf::FieldShape* linear2 = apf::getSBPQuadratic();

  // for linear meshes
  apf::FieldShape* linear2 = apf::getLagrange(1);
  std::cout << "shape name = " << linear2->getName() << std::endl;
  apf::changeMeshShape(m, linear2, true);



  // create all vertices
  for (int k = 0; k < (numElz+1); ++k)
  { 
    for (int j = 0; j < (numEly+1); ++j)  // loop up columns (bottom to top)
    {
      for (int i = 0; i < (numElx+1); ++i) // loop over rows (left to right)
      {
         std::cout << "creating point at x = " << x_i << " , y = " << y_i << std::endl;
         coords_i[0] = x_i;
         coords_i[1] = y_i;
         coords_i[2] = z_i;
         
         vertices[i][j][k] = m->createVert(0);
         m->setPoint(vertices[i][j][k], 0, coords_i);
         x_i = x_i + x_spacing;
      }
      y_i = y_i + y_spacing;  // increment y_i
      x_i = x_0;  // reset x_i to beginning of row
    } 
    // increment z, reset others
    z_i = z_i + z_spacing;
    x_i = x_0;
    y_i = y_0;
  } 

  int elnum = 1;
    // build element from verticies
  for (int k = 0; k < numElz; ++k)
  {
    for (int j = 0; j < numEly; ++j)
    {
      for (int i = 0; i < numElx; ++i)
      {

        std::cout << "creating tets on cube " << elnum << std::endl;
        
  //      std::cout << "creating element i = " << i << " , j = " << j << std::endl;
        // get correct vertices
        tet1[0] = vertices[i][j][k];
        tet1[1] = vertices[i+1][j][k];
        tet1[2] = vertices[i+1][j+1][k+1];
        tet1[3] = vertices[i+1][j][k+1];
        
        tet2[0] = vertices[i][j][k];
        tet2[1] = vertices[i+1][j+1][k+1];
        tet2[2] = vertices[i+1][j][k+1];
        tet2[3] = vertices[i][j][k+1];

        tet3[0] = vertices[i][j][k];
        tet3[1] = vertices[i+1][j][k];
        tet3[2] = vertices[i+1][j+1][k];
        tet3[3] = vertices[i+1][j+1][k+1];

        tet4[0] = vertices[i][j][k];
        tet4[1] = vertices[i][j+1][k];
        tet4[2] = vertices[i+1][j+1][k];
        tet4[3] = vertices[i+1][j+1][k+1];

        tet5[0] = vertices[i][j][k];
        tet5[1] = vertices[i][j+1][k];
        tet5[2] = vertices[i+1][j+1][k+1];
        tet5[3] = vertices[i][j+1][k+1];

        tet6[0] = vertices[i][j][k];
        tet6[1] = vertices[i+1][j+1][k+1];
        tet6[2] = vertices[i][j+1][k+1];
        tet6[3] = vertices[i][j][k+1];

       all_tets[0] = tet1;
       all_tets[1] = tet2;
       all_tets[2] = tet3;
       all_tets[3] = tet4;
       all_tets[4] = tet5;
       all_tets[5] = tet6;
       
       for (int n = 0; n < 6; ++n)
       {
        std::cout << "creating local tet number " << n << std::endl;
        apf::buildElement(m, 0, apf::Mesh::TET, all_tets[n]);
       }




//        std::cout << "Element " << el_num << " has verticies "; 
//        std::cout << i + ((numElx+1)*j) << " , " << i+1 + ((numElx+1)*j) << " , " << i+1 + ((numElx+1)*(j+1));
//        std::cout << " , " << i + ((numElx+1)*(j+1)) << std::endl;
  //      std::cout << "assembled vertices" << std::endl;
  //
        // counterclockwise ordering
//        apf::buildElement(m, 0, apf::Mesh::TET, tet1);

//        vertices_i[0] = vertices[i+1][j];
//        vertices_i[1] = vertices[i+1][j+1];
//        vertices_i[2] = vertices[i][j+1];
//        apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vertices_i);

  /*
        vertices_i[1] = vertices[i][[j+1];
        apf::buildElement(m, 0, apf::Mesh::TRIANGLE, vertices_i);
  */

       ++elnum;
       }
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
  std::cout << "deriving model" << std::endl;
  apf::deriveMdsModel(m);
  std::cout << "finished deriving model" << std::endl;
  m->acceptChanges();
  std::cout << "accepted changes" << std::endl;
  m->verify();
  std::cout << "verified" << std::endl;

 // for quadratic meshes
//  apf::FieldShape* linear2 = apf::getSerendipity();
//  apf::FieldShape* linear2 = apf::getSBPShape(1);
//  apf::FieldShape* linear2 = apf::getLagrange(2);
//  apf::changeMeshShape(m, linear2, true);  // last argument should be true for second order

  std::cout << "changed mesh shape" << std::endl;
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

