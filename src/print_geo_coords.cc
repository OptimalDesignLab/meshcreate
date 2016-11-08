
#include <apf.h>
#include <gmi_mesh.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apfMesh2.h>
#include <PCU.h>
#include <apfShape.h>
#include <cstdlib>
#include <iostream>

double functional(apf::Vector3 coords)
{
  return coords.x()*coords.x() + coords.y()*coords.y();
}

void printGeoCoords(apf::Mesh* m, int geo_dim, int geo_tag)
{


  std::cout << "printing coordinates of nodes on geometric entity " << geo_tag << " of dimension " << geo_tag << std::endl;

  apf::ModelEntity* me;
  int model_dim;
  int model_tag;

  apf::MeshEntity * e;
  int entity_num = 0;

  apf::Vector3 coords;

  apf::FieldShape* fshape = m->getShape();
  int nnodes;  // number of nodes on current element
  int e_type;  // type of current element

  for (int dim = 0; dim <= m->getDimension(); ++dim)
  {
    std::cout << "checking mesh entities of dimension " << dim << std::endl;
    apf::MeshIterator* it = m->begin(dim);
    entity_num = 0;


    while ( (e = m->iterate(it)) )  // loop over entities of current dimension
    {
      me = m->toModel(e);
      model_dim = m->getModelType(me);
      model_tag = m->getModelTag(me);

      if ( model_dim == geo_dim && model_tag == geo_tag)
      {
        e_type = m->getType(e);
        nnodes = fshape->countNodesOn(e_type);

        if (nnodes > 0)
          std::cout << "entity " << entity_num << std::endl;

        for (int i = 0; i < nnodes; ++i)
        {
          m->getPoint(e, i, coords);
          std::cout << "  node " << i << " cordinates = (" << coords.x() << ", " << coords.y() << "), functional = " << functional(coords) << std::endl;

        }  // end loop over nodes


      }  // end if match
        

      entity_num++;

    }  // end loop over entities

  }  // end loop dim

}  // end function

int main(int argc, char** argv)
{
  MPI_Init(&argc,&argv);
  PCU_Comm_Init();
  if ( argc != 5 ) {
    if ( !PCU_Comm_Self() )
      printf("Usage: %s <model> <mesh> <geomtric dimension> <geometric tag>\n", argv[0]);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
  gmi_register_mesh();
  gmi_register_null();
  apf::Mesh2* m = apf::loadMdsMesh(argv[1],argv[2]);

  int geo_dim, geo_tag;
  sscanf(argv[3], "%d", &geo_dim);
  sscanf(argv[4], "%d", &geo_tag);


  printGeoCoords(m, geo_dim, geo_tag);

//  apf::FieldShape* fshape = apf::getLagrange(1);
//  apf::changeMeshShape(m, fshape, false);
//  apf::writeVtkFiles(argv[3], m);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}


