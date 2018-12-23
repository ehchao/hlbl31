#define MAIN_PROGRAM
//#include <sstream>
//#include <cmath>
//#include <vector>
//#include "getL.h"
#include <mpi.h>
#include "observer.h"
//#include "qdp.h"
#include "contraction_3_1.h"

using namespace QDP;
using namespace std;

//fix coordinates wrt the origin
void fix_coord(int crd[4], int origin[4]){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < 4 ; i++){crd[i] = (crd[i]-origin[i]) % lattsize[i];}
}

void fix_coord(multi1d<int> & crd, const multi1d<int> & origin, int length){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < length ; i++){crd[i] = (crd[i]-origin[i]) % lattsize[i];}
}

int main(int argc, char *argv[]){

    observer_init(argc, argv);

    QDPIO::cout<< "observer init done" <<endl;
    read_hlbl31_infile(argc, argv);
    
    QDPIO::cout << "read done" << endl;

    int status[3];
    
    

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    QDPIO::cout<< "Testing the kernel writer ... "

    char* name1= "3L_kernel_hdf5";
    char* name2 = "L_x_kernel_hdf5";

    multi4d<LatticeComplexD> data(4,4,4,4);
    //TODO continue from here
    for(int i =0; i<4; i++){
      for(int j =0; j <4; j++){
        for(int k= 0 ; k<4; k++){
          for(int l =0; l < 4; l++){ 
            data[i][j][k][l] = cmplx(Real(i+j), Real(k*l));
          }
        }
      }
    }

    write_3L_kernel_hdf5(hdf5_writer &writer, name1, const multi4d<LatticeComplexD> & data);
    
    write_L_x_kernel_hdf5(hdf5_writer &writer, char *name, const multi4d<LatticeComplexD> & data);

    QDPIO::cout <<"Hello world" <<endl;
    Layout::destroy();
    QDP_finalize();
  
   }
   return 0;
}

