#include "qdp.h"
#include "hdf5_output.h"
#include "sftmom.h"
#include "qdpinterface.h"
extern "C"
{
#include "observ.h"
#include "apps.h"
}
using namespace QDP;

/*void init_3L_kernel_hdf5(hdf5_writer &writer){
    multi1d<int> lattsize = Layout::lattSize();
    int dims[8] = {6,4,4,4, lattsize[0], lattsize[1], lattsize[2], lattsize[3]}; 
    char* axis_label[8] = {"rhosig", "mu", "nu", "lambda", "x", "y", "z", "t"};
    
    //char* axis_info[8];
    //axis_info[0] = "rhosig = [rho, sigma]";
    //for(int i = 1; i < 8; i++){axis_info[i] = NULL;}
    writer.write_real_init("3L_kernel", 8, dims , axis_label);
}*/

/*void init_L_x_kernel_hdf5(hdf5_writer &writer){
    multi1d<int> lattsize = Layout::lattSize();
    int dims[8] = {6,4,4,4, lattsize[0], lattsize[1], lattsize[2], lattsize[3]}; 
    char* axis_label[8] = {"sigma", "mu", "nu", "lambda", "x", "y", "z", "t"};

    writer.write_real_init("L_x_kernel", 8, dims , axis_label);
}
*/

//y fixed in the init of writer
/*void write_3L_kernel_hdf5(hdf5_writer &writer, int x, int y , int z, int t, const multi4d<RealD> & data){
    int dims_offset[8] = {0,0,0,0, x, y, z, t};
    writer.write_data("3L_kernel", dims_offset, 0, 1, 2, 3, data);
}*/

void write_3L_kernel_hdf5_lattice(hdf5_writer &writer, const char * filename, double  ******** data){
    //StopWatch timer;
    //timer.reset();
    //timer.start();
    int dims_offset[8] = {0,0,0,0, 0, 0, 0, 0};
    writer.write_lattice_data_parallel(filename, "3L_kernel", dims_offset, 6, 4, 4, 4, data);
    //timer.stop();
    //QDPIO::cout<<"writing finished ( "<<timer.getTimeInSeconds() << " s)" <<endl;
}

void write_L_x_kernel_hdf5_lattice(hdf5_writer &writer,const char *filename, double  ******** data){
    //StopWatch timer;
    //timer.reset();
    //timer.start();
    int dims_offset[8] = {0,0,0,0, 0, 0, 0, 0};
    writer.write_lattice_data_parallel(filename, "L_x_kernel", dims_offset, 4, 4, 4, 4, data);
    //timer.stop();
    //QDPIO::cout<<"writing finished ( "<<timer.getTimeInSeconds() << " s)" <<endl;
}

//y fixed in the init
//use sigma instead of rhosig here
/*void write_L_x_kernel_hdf5(hdf5_writer &writer, int x,  int y, int z, int t,  const multi4d<ComplexD> & data){
    int dims_offset[8] = {0,0,0,0,x,y,z,t};
    writer.write_data("L_x_kernel", dims_offset, 0, 1, 2, 3, data);
}*/

//the z*3pt corr is stored as a function of y
void init_int_z_3pt_loop_hdf5(hdf5_writer &writer){
    multi1d<int> lattsize = Layout::lattSize();
    int dims[8] = {4,4,4,4, lattsize[0], lattsize[1], lattsize[2], lattsize[3]}; 
    char* axis_label[8] = {"rho", "nu", "lambda", "sigma", "y0", "y1", "y2", "y3"};

    writer.write_real_init("int_z_3pt_y0z", 8, dims , axis_label);
}


void write_int_z_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const multi4d<ComplexD> & data){
    int dims_offset[8] = {0,0,0,0,y0,y1,y2,y3};
    writer.write_data("int_z_3pt_y0z", dims_offset, 0, 1, 2, 3, data);
}

//the 3pt corr is stored as a function of y
void init_int_3pt_loop_hdf5(hdf5_writer &writer){
    multi1d<int> lattsize = Layout::lattSize();
    int dims[7] = {4,4,4, lattsize[0], lattsize[1], lattsize[2], lattsize[3]}; 
    char* axis_label[7] = {"nu", "lambda", "sigma", "y0", "y1", "y2", "y3"};

    writer.write_real_init("int_3pt_y0z", 7, dims , axis_label);
}


void write_int_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const multi3d<ComplexD> & data){
    int dims_offset[7] = {0,0,0,y0,y1,y2,y3};
    writer.write_data("int_3pt_y0z", dims_offset, 0, 1, 2, data);
}

//function of y
//the only index is rhosig
void init_int_kernel_3pt_loop_hdf5(hdf5_writer & writer){
    multi1d<int> lattsize = Layout::lattSize();
    int dims[5] = {6, lattsize[0], lattsize[1], lattsize[2], lattsize[3]}; 
    char* axis_label[5] = {"rhosig", "y0", "y1", "y2", "y3"};

    writer.write_real_init("int_L_3pt_xy0", 5, dims , axis_label);
}


void write_int_kernel_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const multi1d<ComplexD> & data){
    int dims_offset[5] = {0,y0,y1,y2,y3};
    writer.write_data("int_L_3pt_xy0", dims_offset, 0, data);
}
