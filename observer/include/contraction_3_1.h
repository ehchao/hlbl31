#ifndef CONTRACTION_3_1_H
#define CONTRACTION_3_1_H
#include <vector>
#include <qdp.h>
extern "C"{
  #include "KQED.h"
}

#ifdef MAIN_PROGRAM
extern int Nf=2;
#endif

//typedef for functions which take a double and return a double
typedef double (*funcD)(double a);

void fix_coord(int out[4], const int crd[4], const int origin[4]);

void fix_coord(QDP::multi1d<int> & out, const QDP::multi1d<int> & crd, const QDP::multi1d<int> & origin);

//template<typename T>
extern void multi1d2array(int * out, const QDP::multi1d<int> & input, int length );

extern void array2multi1d(QDP::multi1d<int> & out, int* input, int length);

void callKernelQED(double *x, double* y, double kern[6][4][4][4]);

void callIpihatFermLoop_antisym(double *x, double *y, double ipihat[6][4][4][4]);

void callInitKernelQED(void);

void kernel_3L_array(double out[6][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4]);

void kernel_L_x_array(double out[4][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4]);

void kernel_3L_new(double out[6][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4], struct QED_kernel_temps *kert);

void kernel_L_x_new(double out[4][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4], struct QED_kernel_temps *kert);

void array2multi4dlattice(QDP::multi4d<QDP::LatticeRealD> & out, double in[6][4][4][4], int size0, int size1, int size2, int size3, const QDP::multi1d<int> & coord);

void write_im_3_pt_llc(QDP::multi3d<QDP::LatticeRealD> &out, const QDP::multi1d<QDP::LatticeColorMatrixD> &u, const QDP::LatticePropagatorD &s_0, const QDP::LatticePropagatorD &s_y, int nu, int lambda, int sigma, const QDP::multi1d<int> & coord_y);

void write_3_pt_lll(QDP::multi3d<QDP::LatticeComplexD> &out, const QDP::LatticePropagatorD &s_0, const QDP::LatticePropagatorD &s_y, int mu, int nu, int lambda, const QDP::multi1d<int> & coord_y);

void compute_propagator(QDP::LatticePropagatorD &out, QDP::multi1d<int> &src_pos, int propid);

void write_sin_func(QDP::LatticeRealD &out, int sigma);

//void int_3pt_times_z(double out[4][4][4][4], const QDP::multi3d<QDP::LatticeRealD> & in);

extern void read_hlbl31(void);

extern void read_hlbl31_infile(int argc, char * argv[]);

extern double check_3_pt_llc(const QDP::multi1d<QDP::LatticeColorMatrixD> & u, const QDP::LatticePropagatorD & s_0, const QDP::LatticePropagatorD &s_y, int nu, int lambda, const QDP::multi1d<int> & crd_y);

double check_3_pt_lll(const QDP::multi1d<QDP::LatticeColorMatrixD> & u, const QDP::LatticePropagatorD & s_0, const QDP::LatticePropagatorD & s_y, const QDP::multi1d<int> & crd_x, const QDP::multi1d<int> & crd_y, const QDP::multi1d<int> & origin, double am, double a);

extern double check_3_pt_lll_wi(const QDP::LatticePropagatorD & s_0, const QDP::LatticePropagatorD &s_y, int nu, int lambda, const QDP::multi1d<int> & crd_y);

double id_func(double a);

double sin_func(double a);

void int_3pt_times_func(double out[4][4][4][4], QDP::multi3d<QDP::LatticeRealD> & corr, funcD f, int origin[4]);

void init_3L_kernel_hdf5(hdf5_writer &writer);

void init_L_x_kernel_hdf5(hdf5_writer &writer);

void write_3L_kernel_hdf5(hdf5_writer &writer,  int x, int y , int z, int t, const QDP::multi4d<QDP::RealD> & data);
        
void write_3L_kernel_hdf5_lattice(hdf5_writer &writer, const char * name, double ******** data);
void write_L_x_kernel_hdf5_lattice(hdf5_writer &writer,const char *filename, double  ******** data);

void write_L_x_kernel_hdf5(hdf5_writer &writer,  int x, int y, int z, int t,  const QDP::multi4d<QDP::ComplexD> & data);
void init_int_z_3pt_loop_hdf5(hdf5_writer &writer);
void write_int_z_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const QDP::multi4d<QDP::ComplexD> & data);
void init_int_3pt_loop_hdf5(hdf5_writer &writer);
void write_int_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const QDP::multi3d<QDP::ComplexD> & data);
void init_int_kernel_3pt_loop_hdf5(hdf5_writer & writer);
void write_int_kernel_3pt_loop_hdf5(hdf5_writer &writer, int y0, int y1, int y2, int y3,  const QDP::multi1d<QDP::ComplexD> & data);

void multi4dlattice2multi4d(QDP::multi4d<QDP::ComplexD> & out, const QDP::multi4d<QDP::LatticeRealD> & in, int size0, int size1, int size2, int size3, const QDP::multi1d<int> & coord);
void int_L_3pt_loop(double out[6], const QDP::multi3d<QDP::LatticeRealD> & corr_3pt, const QDP::multi1d<int> & origin, const QDP::multi1d<int> & coord_y, double a);
void kernel_L(QDP::multi4d<QDP::LatticeRealD> &out, const int origin[4], const int crd_y[4], double a);//, double kern[6][4][4][4] );


//template<class T>
//void int_times_z( T & out, const T & in, int rho);


#endif
