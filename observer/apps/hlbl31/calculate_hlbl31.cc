#include <sstream>
#include <cmath>
#include <vector>
#include "hdf5_output.h"
#include "contraction_3_1.h"
#include "observer.h"
#include <mpi.h>
#include <qdp.h>
//#define FERML
extern "C"
{
#include "apps.h"
#include "getL.h"
#include "KQED.h"
}
#include "qdpinterface.h"
//#include "contraction_3_1.h"
#ifndef PI
#define PI 3.14159265358979323846
#endif 

using namespace QDP;

void multi1d2array(int * out, const multi1d<int> & input, int length ){
    for(int i = 0; i < length; i ++){ *(out+i) = input[i]; }
}

//convert array to multi1d
void array2multi1d(multi1d<int> & out, int* input, int length){
    for(int i =0; i < length; i++){out[i] = input[i];}
}

//fix coordinates in the argument of a pure funtion which will be multiplied with a correlation funcion and integrated over
////eg., z-> z if z<L/2-1 and -> L - z else
void fix_coord(int out[4], const int crd[4], const int origin[4]){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < 4 ; i++){
       if(crd[i] - origin[i] > lattsize[i] / 2 ) out[i] = crd[i] - origin[i] - lattsize[i];
       else if(crd[i] - origin[i] == lattsize[i]/2) out[i] = 0;
       else out[i] = crd[i]-origin[i];
    }
}

void fix_coord(multi1d<int> & out, const multi1d<int> & crd, const multi1d<int> & origin){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < 4 ; i++){
      if(crd[i]-origin[i] > lattsize[i]/2) out[i] = crd[i] - origin[i] - lattsize[i];
      else if(crd[i] - origin[i] == lattsize[i]/2) out[i] = 0;
      else out[i] = crd[i] - origin[i];
    }
}

/*
void callInitKernelQED(){
  StopWatch timer;
  timer.reset();
  timer.start();
  QDPIO::cout<< "Call init_kernelQED" << endl; 
  init_kernelQED();
  timer.stop();
  QDPIO::cout<< "call init_kernelQED done ( " << timer.getTimeInSeconds() << " s ) "<<endl;
}
*/

/*
#ifdef FERML
void callIpihatFermLoop_antisym(double * x, double * y , double ipihat[6][4][4][4]){
  StopWatch timer;
  timer.reset();
  timer.start();
  QDPIO::cout<< "Call ipihatFermLoop_antisym" << endl;
  ipihatFermLoop_antisym(x,y, ipihat);
  timer.stop();
  QDPIO::cout<< "call ipihatFermLoop_antisym done ( " << timer.getTimeInSeconds() << " s ) "<<endl;

}
#endif
*/
/*
void callKernelQED(double* x, double * y, double kern[6][4][4][4]){
  kernelQED(x,y,kern);
}
*/

/*
//L(x,y)+L(y,x)+L(-x,y+x)
//in the entries, the coordinates are given on the lattice point, but will be transformed to the right ones which are needed when doing the integration
//output dimension: [6 rhosig][4 mu][4 nu][4 lambda]
//scalx and scaly are the scales of x and y (the kernel knows nothing about the lattice)
void kernel_3L_array(double out[6][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4]){
  int crd_x_fixed[4], crd_y_fixed[4];
  fix_coord(crd_x_fixed, crd_x, origin);
  fix_coord(crd_y_fixed, crd_y, origin);
  double ker1[6][4][4][4];
  double ker2[6][4][4][4];
  double ker3[6][4][4][4];
  double crd_x_d[4];
  double crd_y_d[4];
  for(int i =0; i < 4 ; i ++){crd_x_d[i] = crd_x_fixed[i]*scalx[i]; crd_y_d[i] = crd_y_fixed[i]*scaly[i];}
  callKernelQED(crd_x_d, crd_y_d, ker1);
  callKernelQED(crd_y_d, crd_x_d, ker2);
  double crd_minus_x[4];
  double crd_y_minus_x[4];
  for(int cpt =0; cpt <4; cpt++ ){
      crd_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt];
      crd_y_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt] + crd_y_fixed[cpt]*scaly[cpt];
  }
  callKernelQED(crd_minus_x, crd_y_minus_x, ker3);

  for(int rhosig = 0; rhosig < 6; rhosig ++){
    for(int mu = 0; mu < 4; mu++){
      for(int nu = 0; nu < 4; nu ++){
        for(int lambda = 0; lambda < 4; lambda ++){
          out[rhosig][mu][nu][lambda] = ker1[rhosig][mu][nu][lambda] + ker2[rhosig][nu][mu][lambda] + ker3[rhosig][nu][lambda][mu];
        }
      }
    }
  }    
}


//L(x,y)*x
//sum over rho => the output is of dimension [4(sigma)][4 (mu)][4 (nu) ][4 (lambda)] (the first coordinate is for sigma = {0,1,2,3})
void kernel_L_x_array(double out[4][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4] ){
  int crd_x_fixed[4], crd_y_fixed[4];
  fix_coord(crd_x_fixed, crd_x, origin);
  fix_coord(crd_y_fixed, crd_y, origin);
  double kern[6][4][4][4];
  double crd_minus_x[4];
  double crd_y_minus_x[4];
  for(int j = 0; j <4 ; j++){
      crd_minus_x[j] = -crd_x_fixed[j]*scalx[j];
      crd_y_minus_x[j] = -crd_x_fixed[j]*scalx[j] + crd_y_fixed[j]*scaly[j];
  }
  callKernelQED(crd_minus_x, crd_y_minus_x,kern);
  for(int sigma = 0; sigma < 4; sigma ++){
      for(int mu = 0; mu < 4; mu++){
          for(int nu = 0; nu < 4; nu ++){
              for(int lambda = 0; lambda < 4; lambda ++){
                  switch(sigma){
                    case 0: 
                      out[sigma][mu][nu][lambda] = -kern[0][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] -  kern[1][nu][lambda][mu]*crd_x_fixed[2]*scalx[2] - kern[2][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 1: 
                      out[sigma][mu][nu][lambda] = kern[0][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] -  kern[3][nu][lambda][mu]*crd_x_fixed[2]*scalx[2] - kern[4][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 2: 
                      out[sigma][mu][nu][lambda] = kern[1][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] +  kern[3][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] - kern[5][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 3: 
                      out[sigma][mu][nu][lambda] = kern[2][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] +  kern[4][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] + kern[5][nu][lambda][mu]*crd_x_fixed[2]*scalx[2];
                      break;
                  }
              }
          }
      }
  }    
}
*/

//imaginary part of the 3 point function (which is itself purely imaginary) with 2 local currents and one conserved
//output : function in z (fixed y) [nu][lambda][sigma]
void write_im_3_pt_llc(multi3d<LatticeRealD> &out, const multi1d<LatticeColorMatrixD> &u, const LatticePropagatorD &s_0, const LatticePropagatorD &s_y, int nu, int lambda, int sigma, const multi1d<int> & coord_y){
    
    StopWatch timer;
    timer.reset();
    timer.start();
    QDPIO::cout << "Writing the imaginary part of the three point current function (llc)" << endl;

    LatticePropagatorD s_0_inverse = Gamma(15)*adj(s_0)*Gamma(15);

    out[nu][lambda][sigma] = imag(trace( Gamma(1<<nu)*peekSite(s_0, coord_y)*Gamma(1<<lambda)*0.5*( (shift(s_0_inverse, 1, sigma)*Gamma(1<<sigma) + shift(s_0_inverse,1,sigma)) * adj(u[sigma]) * s_y + (s_0_inverse*Gamma(1<<sigma) - s_0_inverse) * u[sigma] * shift(s_y, 1, sigma)  )));
    
    timer.stop();
    QDPIO::cout << "3-pt llc current correlation function finished ( " << timer.getTimeInSeconds()<<" sec )" << endl;
}

//3pt fct for three local currents (NOT the entire expression ! For VVV, take the imag and iiiiiiiiiiy 2)
//Although only its imaginary part matters in hlbl, we would like to compare the result to the VVA correlator, which is real (it suffices to replace s_y by gamma_5 s_y and s_0 by s_0 gamma5 in the argument, take the real part and multiply by 2)
//output: function of x
void write_3_pt_lll(multi3d<LatticeComplexD> &out, const LatticePropagatorD &s_0, const LatticePropagatorD &s_y, int mu, int nu, int lambda, const multi1d<int> & coord_y){

    StopWatch timer;
    timer.reset();
    timer.start();
    QDPIO::cout << "Writing the three point current function (lll)" <<endl;
    LatticePropagatorD s_0_inverse = Gamma(15)*adj(s_0)*Gamma(15);
        
    out[mu][nu][lambda] = trace( Gamma(1<<mu)*s_y*Gamma(1<<nu)*peekSite(s_0, coord_y)*Gamma(1<<lambda)*s_0_inverse);
    
    timer.stop();
    QDPIO::cout << "3-pt lll current correlation function finished in ( " << timer.getTimeInSeconds() << " sec )" << endl;
}

//function to compute propagators S(.,0 ) and S(.,y)
//use "observer" to generate source
void compute_propagator(LatticePropagatorD &out, multi1d<int> &src_pos, int propid){

    StopWatch timer;
    timer.reset();
    timer.start();
    QDPIO::cout << "Computing propagator ..." << endl;
    
    LatticePropagatorD source;
    create_point_source(src_pos, source);
    
    QDPIO::cout << "Call observer calculate_propagator" <<endl;
    calculate_propagator(propid, source, out, 1);

    timer.stop();
    QDPIO::cout << "Finished computing propagator " << timer.getTimeInSeconds() << " sec )" <<endl; 
}

//sin function used for the checks
void write_sin_func(LatticeRealD &out, int sigma){
    LatticeRealD z_sigma = Layout::latticeCoordinate(sigma);
    int sites_on_node = Layout::sitesOnNode();
    multi1d<int> lattice_size = Layout::lattSize();
    multi1d<RealD> func(sites_on_node);
    QDP_extract(func, z_sigma, all);
    for(int i = 0; i < sites_on_node; i++){
        func[i] = sin(2*PI/lattice_size[sigma]*(toWordType(func[i]) + 0.5) );
    }
 
    QDP_insert(out, func, all);

}

//check the <llc> with Ward identity
//should use the SU(2) isospin symmetry
double check_3_pt_llc(const multi1d<LatticeColorMatrixD> & u, const LatticePropagatorD & s_0, const LatticePropagatorD &s_y, int nu, int lambda, const multi1d<int> & crd_y){

        //StopWatch timer;
        //timer.reset();
        //timer.start();
        QDPIO::cout << "Applying check function for the three point llc current for nu ="<< nu << " lambda = " << lambda  << endl;

        multi3d<LatticeRealD> prop(4,4,4);
        double res = 0;

        for(int sigma = 0; sigma <4; sigma ++){
            write_im_3_pt_llc(prop, u, s_0, s_y, nu, lambda, sigma, crd_y);

            LatticeRealD sin_func;
            write_sin_func(sin_func, sigma);

            res += toWordType(sum(sin_func*prop[nu][lambda][sigma]));
        }
 
        return res;
        //timer.stop();
        //QDPIO::cout << "Check finished ( "<< timer.getTimeInSeconds() << " sec).";
        //QDPIO::cout<< "The Ward identity is verified at the precision " << res << endl;

}

//check the <lll> by inserting a gamma5 at 0, turn off the gauge
//Compare with the result of Harvey in free theory with free Wilson fermions
//In the routine of Harvey, the pole mass is used ( am^0 = am(1-b_m am), where b_m = -1/2 ), as well as the O(a) improved vector and axial currents (where b_V = b_A = 1 at tree level).
//need to put kappa = 1/( 2(am (1-b_m am) + 4 )) in the infile
double check_3_pt_lll(const multi1d<LatticeColorMatrixD> & u, const LatticePropagatorD & s_0, const LatticePropagatorD & s_y, const multi1d<int> & crd_x, const multi1d<int> & crd_y, const multi1d<int> & origin, double am, double a){
  double out = 0;
  multi1d<int> lattsize = Layout::lattSize();
  int L_tmp[4];
  multi1d2array(L_tmp, lattsize, 4);
  double L[4];
  
  int x_tmp[4], y_tmp[4], origin_tmp[4];
  multi1d2array(x_tmp, crd_x,4);
  multi1d2array(y_tmp, crd_y, 4);
  multi1d2array(origin_tmp, origin, 4);
  double xv[4] ;
  double yv[4] ;
  double x_minus_yv[4];
  double minus_yv[4];
  for(int i =0; i < 4; i++){
    //in Harvey's code, the time direction is 0.
    //Lengths are given in lattice unit 
    int ind = (i+3)%4;
    L[i] = L_tmp[ind];
    xv[i] = (x_tmp[ind] - origin_tmp[ind]) % L_tmp[ind];
    yv[i] = (y_tmp[ind] - origin_tmp[ind]) % L_tmp[ind];
    QDPIO::cout << "L, xv, yv [" << i << "] (time first) " <<L[i] << " " << xv[i] << " " << yv[i] <<endl;
    x_minus_yv[i] = ( xv[i] - yv[i]);
    minus_yv[i] = -yv[i];
  }
 
  multi3d<LatticeComplexD> lll_qdp(4,4,4);
  for(int mu = 0 ; mu < 4; mu ++){
    for(int nu = 0; nu < 4; nu ++){
      for(int lambda = 0; lambda <4; lambda ++){
        int muh = (mu+1)%4;
        int nuh = (nu+1)%4;
        int lambdah = (lambda +1)%4;
        double lll_fwf = getavv_floopL(am, L, minus_yv, x_minus_yv, lambdah, muh, nuh);
        //double lll_fwf = getavv_floopL(am, L, y_minus_xv, minus_xv, nuh, lambdah, muh);
        //replace s_0 <- gamma5 s_0, s_y <- gamma5 s_y in the function argument
        write_3_pt_lll(lll_qdp, s_0, s_y, mu, nu, lambda, crd_y);
        //write_3_pt_lll(lll_qdp, s_0*Gamma(15), Gamma(15)*s_y, mu, nu, lambda, crd_y);
        QDPIO::cout << "Comparing the VVA 3 pt function with 3 times Harvey's result with am = " << am << " and (mu,nu,lambda) = " << mu << " " << nu << " " << lambda <<" : "<<lll_fwf*3 <<endl;
        double val = - 2*toWordType(peekSite(real(lll_qdp[mu][nu][lambda]), crd_x)); 
        double incr =  abs(val - 3*lll_fwf ); 
        QDPIO::cout << "Computed by my code : " << val << endl;       
        QDPIO::cout << "Difference : " << incr << endl;
        out += incr;
      }
    }
  }
  
  return out;
}


//Multiply by a function f(z_rho) and do the integration
//while testing, choose f = sin
//when computing, choose f = id
void int_3pt_times_func(double out[4][4][4][4], multi3d<LatticeRealD> & corr, funcD f, int origin[4]){
  multi1d<LatticeRealD> flattice(4);
  multi1d<int> lattsize = Layout::lattSize();
  multi1d<int> temp(4);
  for(int i  = 0; i < lattsize[0]; i ++ ){
   for(int j = 0; j < lattsize[1]; j ++){
     for(int k = 0; k <lattsize[2]; k ++){
       for(int l = 0; l < lattsize[3]; l++){
         for(int rho = 0; rho < 4; rho ++){
           int tmp[4] = {i,j,k,l};
           int tmp_fixed[4];
           fix_coord(tmp_fixed, tmp, origin); 
           array2multi1d(temp, tmp, 4);
           for(int rho = 0; rho < 4; rho ++){
#if 0
//test mode
             //QDPIO::cout<<"Factor 2*PI/L is applied for the int_3pt_times_func test"<<endl;
             pokeSite(flattice[rho], (Real) f(2*PI*tmp_fixed[rho]/lattsize[rho]), temp);
#else
             //QDPIO::cout<<"Normal mode for int_3pt_times_func"<<endl;
             pokeSite(flattice[rho], (Real) f(tmp_fixed[rho]), temp);
#endif
           }
          }
        }
      }
    }
  }
  for(int rho = 0; rho < 4; rho ++){
    for(int nu = 0 ; nu < 4; nu ++){
      for(int lambda = 0 ; lambda < 4; lambda ++){
        for(int sigma = 0 ; sigma < 4 ; sigma ++){
          out[rho][nu][lambda][sigma] = toWordType(sum(flattice[rho]*corr[nu][lambda][sigma]));
        }
      }
    }
  }
}

double id_func(double a){
  double out = (double) a;
  return out;
}

double sin_func(double a){
  double out = sin((double) a);
  return out;
}


/*
//the function shoulde be called by root
//kern[6][4][4][4] is needed just because of the Bcasting
//Here, crd_y should be the vector_y which is relative to the origin
void kernel_L(multi4d<LatticeRealD> &out, const int origin[4], const int crd_y[4], double a){
  double kern[6][4][4][4];

  multi1d<int> lattsize = Layout::lattSize();
  int crd_y_fixed[4];
  fix_coord(crd_y_fixed, crd_y, origin);
  multi1d<int> coord_x(4);

  QDPIO::cout<< "Constructing kernel_L" <<endl;

  RealD tmp;
  for(int i  = 0; i < lattsize[0]; i ++ ){
    for(int j = 0; j < lattsize[1]; j ++){
      for(int k = 0; k <lattsize[2]; k ++){
        for(int l = 0; l < lattsize[3]; l++){
          int crd_x[4] = {i,j,k,l};
          int crd_x_fixed[4];
          fix_coord(crd_x_fixed, crd_x, origin);
          double crd_x_d[4], crd_y_d[4];
          for(int cc =0; cc < 4 ; cc++) {crd_x_d[cc] = a*crd_x_fixed[cc]; crd_y_d[cc] = a*crd_y_fixed[cc];}
          //QDPIO::cout<< "crd_x_d "<<crd_x_d[0] << " " << crd_x_d[1] << " "<< crd_x_d[2] << " " << crd_x_d[3] <<endl;
          //QDPIO::cout<< "crd_y_d "<< crd_y_d[0]<<" " << crd_y_d[1]<< " " << crd_y_d[2] << " " <<crd_y_d[3] << endl;
          if(Layout::nodeNumber() == 0) {kernelQED(crd_x_d, crd_y_d, kern);}
          //QDPIO::cout<< "Bcasting kern for (i,j,k,l) = " <<i << " " << j << " " << k << " " << l <<endl;
          MPI_Bcast(&(kern[0][0][0][0]), 384, MPI_DOUBLE, 0 , MPI_COMM_WORLD);
          for(int rhosig = 0; rhosig < 6; rhosig ++){
            for(int mu = 0; mu < 4; mu++){
              for(int nu = 0; nu < 4; nu ++){
                for(int lambda = 0; lambda < 4; lambda ++){
                  tmp = (Real) (kern[rhosig][mu][nu][lambda]);
                  array2multi1d(coord_x, crd_x, 4);
                  pokeSite(out[rhosig][mu][nu][lambda], tmp, coord_x);
                }
              }
            }
          }
        }
      }
    }
  }
}
*/

/*
//integrate the product of L*3pt(xy0) over x
//kern_aux[6][4][4][4] is needed just because of the Bcasting
void int_L_3pt_loop(double out[6], const multi3d<LatticeRealD> & corr_3pt, const multi1d<int> & origin, const multi1d<int> & coord_y, double a){
  QDPIO::cout << "Computing L_3pt_loop"<<endl;
  int crd_y[4];
  int crd_o[4];
  multi1d2array(crd_o, origin, 4);
  multi1d2array(crd_y, coord_y, 4);
  multi4d<LatticeRealD> kern(6,4,4,4);
  kernel_L(kern, crd_o , crd_y, a);
  for(int rhosig= 0 ; rhosig< 6 ; rhosig ++ ){
    for(int mu = 0; mu < 4; mu++){
      for(int nu = 0; nu<4; nu ++){
        for(int lambda = 0; lambda < 4; lambda ++){
          out[rhosig] = toWordType(sum(kern[rhosig][mu][nu][lambda]*corr_3pt[mu][nu][lambda]) );
        }
      }
    }
  }
}
*/

//peekSite a size0*size1*size2*size3 multi4d<LatticeRealD> object 
void multi4dlattice2multi4d(multi4d<ComplexD> & out, const multi4d<LatticeRealD> & in, int size0, int size1, int size2, int size3, const multi1d<int> & coord){
  for(int i = 0; i < size0; i ++)
    for(int j = 0 ; j < size1; j ++)
      for(int k =0; k < size2; k ++)
        for(int l =0; l < size3; l ++)
          out[i][j][k][l] = (ComplexD) (peekSite(in[i][j][k][l], coord));
}

void array2multi4dlattice(multi4d<LatticeRealD> & out, double in[6][4][4][4], int size0, int size1, int size2, int size3, const multi1d<int> & coord){
  for(int i = 0; i < size0; i ++)
    for(int j = 0 ; j < size1; j ++)
      for(int k =0; k < size2; k ++)
        for(int l =0; l < size3; l ++)
          pokeSite(out[i][j][k][l], (RealD)in[i][j][k][l], coord);
}

//use Jamie's code
void kernel_3L_new(double out[6][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4], struct QED_kernel_temps *kert){
  int crd_x_fixed[4], crd_y_fixed[4];
  fix_coord(crd_x_fixed, crd_x, origin);
  fix_coord(crd_y_fixed, crd_y, origin);
  double ker1[6][4][4][4];
  double ker2[6][4][4][4];
  double ker3[6][4][4][4];
  double crd_x_d[4];
  double crd_y_d[4];
  for(int i =0; i < 4 ; i ++){crd_x_d[i] = crd_x_fixed[i]*scalx[i]; crd_y_d[i] = crd_y_fixed[i]*scaly[i];}
  QED_kernel_L0(crd_x_d, crd_y_d, kert, ker1);
  QED_kernel_L0(crd_y_d, crd_x_d, kert, ker2);
  double crd_minus_x[4];
  double crd_y_minus_x[4];
  for(int cpt =0; cpt <4; cpt++ ){
      crd_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt];
      crd_y_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt] + crd_y_fixed[cpt]*scaly[cpt];
  }
  QED_kernel_L0(crd_minus_x, crd_y_minus_x, kert, ker3);

  for(int rhosig = 0; rhosig < 6; rhosig ++){
    for(int mu = 0; mu < 4; mu++){
      for(int nu = 0; nu < 4; nu ++){
        for(int lambda = 0; lambda < 4; lambda ++){
          out[rhosig][mu][nu][lambda] = ker1[rhosig][mu][nu][lambda] + ker2[rhosig][nu][mu][lambda] + ker3[rhosig][nu][lambda][mu];
        }
      }
    }
  }
}

//use Jamie's code
void kernel_L_x_new(double out[4][4][4][4], const int origin[4], const int crd_x[4], const int crd_y[4], const double scalx[4], const double scaly[4], struct QED_kernel_temps * kert ){
  int crd_x_fixed[4], crd_y_fixed[4];
  fix_coord(crd_x_fixed, crd_x, origin);
  fix_coord(crd_y_fixed, crd_y, origin);
  double kern[6][4][4][4];
  double crd_minus_x[4];
  double crd_y_minus_x[4];
  for(int j = 0; j <4 ; j++){
      crd_minus_x[j] = -crd_x_fixed[j]*scalx[j];
      crd_y_minus_x[j] = -crd_x_fixed[j]*scalx[j] + crd_y_fixed[j]*scaly[j];
  }
  QED_kernel_L0(crd_minus_x, crd_y_minus_x,kert, kern);
  for(int sigma = 0; sigma < 4; sigma ++){
      for(int mu = 0; mu < 4; mu++){
          for(int nu = 0; nu < 4; nu ++){
              for(int lambda = 0; lambda < 4; lambda ++){
                  switch(sigma){
                    case 0:
                      out[sigma][mu][nu][lambda] = -kern[0][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] -  kern[1][nu][lambda][mu]*crd_x_fixed[2]*scalx[2] - kern[2][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 1:
                      out[sigma][mu][nu][lambda] = kern[0][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] -  kern[3][nu][lambda][mu]*crd_x_fixed[2]*scalx[2] - kern[4][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 2:
                      out[sigma][mu][nu][lambda] = kern[1][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] +  kern[3][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] - kern[5][nu][lambda][mu]*crd_x_fixed[3]*scalx[3];
                      break;
                    case 3:
                      out[sigma][mu][nu][lambda] = kern[2][nu][lambda][mu]*crd_x_fixed[0]*scalx[0] +  kern[4][nu][lambda][mu]*crd_x_fixed[1]*scalx[1] + kern[5][nu][lambda][mu]*crd_x_fixed[2]*scalx[2];
                      break;
                  }
              }
          }
      }
  }
}

