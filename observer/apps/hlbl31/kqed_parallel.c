extern "C"{
#include "KQED.h"
}

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
  QED_kernel_L3(crd_x_d, crd_y_d, kert, ker1);
  QED_kernel_L3(crd_y_d, crd_x_d, kert, ker2);
  double crd_minus_x[4];
  double crd_y_minus_x[4];
  for(int cpt =0; cpt <4; cpt++ ){
      crd_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt];
      crd_y_minus_x[cpt] = -crd_x_fixed[cpt]*scalx[cpt] + crd_y_fixed[cpt]*scaly[cpt];
  }
  QED_kernel_L3(crd_minus_x, crd_y_minus_x, kert, ker3);

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
  QED_kernel_L3(crd_minus_x, crd_y_minus_x,kern);
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


