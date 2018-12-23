#ifndef GETL_H
#define GETL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NF 64 
#define N 4400
#define NMAX 100 
// 64
/*
// scan up to 2.18 used:
#define nfx(x) ( 8  + (int) (5.0*x))
// #define nfx(x) ( (int) (5.0*x)>10 ? 18 : 8 + (int) (5.0*x))

// dy should be set to y1-y2
#define interpol3(y,y1,y2,dy,f1,f2,g1,g2)     ((-f1* (y - y2)*(y - y2)* (2* y - 3* y1 + y2) + \
       (y - y1)* (f2* (y - y1)* (2* y + y1 - 3* y2) \
     + (y - y2)*dy* (g1* y + g2* y - g2* y1 - g1* y2)))/(dy*dy*dy))

static const int Idm[24]={0,0,1,0,1,2, -1, -1, -1, -1, 0,0,1,0,1,2, -1, -1, -1, -1, -1, -1, 1,0}; 
static int Iconv[24]={0,1,2,3,4,5, -1,-1,-1,-1, 6, 7, 8, 9, 10, 11, -1, -1, -1, -1, -1, -1, 12, 13};

// double (*Func_sig[30])(int,double,double,double);
// static void (*Func_uch[3])(int,double,double *);
double (*Func_usm[4])(int,double,double *);

// #define READ_ON_THE_FLY

#define XMYSWAP_S
#define XMYSWAP_V
#define XMYSWAP_T

#define SYMG

#define TAYLORX
#define TAYLORY

#ifdef TAYLORX
#define NY_TAY 810
void Tabd_xeq0(double yv[4], double ***tI, double ***tII, double ***tIII);
double Stp_tay, YY_tay[NY_TAY];
double G0dx_tay[NY_TAY], G0dy_tay[NY_TAY];
double Gl2_tay[NY_TAY], Gl2dy_tay[NY_TAY], Gl21_tay[NY_TAY], Gl21dy_tay[NY_TAY], Gl3_tay[NY_TAY], Gl3dy_tay[NY_TAY];
double G21_tay[NY_TAY], G21dy_tay[NY_TAY], G22A_tay[NY_TAY], G22B_tay[NY_TAY], G22Ady_tay[NY_TAY], G22Bdy_tay[NY_TAY];
double  G3A_tay[NY_TAY], G3B_tay[NY_TAY], G3Ady_tay[NY_TAY], G3Bdy_tay[NY_TAY], G31A_tay[NY_TAY], G31B_tay[NY_TAY], G31Ady_tay[NY_TAY], G31Bdy_tay[NY_TAY];
#endif

#ifdef TAYLORY
void Tabd_yeq0(double xv[4], double ***tI, double ***tII, double ***tIII);
#define NX_TAY 100
double alpha0dx_0p_taY[NX_TAY], alpha0_1p_taY[NX_TAY], alpha3_0p_taY[NX_TAY], beta2_1p_taY[NX_TAY], alpha3_1p_taY[NX_TAY];
double  alpha1_0p_taY[NX_TAY], beta4_1p_taY[NX_TAY], alpha1_1p_taY[NX_TAY], alpha1dx_0p_taY[NX_TAY];
double alpha1dx_1p_taY[NX_TAY], alpha3dx_0p_taY[NX_TAY], alpha3dxdx_0p_taY[NX_TAY], beta2dx_1p_taY[NX_TAY], alpha3dx_1p_taY[NX_TAY];
#endif

#define NFFA 14 

double ***Ffp[NFFA], ***Ffm[NFFA];

// #define TESTING
double testff(int is_default, double x, double cb, double y, double *lambdas, double *v);
double testff2(int is_default, double x, double cb, double y, double *lambdas, double *v);

void getff(int nff, int nf, int nm, int mm, int ndy, int ndcb, 
	   double x, double y, double *cbv, double *fm, double *fp, double *ff);

double chebUsum(int nk, double x, double *co);
double dchebUsum(int nk, double x, double *co);
double ddchebUsum(int nk, double x, double *co);
double dddchebUsum(int nk, double x, double *co);

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint(double *xa, double *ya, double *y2a, int n, double x, double *y);
void gaussj(double **a, int n, double **b, int m);
int computeG();

void chnr_dS(double xv[4], double yv[4], double *dxv, double *dyv);
void chnr_dT(double xv[4], double yv[4], double ***dxv, double ***dyv);
void chnr_dV(double xv[4], double yv[4], double ***dv);

void init_ff(double xh, double yh, int nstpx, int nstpy);
double extractff(int nm, int ndy, int ndcb, double x, double cb, double y);
void readf(FILE *fp, int ny, double *xr, double *yr, double *am, double *ap, int verbose);
double accessf(FILE *fp, int flag_hy, int use_y_derivs, int ny, int nm, 
	       int ndy, int ndcb, double *x1, double *y1, double *y2, double x, double *cbv, double y);
double accessv(int flag_hy, int use_y_derivs, int ix, int iy, int nm, 
	       int ndy, int ndcb, double *x1, double *y1, double *y2, double x, double *cbv, double y);

void gI(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4]);
int getgI(int rho, int sig, int mu, int nu, int lda, int alf, int bet, int dta);
void gIII(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4]);
void gII(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4]);
int getgII(int rho, int sig, int mu, int nu, int lda, int alf, int bet, int dta);
int diractrace(int n, int *idx);


void dKpiB(double x, double cb, double y, double kpi[4][4][4], double kpiB[4][4][4]);

//copy from amu_forlattice.c

//#define FERML
*/
//void init_kernelQED(void);
//int kernelQED(double *xv, double *yv, double kerv[6][4][4][4]);
//int kernelQED_xoryeq0(int flagxy, double *qv, double kerv[6][4][4][4]);
//void init_G();
/*
static  double sxv[4], syv[4], ***txv, ***tyv, ***vv;
static  double sxv_xswapy[4], syv_xswapy[4], ***txv_xswapy, ***tyv_xswapy, ***vv_xswapy;
static  double sxv_yswapxmy[4], syv_yswapxmy[4], ***txv_yswapxmy, ***tyv_yswapxmy, ***vv_yswapxmy;
static  double sxv_yswapxmy_xswapy[4], syv_yswapxmy_xswapy[4], ***txv_yswapxmy_xswapy, ***tyv_yswapxmy_xswapy, ***vv_yswapxmy_xswapy;
static  double sxv5[4], syv5[4], ***txv5, ***tyv5, ***vv5;
static  double sxv6[4], syv6[4], ***txv6, ***tyv6, ***vv6;
static int *gIv[4][4][6][4], *gIIv[4][4][6][4], *gIIIv[4][4][6][4];
static int *alfIv[4][4][6][4], *alfIIv[4][4][6][4], *alfIIIv[4][4][6][4], *betIv[4][4][6][4], *betIIv[4][4][6][4], *betIIIv[4][4][6][4],
    *dtaIv[4][4][6][4], *dtaIIv[4][4][6][4], *dtaIIIv[4][4][6][4], nIv[4][4][6][4], nIIv[4][4][6][4], nIIIv[4][4][6][4];

//#define SYMXY
// #define SYMXY0

#define v4norm(x) sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
#define kro(mu,nu) (mu==nu ? 1 : 0)

#ifdef FERML
static int g8[4][4][4][4][4][4][4][4];

void PihatFermionLoop(double xv[4], double yv[4], double pihat[6][4][4][4]);
void pihat1(double x[4], double y[4], double vpihat1[4][4][4][4][4], double vpihatr1[4][4][4][4]);
void ipihatFermLoop_antisym(double x[4], double y[4], double vpihat[6][4][4][4]);
void init_g8();
double amu_heavlept(double mh);
double amu_lightlept(double ml);
void bessk(int n, double x, double *bv);
#endif

#define SCAL_FAC 1.0

#ifdef AMU_FORLATTICE
double Mv = 0.5; // ratio of lepton-loop mass to muon mass 
double AlfQED = 1.0/137.035999;
#endif
*/
//here are functions used by the lll check
//#ifdef VVA
//void bessk(int n, double x, double *bv);  //declared before
void s0prop(double ma, double x, double a[2]);
void s0propL(double ma, double xv[4], double Ll[4], double a[5]);
double getavv_floopL(double ma, double Ll[4], double x[4], double y[4], int mu, int nu, int rho);
void set_epstensor();

void momentum_sum3d(int ns, void (*func)(double,double,double,double *,double *), double *p, double *sum,
            int L1, int L2, int L3, double theta1, double theta2, double theta3);
void thermprop(double px, double py, double pz, double *p, double prop[5]);
void prop_smd(double px, double py, double pz, double *p, double qr[5]);
static int Epst[4][4][4][4];
//#endif


#endif
