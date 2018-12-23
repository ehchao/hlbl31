
// gcc -o amu_forlattice  amu_forlattice.c  getff-new.c   traces.c  -lm

#include "getL.h"
// #include "../../CUBATURE/cubature-1.0.2/cubature.h"

#define PI 3.14159265358979323846
#define PI5 (0.5*PI)
#define pi PI
#define pi2 (PI*2.0)
#define PI8 25.1327412287183459077011
#define ZETA3 1.2020569031595942854

#define FERML

int kernelQED(double *xv, double *yv, double kerv[6][4][4][4]);
int kernelQED_xoryeq0(int flagxy, double *qv, double kerv[6][4][4][4]);
void init_G();

static  double sxv[4], syv[4], ***txv, ***tyv, ***vv;
static  double sxv_xswapy[4], syv_xswapy[4], ***txv_xswapy, ***tyv_xswapy, ***vv_xswapy;
static  double sxv_yswapxmy[4], syv_yswapxmy[4], ***txv_yswapxmy, ***tyv_yswapxmy, ***vv_yswapxmy;
static  double sxv_yswapxmy_xswapy[4], syv_yswapxmy_xswapy[4], ***txv_yswapxmy_xswapy, ***tyv_yswapxmy_xswapy, ***vv_yswapxmy_xswapy;
static  double sxv5[4], syv5[4], ***txv5, ***tyv5, ***vv5;
static  double sxv6[4], syv6[4], ***txv6, ***tyv6, ***vv6;
static int *gIv[4][4][6][4], *gIIv[4][4][6][4], *gIIIv[4][4][6][4];
static int *alfIv[4][4][6][4], *alfIIv[4][4][6][4], *alfIIIv[4][4][6][4], *betIv[4][4][6][4], *betIIv[4][4][6][4], *betIIIv[4][4][6][4],
    *dtaIv[4][4][6][4], *dtaIIv[4][4][6][4], *dtaIIIv[4][4][6][4], nIv[4][4][6][4], nIIv[4][4][6][4], nIIIv[4][4][6][4];

// #define SYMXY
// #define SYMXY0

#ifdef FERML
#define v4norm(x) sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
#define kro(mu,nu) (mu==nu ? 1 : 0)

static int g8[4][4][4][4][4][4][4][4];

void PihatFermionLoop(double xv[4], double yv[4], double pihat[6][4][4][4]);
void pihat1(double x[4], double y[4], double vpihat1[4][4][4][4][4], double vpihatr1[4][4][4][4]);
void ipihatFermLoop_antisym(double x[4], double y[4], double vpihat[6][4][4][4]);
void init_g8();
double amu_heavlept(double mh);
double amu_lightlept(double ml);
void bessk(int n, double x, double *bv);

#define SCAL_FAC 1.0
double Mv = 0.5; // ratio of lepton-loop mass to muon mass 
double AlfQED = 1.0/137.035999;
#endif 

 
//int main(int argc, char *argv[])
void init_kernelQED()
{
  int k, n, i, count, nf, nff, j, jmax, nm, nn, ndcb, l, l1, inv_order, isuc, idx, ndint, ell;
  int mu, nu, lda, rhosig, rho, sig, ix, iy, is_x, is_y;
  double err1, err2, symfact, r, beta, tmp, kerv[6][4][4][4], kerv2[6][4][4][4];
  double u, v, phi, f, g, h0[N], h1[N], zv[N], vg3bs[N], vg3bi[N], xi1, xi2, xi3, xi4, xi5, stp, bstp, sti, stf, amu, sum, ya, yb, yc;
  double x, y, b, z, zup, zdn, val[10],  xmin[10], xmax[10], kpi[4][4][4], kpib[4][4][4];
  double  fnorm, cbv[N], *kptr[4], xymin, xymax, stpx;
  double cf, sf, sf2, cb, sb, xbuf;
  double integrmax, pihat[6][4][4][4], pihatp[6][4][4][4], pihatm[6][4][4][4];
  double reqAbsError, reqRelError;
  double xh, yh, fpi, ker, xsq, arg2, vmin, vmax, hv[4][4][4];
  double yv[4], xv[4], xsv[4], alat, xvMv[4], yvMv[4], xmyv[4], xmysq, xmy, fdat[10], fval[10], err[10], vi[4], co, si;
  double pref, convf, norm_pihatpi0, two_rhosig, cpi0,   norm_fermloop;
  int nstpx, nstpy, ndy, use_y_derivs, nxv[4], nvol, i0, i1, i2, i3, alf, bet, dta;
  unsigned norm, maxEval, dim, fdim;
  FILE *fp;

  //  printf("AlfQED= %lg, Mv= %lg amu_heavlept = %.8lg\n", AlfQED, Mv, amu_heavlept(Mv)); return 1;
  //  printf("AlfQED= %.11lg, 1.0/Mv= %.11lg amu_lightlept = %.11lg\n", AlfQED, 1.0/Mv, amu_lightlept(Mv)); return 1;

#ifdef FERML
  init_g8();
#endif 

  fp=fopen("in/input", "r");
  fscanf(fp, "%d %lg %lg %lg %d %d %d %d %d %lg %lg", &nm, &x, cbv+0, &y, &ndy, &ndcb, &k, &nstpx, &nstpy, &xh, &yh);
  use_y_derivs = k && (ndy==0);
  printf("# nm= %d    x= %lg   cbv[0]= %lg    y= %lg   ndy= %d    ndcb= %d   use_y_derivs= %d\n", 
	 nm, x, cbv[0], y, ndy, ndcb, use_y_derivs);
  printf("# nstpx= %d  nstpy= %d   xh= %lg   yh= %lg\n", nstpx, nstpy, xh, yh);
  fclose(fp);
  init_ff(xh, yh, nstpx, nstpy);
  printf("# Initialized QED kernel FFs...\n");

  txv = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv[i][j] = (double *) malloc(4*sizeof(double));
  tyv = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv[i][j] = (double *) malloc(4*sizeof(double));
  vv = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv[i][j] = (double *) malloc(4*sizeof(double));

  txv_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv_xswapy[i][j] = (double *) malloc(4*sizeof(double));
  tyv_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv_xswapy[i][j] = (double *) malloc(4*sizeof(double));
  vv_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv_xswapy[i][j] = (double *) malloc(4*sizeof(double));

  txv_yswapxmy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv_yswapxmy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv_yswapxmy[i][j] = (double *) malloc(4*sizeof(double));
  tyv_yswapxmy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv_yswapxmy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv_yswapxmy[i][j] = (double *) malloc(4*sizeof(double));
  vv_yswapxmy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv_yswapxmy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv_yswapxmy[i][j] = (double *) malloc(4*sizeof(double));

  txv_yswapxmy_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv_yswapxmy_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv_yswapxmy_xswapy[i][j] = (double *) malloc(4*sizeof(double));
  tyv_yswapxmy_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv_yswapxmy_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv_yswapxmy_xswapy[i][j] = (double *) malloc(4*sizeof(double));
  vv_yswapxmy_xswapy = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv_yswapxmy_xswapy[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv_yswapxmy_xswapy[i][j] = (double *) malloc(4*sizeof(double));

  txv5 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv5[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv5[i][j] = (double *) malloc(4*sizeof(double));
  tyv5 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv5[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv5[i][j] = (double *) malloc(4*sizeof(double));
  vv5 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv5[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv5[i][j] = (double *) malloc(4*sizeof(double));

  txv6 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) txv6[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) txv6[i][j] = (double *) malloc(4*sizeof(double));
  tyv6 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) tyv6[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) tyv6[i][j] = (double *) malloc(4*sizeof(double));
  vv6 = (double ***) malloc(4*sizeof(double **));
  for(i=0;i<4;i++) vv6[i] = (double **) malloc(4*sizeof(double *));
  for(i=0;i<4;i++) for(j=0;j<4;j++) vv6[i][j] = (double *) malloc(4*sizeof(double));

  for(i=0;i<4;i++)  for(j=0;j<4;j++) for(k=0;k<6;k++) for(ell=0;ell<4;ell++) {
      gIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      alfIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      betIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      dtaIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      gIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      alfIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      betIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      dtaIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      gIIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      alfIIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      betIIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
      dtaIIIv[i][j][k][ell] = (int *) malloc(NMAX*sizeof(int));
	}
  printf("# main: Allocated memory for QED kernel..."); fflush(stdout);
  init_G();
  printf("done\n");
}
/*  convf = pow(Mv,7);
  pref = (pow(4.0*pi*AlfQED,3)/3)*(2.0*pi*pi)*(4*pi);
  two_rhosig = 2.0;
#ifdef FERML
  norm_fermloop = pow(1.0,4);
  fnorm = two_rhosig*pref*(convf*norm_fermloop);
#endif

#ifdef SYMG
  symfact = 0.5;
#else 
  symfact = 1.0;
#endif 

  vi[0] = 1.228; // 0.0242223; // 0.0021; // x
  vi[1] = 0.209243; // beta
  vi[2] = 0.0484223; // 0.218; // y
  printf("# vi[0]= %lg   vi[1] = %lg   vi[2] = %lg\n", vi[0], vi[1], vi[2]);
  co = cos(vi[1]);
  si = sin(vi[1]);
  yv[0] = yv[1] = yv[2] = 0.0;  yv[3] = vi[2];
  xv[0] = xv[1] = 0.0;  xv[2] = vi[0]*si;  xv[3] = vi[0]*co;
  xvMv[0] = xv[0]*Mv;  xvMv[1] = xv[1]*Mv;  xvMv[2] = xv[2]*Mv;  xvMv[3] = xv[3]*Mv;
  yvMv[0] = yv[0]*Mv;  yvMv[1] = yv[1]*Mv;  yvMv[2] = yv[2]*Mv;  yvMv[3] = yv[3]*Mv;
  i = kernelQED(xv, yv, kerv);
#if 0
  ipihatFermLoop_antisym(xvMv, yvMv, pihat);
  tmp = 0.0;
  for(i=0;i<6;i++) for(j=0;j<4;j++) for(k=0;k<4;k++) for(l=0;l<4;l++) {
	  tmp += pihat[i][j][k][l]*kerv[i][j][k][l];
	}
  tmp *= fnorm;
  printf("x= %lf beta= %lf  y= %lf  iPihat*L_QED= %.11lg\n", vi[0], vi[1], vi[2], pow(vi[0],3.0)*pow(vi[2],4.0)*si*si*tmp);
#endif
#if 1
  is_x = 1;
  is_y = 0;
  //  i = kernelQED_xoryeq0(is_x, yv, kerv2); 
  i = kernelQED_xoryeq0(is_y, xv, kerv2); 
  //  return 1;
  for(i=0;i<6;i++) for(j=0;j<4;j++) for(k=0;k<4;k++) for(l=0;l<4;l++) {
	  printf("rhosig= %d\tmu= %d\tnu= %d\tlda= %d\t%lg   %lg\n", i, j, k, l, symfact*kerv[i][j][k][l], symfact*kerv2[i][j][k][l]);
	}
#endif
  return 1;


}
*/


#if 0
#define SINGLE_CONTRIB
// only vector contribution:
#define CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv)\
	    ker=0.0;\
          for(i=0;i<nIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIv[mu][nu][rhosig][lda][i];\
	      bet = betIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIv[mu][nu][rhosig][lda][i];\
	      g = gIv[mu][nu][rhosig][lda][i];\
	      ker += g*vv[alf][bet][dta];	\
	    }
#endif
#if 0
#define SINGLE_CONTRIB
// only tensor contribution:
#define CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv)\
              ker=0.0;\
	    for(i=0;i<nIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIv[mu][nu][rhosig][lda][i];\
	      g = gIIv[mu][nu][rhosig][lda][i];\
	      ker += g*txv[bet][alf][dta];\
	    }\
	    for(i=0;i<nIIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIIv[mu][nu][rhosig][lda][i];\
	      g = gIIIv[mu][nu][rhosig][lda][i];\
	      ker += g*(txv[alf][bet][dta]+tyv[alf][bet][dta]);\
	    }
#endif
#if 0
#define SINGLE_CONTRIB
// only scalar contribution:
#define CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv)\
	    for(i=0;i<nIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIv[mu][nu][rhosig][lda][i];\
	      g = gIIv[mu][nu][rhosig][lda][i];\
	      if(bet==dta) ker += 0.25*g*sxv[alf];	\
	    }\
	    for(i=0;i<nIIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIIv[mu][nu][rhosig][lda][i];\
	      g = gIIIv[mu][nu][rhosig][lda][i];\
	      if(alf==dta) ker += 0.25*g*(sxv[bet]+syv[bet]);\
	    }
#endif

#if 0
#define SINGLE_CONTRIB
#define CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv)\
	    ker=0.0;\
          for(i=0;i<nIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIv[mu][nu][rhosig][lda][i];\
	      bet = betIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIv[mu][nu][rhosig][lda][i];\
	      g = gIv[mu][nu][rhosig][lda][i];\
	      ker += g*vv[alf][bet][dta];	\
	    }
#endif 

#ifndef SINGLE_CONTRIB
// normal case: all contributions included
#define CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv)\
	    ker=0.0;\
          for(i=0;i<nIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIv[mu][nu][rhosig][lda][i];\
	      bet = betIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIv[mu][nu][rhosig][lda][i];\
	      g = gIv[mu][nu][rhosig][lda][i];\
	      ker += g*vv[alf][bet][dta];	\
	    }\
	    for(i=0;i<nIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIv[mu][nu][rhosig][lda][i];\
	      g = gIIv[mu][nu][rhosig][lda][i];\
	      ker += g*txv[bet][alf][dta];\
	      if(bet==dta) ker += 0.25*g*sxv[alf];\
	    }\
	    for(i=0;i<nIIIv[mu][nu][rhosig][lda];i++) {\
	      alf = alfIIIv[mu][nu][rhosig][lda][i];\
	      bet = betIIIv[mu][nu][rhosig][lda][i];\
	      dta = dtaIIIv[mu][nu][rhosig][lda][i];\
	      g = gIIIv[mu][nu][rhosig][lda][i];\
	      ker += g*(txv[alf][bet][dta]+tyv[alf][bet][dta]);\
	      if(alf==dta) ker += 0.25*g*(sxv[bet]+syv[bet]);\
	    }
#endif



int kernelQED(double *xv, double *yv, double kerv[6][4][4][4])
{
  double xmyv[4], ymxv[4], co, si, g, xmysq;
  double ker, x, y, xsq, ysq, xovy;
  int mu, nu, lda, rhosig, rho, sig, alf, bet, dta, i;

  for(rhosig=0;rhosig<6;rhosig++) for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) for(lda=0;lda<4;lda++) kerv[rhosig][mu][nu][lda]=0.0;

  //  printf("xv= %lf  %lf  %lf  %lf \n", xv[0], xv[1], xv[2], xv[3]);
  //  printf("yv= %lf  %lf  %lf  %lf \n", yv[0], yv[1], yv[2], yv[3]); 

  // FFs of QED kernel and their derivatives wrt x,cb,y
  chnr_dS(xv, yv, sxv, syv);
  chnr_dT(xv, yv, txv, tyv);
  chnr_dV(xv, yv, vv);
#ifdef SYMXY
  chnr_dS(yv, xv, sxv_xswapy, syv_xswapy);
  chnr_dT(yv, xv, txv_xswapy, tyv_xswapy);
  chnr_dV(yv, xv, vv_xswapy);
#ifdef SYMXY0
  xmyv[0] = xv[0]-yv[0]; xmyv[1] = xv[1]-yv[1];  xmyv[2] = xv[2]-yv[2];  xmyv[3] = xv[3]-yv[3];
  ymxv[0] = -xmyv[0];   ymxv[1] = -xmyv[1];   ymxv[2] = -xmyv[2];   ymxv[3] = -xmyv[3]; 

  chnr_dS(xv, xmyv, sxv_yswapxmy, syv_yswapxmy);
  chnr_dT(xv, xmyv, txv_yswapxmy, tyv_yswapxmy);
  chnr_dV(xv, xmyv, vv_yswapxmy);

  chnr_dS(yv, ymxv, sxv_yswapxmy_xswapy, syv_yswapxmy_xswapy);
  chnr_dT(yv, ymxv, txv_yswapxmy_xswapy, tyv_yswapxmy_xswapy);
  chnr_dV(yv, ymxv, vv_yswapxmy_xswapy);

  chnr_dS(xmyv, xv, sxv5, syv5);
  chnr_dT(xmyv, xv, txv5, tyv5);
  chnr_dV(xmyv, xv, vv5);

  chnr_dS(ymxv, yv, sxv6, syv6);
  chnr_dT(ymxv, yv, txv6, tyv6);
  chnr_dV(ymxv, yv, vv6);
#endif
#endif

  for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) { 
      //      gI(mu,nu,nIv,gIv,alfIv,betIv,dtaIv);
      //      gII(mu,nu,nIIv,gIIv,alfIIv,betIIv,dtaIIv);
      //      gIII(mu,nu,nIIIv,gIIIv,alfIIIv,betIIIv,dtaIIIv); 
      
      rhosig=0;
      for(rho=0;rho<3;rho++) for(sig=rho+1;sig<4;sig++) {
	  for(lda=0;lda<4;lda++) {

	    CONSTRUCT_KERNEL(sxv, syv, vv, txv, tyv);	    
	    kerv[rhosig][mu][nu][lda] += ker;
#ifdef SYMXY
	    CONSTRUCT_KERNEL(sxv_xswapy, syv_xswapy, vv_xswapy, txv_xswapy, tyv_xswapy);
	    kerv[rhosig][nu][mu][lda] += ker;
#ifdef SYMXY0
	    CONSTRUCT_KERNEL(sxv_yswapxmy, syv_yswapxmy, vv_yswapxmy, txv_yswapxmy, tyv_yswapxmy);
	    kerv[rhosig][lda][nu][mu] -= ker;

	    CONSTRUCT_KERNEL(sxv_yswapxmy_xswapy, syv_yswapxmy_xswapy, vv_yswapxmy_xswapy, txv_yswapxmy_xswapy, tyv_yswapxmy_xswapy); 
	    kerv[rhosig][nu][lda][mu] -= ker;

	    CONSTRUCT_KERNEL(sxv5, syv5, vv5, txv5, tyv5);
	    kerv[rhosig][lda][mu][nu] -= ker;

	    CONSTRUCT_KERNEL(sxv6, syv6, vv6, txv6, tyv6); 
	    kerv[rhosig][mu][lda][nu] -= ker;
#endif
#endif	
	  }
	  ++rhosig;
	}
    }

  return(0);
}	  


// flagxy = 1:  kernel at x=0; 2nd argument is vector y
// flagxy = 0:  kernel at y=0; 2nd argument is vector x
int kernelQED_xoryeq0(int flagxy, double *qv, double kerv[6][4][4][4])
{
  double xmyv[4], ymxv[4], co, si, g, xmysq;
  double ker, x, y, xsq, ysq, xovy;
  int mu, nu, lda, rhosig, rho, sig, alf, bet, dta, i;

  for(rhosig=0;rhosig<6;rhosig++) for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) for(lda=0;lda<4;lda++) kerv[rhosig][mu][nu][lda]=0.0;

  // FFs of QED kernel and their derivatives wrt x,cb,y
  //  chnr_dS(xv, yv, sxv, syv);  chnr_dT(xv, yv, txv, tyv);  chnr_dV(xv, yv, vv);
  if(flagxy==1)  Tabd_xeq0(qv, vv, txv, tyv); 
  else           Tabd_yeq0(qv, vv, txv, tyv); 
  // now: vv is T^{(I)}_{\alpha\beta\delta}(0,y)   (the vector contribution)
  //      txv is T^{(II)}_{\alpha\beta\delta}(0,y) (including the scalar contribution)
  //      tyv is T^{(III)}_{\alpha\beta\delta}(0,y) (including the scalar contribution)

  for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) { 
      //      gI(mu,nu,nIv,gIv,alfIv,betIv,dtaIv);
      //      gII(mu,nu,nIIv,gIIv,alfIIv,betIIv,dtaIIv);
      //      gIII(mu,nu,nIIIv,gIIIv,alfIIIv,betIIIv,dtaIIIv); 
      
      rhosig=0;
      for(rho=0;rho<3;rho++) for(sig=rho+1;sig<4;sig++) {
	  for(lda=0;lda<4;lda++) {	    
	    
	    ker=0.0;
	    for(i=0;i<nIv[mu][nu][rhosig][lda];i++) {
	      alf = alfIv[mu][nu][rhosig][lda][i];
	      bet = betIv[mu][nu][rhosig][lda][i];
	      dta = dtaIv[mu][nu][rhosig][lda][i];
	      g = gIv[mu][nu][rhosig][lda][i];
	      ker += g*vv[alf][bet][dta];	  
	    }
	    for(i=0;i<nIIv[mu][nu][rhosig][lda];i++) {
	      alf = alfIIv[mu][nu][rhosig][lda][i];
	      bet = betIIv[mu][nu][rhosig][lda][i];
	      dta = dtaIIv[mu][nu][rhosig][lda][i];
	      g = gIIv[mu][nu][rhosig][lda][i];
	      ker += g*txv[alf][bet][dta]; 
	    }
	    for(i=0;i<nIIIv[mu][nu][rhosig][lda];i++) {
	      alf = alfIIIv[mu][nu][rhosig][lda][i];
	      bet = betIIIv[mu][nu][rhosig][lda][i];
	      dta = dtaIIIv[mu][nu][rhosig][lda][i];
	      g = gIIIv[mu][nu][rhosig][lda][i];
	      ker += g*tyv[alf][bet][dta];
	    }
	    kerv[rhosig][mu][nu][lda] += ker;
	  }
	  ++rhosig;
	}
    }
  
  return(0);
}	  

void init_G()
{
  gI(nIv,gIv,alfIv,betIv,dtaIv);
  gII(nIIv,gIIv,alfIIv,betIIv,dtaIIv);
  gIII(nIIIv,gIIIv,alfIIIv,betIIIv,dtaIIIv);
}

#ifdef FERML

// returns contribution of the lbl contribution due to a heavy lepton; mh = (lepton mass)/(muon mass)
// including the (alpha/pi)^3 factor
// from Remiddi-Laporta PLB 301 (1993) 440 Eq. (4)
double amu_heavlept(double mh)
{
  double t2, t4, t6, t8, t10, mh2, mh4, ln, zeta2;
  zeta2 = pi*pi/6;
  mh2 = mh*mh;
  mh4 = mh2*mh2;
  ln = log(mh);
  t2= (1.5*ZETA3-19.0/16)/mh2;
  t4= (-161.0/810*ln*ln - 16189.0/48600*ln + 13.0/18*ZETA3 - 161.0/9720*pi*pi - 831931.0/972000)/(mh2*mh2);
  t6= (17.0/36*ZETA3-13.0/224*zeta2-1840256147.0/3556224000-4381.0/120960*4.0*ln*ln-24761.0/317520*2.0*ln)/(mh4*mh2);
  t8= (7.0/20*ZETA3 - 2047.0/54000*zeta2 - 453410778211.0/1200225600000 - 5207.0/189000*4.0*ln*ln - 41940853.0/952560000*2.0*ln)/(mh4*mh4);
  t10=(5.0/18*ZETA3 - 1187.0/44550*zeta2 - 86251554753071.0/287550049248000 - 328337.0/14968800*4.0*ln*ln - 640572781.0/23051952000*2.0*ln)/(mh4*mh4*mh2);
  printf("# t2=%.11lg\n# t4=%.11lg\n# t6=%.11lg\n# t8=%.11lg\n# t10=%.11lg\n# sum= %.11lg\n", t2, t4, t6, t8, t10, t2+t4+t6+t8+t10);
  return(AlfQED*AlfQED*AlfQED/(pi*pi*pi)*(t2+t4+t6+t8+t10));
  // return((t2+t4+t6+t8+t10));
}

// contribution of a light lepton in a_mu, with ml=m_lepton/m_mu
// // from Remiddi-Laporta PLB 301 (1993) 440 Eq. (2)
double amu_lightlept(double ml)
{
  double t1, t2, t3, t4, t5, t6, t7, ml2, ml4, ln, l2, a4;
  a4 = 0.51747906167389938633;//PolyLog[4, 1/2]
  ml2 = ml*ml;
  ml4 = ml2*ml2;
  ln = -log(ml);
  l2 = log(2.0);
  t1 = 2.0/3*pi*pi*ln;
  t2 = 59.0/270*pi*pi*pi*pi-3.0*ZETA3-10.0/3*pi*pi+2.0/3;
  t3 = ml*(4.0/3*pi*pi*ln-196.0/3*pi*pi*l2+424.0/9*pi*pi);
  t4 = ml2*(-2.0/3*ln*ln*ln+(pi*pi/9-20.0/3)*ln*ln-(16.0/135*pi*pi*pi*pi+4.0*ZETA3-32.0/9*pi*pi+61.0/3)*ln
	      +4.0/3*ZETA3*pi*pi-61.0/270*pi*pi*pi*pi+3.0*ZETA3+25.0/18*pi*pi-283.0/12);
  t5 = ml2*ml*(10.0/9*pi*pi*ln-11.0/9*pi*pi);
  t6 = ml2*ml2*(7.0/9*ln*ln*ln+41.0/18*ln*ln+13.0/9*pi*pi*ln+517.0/108*ln+ZETA3/2+191.0/216*pi*pi+13283.0/2592);
  printf("# t1=%.11lg\n# t2=%.11lg\n# t3=%.11lg\n# t4=%.11lg\n# t5=%.11lg\n# t6=%.11lg\n# sum= %.11lg\n", t1, t2, t3, t4, t5, t6, t1+t2+t3+t4+t5+t6);
  /* vac pol:
  t1 = 2.0/9*ln*ln;
  t2 = (ZETA3-2.0/3*pi*pi*l2+pi*pi/9+31.0/27)*ln;
  t3 = 11.0/216*pi*pi*pi*pi-2.0/9*pi*pi*l2*l2-8.0/3*a4-1.0/9*l2*l2*l2*l2 -3.0*ZETA3+5.0/3*pi*pi*l2-25.0/18*pi*pi+1075.0/216;
  t4 = ml*(-13.0/18*pi*pi*pi-16.0/9*pi*pi*l2+3199.0/1080*pi*pi);
  t5 = ml2*(10.0/3*ln*ln-11.0/9*ln-14.0/3*pi*pi*ln2-2.0*ZETA3+49.0/12*pi*pi-131.0/54);
  t6 = ml*ml2*(4.0/3*pi*pi*ln+35.0/12*pi*pi*pi-16.0/3*pi*pi*l2-5771.0/1080*pi*pi);
  t7 = ml2*ml2*(-25.0/9*ln*ln*ln-1369.0/180*ln*ln+(-2.0*ZETA3+4*pi*pi*l2-269.0/144*pi*pi-7496.0/675)*ln
		    -43.0/108*pi*pi*pi*pi+8.0/9*pi*pi*l2*l2+80.0/3*a4+10.0/9*l2*l2*l2*l2+411.0/32*ZETA3+89.0/48*pi*pi*l2-1061.0/864*pi*pi-274511.0/54000);
  */
  return(AlfQED*AlfQED*AlfQED/(pi*pi*pi)*(t1+t2+t3+t4+t5+t6));
  // return((t1+t2+t3+t4+t5+t6)); // should give 20.9479242 for ml = 1.0/206.768262.
}

void ipihatFermLoop_antisym(double x[4], double y[4], double vpihat[6][4][4][4])
{
  double vpihat1[4][4][4][4][4], vpihatr1[4][4][4][4];
  double vpihat2[4][4][4][4][4], vpihatr2[4][4][4][4];
  double vpihat3[4][4][4][4][4], vpihatr3[4][4][4][4];
  double ymx[4], mx[4];
  int rho, sig, mu, nu, lda, rhosig;
  
  for(mu=0;mu<4;mu++) {
    mx[mu] = -x[mu];
    ymx[mu] = y[mu] - x[mu];
  }

  pihat1(x,y,vpihat1,vpihatr1);
  pihat1(ymx,mx,vpihat2,vpihatr2);
  pihat1(mx,ymx,vpihat3,vpihatr3);
  rhosig=-1;
  for(rho=0;rho<4;rho++) {
    for(sig=rho+1;sig<4;sig++) {
      ++rhosig;
      for(mu=0;mu<4;mu++) {
	for(nu=0;nu<4;nu++) {
	  for(lda=0;lda<4;lda++) {
	    vpihat[rhosig][mu][nu][lda] = 0.5*(vpihat1[rho][mu][nu][lda][sig]
	      + vpihat2[rho][nu][lda][mu][sig] + x[rho]*vpihatr2[nu][lda][mu][sig]
	      + vpihat3[rho][lda][nu][mu][sig] + x[rho]*vpihatr3[lda][nu][mu][sig]
	      - (vpihat1[sig][mu][nu][lda][rho]
	      + vpihat2[sig][nu][lda][mu][rho] + x[sig]*vpihatr2[nu][lda][mu][rho]
		 + vpihat3[sig][lda][nu][mu][rho] + x[sig]*vpihatr3[lda][nu][mu][rho]));
	  }
	}
      }
    }
  }
}


// returns pihat1(rhosig,mu,nu,lambda)
void pihat1(double x[4], double y[4], double vpihat1[4][4][4][4][4], double vpihatr1[4][4][4][4])
{
  double g[4], h[4][4], fhatv[4][4], fv[4][4][4];
  double p, q[4], ell[4][4];
  double xmy[4], xn, yn, xmyn, norm, t1[4], tr1, bx[3], by[3], bxmy[3];
  int alf, bta, gam, dta, lda, sig, rho, mu, nu;

  norm = 8.2336243479292e-07; // 2/(2*pi)**8

  for(alf=0;alf<4;alf++) {
    xmy[alf] = x[alf] - y[alf];
  }
  xn = v4norm(x);
  yn = v4norm(y);
  xmyn = v4norm(xmy);

  bessk(2, xn, bx); bx[1] /= xn; bx[2] /= (xn*xn);
  bessk(2, yn, by); by[1] /= yn; by[2] /= (yn*yn);
  bessk(2, xmyn, bxmy); bxmy[1] /= xmyn; bxmy[2] /= (xmyn*xmyn);

  for(rho=0;rho<4;rho++) {
    g[rho] = pi*pi*y[rho]*by[0];
    for(dta=0;dta<4;dta++) {
      fhatv[rho][dta] = pi*pi*(y[rho]*y[dta]*by[1]+kro(rho,dta)*by[0]);
      h[rho][dta] = pi*pi*(y[rho]*y[dta]*by[1]-kro(rho,dta)*by[0]);
      for(gam=0;gam<4;gam++) {
	fv[rho][dta][gam] = pi*pi*(by[2]*y[gam]*y[dta]*y[rho]+(kro(rho,dta)*y[gam]-kro(gam,rho)*y[dta]-kro(gam,dta)*y[rho])*by[1]);
      }
    }
  }

  p = 2*pi*pi*by[0];
  for(gam=0;gam<4;gam++) {
    q[gam] = 2*pi*pi*y[gam]*by[1];
    for(dta=0;dta<4;dta++) {
      ell[gam][dta] = 2.0*pi*pi*(y[gam]*y[dta]*by[2]-kro(gam,dta)*by[1]);
    }
  }


  for(sig=0;sig<4;sig++) {
    for(mu=0;mu<4;mu++) {
      for(nu=0;nu<4;nu++) {
	for(lda=0;lda<4;lda++) {

	  for(rho=0;rho<4;rho++) t1[rho] = 0.0;	  
	  tr1 = 0.0;
	  for(alf=0;alf<4;alf++) {
	    for(bta=0;bta<4;bta++) {
	      for(gam=0;gam<4;gam++ ) {
		for(dta=0;dta<4;dta++) {
		  for(rho=0;rho<4;rho++) t1[rho] += -x[alf]*xmy[bta]*bx[2]*bxmy[2]*fv[rho][dta][gam]*g8[alf][mu][bta][nu][gam][sig][dta][lda];
		  tr1 += -x[alf]*xmy[bta]*bx[2]*bxmy[2]*ell[gam][dta] *g8[alf][mu][bta][nu][gam][sig][dta][lda];
		}
	      }
	    }
	  }
	  for(rho=0;rho<4;rho++) t1[rho] += bx[1]*bxmy[1]*g[rho]*g8[0][0][0][0][mu][nu][sig][lda];
	  tr1 += bx[1]*bxmy[1]*p*g8[0][0][0][0][mu][nu][sig][lda];
	  for(alf=0;alf<4;alf++) {
	    for(bta=0;bta<4;bta++) {
	      for(rho=0;rho<4;rho++) t1[rho] +=-x[alf]*xmy[bta]*bx[2]*bxmy[2]*g[rho]*g8[0][0][alf][mu][bta][nu][sig][lda];
	      tr1 += -x[alf]*xmy[bta]*bx[2]*bxmy[2]*p*g8[0][0][alf][mu][bta][nu][sig][lda];
	    }
	  }
	  for(alf=0;alf<4;alf++) {
	    for(gam=0;gam<4;gam++) {
	      for(rho=0;rho<4;rho++) t1[rho] += -x[alf]*bx[2]*bxmy[1]*h[rho][gam]*g8[0][0][alf][mu][nu][gam][sig][lda];
	      tr1 += -x[alf]*bx[2]*bxmy[1]*q[gam]*g8[0][0][alf][mu][nu][gam][sig][lda];
	    }
	  }
	  for(bta=0;bta<4;bta++) {
	    for(gam=0;gam<4;gam++) {
	      for(rho=0;rho<4;rho++) t1[rho] += xmy[bta]*bx[1]*bxmy[2]*h[rho][gam]*g8[0][0][mu][bta][nu][gam][sig][lda];
	      tr1 += xmy[bta]*bx[1]*bxmy[2]*q[gam]*g8[0][0][mu][bta][nu][gam][sig][lda];
	    }
	  }
	  for(alf=0;alf<4;alf++) {
	    for(dta=0;dta<4;dta++) {
	      for(rho=0;rho<4;rho++) t1[rho] += -x[alf]*bx[2]*bxmy[1]*fhatv[rho][dta]*g8[0][0][alf][mu][nu][sig][dta][lda];
	      tr1 += -x[alf]*bx[2]*bxmy[1]*q[dta]*g8[0][0][alf][mu][nu][sig][dta][lda];
	    }
	  }
	  for(bta=0;bta<4;bta++) {
	    for(dta=0;dta<4;dta++) {
	      for(rho=0;rho<4;rho++) t1[rho]+= xmy[bta]*bx[1]*bxmy[2]*fhatv[rho][dta]*g8[0][0][mu][bta][nu][sig][dta][lda];
	      tr1 += xmy[bta]*bx[1]*bxmy[2]*q[dta]*g8[0][0][mu][bta][nu][sig][dta][lda]; // HERE
	    }
	  }
	  for(gam=0;gam<4;gam++ ) {
	    for(dta=0;dta<4;dta++) {
	      for(rho=0;rho<4;rho++) t1[rho] += bx[1]*bxmy[1]*fv[rho][dta][gam]*g8[0][0][mu][nu][gam][sig][dta][lda];
	      tr1 += bx[1]*bxmy[1]*ell[gam][dta]*g8[0][0][mu][nu][gam][sig][dta][lda];
	    }
	  }
	  
	  for(rho=0;rho<4;rho++)    vpihat1[rho][mu][nu][lda][sig] = norm*t1[rho];
	  vpihatr1[mu][nu][lda][sig] = norm*tr1 ;
	}
      }
    }
  }
  

}



void init_g8()
{
  int idx[8], alf, bet, dta, mu, nu, rho, sig, lda;

  for(mu=0;mu<4;mu++) {
    for(nu=0;nu<4;nu++) {
      for(sig=0;sig<4;sig++) {
	for(lda=0;lda<4;lda++) {
	  for(rho=0;rho<4;rho++) {
	    for(alf=0;alf<4;alf++) {
	      for(bet=0;bet<4;bet++) {
		for(dta=0;dta<4;dta++) {
		  idx[0]=dta;
		  idx[1]=rho;
		  idx[2]=sig;
		  idx[3]=mu;
		  idx[4]=alf;
		  idx[5]=nu;
		  idx[6]=bet;
		  idx[7]=lda;
		  g8[dta][rho][sig][mu][alf][nu][bet][lda] = 4*diractrace(8, idx);
		}
	      }
	    }
	  }
	}
      }
    }
  }

}

#endif



/*
void bessk(int n, double x, double *bv)
{
  int j;
  double bk, bkm, bkp, tox, y, bi0, bi1, sq;

  tox=2.0/x;

  if(x<=2.0) {
    y=x/3.75;
    y*=y;
    bi0=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
	 +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    bi1=x*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
       +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));

    y=0.25*x*x;
    sq = log(x/2.0);
    bkm=(-sq*bi0)+(-0.57721566+y*(0.42278420
	+y*(0.23069756+y*(0.3488590e-1+y*(0.262698e-2
	  +y*(0.10750e-3+y*0.74e-5))))));

    bk = (sq*bi1)+(1.0/x)*(1.0+y*(0.15443144
       +y*(-0.67278579+y*(-0.18156897+y*(-0.1919402e-1
	 +y*(-0.110404e-2+y*(-0.4686e-4)))))));
  } else {
    sq = (exp(-x)/sqrt(x));
    bkm=sq*(1.25331414+tox*(-0.7832358e-1
	 +tox*(0.2189568e-1+tox*(-0.1062446e-1+tox*(0.587872e-2
      +tox*(-0.251540e-2+tox*0.53208e-3))))));

    bk = sq*(1.25331414+tox*(0.23498619
	  +tox*(-0.3655620e-1+tox*(0.1504268e-1+tox*(-0.780353e-2
      +tox*(0.325614e-2+tox*(-0.68245e-3)))))));
  }

  bv[0]=bkm;
  bv[1]=bk;
  for(j=1;j<n;j++) {
    bkp=bkm+j*tox*bk;
    bkm=bk;
    bv[j+1]=bk=bkp;
  }
}
*/

