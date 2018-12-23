// gcc -o vva4enhung  vva4enhung.c  -lm

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getL.h"

//#define VVA

#ifndef PI
#define PI 3.14159265358979323846
#define PI5 (0.5*PI)
#define pi PI
#define pi2 (PI*2.0)
#define PI8 25.1327412287183459077011
#define ZETA3 1.2020569031595942854
#endif

#define v4norm(x) sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
#define kro(mu,nu) (mu==nu ? 1 : 0)
#define Iov2pisq 0.050660591821168885722

//  #define CONTI

void bessk(int n, double x, double *bv);
void s0prop(double ma, double x, double a[2]);
void s0propL(double ma, double xv[4], double Ll[4], double a[5]);
double getavv_floopL(double ma, double Ll[4], double x[4], double y[4], int mu, int nu, int rho);
void set_epstensor();

void momentum_sum3d(int ns, void (*func)(double,double,double,double *,double *), double *p, double *sum, 
		    int L1, int L2, int L3, double theta1, double theta2, double theta3);
void thermprop(double px, double py, double pz, double *p, double prop[5]);
void prop_smd(double px, double py, double pz, double *p, double qr[5]);


static int Epst[4][4][4][4];

// this program computes the < V_mu(x) V_nu(y) A_rho(0) > 
// for a free massive Dirac fermion on a 4d torus, antiperiodic bc in time (=direction 0) and periodic in space.
// If the flag CONTI is switched on, the calculation uses continuum propagators;
// if not, Wilson propagators are used. The treelevel O(a) improvement of the mass and of the vector
// and axial-vector currents is used.
/*
int main(int argc, char *argv[])
{
  double ma, al, a9, check, vvp[4], vvm[4], deriv[4], buf, stp, f1, f2, f, ya[4];
  double xmy[4], x[4], y[4], z[4], Ll[4], a[5], xn, yn, sum, sumabs, eps;
  int mu, nu, sig, lda, kmax, rho, i;
  double x1, x2, x3, mx[4], ymx[4], avv, avv2, bA, bV;

  set_epstensor();

  al =  1.0 ; // to take the continuum limit, make al smaller and smaller
  ma=0.076*al;   // the mass in lattice units
  Ll[1] = Ll[2] = Ll[3] = 32/al;
  Ll[0] = 2.0*Ll[1];
#ifdef CONTI
  printf("Continuum calculation\n");
#else
  printf("Lattice calculation\n");
#endif 
  for(i=0;i<4;i++) printf("# Ll[%d]= %lf\n", i, Ll[i]);

  x[0]= 5/al;
  x[1]= 4/al;
  x[2]= 3/al; 
  x[3]= 4/al;
  y[0]= -4/al;
  y[1]= -4/al;
  y[2]= 4/al;
  y[3]= 6/al;

  for(i=0;i<4;i++)  xmy[i] = x[i] - y[i];
  printf("# mass= %lg\n", ma);
  printf("norm(x)= %lg\n", v4norm(x));
  printf("norm(y)= %lg\n", v4norm(y));
  printf("norm(x-y)= %lg\n", v4norm(xmy));
  for(i=0;i<4;i++) printf("# x[%d]= %lf\ty[%d]= %lf\n", i, x[i], i, y[i]);
#ifdef CONTI
  a9 = pow(ma,-9);
#else
  bA=bV=1.0;
  a9 = pow(ma,-9)*(1.0+bV*ma)*(1.0+bV*ma)*(1.0+bA*ma);// O(a) improvement of V and A currents
#endif
  rho=2;
  printf("rho= %d\n", rho);
  printf("mu\tnu\t VVA/m^9\n");
  for(mu=0;mu<4;mu++) {
    for(nu=0;nu<4;nu++) {
      avv= a9*getavv_floopL(ma, Ll, x, y, mu, nu, rho);
#ifdef CHECK_BOSE_SYMM
      avv2= a9*getavv_floopL(ma, Ll, y, x, nu, mu, rho);
      printf("# %d\t%d\t%lg\t%lg\t%lg\n", mu, nu, avv, avv2, 2.0*(avv-avv2)/(avv+avv2));
#else
      printf("# %d\t%d\t%lg\n", mu, nu, avv);
#endif 
    }
  }
  printf("\n");

#if 0
  printf("Check of WI in y\n");
  eps = 1.0e-3;
  for(mu=0;mu<4;mu++) {
    for(nu=0;nu<4;nu++) {
      y[nu] += eps;
      vvp[nu]= a9*getavv_floopL(ma, Ll, x, y, mu, nu, rho);
      y[nu] -= 2*eps;
      vvm[nu]= a9*getavv_floopL(ma, Ll, x, y, mu, nu, rho);
      y[nu] += eps;
      deriv[nu] = (vvp[nu]-vvm[nu])/(2.0*eps);
    }
    sum = deriv[0] + deriv[1] + deriv[2] + deriv[3];
    sumabs = fabs(deriv[0]) + fabs(deriv[1]) + fabs(deriv[2]) + fabs(deriv[3]);
    printf("mu= %d: WI in Vnu(y) gives  %lg (sumabs= %lg, ratio= %lg)\n", mu, sum, sumabs, sum/sumabs);
  }
  printf("\nCheck of WI in x\n");
  for(nu=0;nu<4;nu++) {
    for(mu=0;mu<4;mu++) {
      x[mu] += eps;
      vvp[mu]= a9*getavv_floopL(ma, Ll, x, y, mu, nu, rho);
      x[mu] -= 2*eps;
      vvm[mu]= a9*getavv_floopL(ma, Ll, x, y, mu, nu, rho);
      x[mu] += eps;
      deriv[mu] = (vvp[mu]-vvm[mu])/(2.0*eps);
    }
    sum = deriv[0] + deriv[1] + deriv[2] + deriv[3];
    sumabs = fabs(deriv[0]) + fabs(deriv[1]) + fabs(deriv[2]) + fabs(deriv[3]);
    printf("nu= %d: WI in Vmu(x) gives  %lg (sumabs= %lg, ratio= %lg)\n", nu, sum, sumabs, sum/sumabs);
  }
#endif

  return 1;
}
*/


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


// massive case: a= m^2 K2(m|x|)/ (4*pi^2 x^2)   and   b= m^2 K1(m|x|)/ (4*pi^2|x|)
void s0prop(double ma, double x, double a[2])
{
  double bv[3], tmp, msq;

  if(ma==0.0) {
    a[0] = Iov2pisq/(x*x*x*x);
    a[1]=0.0;
  } else {
    bessk(2, ma*x, bv);
    tmp = 0.5/(x*x);
    msq = ma*ma;
    a[0] = tmp*Iov2pisq*msq*bv[2];
    a[1] = tmp*x*Iov2pisq*msq*bv[1];
  }
}


#define Ns 8
#define N0 Ns
#define N1 Ns
#define N2 Ns
#define N3 Ns

void s0propL(double ma, double xv[4], double Ll[4], double a[5])
{
  double bv[3], tmp, xw[4], x, tmp2, fn, xsq, ph;
  static double ph0=-1.0; // antiperiodic bc in time
  int nx0, nx1, nx2, nx3;

  a[0]=a[1]=a[2]=a[3]=a[4]=0.0;

  if(ma==0.0) {
    for(nx0=-N0;nx0<=N0;nx0++) {
      xw[0] = xv[0] + nx0*Ll[0];
      ph = pow(ph0,nx0)*Iov2pisq;
      for(nx1=-N1;nx1<=N1;nx1++) {
	xw[1] = xv[1] + nx1*Ll[1];
	for(nx2=-N2;nx2<=N2;nx2++) {
	  xw[2] = xv[2] + nx2*Ll[2];
	  for(nx3=-N3;nx3<=N3;nx3++) {
	    xw[3] = xv[3] + nx3*Ll[3];
	    xsq = xw[0]*xw[0]+xw[1]*xw[1]+xw[2]*xw[2]+xw[3]*xw[3];
	    tmp = ph/(xsq*xsq);
	    a[0] += tmp*xw[0]; 
	    a[1] += tmp*xw[1]; 
	    a[2] += tmp*xw[2]; 
	    a[3] += tmp*xw[3]; 
	  }
	}
      }
    }
  }   else {
    fn = 0.5*Iov2pisq*ma*ma;
    for(nx0=-N0;nx0<=N0;nx0++) {
      xw[0] = xv[0] + nx0*Ll[0];
      ph = pow(ph0,nx0);
      for(nx1=-N1;nx1<=N1;nx1++) {
	xw[1] = xv[1] + nx1*Ll[1];
	for(nx2=-N2;nx2<=N2;nx2++) {
	  xw[2] = xv[2] + nx2*Ll[2];
	  for(nx3=-N3;nx3<=N3;nx3++) {
	    xw[3] = xv[3] + nx3*Ll[3];
	    xsq = xw[0]*xw[0]+xw[1]*xw[1]+xw[2]*xw[2]+xw[3]*xw[3];
	    x = sqrt(xsq);
	    bessk(2, ma*x, bv);
	    tmp = ph*fn/(xsq);
	    tmp2 = tmp*bv[2];
	    a[0] += tmp2*xw[0];
	    a[1] += tmp2*xw[1];
	    a[2] += tmp2*xw[2];
	    a[3] += tmp2*xw[3];
	    a[4] += tmp*x*bv[1];
	  }
	}
      }
    }
  }
}
#undef N3
#undef N2
#undef N1
#undef Ns
#undef N0




// computes < V_mu(x) V_nu(y) A_rho(0) >
double getavv_floopL(double ma, double Ll[4], double x[4], double y[4], int mu, int nu, int rho)
{
  const double bm=-0.5;
  int alf, bta, gam, dta;
  int i, ns1, ns2, ns3;
  double xmy[4], ymz[4], xm[4];
  double smx[5], sy[5], sxmy[5], t1, sum1, p[100];

  set_epstensor();

  t1 = 0.0;

  for(alf=0;alf<4;alf++) {
    xmy[alf] = x[alf] - y[alf];
    xm[alf] = -x[alf];
  }

#ifdef CONTI
  s0propL(ma, xm, Ll, smx);
  s0propL(ma, y, Ll, sy);
  s0propL(ma, xmy, Ll, sxmy);
#else
  ns1=Ll[1];
  ns2=Ll[2];
  ns3=Ll[3];
  p[0] = ma*(1.0-bm*ma);// O(a) improvement of the mass
  p[1] = Ll[0];
  for(i=0;i<4;i++)  p[2+i] = xm[i];
  momentum_sum3d(5, prop_smd, p, smx, ns1, ns2, ns3, 0.0, 0.0, 0.0);
  for(i=0;i<4;i++)  p[2+i] = y[i];
  momentum_sum3d(5, prop_smd, p, sy, ns1, ns2, ns3, 0.0, 0.0, 0.0);
  for(i=0;i<4;i++)  p[2+i] = xmy[i];
  momentum_sum3d(5, prop_smd, p, sxmy, ns1, ns2, ns3, 0.0, 0.0, 0.0);
#endif 
  sum1=smx[4]*sxmy[4];
  for(bta=0;bta<4;bta++)  sum1 -= smx[bta]*sxmy[bta];
  for(alf=0;alf<4;alf++) t1 += Epst[mu][nu][rho][alf]*(-sy[4]*smx[4]*sxmy[alf] + sxmy[4]*sy[4]*smx[alf] + sum1*sy[alf]);

  for(alf=0;alf<4;alf++) {
    for(bta=0;bta<4;bta++) {
      t1 += Epst[nu][rho][alf][bta]*sy[bta]*(sxmy[alf]*smx[mu] + sxmy[mu]*smx[alf]);
      t1 -= sxmy[alf]*smx[bta]*(sy[nu]*Epst[mu][rho][alf][bta] + sy[rho]*Epst[mu][nu][alf][bta]);
    }
  }
  if(rho==nu) for(alf=0;alf<4;alf++) {
      for(bta=0;bta<4;bta++) for(gam=0;gam<4;gam++) {
	  t1 += sxmy[alf]*smx[bta]*Epst[mu][alf][bta][gam]*sy[gam];
	}
    }

  return(8.0*t1); 
}


#define epst3(mu,nu,rho,sig,s) {\
    Epst[mu][nu][rho][sig] = s;\
    Epst[nu][mu][rho][sig] = -s;\
    Epst[rho][nu][mu][sig] = -s;\
    Epst[mu][rho][nu][sig] = -s;\
    Epst[nu][rho][mu][sig] = s;\
    Epst[rho][mu][nu][sig] = s;}

void set_epstensor()
{
  int mu, nu, rho, sig, s;

  for(mu=0;mu<4;mu++)   for(nu=0;nu<4;nu++) 
			  for(rho=0;rho<4;rho++)   for(sig=0;sig<4;sig++)  Epst[mu][nu][rho][sig]=0;

  mu=0; nu=1; rho=2; sig=3; s=1;
  epst3(mu,nu,rho,sig,s);
  mu=0; nu=1; rho=3; sig=2; s=-1;
  epst3(mu,nu,rho,sig,s);
  mu=0; nu=3; rho=2; sig=1; s=-1;
  epst3(mu,nu,rho,sig,s);
  mu=3; nu=1; rho=2; sig=0; s=-1;
  epst3(mu,nu,rho,sig,s);

}


// fermion propagator in the time-momentum representation,
// for antiperiodic bc in time
void thermprop(double px, double py, double pz, double *p, double prop[5])
{  
  double Ap, Bp, C, Cs, phx, phy, phz, pox, poy, poz, phx2, phy2, phz2;
  double phat2, c33, denom, brack, th, wp, am, L0, x0, Dp, chL0, fs, sh;

  am = p[0];
  L0 = p[1];
  x0 = p[2];

  phx = sin(0.5*px);
  phy = sin(0.5*py);
  phz = sin(0.5*pz);
  phx2 = phx*phx;
  phy2 = phy*phy;
  phz2 = phz*phz;
  phat2 = 4.0*(phx2 + phy2 + phz2);
  pox = sin(px);
  poy = sin(py);
  poz = sin(pz);
  Cs = p[0]+0.5*phat2;
  Ap = 1.0+Cs;
  Bp = p[0]*p[0]+(1.0+p[0])*phat2+8.0*(phx2*phy2+phx2*phz2+phy2*phz2);
  C = Cs - 0.5*Bp/Ap;
  wp = 2.0*asinh(0.5*sqrt(Bp/Ap));
  sh = sinh(wp);
  Dp = 2.0*Ap*sinh(wp);
  chL0 = cosh(0.5*wp*L0);
  fs = sinh(wp*(0.5*p[1]-p[2]))/(Dp*chL0);
  prop[1] = -pox*fs;
  prop[2] = -poy*fs;
  prop[3] = -poz*fs;
  if(x0==0.0) {
    prop[0] = 0.0;
    prop[4] = fs*C + sinh(wp)/Dp;
  } else {
    prop[0] = sinh(wp)*cosh(wp*(0.5*p[1]-p[2]))/(Dp*chL0);
    prop[4] = fs*C;
  }
}

// computes coeff. of gamma{0,1,2,3} and of the unit matrix in the summand 
// for the position-space quark propagator
void prop_smd(double px, double py, double pz, double *p, double qr[5])
{
  double prop[5], x, y, z, co, si, xdotp;
  
  x = p[3];
  y = p[4];
  z = p[5];

  thermprop(px, py, pz, p, prop);
  xdotp = x*px + y*py + z*pz;
  co = cos(xdotp);
  si = sin(xdotp);
  qr[0] +=  co * prop[0];
  qr[1] += -si * prop[1];
  qr[2] += -si * prop[2];
  qr[3] += -si * prop[3];
  qr[4] +=  co * prop[4];

}

// output is sum[0..(ns-1)]
void momentum_sum3d(int ns, void (*func)(double,double,double,double *,double *), double *p, double *sum, 
		      int L1, int L2, int L3, double theta1, double theta2, double theta3)
{
  double p1, p2, p3, mom_unit1, mom_unit2, mom_unit3;
  double mom_offset1, mom_offset2, mom_offset3, fn;
  int i1, i2, i3, i, iflag;

  mom_unit1 = pi2/L1;
  mom_unit2 = pi2/L2;
  mom_unit3 = pi2/L3;

  mom_offset1 = theta1/L1;
  mom_offset2 = theta2/L2;
  mom_offset3 = theta3/L3;

  if(p[2]<0.0) {
    iflag = 1;
    p[2] = -p[2];
  }  else iflag = 0;
  for(i=0;i<ns;i++)  sum[i]=0.0;

  for(i1=0;i1<L1;i1++) {
    p1 = mom_offset1 + i1*mom_unit1;
    for(i2=0;i2<L2;i2++) {
      p2 = mom_offset2 + i2*mom_unit2;
      for(i3=0;i3<L3;i3++) {
	p3 = mom_offset3 + i3*mom_unit3;
	func(p1,p2,p3,p,sum);
      }
    }
  }
  fn = 1.0/(L1*L2*L3);
  for(i=0;i<ns;i++)   sum[i] *= fn;
  if(iflag) {
    sum[0]=-sum[0];
    p[2] = -p[2];
  }
}

