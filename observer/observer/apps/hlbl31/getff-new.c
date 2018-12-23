#include "getL.h"


// returns the form factor, given the coefficients fm[0..(nf-1)] and fp[0..(nf-1)]
// e.g. fm = alpha^{(3)}_{m-}  and fp = alpha^{(3)}_{m+}
// nf = length of vectors fm and fp
// nff = length of vectors cbv and ff
// nm = index of the sum sigma that appeared in the integrand
// mm = whether the form factor is expanded in ChebU, ChebU' or ChebU'' (mm=0,1,2)
// ndy and ndcb = # derivatives with respect to y and cos(beta) respectively
void getff(int nff, int nf, int nm, int mm, int ndy, int ndcb, 
	   double x, double y, double *cbv, double *fm, double *fp, double *ff)
{
  double y1, yp, facm, facp, fval[N];
  int j, k, nmr, mshm, mshp;

  if(ndy<0 || ndy>1 ) { printf("getff: Unexpected value of ndy=%d \n", ndy); exit(0);}
  if(mm+ndcb>3) { printf("getff: Unexpected value of mm=%d and ndcb=%d\n", mm, ndcb); exit(0);}
  if(nf+mm>N-2) { printf("getff: Too many coeffs in routine getff. Allocated more memory.\n"); exit(0);}

  y1 = 1.0/y;
  nmr = nm%10;
  switch(nmr) {
  case 0: yp = 1.0; mshm=2; mshp=0; break;
  case 1: yp = 1.0; mshm=2; mshp=0; break;
  case 2: yp =  y1; mshm=3; mshp=-1; break;
  case 3: yp = 1.0; mshm=2; mshp=0; break;
  case 4: yp =   y; mshm=1; mshp=1; break;
  case 5: yp =y1*y1; mshm=4; mshp=-2; break;
  default: printf("Unexpected value of nm= %d\n", nm); exit(0); break;
  }
  yp *= (ndy==1 ? y1 : 1.0);
  facm= y1*y1*yp;
  facp= yp;
  for(j=0;j<mm;j++) {
    fval[j] = 0.0;
    facm *= y1;
    facp *= y;
  }
  for(j=mm;j<nf+mm;j++) {
    fval[j] = (ndy==1 ?  -(j+mshm)*facm*fm[j-mm] + (j+mshp)*facp*fp[j-mm] : facm*fm[j-mm] + facp*fp[j-mm]);
    facm *= y1;
    facp *= y;
  }
  for(k=0;k<nff;k++) ff[k] = Func_usm[mm + ndcb](nf+mm, cbv[k], fval); 
}


// returns Sum[co[k]*ChebyshevU[k,x],{k,0,nk-1}]
// input: array co[0...(nk-1)], not modified
// assumptions: nk>=1
// implements the Clenshaw Recurrence Formula
// based on Numerical Recipes in C 1992, section 5.5
// U_{n+1}(x) = 2x U_n(x) - U_{n-1}(x)
// i.e.  alpha(n,x)=2*x     beta(n,x)=-1
double chebUsum(int nk, double x, double *co)
{
  int n, j;
  double ya, yb, two, ytmp, twox;

  n = nk-1;
  ya=0.0;
  yb=0.0;
  twox = 2.0*x;
  for(j=n;j>=1;j--) {
    ytmp = yb;
    yb = twox*yb - ya + co[j];
    ya = ytmp;
  }
  return(-ya + twox*yb + co[0]);
}


// Clenshaw for Sum[co[k]*Derivative[0,1][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/n       beta(n,z) = -(n+2)/n
// y_k = alf(k,x)*y_{k+1} + beta(k+1,x)*y_{k+2} + co[k]
double dchebUsum(int nk, double x, double *co)
{
  int n, j;
  double ya, yb, two, ytmp, twox;

  n = nk-1;
  ya=0.0;
  yb=0.0;
  twox = 2.0*x;

  for(j=n;j>=1;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/j - (j+3)*ya/(j+1) + co[j]; 
    ya = ytmp;
  }
  return(2.0*yb);
}


// Clenshaw for Sum[co[k]*Derivative[0,2][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/(n-1)       beta(n,z) = -(n+3)/(n-1)
double ddchebUsum(int nk, double x, double *co)
{
  int n, j;
  double ya, yb, two, ytmp, twox;

  n = nk-1;
  ya=0.0;
  yb=0.0;
  twox = 2.0*x;
  for(j=n;j>=2;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/(j-1) - (j+4)*ya/j + co[j];
    ya = ytmp;
  }
  return(8.0*yb);
}


// Clenshaw for Sum[co[k]*Derivative[0,3][ChebyshevU][k,x],{k,0,nk-1}]
// alf(n,z) = 2z(n+1)/(n-2)       beta(n,z) = -(n+4)/(n-2)
double dddchebUsum(int nk, double x, double *co)
{
  int n, j;
  double ya, yb, two, ytmp, twox;

  n = nk-1;
  ya=0.0;
  yb=0.0;
  twox = 2.0*x;
  for(j=n;j>=3;j--) {
    ytmp = yb;
    yb = twox*(j+1)*yb/(j-2) - (j+5)*ya/(j-1) + co[j];
    ya = ytmp;
  }
  return(48.0*yb);
}


void spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
	int i,k;
	double p,qn,sig,un,*u;

	u = (double *) malloc((n-1)*sizeof(double));
	if (yp1 > 0.99e30)
		y2[0]=u[0]=0.0;
	else {
		y2[0] = -0.5;
		u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=1;i<n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free(u);
}


void splint(double *xa, double *ya, double *y2a, int n, double x, double *y)
{
	int klo,khi,k;
	double h,b,a;

	klo=0;
	khi=n-1;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) { printf("Bad XA input to routine SPLINT"); exit(0); }
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	y[0]=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	y[1]=(ya[khi]-ya[klo])/h - (3.0*a*a-1.0)/6.0*h*y2a[klo] + (3.0*b*b-1.0)/6.0*h*y2a[khi];
}


#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m)
{
  int *indxc,*indxr,*ipiv;
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,temp;
  
  indxc= (int *) malloc(n* sizeof(int));
  indxr= (int *) malloc(n* sizeof(int));
  ipiv = (int *) malloc(n* sizeof(int));


  for (j=0;j<n;j++) ipiv[j]=0;
  for (i=0;i<n;i++) {
    big=0.0;
   
    for (j=0;j<n;j++)
      if (ipiv[j] != 1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big=fabs(a[j][k]);
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] > 1) { printf("gaussj: Singular Matrix-1"); exit(0);}
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) { printf("gaussj: Singular Matrix-2"); exit(0);}
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
    for (ll=0;ll<n;ll++)
      if (ll != icol) {
	dum=a[ll][icol];
	a[ll][icol]=0.0;
	for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
	for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
  }
  for (l=n-1;l>=0;l--) {
    if (indxr[l] != indxc[l])
      for (k=0;k<n;k++)
	SWAP(a[k][indxr[l]],a[k][indxc[l]]);
  }
  free(ipiv);
  free(indxr);
  free(indxc);
}

#define NFF 512
static double XX[NFF], YY[NFF], xb_ff, xstp_ff, ystp_ff, xstp3_ff;
// static   double xstp, ystp, xh, yh;
static   int nstpx_ff, nstpy_ff, mul_ff, ib_ff;

void init_ff(double xh, double yh, int nstpx, int nstpy)
{
  int i, nx, j, k, kq, nmc;
  double xstp, ystp, x, y, stmp, xbuf;
  char c_char[81];
  FILE *fp;

  Func_usm[0] = chebUsum;
  Func_usm[1] = dchebUsum;
  Func_usm[2] = ddchebUsum;
  Func_usm[3] = dddchebUsum;

  nstpx_ff = nstpx;
  nstpy_ff = nstpy;
  ystp_ff = ystp = yh/nstpy;
  xstp_ff =  xstp = ystp; // xh/nstpx;
  xb_ff = 0.363333333334;
  ib_ff = 15;
  mul_ff = 3;
  xstp3_ff = mul_ff*xstp;

#ifdef ECHO
  printf("xstp= %lg   ystp= %lg\n", xstp, ystp);
#endif
  for(i=0;i<ib_ff;i++) XX[i] = (i+1)*xstp;
  for(i=ib_ff;i<nstpx;i++) XX[i] = (i - ib_ff + 6)*xstp3_ff;
  for(i=0;i<nstpy;i++) YY[i] = (i+1)*ystp;

#ifndef READ_ON_THE_FLY
  for(nmc=0;nmc<NFFA;nmc++) {  
    Ffp[nmc] = (double ***) malloc(nstpx*sizeof(double **));
    for(i=0;i<nstpx;i++) Ffp[nmc][i] = (double **) malloc(nstpy*sizeof(double *));
    for(i=0;i<nstpx;i++) {
      nx= nfx(XX[i]);
      for(j=0;j<nstpy;j++) Ffp[nmc][i][j] = (double *) malloc(nx*sizeof(double));
    }
    Ffm[nmc] = (double ***) malloc(nstpx*sizeof(double **));
    for(i=0;i<nstpx;i++) Ffm[nmc][i] = (double **) malloc(nstpy*sizeof(double *));
    for(i=0;i<nstpx;i++) {
      nx= nfx(XX[i]);
      for(j=0;j<nstpy;j++) Ffm[nmc][i][j] = (double *) malloc(nx*sizeof(double));
    }
    
    sprintf(c_char, "in/ffxyh%d.txt", nmc); 
    //      sprintf(c_char, "ffxy%d.txt", nmc); 
    //    sprintf(c_char, "ffy%d.txt", nmc); 
    //    sprintf(c_char, "ff%d.txt", nmc);
    // printf("Reading from file %s\n", c_char);
    fp = fopen(c_char, "r");
    for(i=0;i<nstpx;i++) {
      nx= nfx(XX[i]);
      for(j=0;j<nstpy;j++) {
	fscanf(fp, "%lg %lg", &x, &y);
	// fscanf(fp, "%lg %lg", XX+i, YY+j);
	if(fabs(XX[i]- x)>1.e-6 || fabs(YY[j]-y)>1.e-6) { printf("x= %lg vs. XX[%d]= %lg   y= %lg vs. YY[%d]= %lg\n", x, i, XX[i], y, j, YY[j]); exit(0); }
	for(k=0;k<nx;k++)  {
	  fscanf(fp, "%d %lg %lg", &kq, Ffm[nmc][i][j]+k, Ffp[nmc][i][j]+k);
	  if(kq!=k)  { printf("Problem at: nmc= %d x= %lg y= %lg kq= %d .neq. k= %d\n", nmc, x, y, kq, k); exit(0); }
	}
      }
    }
    fclose(fp);
  }

  //  for(k=0;k<nfx(XX[0]);k++) printf("%d\t%.11lg\t%.11lg\n", k,  Ffm[6][0][1][k],   Ffp[6][0][1][k]);
#endif

#ifdef TAYLORX
  //  fp=fopen("g_tay.txt", "r");
  fp=fopen("in/g_tay810.txt", "r");
  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, G0dy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dy_tay[j]); }
  Stp_tay = stmp= YY_tay[NY_TAY-1]/NY_TAY;
  printf("# Stp_tay= %.11lg\n", Stp_tay);
  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, G0dx_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl2_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl2dy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl21_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl21dy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl3_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, Gl3dy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, G21_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg", YY_tay+j, G21dy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G22A_tay+j, G22B_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G22Ady_tay+j, G22Bdy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G3A_tay+j, G3B_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G3Ady_tay+j, G3Bdy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G31A_tay+j, G31B_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  for(j=0;j<NY_TAY;j++) { fscanf(fp, "%lg %lg %lg", YY_tay+j, G31Ady_tay+j, G31Bdy_tay+j); } // {printf("%lf\t%.11lg\n", YY_tay[j], G0dx_tay[j]); }
  Stp_tay = YY_tay[NY_TAY-1]/NY_TAY;
  if(Stp_tay!=stmp) {printf("Trouble reading in Taylor expansion.\n"); exit(0);}
  stmp = Stp_tay;

  fclose(fp);

#endif

#ifdef TAYLORY
  fp=fopen("in/taylorY2.txt", "r");  
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha0dx_0p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha0_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha3_0p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, beta2_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha3_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha1_0p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, beta4_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha1_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha1dx_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha3dx_0p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha3dxdx_0p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, beta2dx_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha3dx_1p_taY+j); }
  for(j=0;j<NX_TAY;j++) { fscanf(fp, "%lg %lg", &xbuf, alpha1dx_0p_taY+j); }
  fclose(fp);

#if 0
  for(j=0;j<5;j++)  printf("### x= %lf\talpha0dx_0p_taY[%d]= %lg\n", XX[j], j, alpha0dx_0p_taY[j]);
  for(j=0;j<5;j++)  printf("### x= %lf\talpha0_1p_taY[%d]= %lg\n", XX[j], j, alpha0_1p_taY[j]);
  for(j=0;j<5;j++)  printf("### x= %lf\talpha1dx_0p_taY[%d]= %lg\n", XX[j], j, alpha1dx_0p_taY[j]);
  j=37; 
  printf("alpha3dxdx_0p_taY[%d]= %lg  discrete derivative = %lg \n", j,
		alpha3dxdx_0p_taY[j], (alpha3dx_0p_taY[j+1]-alpha3dx_0p_taY[j-1])/(XX[j+1]-XX[j-1]));
  printf("beta2dx_1p_taY[%d]= %lg  discrete derivative = %lg \n", j,
		beta2dx_1p_taY[j], (beta2_1p_taY[j+1]-beta2_1p_taY[j-1])/(XX[j+1]-XX[j-1]));
#endif

#endif
}


double extractff(int nm, int ndy, int ndcb, double x, double cb, double y)
{
  int ix1;
  int ix, iy, nmc, nmd, use_x_derivs, use_y_derivs;
  int flag_hx, flag_hy, nstpx, nstpy;
  FILE *fp, *gp;
  char c_char[81], d_char[81];
  double f1iy, f2iy, g1iy, g2iy, xbuf, x1, x2, y1, y2, a, fi, dx, stpx, stpy;
  use_y_derivs = (ndy==0);
  nstpx = nstpx_ff;
  nstpy = nstpy_ff;
  stpx = xstp_ff; // XX[1]-XX[0];
  stpy = ystp_ff; // YY[1]-YY[0];
  //  for(ix=0;ix<nstpx;ix++) if(XX[ix]>x) break;   ix1 = ix;
  ix = (x<xb_ff ? x/stpx : ib_ff + (x-xb_ff)/xstp3_ff);
  //  if(ix1!=ix) {printf("x= %lf : ix1=%d not equal ix=%d.\n", x, ix1, ix); exit(0); }
  if(x>XX[nstpx-1]) {
    printf("# Warning: x is higher than upper edge.\n");
    flag_hx = 1;
  } else flag_hx = 0;
  ix = (ix>=1 ? ix-1 : ix);
  //if(ix>=1 && ix<nstpx-1 && (x<XX[ix] || x>XX[ix+1])) { printf("Sthg wrong in extractff.\n"); exit(0);} // just a check; can remove later 

  if(y>YY[nstpy-1]) {
    printf("# Warning: y is higher than upper edge.\n");
    flag_hy = 1;
  } else  flag_hy = 0;
  //  for(iy=0;iy<nstpy;iy++) if(YY[iy]>y) break; 
  iy = y/stpy;
  iy = (iy>=1 ? iy-1 : iy);
  //  printf("ix= %d   iy= %d\n", ix, iy);
  // printf("x= %lg  y= %lg\n", x, y);  printf("XX[0]= %lg XX[1]= %lg XX[2]= %lg\n", XX[0], XX[1], XX[2]);  return 0; 

  nmc = Iconv[nm];
#ifdef READ_ON_THE_FLY
  sprintf(c_char, "in/ff%d.txt", nmc);
  //  printf("# extractff: Opening file %s ...\n", c_char);
  fp = fopen(c_char, "r");
  if(nm<10 || nm==12 || nm==13) {
    use_x_derivs = 1;
    nmd = Iconv[nm+10];
    sprintf(d_char, "in/ff%d.txt", nmd);
    //    printf("# extractff: Opening file %s ...\n", d_char);
    gp = fopen(d_char, "r");
  } else use_x_derivs = 0;

  f1iy = accessf(fp, flag_hy, use_y_derivs, ix*nstpy+iy, nm, ndy, ndcb, &x1, &y1, &y2, x, &cb, y);
  //  printf("f1iy= %lg\n", f1iy);
  if(use_x_derivs) {
    g1iy = accessf(gp, flag_hy, use_y_derivs, ix*nstpy+iy, nm, ndy, ndcb, &x1, &y1, &y2, x, &cb, y);
    //    printf("g1iy= %lg\n", g1iy);
  }
  if(!flag_hx) {
    f2iy = accessf(fp, flag_hy, use_y_derivs, nstpy-2, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
    //    printf("f2iy= %lg\n", f2iy);
    if(use_x_derivs) {
      g2iy = accessf(gp, flag_hy, use_y_derivs, nstpy-2, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
      //  printf("g2iy= %lg\n", g2iy);
      dx = x1-x2;
      fi = interpol3(x,x1,x2,dx,f1iy,f2iy,g1iy,g2iy);
    } else {
      a = (x2-x)/(x2-x1);
      fi = a*f1iy + (1.0-a)*f2iy;
    }
  } else fi = f1iy;
  //  printf("Final result: fi= %lg\n", fi);

  fclose(fp);
  if(use_x_derivs) fclose(gp);
#else 
  if(nm<10 || nm==12 || nm==13) {
    use_x_derivs = 1;
    //    nmd = Iconv[nm+10];
  } else use_x_derivs = 0;

  f1iy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x1, &y1, &y2, x, &cb, y);
  //  printf("f1iy= %lg\n", f1iy);
  if(use_x_derivs) {
    g1iy = accessv(flag_hy, use_y_derivs, ix, iy, nm+10, ndy, ndcb, &x1, &y1, &y2, x, &cb, y);
    //    printf("g1iy= %lg\n", g1iy);
  }

  if(!flag_hx) {
    f2iy = accessv(flag_hy, use_y_derivs, ix+1, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
    //    printf("f2iy= %lg\n", f2iy);
    if(use_x_derivs) {
      g2iy = accessv(flag_hy, use_y_derivs, ix+1, iy, nm+10, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
      //  printf("g2iy= %lg\n", g2iy);
      dx = x1-x2;
      fi = interpol3(x,x1,x2,dx,f1iy,f2iy,g1iy,g2iy);
    } else {
      a = (x2-x)/(x2-x1);
      fi = a*f1iy + (1.0-a)*f2iy;
    }
  } else fi = f1iy;
  //  printf("Final result: fi= %lg\n", fi);

#endif
  return(fi);
}

// #define VERBOSE 

void readf(FILE *fp, int ny, double *xr, double *yr, double *am, double *ap, int verbose)
{
  int ij, j, k;

  for(ij=0;ij<ny;ij++) {
    fscanf(fp, "%lg %lg", xr, yr);
    if(verbose) {
#ifdef VERBOSE
      printf("# Using values x1= %lg    y1= %lg\n", xr[0], yr[0]);
#endif
      for(j=0;j<nfx(xr[0]);j++) {
	fscanf(fp, "%d %lg %lg", &k, am+j, ap+j);
#ifdef VERBOSE
	printf( "# %d %lg %lg\n", k, am[j], ap[j]);  
#endif 
      }
    } else for(j=0;j<nfx(xr[0]);j++) fscanf(fp, "%d %lg %lg", &k, am+j, ap+j);
  }  
}

#ifdef READ_ON_THE_FLY
double accessf(FILE *fp, int flag_hy, int use_y_derivs, int ny, int nm, 
	       int ndy, int ndcb, double *x1, double *y1, double *y2, double x, double *cbv, double y)
{
  double f1iy, a, f11, f12, g11, g12, dy, xr, yr, bufm[N], bufp[N], am[N], ap[N];
  int nChterms;

  readf(fp, ny, &xr, &yr, bufm, bufp, 0);
  readf(fp, 1, x1, y1, am, ap, 1);
  nChterms = nfx(*x1);
    {
    getff(1,nChterms,nm,Idm[nm],ndy,ndcb,*x1,*y1,cbv,am,ap,&f11);
#ifdef VERBOSE
    printf("%d\t%lf\t%lf\t%lf\t%.11lg\n", nChterms, x, cbv[0], y, f11);
#endif
  }
  if(use_y_derivs) getff(1,nChterms,nm,Idm[nm],ndy+1,ndcb,*x1,*y1,cbv,am,ap,&g11); // get the derivative wrt y
  if(!flag_hy) {
    readf(fp, 1, x1, y2, am, ap, 1);
    {
      getff(1,nChterms,nm,Idm[nm],ndy,ndcb,*x1,*y2,cbv,am,ap,&f12);
#ifdef VERBOSE
      printf("%d\t%lf\t%lf\t%lf\t%.11lg\n", nChterms, x, cbv[0], y, f12);
#endif
    }
    if(use_y_derivs) {
      getff(1,nChterms,nm,Idm[nm],ndy+1,ndcb,*x1,y2[0],cbv,am,ap,&g12);
      dy = y1[0]-y2[0];
      f1iy = interpol3(y,y1[0],y2[0],dy,f11,f12,g11,g12);
    } else {
      a = (y2[0]-y)/(y2[0]-y1[0]);
      f1iy = a*f11 + (1.0-a)*f12; // value interpolated to target y, at x=x1[0];
    }
  } else f1iy = f11;

  return(f1iy);
}

#else 
// case where you have read in the weight functions upon initialization
// interpolates the form factor[nm] to the target point y using the grid
double accessv(int flag_hy, int use_y_derivs, int ix, int iy, int nm, 
	       int ndy, int ndcb, double *x1, double *y1, double *y2, double x, double *cbv, double y)
{
  double f1iy, a, f11, f12, g11, g12, dy, xr, yr, bufm[N], bufp[N], am[N], ap[N];
  double xtmp, xstp, ystp;
  int k, nmc, iy2, nx;

  //  readf(fp, ny, &xr, &yr, bufm, bufp, 0);
  //  readf(fp, 1, x1, y1, am, ap, 1);

  nmc = Iconv[nm];
  xstp = xstp_ff; // XX[0];
  ystp = ystp_ff; // YY[0];
  nx =  nfx(XX[ix]) ; // (XX[ix]>1.0 ?  nfx(XX[ix])-2 : nfx(XX[ix])); //    nfx(XX[ix]); // HERE 
  iy2 = iy+1;
  *x1 = XX[ix];
  *y1 = YY[iy];
  *y2 = YY[iy2];
  for(k=0;k<nx;k++) { am[k] = Ffm[nmc][ix][iy][k]; ap[k] = Ffp[nmc][ix][iy][k]; /*printf("# %d %lg %lg\n", k, am[k], ap[k]);*/}
  getff(1,nx,nm,Idm[nm],ndy,ndcb,*x1,*y1,cbv,am,ap,&f11);
#ifdef VERBOSE
  printf("%d\t%lf\t%lf\t%lf\t%.11lg\n", nx, x, cbv[0], y, f11);
#endif

  if(use_y_derivs) getff(1,nx,nm,Idm[nm],ndy+1,ndcb,*x1,*y1,cbv,am,ap,&g11); // get the derivative wrt y
  if(!flag_hy) {
    for(k=0;k<nx;k++) { am[k] = Ffm[nmc][ix][iy2][k]; ap[k] = Ffp[nmc][ix][iy2][k]; }
    getff(1,nx,nm,Idm[nm],ndy,ndcb,*x1,*y2,cbv,am,ap,&f12);
#ifdef VERBOSE
    printf("%d\t%lf\t%lf\t%lf\t%.11lg\n", nx, x, cbv[0], y, f12);
#endif
    if(use_y_derivs) {
      getff(1,nx,nm,Idm[nm],ndy+1,ndcb,*x1,y2[0],cbv,am,ap,&g12);
      dy = YY[iy]-YY[iy2];
      f1iy = interpol3(y,y1[0],y2[0],dy,f11,f12,g11,g12);
    } else {
      a = (YY[iy2]-y)/(YY[iy2]-YY[iy]);
       f1iy = a*f11 + (1.0-a)*f12; // value interpolated to target y, at x=x1[0];
     }
   } else f1iy = f11;

   return(f1iy);
 }
#endif 


 // #undef VERBOSE

#define SCALPROD(xv,yv) ((xv[0])*(yv[0])+(xv[1])*(yv[1])+(xv[2])*(yv[2])+(xv[3])*(yv[3]))
#define EPSIN2  1.0e-12
#define SETINVARIANTS {  xsq = SCALPROD(xv,xv); xsq = (xsq<EPSIN2 ? EPSIN2 : xsq);  x = sqrt(xsq); \
    ysq = SCALPROD(yv,yv);  ysq = (ysq<EPSIN2 ? EPSIN2 : ysq);  y = sqrt(ysq); \
     xdoty = SCALPROD(xv,yv); cb = xdoty/(x*y);}


void Tabd_xeq0(double yv[4], double ***tI, double ***tII, double ***tIII)
{
  double dg0dy[4], dg0dx[4], yhat[4], fx, fy, y, ysq, ell2, ell3a, dell2adx, dell2ady;
  int iy_tay, iy2_tay, alf, bet, dta, dlta[4][4];
  double ay, xtmp, dg1dx, dg1dy, dg2dx, ddg2dxdy, phi1, phi2, yad;

  ysq = SCALPROD(yv,yv);   y = sqrt(ysq); 

  iy_tay = y/Stp_tay;
  iy_tay = (iy_tay>=1 ? iy_tay-1 : iy_tay);
  iy2_tay=iy_tay+1;
  ay = (YY_tay[iy2_tay]-y)/(YY_tay[iy2_tay]-YY_tay[iy_tay]);

  fx = ay*G0dx_tay[iy_tay] + (1.0-ay)*G0dx_tay[iy2_tay];
  fy = (ay*G0dy_tay[iy_tay] + (1.0-ay)*G0dy_tay[iy2_tay]);
  for(bet=0;bet<4;bet++) {
    yhat[bet] = yv[bet]/y;
    dg0dy[bet] = yhat[bet]*fy;
    dg0dx[bet] = yhat[bet]*fx;
  }

  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++)   dlta[alf][dta] = 0.0;
  for(alf=0;alf<4;alf++)   dlta[alf][alf] = 1.0;

  ell2 = (ay*Gl2_tay[iy_tay] + (1.0-ay)*Gl2_tay[iy2_tay]);
  dell2adx = ay*Gl21_tay[iy_tay] + (1.0-ay)*Gl21_tay[iy2_tay];
  dell2ady = ay*Gl2dy_tay[iy_tay] + (1.0-ay)*Gl2dy_tay[iy2_tay];
  ell3a = ay*Gl3_tay[iy_tay] + (1.0-ay)*Gl3_tay[iy2_tay];

  //  ell2= dell2adx= dell2ady = ell3a = 0.0; 

  // (d/dxbeta+d/dybeta) T_{alpha delta}
  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      yad = (yv[alf]*yv[dta]-0.25*ysq*dlta[alf][dta]);
      for(bet=0;bet<4;bet++) {
	xtmp = yad*yhat[bet]*(dell2adx + dell2ady);
	xtmp += (dlta[alf][bet]*yv[dta]+dlta[bet][dta]*yv[alf]-0.5*yv[bet]*dlta[alf][dta])*(ell3a + ell2);
	tIII[alf][bet][dta] = xtmp + 0.25*dlta[alf][dta]*(dg0dx[bet]+dg0dy[bet]);
      }
    }

  // d/dxalpha T_{beta delta}
  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      for(bet=0;bet<4;bet++) {
	xtmp = (yv[bet]*yv[dta]-0.25*ysq*dlta[bet][dta])*yhat[alf]*(dell2adx );
	xtmp += (dlta[bet][alf]*yv[dta]+dlta[alf][dta]*yv[bet]-0.5*yv[alf]*dlta[bet][dta])*(ell3a);
	tII[alf][bet][dta] = xtmp + 0.25*dlta[bet][dta]*dg0dx[alf];
      }
    }

  dg2dx = ay*G21_tay[iy_tay] + (1.0-ay)*G21_tay[iy2_tay];
  ddg2dxdy = ay*G21dy_tay[iy_tay] + (1.0-ay)*G21dy_tay[iy2_tay];
  dg1dx = ay*G31A_tay[iy_tay]+(1.0-ay)*G31A_tay[iy2_tay]-(2.0/3)*(ay*G31B_tay[iy_tay]+(1.0-ay)*G31B_tay[iy2_tay])
	    -y*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]);
  dg1dy = ay*G3Ady_tay[iy_tay]+(1.0-ay)*G3Ady_tay[iy2_tay]-(ay*G3Bdy_tay[iy_tay]+(1.0-ay)*G3Bdy_tay[iy2_tay]);
  //  ddg2adxdx = 2.0*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]+(ay*G22B_tay[iy_tay]+(1.0-ay)*G22B_tay[iy2_tay])*(0.4+0.6*(2.0*cb*cb-1)));
  phi1 = 1.0*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]+(ay*G22B_tay[iy_tay]+(1.0-ay)*G22B_tay[iy2_tay])*(-0.2));
  phi2 = 1.0*((ay*G22B_tay[iy_tay]+(1.0-ay)*G22B_tay[iy2_tay])*(1.2));

  //  dg2dx = ddg2dxdy = phi1 = phi2 = 0.0; 
  //  dg1dx = dg1dy = 0.0;

  //   printf("## dg2dx= %.11lg ddg2dxdy= %.11lg  phi1= %.11lg  phi2= %.11lg\n", dg2dx, ddg2dxdy, phi1, phi2);
  //   printf("## dg1dx= %.11lg   dg1dy= %.11lg\n", dg1dx, dg1dy);

  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      for(bet=0;bet<4;bet++) {
	// first the d/dxalpha d/dybeta V_delta terms:  
        xtmp = dlta[alf][dta]*yhat[bet]*dg1dy;
	xtmp += (dlta[bet][dta]*yhat[alf]+yhat[dta]*dlta[alf][bet])*dg2dx;
        xtmp += yhat[bet]*yhat[alf]*yhat[dta]*(y*ddg2dxdy-dg2dx);
	// now the d/dxalpha d/dxbeta V_delta terms:
	xtmp += (dlta[alf][dta]*yhat[bet]+dlta[bet][dta]*yhat[alf])*dg1dx;
	xtmp += 2.0*yv[dta]*(dlta[alf][bet]*phi1 + yhat[alf]*yhat[bet]*phi2);
        tI[alf][bet][dta] = xtmp;
      }
    }

}


void Tabd_yeq0(double xv[4], double ***tI, double ***tII, double ***tIII)
{
  double dg0dy[4], dg0dx[4], xhat[4], fx, fy, x, xsq, xad, xa, xb;
  double ell1, ell3, dell1dx, dell1dycb, stpx;
  int ix_tay, ix2_tay, alf, bet, dta, dlta[4][4], nstpx, flag_hx;
  double ax, xtmp, dg1dx, dg1dycb, dg2dx, ddg1dxdx, ddg1dxdycb;

  xsq = SCALPROD(xv,xv);   x = sqrt(xsq); 

  nstpx = nstpx_ff;
  stpx = xstp_ff;
  ix_tay = (x<xb_ff ? x/stpx : ib_ff + (x-xb_ff)/xstp3_ff);
  if(x>XX[nstpx-1]) {
    printf("# Warning in Tabd_yeq0: x is higher than upper edge.\n");
    flag_hx = 1;
  } else flag_hx = 0;
  ix_tay = (ix_tay>=1 ? ix_tay-1 : ix_tay);
  ix2_tay=ix_tay+1;
  ax = (XX[ix2_tay]-x)/(XX[ix2_tay]-XX[ix_tay]);
  xa = XX[ix_tay];
  xb = XX[ix2_tay];
  //  printf("## ix_tay= %d ax= %lg xa= %lg xb= %lg\n", ix_tay, ax, xa, xb);
  
  fx = ax*alpha0dx_0p_taY[ix_tay] + (1.0-ax)*alpha0dx_0p_taY[ix2_tay];
  fy = 2.0*(ax*alpha0_1p_taY[ix_tay] + (1.0-ax)*alpha0_1p_taY[ix2_tay]);
  for(bet=0;bet<4;bet++) {
    xhat[bet] = xv[bet]/x;
    dg0dy[bet] = xhat[bet]*fy;
    dg0dx[bet] = xhat[bet]*fx;
  }

  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++)   dlta[alf][dta] = 0.0;
  for(alf=0;alf<4;alf++)   dlta[alf][alf] = 1.0;

  ell1 = (4.0/3.0)*(ax*alpha1_0p_taY[ix_tay]/(xa*xa*xa*xa)+(1.0-ax)*alpha1_0p_taY[ix2_tay]/(xb*xb*xb*xb));
  ell3 = ax*beta4_1p_taY[ix_tay]/(xa*xa)+(1.0-ax)*beta4_1p_taY[ix2_tay]/(xb*xb);
  //  dell1dx = -4.0/x*ell1 + (4.0/(3.0*xsq*xsq))*(ax*alpha1dx_0p_taY[ix_tay]+(1.0-ax)*alpha1dx_0p_taY[ix2_tay]);
  dell1dx = (4.0/3)*(ax*(-4.0*alpha1_0p_taY[ix_tay]/(xa*xa*xa*xa*xa)+alpha1dx_0p_taY[ix_tay]/(xa*xa*xa*xa))
		     +(1.0-ax)*(-4.0*alpha1_0p_taY[ix2_tay]/(xb*xb*xb*xb*xb)+alpha1dx_0p_taY[ix2_tay]/(xb*xb*xb*xb)));
  dell1dycb = 2.0*(ax*(4.0/3*alpha1_1p_taY[ix_tay]/(xa*xa*xa*xa)-beta4_1p_taY[ix_tay]/(xa*xa*xa))+
		    (1.0-ax)*(4.0/3*alpha1_1p_taY[ix2_tay]/(xb*xb*xb*xb)-beta4_1p_taY[ix2_tay]/(xb*xb*xb)));

    // (d/dxbeta+d/dybeta) T_{alpha delta}
  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      xad = (xv[alf]*xv[dta]-0.25*xsq*dlta[alf][dta]);
      for(bet=0;bet<4;bet++) {
	xtmp = xad*xhat[bet]*(dell1dx + dell1dycb);
	xtmp += (dlta[alf][bet]*xv[dta]+dlta[bet][dta]*xv[alf]-0.5*xv[bet]*dlta[alf][dta])*(ell1 + ell3);
	tIII[alf][bet][dta] = xtmp + 0.25*dlta[alf][dta]*(dg0dx[bet]+dg0dy[bet]);
      }
    }

  // d/dxalpha T_{beta delta}
  for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      for(bet=0;bet<4;bet++) {
	xtmp = (xv[bet]*xv[dta]-0.25*xsq*dlta[bet][dta])*xhat[alf]*(dell1dx);
	xtmp += (dlta[bet][alf]*xv[dta]+dlta[alf][dta]*xv[bet]-0.5*xv[alf]*dlta[bet][dta])*(ell1);
	tII[alf][bet][dta] = xtmp + 0.25*dlta[bet][dta]*dg0dx[alf];
      }
    }

  dg1dx = ax*alpha3dx_0p_taY[ix_tay] + (1.0-ax)*alpha3dx_0p_taY[ix2_tay];
  ddg1dxdx = ax*alpha3dxdx_0p_taY[ix_tay] + (1.0-ax)*alpha3dxdx_0p_taY[ix2_tay];
  dg2dx = 2.0*(ax*beta2dx_1p_taY[ix_tay] + (1.0-ax)*beta2dx_1p_taY[ix2_tay]);
  dg1dycb = 2.0*(ax*(-beta2_1p_taY[ix_tay]/xa + alpha3_1p_taY[ix_tay]) + (1.0-ax)*(-beta2_1p_taY[ix2_tay]/xb + alpha3_1p_taY[ix2_tay]));
  ddg1dxdycb = 2.0*(ax*(1.0/(xa*xa)*beta2_1p_taY[ix_tay]-1.0/xa*beta2dx_1p_taY[ix_tay]+alpha3dx_1p_taY[ix_tay])
		    +(1.0-ax)*(1.0/(xb*xb)*beta2_1p_taY[ix2_tay]-1.0/xb*beta2dx_1p_taY[ix2_tay]+alpha3dx_1p_taY[ix2_tay]));

  
  //  printf("ddg1dxdx= %lg  discrete deriv= %lg\n", ddg1dxdx, (alpha3dx_0p_taY[ix2_tay]-alpha3dx_0p_taY[ix_tay])/(XX[ix2_tay]-XX[ix_tay]));
  //  printf("ddg1dxdycb= %lg discrete deriv= %lg\n", ddg1dxdycb, 
  //	 ((beta2_1p_taY[ix2_tay]/xb + alpha3_1p_taY[ix2_tay]) - (beta2_1p_taY[ix_tay]/xa + alpha3_1p_taY[ix_tay]))*2.0/(XX[ix2_tay]-XX[ix_tay]));
  //  printf("dg1dx= %lg  dg2dx= %lg  dg1dycb= %lg\n", dg1dx, dg2dx, dg1dycb);
  //  dg1dx = ddg1dxdx = dg1dycb = ddg1dxdycb = 0.0; 
  //  dg2dx = 0.0; 

    for(alf=0;alf<4;alf++) for(dta=0;dta<4;dta++) {
      for(bet=0;bet<4;bet++) {
	// first the d/dxalpha d/dybeta V_delta terms:  
        xtmp = (dlta[alf][dta]*xhat[bet]+dlta[bet][dta]*xhat[alf]+dlta[alf][bet]*xhat[dta]-xhat[alf]*xhat[bet]*xhat[dta])*dg1dx;
	xtmp += (xhat[alf]*xhat[bet]*xv[dta])*ddg1dxdx;
	xtmp += (dlta[alf][dta]*xhat[bet]+dlta[alf][bet]*xhat[dta]-xhat[alf]*xhat[bet]*xhat[dta])*(dg1dycb);
	// now the d/dxalpha d/dxbeta V_delta terms:
	xtmp += xhat[alf]*xhat[bet]*xv[dta]*ddg1dxdycb;
	xtmp += dlta[bet][dta]*xhat[alf]*dg2dx;
        tI[alf][bet][dta] = xtmp;
      }
    }

}

 // return dyv[bet]= \partial^(y)_\beta < I >_epsilon
 // return dxv[bet]= \partial^(x)_\beta < I >_epsilon
 void chnr_dS(double xv[4], double yv[4], double *dxv, double *dyv)
 {
   int bet;
   double x, cb, y, yhat, c2, xhat, c4;
   double dg0dy, dg0dx, dg0dcb;
   double xsq, ysq, xdoty;
#ifdef TAYLORX
   int ix, iy, iy2_tay, iy_tay, flag_hy, use_y_derivs, nm, ndy, ndcb, nstpy;
   double ax, ay, f1, f2, x1, x2, y1, y2, stpy;
#endif
#ifdef XMYSWAP_S
   int flag, flag2;
   double rx, rxy, xmysq, xmy, cborig, yorig, ysqorig, dgs0dy, dgs0dx, dgs0dcb, cd, ca;
#endif 
   SETINVARIANTS;

#ifdef TESTING
   double vtmp[10];
   testff(1, x, cb, y, &c2, vtmp);
   dg0dx = vtmp[0];
   dg0dy = vtmp[2];
   dg0dcb = vtmp[1];
#else 
#ifdef XMYSWAP_S
   xmysq=(xsq+ysq)*(1.00000000000001)-2.0*xdoty;
   xmy=sqrt(xmysq);
   rx = (x>xmy ? xmy/x : x/xmy);
   rxy = (x>y ? y/x : x/y);
   if(rx<rxy && xmy<YY[nstpy_ff-1])  { // && xmy>=0.0436
     // swap y for x-y
     yorig = y;
     cborig = cb;
     ysqorig = ysq;
     cb = (x-y*cb)/xmy;
     y = xmy;
     flag2=1;
   } else flag2=0;
#endif
#ifdef TAYLORX 
   if(x<XX[0]) {
     nstpy = nstpy_ff;
     stpy = YY[1]-YY[0];
     iy_tay = y/Stp_tay;
     iy_tay = (iy_tay>=1 ? iy_tay-1 : iy_tay);
     iy2_tay=iy_tay+1;
     ay = (YY_tay[iy2_tay]-y)/(YY_tay[iy2_tay]-YY_tay[iy_tay]);
     ax = (XX[0]-x)/(XX[0]);    
     //     for(iy=0;iy<nstpy;iy++) if(YY[iy]>y) break; 
     iy = y/stpy;
     iy = (iy>=1 ? iy-1 : iy);
     flag_hy=0;
     ix=0;
     //     printf("iy_tay= %d  iy2_tay= %d\n", iy_tay, iy2_tay);
     nm=10; ndy=0; ndcb=0;
     use_y_derivs=1;
     f1 = (ay*G0dx_tay[iy_tay] + (1.0-ay)*G0dx_tay[iy2_tay])*cb;
     f2 = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg0dx = ax*f1 + (1.0-ax)*f2;

     nm=0; ndy=1; ndcb=0;
     use_y_derivs=0;
     f1 = (ay*G0dy_tay[iy_tay] + (1.0-ay)*G0dy_tay[iy2_tay]);
     f2 = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg0dy = ax*f1 + (1.0-ax)*f2;

     nm=0; ndy=0; ndcb=1;
     use_y_derivs=1;
     f1 = 0.0;
     f2 = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg0dcb = ax*f1 + (1.0-ax)*f2;     
  } else 
#endif
  {
    dg0dx = extractff(10, 0, 0, x, cb, y);
    dg0dy = extractff(0, 1, 0, x, cb, y);
    dg0dcb = extractff(0, 0, 1, x, cb, y);
  }
#ifdef XMYSWAP_S
   if(flag2) {
     dgs0dx  = dg0dx;
     dgs0dcb = dg0dcb;
     dgs0dy  = dg0dy;
     ca = cb; 
     y = yorig;
     cb = cborig;
     ysq = ysqorig;
     cd = (y-x*cb)/xmy;
     dg0dx = dgs0dx + (1.0-ca*ca)/xmy*dgs0dcb + ca * dgs0dy;
     dg0dcb = -(ysq*cd*dgs0dcb + x*y*xmy*dgs0dy)/xmysq; 
     dg0dy = -(cb+ca*cd)/xmy*dgs0dcb+cd*dgs0dy;
   } 
   //printf("flag2= %d YY[nstpy_ff-1]= %lg  xmy= %lg rx= %lg  rxy= %lg\n", flag2, YY[nstpy_ff-1], xmy, rx, rxy);
   //printf("extractff(0, 0, 0, x, cb, y)= %lg   extractff(0, 0, 0, x, ca, xmy)= %lg\n", extractff(0, 0, 0, x, cb, y), extractff(0, 0, 0, x, ca, xmy));
#endif
   //   printf("x= %.11lg   cb= %.11lg   y= %.11lg\n", x, cb, y);
   //   printf("dg0dx= %.11lg   dg0dcb= %.11lg   dg0dy= %.11lg\n", dg0dx, dg0dcb, dg0dy); // exit(0);
#endif 

  for(bet=0;bet<=3;bet++)  {
    yhat = yv[bet]/y;
    xhat = xv[bet]/x;
    c2 = (xhat-cb*yhat)/y;
    c4 = (yhat-cb*xhat)/x;
    dyv[bet] = yhat*dg0dy+c2*dg0dcb;
    dxv[bet] = xhat*dg0dx+c4*dg0dcb;
  }
}

// return dyv[alf][bet][dta]= \partial^(y)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
// and    dxv[alf][bet][dta]= \partial^(x)_\beta <(epsilon_alf epsilon_dta - 1/4 delta_{alf dta}) I>_epsilon
void chnr_dT(double xv[4], double yv[4], double ***dxv, double ***dyv)
{
  int alf, bet, dta, d_ab, d_bd, d_ad;
  double x, cb, y, xsq, ysq, xdoty;
  double t1, t2, t3, t4, t5, c1, c2, c3, c4;
  double ell1, ell2, ell3, v1, ell4, dv1dx, dv1dy, dv1dcb, dell2dx, dell2dy, dell2dcb;
  double dell1dx, dell1dy, dell1dcb, dell3dx, dell3dy, dell3dcb, dell4dy, dell4dx, dell4dcb;
#ifdef TAYLORX
  int ix, iy, iy2_tay, iy_tay, flag_hy, use_y_derivs, nm, ndy, ndcb, nstpy;
  double x1, x2, y1, y2, ax, ay, xb, xbsq, ell2a, ell2b, dell2adx, dell2bdx, dell2ady, dell2bdy, dell2adcb, dell2bdcb, stpy;
  double v1b, dv1bdx, dv1bdy, dv1bdcb, ell4b, dell4bdx, dell4bdy, dell4bdcb;
  double ell3a, ell3b, ell1b, dell3ady, dell3bdy, dell3adcb, dell3bdcb, dell3bdx;
#endif
#ifdef XMYSWAP_T
   int flag, flag2;
   double rx, rxy, xmysq, xmy, cborig, yorig, ysqorig, ell2s, dells2dy, dells2dx, dells2dcb, cd, ca;
   double ell1s, dells1dy, dells1dx, dells1dcb, ell3s, dells3dy, dells3dx, dells3dcb;
#endif 

  SETINVARIANTS;

#ifdef TESTING
  double vtmp[10];
  testff2(1, x, cb, y, &c2, vtmp);
#if 1
  ell1 = vtmp[0];
  dell1dx = vtmp[1];
  dell1dcb = vtmp[2];
  dell1dy = vtmp[3];
  ell2 = 0;
  dell2dx = 0;
  dell2dcb = 0;
  dell2dy = 0;
  ell3 = 0;
  dell3dx = 0;
  dell3dcb = 0;
  dell3dy = 0;
#endif
#if 0
  ell2 = vtmp[0];
  dell2dx = vtmp[1];
  dell2dcb = vtmp[2];
  dell2dy = vtmp[3];
  ell1 = 0;
  dell1dx = 0;
  dell1dcb = 0;
  dell1dy = 0;
  ell3 = 0;
  dell3dx = 0;
  dell3dcb = 0;
  dell3dy = 0;
#endif
#if 0
  ell3 = vtmp[0];
  dell3dx = vtmp[1];
  dell3dcb = vtmp[2];
  dell3dy = vtmp[3];
  ell1 = 0;
  dell1dx = 0;
  dell1dcb = 0;
  dell1dy = 0;
  ell2 = 0;
  dell2dx = 0;
  dell2dcb = 0;
  dell2dy = 0;
#endif
#else
#ifdef XMYSWAP_T
   xmysq=(xsq+ysq)*(1.00000000000001)-2.0*xdoty;
   xmy=sqrt(xmysq);
   rx = (x>xmy ? xmy/x : x/xmy);
   rxy = (x>y ? y/x : x/y);
   if(rx<rxy && xmy<YY[nstpy_ff-1])  { //  && xmy>=0.0436
     // swap y for x-y
     yorig = y;
     cborig = cb;
     ysqorig = ysq;
     cb = (x-y*cb)/xmy;
     y = xmy;
     ysq = xmysq;
     flag2=1;
   } else flag2=0;
#endif
#ifdef TAYLORX 
   if(x<XX[0]) {
     nstpy = nstpy_ff;
     stpy = YY[1]-YY[0];
     stpy = YY[1]-YY[0];
     iy_tay = y/Stp_tay;
     iy_tay = (iy_tay>=1 ? iy_tay-1 : iy_tay);
     iy2_tay=iy_tay+1;
     ay = (YY_tay[iy2_tay]-y)/(YY_tay[iy2_tay]-YY_tay[iy_tay]);
     ax = (XX[0]-x)/(XX[0]);    
     //     for(iy=0;iy<nstpy;iy++) if(YY[iy]>y) break; 
     iy = y/stpy;
     iy = (iy>=1 ? iy-1 : iy);
     flag_hy=0;
     ix=0;
     xb=XX[0];
     xbsq=xb*xb;
     //     printf("iy_tay= %d  iy2_tay= %d\n", iy_tay, iy2_tay);
     nm=5; ndy=0; ndcb=0;
     use_y_derivs=1;
     ell2a = (ay*Gl2_tay[iy_tay] + (1.0-ay)*Gl2_tay[iy2_tay]);
     ell2b = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ell2 = ax*ell2a + (1.0-ax)*ell2b;

     nm=15; ndy=0; ndcb=0;
     use_y_derivs=1;
     dell2adx = (ay*Gl21_tay[iy_tay] + (1.0-ay)*Gl21_tay[iy2_tay])*cb;
     dell2bdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dell2dx = ax*dell2adx + (1.0-ax)*dell2bdx;

     nm=5; ndy=1; ndcb=0;
     use_y_derivs=0;
     dell2ady = (ay*Gl2dy_tay[iy_tay] + (1.0-ay)*Gl2dy_tay[iy2_tay]);
     dell2bdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dell2dy = ax*dell2ady + (1.0-ax)*dell2bdy;

     nm=5; ndy=0; ndcb=1;
     use_y_derivs=1;
     dell2adcb = 0.0;
     dell2bdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dell2dcb = ax*dell2adcb + (1.0-ax)*dell2bdcb;

     nm=1; ndy=0; ndcb=0;
     use_y_derivs=1;
     v1b = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=11; ndy=0; ndcb=0;
     use_y_derivs=1;
     dv1bdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=1; ndy=1; ndcb=0;
     use_y_derivs=0;
     dv1bdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=1; ndy=0; ndcb=1;
     use_y_derivs=1;
     dv1bdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);

     nm=4; ndy=0; ndcb=0;
     use_y_derivs=1;
     ell4b = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=14; ndy=0; ndcb=0;
     use_y_derivs=1;
     dell4bdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=4; ndy=1; ndcb=0;
     use_y_derivs=0;
     dell4bdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=4; ndy=0; ndcb=1;
     use_y_derivs=1;
     dell4bdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);

     ell3a = ay*Gl3_tay[iy_tay] + (1.0-ay)*Gl3_tay[iy2_tay];
     ell3b = ell4b/(2.0*xbsq*ysq) - y*cb*ell2b/xb;
     ell3 = ax*ell3a + (1.0-ax)*ell3b;

     ell1b = 4.0/(3.0*xbsq*xbsq)*(v1b - xbsq*ysq*(cb*cb-0.25)*ell2b - 1.5*xbsq*xb*y*cb*ell3b);
     ell1 = ell1b;

     dell3ady = ay*Gl3dy_tay[iy_tay] + (1.0-ay)*Gl3dy_tay[iy2_tay];
     dell3bdy = (-ell4b/(ysq*y) + dell4bdy/(2.0*ysq) - cb*xb*ell2b - y*xb*cb*dell2bdy)/(xbsq);
     dell3dy = ax*dell3ady + (1.0-ax)*dell3bdy;

     dell3adcb = 0.0;
     dell3bdcb = (dell4bdcb/(2.0*xb*ysq) - y*ell2b - y*cb*dell2bdcb)/xb;
     dell3dcb = ax*dell3adcb + (1.0-ax)*dell3bdcb;

     dell3bdx = (-ell4b/(xb*ysq) + dell4bdx/(2.0*ysq) + y*cb*ell2b - y*xb*cb*dell2bdx)/xbsq;
     dell3dx = dell3bdx;

     dell1dy = 4.0/(3.0*xbsq*xbsq)*(dv1bdy-2.0*xbsq*y*(cb*cb-0.25)*ell2b-xbsq*ysq*(cb*cb-0.25)*dell2bdy
				  -1.5*xbsq*xb*cb*ell3b-1.5*xbsq*xb*y*cb*dell3bdy);
     dell1dcb = 4.0/(3.0*xbsq*xbsq)*(dv1bdcb - 2.0*xbsq*ysq*cb*ell2b-xbsq*ysq*(cb*cb-0.25)*dell2bdcb
				   -1.5*xbsq*xb*y*ell3b-1.5*xbsq*xb*y*cb*dell3bdcb);
     dell1dx = -4.0*ell1b/xb + 4.0/(3.0*xbsq*xbsq)*(dv1bdx - 2.0*xb*ysq*(cb*cb-0.25)*ell2b-xbsq*ysq*(cb*cb-0.25)*dell2bdx
						-4.5*xbsq*y*cb*ell3b - 1.5*xbsq*xb*y*cb*dell3bdx);
  } else 
#endif
     { 
  ell2 = extractff(5, 0, 0, x, cb, y);
  ell4 = extractff(4, 0, 0, x, cb, y);
  v1 = extractff(1, 0, 0, x, cb, y);
  ell3 = ell4/(2.0*xsq*ysq) - y*cb*ell2/x;
  ell1 = 4.0/(3.0*xsq*xsq)*(v1 - xsq*ysq*(cb*cb-0.25)*ell2 - 1.5*xsq*x*y*cb*ell3);

  dv1dx = extractff(11, 0, 0, x, cb, y);
  dv1dy = extractff(1, 1, 0, x, cb, y); 
  dv1dcb = extractff(1, 0, 1, x, cb, y);
  dell2dx = extractff(15, 0, 0, x, cb, y);
  dell2dy = extractff(5, 1, 0, x, cb, y);
  dell2dcb = extractff(5, 0, 1, x, cb, y);
  dell4dy = extractff(4, 1, 0, x, cb, y);
  dell4dcb = extractff(4, 0, 1, x, cb, y);
  dell4dx = extractff(14, 0, 0, x, cb, y);


  dell3dy = (-ell4/(ysq*y) + dell4dy/(2.0*ysq) - cb*x*ell2 - y*x*cb*dell2dy)/(xsq);
  dell3dcb = (dell4dcb/(2.0*x*ysq) - y*ell2 - y*cb*dell2dcb)/x;
  dell3dx = (-ell4/(x*ysq) + dell4dx/(2.0*ysq) + y*cb*ell2 - y*x*cb*dell2dx)/xsq;
  dell1dy = 4.0/(3.0*xsq*xsq)*(dv1dy-2.0*xsq*y*(cb*cb-0.25)*ell2-xsq*ysq*(cb*cb-0.25)*dell2dy
			       -1.5*xsq*x*cb*ell3-1.5*xsq*x*y*cb*dell3dy);
  dell1dcb = 4.0/(3.0*xsq*xsq)*(dv1dcb - 2.0*xsq*ysq*cb*ell2-xsq*ysq*(cb*cb-0.25)*dell2dcb
				-1.5*xsq*x*y*ell3-1.5*xsq*x*y*cb*dell3dcb);
  dell1dx = -4.0*ell1/x + 4.0/(3.0*xsq*xsq)*(dv1dx - 2.0*x*ysq*(cb*cb-0.25)*ell2-xsq*ysq*(cb*cb-0.25)*dell2dx
					     -4.5*xsq*y*cb*ell3 - 1.5*xsq*x*y*cb*dell3dx);
     }

#ifdef XMYSWAP_T
   if(flag2) {
     ell2s = ell2;
     dells2dx  = dell2dx;
     dells2dcb = dell2dcb;
     dells2dy  = dell2dy;
     ca = cb; 
     y = yorig;
     cb = cborig;
     ysq = ysqorig;
     cd = (y-x*cb)/xmy;
     ell2 = ell2s;
     dell2dx = dells2dx + (1.0-ca*ca)/xmy*dells2dcb + ca * dells2dy;
     dell2dcb = -(ysq*cd*dells2dcb + x*y*xmy*dells2dy)/xmysq;
     dell2dy = -(cb+ca*cd)/xmy*dells2dcb+cd*dells2dy;
     ell3s = -ell3;
     dells3dx  = -dell3dx;
     dells3dcb = -dell3dcb;
     dells3dy  = -dell3dy;
     ell3 =  ell3s - ell2; 
     dell3dx =  dells3dx + (1.0-ca*ca)/xmy*dells3dcb + ca * dells3dy - dell2dx;
     dell3dcb =  -(ysq*cd*dells3dcb + x*y*xmy*dells3dy)/xmysq - dell2dcb;
     dell3dy =  -(cb+ca*cd)/xmy*dells3dcb+cd*dells3dy - dell2dy;
     ell1s = ell1;
     dells1dx  = dell1dx;
     dells1dcb = dell1dcb;
     dells1dy  = dell1dy;
     ell1 = ell1s - ell2 -2.0*ell3;
     dell1dx = dells1dx + (1.0-ca*ca)/xmy*dells1dcb + ca * dells1dy - dell2dx - 2.0*dell3dx;
     dell1dcb = -(ysq*cd*dells1dcb + x*y*xmy*dells1dy)/xmysq - dell2dcb -2.0*dell3dcb;
     dell1dy = -(cb+ca*cd)/xmy*dells1dcb+cd*dells1dy - dell2dy -2.0*dell3dy;
   } 
   //   printf("flag2= %d YY[nstpy_ff-1]= %lg  xmy= %lg rx= %lg  rxy= %lg\n", flag2, YY[nstpy_ff-1], xmy, rx, rxy);
#endif

  //  ell2=dell2dx=dell2dcb=dell2dy=0.0;  ell3=dell3dx=dell3dcb=dell3dy=0.0; 
   //      printf("x= %.11lg   cb= %.11lg   y= %.11lg\n", x, cb, y);
   //     printf("ell1= %.11lg  dell1dx= %.11lg  dell1dcb= %.11lg  dell1dy= %.11lg\n", ell1, dell1dx, dell1dcb, dell1dy);
   //     printf("ell2= %.11lg  dell2dx= %.11lg  dell2dcb= %.11lg  dell2dy= %.11lg\n", ell2, dell2dx, dell2dcb, dell2dy);
   //     printf("ell3= %.11lg  dell3dx= %.11lg  dell3dcb= %.11lg  dell3dy= %.11lg\n", ell3, dell3dx, dell3dcb, dell3dy); exit(0);
#endif

  for(bet=0;bet<=3;bet++)  {  
    c1 = yv[bet]/y;
    c3 = xv[bet]/x;
    c2 = (c3-cb*c1)/y;
    c4 = (c1-cb*c3)/x;
    for(alf=0;alf<=3;alf++)  {
      d_ab = (alf==bet ? 1 : 0);
      for(dta=0;dta<=3;dta++) {
	d_bd = (bet==dta ? 1 : 0);
	d_ad = (alf==dta ? 1 : 0);
	t1 = d_ab*yv[dta]+d_bd*yv[alf]-0.5*d_ad*yv[bet];
	t2 = d_bd*xv[alf]+d_ab*xv[dta]-0.5*d_ad*xv[bet];
	t3 = xv[alf]*xv[dta]-0.25*xsq*d_ad;
	t4 = yv[alf]*yv[dta]-0.25*ysq*d_ad;
	t5 = xv[alf]*yv[dta]+yv[alf]*xv[dta]-0.5*xdoty*d_ad;
	dyv[alf][bet][dta] = t1*ell2 + t2*ell3 + t3*(c1*dell1dy + c2*dell1dcb) 
	  + t4*(c1*dell2dy + c2*dell2dcb) + t5*(c1*dell3dy + c2*dell3dcb);
	dxv[alf][bet][dta] = t2*ell1 + t1*ell3 + t3*(c3*dell1dx+c4*dell1dcb)
	  + t4*(c3*dell2dx+c4*dell2dcb) + t5*(c3*dell3dx+c4*dell3dcb);
      }
    }
  }
  /*
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][0][1]= %.11lg\n", x, cb, y, dxv[0][0][1]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][1][1]= %.11lg\n", x, cb, y, dxv[0][1][1]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][2][1]= %.11lg\n", x, cb, y, dxv[0][2][1]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][3][1]= %.11lg\n", x, cb, y, dxv[0][3][1]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][0][0]= %.11lg\n", x, cb, y, dxv[0][0][0]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][1][0]= %.11lg\n", x, cb, y, dxv[0][1][0]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][2][0]= %.11lg\n", x, cb, y, dxv[0][2][0]);
  printf("# x= %lg cb= %lg  y= %lg  dxv[0][3][0]= %.11lg\n", x, cb, y, dxv[0][3][0]);
  */
}

// returns dv[alf][bet][dta] = \partial^(x)_\alpha (\partial^(x)_\beta + \partial^(y)_\beta) < \epsilon_\delta I>_\epsilon
void chnr_dV(double xv[4], double yv[4], double ***dv)
{
  int alf, bet, dta, d_ab, d_ad, d_bd;
  double xhat[4], c2v[4], yhat[4], c4v[4];
  double x, cb, y, xsq, ysq, xdoty, f;
  double g2, dg1dx, dg1dy, dg1dcb, dg2dx, dg2dy, dg2dcb, dg3dx, dg3dy, dg3dcb;
  double ddg2dxdx, ddg2dxdy, ddg2dxdcb, ddg2dydcb, ddg2dcbdcb;
  double ddg3dxdx, ddg3dxdy, ddg3dxdcb, ddg3dydcb, ddg3dcbdcb;
  double ddg1dxdx, ddg1dxdy, ddg1dxdcb, ddg1dydcb, ddg1dcbdcb;	
  double stpy;
#ifdef TAYLORX
  int ix, iy, iy2_tay, iy_tay, flag_hy, use_y_derivs, nm, ndy, ndcb, nstpy;
  double x1, x2, y1, y2, ax, ay, xb, xbsq, g2a, g2b, dg2adx, dg2bdx, dg2ady, dg2bdy, dg2adcb, dg2bdcb;
  double ddg2adxdx, ddg2bdxdx, ddg2adxdy, ddg2bdxdy, ddg2adxdcb, ddg2bdxdcb, ddg2adydcb, ddg2bdydcb, ddg2adcbdcb, ddg2bdcbdcb;
  double dg1adx, dg1bdx, dg1ady, dg1bdy, dg1adcb, dg1bdcb;
  double ddg1adxdx, ddg1bdxdx, ddg1adxdy, ddg1bdxdy, ddg1adxdcb, ddg1bdxdcb, ddg1adydcb, ddg1bdydcb, ddg1adcbdcb, ddg1bdcbdcb;
  double dg3bdx, dg3bdy, dg3bdcb;
  double ddg3bdxdx, ddg3bdxdy, ddg3bdxdcb, ddg3bdydcb, ddg3bdcbdcb;
#endif
#ifdef XMYSWAP_V
  int flag, flag2, iy1;
   double rx, rxy, xmysq, xmy, cborig, yorig, ysqorig, g2s, dgs2dy, dgs2dx, dgs2dcb, cd, ca;
   double g1s, dgs1dy, dgs1dx, dgs1dcb;
   double ya, yb, ddg2dydy, ddg3dydy, ddg1dydy, fq1, fq2, ddgs2dxdx, ddgs2dxdcb, ddgs2dxdy, ddgs2dydcb, ddgs2dcbdcb, ddgs2dydy;
   double ddgs1dxdx, ddgs1dxdcb, ddgs1dxdy, ddgs1dydcb, ddgs1dcbdcb, ddgs1dydy;
#ifdef TAYLORX
   double ddg2adydy, ddg2bdydy, ddg1adydy, ddg1bdydy, dg2bdy_ya, dg2bdy_yb, dg3bdy_ya, dg3bdy_yb, ddg3bdydy;
#endif 
#endif 

  SETINVARIANTS;

#ifdef TESTING
  double vtmp[10], gtmp1, gtmp2, gtmp3, gtmp4;
#if (1==0)
  g2=0.0;
  dg2dx=0.0;
  dg2dy=0.0;
  dg2dcb=0.0;
  ddg2dxdx=0.0;
  ddg2dxdy=0.0;
  ddg2dxdcb=0.0;
  ddg2dydcb=0.0;
  ddg2dcbdcb=0.0;
  testff(1, x, cb, y, &f, vtmp);
  dg1dx = vtmp[0];
  dg1dcb = vtmp[1];
  dg1dy = vtmp[2];
  ddg1dxdx = vtmp[3];
  ddg1dxdcb = vtmp[4];
  ddg1dxdy = vtmp[5];
  ddg1dydcb = vtmp[6];
  ddg1dcbdcb = vtmp[7];
#else
  dg1dx=0.0;
  dg1dy=0.0;
  dg1dcb=0.0;
  ddg1dxdx=0.0;
  ddg1dxdy=0.0;
  ddg1dxdcb=0.0;
  ddg1dydcb=0.0;
  ddg1dcbdcb=0.0;
  testff(1, x, cb, y, &f, vtmp);
  dg2dx = vtmp[0];
  dg2dcb = vtmp[1];
  dg2dy = vtmp[2];
  ddg2dxdx = vtmp[3];
  ddg2dxdcb = vtmp[4];
  ddg2dxdy = vtmp[5];
  ddg2dydcb = vtmp[6];
  ddg2dcbdcb = vtmp[7];
  //  printf("dg2dx= %lf  dg2dcb= %lf dg2dy = %lf ddg2dxdx= %lf\n", dg2dx, dg2dcb, dg2dy, ddg2dxdx);
  //  printf("dg1dcb= %lf  dg2dcb= %lf dg1dx = %lf\n", dg1dcb, dg2dcb, dg1dx);
#endif 
#else
#ifdef XMYSWAP_V
   xmysq=(xsq+ysq)*(1.00000000000001)-2.0*xdoty;
   xmy=sqrt(xmysq);
   rx = (x>xmy ? xmy/x : x/xmy);
   rxy = (x>y ? y/x : x/y);
   if(rx<rxy && xmy<YY[nstpy_ff-2])  { // && xmy>=0.0436
     // swap y for x-y
     yorig = y;
     cborig = cb;
     ysqorig = ysq;
     cb = (x-y*cb)/xmy;
     y = xmy;
     flag2=1;
     //     printf("flag2= %d : SWAP!\n", flag2);
   } else { flag2=0;   /* printf("flag2= %d : no swap.\n", flag2); */}

#endif
#ifdef TAYLORX 
   if(x<XX[0]) {
     //     printf("Taylor expansion!\n");
     nstpy = nstpy_ff;
     stpy = YY[1]-YY[0];
     iy_tay = y/Stp_tay;
     iy_tay = (iy_tay>=1 ? iy_tay-1 : iy_tay);
     iy2_tay=iy_tay+1;
     ay = (YY_tay[iy2_tay]-y)/(YY_tay[iy2_tay]-YY_tay[iy_tay]);
     ax = (XX[0]-x)/(XX[0]);    
     //     for(iy=0;iy<nstpy;iy++) if(YY[iy]>y) break; 
     iy = y/stpy;
     iy = (iy>=1 ? iy-1 : iy);
     flag_hy=0;
     ix=0;
     xb=XX[0];
     xbsq=xb*xb;

     nm=2; ndy=0; ndcb=0;
     use_y_derivs=1;
     g2a = 0.0; // (ay*G2_tay[iy_tay] + (1.0-ay)*G2_tay[iy2_tay]);
     g2b = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     g2 = ax*g2a + (1.0-ax)*g2b;

     nm=12; ndy=0; ndcb=1;
     use_y_derivs=1;
     ddg2adxdcb = (ay*G21_tay[iy_tay] + (1.0-ay)*G21_tay[iy2_tay]);
     ddg2bdxdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ddg2dxdcb = ax*ddg2adxdcb + (1.0-ax)*ddg2bdxdcb;
     //printf("ddg2adxdcb= %.11lg  ddg2bdxdcb= %.11lg\n", ddg2adxdcb, ddg2bdxdcb);

     nm=12; ndy=0; ndcb=0;
     use_y_derivs=1;
     dg2adx = ddg2adxdcb*cb;
     dg2bdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg2dx = ddg2dxdcb*cb; // ax*dg2adx + (1.0-ax)*dg2bdx;
     //printf("dg2adx= %.11lg  dg2bdx= %.11lg\n", dg2adx, dg2bdx); 

     nm=2; ndy=0; ndcb=1;
     use_y_derivs=1;
     dg2adcb = 0.0; // (ay*G21_tay[iy_tay] + (1.0-ay)*G21_tay[iy2_tay])*x; 
     dg2bdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg2dcb = ddg2dxdcb*x; // dg2dx*x/cb; // ax*dg2adcb + (1.0-ax)*dg2bdcb;
     //printf("dg2adcb= %.11lg  dg2bdcb= %.11lg\n", dg2adcb, dg2bdcb); 

     nm=2; ndy=1; ndcb=0;
     use_y_derivs=0;
     dg2ady = 0.0; // (ay*G2dy_tay[iy_tay] + (1.0-ay)*G2dy_tay[iy2_tay])*cb;
     dg2bdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     dg2dy = ax*dg2ady + (1.0-ax)*dg2bdy;
     //printf("dg2ady= %.11lg  dg2bdy= %.11lg\n", dg2ady, dg2bdy); 

     nm=12; ndy=1; ndcb=0;
     use_y_derivs=0;
     ddg2adxdy = (ay*G21dy_tay[iy_tay] + (1.0-ay)*G21dy_tay[iy2_tay])*cb;
     ddg2bdxdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ddg2dxdy = ax*ddg2adxdy + (1.0-ax)*ddg2bdxdy;
     //printf("ddg2adxdy= %.11lg  ddg2bdxdy= %.11lg\n", ddg2adxdy, ddg2bdxdy);

     nm=22; ndy=0; ndcb=0;
     use_y_derivs=1;
     ddg2adxdx = 2.0*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]+(ay*G22B_tay[iy_tay]+(1.0-ay)*G22B_tay[iy2_tay])*(0.4+0.6*(2.0*cb*cb-1)));
     ddg2bdxdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ddg2dxdx = ax*ddg2adxdx + (1.0-ax)*ddg2bdxdx;
     //printf("ddg2adxdx= %.11lg  ddg2bdxdx= %.11lg\n", ddg2adxdx, ddg2bdxdx);

     nm=2; ndy=0; ndcb=2;
     use_y_derivs=1;
     ddg2adcbdcb = 0.0;
     ddg2bdcbdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ddg2dcbdcb = ddg2bdcbdcb*xsq/xbsq; // the function starts quadratically // ax*dg2adcbdcb + (1.0-ax)*dg2bdcbdcb;
     //printf("ddg2adcbdcb= %.11lg  ddg2bdcbdcb= %.11lg\n", ddg2adcbdcb, ddg2bdcbdcb);

     nm=2; ndy=1; ndcb=1;
     use_y_derivs=0;
     ddg2adydcb = 0.0;
     ddg2bdydcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     ddg2dydcb = ax*ddg2adydcb + (1.0-ax)*ddg2bdydcb;
     //printf("ddg2adydcb= %.11lg  ddg2bdydcb= %.11lg\n", ddg2adydcb, ddg2bdydcb);

     nm=13; ndy=0; ndcb=0;
     use_y_derivs=1;
     dg3bdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=3; ndy=1; ndcb=0;
     use_y_derivs=0;
     dg3bdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=3; ndy=0; ndcb=1;
     use_y_derivs=1;
     dg3bdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=13; ndy=1; ndcb=0;
     use_y_derivs=0;
     ddg3bdxdy = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=13; ndy=0; ndcb=1;
     use_y_derivs=1;
     ddg3bdxdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=23; ndy=0; ndcb=0;
     use_y_derivs=1;
     ddg3bdxdx = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=3; ndy=0; ndcb=2;
     use_y_derivs=1;
     ddg3bdcbdcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);
     nm=3; ndy=1; ndcb=1;
     use_y_derivs=0;
     ddg3bdydcb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, y);

     dg1bdx = dg3bdx + y/xbsq*cb*g2b - y/xb*cb*dg2bdx;
     dg1bdy = dg3bdy - cb/xb*g2b -(y/xb)*cb*dg2bdy;
     dg1bdcb = dg3bdcb - y/xb*g2b - y/xb*cb*dg2bdcb;
     ddg1bdxdx =  ddg3bdxdx -2.0*y/(xbsq*xb)*cb*g2b+y/xbsq*cb*dg2bdx + y/xbsq*cb*dg2bdx - y/xb*cb*ddg2bdxdx;
     ddg1bdxdy = ddg3bdxdy +cb/xbsq*g2b + y*cb/xbsq*dg2bdy -cb/xb*dg2bdx - y/xb*cb*ddg2bdxdy;
     ddg1bdxdcb = ddg3bdxdcb +y/xbsq*g2b +y/xbsq*cb*dg2bdcb -y/xb*dg2bdx -y/xb*cb*ddg2bdxdcb;
     ddg1bdydcb = ddg3bdydcb - g2b/xb -cb/xb*dg2bdcb -y/xb*dg2bdy - y/xb*cb*ddg2bdydcb;
     ddg1bdcbdcb = ddg3bdcbdcb - 2.0*y/xb*dg2bdcb - y/xb*cb*ddg2bdcbdcb;
  
     dg1adx = (ay*G31A_tay[iy_tay]+(1.0-ay)*G31A_tay[iy2_tay]-(2.0/3)*(ay*G31B_tay[iy_tay]+(1.0-ay)*G31B_tay[iy2_tay])
	       -y*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]))*cb;
     dg1ady = ay*G3Ady_tay[iy_tay]+(1.0-ay)*G3Ady_tay[iy2_tay]-(ay*G3Bdy_tay[iy_tay]+(1.0-ay)*G3Bdy_tay[iy2_tay]);
     dg1adcb = 0.0;

     ddg1adxdy = (ay*G31Ady_tay[iy_tay]+(1.0-ay)*G31Ady_tay[iy2_tay]-(2.0/3)*(ay*G31Bdy_tay[iy_tay]+(1.0-ay)*G31Bdy_tay[iy2_tay])
		  -(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay])-y*(ay*G22Ady_tay[iy_tay]+(1.0-ay)*G22Ady_tay[iy2_tay]))*cb;
     ddg1adxdcb = (ay*G31A_tay[iy_tay]+(1.0-ay)*G31A_tay[iy2_tay]-(2.0/3)*(ay*G31B_tay[iy_tay]+(1.0-ay)*G31B_tay[iy2_tay])
	       -y*(ay*G22A_tay[iy_tay]+(1.0-ay)*G22A_tay[iy2_tay]));
     ddg1adxdx = ddg1bdxdx; 
     ddg1adcbdcb = 0.0;
     ddg1adydcb = 0.0;

     ddg1dxdcb = ax*ddg1adxdcb + (1.0-ax)*ddg1bdxdcb;
     dg1dx = ddg1dxdcb*cb; // ax*dg1adx + (1.0-ax)*dg1bdx;     
     dg1dcb = ddg1dxdcb*x; // ax*dg1adcb + (1.0-ax)*dg1bdcb;     
     dg1dy = ax*dg1ady + (1.0-ax)*dg1bdy;
     ddg1dxdy = ax*ddg1adxdy + (1.0-ax)*ddg1bdxdy;
     ddg1dxdx = ax*ddg1adxdx + (1.0-ax)*ddg1bdxdx;
     ddg1dcbdcb = ddg1bdcbdcb*xsq/xbsq; // the function starts quadratically 
     ddg1dydcb = ax*ddg1adydcb + (1.0-ax)*ddg1bdydcb;

     //printf("dg1adx= %.11lg  dg1bdx= %.11lg\n", dg1adx, dg1bdx);
     //printf("dg1ady= %.11lg  dg1bdy= %.11lg\n", dg1ady, dg1bdy);
     //printf("dg1adcb= %.11lg  dg1bdcb= %.11lg\n", dg1adcb, dg1bdcb);
     //printf("ddg1adxdy= %.11lg ddg1bdxdy= %.11lg\n", ddg1adxdy, ddg1bdxdy);
     //printf("ddg1adxdcb= %.11lg ddg1bdxdcb= %.11lg\n", ddg1adxdcb, ddg1bdxdcb);
     //printf("ddg1adxdx= %.11lg ddg1bdxdx= %.11lg\n", ddg1adxdx, ddg1bdxdx);
     //printf("ddg1adcbdcb= %.11lg ddg1bdcbdcb= %.11lg\n", ddg1adcbdcb, ddg1bdcbdcb);
     //printf("ddg1adydcb= %.11lg ddg1bdydcb= %.11lg\n", ddg1adydcb, ddg1bdydcb);

#ifdef XMYSWAP_V
     if(flag2) {
       // in the case of a swap, you need the second derivative wrt y
       ya = YY[iy];
       yb = ya+stpy;
       ddg2adydy= 0.0;
       nm=2; ndy=1; ndcb=0;
       use_y_derivs=0;
       dg2bdy_ya = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, ya);
       dg2bdy_yb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, yb);
       ddg2bdydy = (dg2bdy_yb-dg2bdy_ya)/stpy;
       ddg2dydy = ax*ddg2adydy + (1.0-ax)*ddg2bdydy;
       //       printf("Doing x-Taylor expansion after swap: iy= %d  y= %lg  ya= %lg  yb= %lg  fq1= %lg  fq2= %lg  ddg2dydy= %lg \n", iy, y, ya, yb, dg2bdy_ya, dg2bdy_yb, ddg2dydy);
       
       ddg1adydy = (G3Ady_tay[iy2_tay]-G3Bdy_tay[iy2_tay]-(G3Ady_tay[iy_tay]-G3Bdy_tay[iy_tay]))/Stp_tay; 
       // ay*G3Ady_tay[iy_tay]+(1.0-ay)*G3Ady_tay[iy2_tay]-(ay*G3Bdy_tay[iy_tay]+(1.0-ay)*G3Bdy_tay[iy2_tay]);
       nm=3; ndy=1; ndcb=0;
       use_y_derivs=0;
       dg3bdy_ya = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, ya);
       dg3bdy_yb = accessv(flag_hy, use_y_derivs, ix, iy, nm, ndy, ndcb, &x2, &y1, &y2, x, &cb, yb);
       //     dg3bdy_ya = extractff(3, 1, 0, x, cb, ya);     dg3bdy_yb = extractff(3, 1, 0, x, cb, yb);
       ddg3bdydy = (dg3bdy_yb-dg3bdy_ya)/stpy;
       ddg1bdydy = ddg3bdydy -(y*ddg2bdydy+2.0*dg2bdy)*cb/x;
       ddg1dydy = ax*ddg1adydy + (1.0-ax)*ddg1bdydy;
     }
#endif
   } else 
#endif
     {
  g2 = extractff(2, 0, 0, x, cb, y);
  //  g3 = extractff(3, 0, 0, x, cb, y);
  //  g1 = g3 -(y/x)*cb*g2;

  dg2dx = extractff(12, 0, 0, x, cb, y);
  dg2dy = extractff(2, 1, 0, x, cb, y);
  dg2dcb = extractff(2, 0, 1, x, cb, y);
  dg3dx = extractff(13, 0, 0, x, cb, y);
  dg3dy = extractff(3, 1, 0, x, cb, y);
  dg3dcb = extractff(3, 0, 1, x, cb, y);
  dg1dx = dg3dx + y/xsq*cb*g2 - y/x*cb*dg2dx;
  dg1dy = dg3dy - cb/x*g2 -(y/x)*cb*dg2dy;
  dg1dcb = dg3dcb - y/x*g2 - y/x*cb*dg2dcb;

  ddg2dxdx = extractff(22, 0, 0, x, cb, y);
  ddg2dxdy = extractff(12, 1, 0, x, cb, y);
  ddg2dxdcb = extractff(12, 0, 1, x, cb, y);
  ddg2dcbdcb = extractff(2, 0, 2, x, cb, y);
  ddg2dydcb = extractff(2, 1, 1, x, cb, y);

  ddg3dxdx = extractff(23, 0, 0, x, cb, y);
  ddg3dxdy = extractff(13, 1, 0, x, cb, y);
  ddg3dxdcb = extractff(13, 0, 1, x, cb, y);
  ddg3dcbdcb = extractff(3, 0, 2, x, cb, y);
  ddg3dydcb = extractff(3, 1, 1, x, cb, y);

  ddg1dxdx =  ddg3dxdx -2.0*y/(xsq*x)*cb*g2+y/xsq*cb*dg2dx + y/xsq*cb*dg2dx - y/x*cb*ddg2dxdx;
  ddg1dxdy = ddg3dxdy +cb/xsq*g2 + y*cb/xsq*dg2dy -cb/x*dg2dx - y/x*cb*ddg2dxdy;
  ddg1dxdcb = ddg3dxdcb +y/xsq*g2 +y/xsq*cb*dg2dcb -y/x*dg2dx -y/x*cb*ddg2dxdcb;
  ddg1dydcb = ddg3dydcb - g2/x -cb/x*dg2dcb -y/x*dg2dy - y/x*cb*ddg2dydcb;
  ddg1dcbdcb = ddg3dcbdcb - 2.0*y/x*dg2dcb - y/x*cb*ddg2dcbdcb;

#ifdef XMYSWAP_V
  if(flag2) {
    // in case of swap, you need the second derivative wrt y
    stpy=YY[1]-YY[0];
    iy1 = y/stpy;
    //   iy = (iy>=1 ? iy-1 : iy);
    ya = iy1*stpy;
    yb = ya+stpy;
    fq1 = extractff(2, 1, 0, x, cb, ya);
    fq2 = extractff(2, 1, 0, x, cb, yb);
    ddg2dydy = (fq2-fq1)/stpy;
    //    printf("flag2= %d  y= %lg  ya= %lg  yb= %lg  fq1= %lg  fq2= %lg  ddg2dydy= %lg \n", flag2, y, ya, yb, fq1, fq2, ddg2dydy);
    fq1 = extractff(3, 1, 0, x, cb, ya);
    fq2 = extractff(3, 1, 0, x, cb, yb);
    ddg3dydy = (fq2-fq1)/stpy; 
    ddg1dydy = ddg3dydy -(y*ddg2dydy+2.0*dg2dy)*cb/x;
  }
#endif
     }

#ifdef XMYSWAP_V
   if(flag2) {
     // here, convert the derivatives of g2(x,ca,xmy) into the derivatives of g2(x,cb,y)
     g2s = -g2;
     dgs2dx  = -dg2dx;
     dgs2dcb = -dg2dcb;
     dgs2dy  = -dg2dy;
     ca = cb; 
     y = yorig;
     cb = cborig;
     ysq = ysqorig;
     cd = (y-x*cb)/xmy;
     g2 = g2s;
     dg2dx = dgs2dx + (1.0-ca*ca)/xmy*dgs2dcb + ca * dgs2dy;
     dg2dcb = -(ysq*cd*dgs2dcb + x*y*xmy*dgs2dy)/xmysq;
     dg2dy = -(cb+ca*cd)/xmy*dgs2dcb+cd*dgs2dy;
     //     g1s = g1;
     dgs1dx  = dg1dx;
     dgs1dcb = dg1dcb;
     dgs1dy  = dg1dy;
     //     g1 = g1s - g2;
     dg1dx = dgs1dx + (1.0-ca*ca)/xmy*dgs1dcb + ca * dgs1dy - dg2dx;
     dg1dcb = -(ysq*cd*dgs1dcb + x*y*xmy*dgs1dy)/xmysq - dg2dcb;
     dg1dy = -(cb+ca*cd)/xmy*dgs1dcb+cd*dgs1dy - dg2dy;

     ddgs2dxdx  = -ddg2dxdx;
     ddgs2dxdy  = -ddg2dxdy;
     ddgs2dxdcb  = -ddg2dxdcb;
     ddgs2dydcb  = -ddg2dydcb;
     ddgs2dcbdcb  = -ddg2dcbdcb;
     ddgs2dydy  = -ddg2dydy;

     ddgs1dxdx  = ddg1dxdx;
     ddgs1dxdy  = ddg1dxdy;
     ddgs1dxdcb  = ddg1dxdcb;
     ddgs1dydcb  = ddg1dydcb;
     ddgs1dcbdcb  = ddg1dcbdcb;
     ddgs1dydy  = ddg1dydy;

     ddg2dxdx = (ddgs2dxdx + (1.0-ca*ca)/xmy*ddgs2dxdcb + ca * ddgs2dxdy)
       - 3.0*ca*(1.0-ca*ca)/xmysq*dgs2dcb + (1.0-ca*ca)/xmy*(ddgs2dxdcb + (1.0-ca*ca)/xmy*ddgs2dcbdcb + ca * ddgs2dydcb)
       + (1.0-ca*ca)/xmy*dgs2dy + ca*(ddgs2dxdy + (1.0-ca*ca)/xmy*ddgs2dydcb + ca * ddgs2dydy);
     ddg2dxdcb = ysq/(xmysq*xmy)*(cb+3.0*ca*cd)*dgs2dcb-ysq/xmysq*cd*(ddgs2dxdcb + (1.0-ca*ca)/xmy*ddgs2dcbdcb + ca * ddgs2dydcb)
       + y/xmysq*(x*ca-xmy)*dgs2dy-x*y/xmy*(ddgs2dxdy + (1.0-ca*ca)/xmy*ddgs2dydcb + ca * ddgs2dydy);
     ddg2dxdy = (-(cb+ca*cd)/xmy*ddgs2dxdcb+cd*ddgs2dxdy)+(2.0*ca*cb+cd*(3.0*ca*ca-1.0))/xmysq*dgs2dcb
       +(1.0 - ca*ca)/xmy*(-(cb + ca*cd)/xmy*ddgs2dcbdcb + cd*ddgs2dydcb)
       -(cb+ca*cd)/xmy*dgs2dy + ca*(-(cb+ca*cd)/xmy*ddgs2dydcb+cd*ddgs2dydy);
     ddg2dydcb = y/(xmysq*xmy)*((3.0*cd*cd-1.0)*y-2*xmy*cd)*dgs2dcb-ysq/(xmysq)*cd*(-(cb+ca*cd)/xmy*ddgs2dcbdcb+cd*ddgs2dydcb)
       + x/xmysq*(y*cd-xmy)*dgs2dy-x*y/xmy*(-(cb+ca*cd)/xmy*ddgs2dydcb+cd*ddgs2dydy);
     ddg2dcbdcb = x*ysq/(xmysq*xmysq)*(xmy-3.0*y*cd)*dgs2dcb - ysq/(xmysq)*cd*(-(ysq*cd*ddgs2dcbdcb + x*y*xmy*ddgs2dydcb)/xmysq)
       -xsq*ysq/(xmysq*xmy)*dgs2dy - x*y/xmy*(-(ysq*cd*ddgs2dydcb + x*y*xmy*ddgs2dydy)/xmysq);

     ddg1dxdx = (ddgs1dxdx + (1.0-ca*ca)/xmy*ddgs1dxdcb + ca * ddgs1dxdy)
       - 3.0*ca*(1.0-ca*ca)/xmysq*dgs1dcb + (1.0-ca*ca)/xmy*(ddgs1dxdcb + (1.0-ca*ca)/xmy*ddgs1dcbdcb + ca * ddgs1dydcb)
       + (1.0-ca*ca)/xmy*dgs1dy + ca*(ddgs1dxdy + (1.0-ca*ca)/xmy*ddgs1dydcb + ca * ddgs1dydy) - ddg2dxdx;
     ddg1dxdcb = ysq/(xmysq*xmy)*(cb+3.0*ca*cd)*dgs1dcb-ysq/xmysq*cd*(ddgs1dxdcb + (1.0-ca*ca)/xmy*ddgs1dcbdcb + ca * ddgs1dydcb)
       + y/xmysq*(x*ca-xmy)*dgs1dy-x*y/xmy*(ddgs1dxdy + (1.0-ca*ca)/xmy*ddgs1dydcb + ca * ddgs1dydy) - ddg2dxdcb;
     ddg1dxdy = (-(cb+ca*cd)/xmy*ddgs1dxdcb+cd*ddgs1dxdy)+(2.0*ca*cb+cd*(3.0*ca*ca-1.0))/xmysq*dgs1dcb
       +(1.0 - ca*ca)/xmy*(-(cb + ca*cd)/xmy*ddgs1dcbdcb + cd*ddgs1dydcb)
       -(cb+ca*cd)/xmy*dgs1dy + ca*(-(cb+ca*cd)/xmy*ddgs1dydcb+cd*ddgs1dydy) - ddg2dxdy;
     ddg1dydcb = y/(xmysq*xmy)*((3.0*cd*cd-1.0)*y-2*xmy*cd)*dgs1dcb-ysq/(xmysq)*cd*(-(cb+ca*cd)/xmy*ddgs1dcbdcb+cd*ddgs1dydcb)
       + x/xmysq*(y*cd-xmy)*dgs1dy-x*y/xmy*(-(cb+ca*cd)/xmy*ddgs1dydcb+cd*ddgs1dydy) - ddg2dydcb;
     ddg1dcbdcb = x*ysq/(xmysq*xmysq)*(xmy-3.0*y*cd)*dgs1dcb - ysq/(xmysq)*cd*(-(ysq*cd*ddgs1dcbdcb + x*y*xmy*ddgs1dydcb)/xmysq)
       -xsq*ysq/(xmysq*xmy)*dgs1dy - x*y/xmy*(-(ysq*cd*ddgs1dydcb + x*y*xmy*ddgs1dydy)/xmysq) - ddg2dcbdcb;

     /*     ddg1dxdx = (ddgs1dxdx + (1.0-ca*ca)/xmy*ddgs1dxdcb + ca * ddgs1dxdy - ddg2dxdx)
       - ca*(1.0-ca*ca)/xmysq*dgs1dcb + (1.0-ca*ca)/xmy*(ddgs1dxdcb + (1.0-ca*ca)/xmy*ddgs1dcbdcb + ca * ddgs1dydcb - ddg2dxdcb)
       + (1.0-ca*ca)/xmy*dgs1dy + ca*(ddgs1dxdy + (1.0-ca*ca)/xmy*ddgs1dydcb + ca * ddgs1dydy - ddg2dxdy) - ddg2dxdx;*/

     //   printf("flag2= %d YY[nstpy_ff-1]= %lg  xmy= %lg rx= %lg  rxy= %lg\n", flag2, YY[nstpy_ff-1], xmy, rx, rxy);
   } 
#endif
   /*
  printf("x= %.11lg   cb= %.11lg   y= %.11lg \n", x, cb, y);
  printf("g2= %.11lg\n", g2);
  printf("dg2dx= %.11lg   dg2dy= %.11lg   dg2dcb= %.11lg\n", dg2dx, dg2dy, dg2dcb); 
  printf("ddg2dxdx= %.11lg   ddg2dxdy= %.11lg   ddg2dxdcb= %.11lg  ddg2dcbdcb= %.11lg  ddg2dydcb= %.11lg\n", 
	 ddg2dxdx, ddg2dxdy, ddg2dxdcb, ddg2dcbdcb, ddg2dydcb);
  printf("dg1dx= %.11lg   dg1dy= %.11lg   dg1dcb= %.11lg\n", dg1dx, dg1dy, dg1dcb);
    printf("ddg1dxdx= %.11lg   ddg1dxdy= %.11lg   ddg1dxdcb= %.11lg  ddg1dcbdcb= %.11lg  ddg1dydcb= %.11lg\n", 
	 ddg1dxdx, ddg1dxdy, ddg1dxdcb, ddg1dcbdcb, ddg1dydcb);   exit(0); 
   */
#endif 

  for(bet=0;bet<=3;bet++)  {  
    yhat[bet] = yv[bet]/y;
    xhat[bet] = xv[bet]/x;
    c2v[bet] = (xhat[bet]-cb*yhat[bet])/y;
    c4v[bet] = (yhat[bet]-cb*xhat[bet])/x;
  }

  //  printf("dg1dx= %lg dg1dy= %lg dg1dcb= %lg ddg1dydcb= %lg\n", dg1dx, dg1dy, dg1dcb, ddg1dydcb);
  //  printf("ddg1dxdx= %lg  ddg1dxdy= %lg   ddg1dxdcb= %lg ddg1dcbdcb= %lg\n", ddg1dxdx, ddg1dxdy, ddg1dxdcb, ddg1dcbdcb);
  // dg1dx = dg1dcb = dg1dy = ddg1dxdx = ddg1dxdy = ddg1dxdcb = ddg1dydcb = ddg1dcbdcb = 0.0; 
  //  dg2dx=dg2dcb=ddg2dxdx=ddg2dxdy=ddg2dxdcb=ddg2dydcb=ddg2dcbdcb=0.0; 

  for(bet=0;bet<=3;bet++)  {  
    for(alf=0;alf<=3;alf++)  {
      d_ab = (alf==bet ? 1 : 0);
      for(dta=0;dta<=3;dta++) {
	d_bd = (bet==dta ? 1 : 0);
	d_ad = (alf==dta ? 1 : 0);
	f = d_bd*((dg1dx+dg2dx)*xhat[alf]+(dg1dcb+dg2dcb)*c4v[alf]); 
	f+= d_ad*(dg1dx*xhat[bet]+dg1dy*yhat[bet]+dg1dcb*(xhat[bet]/y+yhat[bet]/x-cb*(xhat[bet]/x+yhat[bet]/y)));
	f+= (xv[dta]*ddg1dxdx+yv[dta]*ddg2dxdx)*xhat[alf]*xhat[bet]+(xv[dta]*ddg1dxdy+yv[dta]*ddg2dxdy)*xhat[alf]*yhat[bet];
	f+=((xv[dta]*ddg1dxdcb+yv[dta]*ddg2dxdcb)*xhat[bet]+(xv[dta]*ddg1dydcb+yv[dta]*ddg2dydcb)*yhat[bet])*c4v[alf];
	f+=((xv[dta]*ddg1dxdcb+yv[dta]*ddg2dxdcb)*xhat[alf]+(xv[dta]*ddg1dcbdcb+yv[dta]*ddg2dcbdcb)*c4v[alf])*(c2v[bet]+c4v[bet]);
	f+=(xv[dta]*dg1dx+yv[dta]*dg2dx)*(d_ab-xhat[alf]*xhat[bet])/x;
	dv[alf][bet][dta]=f+(xv[dta]*dg1dcb+yv[dta]*dg2dcb)/x*
	  ((cb*(3.0*xhat[alf]*xhat[bet]-d_ab)-xhat[alf]*yhat[bet]-yhat[alf]*xhat[bet])/x
	   +(d_ab+xhat[alf]*yhat[bet]*cb-xhat[alf]*xhat[bet]-yhat[alf]*yhat[bet])/y);
      }
    }
  }

#if (1==0)
#define EPSI 1.0e-5
  for(bet=0;bet<=3;bet++)  {
    for(alf=0;alf<=3;alf++)  {
      for(dta=0;dta<=3;dta++) {
	if(alf==bet) {
	  SETINVARIANTS;
	  gtmp1 = yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  xv[alf] += EPSI;
	  SETINVARIANTS;
	  gtmp2 = yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  xv[alf] -= 2*EPSI;
	  SETINVARIANTS;
	  gtmp3 = yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  gtmp4 = (gtmp2-2.0*gtmp1+gtmp3)/(EPSI*EPSI);
	  xv[alf] += EPSI;
	} else {
	  xv[alf] += EPSI; xv[bet] += EPSI;
	  SETINVARIANTS;
	  gtmp4 = yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  xv[bet] -= 2*EPSI;
	  SETINVARIANTS;
	  gtmp4 -= yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  xv[alf] -= 2*EPSI;
	  SETINVARIANTS;
	  gtmp4 += yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  xv[bet] += 2*EPSI;  // printf("xv[%d]= %lg  yv[%d]= %lg\n", alf, xv[alf], bet, yv[bet]);
	  SETINVARIANTS;
	  gtmp4 -= yv[dta]*testff(1, x, cb, y, &f, vtmp);
	  gtmp4 /= (4.0*EPSI*EPSI);
	  xv[bet] -= EPSI;
	  xv[alf] += EPSI;
	}
	xv[alf] += EPSI; yv[bet] += EPSI;  // printf("xv[%d]= %lg  yv[%d]= %lg\n", alf, xv[alf], bet, yv[bet]);
	SETINVARIANTS;
	gtmp1 = yv[dta]*testff(1, x, cb, y, &f, vtmp);
	yv[bet] -= 2*EPSI;  // printf("xv[%d]= %lg  yv[%d]= %lg\n", alf, xv[alf], bet, yv[bet]);
	SETINVARIANTS;
	gtmp1 -= yv[dta]*testff(1, x, cb, y, &f, vtmp);
	xv[alf] -= 2*EPSI;  // printf("xv[%d]= %lg  yv[%d]= %lg\n", alf, xv[alf], bet, yv[bet]);
	SETINVARIANTS;
	gtmp1 += yv[dta]*testff(1, x, cb, y, &f, vtmp);
	yv[bet] += 2*EPSI;  // printf("xv[%d]= %lg  yv[%d]= %lg\n", alf, xv[alf], bet, yv[bet]);
	SETINVARIANTS;
	gtmp1 -= yv[dta]*testff(1, x, cb, y, &f, vtmp);
	gtmp1 /= (4.0*EPSI*EPSI);
	yv[bet] -= EPSI;
	xv[alf] += EPSI;
	if(fabs((gtmp1+gtmp4)/dv[alf][bet][dta]-1.0)>5.0e-4) 
	  printf("!!! alf= %d bet= %d  dta= %d  chainrule: %lg  direct: %lg\n", alf, bet, dta, dv[alf][bet][dta], (gtmp1+gtmp4));
      }
    }
  }
  exit(0); // HACK  
#undef EPSI
#endif
}

/*
// this routine is not needed,
//  it just shows how to use gI(..), gII(..) and gIII(..) 
#define NMAX 4096

int computeG()
{
  int i, j, n, mu, nu, rho, sig, rhosig, lda, idx[20], irs, dta, sgn, cnt;
  int *gt[6][4], *alfv[6][4], *betv[6][4], *dtav[6][4], nv[6][4];

  for(i=0;i<6;i++)  for(j=0;j<4;j++)  {
      gt[i][j] = (int *) malloc(NMAX*sizeof(int));
      alfv[i][j] = (int *) malloc(NMAX*sizeof(int));
      betv[i][j] = (int *) malloc(NMAX*sizeof(int));
      dtav[i][j] = (int *) malloc(NMAX*sizeof(int));
    }


  cnt=0;
  for(mu=0;mu<=3;mu++) for(nu=0;nu<=3;nu++) {
      gIII(mu, nu, nv, gt, alfv, betv, dtav);
      for(irs=0;irs<12;irs++) {
	rhosig=rhosigw[irs];
	sgn = sgnw[irs];
	rho = rhow[irs];
	sig = sigw[irs];
	for(lda=0;lda<=3;lda++) {
	  for(i=0;i<nv[rhosig][lda];i++) {
	    printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
		   mu, nu, rho, sig, lda, alfv[rhosig][lda][i], betv[rhosig][lda][i], 
		   dtav[rhosig][lda][i], sgn*gt[rhosig][lda][i]);
	    ++cnt;
	  }
	}
      }
    }

}
#undef NMAX
  */
