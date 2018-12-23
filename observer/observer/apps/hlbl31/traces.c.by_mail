

#include "getL.h"



static int rhosigv[4][4], rhosigw[12], sgnw[12], rhow[12], sigw[12];


void init_dirac()
{
  rhosigw[0] = 0;  sgnw[0] = 1;   rhow[0]=0;  sigw[0]=1;
  rhosigw[1] = 1;  sgnw[1] = 1;   rhow[1]=0;  sigw[1]=2;
  rhosigw[2] = 2;  sgnw[2] = 1;   rhow[2]=0;  sigw[2]=3;
  rhosigw[3] = 0;  sgnw[3] =-1;   rhow[3]=1;  sigw[3]=0;

  rhosigw[4] = 3;  sgnw[4] = 1;   rhow[4]=1;  sigw[4]=2;
  rhosigw[5] = 4;  sgnw[5] = 1;   rhow[5]=1;  sigw[5]=3;
  rhosigw[6] = 1;  sgnw[6] =-1;   rhow[6]=2;  sigw[6]=0;
  rhosigw[7] = 3;  sgnw[7] =-1;   rhow[7]=2;  sigw[7]=1;

  rhosigw[8] = 5;  sgnw[8] = 1;   rhow[8]=2;  sigw[8]=3;
  rhosigw[9] = 2;  sgnw[9] = -1;   rhow[9]=3;  sigw[9]=0;
  rhosigw[10]= 4;  sgnw[10] = -1;   rhow[10]=3;  sigw[10]=1;
  rhosigw[11]= 5;  sgnw[11] = -1;   rhow[11]=3;  sigw[11]=2;

  rhosigv[0][1] = 0;
  rhosigv[0][2] = 1;
  rhosigv[0][3] = 2;
  rhosigv[1][2] = 3;
  rhosigv[1][3] = 4;
  rhosigv[2][3] = 5;
}



void gI(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4])
{
  int rho, sig, lda, mu, nu;
  int alf, bet, dta;
  int i, rhosig, g;

  for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) {
      rhosig=-1;
      for(rho=0;rho<=2;rho++)   for(sig=rho+1;sig<=3;sig++)   {
	  ++rhosig;
	  for(lda=0;lda<=3;lda++)  {
	    i=0;
	    for(alf=0;alf<=3;alf++)   for(bet=0;bet<=3;bet++)   for(dta=0;dta<=3;dta++) {	  
		  g=getgI(rho,sig,mu,nu,lda,alf,bet,dta);
#ifdef SYMG
		  g-=getgI(rho,sig,lda,nu,mu,bet,alf,dta);
#endif		 
		  if(g!=0) {
		    gt[mu][nu][rhosig][lda][i] = g;
		    alfv[mu][nu][rhosig][lda][i]= alf;
		    betv[mu][nu][rhosig][lda][i]= bet;
		    dtav[mu][nu][rhosig][lda][i]= dta;
		    ++i;
		  }
		}
	    nv[mu][nu][rhosig][lda]=i;
	    if(i>NMAX) { printf("Exceeding array size.\n"); exit(0); }
	  }
	}
    }
}


void gII(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4])
{
  int rho, sig, lda, mu, nu;
  int alf, bet, dta;
  int i, rhosig, g, h;

  for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) {
      rhosig=-1;
      for(rho=0;rho<=2;rho++)   for(sig=rho+1;sig<=3;sig++)   {
	  ++rhosig;
	  for(lda=0;lda<=3;lda++)  {
	    i=0;
	    for(bet=lda, alf=0;alf<=3;alf++)   for(dta=0;dta<=3;dta++) {	  
		g = getgII(rho,sig,mu,nu,lda,alf,bet,dta);
#ifdef SYMG
		h = getgII(rho,sig,lda,alf,mu,nu,bet,dta);
		h-= getgII(rho,sig,lda,mu,bet,nu,alf,dta);
		h+= getgII(rho,sig,lda,mu,bet,alf,nu,dta);		
		g -= h;
#endif
		if(g!=0) {
		  gt[mu][nu][rhosig][lda][i] = g;
		  alfv[mu][nu][rhosig][lda][i]= alf;
		  betv[mu][nu][rhosig][lda][i]= bet;
		  dtav[mu][nu][rhosig][lda][i]= dta;
		  ++i;
		}
	      }
	    nv[mu][nu][rhosig][lda]=i;
	    if(i>NMAX) { printf("Exceeding array size.\n"); exit(0); }
	  }
	}
    }
}


void gIII(int nv[4][4][6][4], int *gt[4][4][6][4], int *alfv[4][4][6][4], int *betv[4][4][6][4], int *dtav[4][4][6][4])
{
  int rho, sig, lda, mu, nu;
  int alf, bet, dta;
  int i, rhosig, g;

  for(mu=0;mu<4;mu++) for(nu=0;nu<4;nu++) {
      rhosig=-1;
      for(rho=0;rho<=2;rho++)   for(sig=rho+1;sig<=3;sig++)   {
	  ++rhosig;
	  for(lda=0;lda<=3;lda++) {
	    i=0;
	    for(alf=0;alf<=3;alf++)   for(bet=0;bet<=3;bet++)  for(dta=0;dta<=3;dta++) {	  
		  g = getgII(rho,sig,mu,bet,lda,nu,alf,dta);
		  g-= getgII(rho,sig,mu,lda,alf,nu,bet,dta);
		  g+= getgII(rho,sig,mu,lda,alf,bet,nu,dta);
#ifdef SYMG
		  g -= getgII(rho,sig,lda,nu,mu,bet,alf,dta);
#endif
		  if(g!=0) {
		    gt[mu][nu][rhosig][lda][i] = g;
		    alfv[mu][nu][rhosig][lda][i]= alf;
		    betv[mu][nu][rhosig][lda][i]= bet;
		    dtav[mu][nu][rhosig][lda][i]= dta;
		    ++i;
		  }
		}
	    nv[mu][nu][rhosig][lda]=i;	 
	    if(i>NMAX) { printf("Exceeding array size.\n"); exit(0); }
	  }
	}
    }
}


int getgI(int rho, int sig, int mu, int nu, int lda, int alf, int bet, int dta)
{
  int idx[8];
  int g;

  if(rho==sig) return(0);
  idx[0]=dta;
  idx[1]=rho;
  idx[2]=sig;
  idx[3]=mu;
  idx[4]=alf;
  idx[5]=nu;
  idx[6]=bet;
  idx[7]=lda;
  g = diractrace(8,idx);
  if(dta==sig) {
    idx[0]=rho;
    idx[1]=mu;
    idx[2]=alf;
    idx[3]=nu;
    idx[4]=bet;
    idx[5]=lda;
    g += diractrace(6,idx);
  }
  if(dta==rho) {
    idx[0]=sig;
    idx[1]=mu;
    idx[2]=alf;
    idx[3]=nu;
    idx[4]=bet;
    idx[5]=lda;
    g -= diractrace(6,idx);
  }
  return(g);
}

int getgII(int rho, int sig, int mu, int nu, int lda, int alf, int bet, int dta)
{
  int idx[6], g;

  if(rho==sig || bet!=lda) return(0);
  idx[0]=dta;
  idx[1]=rho;
  idx[2]=sig;
  idx[3]=mu;
  idx[4]=alf;
  idx[5]=nu;
  g = diractrace(6,idx);
  if(dta==sig) {
    idx[0]=rho;
    idx[1]=mu;
    idx[2]=alf;
    idx[3]=nu;
    g += diractrace(4,idx);
  }
  if(dta==rho) {
    idx[0]=sig;
    idx[1]=mu;
    idx[2]=alf;
    idx[3]=nu;
    g -= diractrace(4,idx);
  }
  return(-2*g);
}

// returns 1/4*tr(gamma_{idx[0]}*gamma_{idx[1]}...*gamma_{idx[n-1]})
int diractrace(int n, int *idx)
{
  int j, k, i, it, *idxred;

  if(n%2==1)  return 0;
  if(n==0)  return 1;
  if(n==2) return((idx[0]==idx[1] ? 1 : 0));

  idxred = malloc((n-2)*sizeof(int));
  for(it=0, i=0;i<n-1;i++) if(idx[i]==idx[n-1]) {
      for(j=k=0;j<n-1;j++) if(j!=i) idxred[k++] = idx[j];
      it += (i%2==0 ? 1 : -1)*diractrace(n-2,idxred);
    }
  free(idxred);
  return(it);
}


