#ifndef PI_PERT_H
#define PI_PERT_H

// returns contribution of the lbl contribution due to a heavy lepton; mh = (lepton mass)/(muon mass)
// including the (alpha/pi)^3 factor
// from Remiddi-Laporta PLB 301 (1993) 440 Eq. (4)
double
amu_heavlept( const double mh ) ;

// contribution of a light lepton in a_mu, with ml=m_lepton/m_mu
// from Remiddi-Laporta PLB 301 (1993) 440 Eq. (2)
double
amu_lightlept( const double ml ) ;

// bessel function
void
bessk( const int n,
       const double x,
       double *bv ) ;

// initialises 8-index gamma thing
void
init_g8( struct QED_kernel_temps *t ) ;

//
void
ipihatFermLoop_antisym( const double x[4] ,
			const double y[4] ,
			const struct QED_kernel_temps t ,
			double vpihat[6][4][4][4] ) ;

// returns pihat1(rhosig,mu,nu,lambda)
void
pihat1( const double x[4] ,
	const double y[4] ,
	const struct QED_kernel_temps t ,
	double vpihat1[4][4][4][4][4],
	double vpihatr1[4][4][4][4] ) ;

// vacuum polarisation
double
vac_pol( const double ml ) ;

#endif
