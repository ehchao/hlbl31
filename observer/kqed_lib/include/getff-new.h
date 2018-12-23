#ifndef GETFF_NEW_H
#define GETFF_NEW_H

// accesses V apparently
double
accessv( const bool flag_hy , const bool use_y_derivs,
	 const int ix, const int iy,
	 const FFidx nm, const bool ndy , const NDCB ndcb ,
	 const double cb , const double y ,
	 const struct Grid_coeffs Grid ) ;

double
extractff( const FFidx nm, const bool ndy, const NDCB ndcb ,
	   const struct invariants Inv ,
	   const struct Grid_coeffs Grid ) ;

// little helper function
static inline double
lerp( const double a ,
      const double T1 ,
      const double T2 )
{
  return a*T1 + (1.0-a)*T2 ;
} ;

// scan up to 2.18 used:
int
nfx( const double x ) ;

// scalar product of two four-vectors
double
SCALPROD( const double xv[4] ,
	  const double yv[4] ) ;

// initialises the invariants struct
struct invariants
set_invariants( const double xv[4] ,
		const double yv[4] ,
		const struct Grid_coeffs Grid ) ;

#endif
