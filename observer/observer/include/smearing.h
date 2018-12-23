/*
 * Copyright 2014 Jeremy Green
 */
#include <qdp.h>

template <typename T> 
void wuppertal_smear(T &out, const T &in,
		     const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
		     double alpha, int n, const QDP::Subset &subs=QDP::all,
		     double scalefac=1.0)
{
  double c1 = scalefac/(1.0 + 6.0*alpha);
  double c2 = alpha*c1;

  T tmp;
  T *t1, *t2, *t3;

  if (n % 2 == 0) {
    t1 = &out;
    t2 = &tmp;
  } else {
    t1 = &tmp;
    t2 = &out;
  }

  (*t1) = in;
  (*t2) = QDP::zero;
  for (int i = 0; i < n; i++) {
    (*t2)[subs] = (*t1) * c1;
    for (int j = 0; j < 3; j++) { // only spatial directions
      (*t2)[subs] += u[j] * QDP::shift(*t1,1,j) * c2;
      (*t2)[subs] += QDP::shift(adj(u[j]) * (*t1), -1, j) * c2;
    }

    t3 = t1;
    t1 = t2;
    t2 = t3;
  }
}

/* for compatibility with measure code */
template <typename T>
void jacobi_smear(T &out, const T &in,
		  const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
		  double kappa, int n, const QDP::Subset &subs=QDP::all)
{
  double alpha = kappa/(1.0+kappa);
  double scalefac = (1.0+6.0*alpha)*(1.0+kappa)/(1.0+6.0*kappa);
  wuppertal_smear(out, in, u, alpha, n, subs, scalefac);
  out[subs] *= 1.0/QDP::sqrt(n*QDP::Layout::vol());
}

void ape_smear(QDP::multi1d<QDP::LatticeColorMatrixD> &out,
	       const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
	       double alpha, int n=1, int skipdir=3);
