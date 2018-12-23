/**
* @file smearing.cc
* \brief  Smearing functions
* @author Jeremy Green
* @version 
* @date 2014-12-01
*/
#include <qdp.h>
#include "smearing.h"


/* both forward and backward staples */
void staple(QDP::LatticeColorMatrixD & out,
            const QDP::LatticeColorMatrixD & links_parallel,
            const QDP::LatticeColorMatrixD & links_perp,
            int dir_parallel, int dir_perp)
{
  QDP::LatticeColorMatrixD sperp = QDP::shift(links_perp, 1, dir_parallel);
  out = links_perp * QDP::shift(links_parallel, 1, dir_perp) * adj(sperp)
    + QDP::shift(adj(links_perp) * links_parallel * sperp, -1, dir_perp);
}

static inline QDP::LatticeRealD det_im(const QDP::LatticeColorMatrixD & u)
{
  QDP::LatticeComplexD det1, det2, det3;
  det1 = QDP::peekColor(u, 1, 1) * QDP::peekColor(u, 2, 2)
    - QDP::peekColor(u, 1, 2) * QDP::peekColor(u, 2, 1);
  det2 = QDP::peekColor(u, 1, 0) * QDP::peekColor(u, 2, 2)
    - QDP::peekColor(u, 1, 2) * QDP::peekColor(u, 2, 0);
  det3 = QDP::peekColor(u, 1, 0) * QDP::peekColor(u, 2, 1)
    - QDP::peekColor(u, 1, 1) * QDP::peekColor(u, 2, 0);

  return QDP::imag(QDP::peekColor(u, 0, 0) * det1
                   - QDP::peekColor(u, 0, 1) * det2
                   + QDP::peekColor(u, 0, 2) * det3);
}

static inline void approx_project_to_su3(QDP::LatticeColorMatrixD & u, int iter)
{
  QDP::LatticeRealD d = 1.0 / QDP::sqrt(QDP::localNorm2(u) / 3.0);
  u *= d;
  QDP::LatticeColorMatrixD x;
  for (int i = 0; i < iter; i++)
  {
    x = u * (QDP::adj(u) * u * (-0.5) + 1.5);
    u = QDP::cmplx(1.0, -det_im(x) / 3.0) * x;
  }
}

/* designed to be consistent with measure code
 * may want to be more flexible with projection */
void ape_smear(QDP::multi1d < QDP::LatticeColorMatrixD > &out,
               const QDP::multi1d < QDP::LatticeColorMatrixD > &u,
               double alpha, int n, int skipdir)
{
  QDP::LatticeColorMatrixD tmp;
  QDP::multi1d < QDP::LatticeColorMatrixD > tmp_u(4);
  QDP::multi1d < QDP::LatticeColorMatrixD > *t1, *t2, *t3;

  if (skipdir >= 0 && skipdir < 4)
    out[skipdir] = u[skipdir];

  if (n % 2 == 0)
  {
    t1 = &out;
    t2 = &tmp_u;
  }
  else
  {
    t1 = &tmp_u;
    t2 = &out;
  }

  (*t1) = u;
  for (int step = 0; step < n; step++)
  {
    for (int i = 0; i < 4; i++)
    {
      if (i == skipdir)
        continue;
      (*t2)[i] = (*t1)[i] * (1.0 - alpha);
      for (int j = 0; j < 4; j++)
      {
        if (j == i || j == skipdir)
          continue;
        staple(tmp, (*t1)[i], (*t1)[j], i, j);
        (*t2)[i] += (alpha / 6.0) * tmp;
      }
      approx_project_to_su3((*t2)[i], 7);
    }

    t3 = t1;
    t1 = t2;
    t2 = t3;
  }
}
