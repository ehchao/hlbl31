/**
* @file twistedbc.cc
* \brief  Applies twisted boundary conditions to gauge links
* @author Dalibor Djukanovic
* @version 0.1
* @date 2015-02-06
*/

#include "qdp.h"
#include "observ.h"
#include "qdpinterface.h"


namespace QDP
{
/* --------------------------------------------------------------------------*/
/**
* \brief  Apply twist angle to gauge link, i.e.
* \exp[i \theta_\mu/L_\mu]
*
* \param mom
* \param u
* \param 
*/
/* ----------------------------------------------------------------------------*/
  void twistedbc(multi1d < Double > &theta, multi1d < LatticeColorMatrixD > &u,
                 multi1d < LatticeColorMatrixD > &unew)
  {

    ComplexD phase;

    for (int i = 0; i < 4; i++)
    {
      phase =
        cmplx(RealD(cos(theta[i] / lat_geom[i])),
              RealD(sin(theta[i] / lat_geom[i])));
      unew[i] = phase * u[i];
    }
  }
}
