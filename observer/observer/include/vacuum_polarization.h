#include <qdp.h>

void c2_src_vecloc_snk_veccons(QDP::multi2d<QDP::LatticeRealD> &corr,
                               const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
                               const QDP::LatticePropagatorD &prop_l,
                               const QDP::LatticePropagatorD &prop_r);

void c2_VP_mom0(QDP::multi3d<QDP::RealD> &out,
                const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
                const QDP::LatticePropagatorD &prop);
