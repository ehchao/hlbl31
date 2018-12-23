#ifndef MESON_3PT_H
#define MESON_3PT_H

namespace QDP{

void calculate_meson_3pt(LatticePropagatorD & prop1, SpinMatrixD & Ga, SpinMatrixD &O,int ts, int propid, multi1d<int> & mom , multi1d<int> & pos);
}
#endif
