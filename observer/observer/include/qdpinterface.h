#ifndef QDP_INTERFACE_H
#define QDP_INTERFACE_H
#include "qdp.h"
#include "sftmom.h"
#include "hdf5_output.h"
extern "C"
{
#include "global.h"
#include "su3.h"
}

#if defined MAIN_PROGRAM
#define EXTERN
#else
#define EXTERN extern
#endif
namespace QDP
{



  EXTERN multi1d < int >lat_geom;
    QDP::SpinMatrixD cnvrt_gamma_to_qdp(int gamma);
    QDP::LatticeComplexD calculate_meson_2pt(const QDP::
                                             LatticePropagatorD & prop1,
                                             const QDP::
                                             LatticePropagatorD & prop2,
                                             int gamma_source, int gamma_sink);
  void calculate_meson_2pt_full(const QDP::LatticePropagatorD & prop1,
                                const QDP::LatticePropagatorD & prop2,
                                const QDP::SftMom & fmom, hdf5_writer & writer,
                                const char *group_name, int src_pos_ind,
                                int prop_combination);

  void observer_init(int argc, char *argv[]);
  void copy_from_openQCD(LatticeFermionD & u, const spinor_dble * psi);
  void copy_to_openQCD(const LatticeFermionD & u, spinor_dble * psi);
  void copy_to_openQCD(const multi1d < LatticeFermionD > u, spinor_dble ** psi);
  void copy_from_openQCD(LatticePropagatorD & u, const spinor_dble *const * psi);
  void copy_to_openQCD(const LatticePropagatorD & u, spinor_dble ** psi);
  void copy_from_openQCD(multi1d < LatticeColorMatrixD > &u, const su3_dble * ud);
  void copy_to_openQCD(const multi1d < LatticeColorMatrixD > &u, su3_dble * ud);
  void create_point_source(int *src_pos, LatticePropagatorD & ferm);
  void create_point_source(multi1d < int > & src_pos, LatticePropagatorD & ferm);
  void create_wuppertal_source(multi1d < int > &src_pos, LatticePropagatorD & ferm, int propid);
  void create_jacobi_source(multi1d < int > &src_pos, LatticePropagatorD & ferm, int propid);
  void create_momentum_source(multi1d < int > &src_pos, multi1d < int > &mom,
                              LatticePropagatorD & ferm);
  void create_momentum_vol_source(multi1d < int >&src_pos, multi1d < int >&mom,
                                  LatticePropagatorD & ferm);
  void calculate_propagator(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & prop);
  void calculate_propagator(int propid, int *source, LatticePropagatorD & prop);
  void calculate_propagator(int propid, multi1d < int >&source,
                            LatticePropagatorD & prop);
  void calculate_propagator(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & prop, bool);
  void calculate_propagator_once(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & prop);
  void calculate_propagator(int propid, multi1d < int >&source,
                            LatticePropagatorD & prop, bool);

  void calculate_propagators(multi1d < int > &source,
                             multi1d < LatticePropagatorD > &quark_prop,
                             bool recompute);
  void calculate_propagators(int * source,
                             multi1d < LatticePropagatorD > &quark_prop,
                             bool recompute=false);
  void apply_link_smearing(int link_smearing_idx);
  void reset_link_smearing(void);



  void twistedbc(multi1d < Double > &theta, multi1d < LatticeColorMatrixD > &uu,
                 multi1d < LatticeColorMatrixD > &uunew);




    //! Function object used for constructing the time-slice set
    class TimeSliceFunc : public SetFunc
    {
    public:
      TimeSliceFunc(int dir): dir_decay(dir) {}

      int operator() (const multi1d<int>& coordinate) const 
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 0 ;
    } else {
      return coordinate[dir_decay] ;
    }
  }

      int numSubsets() const 
  {
    if ((dir_decay<0)||(dir_decay>=Nd)) {
      return 1 ;
    } else {
      return Layout::lattSize()[dir_decay] ;
    }
  }

    private:
      TimeSliceFunc() {}  // hide default constructor

      int dir_decay;
    };

}



#endif
