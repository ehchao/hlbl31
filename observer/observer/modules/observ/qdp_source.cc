/**
* @file qdp_source.c
* \brief  Construction of source with openQCD
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-11-19
*/

#include "observer.h"

namespace QDP
{


/* --------------------------------------------------------------------------*/
/**
* \brief  Creates point source
*
* \param[in] src_pos: Cartesian coordinates for the src position: (x, y, z, t)
* \param[out] eta: Source field for a point source, i.e. 12x12 unity matrix.
*/
/* ----------------------------------------------------------------------------*/


  void create_point_source(multi1d < int > &src_pos, LatticePropagatorD & eta)
  {
    PropagatorD tmpFermion;
    ColorMatrixD tmpColorVector = zero;
    ComplexD one = QDP::cmplx(Real(1), 0);

      eta = zero;


      QDPIO::cout << "Creating point source [x,y,z,t] :\
      [" << src_pos[0] << "," << src_pos[1] << "," << src_pos[2] << "," << src_pos[3] << "]" << endl;

      start_time(create_point_source);
    for (int color = 0; color < 3; color++)
    {
      tmpFermion = zero;
      pokeColor(tmpColorVector, one, color, color);
    }
    for (int spin = 0; spin < 4; spin++)
    {
      pokeSpin(tmpFermion, tmpColorVector, spin, spin);
    }
    pokeSite(eta, tmpFermion, src_pos);
    stop_time(create_point_source);


  }

  void create_point_source(int *src_pos, LatticePropagatorD & eta)
  {
    multi1d < int >srcp(4);

    srcp[0] = src_pos[0];
    srcp[1] = src_pos[1];
    srcp[2] = src_pos[2];
    srcp[3] = src_pos[3];

    create_point_source(srcp, eta);
  }



  void create_wuppertal_source(multi1d < int > &src_pos, LatticePropagatorD & eta,
                               int propid)
  {
    PropagatorD tmpFermion;
    ColorMatrixD tmpColorVector = zero;
    ComplexD one = QDP::cmplx(Real(1), 0);
    multi1d < LatticeColorMatrixD > u(4);
    su3_dble *upoint;
    Set singlets;


    singlets.make(TimeSliceFunc(3));
    eta = zero;


    QDPIO::cout << "Creating Wuppertal source [x,y,z,t] :\
      [" << src_pos[0] << "," << src_pos[1] << "," << src_pos[2] << "," << src_pos[3] << "]" << endl;


    if (prop[propid].link_smearing_idx >= 0)
      QDPIO::cout << "Using link smearing index: " << prop[propid].
        link_smearing_idx << endl;
    if (prop[propid].smearing_steps_source > 0)
    {
      QDPIO::cout << "Number of source smearing steps: " << prop[propid].
        smearing_steps_source << endl;
      QDPIO::cout << "Source smearing parameter: " << prop[propid].
        smearing_parm_source << endl;
    }

    if (prop[propid].smearing_steps_source < 0)
      error(1,1,"qdp_sources.cc","Missing smearing parameters in infile\n");
    start_time(create_wuppertal_source);
    for (int color = 0; color < 3; color++)
    {
      tmpFermion = zero;
      pokeColor(tmpColorVector, one, color, color);
    }
    for (int spin = 0; spin < 4; spin++)
    {
      pokeSpin(tmpFermion, tmpColorVector, spin, spin);
    }

    pokeSite(eta, tmpFermion, src_pos);
    /* Smearing Source */
    if (prop[propid].link_smearing_idx >= 0)
    {
      apply_link_smearing(prop[propid].link_smearing_idx);
      u = lsmearing_stat.u_field;
    }
    else
    {
      upoint = udfld();
      copy_from_openQCD(u, upoint);
    }
    wuppertal_smear(eta, eta, u, prop[propid].smearing_parm_source,
                    prop[propid].smearing_steps_source,singlets[src_pos[3]]);

    stop_time(create_wuppertal_source);


  }

  void create_jacobi_source(multi1d < int > &src_pos, LatticePropagatorD & eta,
                               int propid)
  {
    PropagatorD tmpFermion;
    ColorMatrixD tmpColorVector = zero;
    ComplexD one = QDP::cmplx(Real(1), 0);
    multi1d < LatticeColorMatrixD > u(4);
    su3_dble *upoint;
    Set singlets;


    singlets.make(TimeSliceFunc(3));
    eta = zero;


    QDPIO::cout << "Creating Jacobi smeared source [x,y,z,t] :\
      [" << src_pos[0] << "," << src_pos[1] << "," << src_pos[2] << "," << src_pos[3] << "]" << endl;


    if (prop[propid].link_smearing_idx >= 0)
      QDPIO::cout << "Using link smearing index: " << prop[propid].
        link_smearing_idx << endl;
    if (prop[propid].smearing_steps_source > 0)
    {
      QDPIO::cout << "Number of source smearing steps: " << prop[propid].
        smearing_steps_source << endl;
      QDPIO::cout << "Source smearing parameter: " << prop[propid].
        smearing_parm_source << endl;
    }

    if (prop[propid].smearing_steps_source < 0)
      error(1,1,"qdp_sources.cc","Missing smearing parameters in infile\n");
    start_time(create_jacobi_source);
    for (int color = 0; color < 3; color++)
    {
      tmpFermion = zero;
      pokeColor(tmpColorVector, one, color, color);
    }
    for (int spin = 0; spin < 4; spin++)
    {
      pokeSpin(tmpFermion, tmpColorVector, spin, spin);
    }

    pokeSite(eta, tmpFermion, src_pos);
    /* Smearing Source */
    if (prop[propid].link_smearing_idx >= 0)
    {
      apply_link_smearing(prop[propid].link_smearing_idx);
      u = lsmearing_stat.u_field;
    }
    else
    {
      upoint = udfld();
      copy_from_openQCD(u, upoint);
    }
    jacobi_smear(eta, eta, u, prop[propid].smearing_parm_source,
                    prop[propid].smearing_steps_source,singlets[src_pos[3]]);

    stop_time(create_jacobi_source);


  }


/* --------------------------------------------------------------------------*/
/**
* \brief  Create momentum source, where the phase is \exp^(i mom \cdot x)
*
* \param src_pos Source position [x,y,z,t]
* \param mom     Momentum for momentum source [px,py,px]
* \param eta     Returns momentum source
*/
/* ----------------------------------------------------------------------------*/
  void create_momentum_source(multi1d < int >& src_pos, multi1d < int >&mom,
                              LatticePropagatorD & eta)
  {

    SftMom mom_phase(0, src_pos, mom);
    QDPIO::cout << "Creating momentum source [x,y,z,t] :\
      [" << src_pos[0] << "," << src_pos[1] << "," << src_pos[2] << "," << src_pos[3] << "]" << endl;
    QDPIO::cout << "with momentum [p_x,p_y,p_z] :\
      [" << mom[0] << "," << mom[1] << "," << mom[2] << "]" << endl;

    create_point_source(src_pos, eta);
    /* Multiplying lots of zeros butmost compact way of implementation */
    eta *= mom_phase[0];


  }
/* --------------------------------------------------------------------------*/
/**
* \brief  Create momentum volume source, where the phase is \exp^(i mom \cdot x)
*
* \param src_pos Source position
* \param mom     Momentum for momentum source [px,py,px]
* \param eta     Returns momentum source
*/
/* ----------------------------------------------------------------------------*/
  void create_momentum_vol_source(multi1d < int >& src_pos, multi1d < int >&mom,
                                  LatticePropagatorD & eta)
  {

    PropagatorD tmpFermion;
    ColorMatrixD tmpColorVector = zero;
    ComplexD one = QDP::cmplx(Real(1), 0);
    SftMom mom_phase(0, src_pos, mom);

    eta = zero;


    QDPIO::cout << "Creating momentum wall source with momentum [p_x,p_y,p_z] :\
      [" << mom[0] << "," << mom[1] << "," << mom[2] << "]" << endl;

    for (int color = 0; color < 3; color++)
    {
      tmpFermion = zero;
      pokeColor(tmpColorVector, one, color, color);
    }
    for (int spin = 0; spin < 4; spin++)
    {
      pokeSpin(tmpFermion, tmpColorVector, spin, spin);
    }

    eta = tmpFermion;
    eta *= mom_phase[0];


  }
}
