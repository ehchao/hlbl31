/*#include <qdp.h>
#include <qdpinterface.h>
extern "C"
{
#include "global.h"
#include "flags.h"
#include "lattice.h"
#include "su3.h"
#include "apps.h"
#include "observ.h"
}*/

#include "observer.h"

#define PROP_DEBUG 0

namespace QDP
{
  static multi1d < int >prop_map;

  void observer_init(int argc, char *argv[])
  {


    QDP_initialize(&argc, &argv);       /* Initializes QDP */

    lat_geom.resize(4);
    init_arrays(argc, argv);
    read_infile(argc, argv);
    geometry();
    alloc_wspace(0, 12, 0, 0);

    multi1d < int >nrow(Nd);
      lat_geom[0] = NPROC1 * L1;
      lat_geom[1] = NPROC2 * L2;
      lat_geom[2] = NPROC3 * L3;
      lat_geom[3] = NPROC0 * L0;
      print_info();
      Layout::setLattSize(lat_geom);
      Layout::create();
    prop_map.resize(observer.no_prop);
    for(int i=0;i<observer.no_prop;i++)
        prop_map[i]=-1;
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Converts \gamma matrices from old measure code convention to QDP.
* Conversion used from openQCD -> QDP
* \gamma_0 -> \gamma_3
* \gamma_1 -> \gamma_0
* \gamma_2 -> \gamma_1
* \gamma_3 -> \gamma_2
*
* \param gamma: Integer representing \gamma in measurecode convention
*
* \return   Spinmatrix object representing \gamma in QDP 
*/
/* ----------------------------------------------------------------------------*/
  QDP::SpinMatrixD cnvrt_gamma_to_qdp(int gamma)
  {

    SpinMatrixD g;
    SpinMatrixD mone = -1.;
    SpinMatrixD one = 1.;

    switch (gamma)
    {
    case G0:
      g = one * Gamma(8);
      break;
    case G1:
      g = one * Gamma(1);
      break;
    case G2:
      g = one * Gamma(2);
      break;
    case G3:
      g = one * Gamma(4);
      break;
    case G5:
      g = mone * Gamma(15);
      break;
    case G0G5:
      g = one * Gamma(7);
      break;
    case G1G5:
      g = mone * Gamma(14);
      break;
    case G2G5:
      g = one * Gamma(13);
      break;
    case G3G5:
      g = mone * Gamma(11);
      break;
    case G0G1:
      g = mone * Gamma(9);
      break;
    case G0G2:
      g = mone * Gamma(10);
      break;
    case G0G3:
      g = mone * Gamma(12);
      break;
    case G1G2:
      g = one * Gamma(3);
      break;
    case G1G3:
      g = one * Gamma(5);
      break;
    case G2G3:
      g = one * Gamma(6);
      break;
    case Id:
      g = one;
      break;
    default:
      error(1, 1, "qdp_utils.c", "Unkown Gamma Matrix");
    }

    return g;
  }


  void reset_link_smearing(int link_smearing_idx)
  {
    if (smearing_parms[link_smearing_idx].stat.status == COMPUTED)
    {
      smearing_parms[link_smearing_idx].stat.status = INIT;
    }

  }


  static void calculate_link_smearing(int link_smearing_idx,
                                      multi1d < LatticeColorMatrixD > &up)
  {

    switch (cnvrt_smear(smearing_parms[link_smearing_idx].smearing_t))
    {
    case APE:
      ape_smear(lsmearing_stat.u_field, up,
                smearing_parms[link_smearing_idx].smearing_parm,
                smearing_parms[link_smearing_idx].smearing_steps,
                smearing_parms[link_smearing_idx].skipdir);

      break;
    default:
      error(1, 1, "qdp_utils.c", "Unknown smearing type");
    }

  }

  void apply_link_smearing(int link_smearing_idx)
  {
    multi1d < LatticeColorMatrixD > up(4);
    su3_dble *upoint;
    int sign;

    start_time(apply_link_smearing);
    upoint = udfld();

    /* Change to periodic boundary */
    sign = chs_ubnd(1);

    if (lsmearing_stat.idx != link_smearing_idx)
    {
      if (lsmearing_stat.u_field.size() != 4)
        lsmearing_stat.u_field.resize(4);
      copy_from_openQCD(up, upoint);
      calculate_link_smearing(link_smearing_idx, up);
      smearing_parms[link_smearing_idx].stat.status = COMPUTED;
      lsmearing_stat.idx = link_smearing_idx;
    }
    else
    {
      if (smearing_parms[link_smearing_idx].stat.status == COMPUTED)
      {
        QDPIO::cout << "Link smearing already applied not doing anything!" <<
          endl;
      }
      else
      {
        QDPIO::cout << "Link smearing already applied not doing anything!" <<
          endl;
        copy_from_openQCD(up, upoint);
        calculate_link_smearing(link_smearing_idx, up);
      }
    }

    /* Change to old value */
    if (sign == 1)
      chs_ubnd(-1);
    stop_time(apply_link_smearing);
  }



  void calculate_propagator(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & quark_prop, bool recompute)
  {
    if (recompute)
    {
      recompute_propagators();
      calculate_propagator(propid, source, quark_prop);
    }
    else
    {
      calculate_propagator(propid, source, quark_prop);
    }

  }

  void calculate_propagator(int propid, multi1d < int >&src_pos,
                            LatticePropagatorD & quark_prop)
  {
    calculate_propagator(propid, src_pos, quark_prop, false);
  };

  void calculate_propagator(int propid, int *src_pos,
                            LatticePropagatorD & quark_prop)
  {
    multi1d < int >pos(4);

    pos[0] = src_pos[0];
    pos[1] = src_pos[1];
    pos[2] = src_pos[2];
    pos[3] = src_pos[3];

    calculate_propagator(propid, pos, quark_prop, false);
  };

  void calculate_propagator(int propid, multi1d < int >&src_pos,
                            LatticePropagatorD & quark_prop, bool recompute)
  {
    LatticePropagatorD source;
    su3_dble *upoint;
    multi1d < LatticeColorMatrixD > u(4);
    int noc = 0;
    int i;
    spinor_dble *psi[12];

    start_time(calculate_propagator);


    if (prop_map[propid] > -1)
    {
      QDPIO::cout << "Reusing solution of [Propagator " << prop_map[propid] <<
        "] for [Propagator " << propid << "]" << endl;
      /* Check if it has been calculated already */
      /* This check should never fail */
      for (i = 0; i < 12; i++)
      {
        alloc_prop(propid, i);
        psi[i] = observer.prop_list[prop_map[propid]].sd[i];
        if ((observer.prop_list[prop_map[propid]].spin_idx_status[i] ==
             COMPUTED))
        {
          noc++;
        }
      }

      //error(noc != 12, 1, "qdp_utils.cc",
      //      "Wanted to reuse uncomputed propagator");

      start_time(copy_solution);
      copy_from_openQCD(quark_prop, psi);
      stop_time(copy_solution);
    }
    else
    {
      /* Prepare Source */
      switch (cnvrt_src(prop[propid].src_t))
      {
      case POINT:
        create_point_source(src_pos, source);
        if (prop[propid].link_smearing_idx >= 0)
        {
          message("Link smearing for point source is ignored!\n");
        }
        break;
      case WUPPERTAL:
        create_wuppertal_source(src_pos, source, propid);
        break;
      case JACOBI:
        create_jacobi_source(src_pos, source, propid);
        break;
      default:
        error(1, 1, "qdp_utils.cc", "Source type %s not defined",
              prop[propid].src_t);
        break;
      }

      if (recompute)
      {
        recompute_propagators();
        calculate_propagator(propid, source, quark_prop);
      }
      else
      {
        calculate_propagator(propid, source, quark_prop);
      }
    }

    /* Link smearing */
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
    /* Smear the sink */
    switch (cnvrt_src(prop[propid].src_t))
    {
    case POINT:
      break;
    case WUPPERTAL:
      if (prop[propid].smearing_steps_sink < 0)
      {
        /*  error(1, 1, "qdp_utils.cc",
           "Smearing parameters for sink are missing!\n");
         */
        QDPIO::cout << "No sink smearing applied." << endl;
      }
      else
      {
        QDPIO::cout << "Number of sink smearing steps: " <<
          prop[propid].smearing_steps_sink << endl;
        QDPIO::cout << "Sink smearing parameter: " << prop[propid].
          smearing_parm_sink << endl;
        wuppertal_smear(quark_prop, quark_prop, u,
                        prop[propid].smearing_parm_sink,
                        prop[propid].smearing_steps_sink);
      }

      break;
    case JACOBI:
      if (prop[propid].smearing_steps_sink < 0)
      {
        /*  error(1, 1, "qdp_utils.cc",
           "Smearing parameters for sink are missing!\n");
         */
        QDPIO::cout << "No sink smearing applied." << endl;
      }
      else
      {
        QDPIO::cout << "Number of sink smearing steps: " <<
          prop[propid].smearing_steps_sink << endl;
        QDPIO::cout << "Sink smearing parameter: " << prop[propid].
          smearing_parm_sink << endl;
        jacobi_smear(quark_prop, quark_prop, u,
                        prop[propid].smearing_parm_sink,
                        prop[propid].smearing_steps_sink);
      }

      break;
    default:
      error(1, 1, "qdp_utils.c", "Source type %s not defined",
            prop[propid].src_t);
      break;
    }

    stop_time(calculate_propagator);

  }

  void calculate_propagator(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & quark_prop)
  {
    spinor_dble *psi[12];
    spinor_dble **eta;
    int i, l;
    int noc;
    int status[3];
    StopWatch swatch;

    noc = 0;
    eta = reserve_wsd(12);
    start_time(full_solve);
    for (i = 0; i < 12; i++)
    {
      alloc_prop(propid, i);
      psi[i] = observer.prop_list[propid].sd[i];
      if ((observer.prop_list[propid].spin_idx_status[i] == COMPUTED))
      {
        noc++;
      }
    }

    if (noc == 12)
    {
      message("Propagator already computed. Using old value\n");

    }
    else
    {
      start_time(copy_source);
      copy_to_openQCD(source, eta);
      stop_time(copy_source);

      int chs_query = chs_ubnd(-1);
      if (chs_query)
        QDPIO::cout << "Changed sign of gauge fields on boundary to -1 to implement antiperiodic bcs for solve." << endl;

      if ((query_flags(SWD_UP2DATE) != 1))
      {
        QDPIO::cout << "Recompute propagators SWD not up to date." << endl;
        recompute_propagators();
      }

      start_time(calculate_full_propagator);
      for (i = 0; i < 12; i++)
      {
        if ((observer.prop_list[propid].spin_idx_status[i] != COMPUTED))
        {
          for (l = 0; l < 3; l++)
          {
            status[l] = 0;
          }
          start_time(propagator_component);
          solve_dirac_comp(eta[i], psi[i], status, propid,i);
          stop_time(propagator_component);
          observer.prop_list[propid].spin_idx_status[i] = COMPUTED;
        }
        else if ((observer.prop_list[propid].spin_idx_status[i] == COMPUTED))
        {
          if (PROP_DEBUG > 0)
            message("Propagator already computed. Using old value\n");
        }
      }

      
      if (chs_query == 1)
      {
        chs_query = chs_ubnd(1);
        QDPIO::cout << "Returned sign of gauge fields to +1 after solve." << endl;
      }
      stop_time(calculate_full_propagator);
    }

    start_time(copy_solution);
    copy_from_openQCD(quark_prop, psi);
    stop_time(copy_solution);
    release_wsd();
    stop_time(full_solve);

  }

  void calculate_propagator_once(int propid, LatticePropagatorD & source,
                            LatticePropagatorD & quark_prop)
  {
    spinor_dble **psi;
    spinor_dble **eta;
    int i;
    int status[3];
    StopWatch swatch;
    int sign;


    /* Change to antiperiodic boundary */
    sign = chs_ubnd(-1);

    psi = reserve_wsd(12);
    eta = reserve_wsd(12);

    start_time(full_solve);
    start_time(copy_source);
    copy_to_openQCD(source, eta);
    stop_time(copy_source);


      start_time(calculate_full_propagator);
      for(i=0;i<12;i++)
      {
      start_time(propagator_component);
      solve_dirac_comp(eta[i], psi[i], status, propid,i);
      stop_time(propagator_component);
      }
      stop_time(calculate_full_propagator);

    start_time(copy_solution);
    copy_from_openQCD(quark_prop, psi);
    stop_time(copy_solution);
    release_wsd();
    release_wsd();
    if (sign == 1)
      chs_ubnd(1);
    stop_time(full_solve);

  }

  static int max_int(void)
  {
    int lmax = 0;
    for (int i = 0; i < observer.no_prop; i++)
    {
      if (observer.idx_prop[i] > lmax)
        lmax = observer.idx_prop[i];
    }

    return lmax;
  }

  void calculate_propagators(int *source,
                             multi1d < LatticePropagatorD > &quark_prop,
                             bool recompute)
  {
    multi1d < int >src(4);

    src[0] = source[0];
    src[1] = source[1];
    src[2] = source[2];
    src[3] = source[3];

    calculate_propagators(src, quark_prop, recompute);

  }

  int wronce = 0;

  int identify_props(prop_parms_t * test, multi1d < int >&prop_map)
  {
    int i, j;
    int num = 0;
    int loop;

    loop = prop_map.size();
    for (i = 0; i < loop; i++)
    {
      prop_map[i] = -1;
    }
    for (i = 0; i < loop; i++)
    {
      for (j = i + 1; j < loop; j++)
      {
        if (!compare_structs_simple(test + i, test + j) && (prop_map[j] < 0))
        {
          prop_map[j] = i;
          num++;
          QDPIO::cout<<"prop_map["<<j<<"]->"<<prop_map[j]<<endl;
        }
      }

    }

    if (wronce == 0)
    {
      QDPIO::cout << "Found " << num << " propagators to differ after solve." <<
        endl;
      for (i = 0; i < loop; i++)
      {
        if (prop_map[i] > -1)
        {
          QDPIO::cout << "[Propagator " << i << "] is identical to [Propagator"
            << prop_map[i] << "] before solve" << endl;
        }
      }
      wronce++;

    }
    return num;
  }

  int identify_props_observer(prop_parms_t * t, multi1d < int >&prop_map)
  {

    prop_parms_t *tmp_prop_list;
    int i, j;
    int num = 0;
    int loop;

    tmp_prop_list =
      (prop_parms_t *) malloc(observer.no_prop * sizeof(prop_parms_t));

    for (int propid = 0; propid < observer.no_prop; propid++)
    {
      tmp_prop_list[propid] = prop[observer.idx_prop[propid]];
    }

    loop = prop_map.size();
    for (i = 0; i < loop; i++)
    {
      prop_map[i] = -1;
    }
    for (i = 0; i < loop; i++)
    {
      for (j = i + 1; j < loop; j++)
      {
        if (!compare_structs_simple(tmp_prop_list + i, tmp_prop_list + j)
            && (prop_map[j] < 0))
        {
          prop_map[observer.idx_prop[j]] = observer.idx_prop[i];
          num++;
        }
      }

    }

    if (wronce == 0)
    {
      QDPIO::cout << "Found " << num <<
        " propagators to differ only after solve." << endl;
      for (i = 0; i < loop; i++)
      {
        if (prop_map[i] > -1)
        {
          QDPIO::cout << "[Propagator " << i << "] is identical to [Propagator "
            << prop_map[i] << "] before solve." << endl;
        }
      }
      wronce++;

    }
    return num;

  }

  void calculate_propagators(multi1d < int >&source,
                             multi1d < LatticePropagatorD > &quark_prop,
                             bool recompute = false)
  {

    int lmax;
    lmax = max_int();
    if ((lmax + 1) != observer.no_prop)
    {
      error(1, 1, "qdp_utils.cc",
            "Maximum propid must be equal to( observer.no_prop-1) [propid,no_prop]: [%d,%d]\n",
            lmax, observer.no_prop);
    }

    /* Check which propagators differ only after the propagtors */
    //if (QDP::Layout::nodeNumber() == 0)
    //{
    identify_props_observer(prop, prop_map);
    //}

    //QDPInternal::broadcast(prop_map);
    if (quark_prop.size() != observer.no_prop)
    {
      quark_prop.resize(observer.no_prop);
    }
    for (int propid = 0; propid < observer.no_prop; propid++)
    {
      if (recompute)
      {
        recompute_propagators();
        calculate_propagator(observer.idx_prop[propid], source,
                             quark_prop[observer.idx_prop[propid]]);
      }
      else
      {
        calculate_propagator(observer.idx_prop[propid], source,
                             quark_prop[observer.idx_prop[propid]]);
      }
    }
  }

  void import_config(int cnfg_no)
  {
    char import_cnfg_cname[NAME_SIZE];
    snprintf(import_cnfg_cname, sizeof(import_cnfg_cname), "%s/%sn%d", cnfg_dir,
             nbase, cnfg_no);

    start_time(import_cnfg);
    import_cnfg(import_cnfg_cname);
    stop_time(import_cnfg);
  }

}
