/**
* @file propagators.c
* \brief  Calculate the quark propagators as provided in the infile
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-05-14
*
* v0.1: Initial Version.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#include "flags.h"
#include "linalg.h"
#include "sap.h"
#include "dfl.h"
#include "dirac.h"
#include "forces.h"
#include "global.h"
#include "observ.h"
#include "archive.h"

#define PROP_DEBUG 0

static spinor_dble *psi, *phi[12];

void solve_dirac_comp(spinor_dble * eta, spinor_dble * wpsi, int *status,
                 int idx_prop, int comp)
{
  solver_parms_t sp;
  sap_parms_t sap;
  double mu;
  double residue;
  int mins0;
/*  double norm_rho;
  spinor_dble *re;
  re = reserve_wsd(1)[0];*/

  if (prop[idx_prop].truncated_solver)
    mins0 = -1;
  else
    mins0 = 0;
  mu = prop[idx_prop].mu;
  set_sw_parms(1 / (2. * prop[idx_prop].kappa) - 4.);
  sp = solver_parms(prop[idx_prop].idx_solver);
  if (sp.solver == CGNE)
  {
    mulg5_dble(VOLUME, eta);

    tmcg(sp.nmx, sp.res, mus, eta, eta, status);

    error_root(status[0] < mins0, 1, "solve_dirac [ms4.c]",
               "CGNE solver failed (status = %d)", status[0]);

    Dw_dble(-mus, eta, wpsi);
    mulg5_dble(VOLUME, wpsi);
  }
  else if (sp.solver == SAP_GCR)
  {
    sap = sap_parms();
    set_sap_parms(sap.bs, sp.isolv, sp.nmr, sp.ncy);

    sap_gcr(sp.nkv, sp.nmx, sp.res, mus, eta, wpsi, status);

    error_root(status[0] < mins0, 1, "solve_dirac [ms4.c]",
               "SAP_GCR solver failed (status = %d)", status[0]);
  }
  else if (sp.solver == DFL_SAP_GCR)
  {
    sap = sap_parms();
    set_sap_parms(sap.bs, sp.isolv, sp.nmr, sp.ncy);

    residue=dfl_sap_gcr2(sp.nkv, sp.nmx, sp.res, mus, eta, wpsi, status);
    message("Propagator %d component %d : Status [ %d , %d , %d ] Residue %.15e Norm square of solution is %.15e Relative residue %.15e\n",idx_prop,comp,status[0],status[1],status[2],residue,norm_square_dble(VOLUME,1,wpsi),residue/sqrt(norm_square_dble(VOLUME,1,eta)));

    error_root((status[0] < mins0) || (status[1] < 0), 1,
               "solve_dirac [ms4.c]", "DFL_SAP_GCR solver failed "
               "(status = %d,%d,%d)", status[0], status[1], status[2]);
  }
  else
    error_root(1, 1, "solve_dirac [ms4.c] ", "Unknown or unsupported solver");


  /*Dw_dble(prop[idx_prop].mu, wpsi, re);
  mulr_spinor_add_dble(VOLUME, re, eta, (double) -1.f);
  residue = norm_square_dble(VOLUME, 1, re);
  norm_rho=norm_square_dble(VOLUME,1,eta);
  
  message("Residue is %e \n", sqrt(residue/norm_rho));
  release_wsd();*/
}

/* --------------------------------------------------------------------------*/
/**
* \brief  This will solve the Dirac equation.
*
* \param eta
* \param wpsi
* \param status
* \param idx_prop
*/
/* ----------------------------------------------------------------------------*/
void solve_dirac(spinor_dble * eta, spinor_dble * wpsi, int *status,
                 int idx_prop)
{
  solver_parms_t sp;
  sap_parms_t sap;
  double mu;
  double residue;
  int mins0;
/*  double norm_rho;
  spinor_dble *re;
  re = reserve_wsd(1)[0];*/

  if (prop[idx_prop].truncated_solver)
    mins0 = -1;
  else
    mins0 = 0;
  mu = prop[idx_prop].mu;
  set_sw_parms(1 / (2. * prop[idx_prop].kappa) - 4.);
  sp = solver_parms(prop[idx_prop].idx_solver);
  if (sp.solver == CGNE)
  {
    mulg5_dble(VOLUME, eta);

    tmcg(sp.nmx, sp.res, mus, eta, eta, status);

    error_root(status[0] < mins0, 1, "solve_dirac [ms4.c]",
               "CGNE solver failed (status = %d)", status[0]);

    Dw_dble(-mus, eta, wpsi);
    mulg5_dble(VOLUME, wpsi);
  }
  else if (sp.solver == SAP_GCR)
  {
    sap = sap_parms();
    set_sap_parms(sap.bs, sp.isolv, sp.nmr, sp.ncy);

    sap_gcr(sp.nkv, sp.nmx, sp.res, mus, eta, wpsi, status);

    error_root(status[0] < mins0, 1, "solve_dirac [ms4.c]",
               "SAP_GCR solver failed (status = %d)", status[0]);
  }
  else if (sp.solver == DFL_SAP_GCR)
  {
    sap = sap_parms();
    set_sap_parms(sap.bs, sp.isolv, sp.nmr, sp.ncy);

    residue=dfl_sap_gcr2(sp.nkv, sp.nmx, sp.res, mus, eta, wpsi, status);
    message("Propagator %d : Status [ %d , %d , %d ] Residue %e\n",idx_prop,status[0],status[1],status[2],residue);

    error_root((status[0] < mins0) || (status[1] < 0), 1,
               "solve_dirac [ms4.c]", "DFL_SAP_GCR solver failed "
               "(status = %d,%d,%d)", status[0], status[1], status[2]);
  }
  else
    error_root(1, 1, "solve_dirac [ms4.c] ", "Unknown or unsupported solver");


  /*Dw_dble(prop[idx_prop].mu, wpsi, re);
  mulr_spinor_add_dble(VOLUME, re, eta, (double) -1.f);
  residue = norm_square_dble(VOLUME, 1, re);
  norm_rho=norm_square_dble(VOLUME,1,eta);
  
  message("Residue is %e \n", sqrt(residue/norm_rho));
  release_wsd();*/
}

spinor_dble *propagator(int idx_prop, int isrc)
{
  int l, stat[3];
  double wt1, wt2;
  double residue;
  double norm_rho;
  spinor_dble *eta;
  spinor_dble *re;

  set_sw_parms(1 / (2. * prop[idx_prop].kappa) - 4.);
  alloc_prop(idx_prop, isrc);
  psi = observer.prop_list[idx_prop].sd[isrc];

  MPI_Barrier(MPI_COMM_WORLD);
  wt1 = MPI_Wtime();
  eta = reserve_wsd(1)[0];
  re = reserve_wsd(1)[0];

  if ((observer.prop_list[idx_prop].spin_idx_status[isrc] != COMPUTED))
  {
    for (l = 0; l < 3; l++)
    {
      stat[l] = 0;
    }
    create_source(prop[idx_prop].src_t, prop[idx_prop].src_pos, isrc / 3,
                  isrc % 3, eta);
    solve_dirac(eta, psi, stat, idx_prop);

    observer.prop_list[idx_prop].spin_idx_status[isrc] = COMPUTED;
  }
  else if ((observer.prop_list[idx_prop].spin_idx_status[isrc] == COMPUTED))
  {
    if (PROP_DEBUG > 0)
      message("Propagator already computed. Using old value\n");
  }
  else
  {
    error(1, 1, "[propagator.c]", "Not implemented yet\n");
  }

  wt2 = MPI_Wtime();

  Dw_dble(prop[idx_prop].mu, psi, re);
  mulr_spinor_add_dble(VOLUME, re, eta, (double) -1.f);
  residue = norm_square_dble(VOLUME, 1, re);
  norm_rho=norm_square_dble(VOLUME,1,eta);
  MPI_Barrier(MPI_COMM_WORLD);
  if (PROP_DEBUG >= 0)
  {
    message
      ("Calcuation of propagator %d , src pos = (%d,%d,%d,%d), component (spin,color) = (%d,%d)\n",
       idx_prop, prop[idx_prop].src_pos[0], prop[idx_prop].src_pos[1],
       prop[idx_prop].src_pos[2], prop[idx_prop].src_pos[3], isrc / 3,
       isrc % 3);
    mpi_print_stats(wt1, wt2, "Propagator Calculation");
    fflush(stdout);
  }
  release_wsd();
  release_wsd();
  return psi;
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Returns pointer to the propgator described in the infile by index idx_prop.
*
* \param idx_prop Index of propgator as given in the infile
*
* \return   Retunrs pointer to propgator.
*/
/* ----------------------------------------------------------------------------*/
spinor_dble **full_propagator(int idx_prop)
{
  int isrc;
  double wt1, wt2;

  wt1 = MPI_Wtime();

  for (isrc = 0; isrc < 12; isrc++)
  {
    phi[isrc] = propagator(idx_prop, isrc);
  }

  wt2 = MPI_Wtime();

  mpi_print_stats(wt1, wt2, "Calculation Fullpropagator");
  return phi;
}


void set_source_position(int src_pos[4])
{
  int i;

  for (i = 0; i < observer.no_prop; i++)
  {
    prop[i].src_pos[0] = src_pos[0];
    prop[i].src_pos[1] = src_pos[1];
    prop[i].src_pos[2] = src_pos[2];
    prop[i].src_pos[3] = src_pos[3];
  }


}

void set_isource_position(int src_pos[4], int i)
{
  prop[i].src_pos[0] = src_pos[0];
  prop[i].src_pos[1] = src_pos[1];
  prop[i].src_pos[2] = src_pos[2];
  prop[i].src_pos[3] = src_pos[3];
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Sets flag of all propagtor to be recomputed on next call.
*/
/* ----------------------------------------------------------------------------*/
void recompute_propagators(void)
{
  int i, j;

  for (i = 0; i < observer.no_prop; i++)
  {
    for (j = 0; j < 12; j++)
    {
      if(observer.prop_list[i].spin_idx_status[j] == INIT)
        continue;
      observer.prop_list[i].spin_idx_status[j] = ALLOCATED;
    }
  }
}

static void poke_colour_column(su3_dble * r, su3_vector_dble * l, int column)
{
  switch (column)
  {
  case 0:
    (*r).c11.re = (*l).c1.re;
    (*r).c11.im = (*l).c1.im;
    (*r).c21.re = (*l).c2.re;
    (*r).c21.im = (*l).c2.im;
    (*r).c31.re = (*l).c3.re;
    (*r).c31.im = (*l).c3.im;
    break;
  case 1:
    (*r).c12.re = (*l).c1.re;
    (*r).c12.im = (*l).c1.im;
    (*r).c22.re = (*l).c2.re;
    (*r).c22.im = (*l).c2.im;
    (*r).c32.re = (*l).c3.re;
    (*r).c32.im = (*l).c3.im;
    break;
  case 2:
    (*r).c13.re = (*l).c1.re;
    (*r).c13.im = (*l).c1.im;
    (*r).c23.re = (*l).c2.re;
    (*r).c23.im = (*l).c2.im;
    (*r).c33.re = (*l).c3.re;
    (*r).c33.im = (*l).c3.im;
    break;

  }
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Legacy conversion into format of old measure code.
*
* \param volume size of the Volume, used to determine the number of spinors
* \param eta Source spinor, i.e. has to be an array of 12 pointer to the 
* solution of the Dirac Equation
* \param fpsi Full_spinor_dble in the old format.
*/
/* ----------------------------------------------------------------------------*/
void cnvrt_spinor_to_full_spinor(int volume, spinor_dble ** eta,
                                 full_spinor_dble * fpsi)
{
  int id, ic, is, i = 0;
  spinor_dble *rpk, *rpe;


  for (is = 0; is < 12; is++)
  {
    id = is / 3;
    ic = is % 3;
    rpk = eta[is];
    rpe = eta[is] + volume;
    i = 0;
    for (; rpk < rpe; rpk++)
    {
      switch (id)
      {
      case 0:
        poke_colour_column(&((*(fpsi + i)).c11), &((*rpk).c1), ic);
        poke_colour_column(&((*(fpsi + i)).c21), &((*rpk).c2), ic);
        poke_colour_column(&((*(fpsi + i)).c31), &((*rpk).c3), ic);
        poke_colour_column(&((*(fpsi + i)).c41), &((*rpk).c4), ic);
        break;
      case 1:
        poke_colour_column(&((*(fpsi + i)).c12), &((*rpk).c1), ic);
        poke_colour_column(&((*(fpsi + i)).c22), &((*rpk).c2), ic);
        poke_colour_column(&((*(fpsi + i)).c32), &((*rpk).c3), ic);
        poke_colour_column(&((*(fpsi + i)).c42), &((*rpk).c4), ic);
        break;
      case 2:
        poke_colour_column(&((*(fpsi + i)).c13), &((*rpk).c1), ic);
        poke_colour_column(&((*(fpsi + i)).c23), &((*rpk).c2), ic);
        poke_colour_column(&((*(fpsi + i)).c33), &((*rpk).c3), ic);
        poke_colour_column(&((*(fpsi + i)).c43), &((*rpk).c4), ic);
        break;
      case 3:
        poke_colour_column(&((*(fpsi + i)).c14), &((*rpk).c1), ic);
        poke_colour_column(&((*(fpsi + i)).c24), &((*rpk).c2), ic);
        poke_colour_column(&((*(fpsi + i)).c34), &((*rpk).c3), ic);
        poke_colour_column(&((*(fpsi + i)).c44), &((*rpk).c4), ic);
        break;
      }
      i += 1;
    }
    eta[is] = rpe - volume;
  }
}
