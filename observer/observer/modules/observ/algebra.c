#define ALGEBRA_C

#include "su3.h"
#include "global.h"
#include "observ.h"


#define VOL_SPINOR_SIZE VOLUME


static void sign_flip(spinor_dble * eta, spinor_dble * rho)
{
  rho->c1.c1.re = -eta->c1.c1.re;
  rho->c1.c1.im = -eta->c1.c1.im;
  rho->c1.c2.re = -eta->c1.c2.re;
  rho->c1.c2.im = -eta->c1.c2.im;
  rho->c1.c3.re = -eta->c1.c3.re;
  rho->c1.c3.im = -eta->c1.c3.im;

  rho->c2.c1.re = -eta->c2.c1.re;
  rho->c2.c1.im = -eta->c2.c1.im;
  rho->c2.c2.re = -eta->c2.c2.re;
  rho->c2.c2.im = -eta->c2.c2.im;
  rho->c2.c3.re = -eta->c2.c3.re;
  rho->c2.c3.im = -eta->c2.c3.im;

  rho->c3.c1.re = -eta->c3.c1.re;
  rho->c3.c1.im = -eta->c3.c1.im;
  rho->c3.c2.re = -eta->c3.c2.re;
  rho->c3.c2.im = -eta->c3.c2.im;
  rho->c3.c3.re = -eta->c3.c3.re;
  rho->c3.c3.im = -eta->c3.c3.im;

  rho->c4.c1.re = -eta->c4.c1.re;
  rho->c4.c1.im = -eta->c4.c1.im;
  rho->c4.c2.re = -eta->c4.c2.re;
  rho->c4.c2.im = -eta->c4.c2.im;
  rho->c4.c3.re = -eta->c4.c3.re;
  rho->c4.c3.im = -eta->c4.c3.im;
}

static void sp_assign(spinor_dble * eta, spinor_dble * rho)
{
  rho->c1.c1.re = eta->c1.c1.re;
  rho->c1.c1.im = eta->c1.c1.im;
  rho->c1.c2.re = eta->c1.c2.re;
  rho->c1.c2.im = eta->c1.c2.im;
  rho->c1.c3.re = eta->c1.c3.re;
  rho->c1.c3.im = eta->c1.c3.im;

  rho->c2.c1.re = eta->c2.c1.re;
  rho->c2.c1.im = eta->c2.c1.im;
  rho->c2.c2.re = eta->c2.c2.re;
  rho->c2.c2.im = eta->c2.c2.im;
  rho->c2.c3.re = eta->c2.c3.re;
  rho->c2.c3.im = eta->c2.c3.im;

  rho->c3.c1.re = eta->c3.c1.re;
  rho->c3.c1.im = eta->c3.c1.im;
  rho->c3.c2.re = eta->c3.c2.re;
  rho->c3.c2.im = eta->c3.c2.im;
  rho->c3.c3.re = eta->c3.c3.re;
  rho->c3.c3.im = eta->c3.c3.im;

  rho->c4.c1.re = eta->c4.c1.re;
  rho->c4.c1.im = eta->c4.c1.im;
  rho->c4.c2.re = eta->c4.c2.re;
  rho->c4.c2.im = eta->c4.c2.im;
  rho->c4.c3.re = eta->c4.c3.re;
  rho->c4.c3.im = eta->c4.c3.im;
}

static void sp_cmplx_massign(spinor_dble * eta, spinor_dble * rho)
{
  rho->c1.c1.re = eta->c1.c1.im;
  rho->c1.c1.im = -eta->c1.c1.re;
  rho->c1.c2.re = eta->c1.c2.im;
  rho->c1.c2.im = -eta->c1.c2.re;
  rho->c1.c3.re = eta->c1.c3.im;
  rho->c1.c3.im = -eta->c1.c3.re;

  rho->c2.c1.re = eta->c2.c1.im;
  rho->c2.c1.im = -eta->c2.c1.re;
  rho->c2.c2.re = eta->c2.c2.im;
  rho->c2.c2.im = -eta->c2.c2.re;
  rho->c2.c3.re = eta->c2.c3.im;
  rho->c2.c3.im = -eta->c2.c3.re;

  rho->c3.c1.re = eta->c3.c1.im;
  rho->c3.c1.im = -eta->c3.c1.re;
  rho->c3.c2.re = eta->c3.c2.im;
  rho->c3.c2.im = -eta->c3.c2.re;
  rho->c3.c3.re = eta->c3.c3.im;
  rho->c3.c3.im = -eta->c3.c3.re;

  rho->c4.c1.re = eta->c4.c1.im;
  rho->c4.c1.im = -eta->c4.c1.re;
  rho->c4.c2.re = eta->c4.c2.im;
  rho->c4.c2.im = -eta->c4.c2.re;
  rho->c4.c3.re = eta->c4.c3.im;
  rho->c4.c3.im = -eta->c4.c3.re;
}

static void sp_cmplx_assign(spinor_dble * eta, spinor_dble * rho)
{
  rho->c1.c1.re = -eta->c1.c1.im;
  rho->c1.c1.im = eta->c1.c1.re;
  rho->c1.c2.re = -eta->c1.c2.im;
  rho->c1.c2.im = eta->c1.c2.re;
  rho->c1.c3.re = -eta->c1.c3.im;
  rho->c1.c3.im = eta->c1.c3.re;

  rho->c2.c1.re = -eta->c2.c1.im;
  rho->c2.c1.im = eta->c2.c1.re;
  rho->c2.c2.re = -eta->c2.c2.im;
  rho->c2.c2.im = eta->c2.c2.re;
  rho->c2.c3.re = -eta->c2.c3.im;
  rho->c2.c3.im = eta->c2.c3.re;

  rho->c3.c1.re = -eta->c3.c1.im;
  rho->c3.c1.im = eta->c3.c1.re;
  rho->c3.c2.re = -eta->c3.c2.im;
  rho->c3.c2.im = eta->c3.c2.re;
  rho->c3.c3.re = -eta->c3.c3.im;
  rho->c3.c3.im = eta->c3.c3.re;

  rho->c4.c1.re = -eta->c4.c1.im;
  rho->c4.c1.im = eta->c4.c1.re;
  rho->c4.c2.re = -eta->c4.c2.im;
  rho->c4.c2.im = eta->c4.c2.re;
  rho->c4.c3.re = -eta->c4.c3.im;
  rho->c4.c3.im = eta->c4.c3.re;
}

static void spinor_assign(int vol, spinor_dble * eta, spinor_dble * rho)
{
  spinor_dble *sd;
  int i = 0;
  sd = eta + vol;

  for (; eta < sd; eta++)
  {
    sp_assign(eta, rho + i);
    i++;
  }
  eta = sd - vol;

}

static void spinor_cmplx_massign(int vol, spinor_dble * eta,
                                 spinor_dble * rho)
{
  spinor_dble *sd;
  int i = 0;
  sd = eta + vol;

  for (; eta < sd; eta++)
  {
    sp_cmplx_massign(eta, rho + i);
    i++;
  }
  eta = sd - vol;

}

static void spinor_cmplx_assign(int vol, spinor_dble * eta, spinor_dble * rho)
{
  spinor_dble *sd;
  int i = 0;

  sd = eta + vol;

  for (; eta < sd; eta++)
  {
    sp_cmplx_assign(eta, rho + i);
    i++;
  }
  eta = sd - vol;

}

static void spinor_sign_flip(int vol, spinor_dble * eta, spinor_dble * rho)
{
  spinor_dble *sd;
  int i = 0;
  sd = eta + vol;

  for (; eta < sd; eta++)
  {
    sign_flip(eta, rho + i);
    i++;
  }
  eta = sd - vol;
}




#define sp_weyl_assign(eta, rho)\
  rho->c1.c1.re= (eta->c1.c1.re); rho->c1.c1.im= (eta->c1.c1.im);\
  rho->c1.c2.re= (eta->c1.c2.re); rho->c1.c2.im= (eta->c1.c2.im);\
  rho->c1.c3.re= (eta->c1.c3.re); rho->c1.c3.im= (eta->c1.c3.im);\
  rho->c2.c1.re= (eta->c2.c1.re); rho->c2.c1.im= (eta->c2.c1.im);\
  rho->c2.c2.re= (eta->c2.c2.re); rho->c2.c2.im= (eta->c2.c2.im);\
  rho->c2.c3.re= (eta->c2.c3.re); rho->c2.c3.im= (eta->c2.c3.im);

#define sp_weyl_x_assign(eta, rho)\
  rho->c1.c1.re= (eta->c2.c1.re); rho->c1.c1.im= (eta->c2.c1.im);\
  rho->c1.c2.re= (eta->c2.c2.re); rho->c1.c2.im= (eta->c2.c2.im);\
  rho->c1.c3.re= (eta->c2.c3.re); rho->c1.c3.im= (eta->c2.c3.im);\
  rho->c2.c1.re= (eta->c1.c1.re); rho->c2.c1.im= (eta->c1.c1.im);\
  rho->c2.c2.re= (eta->c1.c2.re); rho->c2.c2.im= (eta->c1.c2.im);\
  rho->c2.c3.re= (eta->c1.c3.re); rho->c2.c3.im= (eta->c1.c3.im);

#define sp_weyl_massign( eta,  rho)\
  rho->c1.c1.re=-(eta->c1.c1.re); rho->c1.c1.im=-(eta->c1.c1.im);\
  rho->c1.c2.re=-(eta->c1.c2.re); rho->c1.c2.im=-(eta->c1.c2.im);\
  rho->c1.c3.re=-(eta->c1.c3.re); rho->c1.c3.im=-(eta->c1.c3.im);\
  rho->c2.c1.re=-(eta->c2.c1.re); rho->c2.c1.im=-(eta->c2.c1.im);\
  rho->c2.c2.re=-(eta->c2.c2.re); rho->c2.c2.im=-(eta->c2.c2.im);\
  rho->c2.c3.re=-(eta->c2.c3.re); rho->c2.c3.im=-(eta->c2.c3.im);

#define sp_weyl_x_massign( eta,  rho)\
  rho->c1.c1.re=-(eta->c2.c1.re); rho->c1.c1.im=-(eta->c2.c1.im);\
  rho->c1.c2.re=-(eta->c2.c2.re); rho->c1.c2.im=-(eta->c2.c2.im);\
  rho->c1.c3.re=-(eta->c2.c3.re); rho->c1.c3.im=-(eta->c2.c3.im);\
  rho->c2.c1.re=-(eta->c1.c1.re); rho->c2.c1.im=-(eta->c1.c1.im);\
  rho->c2.c2.re=-(eta->c1.c2.re); rho->c2.c2.im=-(eta->c1.c2.im);\
  rho->c2.c3.re=-(eta->c1.c3.re); rho->c2.c3.im=-(eta->c1.c3.im);

#define sp_weyl_x_muassign( eta,  rho)\
  rho->c1.c1.re=-(eta->c2.c1.re); rho->c1.c1.im=-(eta->c2.c1.im);\
  rho->c1.c2.re=-(eta->c2.c2.re); rho->c1.c2.im=-(eta->c2.c2.im);\
  rho->c1.c3.re=-(eta->c2.c3.re); rho->c1.c3.im=-(eta->c2.c3.im);\
  rho->c2.c1.re=(eta->c1.c1.re); rho->c2.c1.im=(eta->c1.c1.im);\
  rho->c2.c2.re=(eta->c1.c2.re); rho->c2.c2.im=(eta->c1.c2.im);\
  rho->c2.c3.re=(eta->c1.c3.re); rho->c2.c3.im=(eta->c1.c3.im);

#define sp_weyl_x_mlassign( eta,  rho)\
  rho->c1.c1.re=(eta->c2.c1.re); rho->c1.c1.im=(eta->c2.c1.im);\
  rho->c1.c2.re=(eta->c2.c2.re); rho->c1.c2.im=(eta->c2.c2.im);\
  rho->c1.c3.re=(eta->c2.c3.re); rho->c1.c3.im=(eta->c2.c3.im);\
  rho->c2.c1.re=-(eta->c1.c1.re); rho->c2.c1.im=-(eta->c1.c1.im);\
  rho->c2.c2.re=-(eta->c1.c2.re); rho->c2.c2.im=-(eta->c1.c2.im);\
  rho->c2.c3.re=-(eta->c1.c3.re); rho->c2.c3.im=-(eta->c1.c3.im);

#define sp_weyl_iassign( eta,  rho)\
  rho->c1.c1.re=-(eta->c1.c1.im); rho->c1.c1.im= (eta->c1.c1.re);\
  rho->c1.c2.re=-(eta->c1.c2.im); rho->c1.c2.im= (eta->c1.c2.re);\
  rho->c1.c3.re=-(eta->c1.c3.im); rho->c1.c3.im= (eta->c1.c3.re);\
  rho->c2.c1.re=-(eta->c2.c1.im); rho->c2.c1.im= (eta->c2.c1.re);\
  rho->c2.c2.re=-(eta->c2.c2.im); rho->c2.c2.im= (eta->c2.c2.re);\
  rho->c2.c3.re=-(eta->c2.c3.im); rho->c2.c3.im= (eta->c2.c3.re);

#define sp_weyl_x_iassign( eta,  rho)\
  rho->c1.c1.re=-(eta->c2.c1.im); rho->c1.c1.im= (eta->c2.c1.re);\
  rho->c1.c2.re=-(eta->c2.c2.im); rho->c1.c2.im= (eta->c2.c2.re);\
  rho->c1.c3.re=-(eta->c2.c3.im); rho->c1.c3.im= (eta->c2.c3.re);\
  rho->c2.c1.re=-(eta->c1.c1.im); rho->c2.c1.im= (eta->c1.c1.re);\
  rho->c2.c2.re=-(eta->c1.c2.im); rho->c2.c2.im= (eta->c1.c2.re);\
  rho->c2.c3.re=-(eta->c1.c3.im); rho->c2.c3.im= (eta->c1.c3.re);

#define sp_weyl_miassign( eta,  rho)\
  rho->c1.c1.re= (eta->c1.c1.im); rho->c1.c1.im=-(eta->c1.c1.re);\
  rho->c1.c2.re= (eta->c1.c2.im); rho->c1.c2.im=-(eta->c1.c2.re);\
  rho->c1.c3.re= (eta->c1.c3.im); rho->c1.c3.im=-(eta->c1.c3.re);\
  rho->c2.c1.re= (eta->c2.c1.im); rho->c2.c1.im=-(eta->c2.c1.re);\
  rho->c2.c2.re= (eta->c2.c2.im); rho->c2.c2.im=-(eta->c2.c2.re);\
  rho->c2.c3.re= (eta->c2.c3.im); rho->c2.c3.im=-(eta->c2.c3.re);

#define sp_weyl_miuassign( eta,  rho)\
  rho->c1.c1.re= (eta->c1.c1.im); rho->c1.c1.im=-(eta->c1.c1.re);\
  rho->c1.c2.re= (eta->c1.c2.im); rho->c1.c2.im=-(eta->c1.c2.re);\
  rho->c1.c3.re= (eta->c1.c3.im); rho->c1.c3.im=-(eta->c1.c3.re);\
  rho->c2.c1.re=-(eta->c2.c1.im); rho->c2.c1.im=(eta->c2.c1.re);\
  rho->c2.c2.re=-(eta->c2.c2.im); rho->c2.c2.im=(eta->c2.c2.re);\
  rho->c2.c3.re=-(eta->c2.c3.im); rho->c2.c3.im=(eta->c2.c3.re);

#define sp_weyl_milassign( eta,  rho)\
  rho->c1.c1.re=-(eta->c1.c1.im); rho->c1.c1.im= (eta->c1.c1.re);\
  rho->c1.c2.re=-(eta->c1.c2.im); rho->c1.c2.im= (eta->c1.c2.re);\
  rho->c1.c3.re=-(eta->c1.c3.im); rho->c1.c3.im= (eta->c1.c3.re);\
  rho->c2.c1.re= (eta->c2.c1.im); rho->c2.c1.im=-(eta->c2.c1.re);\
  rho->c2.c2.re= (eta->c2.c2.im); rho->c2.c2.im=-(eta->c2.c2.re);\
  rho->c2.c3.re= (eta->c2.c3.im); rho->c2.c3.im=-(eta->c2.c3.re);

#define sp_weyl_x_miassign( eta,  rho)\
  rho->c1.c1.re= (eta->c2.c1.im); rho->c1.c1.im=-(eta->c2.c1.re);\
  rho->c1.c2.re= (eta->c2.c2.im); rho->c1.c2.im=-(eta->c2.c2.re);\
  rho->c1.c3.re= (eta->c2.c3.im); rho->c1.c3.im=-(eta->c2.c3.re);\
  rho->c2.c1.re= (eta->c1.c1.im); rho->c2.c1.im=-(eta->c1.c1.re);\
  rho->c2.c2.re= (eta->c1.c2.im); rho->c2.c2.im=-(eta->c1.c2.re);\
  rho->c2.c3.re= (eta->c1.c3.im); rho->c2.c3.im=-(eta->c1.c3.re);


void adjoint(spinor_dble ** psi, spinor_dble * phi)
{
  int i, j;
  complex_dble *val, *res;
  for (j = 0; j < 12; j++)
  {
    for (i = 0; i < 12; i++)
    {
      res = (complex_dble *) (phi + i);
      val = (complex_dble *) (psi + j * VOL_SPINOR_SIZE);
      (res + j)->re = (val + i)->re;
      (res + j)->im = -(val + i)->im;
    }
  }
}

void full_adjoint(int vol, spinor_dble ** psi, spinor_dble ** phi)
{
  int i, j, k;
  complex_dble *val, *res;

  for (k = 0; k < vol; k++)
  {
    for (j = 0; j < 12; j++)
    {
      for (i = 0; i < 12; i++)
      {
        res = (complex_dble *) (phi[i] + k);
        val = (complex_dble *) (psi[j] + k);
        (res + j)->re = (val + i)->re;
        (res + j)->im = -(val + i)->im;
      }
    }
  }
}

void transpose(spinor_dble ** psi, spinor_dble * phi)
{
  int i, j;
  complex_dble *val, *res;

  for (j = 0; j < 12; j++)
  {
    for (i = 0; i < 12; i++)
    {
      res = (complex_dble *) (phi + i);
      val = (complex_dble *) (psi + j * VOL_SPINOR_SIZE);
      (res + j)->re = (val + i)->re;
      (res + j)->im = (val + i)->im;
    }
  }
}

void full_transpose(int vol, spinor_dble ** psi, spinor_dble ** phi)
{
  int i, j, k;
  complex_dble *val, *res;

  for (k = 0; k < vol; k++)
  {
    for (j = 0; j < 12; j++)
    {
      for (i = 0; i < 12; i++)
      {
        res = (complex_dble *) (phi[i] + k);
        val = (complex_dble *) (psi[j] + k);
        (res + j)->re = (val + i)->re;
        (res + j)->im = (val + i)->im;
      }
    }
  }
}


void lmult_gamma(int vol, spinor_dble ** psi, int gamma, spinor_dble ** phi)
{
  int i, k;
  spinor_dble *eta, *rho;
  switch (gamma)
  {
  case Id:
    for (i = 0; i < 12; i++)
    {
      spinor_assign(vol, psi[i], phi[i]);
    }
    break;
  case G0:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_massign(((weyl_dble *) & (rho->c3)),
                        ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_massign(((weyl_dble *) & (rho->c1)),
                        ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G1:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_miassign(((weyl_dble *) & (rho->c3)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G2:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_muassign(((weyl_dble *) & (rho->c3)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_mlassign(((weyl_dble *) & (rho->c1)),
                           ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G3:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_miuassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_milassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G5:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_assign(((weyl_dble *) & (rho->c1)),
                       ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_massign(((weyl_dble *) & (rho->c3)),
                        ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G0G5:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_assign(((weyl_dble *) & (rho->c3)),
                       ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_massign(((weyl_dble *) & (rho->c1)),
                        ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G1G5:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G2G5:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_mlassign(((weyl_dble *) & (rho->c3)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_mlassign(((weyl_dble *) & (rho->c1)),
                           ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G3G5:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_milassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_milassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G0G1:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_miassign(((weyl_dble *) & (rho->c1)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G0G2:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_muassign(((weyl_dble *) & (rho->c1)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_mlassign(((weyl_dble *) & (rho->c3)),
                           ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G0G3:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_miuassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_milassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G1G2:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_milassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_milassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G1G3:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_muassign(((weyl_dble *) & (rho->c1)),
                           ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_muassign(((weyl_dble *) & (rho->c3)),
                           ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case G2G3:
    for (i = 0; i < 12; i++)
    {
      rho = psi[i];
      eta = psi[i] + vol;
      k = 0;
      for (; rho < eta; rho++)
      {
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c1)),
                          ((weyl_dble *) & ((phi[i] + k)->c1)));
        sp_weyl_x_iassign(((weyl_dble *) & (rho->c3)),
                          ((weyl_dble *) & ((phi[i] + k)->c3)));
        k++;
      }
      psi[i] = eta - vol;
    }
    break;
  case unknown:
    error(1, 1, "[algebra.c]", "Unknown gamma structure\n");
    break;
  }
}

void rmult_gamma(int vol, spinor_dble ** psi, int gamma, spinor_dble ** phi)
{


  switch (gamma)
  {
  case Id:

    spinor_assign(vol, psi[0], phi[0]);
    spinor_assign(vol, psi[1], phi[1]);
    spinor_assign(vol, psi[2], phi[2]);

    spinor_assign(vol, psi[3], phi[3]);
    spinor_assign(vol, psi[4], phi[4]);
    spinor_assign(vol, psi[5], phi[5]);

    spinor_assign(vol, psi[6], phi[6]);
    spinor_assign(vol, psi[7], phi[7]);
    spinor_assign(vol, psi[8], phi[8]);

    spinor_assign(vol, psi[9], phi[9]);
    spinor_assign(vol, psi[10], phi[10]);
    spinor_assign(vol, psi[11], phi[11]);
    break;
  case G0:
    spinor_sign_flip(vol, psi[6], phi[0]);
    spinor_sign_flip(vol, psi[7], phi[1]);
    spinor_sign_flip(vol, psi[8], phi[2]);

    spinor_sign_flip(vol, psi[9], phi[3]);
    spinor_sign_flip(vol, psi[10], phi[4]);
    spinor_sign_flip(vol, psi[11], phi[5]);

    spinor_sign_flip(vol, psi[0], phi[6]);
    spinor_sign_flip(vol, psi[1], phi[7]);
    spinor_sign_flip(vol, psi[2], phi[8]);

    spinor_sign_flip(vol, psi[3], phi[9]);
    spinor_sign_flip(vol, psi[4], phi[10]);
    spinor_sign_flip(vol, psi[5], phi[11]);
    break;
  case G2:

    spinor_sign_flip(vol, psi[9], phi[0]);
    spinor_sign_flip(vol, psi[10], phi[1]);
    spinor_sign_flip(vol, psi[11], phi[2]);

    spinor_assign(vol, psi[6], phi[3]);
    spinor_assign(vol, psi[7], phi[4]);
    spinor_assign(vol, psi[8], phi[5]);

    spinor_assign(vol, psi[3], phi[6]);
    spinor_assign(vol, psi[4], phi[7]);
    spinor_assign(vol, psi[5], phi[8]);

    spinor_sign_flip(vol, psi[0], phi[9]);
    spinor_sign_flip(vol, psi[1], phi[10]);
    spinor_sign_flip(vol, psi[2], phi[11]);
    break;
  case G5:

    spinor_assign(vol, psi[0], phi[0]);
    spinor_assign(vol, psi[1], phi[1]);
    spinor_assign(vol, psi[2], phi[2]);

    spinor_assign(vol, psi[3], phi[3]);
    spinor_assign(vol, psi[4], phi[4]);
    spinor_assign(vol, psi[5], phi[5]);

    spinor_sign_flip(vol, psi[6], phi[6]);
    spinor_sign_flip(vol, psi[7], phi[7]);
    spinor_sign_flip(vol, psi[8], phi[8]);

    spinor_sign_flip(vol, psi[9], phi[9]);
    spinor_sign_flip(vol, psi[10], phi[10]);
    spinor_sign_flip(vol, psi[11], phi[11]);
    break;
  case G0G5:
    spinor_sign_flip(vol, psi[6], phi[0]);
    spinor_sign_flip(vol, psi[7], phi[1]);
    spinor_sign_flip(vol, psi[8], phi[2]);

    spinor_sign_flip(vol, psi[9], phi[3]);
    spinor_sign_flip(vol, psi[10], phi[4]);
    spinor_sign_flip(vol, psi[11], phi[5]);

    spinor_assign(vol, psi[0], phi[6]);
    spinor_assign(vol, psi[1], phi[7]);
    spinor_assign(vol, psi[2], phi[8]);

    spinor_assign(vol, psi[3], phi[9]);
    spinor_assign(vol, psi[4], phi[10]);
    spinor_assign(vol, psi[5], phi[11]);
    break;
  case G2G5:

    spinor_sign_flip(vol, psi[9], phi[0]);
    spinor_sign_flip(vol, psi[10], phi[1]);
    spinor_sign_flip(vol, psi[11], phi[2]);

    spinor_assign(vol, psi[6], phi[3]);
    spinor_assign(vol, psi[7], phi[4]);
    spinor_assign(vol, psi[8], phi[5]);

    spinor_sign_flip(vol, psi[3], phi[6]);
    spinor_sign_flip(vol, psi[4], phi[7]);
    spinor_sign_flip(vol, psi[5], phi[8]);

    spinor_assign(vol, psi[0], phi[9]);
    spinor_assign(vol, psi[1], phi[10]);
    spinor_assign(vol, psi[2], phi[11]);
    break;
  case G0G2:

    spinor_assign(vol, psi[3], phi[0]);
    spinor_assign(vol, psi[4], phi[1]);
    spinor_assign(vol, psi[5], phi[2]);

    spinor_sign_flip(vol, psi[0], phi[3]);
    spinor_sign_flip(vol, psi[1], phi[4]);
    spinor_sign_flip(vol, psi[2], phi[5]);

    spinor_sign_flip(vol, psi[9], phi[6]);
    spinor_sign_flip(vol, psi[10], phi[7]);
    spinor_sign_flip(vol, psi[11], phi[8]);

    spinor_assign(vol, psi[6], phi[9]);
    spinor_assign(vol, psi[7], phi[10]);
    spinor_assign(vol, psi[8], phi[11]);
    break;
  case G1G3:

    spinor_assign(vol, psi[3], phi[0]);
    spinor_assign(vol, psi[4], phi[1]);
    spinor_assign(vol, psi[5], phi[2]);

    spinor_sign_flip(vol, psi[0], phi[3]);
    spinor_sign_flip(vol, psi[1], phi[4]);
    spinor_sign_flip(vol, psi[2], phi[5]);

    spinor_assign(vol, psi[9], phi[6]);
    spinor_assign(vol, psi[10], phi[7]);
    spinor_assign(vol, psi[11], phi[8]);

    spinor_sign_flip(vol, psi[6], phi[9]);
    spinor_sign_flip(vol, psi[7], phi[10]);
    spinor_sign_flip(vol, psi[8], phi[11]);
    break;
  case G1:

    spinor_cmplx_assign(vol, psi[9], phi[0]);
    spinor_cmplx_assign(vol, psi[10], phi[1]);
    spinor_cmplx_assign(vol, psi[11], phi[2]);

    spinor_cmplx_assign(vol, psi[6], phi[3]);
    spinor_cmplx_assign(vol, psi[7], phi[4]);
    spinor_cmplx_assign(vol, psi[8], phi[5]);

    spinor_cmplx_massign(vol, psi[3], phi[6]);
    spinor_cmplx_massign(vol, psi[4], phi[7]);
    spinor_cmplx_massign(vol, psi[5], phi[8]);

    spinor_cmplx_massign(vol, psi[0], phi[9]);
    spinor_cmplx_massign(vol, psi[1], phi[10]);
    spinor_cmplx_massign(vol, psi[2], phi[11]);
    break;

  case G3:

    spinor_cmplx_assign(vol, psi[6], phi[0]);
    spinor_cmplx_assign(vol, psi[7], phi[1]);
    spinor_cmplx_assign(vol, psi[8], phi[2]);

    spinor_cmplx_massign(vol, psi[9], phi[3]);
    spinor_cmplx_massign(vol, psi[10], phi[4]);
    spinor_cmplx_massign(vol, psi[11], phi[5]);

    spinor_cmplx_massign(vol, psi[0], phi[6]);
    spinor_cmplx_massign(vol, psi[1], phi[7]);
    spinor_cmplx_massign(vol, psi[2], phi[8]);

    spinor_cmplx_assign(vol, psi[3], phi[9]);
    spinor_cmplx_assign(vol, psi[4], phi[10]);
    spinor_cmplx_assign(vol, psi[5], phi[11]);
    break;
  case G1G5:

    spinor_cmplx_assign(vol, psi[9], phi[0]);
    spinor_cmplx_assign(vol, psi[10], phi[1]);
    spinor_cmplx_assign(vol, psi[11], phi[2]);

    spinor_cmplx_assign(vol, psi[6], phi[3]);
    spinor_cmplx_assign(vol, psi[7], phi[4]);
    spinor_cmplx_assign(vol, psi[8], phi[5]);

    spinor_cmplx_assign(vol, psi[3], phi[6]);
    spinor_cmplx_assign(vol, psi[4], phi[7]);
    spinor_cmplx_assign(vol, psi[5], phi[8]);

    spinor_cmplx_assign(vol, psi[0], phi[9]);
    spinor_cmplx_assign(vol, psi[1], phi[10]);
    spinor_cmplx_assign(vol, psi[2], phi[11]);
    break;
  case G3G5:

    spinor_cmplx_assign(vol, psi[6], phi[0]);
    spinor_cmplx_assign(vol, psi[7], phi[1]);
    spinor_cmplx_assign(vol, psi[8], phi[2]);

    spinor_cmplx_massign(vol, psi[9], phi[3]);
    spinor_cmplx_massign(vol, psi[10], phi[4]);
    spinor_cmplx_massign(vol, psi[11], phi[5]);

    spinor_cmplx_assign(vol, psi[0], phi[6]);
    spinor_cmplx_assign(vol, psi[1], phi[7]);
    spinor_cmplx_assign(vol, psi[2], phi[8]);

    spinor_cmplx_massign(vol, psi[3], phi[9]);
    spinor_cmplx_massign(vol, psi[4], phi[10]);
    spinor_cmplx_massign(vol, psi[5], phi[11]);
    break;
  case G0G1:

    spinor_cmplx_massign(vol, psi[3], phi[0]);
    spinor_cmplx_massign(vol, psi[4], phi[1]);
    spinor_cmplx_massign(vol, psi[5], phi[2]);

    spinor_cmplx_massign(vol, psi[0], phi[3]);
    spinor_cmplx_massign(vol, psi[1], phi[4]);
    spinor_cmplx_massign(vol, psi[2], phi[5]);

    spinor_cmplx_assign(vol, psi[9], phi[6]);
    spinor_cmplx_assign(vol, psi[10], phi[7]);
    spinor_cmplx_assign(vol, psi[11], phi[8]);

    spinor_cmplx_assign(vol, psi[6], phi[9]);
    spinor_cmplx_assign(vol, psi[7], phi[10]);
    spinor_cmplx_assign(vol, psi[8], phi[11]);
    break;
  case G0G3:

    spinor_cmplx_massign(vol, psi[0], phi[0]);
    spinor_cmplx_massign(vol, psi[1], phi[1]);
    spinor_cmplx_massign(vol, psi[2], phi[2]);

    spinor_cmplx_assign(vol, psi[3], phi[3]);
    spinor_cmplx_assign(vol, psi[4], phi[4]);
    spinor_cmplx_assign(vol, psi[5], phi[5]);

    spinor_cmplx_assign(vol, psi[6], phi[6]);
    spinor_cmplx_assign(vol, psi[7], phi[7]);
    spinor_cmplx_assign(vol, psi[8], phi[8]);

    spinor_cmplx_massign(vol, psi[9], phi[9]);
    spinor_cmplx_massign(vol, psi[10], phi[10]);
    spinor_cmplx_massign(vol, psi[11], phi[11]);
    break;
  case G1G2:

    spinor_cmplx_assign(vol, psi[0], phi[0]);
    spinor_cmplx_assign(vol, psi[1], phi[1]);
    spinor_cmplx_assign(vol, psi[2], phi[2]);

    spinor_cmplx_massign(vol, psi[3], phi[3]);
    spinor_cmplx_massign(vol, psi[4], phi[4]);
    spinor_cmplx_massign(vol, psi[5], phi[5]);

    spinor_cmplx_assign(vol, psi[6], phi[6]);
    spinor_cmplx_assign(vol, psi[7], phi[7]);
    spinor_cmplx_assign(vol, psi[8], phi[8]);

    spinor_cmplx_massign(vol, psi[9], phi[9]);
    spinor_cmplx_massign(vol, psi[10], phi[10]);
    spinor_cmplx_massign(vol, psi[11], phi[11]);
    break;
  case G2G3:

    spinor_cmplx_assign(vol, psi[3], phi[0]);
    spinor_cmplx_assign(vol, psi[4], phi[1]);
    spinor_cmplx_assign(vol, psi[5], phi[2]);

    spinor_cmplx_assign(vol, psi[0], phi[3]);
    spinor_cmplx_assign(vol, psi[1], phi[4]);
    spinor_cmplx_assign(vol, psi[2], phi[5]);

    spinor_cmplx_assign(vol, psi[9], phi[6]);
    spinor_cmplx_assign(vol, psi[10], phi[7]);
    spinor_cmplx_assign(vol, psi[11], phi[8]);

    spinor_cmplx_assign(vol, psi[6], phi[9]);
    spinor_cmplx_assign(vol, psi[7], phi[10]);
    spinor_cmplx_assign(vol, psi[8], phi[11]);
    break;
  case unknown:
    error(1, 1, "[algebra.c]", "Unknown gamma structure\n");
    break;
  }

}
