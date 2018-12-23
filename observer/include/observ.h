/**
 * @mainpage
 * This is the observer package.
 * */


/**
* @file observ.h
* \brief  Global defintions for the observer package.
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-06-04
*/

#ifndef OBSERV_H
#define OBSERV_H
/*
#ifndef PI 
#define PI 3.141592653589793
#endif
*/
#include "global.h"
#include "stdio.h"
#include "flags.h"
#include "observer_version.h"

#ifndef MAIN_PROGRAM
#define EXTERN extern
#else
#define EXTERN
#endif

#define   SPATIAL_VOLUME_LOOP_START(x1,x2,x3)\
  for(x1=0;x1<L1;x1++){\
  for(x2=0;x2<L2;x2++){\
  for(x3=0;x3<L3;x3++){

#define  SPATIAL_VOLUME_LOOP_END\
  }}}

#define   VOLUME_LOOP_START(x0,x1,x2,x3)\
  for(x0=0;x0<L0;x1++){\
  for(x1=0;x1<L1;x1++){\
  for(x2=0;x2<L2;x2++){\
  for(x3=0;x3<L3;x3++){

#define  VOLUME_LOOP_END\
  }}}}

typedef enum
{
  G0,
  G1,
  G2,
  G3,
  Id,
  G5,
  G0G5,
  G1G5,
  G2G5,
  G3G5,
  G0G1,
  G0G2,
  G0G3,
  G1G2,
  G1G3,
  G2G3,
  unknown = 999
} gammas;

struct gammatypes
{
  gammas gamma;
  const char *str;
};


typedef enum
{
  POINT, WUPPERTAL,JACOBI
} src_t;

typedef enum
{
  APE
} smear_t;

struct smearingtypes
{
  smear_t src;
  const char *str;
};

struct sourcetypes
{
  src_t src;
  const char *str;
};

typedef enum
{
  INIT, ALLOCATED, COMPUTED, RELEASED, DUMPED
} prop_list_status_t;

typedef struct
{
  int spin_idx_status[12];
  int cnfg;
  spinor_dble *sd[12];
} prop_list_t;

typedef struct
{
  int src_pos[4];
}
src_position;

typedef struct
{
  int no_prop;
  int no_link_smearings;
  int no_src_pos;
  int first_cnfg;
  int last_cnfg;
  int step_cnfg;
  char *infile_contents;
  int level;
  int seed;
  int **src_pos;
  int *idx_prop;
  prop_list_t *prop_list;
} observer_parms_t;


typedef struct
{
  char src_t[NAME_SIZE];
  int src_pos[4];
  double kappa;
  double mu;
  int link_smearing_idx;
  int smearing_steps_source;
  double smearing_parm_source;
  int idx_solver;
  int smearing_steps_sink;
  double smearing_parm_sink;
  int dump;
  int truncated_solver;
} prop_parms_t;


typedef struct
{
  int status;
  void *u_field;
}
link_smearing_stat_t;
 
typedef struct
{
  char smearing_t[NAME_SIZE];
  int smearing_steps;
  int skipdir;
  double smearing_parm;
  link_smearing_stat_t stat;
} smearing_parms_t;



typedef struct
{
  su3_dble c11, c12, c13, c14, c21, c22, c23, c24, c31, c32, c33, c34,
    c41, c42, c43, c44;
} full_spinor_dble;

typedef struct
{
  int no_template;
  int no_meson_2pt;
} app_parms_t;

EXTERN char log_dir[NAME_SIZE], loc_dir[NAME_SIZE];
EXTERN char dat_dir[NAME_SIZE];
EXTERN char cnfg_dir[NAME_SIZE], sfld_dir[NAME_SIZE];
EXTERN char log_file[NAME_SIZE], log_save[NAME_SIZE], end_file[NAME_SIZE];
EXTERN char cnfg_file[NAME_SIZE], sfld_file[NAME_SIZE], nbase[NAME_SIZE];
EXTERN lat_parms_t lat;
EXTERN bc_parms_t bcp;

EXTERN prop_parms_t *prop;
EXTERN smearing_parms_t *smearing_parms;
EXTERN solver_parms_t sp;
EXTERN observer_parms_t observer;
EXTERN app_parms_t apps;
EXTERN double mus;
EXTERN int level;
EXTERN int noexp, endian;

EXTERN int compare_structs_simple(prop_parms_t *, prop_parms_t *);
EXTERN void create_source(char *name, int pos[4], int dirac_index,
                          int color_index, spinor_dble * src_field);
EXTERN void point_source(int pos[4], int dirac_index, int color_index,
                         spinor_dble * src_field);
EXTERN void mpi_print_stats(double t0, double t1, char *message);
EXTERN void read_infile(int argc, char *argv[]);
EXTERN void solve_dirac(spinor_dble * eta, spinor_dble * psi, int *status,
                        int idx_prop);
EXTERN void solve_dirac_comp(spinor_dble * eta, spinor_dble * psi, int *status,
                        int idx_prop, int comp);
EXTERN spinor_dble *propagator(int idx_prop, int isrc);
EXTERN spinor_dble **full_propagator(int idx_prop);
EXTERN void print_info(void);
EXTERN void alloc_prop(int, int);
EXTERN void release_prop(int, int);
EXTERN void release_full_prop(int);
EXTERN void cnvrt_spinor_to_full_spinor(int volume, spinor_dble ** eta,
                                        full_spinor_dble * fpsi);
EXTERN void alloc_wspace(int, int, int, int);
EXTERN int cnvrt_gamma(char *);
EXTERN int cnvrt_src(char *);
EXTERN int cnvrt_smear(char *);
EXTERN const char *cnvrt_smear_str(int);
EXTERN const char *cnvrt_gamma_str(int);
EXTERN const char *cnvrt_src_str(int);
EXTERN void rmult_gamma(int vol, spinor_dble **, int, spinor_dble **);
EXTERN void lmult_gamma(int vol, spinor_dble **, int, spinor_dble **);
EXTERN void full_adjoint(int vol, spinor_dble **, spinor_dble **);
EXTERN void adjoint(spinor_dble **, spinor_dble *);
EXTERN void transpose(spinor_dble **, spinor_dble *);
EXTERN void full_transpose(int vol, spinor_dble **, spinor_dble **);
EXTERN void recompute_propagators(void);
EXTERN void write_all_gammas(char *);
EXTERN void set_source_position(int src_pos[4]);
EXTERN void set_isource_position(int src_pos[4], int index);
EXTERN long read_line_optional(char *tag, char *format, ...);

/* Legacy algebra modules */
EXTERN void mul_gamma_r(int mu, full_spinor_dble * in, full_spinor_dble * out);
EXTERN void mul_gamma_l(int mu, full_spinor_dble * in, full_spinor_dble * out);
EXTERN void copy_fs_sv(int vol, full_spinor_dble * sv, spinor_dble * s, int id,
                       int ic);
EXTERN void adj_full_spinor(full_spinor_dble * in, full_spinor_dble * out);
EXTERN void meson_trace(full_spinor_dble * in1, full_spinor_dble * in2,
                        complex_dble * out);
EXTERN void init_arrays(int argc, char *argv[]);


#endif
