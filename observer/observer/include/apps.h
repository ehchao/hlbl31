#ifndef APPS_H
#define APPS_H

#ifndef MAIN_PROGRAM
#define EXTERN extern
#else
#define EXTERN
#endif
#include "observ.h"
#define GAMMA_CHAR_SIZE 16

typedef struct
{
  int (*idx_prop)[2];           /*Need 2 quark propagators */
  char (*gamma_insertion_source)[GAMMA_CHAR_SIZE];
  char (*gamma_insertion_sink)[GAMMA_CHAR_SIZE];
  int max_momentum;
  int no_of_gso;
  int no_of_gsi;
  int no_meson_2pt;
} meson_2pt_meas_parm_t;

EXTERN meson_2pt_meas_parm_t meson_2pt_meas;

typedef struct
{
  int (*idx_prop)[2];           /*Need 2 quark propagators */
  int max_momentum;
  int no_baryon_2pt;
  char polarization[NAME_SIZE];
} baryon_2pt_meas_parm_t;


EXTERN baryon_2pt_meas_parm_t baryon_2pt_meas;

typedef struct
{
  int (*idx_prop)[2];           /*Need 2 quark propagators */
  // char (*gamma_insertion)[GAMMA_CHAR_SIZE];
  int max_momentum;
  int no_of_ts;
  int *ts;
  int no_baryon_3pt;
} baryon_3pt_meas_parm_t;

EXTERN baryon_3pt_meas_parm_t baryon_3pt_meas;


typedef struct
{
  int random_source_shift;
  int idx_prop;
  int forward_scattering;
  int momentum_0[4];
  int momentum_1[4];
  int max_momentum[4];
} light_by_light_parm_t;

extern light_by_light_parm_t light_by_light_meas;

typedef struct
{
  int lattice_region;
  int random_source_shift;
  int idx_prop;
  char typeofVwVyVxVz[5];
  int num_kernel;
  char (*kernel)[NAME_SIZE];
  int num_vector_y;
  int (*vector_y)[4];
  double atimesmuonmass;
} hlbl_g_minus_2_param_t;

EXTERN hlbl_g_minus_2_param_t hlbl_g_minus_2_meas;

/* Necessary modules to be accessed from other programs */
/*void print_meson_2pt_info(void);
void read_meson_2pt_infile(int argc, char *argv[]);
void calculate_meson_2pt(meson_2pt_meas_parm_t meson,int, int src_pos[4]); 
void write_meson_2pt_data(int config_no,int src_pos[4], meson_2pt_meas_parm_t meson, int *momenta,
                     complex_dble * result);

*/

typedef struct{
    int nprop;
    int* id_prop;
    int vector_y[4];

    //check_mode = 1 : only check; check_mode = 0: only calculation 
    int check_mode;
   //random_offset = 1 : add random offset to the positions 
    int random_offset;
    
    int kernel_type;
    // a times the fermion mass
    double am;
    double a;
} contraction_3_1_param_t;
    
EXTERN contraction_3_1_param_t contraction_3_1_meas;
    
/* Legacy functions */

void calculate_vacuum_polarisation(void);


#endif
