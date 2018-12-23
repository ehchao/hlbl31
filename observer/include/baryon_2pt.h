#include <qdp.h>

void read_baryon_2pt_infile(int argc, char *argv[]);
void print_baryon_2pt_info(void);

namespace QDP
{

  void calculate_baryons(const LatticePropagatorD & light_quark_prop,
                         const LatticePropagatorD & strange_quark_prop,
                         const SftMom & fmom,
                         hdf5_writer & writer,
                         const char *group_name, int isrc, int prop_comb);

  void contract_octet_baryons(const LatticePropagatorD & light_quark_prop,
                              const LatticePropagatorD & strange_quark_prop,
                              const SpinMatrix Pol, const SftMom & fmom,
                              hdf5_writer & writer, const char *group_name,
                              int src_pos_id, int prop_comb);
  void write_baryon_2pt_data_init(hdf5_writer & writer, int ndsets, ...);


}
