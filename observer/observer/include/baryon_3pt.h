#ifndef BARYON_3PT_H
#define BARYON_3PT_H
#include <qdp.h>
#include "sftmom.h"
#include "hdf5_output.h"

#ifndef MAIN_PROGRAM
#define EXTERN extern
#else
#define EXTERN
#endif

extern void read_baryon_3pt_infile(int argc, char* argv[]);
extern void print_baryon_3pt_info(void);
extern void test(void);
/* baryon 3pt funnction*/
namespace QDP
{

EXTERN void calculate_baryon_3pt_full(const LatticePropagatorD & prop1,
                                 const LatticePropagatorD & prop2,
                                 const multi1d<LatticeColorMatrix> &,
                                 const SftMom & fmom,
                                 hdf5_writer & write,
                                 const char *group_name, int src_pos_ind,
                                 int prop_combination,
                                 int ts);

EXTERN void baryon_3pt(const LatticePropagatorD &  quark_prop_0,
                LatticePropagatorD & SeqPropU,
                LatticePropagatorD & SeqPropD, SpinMatrix & OA,
                SpinMatrix & GAB, SpinMatrix & GammaPol, int propid, int ts);

EXTERN void write_baryon_3pt_data_init(hdf5_writer & writer, int ndsets, ...);
}

#endif
