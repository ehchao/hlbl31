#ifndef MESON_2PT_H
#define MESON_2PT_H

#include "qdp.h"

#ifndef MAIN_PROGRAM
#define EXTERN extern
#else
#define EXTERN
#endif

void read_meson_2pt_infile(int argc, char *argv[]);
void print_meson_2pt_info(void);
namespace QDP
{

  void write_meson_2pt_data_init(hdf5_writer & writer, int ndset, ...);
}
#endif
