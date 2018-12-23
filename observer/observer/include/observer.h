#ifndef OBSERVER_H
#define OBSERVER_H

#include "qmp.h"
#include <qdp.h>
#include "qdp_util.h"
#include <sstream>
#include <string>
#include <iomanip>
#include "qdpinterface.h"
#include "mpi.h"
#include "smearing.h"
#include "hdf5_output.h"
#include "sftmom.h"
extern "C"
{
#include "su3.h"
#include "random.h"
#include "su3fcts.h"
#include "flags.h"
#include "utils.h"
#include "lattice.h"
#include "uflds.h"
#include "sflds.h"
#include "linalg.h"
#include "sw_term.h"
#include "dirac.h"
#include "global.h"
#include "block.h"
#include "archive.h"
#include "observ.h"
#include "apps.h"
#include "dfl.h"
}

/* Define timing routines */




#ifdef MAIN_PROGRAM
int global_indent = 0;
std::stringstream time_line_out;
#else
extern int global_indent;
extern std::stringstream time_line_out;
#endif


#define DETAILLEVEL 1

#define start_time(NAME) global_indent++;QDP::StopWatch stopwatch_timer##NAME;stopwatch_timer##NAME.reset();stopwatch_timer##NAME.start();

/* Simple output with indentation */
/*#define stop_time(NAME) stopwatch_timer##NAME.stop(); { std::stringstream tmp_out_buf;  tmp_out_buf<<time_line_out.rdbuf(); time_line_out.str("");time_line_out.clear(); time_line_out<<std::setw(2*global_indent)<<std::setfill('*')<<"*" <<" Time for " #NAME " \
is " << stopwatch_timer##NAME.getTimeInSeconds()<<" [sec]"<<endl<<tmp_out_buf.str();};global_indent--;if(global_indent<=DETAILLEVEL){QDPIO::cout<<time_line_out.str()<<endl;time_line_out.str("");time_line_out.clear();}*/

#define stop_time(NAME) stopwatch_timer##NAME.stop(); { std::stringstream tmp_out_buf;  tmp_out_buf<<std::setw(2*global_indent)<<std::setfill('*')<<"*" <<" Time for " #NAME " \
is " << stopwatch_timer##NAME.getTimeInSeconds()<<" [sec]"<<endl;QDPIO::cout<<tmp_out_buf.str();};global_indent--;



namespace QDP
{

typedef struct
{
int idx;
multi1d<LatticeColorMatrixD> u_field;
}
link_smearing_stat_t;

#ifdef MAIN_PROGRAM
link_smearing_stat_t lsmearing_stat={-1};
#else
extern link_smearing_stat_t lsmearing_stat;
#endif 
void import_config(int cnfg_no);
}

#endif
