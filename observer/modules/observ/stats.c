/**
 * @file stats.c
 * \brief  Methods to gather MPI-Statistics.
 * @author Dalibor Djukanovic
 * @version 0.1
 * @date 2014-02-17
 *
 * v0.1: Initial Version.
 */

#define STATS_C

#include <mpi.h>
#include "utils.h"

typedef struct
{
  double value;
  int rank;
} buf_stat;


/* --------------------------------------------------------------------------*/
/**
 * \brief  Prints mpi timing statistics. The ouput contains
 * maximum, minimum and average wall times.
 *
 * \param t0 Start time, i.e. t0<t1
 * \param t1 End time
 */
/* ----------------------------------------------------------------------------*/
void mpi_print_stats(double t0, double t1, char *message_name)
{
  buf_stat max;
  buf_stat min;
  buf_stat val;
  int myrank;
  double mean;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  val.value = t1 - t0;
  val.rank = myrank;
  MPI_Reduce(&val, &max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(&val, &min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
  MPI_Reduce(&(val.value), &mean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /*message("Timing Statistics - %s : Average %e [sec]\n\
     Max %e [sec] on rank %d -> Deviation is %2d %\n\
     Min %e [sec] on rank %d -> Deviation is %2d %\n",\
     message_name,mean/size,\
     max.value,max.rank,(int)((1-max.value/(mean/size))*100),\
     min.value,min.rank,(int)((1-min.value/(mean/size))*100)); */
  message
    ("Timing Statistics - %s : (Average,Min,Max,MinDev,MaxDev,Rankmin,Rankmax): (%e,%e,%e,%2d,%2d,%d,%d) \n",
     message_name, mean / size, min.value, max.value,
     (int) ((1 - max.value / (mean / size)) * 100),
     (int) ((1 - min.value / (mean / size)) * 100), min.rank, max.rank);
}
