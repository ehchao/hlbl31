/**
* @file sources.c
* \brief  Implementation of various sources.
*
* Currently implemented sources:
* -# point source
* @author Dalibor Djukanovic
* @version 0.1
* @date 2013-12-31
*
* v0.1: Initial Version.
*/


#define SOURCES_C

#include <string.h>
#include "mpi.h"
#include "su3.h"
#include "observ.h"
#include "lattice.h"


/* --------------------------------------------------------------------------*/
/**
* \brief  Generates point source
*
*		  Creates a point source at source position pos[]. The point source 
*		  is defined as
*		  \f[
*		  \eta_{x,\alpha,j}^{(y,\beta,k)} = \delta_{xy}\delta_{\beta\gamma}
*		  	\delta_{jk}
*		  \f]
*
* \param pos[]	Source position given as a four-vector \f$x_\mu\f$ indicating position as:\n
* 		pos[0] = t \n
* 		pos[1] = x \n
* 		pos[2] = y \n
* 		pos[3] = z \n
* \param dirac_index \f$\alpha\f$ of source field spinor_dble (C-style indexing, i.e. first component starts with 0)
* \param color_index j is color index of source field spinor_dble (C-style indexing, i.e. first component starts with 0)
*
* \return Returns global spinor_dble field of point source.
*/
/* ----------------------------------------------------------------------------*/

void point_source(int pos[4], int dirac_index, int color_index, spinor_dble *source_fld)
{
  int process_rank;
  int local_index;
  int my_rank;
  double *spinor_comp;
  spinor_dble sd0={{{0.f}}};
  spinor_dble *sd,*sm;
 
        
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  ipt_global(pos,&process_rank,&local_index);
/*  
    message("Putting point source @ index %d  on %d @ spinor index %d @ color index %d\n",\
    local_index,process_rank,dirac_index,color_index);
*/
  sd=source_fld;
  sm=sd+VOLUME;
  for(;sd < sm;sd++)
  { 
    (*sd)=sd0;
    if(my_rank == process_rank)
    {
      if(sd == (source_fld+local_index))
      {
        spinor_comp=(double*)(sd);
        spinor_comp[6*dirac_index+2*color_index]=1.0;
      }
    }
  }  
}

/* --------------------------------------------------------------------------*/
/**
* \brief Create sources with a given source type and source position
*
* \param src Type of source possible values:\n
* 	POINT
* \param pos[4]  Source position given as a four-vector \f$x_\mu\f$ indicating position as:\n
*        pos[0] = t \n
*        pos[1] = x \n
*        pos[2] = y \n
*        pos[3] = z \n
* \param dirac_index \f$\alpha\f$ of source field spinor_dble (C-style indexing, i.e. first component starts with 0)
* \param color_index j is color index of source field spinor_dble (C-style indexing, i.e. first component starts with 0)
* \param src_field Spinor to which the source is saved.
*/
/* ----------------------------------------------------------------------------*/
void create_source(char *name,int pos[4], int dirac_index, int color_index,spinor_dble * src_field)
{
  if (strcmp(name,"POINT") == 0)
  {
    point_source(pos,dirac_index,color_index,src_field);
  }
  else
  {
    error(1,1,"create_source","Source %s unkown",name);
  }
}
