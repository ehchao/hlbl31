/**
* @file alloc_prop.c
* \brief  Allocate and keep track of allocated propagators.
* Uses the routines from openQCD modules/utils/wspace.c
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-06-05
*
*
*/

#define ALLOC_PROP_C
#include <stdlib.h>
#include "observ.h"

static spinor_dble **wsd0,sd0={{{0.0}}};


void alloc_prop(int idx_prop,int spin_idx)
{
   spinor_dble *sd,*sm;

   if(observer.prop_list[idx_prop].sd[spin_idx]==NULL)
   {

      wsd0=malloc(sizeof(*wsd0));
 
      error_loc((wsd0==NULL),1,"alloc_wsd [wspace.c]",
            "Unable to allocate index arrays");

      wsd0[0]=amalloc(NSPIN*sizeof(**wsd0),ALIGN);

      error_loc(wsd0[0]==NULL,1,"alloc_prop [alloc_prop.c]",
            "Unable to allocate workspace");
     
      observer.prop_list[idx_prop].spin_idx_status[spin_idx]=ALLOCATED;
      observer.prop_list[idx_prop].sd[spin_idx]=wsd0[0];
      sd=wsd0[0];
      sm=sd+NSPIN;

      for (;sd<sm;sd++)
        (*sd)=sd0;
   }
     error_chk();
}


void release_prop(int idx_prop,int spin_idx)
{
    
   if(observer.prop_list[idx_prop].sd[spin_idx]!=NULL)
   {
     afree(observer.prop_list[idx_prop].sd[spin_idx]);
     observer.prop_list[idx_prop].sd[spin_idx]=NULL;
     observer.prop_list[idx_prop].spin_idx_status[spin_idx]=RELEASED;
   }
}

void release_full_prop(int idx_prop)
{
  int i;

  for(i=0;i<12;i++)
  {
    release_prop(idx_prop,i);
  }    
}


