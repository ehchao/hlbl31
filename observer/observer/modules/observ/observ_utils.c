/**
* @file observ_utils.c
* \brief  Utility functions copied from ms4.c from Martin Luescher.
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-06-13
*/
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mpi.h"
#include "observ.h"

#define OBSERV_UTILS_C

#define MAX(n,m) \
   if ((n)<(m)) \
      (n)=(m)

static void dfl_wsize(int *nws,int *nwv,int *nwvd)
{
   dfl_parms_t dp;
   dfl_pro_parms_t dpp;

   dp=dfl_parms();
   dpp=dfl_pro_parms();

   MAX(*nws,dp.Ns+2);
   MAX(*nwv,2*dpp.nkv+2);
   MAX(*nwvd,4);
}

static void wsize(int *nws,int *nwsd,int *nwv,int *nwvd)
{
  int nsd, k;
   (*nws)=0;
   (*nwsd)=0;
   (*nwv)=0;
   (*nwvd)=0;

   nsd=2;

   for (k = 0; k < observer.no_prop; k++) {
     MAX(*nws, 2*solver_parms(prop[k].idx_solver).nkv + 2);
   }
      MAX(*nwsd,1+4);
      dfl_wsize(nws,nwv,nwvd);
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Allocates the workspace necessary given the solver parameters in 
* the infile in addition one cans specifiy the amount of additional workspaces
* to be allocated.
*
* \param nws_loc  allocate nws_loc additional single precision spinor fields.
* \param nwsd_loc allocate nwsd_loc additional double precision spinor fields.
* \param nwv_loc  allocate nvw_loc additional single precision vector fields.
* \param nwvd_loc allocate nwvd_loc additional double precision vector fields.
*/
/* ----------------------------------------------------------------------------*/
void alloc_wspace(int nws_loc,int nwsd_loc, int nwv_loc, int nwvd_loc) 
{
  int nws,nwsd,nwv,nwvd;
  wsize(&nws,&nwsd,&nwv,&nwvd);
  alloc_ws(nws   + nws_loc);
  alloc_wsd(nwsd + nwsd_loc); /*+1 because of the source fields below*/
  alloc_wv(nwv   + nwv_loc);
  alloc_wvd(nwvd + nwvd_loc);
}


static struct gammatypes gamma_array[] = {
{G0,"G0"},
{G1,"G1"},
{G2,"G2"},
{G3,"G3"},
{Id,"Id"},
{G5,"G5"},
{G0G5,"G0G5"},
{G1G5,"G1G5"},
{G2G5,"G2G5"},
{G3G5,"G3G5"},
{G0G1,"G0G1"},
{G0G2,"G0G2"},
{G0G3,"G0G3"},
{G1G2,"G1G2"},
{G1G3,"G1G3"},
{G2G3,"G2G3"}
};

/* --------------------------------------------------------------------------*/
/**
* \brief  Converts string to enum type of the gamma structure
*
* \param str Gamma matrix as string is converted to enum. Currently the
* following Gamma encoding is used (this is the saem as in the old measure code)
*{G0,"G0"},\n
*{G1,"G1"},\n
*{G2,"G2"},\n
*{G3,"G3"},\n
*{Id,"Id"},\n
*{G5,"G5"},\n
*{G0G5,"G0G5"},\n
*{G1G5,"G1G5"},\n
*{G2G5,"G2G5"},\n
*{G3G5,"G3G5"},\n
*{G0G1,"G0G1"},\n
*{G0G2,"G0G2"},\n
*{G0G3,"G0G3"},\n
*{G1G2,"G1G2"},\n
*{G1G3,"G1G3"},\n
*{G2G3,"G2G3"}\n
*
* \return   Integer corresponding to the gamma structure.
*/
/* ----------------------------------------------------------------------------*/



void write_all_gammas(char *all_gammas)
{
  int i;
  const int size = sizeof(gamma_array) / sizeof(gamma_array[0]);
  sprintf(all_gammas,"\n");
  for(i=0;i< size; i++)
  {
    sprintf(all_gammas,"%s\n%s",all_gammas,gamma_array[i].str);
  }

}

int cnvrt_gamma(char* str)
{
 const int size = sizeof(gamma_array) / sizeof(gamma_array[0]);
 int i;
 int unknown=-1;
 for(i = 0; i < size; i++)
 {
   if(strcmp(gamma_array[i].str, str) == 0)
   return gamma_array[i].gamma;
 }
 return unknown;
}

const char *cnvrt_gamma_str(int i)
{
return (const char*) gamma_array[i].str;
}
static int mapa=0,iupa=0,idna=0,ipta=0;




static struct sourcetypes src_array[] = {
{POINT,"POINT"},
{WUPPERTAL,"WUPPERTAL"},
{JACOBI,"JACOBI"},
};


/* --------------------------------------------------------------------------*/
/**
* \brief  Converts string to enum of source type.
*
* \param str: Name of the source type
*
* \return   
*/
/* ----------------------------------------------------------------------------*/
int cnvrt_src(char* str)
{
 const int size = sizeof(src_array) / sizeof(src_array[0]);
 int i;
 int unknown=-1;
 for(i = 0; i < size; i++)
 {
   if(strcmp(src_array[i].str, str) == 0)
   return src_array[i].src;
 }
 return unknown;
}

const char *cnvrt_src_str(int i)
{
return (const char*) src_array[i].str;
}

static struct smearingtypes smear_array[] = {
{APE,"APE"},
};


/* --------------------------------------------------------------------------*/
/**
* \brief  Converts string to enum of source type.
*
* \param str: Name of the source type
*
* \return   
*/
/* ----------------------------------------------------------------------------*/
int cnvrt_smear(char* str)
{
 const int size = sizeof(smear_array) / sizeof(smear_array[0]);
 int i;
 int unknown=-1;
 for(i = 0; i < size; i++)
 {
   if(strcmp(smear_array[i].str, str) == 0)
   return smear_array[i].src;
 }
 return unknown;
}

const char *cnvrt_smear_str(int i)
{
return (const char*) smear_array[i].str;
}

void init_arrays(int argc, char* argv[])
{
  /*int i,j;*/
  int my_rank;
  int ifile;
    
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {
    ifile = find_opt(argc, argv, "-ogeom");
    if(ifile!=0)
    {
      sscanf(argv[1+ifile],"%d",&NPROC0); 
      sscanf(argv[2+ifile],"%d",&NPROC1); 
      sscanf(argv[3+ifile],"%d",&NPROC2); 
      sscanf(argv[4+ifile],"%d",&NPROC3); 
    }
      ifile = find_opt(argc, argv, "-lgeom");
    if(ifile!=0)
    {
      sscanf(argv[1+ifile],"%d",&L0); 
      sscanf(argv[2+ifile],"%d",&L1); 
      sscanf(argv[3+ifile],"%d",&L2); 
      sscanf(argv[4+ifile],"%d",&L3); 
    }
}
  MPI_Bcast(&L0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&L1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&L2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&L3, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&NPROC0, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NPROC1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NPROC2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&NPROC3, 1, MPI_INT, 0, MPI_COMM_WORLD);

  if(mapa==0)
  {
     map=malloc((BNDRY+NPROC%2)*sizeof(int));
     mapa=1;
  }
  if(iupa==0)
  {
     iup=malloc(4*VOLUME*sizeof(int));
     iupa=1;
  }
  if(idna==0)
  {
     idn=malloc(4*VOLUME*sizeof(int));
     idna=1;
  }
  
  if(ipta==0)
  {
    ipt=malloc(VOLUME*sizeof(int));
    ipta=1;
  }
  
}



static char line[NAME_SIZE+1];

static void check_tag(char *tag)
{
   if (tag[0]=='\0')
      return;

   error_root((strspn(tag," 0123456789.")!=0L)||
              (strcspn(tag," \n")!=strlen(tag)),1,
              "check_tag [mutils.c]","Improper tag %s",tag);
}
static char *get_line(void)
{
   char *s,*c;

   s=fgets(line,NAME_SIZE+1,stdin);

   if (s!=NULL)
   {
      error_root(strlen(line)==NAME_SIZE,1,"get_line [mutils.c]",
                 "Input line is longer than NAME_SIZE-1");   

      c=strchr(line,'#');
      if (c!=NULL)
         c[0]='\0';
   }
   
   return s;
}

static long find_tag_optional(char *tag)
{
   int ie;
   long tofs,lofs,ofs;
   char *s,*pl,*pr;

   ie=0;
   tofs=-1L;   
   lofs=ftell(stdin);
   rewind(stdin);
   ofs=ftell(stdin);
   s=get_line();
   
   while (s!=NULL)
   {
      pl=strchr(line,'[');
      pr=strchr(line,']');
         
      if ((pl==(line+strspn(line," \t")))&&(pr>pl))
      {
         if (ofs<lofs)
         {
            ie=0;
            tofs=-1L;
         }
         else
            break;
      }
      else
      {
         pl=line+strspn(line," \t");
         pr=pl+strcspn(pl," \t\n");
         pr[0]='\0';
         
         if (strcmp(pl,tag)==0)
         {
            if (tofs!=-1L)
               ie=1;
            tofs=ofs;
         }
      }

      ofs=ftell(stdin);
      s=get_line();
   }

//   error_root(tofs==-1L,1,"find_tag [mutils.c]","Tag %s not found",tag);   
   error_root(ie!=0,1,"find_tag [mutils.c]",
              "Tag %s occurs more than once in the current section",tag);   
   if(tofs!=-1L)
   {
   ie=fseek(stdin,tofs,SEEK_SET);
   error_root(ie!=0,1,"find_tag [mutils.c]",
              "Unable to go to line with tag %s",tag);
   }
   else
   {
   ie=fseek(stdin,lofs,SEEK_SET);
   error_root(ie!=0,1,"find_tag [mutils.c]",
              "Unable to go to line with tag %s",tag);
   }
   return tofs;
}

long read_line_optional(char *tag,char *format,...)
{
   int my_rank,is,ic;
   long tofs;
   char *pl,*p;
   va_list args;

   MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

   if (my_rank==0)
   {
      check_tag(tag);
      
      if (tag[0]!='\0')
      {
         tofs=find_tag_optional(tag);
         if(tofs<0) return 0;
         get_line();   
         pl=line+strspn(line," \t");
         pl+=strcspn(pl," \t\n");
      }
      else
      {
         p=format;
         p+=strspn(p," ");
         error_root(strstr(p,"%s")==p,1,"read_line [mutils.c]",
                    "String data after empty tag");
         tofs=ftell(stdin);
         pl=get_line();
      }
      
      va_start(args,format);      

      for (p=format;;)
      {
         p+=strspn(p," ");
         ic=0;
         is=2;

         if ((p[0]=='\0')||(p[0]=='\n'))
            break;
         else if (p==strstr(p,"%s"))
            ic=sscanf(pl,"%s",va_arg(args,char*));
         else if (p==strstr(p,"%d"))
            ic=sscanf(pl,"%d",va_arg(args,int*));
         else if (p==strstr(p,"%f"))
            ic=sscanf(pl,"%f",va_arg(args,float*));
         else if (p==strstr(p,"%lf"))
         {
            is=3;
            ic=sscanf(pl,"%lf",va_arg(args,double*));
         }
         else
            error_root(1,1,"read_line [mutils.c]",
                       "Incorrect format string %s on line with tag %s",
                       format,tag);
         
         error_root(ic!=1,1,"read_line [mutils.c]",
                    "Missing data item(s) on line with tag %s",tag);

         p+=is;
         pl+=strspn(pl," \t");
         pl+=strcspn(pl," \t\n");
      }

      va_end(args);

      return tofs;
   }
   else
      return -1L;
}

int compare_structs_simple(prop_parms_t * t1, prop_parms_t * t2)
{
  int comp;

  /* Compare up to solve */
  comp = 0;
  comp = !(t1->dump == t2->dump);
  comp += !(strcmp(t1->src_t, t2->src_t) == 0);
  comp += !(t1->kappa == t2->kappa);
  comp += !(t1->mu == t2->mu);
  comp += !(t1->link_smearing_idx == t2->link_smearing_idx);
  comp += !(t1->smearing_steps_source == t2->smearing_steps_source);
  comp += !(t1->smearing_parm_source == t2->smearing_parm_source);
  comp += !(t1->idx_solver == t2->idx_solver);

  return comp;
}
