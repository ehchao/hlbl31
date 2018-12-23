#include "mpi.h"
extern "C"
{
#include "global.h"
#include "observ.h"
#include "apps.h"
//parameters param;
}
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
//#include "3_1_contraction.h"
  
void read_hlbl31(){

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //int nprop;
    //broadcast if root
    if(rank == 0){
        find_section("hlbl g minus 2");
        //suppose here that they are called like "idx_prop n " and there are exactly Nf prop
   
        /*for(int i =0; i< Nf; i++ ){
            char *propid;
            sprintf(propid, "idx_prop%d", i );
            read_iprms(propid, 1 , param.prop_id);
        }*/
        read_iprms("nprop", 1, & contraction_3_1_meas.nprop);
        contraction_3_1_meas.id_prop = (int*) malloc(contraction_3_1_meas.nprop*sizeof(int));
        read_iprms("idx_prop",contraction_3_1_meas.nprop, contraction_3_1_meas.id_prop);
        read_iprms("vector_y" ,4 ,contraction_3_1_meas.vector_y);
        read_iprms("check_mode", 1, & contraction_3_1_meas.check_mode);
        read_iprms("random_offset", 1 , & contraction_3_1_meas.random_offset);
        read_iprms("kernel_type",1, & contraction_3_1_meas.kernel_type);
        read_dprms("am",1, & contraction_3_1_meas.am);
        read_dprms("a", 1, & contraction_3_1_meas.a);
    }
    //MPI_Bcast(& nprop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.nprop, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int npr = contraction_3_1_meas.nprop;
    if(rank != 0 ){contraction_3_1_meas.id_prop = (int*) malloc(npr*sizeof(int));}
    MPI_Bcast(contraction_3_1_meas.id_prop, npr, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(contraction_3_1_meas.vector_y, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.check_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.random_offset, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.kernel_type, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.am, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(& contraction_3_1_meas.a, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


//read file from arguments
//copy from Nils
//the syntax should be <app_name> -i <input file>
void read_hlbl31_infile(int argc, char *argv[])
{
  int ifile;
  int my_rank;
  FILE *fi;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {

    ifile = find_opt(argc, argv, "-i");

    error_root((ifile == 0)
               || (ifile == (argc - 1)), 1,
               "read_param_infile [read_hlbl_g_minus_2_parms.c]",
               "Syntax: <app_name> -i <input file> ");

    fi = freopen(argv[ifile + 1], "r", stdin);
    error_root(fi == NULL, 1,
               "read_param_infile [read_hlbl_g_minus_2_parms.c]",
               "Unable to open input file");
  }

  read_hlbl31();
  if (my_rank == 0)
    fclose(fi);
}

