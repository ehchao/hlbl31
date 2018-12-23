#define MAIN_PROGRAM
//#include <sstream>
//#include <cmath>
//#include <vector>
//#include "getL.h"
//#include <mpi.h>
#include "observer.h"
#include "contraction_3_1.h"
//extern 'C'{
//#include "KQED.h"
//}

using namespace QDP;
using namespace std;

//1 to use unit gauge config
#define UNIT_GAUGE 1
//1 to check lll, the UNIT_GAUGE must also be set to 1
#define CHECK_LLL 0
//1 to check llc
#define CHECK_LLC 0
//1 to check the Ward identity for lll
#define LLL_WI 0
//1 if we want to compute the relevant kernels; 0 if we want to compute the integrated 3 pt functions
#define COMPUTE_KERNELS 1
//1 to test the integration after mulitplied by z
#define TEST_Z_INTEGRAL 0
//1 to test the integration after mulitplied by the kernel L
#define TEST_L_INTEGRAL 0
//1 for production
#define COMPUTE_INTEGRALS 0

#define TEST_MODE CHECK_LLL || CHECK_LLC || LLL_WI || TEST_Z_INTEGRAL || TEST_L_INTEGRAL

//fix coordinates in the argument of a pure funtion which will be multiplied with a correlation funcion and integrated over
//eg., z-> z if z<L/2-1 and -> L - z else
/*void fix_coord(int out[4], const int crd[4]){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < 4 ; i++){
       if(crd[i] >= lattsize[i] / 2 ) out[i] = crd[i] - lattsize[i];
       else out[i] = crd[i];
    }
}

//assume that the length of the argument is 4
void fix_coord(multi1d<int> & out, const multi1d<int> & crd){
    multi1d<int> lattsize = Layout::lattSize();
    for(int i =0; i < 4 ; i++){ 
      if(crd[i] >= lattsize[i]/2) out[i] = crd[i] - lattsize[i];
      else out[i] = crd[i];
    }
}
*/

void array4d2multi4d(multi4d<RealD> &out, double in[6][4][4][4] , int i, int j, int k , int l){
  for(int a = 0; a < i; a ++){
    for(int b = 0; b < j; b ++){
      for(int c = 0; c < k; c ++){
        for(int d = 0; d < l; d ++){
          out[a][b][c][d] = Real(in[a][b][c][d]);
        }
      }
    }
  }
}

int main(int argc, char *argv[]){

    observer_init(argc, argv);
    read_hlbl31_infile(argc, argv);

    int status[3];

    StopWatch timer;
    multi1d<LatticeColorMatrixD> u(4);
    su3_dble *upoint;
	int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //for the check mode
    double check = 0;

    //flag for init kernelQED
    int flagkerqed = 0;
    
    //flag for kernel computation, needed for the production mode so that it will just be computed once
    int flagcompker = 0;
  
    multi1d<int> lattsize = Layout::lattSize();
    int pos_y[4], origin[4];

    //NAME_SIZE is defined in global.h
    char data_file[NAME_SIZE];
    //snprintf(data_file, NAME_SIZE, "%s/kernels.hlbl31.hdf5", dat_dir);
    //QDPIO::cout <<"Initializing the HDF5 file "<< data_file <<endl;
    //hdf5_writer writer(data_file);
    //hdf5_writer writer(data_file, true);


#if COMPUTE_KERNELS
    QDPIO::cout<<"[Kernel computing mode]" <<endl;
    snprintf(data_file, NAME_SIZE, "%s/kernels.hlbl31.hdf5", dat_dir);
    QDPIO::cout <<"Initializing the HDF5 file "<< data_file <<endl;
    hdf5_writer writer(data_file);
    QDPIO::cout << "Initilalization done" << endl ;

    if(contraction_3_1_meas.random_offset){
        if(rank == 0){
           double dr[4];
           ranlxd(dr,4);
           int tmp, tmpy;
           tmp = (int)(dr[3]*lattsize[3]) % lattsize[3];
           tmpy = (int)(dr[3]*lattsize[3] + contraction_3_1_meas.vector_y[3]) % lattsize[3];
           //we use open boundary condition, the inversion of the Dirac operator is impossible for a point source on the boundary
           while(tmp== 0 || tmp == lattsize[3]-1 || tmpy == 0 || tmpy == lattsize[3]-1){
               tmp = (tmp+1)%lattsize[3];
               tmpy = (tmpy+1)%lattsize[3];
           }
           origin[3] = tmp;
           pos_y[3] = tmpy;                    
           for(int i = 0; i < 3; i++){
               origin[i] = ( (int)(dr[i]*lattsize[i]) )% lattsize[i];
               pos_y[i] =( (int)(contraction_3_1_meas.vector_y[i] + dr[i]*lattsize[i]) )% lattsize[i];
           }
        }
       MPI_Bcast(&origin[0], 4, MPI_INT, 0 , MPI_COMM_WORLD);
       MPI_Bcast(&pos_y[0], 4, MPI_INT, 0 , MPI_COMM_WORLD);
    }
    else{
        for(int i =0; i< 4; i++){
            //put the origin in the middle in order not to be sensitive to the bc
            origin[i] = lattsize[i]/2;
            pos_y[i] = (contraction_3_1_meas.vector_y[i] + lattsize[i]/2) % lattsize[i];       
        }
    }
       
    timer.reset();
    timer.start();
    QDPIO::cout<<"Initializing QED kernel ... "<<endl;
    QED_kernel_temps *kert = (QED_kernel_temps*)malloc(3*sizeof(QED_kernel_temps));
    int ini = initialise(kert);

    multi1d<int> subgrid = QDP::Layout::subgridLattSize();
    multi1d<int> node_coord   = QDP::Layout::nodeCoord();
 
    QDPIO::cout<<"Initializing local kernel value storage ..."<<endl;
    double ker3L[6][4][4][4];
    double kerLx[4][4][4][4];
    double******** ker3L_loc = (double ********)malloc(6*sizeof(double*******));
    double******** kerLx_loc = (double ********)malloc(4*sizeof(double*******));
    for(int a = 0; a < 6; a++){
      ker3L_loc[a] = (double*******)malloc(4*sizeof(double******));
      for(int b = 0; b<4; b++){
        ker3L_loc[a][b] = (double******)malloc(4*sizeof(double*****));
        for(int c = 0; c<4; c++){
          ker3L_loc[a][b][c] = (double*****)malloc(4*sizeof(double****));
          for(int d = 0; d < 4; d++){
            ker3L_loc[a][b][c][d] = (double****)malloc(subgrid[0]*sizeof(double***));
            for(int e = 0; e < subgrid[0]; e++){
              ker3L_loc[a][b][c][d][e] = (double***)malloc(subgrid[1]*sizeof(double**));
              for(int f = 0 ; f< subgrid[1]; f++){
                ker3L_loc[a][b][c][d][e][f] = (double**)malloc(subgrid[2]*sizeof(double*));
                for(int g = 0; g< subgrid[2]; g++){
                  ker3L_loc[a][b][c][d][e][f][g] = (double*)malloc(subgrid[3]*sizeof(double));
                }
              }
            }
          }
        }
      }
    }
    for(int a = 0; a < 4; a++){
      kerLx_loc[a] = (double*******)malloc(4*sizeof(double******));
      for(int b = 0; b<4; b++){
        kerLx_loc[a][b] = (double******)malloc(4*sizeof(double*****));
        for(int c = 0; c<4; c++){
          kerLx_loc[a][b][c] = (double*****)malloc(4*sizeof(double****));
          for(int d = 0; d < 4; d++){
            kerLx_loc[a][b][c][d] = (double****)malloc(subgrid[0]*sizeof(double***));
            for(int e = 0; e < subgrid[0]; e++){
              kerLx_loc[a][b][c][d][e] = (double***)malloc(subgrid[1]*sizeof(double**));
              for(int f = 0 ; f< subgrid[1]; f++){
                kerLx_loc[a][b][c][d][e][f] = (double**)malloc(subgrid[2]*sizeof(double*));
                for(int g = 0; g< subgrid[2]; g++){
                  kerLx_loc[a][b][c][d][e][f][g] = (double*)malloc(subgrid[3]*sizeof(double));
                }
              }
            }
          }
        }
      }
    }

    double scalx[4], scaly[4];
    for(int i = 0 ; i < 4 ; i ++){ scalx[i] = contraction_3_1_meas.a; scaly[i] = contraction_3_1_meas.a;}

    QDPIO::cout<<"writing kernel"<<endl;
    int ss[4];
    for(int s = 0; s < 4; s++){ss[s] = node_coord[s]*subgrid[s];}
    for(int i = 0 ; i < subgrid[0]; i++){
      for(int j = 0 ; j < subgrid[1]; j++){
        for(int k = 0 ; k < subgrid[2]; k++){
          for(int l = 0 ; l < subgrid[3]; l++){
            int coord_x[4] = {i + ss[0],j + ss[1],k + ss[2],l +ss[3]};

            kernel_3L_new(ker3L, origin, coord_x, pos_y, scalx, scaly, kert);
            kernel_L_x_new(kerLx, origin, coord_x, pos_y, scalx, scaly, kert);

            for(int mu = 0 ; mu < 4; mu ++){
              for(int nu = 0 ; nu < 4; nu ++){
                for(int lambda = 0; lambda < 4; lambda ++){
                  for(int rhosig = 0 ; rhosig < 6; rhosig ++)
                    ker3L_loc[rhosig][mu][nu][lambda][i][j][k][l] = ker3L[rhosig][mu][nu][lambda];
                  for(int sigma = 0 ; sigma < 4; sigma ++)
                    kerLx_loc[sigma][mu][nu][lambda][i][j][k][l] = kerLx[sigma][mu][nu][lambda];
                }
              }
            }
          }
        }
      }
    }

    QDPIO::cout <<"Wrtiting 3L_kernel hdf5 ...";
    char kernel3L_file[NAME_SIZE], kernelLx_file[NAME_SIZE]; 
    snprintf(kernel3L_file, NAME_SIZE, "%s/3L_kernel_%d_%d_%d_%d.hdf5", dat_dir, contraction_3_1_meas.vector_y[0], contraction_3_1_meas.vector_y[1], contraction_3_1_meas.vector_y[2], contraction_3_1_meas.vector_y[3]);
    write_3L_kernel_hdf5_lattice(writer, kernel3L_file, ker3L_loc);
    QDPIO::cout << "done"<<endl; 
    QDPIO::cout << "Writing L_x_kernel hdf5 ...";
    snprintf(kernelLx_file, NAME_SIZE, "%s/L_x_kernel_%d_%d_%d_%d.hdf5", dat_dir, contraction_3_1_meas.vector_y[0], contraction_3_1_meas.vector_y[1], contraction_3_1_meas.vector_y[2], contraction_3_1_meas.vector_y[3]);
    write_L_x_kernel_hdf5_lattice(writer, kernelLx_file, kerLx_loc);
    QDPIO::cout << "done" <<endl;   

    timer.stop();
    QDPIO::cout << "3L_kernel and L_x_kernel writing finished ( "<< timer.getTimeInSeconds() << " s)"<<endl;

#else
#if COMPUTE_INTEGRALS
    //Compute the kernel once for all the configs     
    if(rank == 0 && flagkerqed == 0){
      callInitKernelQED(); 
      flagkerqed = 1;
    }
    multi4d<LatticeRealD> kernqed(6,4,4,4);
#endif

    int num_cnfg = observer.last_cnfg - observer.first_cnfg + 1;

    //iterate over configurations
    for(int icnfg = observer.first_cnfg; icnfg <= observer.last_cnfg; icnfg += observer.step_cnfg){
        snprintf(data_file, NAME_SIZE, "%s/%sn%d.hlbl_g_minus_2.hdf5", dat_dir, nbase, icnfg);
        QDPIO::cout <<"Initializing the HDF5 file "<< data_file <<endl;
        hdf5_writer writer(data_file);

        start_ranlux(observer.level, observer.seed ^ icnfg);
        sprintf(cnfg_file, "%s/%sn%d", cnfg_dir, nbase, icnfg);
#if !UNIT_GAUGE
        QDPIO::cout << "Importing config "<< cnfg_file << endl;
        import_cnfg(cnfg_file);
#else   
        QDPIO::cout << "Using unit gauge configuration"<<endl;
#endif

        //copy gauge links, copy from Nils
        timer.reset();
        timer.start();
        QDPIO::cout<< "Loading gauge configuration" <<endl;

        chs_ubnd(-1);
        copy_bnd_ud();
        dfl_modes(status);
        upoint = udfld();
        copy_from_openQCD(u, upoint);
 
        timer.stop();
        QDPIO::cout << "Finished gauge configuration loading ( " << timer.getTimeInSeconds() << " sec )" << endl;
        int pos_y[4];
        int origin[4];
 
        //origin with random offset
        if(contraction_3_1_meas.random_offset){
            if(rank == 0){
               double dr[4];
               ranlxd(dr,4);
               int tmp, tmpy;
               tmp = (int)(dr[3]*lattsize[3]) % lattsize[3];
               tmpy = (int)(dr[3]*lattsize[3] + contraction_3_1_meas.vector_y[3]) % lattsize[3];
               //we use open boundary condition, the inversion of the Dirac operator is impossible for a point source on the boundary
               while(tmp== 0 || tmp == lattsize[3]-1 || tmpy == 0 || tmpy == lattsize[3]-1){
                   tmp = (tmp+1)%lattsize[3];
                   tmpy = (tmpy+1)%lattsize[3];
               }
               origin[3] = tmp;
               pos_y[3] = tmpy;                    
               for(int i = 0; i < 3; i++){
                   origin[i] = ( (int)(dr[i]*lattsize[i]) )% lattsize[i];
                   pos_y[i] =( (int)(contraction_3_1_meas.vector_y[i] + dr[i]*lattsize[i]) )% lattsize[i];
               }
            }
           MPI_Bcast(&origin[0], 4, MPI_INT, 0 , MPI_COMM_WORLD);
           MPI_Bcast(&pos_y[0], 4, MPI_INT, 0 , MPI_COMM_WORLD);
        }
        else{
            for(int i =0; i< 4; i++){
                //put the origin in the middle in order not to be sensitive to the bc
                origin[i] = lattsize[i]/2;
                pos_y[i] = (contraction_3_1_meas.vector_y[i] + lattsize[i]/2) % lattsize[i];       
            }
        }
       
        multi1d<int> origin_multi1d(4);
        multi1d<int> pos_y_multi1d(4); 
        array2multi1d(origin_multi1d, & origin[0], 4);
        array2multi1d(pos_y_multi1d, & pos_y[0], 4);    
 
        //check mode on
        //if(contraction_3_1_meas.check_mode){
        
#if TEST_MODE
        //use only one flavor for test mode
        LatticePropagatorD s_0;
        LatticePropagatorD s_y;
        compute_propagator(s_0, origin_multi1d , 0);
        compute_propagator(s_y, pos_y_multi1d, 0);
#else  
        int nprop = contraction_3_1_meas.nprop;
        multi1d<LatticePropagatorD> s_0(nprop);
        multi1d<LatticePropagatorD> s_y(nprop);
        QDPIO::cout << "number of props : "<< nprop <<endl;
        for(int i = 0; i < nprop; i ++){
          compute_propagator(s_0[i], origin_multi1d , i);
          compute_propagator(s_y[i], pos_y_multi1d, i);
        }
        
#endif


        //QDPIO::cout <<"[Check mode] " << endl;
        timer.reset();
        timer.start();
        double increment = 0;
#if CHECK_LLC
        for(int nu =0; nu < 4; nu ++){
            for(int lambda = 0; lambda < 4; lambda ++){
                increment += check_3_pt_llc(u, s_0, s_y,  nu, lambda,  pos_y_multi1d);
            }
        }
        check += increment;
        timer.stop();
        QDPIO::cout << cnfg_file <<  " : residue for the Ward Identity : " << increment << " ( " << timer.getTimeInSeconds() << " sec)." <<endl;
#endif

#if CHECK_LLL && UNIT_GAUGE
        //in the routine of Harvey, the pole mass is used ( am^0 = am(1-b_m am), where b_m = -1/2 ), as well as the O(a) improved vector and axial currents (where b_V = b_A = 1 at tree level).
        //need to put kappa = 1/( 2(am (1-b_m am) + 4 )) in the infile
        multi1d<int> crd_x(4);
        //arbitrary x far from the boundary
        int crd_x_tmp[4] = {17,16,17,50};
		array2multi1d(crd_x, crd_x_tmp, 4);
        LatticePropagatorD gamma5_s0 = Gamma(15)*s_0;
        LatticePropagatorD gamma5_sy = Gamma(15)*s_y;
		check += check_3_pt_lll(u, gamma5_s0, gamma5_sy, crd_x, pos_y_multi1d,origin_multi1d , contraction_3_1_meas.am, contraction_3_1_meas.a);
        QDPIO::cout<< "Computation for lll check finished " <<endl;
        timer.stop();
        QDPIO::cout << "3 point lll function compared to the C code of Harvey : " << check << " ( " << timer.getTimeInSeconds() << " sec)." <<endl;
#endif

#if (TEST_Z_INTEGRAL || TEST_L_INTEGRAL) && UNIT_GAUGE
        multi1d<int> crd_x(4);
        LatticePropagatorD gamma5_s0 = Gamma(15)*s_0;
        LatticePropagatorD gamma5_sy = Gamma(15)*s_y;
        multi3d<LatticeComplexD> corr_c(4,4,4);
        multi3d<LatticeRealD> corr(4,4,4);
        QDPIO::cout << "writing the 3 pt (lll) correlation function" <<endl;
        for(int mu = 0; mu < 4 ; mu ++){
          for(int nu = 0; nu < 4 ; nu ++){
            for(int lambda = 0 ; lambda < 4; lambda ++){
              write_3_pt_lll(corr_c, gamma5_s0, gamma5_sy, mu, nu, lambda, pos_y_multi1d);
              corr[mu][nu][lambda] = real(corr_c[mu][nu][lambda]);
            }
          }
        }
#endif
#if TEST_Z_INTEGRAL && UNIT_GAUGE
        double z3pt[4][4][4][4];
        timer.reset();
        timer.start();
        QDPIO::cout<<"computing z integral (lll)" << endl;
        int_3pt_times_func(z3pt, corr, sin_func, origin);
        QDPIO::cout<<"z3pt values for (rho,mu,nu,lambda): " <<endl;
        for(int rho = 0 ; rho < 4; rho ++)
          for(int mu = 0 ; mu < 4; mu ++)
            for(int nu = 0 ; nu < 4; nu ++)
              for(int sigma = 0 ; sigma < 4; sigma ++)
                QDPIO::cout << rho<<mu<<nu<<sigma<<" " << z3pt[rho][mu][nu][sigma]<<endl; 
        timer.stop();
        QDPIO::cout<< "time cost: "<<timer.getTimeInSeconds()<<endl;
#endif

#if TEST_L_INTEGRAL && UNIT_GAUGE
        timer.reset();
        timer.start();
        double testlint[6];
        QDPIO::cout<< " Computing int_L_3pt"<<endl;
        if(rank == 0 && flagkerqed ==0) { 
          callInitKernelQED(); 
          flagkerqed = 1;
        }
        int_L_3pt_loop( testlint, corr, origin_multi1d, pos_y_multi1d, contraction_3_1_meas.a);
      
        QDPIO::cout<<"Values with rhosig " <<endl;
        for(int rho = 0; rho <6 ; rho++) QDPIO::cout << rho<<" " << testlint[rho]<<endl; 
        timer.stop();
        QDPIO::cout<< "time cost: "<<timer.getTimeInSeconds()<<endl;

#endif

/*---from here, codes for production--------------*/
#if COMPUTE_INTEGRALS
        double z3pt[4][4][4][4];
        multi3d<LatticeRealD> corr(4,4,4);
        for(int iprop = 0; iprop < nprop; iprop ++){
          QDPIO::cout<< "[Data for propagator "<< iprop <<"]"<<endl;
          timer.reset();
          timer.start();
          QDPIO::cout << "writing the 3 pt (llc) correlation function" <<endl;
          for(int nu = 0; nu < 4 ; nu ++){
            for(int lambda = 0; lambda < 4 ; lambda ++){
              for(int sigma = 0 ; sigma < 4; sigma ++){
              write_im_3_pt_llc(corr, u, s_0[iprop], s_y[iprop], nu, lambda, sigma, pos_y_multi1d);
              }
            }
          }
         
          QDPIO::cout<<"Computing z integral (llc) ... " << endl;
          int_3pt_times_func(z3pt, corr, id_func, origin);
         
          QDPIO::cout<<"Computing integral (llc) ..."<<endl;
          double intloop3pt[4][4][4];
          for(int nu = 0; nu < 4; nu ++)
            for(int lambda = 0 ; lambda < 4; lambda ++)
              for(int sigma = 0; sigma < 4; sigma ++)
                intloop3pt[nu][lambda][sigma] = toWordType(sum(corr[nu][lambda][sigma]));
          
         
          QDPIO::cout<< "Computing L-3pt (lll) integral ... "<<endl;
          double intL3pt[6];
          multi3d<LatticeComplexD> corlll_c(4,4,4);
          multi3d<LatticeRealD> corlll(4,4,4);
          for(int mu =0; mu < 4 ; mu ++){
            for(int nu = 0 ; nu < 4; nu ++){
              for(int lambda = 0; lambda < 4; lambda ++){
                write_3_pt_lll(corlll_c, s_0[iprop], s_y[iprop],  mu,  nu,  lambda, pos_y_multi1d);
                corlll[mu][nu][lambda] = imag(corlll_c[mu][nu][lambda]);
              }
            }
          }

          //Compute kernel, just called once
          if(flagcompker == 0){
            kernel_L(kernqed, origin, pos_y, contraction_3_1_meas.a);
            flagcompker = 1;
          }
          for(int rhosig = 0; rhosig <6; rhosig ++){
            for(int mu = 0; mu < 4; mu ++){
              for(int nu = 0; nu < 4; nu ++){
                for(int lambda = 0 ; lambda < 4; lambda ++){
                  intL3pt[rhosig] = toWordType(sum(kernqed[rhosig][mu][nu][lambda]*corlll[mu][nu][lambda])); 
                }
              }
            }
          }
          //int_L_3pt_loop(intL3pt, corlll, origin_multi1d, pos_y_multi1d, contraction_3_1_meas.a);
         
         
          QDPIO::cout<< "***************************" << endl;
          QDPIO::cout<<"[Result for int_z_3pt]"<<endl;
          QDPIO::cout<<"origin ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< origin_multi1d[i] << " ";}
          QDPIO::cout<<"vector_y ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< pos_y_multi1d[i] << " ";}
          QDPIO::cout<<endl; 
          QDPIO::cout<<"z3pt values for (rho,nu,lambda,sigma): " <<endl;
          for(int rho = 0 ; rho < 4; rho ++)
            for(int mu = 0 ; mu < 4; mu ++)
              for(int nu = 0 ; nu < 4; nu ++)
                for(int sigma = 0 ; sigma < 4; sigma ++)
                  QDPIO::cout << rho<<mu<<nu<<sigma<<" " << z3pt[rho][mu][nu][sigma]<<endl; 
         
          QDPIO::cout<< "***************************" << endl;
          QDPIO::cout<<"[Result for int_3pt]"<<endl;
          QDPIO::cout<<"origin ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< origin_multi1d[i] << " ";}
          QDPIO::cout<<"vector_y ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< pos_y_multi1d[i] << " ";}
          QDPIO::cout<<endl; 
          QDPIO::cout<<"int_3pt values for (nu,lambda,sigma): " <<endl;
          for(int nu = 0; nu < 4; nu ++)
            for(int lambda = 0 ; lambda < 4; lambda ++)
              for(int sigma = 0; sigma < 4; sigma ++)
                QDPIO::cout << nu << lambda << sigma << " " << intloop3pt[nu][lambda][sigma]<<endl;
         
         
          QDPIO::cout<< "***************************" << endl;
          QDPIO::cout<<"[Result for int_L_3pt]"<<endl;
          QDPIO::cout<<"origin ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< origin_multi1d[i] << " ";}
          QDPIO::cout<<"vector_y ";
          for(int i = 0; i < 4; i ++){QDPIO::cout<< pos_y_multi1d[i] << " ";}
          QDPIO::cout<<endl; 
          QDPIO::cout<<"int_L_3pt values for rhosig : " <<endl;
          for(int rho = 0 ; rho < 6; rho ++)
            QDPIO::cout << rho<<" " << intL3pt[rho]<<endl; 
         
         
          timer.stop();
          QDPIO::cout<< "time cost: "<<timer.getTimeInSeconds()<<endl;
          QDPIO::cout<< "***************************" << endl;
        }
#endif        

        //}
#if 0
        //computing mode
        else{
          multi3d<LatticeRealD> corr_3pt_llc(4,4,4);
          multi3d<LatticeComplexD> corr_3pt_lll(4,4,4);
          multi3d<LatticeRealD> corr_3pt_lll_imag(4,4,4);
          multi4d<RealD> z_3pt_int(4,4,4,4);
          multi3d<RealD> int_3pt(4,4,4);
          multi1d<RealD> kernel_3pt_int(6);

          for(int nu = 0; nu < 4; nu ++){
            for(int lambda = 0; lambda < 4; lambda ++){
              for(int sigma = 0; sigma < 4; sigma ++){
                write_im_3_pt_llc(corr_3pt_llc, u, s_0, s_y,  nu, lambda, sigma, pos_y_multi1d);

              }
            }
          }
          for(int mu = 0; mu < 4 ; mu++){
            for(int nu = 0; nu < 4; nu ++){
              for(int lambda = 0; lambda < 4; lambda ++){
                write_3_pt_lll(corr_3pt_lll, s_0, s_y, mu, nu, lambda, pos_y_multi1d );
                corr_3pt_lll_imag[mu][nu][lambda] = real(corr_3pt_lll[mu][nu][lambda]);
              }
            }
          }
          int_L_3pt_loop(kernel_3pt_int, corr_3pt_lll_imag, pos_y_multi1d );


          for(int rho = 0; rho < 4; rho ++){
            LatticeRealD z_rho = Layout::latticeCoordinate(rho);
            for(int nu = 0 ; nu < 4; nu ++){
              for(int lambda = 0 ; lambda < 4 ; lambda ++){
                for(int sigma = 0 ; sigma < 4; sigma ++){
                  z_3pt_int[rho][nu][lambda][sigma] = sum(corr_3pt_llc[nu][lambda][sigma]*z_rho);
                }
              }
            } 
          }

          QDPIO::cout << "Writing HDF5 files"<<endl;
          write_int_z_3pt_loop_hdf5(writer, pos_y[0], pos_y[1], pos_y[2], pos_y[3], 3pt_z_int); 
          3pt_int[nu][lambda][sigma] = sum(3pt_corr[nu][lambda][sigma]);
          write_int_3pt_loop_hdf5(writer, pos_y[0], pos_y[1], pos_y[2], pos_y[3], 3pt_int);
          write_int_kernel_3pt_loop_hdf5(writer, pos_y[0], pos_y[1], pos_y[2], pos_y[3], 3pt_kernel_int );
       }
#endif

    }
#endif

    if(contraction_3_1_meas.check_mode){
#if CHECK_LLC
      QDPIO::cout<< "The Ward identity is verified at precision : "<< check <<endl;
#endif

    }

  
    Layout::destroy();
    QDPIO::cout<<"finalizing qdp"<<endl;
    QDP_finalize();
    QDPIO::cout <<"End" <<endl;
    return 0;
}
