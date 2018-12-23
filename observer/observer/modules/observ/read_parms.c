/**
* @file read_parms.c
* \brief  Reading parameters for the observere package.
* @author Dalibor Djukanovic
* @version 0.1
* @date 2014-06-04
*
* v0.1: Initial Version. Based on ms4 propgram in the openQCD package.
*
*/
#include <stdlib.h>
#include <time.h>
#include "utils.h"
#include "mpi.h"
#include "global.h"
#include "version.h"
#include "observ.h"



static FILE *fin = NULL, *flog = NULL;



static void read_dirs(void)
{
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0)
  {
    find_section("Run name");
    read_line("name", "%s", nbase);

    find_section("Directories");
    read_line("log_dir", "%s", log_dir);
    read_line("dat_dir", "%s", dat_dir);

    if (noexp)
    {
      read_line("loc_dir", "%s", loc_dir);
      cnfg_dir[0] = '\0';
    }
    else
    {
      read_line("cnfg_dir", "%s", cnfg_dir);
      loc_dir[0] = '\0';
    }

    read_line("sfld_dir", "%s", sfld_dir);

    find_section("Configurations");
    read_line("first", "%d", &(observer.first_cnfg));
    read_line("last", "%d", &(observer.last_cnfg));
    read_line("step", "%d", &(observer.step_cnfg));

    find_section("Random number generator");
    read_line("level", "%d", &(observer.level));
    read_line("seed", "%d", &(observer.seed));

    error_root((observer.last_cnfg < observer.first_cnfg)
               || (observer.step_cnfg < 1)
               ||
               (((observer.last_cnfg -
                  observer.first_cnfg) % observer.step_cnfg) != 0), 1,
               "read_dirs [read_infile.c]", "Improper configuration range");
  }

  MPI_Bcast(nbase, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(log_dir, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(loc_dir, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(cnfg_dir, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(dat_dir, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(sfld_dir, NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);

  MPI_Bcast(&observer.first_cnfg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&observer.last_cnfg, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&observer.step_cnfg, 1, MPI_INT, 0, MPI_COMM_WORLD);

  MPI_Bcast(&observer.level, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&observer.seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
}


static void setup_files(void)
{
  time_t now;
  struct tm *lcltime;
  char date[128];
  now = time(NULL);
  lcltime = localtime(&now);
  if (noexp)
    error_root(name_size
               ("%s/%sn%d_%d", loc_dir, nbase, observer.last_cnfg,
                NPROC - 1) >= NAME_SIZE, 1, "setup_files [read_parms.c]",
               "loc_dir name is too long");
  else
    error_root(name_size("%s/%sn%d", cnfg_dir, nbase, observer.last_cnfg) >=
               NAME_SIZE, 1, "setup_files [read_parms.c]",
               "cnfg_dir name is too long");

  /*check_dir_root(sfld_dir);
     error_root(name_size("%s/%sn%d.s%d",sfld_dir,nbase,last,nsrc-1)>=NAME_SIZE,
     1,"setup_files [ms4.c]","sfld_dir name is too long");
   */
  check_dir_root(log_dir);
  error_root(name_size("%s/%s.observ.log~", log_dir, nbase) >= NAME_SIZE,
             1, "setup_files [read_parms.c]", "log_dir name is too long");
  strftime(date, sizeof(date), "%Y-%m-%d_%H_%M_%S", lcltime);
  sprintf(log_file, "%s/%s.observ_%s.log", log_dir, nbase, date);
  sprintf(end_file, "%s/%s.observ.end", log_dir, nbase);
  sprintf(log_save, "%s~", log_file);
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Read lattice parameters from infile\n
* beta \t <DOUBLE> gauge coupling\n
* kappa \t <DOUBLE*> hopping parameters
* csw \t <DOUBLE> improvement coefficient
*/
/* ----------------------------------------------------------------------------*/
static void read_lat_parms(void)
{
  int nk;
  int my_rank;
  double beta, c0, csw, *kappa;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0)
  {
    find_section("Lattice parameters");
    read_line("beta", "%lf", &beta);
    read_line("c0", "%lf", &c0);
    nk = count_tokens("kappa");
    read_line("csw", "%lf", &csw);
  }

  MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&c0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nk, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&csw, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (nk > 0)
  {
    kappa = malloc(nk * sizeof(*kappa));
    error(kappa == NULL, 1, "read_lat_parms [qcd1.c]",
          "Unable to allocate parameter array");
    if (my_rank == 0)
      read_dprms("kappa", nk, kappa);
    MPI_Bcast(kappa, nk, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else
    kappa = NULL;

  set_lat_parms(beta, c0, nk, kappa, csw);

  if (nk > 0)
    free(kappa);

  /*if (append)
     check_lat_parms(fdat);
     else
     write_lat_parms(fdat);
   */
}



/* --------------------------------------------------------------------------*/
/**
* \brief  Read type of boundary condition form infile\n
*  - type: 
*   0: open boundary conditions. Needs\n
*    - cG
*    - cF
*   1: SF boundary conditions. Needs\n
*    - phi  
*    - phi'
*    - cG
*    - cF
*   2: open SF boundary conditions. Needs\n
*    - phi'
*    - cG
*    - cG'
*    - cF
*    - cF'
*   3: periodic boundary conditions.
*
*
*/
/* ----------------------------------------------------------------------------*/
static void read_bc_parms(void)
{
  int bc;
  double cF, cF_prime;
  double phi[2], phi_prime[2];
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0)
  {
    find_section("Boundary conditions");
    read_line("type", "%d", &bc);


    phi[0] = 0.0;
    phi[1] = 0.0;
    phi_prime[0] = 0.0;
    phi_prime[1] = 0.0;
    cF = 1.0;
    cF_prime = 1.0;

    if (bc == 1)
      read_dprms("phi", 2, phi);

    if ((bc == 1) || (bc == 2))
      read_dprms("phi'", 2, phi_prime);

    if (bc != 3)
      read_line("cF", "%lf", &cF);

    if (bc == 2)
      read_line("cF'", "%lf", &cF_prime);
  }

  MPI_Bcast(&bc, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(phi_prime, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cF, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cF_prime, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  bcp = set_bc_parms(bc, 1.0, 1.0, cF, cF_prime, phi, phi_prime);
}


static void read_sap_parms(void)
{
  int bs[4];
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {
    find_section("SAP");
    read_line("bs", "%d %d %d %d", bs, bs + 1, bs + 2, bs + 3);
  }

  MPI_Bcast(bs, 4, MPI_INT, 0, MPI_COMM_WORLD);
  set_sap_parms(bs, 1, 4, 5);
}


static void read_dfl_parms(void)
{
  int bs[4], Ns;
  int ninv, nmr, ncy, nkv, nmx;
  double kappa, mu, res;
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {
    find_section("Deflation subspace");
    read_line("bs", "%d %d %d %d", bs, bs + 1, bs + 2, bs + 3);
    read_line("Ns", "%d", &Ns);
  }

  MPI_Bcast(bs, 4, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Ns, 1, MPI_INT, 0, MPI_COMM_WORLD);
  set_dfl_parms(bs, Ns);

  if (my_rank == 0)
  {
    find_section("Deflation subspace generation");
    read_line("kappa", "%lf", &kappa);
    read_line("mu", "%lf", &mu);
    read_line("ninv", "%d", &ninv);
    read_line("nmr", "%d", &nmr);
    read_line("ncy", "%d", &ncy);
  }

  MPI_Bcast(&kappa, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&mu, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ninv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nmr, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ncy, 1, MPI_INT, 0, MPI_COMM_WORLD);
  set_dfl_gen_parms(kappa, mu, ninv, nmr, ncy);

  if (my_rank == 0)
  {
    find_section("Deflation projection");
    read_line("nkv", "%d", &nkv);
    read_line("nmx", "%d", &nmx);
    read_line("res", "%lf", &res);
  }

  MPI_Bcast(&nkv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nmx, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  set_dfl_pro_parms(nkv, nmx, res);
}



static void read_solver(void)
{
  solver_parms_t sp;
  int i, k;
  int isap, idfl;

  for (i = 0; i < observer.no_prop; i++)
  {
    k = prop[i].idx_solver;
    sp = solver_parms(k);
    if (sp.solver == SOLVERS)
    {
      read_solver_parms(k);
      sp = solver_parms(k);

      if (sp.solver == SAP_GCR)
        isap = 1;
      else if (sp.solver == DFL_SAP_GCR)
      {
        isap = 1;
        idfl = 1;
      }
    }
  }

  if (isap)
    read_sap_parms();

  if (idfl)
    read_dfl_parms();
}


static int max_int(int *arr, int size)
{
  int i;
  int tmp;

  tmp = size;
  for (i = 0; i < size; i++)
  {
    if (arr[i] > tmp)
      tmp = arr[i];
  }

  return tmp + 1;
}

static void read_obs_parms(void)
{
  int myrank;
  int i, j, max, id;
  char prop_sect[NAME_SIZE];
  char src_line[NAME_SIZE];

  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
  {
    observer.no_link_smearings=0;
    find_section("Observer");
    read_line("no_prop", "%d", &(observer.no_prop));
    read_line_optional("no_link_smearings", "%d", &(observer.no_link_smearings));
  }
  MPI_Bcast(&(observer.no_prop), 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&(observer.no_link_smearings), 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(observer.no_link_smearings > 0 )
  {
    /* error(1,1,"read_parms.c","For now a maximum of 1 link smearing is allowed");  */
    smearing_parms = amalloc(observer.no_link_smearings * sizeof(smearing_parms_t), ALIGN);
  }
  
  observer.idx_prop = amalloc(observer.no_prop * sizeof(int), ALIGN);
  if (myrank == 0)
  {
    read_iprms("idx_prop", observer.no_prop, &(observer.idx_prop[0]));
  }
  MPI_Bcast(&(observer.idx_prop[0]), observer.no_prop, MPI_INT, 0,
            MPI_COMM_WORLD);
  max = max_int(observer.idx_prop, observer.no_prop);
  observer.prop_list = amalloc(max * sizeof(prop_list_t), ALIGN);

  for (i = 0; i < max; i++)
  {
    for (j = 0; j < 12; j++)
    {
      observer.prop_list[i].spin_idx_status[j] = INIT;
      observer.prop_list[i].sd[j] = NULL;
    }
  }
  prop = amalloc(max * sizeof(prop_parms_t), ALIGN);

  if (myrank == 0)
  {
    find_section("Source positions");
    read_line("no_src_pos", "%d", &(observer.no_src_pos));
  }
  MPI_Bcast(&(observer.no_src_pos), 1, MPI_INT, 0, MPI_COMM_WORLD);
  


  observer.src_pos =
    (int **) amalloc(observer.no_src_pos * sizeof(int *), ALIGN);
  for (i = 0; i < observer.no_src_pos; i++)
  {
    observer.src_pos[i] = (int *) amalloc(4 * sizeof(int), ALIGN);
  }
  if (myrank == 0)
  {
    for (i = 0; i < observer.no_src_pos; i++)
    {
      sprintf(src_line, "pos%d", i);
      read_iprms(src_line, 4, observer.src_pos[i]);
    }
  }
  for (i = 0; i < observer.no_src_pos; i++)
  {
    MPI_Bcast(observer.src_pos[i], 4, MPI_INT, 0, MPI_COMM_WORLD);
  }

  for (i = 0; i < observer.no_link_smearings; i++)
  {
    id = i;
    /* Default is to smear in time direction */
    smearing_parms[id].skipdir=3;
    if (myrank == 0)
    {
      sprintf(prop_sect, "Linksmearing %d", id);
      find_section(prop_sect);
      read_line("smearing_t", "%s", &(smearing_parms[id].smearing_t));
      read_line_optional("smearing_steps", "%d", &(smearing_parms[id].smearing_steps));
      read_line_optional("smearing_parm", "%lf", &(smearing_parms[id].smearing_parm));
      read_line_optional("skipdir", "%d", &(smearing_parms[id].skipdir));
     }
     MPI_Bcast(&(smearing_parms[id].smearing_t), NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
     MPI_Bcast(&(smearing_parms[id].smearing_steps),1 , MPI_INT, 0 , MPI_COMM_WORLD);
     MPI_Bcast(&(smearing_parms[id].skipdir),1 , MPI_INT, 0 , MPI_COMM_WORLD);
     MPI_Bcast(&(smearing_parms[id].smearing_parm),1 , MPI_DOUBLE, 0 , MPI_COMM_WORLD);
  }
  for (i = 0; i < observer.no_prop; i++)
  {
    id = observer.idx_prop[i];
    if (myrank == 0)
    {
      sprintf(prop_sect, "Propagator %d", id);
      find_section(prop_sect);
      read_line("src_t", "%s", &(prop[id].src_t));
      /*read_iprms("src_pos",4,&(prop[id].src_pos[0])); */
      prop[id].src_pos[0]=0;
      prop[id].src_pos[1]=0;
      prop[id].src_pos[2]=0;
      prop[id].src_pos[3]=0;
      read_line("kappa", "%lf", &(prop[id].kappa));
      read_line("mu", "%lf", &(prop[id].mu));
      read_line("idx_solver", "%d", &(prop[id].idx_solver));
      /* Optional Parameters */
      prop[id].truncated_solver=0;
      prop[id].dump=-1;
      prop[id].smearing_steps_sink=-1;
      prop[id].smearing_steps_source=-1;
      prop[id].smearing_parm_sink=0.;
      prop[id].smearing_parm_source=0.;
      prop[id].link_smearing_idx=-1;
      read_line_optional("truncated_solver", "%d", &prop[id].truncated_solver);
      read_line_optional("dump","%d" , &(prop[id].dump));
      read_line_optional("link_smearing_idx","%d" , &(prop[id].link_smearing_idx));
      if ( prop[id].link_smearing_idx > observer.no_link_smearings )
        error(1,1,"read_parms.c","Link smearing %d not defined.", prop[id].link_smearing_idx);
      read_line_optional("smearing_steps_sink","%d" , &(prop[id].smearing_steps_sink));
      read_line_optional("smearing_parameter_sink","%lf",  &(prop[id].smearing_parm_sink));
      read_line_optional("smearing_steps_source","%d" , &(prop[id].smearing_steps_source));
      read_line_optional("smearing_parameter_source","%lf",  &(prop[id].smearing_parm_source));
    }
    MPI_Bcast(&(prop[id].src_t), NAME_SIZE, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(prop[id].src_pos, 4, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].kappa), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].mu), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].truncated_solver), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].dump), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].idx_solver), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].link_smearing_idx), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].smearing_steps_sink), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].smearing_steps_source), 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].smearing_parm_sink), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&(prop[id].smearing_parm_source), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
}


static void print_observer_parms(void)
{
  int myrank, i, id;
  char buffer[NAME_SIZE];
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (myrank == 0)
  {

    printf("--------------------------------\n");
    printf("Caluclation uses %d source positions\n", observer.no_src_pos);
    printf("++++++++++++++++++++++++++++++++\n");
    for (i = 0; i < observer.no_src_pos; i++)
    {
      printf("pos%d (t,x,y,z): ( %d, %d, %d, %d)\n", i, observer.src_pos[i][0],
             observer.src_pos[i][1], observer.src_pos[i][2],
             observer.src_pos[i][3]);
    }
    printf("++++++++++++++++++++++++++++++++\n");
    printf("--------------------------------\n");
    printf("Propagators:\n\n");
    for (i = 0; i < observer.no_prop; i++)
    {
      id = observer.idx_prop[i];
      printf("++++++++++++++++++++++++++++++++\n");
      printf("Propagator %d:\n", id);
      printf("++++++++++++++++++++++++++++++++\n");
      printf("Kappa: %e\n", prop[id].kappa);
      printf("Solver: %d\n", prop[id].idx_solver);
      if (prop[id].link_smearing_idx >= 0 )
        sprintf(buffer, "%d: (type %s, steps = %d, alpha = %lf, exclude mu = %d)",
            prop[id].link_smearing_idx,
            smearing_parms[prop[id].link_smearing_idx].smearing_t,
            smearing_parms[prop[id].link_smearing_idx].smearing_steps,
            smearing_parms[prop[id].link_smearing_idx].smearing_parm,
            smearing_parms[prop[id].link_smearing_idx].skipdir);
      else
        sprintf(buffer, "none");
      printf("Link smearing %s\n", buffer);
      if ( prop[id].smearing_steps_source > 0 )
        sprintf(buffer, "(type %s, steps = %d, alpha = %.10lf)",
            prop[id].src_t,
            prop[id].smearing_steps_source,
            prop[id].smearing_parm_source);
      else
        sprintf(buffer, "none");
      printf("Fermion field source smearing: %s\n", buffer);
      if ( prop[id].smearing_steps_sink > 0 )
        sprintf(buffer, "(type %s, steps = %d, alpha = %.10lf)",
            prop[id].src_t,
            prop[id].smearing_steps_sink,
            prop[id].smearing_parm_sink);
      else
        sprintf(buffer, "none");
      printf("Fermion field sink smearing: %s\n", buffer);
    }
    printf("--------------------------------\n");
  }
}

/* --------------------------------------------------------------------------*/
/**
* \brief  Read parameters from the infile and save them into global variables\n
* [Run name]\n
* <string> name: Run name\n
* [Directory]\n
* <string> log_dir Directory name where to put logfile with filename <cnfg_dir>/<nbase>.observ.log\n
* <string> cnfg_dir Directory name where configs are\n
* <string> dat_dir Directory name where to put data\n
* <string> loc_dir Directory name where to put local files\n
* <string> sfld_dir Directory name where to put local dumps of propagators.\n
* [Configurations]\n
* <int> first First config number used in <cnfg_dir>/<nbase>n<first>\n
* <int> last Last config number used in <cnfg_dir>/<nbase>n<last>\n
* <int> step Delta in the config number \n
* [Random number generator]\n
* <int> level Random nuber generator level\n
* <int> seed Seed for the random number generator\n
*
* \param argc
* \param argv[]
*/
/* ----------------------------------------------------------------------------*/
void read_infile(int argc, char *argv[])
{
  int ifile;
  /*int geom; */
  /*int lgeom; */
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0)
  {
    flog = freopen("STARTUP_ERROR", "w", stdout);

    ifile = find_opt(argc, argv, "-i");
    endian = endianness();



    error_root((ifile == 0) || (ifile == (argc - 1)), 1, "read_infile [ms4.c]",
               "Syntax: ms4 -i <input file> [-noexp]");

    error_root(endian == UNKNOWN_ENDIAN, 1, "read_infile [ms4.c]",
               "Machine has unknown endianness");

    noexp = find_opt(argc, argv, "-noexp");

    fin = freopen(argv[ifile + 1], "r", stdin);
    error_root(fin == NULL, 1, "read_infile [ms4.c]",
               "Unable to open input file");

    fseek(fin, 0, SEEK_END);
    long fsize = ftell(fin);
    fseek(fin, 0, SEEK_SET);

    observer.infile_contents = malloc((fsize + 1));
    fread(observer.infile_contents, fsize, 1, fin);
    fseek(fin, 0, SEEK_SET);
    observer.infile_contents[fsize] = 0;
  }


  MPI_Bcast(&endian, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&noexp, 1, MPI_INT, 0, MPI_COMM_WORLD);

  read_dirs();
  setup_files();
  read_lat_parms();
  read_bc_parms();
  read_obs_parms();
  read_solver();

  if (my_rank == 0)
    fclose(fin);
}


/* --------------------------------------------------------------------------*/
/**
* \brief  Print information of all read parameters
*/
/* ----------------------------------------------------------------------------*/
void print_info(void)
{
  int isap, idfl;
  /* int n[3]; */
  int n[3];
  long ip;
  int my_rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  if (my_rank == 0)
  {
    ip = ftell(flog);
    fclose(flog);

    if (ip == 0L)
      remove("STARTUP_ERROR");

    flog = freopen(log_file, "w", stdout);
    error_root(flog == NULL, 1, "print_info", "Unable to open log file");
    printf("\n");

    printf("Computation of observables\n");
    printf("--------------------------------\n\n");

    printf("OpenQCD version %s\n", openQCD_RELEASE);
    printf("Observer version %s\n", observer_RELEASE);
    printf("QDP version %s\n", qdp_RELEASE);
    printf("QMP version %s\n", qmp_RELEASE);

    if (endian == LITTLE_ENDIAN)
      printf("The machine is little endian\n");
    else
      printf("The machine is big endian\n");
    if (noexp)
      printf("Configurations are read in imported file format\n\n");
    else
      printf("Configurations are read in exported file format\n\n");

    printf("%dx%dx%dx%d lattice, ", NPROC0 * L0, NPROC1 * L1, NPROC2 * L2,
           NPROC3 * L3);
    printf("%dx%dx%dx%d local lattice\n", L0, L1, L2, L3);
    printf("%dx%dx%dx%d process grid, ", NPROC0, NPROC1, NPROC2, NPROC3);
    printf("%dx%dx%dx%d process block size\n",
           NPROC0_BLK, NPROC1_BLK, NPROC2_BLK, NPROC3_BLK);
    printf("SF boundary conditions on the quark fields\n\n");

    printf("Random number generator:\n");
    printf("level = %d, seed = %d\n\n", observer.level, observer.seed);


    if (bcp.type == 0)
    {
      printf("Open boundary conditions\n");

      n[0] = fdigits(bcp.cF[0]);
      printf("cF = %.*f\n\n", IMAX(n[0], 1), bcp.cF[0]);
    }
    else if (bcp.type == 1)
    {
      printf("SF boundary conditions\n");

      n[0] = fdigits(bcp.cF[0]);
      printf("cF = %.*f\n", IMAX(n[0], 1), bcp.cF[0]);

      n[0] = fdigits(bcp.phi[0][0]);
      n[1] = fdigits(bcp.phi[0][1]);
      n[2] = fdigits(bcp.phi[0][2]);
      printf("phi = %.*f,%.*f,%.*f\n", IMAX(n[0], 1), bcp.phi[0][0],
             IMAX(n[1], 1), bcp.phi[0][1], IMAX(n[2], 1), bcp.phi[0][2]);

      n[0] = fdigits(bcp.phi[1][0]);
      n[1] = fdigits(bcp.phi[1][1]);
      n[2] = fdigits(bcp.phi[1][2]);
      printf("phi' = %.*f,%.*f,%.*f\n\n", IMAX(n[0], 1), bcp.phi[1][0],
             IMAX(n[1], 1), bcp.phi[1][1], IMAX(n[2], 1), bcp.phi[1][2]);
    }
    else if (bcp.type == 2)
    {
      printf("Open-SF boundary conditions\n");

      n[0] = fdigits(bcp.cF[0]);
      printf("cF = %.*f\n", IMAX(n[0], 1), bcp.cF[0]);
      n[1] = fdigits(bcp.cF[1]);
      printf("cF' = %.*f\n", IMAX(n[1], 1), bcp.cF[1]);

      n[0] = fdigits(bcp.phi[1][0]);
      n[1] = fdigits(bcp.phi[1][1]);
      n[2] = fdigits(bcp.phi[1][2]);
      printf("phi' = %.*f,%.*f,%.*f\n\n", IMAX(n[0], 1), bcp.phi[1][0],
             IMAX(n[1], 1), bcp.phi[1][1], IMAX(n[2], 1), bcp.phi[1][2]);
    }
    else
      printf("Periodic boundary conditions\n\n");


    printf("Infile\n\
======\n%s\n", observer.infile_contents);
    print_solver_parms(&isap, &idfl);

    print_observer_parms();

    if (isap)
      print_sap_parms(0);

    if (idfl)
      print_dfl_parms(0);

    printf("Configurations no %d -> %d in steps of %d\n\n",
    observer.first_cnfg, observer.last_cnfg, observer.step_cnfg);
 
    fflush(flog);
  }
}

