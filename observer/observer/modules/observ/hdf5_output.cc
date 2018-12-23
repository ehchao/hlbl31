#include <cassert>
#include <mpi.h>
#include <hdf5.h>
#include <hdf5_hl.h>
/*#include "hdf5_output.h"
extern "C"
{
#include "observ.h"
}
*/
#include "observer.h"
/* assumes complex numbers are stored in a particular way
 * (the only reasonable way of doing so) */
static hid_t complex_type()
{
  hid_t complex_t;

  complex_t = H5Tcreate(H5T_COMPOUND, 2 * sizeof(double));
  H5Tinsert(complex_t, "r", 0, H5T_NATIVE_DOUBLE);
  H5Tinsert(complex_t, "i", sizeof(double), H5T_NATIVE_DOUBLE);
  return complex_t;
}

void hdf5_writer::init(const char *filename, bool mpi)
{
  use_mpi = mpi;
  if (!mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t fapl = H5P_DEFAULT;
  if (use_mpi) {
    fapl = H5Pcreate(H5P_FILE_ACCESS);
    // see main page for MPI_File_open for "info" options */
    H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  }

  file = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  if (file < 0)
  {
    error_root(1, 1, "[hdf5_output.cc]",
               "File %s could not be created.\nFile may already exist!",
               filename);
  }

  complex_t = complex_type();
  axes_group = H5Gcreate(file, "/axes",H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  //H5LTmake_dataset_string(file, "/infile", observer.infile_contents);
  //H5LTmake_dataset_string(file, "/observer_version", observer_RELEASE);
  //H5LTmake_dataset_string(file, "/qdp_version", qdp_RELEASE);
  //H5LTmake_dataset_string(file, "/qmp_version", qmp_RELEASE);

  if (use_mpi) {
    H5Pclose(fapl);

    int mynode = QDP::Layout::nodeNumber();
    QDP::multi1d<int> subgrid = QDP::Layout::subgridLattSize();
    int nsites = QDP::Layout::sitesOnNode();
    mapping = new hsize_t[nsites];
    for (int i = 0; i < nsites; i++) {
      QDP::multi1d<int> coord = QDP::Layout::siteCoords(mynode, i);
      for (int j = 0; j < 4; j++)
	coord[j] %= subgrid[j];
      int ix = 0;
      for (int j = 0; j < 4; j++)
	ix = ix * subgrid[j] + coord[j];
      mapping[ix] = i;
    }
  } else
    mapping = NULL;
}

hdf5_writer::hdf5_writer(const char *filename, bool mpi)
{
  init(filename, mpi);
}

hdf5_writer::hdf5_writer(const char *prefix, int cnfg_no, bool mpi)
{
  char export_corr_name[NAME_SIZE];
  snprintf(export_corr_name, sizeof(export_corr_name), "%s/%s_%sn%d",
           dat_dir, prefix, nbase, cnfg_no);
  init(export_corr_name, mpi);
}

hdf5_writer::~hdf5_writer()
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  H5Gclose(axes_group);
  H5Tclose(complex_t);
  H5Fclose(file);
  if (mapping != NULL)
    delete [] mapping;
}

void hdf5_writer::write_axis_info(const char *name, int count, int size,
                                  const int *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t dataset, dataspace, dataset_parms;
  hsize_t dims[2];
  dims[0] = count;
  dims[1] = size;

  dataset_parms = H5Pcreate(H5P_DATASET_CREATE);
  // write axis info as a single chunk
  if (size > 1) {
    dataspace = H5Screate_simple(2, dims, NULL);
    H5Pset_chunk(dataset_parms, 2, dims);
  } else {
    dataspace = H5Screate_simple(1, dims, NULL);
    H5Pset_chunk(dataset_parms, 1, dims);
  }
  H5Pset_fletcher32(dataset_parms);
  dataset =
    H5Dcreate(axes_group, name, H5T_NATIVE_INT, dataspace, H5P_DEFAULT,
              dataset_parms, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5DSset_scale(dataset, NULL);

  H5Pclose(dataset_parms);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void hdf5_writer::write_axis_info(const char *name, int count, int size,
                                  const double *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t dataset, dataspace, dataset_parms;
  hsize_t dims[2];
  dims[0] = count;
  dims[1] = size;

  dataset_parms = H5Pcreate(H5P_DATASET_CREATE);
  // write axis info as a single chunk
  if (size > 1) {
    dataspace = H5Screate_simple(2, dims, NULL);
    H5Pset_chunk(dataset_parms, 2, dims);
  } else {
    dataspace = H5Screate_simple(1, dims, NULL);
    H5Pset_chunk(dataset_parms, 1, dims);
  }
  H5Pset_fletcher32(dataset_parms);
  dataset =
    H5Dcreate(axes_group, name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT,
              dataset_parms, H5P_DEFAULT);

  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5DSset_scale(dataset, NULL);

  H5Pclose(dataset_parms);
  H5Dclose(dataset);
  H5Sclose(dataspace);
}

void hdf5_writer::write_axis_info(const char *name, int count,
                                  const char *const *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t dataset, dataspace, datatype;
  hsize_t dim = count;

  dataspace = H5Screate_simple(1, &dim, NULL);
  datatype = H5Tcopy(H5T_C_S1);
  H5Tset_size(datatype, H5T_VARIABLE);
  dataset =
    H5Dcreate(axes_group, name, datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT,
              H5P_DEFAULT);

  H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5DSset_scale(dataset, NULL);

  H5Dclose(dataset);
  H5Tclose(datatype);
  H5Sclose(dataspace);
}

hid_t hdf5_writer::write_generic_init_0(const char *name, int ndim,
					const int *dims,
					const char *const *axis_labels,
					const char *const *axis_info, hid_t type)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return -1;
  hid_t dataset, dataspace, dataset_parms;
  hsize_t *hdims = new hsize_t[ndim];
  hsize_t *cdims = new hsize_t[ndim];

  for (int i = 0; i < ndim; i++)
    cdims[i] = hdims[i] = dims[i];
  if (use_mpi) {
    QDP::multi1d<int> subgrid = QDP::Layout::subgridLattSize();
    int j;
    for (j = 0; j < ndim-4; j++) cdims[j]=1;
    for (int i = 0; j < ndim; i++, j++) cdims[j] = subgrid[i];
  } else {
    // combine the last two dimensions into a single chunk
    for (int i = 0; i < ndim-2; i++)
      cdims[i] = 1;
  }

  dataspace = H5Screate_simple(ndim, hdims, NULL);
  dataset_parms = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dataset_parms, ndim, cdims);
  H5Pset_fletcher32(dataset_parms);
  dataset =
    H5Dcreate(file, name, type, dataspace, H5P_DEFAULT, dataset_parms,
              H5P_DEFAULT);

  if (axis_labels != NULL)
    for (int i = 0; i < ndim; i++)
      if (axis_labels[i] != NULL)
        H5DSset_label(dataset, i, axis_labels[i]);
  if (axis_info != NULL)
    for (int i = 0; i < ndim; i++)
      if (axis_info[i] != NULL)
      {
        hid_t info = H5Dopen(axes_group, axis_info[i], H5P_DEFAULT);
        H5DSattach_scale(dataset, info, i);
        H5Dclose(info);
      }
  H5Pclose(dataset_parms);
  H5Sclose(dataspace);
  delete[]hdims;
  delete[]cdims;
  return dataset;
}

void hdf5_writer::write_generic_init(const char *name, int ndim,
				     const int *dims,
				     const char *const *axis_labels,
				     const char *const *axis_info, hid_t type)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t dset = write_generic_init_0(name, ndim, dims,
				    axis_labels, axis_info, type);
  H5Dclose(dset);  
}

void hdf5_writer::write_real_init(const char *name, int ndim, const int *dims,
                                  const char *const *axis_labels,
                                  const char *const *axis_info)
{
  write_generic_init(name, ndim, dims, axis_labels, axis_info,
                     H5T_NATIVE_DOUBLE);
}

void hdf5_writer::write_complex_init(const char *name, int ndim,
                                     const int *dims,
                                     const char *const *axis_labels,
                                     const char *const *axis_info)
{
  write_generic_init(name, ndim, dims, axis_labels, axis_info, complex_t);
}

void hdf5_writer::write_data_generic_0(hid_t dataset, const int *dims_offset,
				       const int *dims_data, int ndim_data,
				       hid_t type, const void *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;

  hid_t dataspace = H5Dget_space(dataset);
  int ndim = H5Sget_simple_extent_ndims(dataspace);

  hsize_t *start = new hsize_t[ndim];
  for (int i = 0; i < ndim; i++)
    start[i] = dims_offset[i];

  hsize_t *spacesize = new hsize_t[ndim];
  H5Sget_simple_extent_dims(dataspace, spacesize, NULL);
  hsize_t *count = new hsize_t[ndim];
  for (int i = 0; i < ndim; i++)
    count[i] = 1;
  for (int i = 0; i < ndim_data; i++)
    count[dims_data[i]] = spacesize[dims_data[i]];
  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
  hid_t memspace = H5Screate_simple(ndim, count, NULL);

  H5Dwrite(dataset, type, memspace, dataspace, H5P_DEFAULT, data);

  H5Sclose(memspace);
  H5Sclose(dataspace);
  delete[]start;
  delete[]spacesize;
  delete[]count;
}

void hdf5_writer::write_data_generic(const char *name, const int *dims_offset,
				     const int *dims_data, int ndim_data,
				     hid_t type, const void *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;

  hid_t dataset = H5Dopen(file, name, H5P_DEFAULT);
  write_data_generic_0(dataset, dims_offset, dims_data, ndim_data, type, data);
  H5Dclose(dataset);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim,
                             const QDP::multi1d < QDP::ComplexD > &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  write_data_generic(name, dims_offset, &dim, 1, complex_t, &data[0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0,
                             int dim1,
                             const QDP::multi2d < QDP::ComplexD > &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0<dim1);
  int dims_data[] = {dim0, dim1};
  write_data_generic(name, dims_offset, dims_data, 2, complex_t, &data[0][0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0,
                             int dim1, int dim2,
                             const QDP::multi3d < QDP::ComplexD > &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0 < dim1);
  assert(dim1 < dim2);
  int dims_data[] = { dim0, dim1, dim2 };
  write_data_generic(name, dims_offset, dims_data, 3, complex_t,
                     &data[0][0][0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0, int dim1, int dim2, int dim3, const QDP::multi4d<QDP::ComplexD> &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0<dim1);
  assert(dim1<dim2);
  assert(dim2<dim3);
  int dims_data[] = {dim0, dim1, dim2, dim3};
  write_data_generic(name, dims_offset, dims_data, 4, complex_t, &data[0][0][0][0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim, const QDP::multi1d<QDP::RealD> &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  write_data_generic(name, dims_offset, &dim, 1, H5T_NATIVE_DOUBLE, &data[0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0,
                             int dim1, const QDP::multi2d < QDP::RealD > &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0 < dim1);
  int dims_data[] = { dim0, dim1 };
  write_data_generic(name, dims_offset, dims_data, 2, H5T_NATIVE_DOUBLE,
                     &data[0][0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0,
                             int dim1, int dim2,
                             const QDP::multi3d < QDP::RealD > &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0 < dim1);
  assert(dim1 < dim2);
  int dims_data[] = { dim0, dim1, dim2 };
  write_data_generic(name, dims_offset, dims_data, 3, H5T_NATIVE_DOUBLE,
                     &data[0][0][0]);
}

void hdf5_writer::write_data(const char *name, const int *dims_offset, int dim0, int dim1, int dim2, int dim3, const QDP::multi4d<QDP::RealD> &data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  assert(dim0<dim1);
  assert(dim1<dim2);
  assert(dim2<dim3);
  int dims_data[] = {dim0, dim1, dim2, dim3};
  write_data_generic(name, dims_offset, dims_data, 4, H5T_NATIVE_DOUBLE, &data[0][0][0][0]);
}


//write multi4d lattice data on node 0
//size"i" is the i-th size of the multi4d data
void hdf5_writer::write_lattice_data_locally(const char *name, int *dims_offset, int size0, int size1, int size2, int size3, const QDP::multi4d<QDP::LatticeRealD> &data){
  //if (QDP::Layout::nodeNumer()!=0)
  //  return;
  int rank, ndim;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  QDP::multi1d<int> lattsize = QDP::Layout::lattSize();
  hid_t dataset, dataspace, memspace;
  int dims_data[] = {0, 1, 2, 3, 4, 5, 6, 7 };

  if(rank == 0){
    dataset = H5Dopen(file, name, H5P_DEFAULT);
    dataspace = H5Dget_space(dataset);
    ndim = H5Sget_simple_extent_ndims(dataspace);
  }
  MPI_Bcast(&ndim, 1, MPI_INT, 0 , MPI_COMM_WORLD);
  hsize_t *start = new hsize_t[ndim];
  hsize_t *spacesize = new hsize_t[ndim];
  hsize_t *count = new hsize_t[ndim];
  if(rank == 0){
    for (int i = 0; i < ndim; i++)
      start[i] = dims_offset[i];
    H5Sget_simple_extent_dims(dataspace, spacesize, NULL);
    for (int i = 0; i < 4; i++)
      count[i] = 1;
    for (int i = 4; i < 8; i++)
      count[dims_data[i]] = spacesize[dims_data[i]];
    //H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
    memspace = H5Screate_simple(ndim, count, NULL);
  }
  double data_d[lattsize[0]][lattsize[1]][lattsize[2]][lattsize[3]];

  for(int a = 0; a < size0; a ++){
    for(int b = 0; b < size1; b ++){
      for(int c = 0; c < size2; c ++){
        for(int d = 0; d < size3; d ++){
          QDP::QDPIO::cout << "writing hdf5 data " <<a << b << c<< d<< endl;
          for(int x0 =0; x0 < lattsize[0]; x0 ++ ){
            for(int x1 = 0; x1<lattsize[1]; x1++){
              for(int x2 =0; x2< lattsize[2]; x2++){
                for(int x3= 0; x3< lattsize[3]; x3++){
                  int tmp[4] = {x0, x1, x2, x3};
                  QDP::multi1d<int> tmpm(4);
                  for(int i = 0; i < 4; i ++){tmpm[i] = tmp[i];}
                  data_d[x0][x1][x2][x3] = QDP::toWordType(QDP::peekSite(data[a][b][c][d], tmpm));
                }
              }
            }
          }
          if(rank == 0){
            start[0] = a;
            start[1] = b;
            start[2] = c;
            start[3] = d;
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(data_d[0][0][0][0]));
          }
        }
      }
    }
  }
  if(rank == 0){ 
    H5Sclose(memspace);
    H5Sclose(dataspace);
    H5Dclose(dataset);
  }
  delete[]start;
  delete[]spacesize;
  delete[]count;

  //delete[]data_d;
}

//write lattice data with multiple processor
//4d array data with the four last positions reserved for the lattice coordinates
//the data should be attribute to each node in a particular way
void hdf5_writer::write_lattice_data_parallel(const char *file_name, const char *dset_name, const int *dims_offset, int size0, int size1, int size2, int size3,  double ******** ker_loc)
{ 
  int rank, ndim;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  QDP::multi1d<int> lattsize = QDP::Layout::lattSize();
  QDP::multi1d<int> subgrid = QDP::Layout::subgridLattSize();
  QDP::multi1d<int> coord   = QDP::Layout::nodeCoord();
  hid_t dataset, dataspace, memspace;
  int dims_data[] = {0, 1, 2, 3, 4, 5, 6, 7 };
  hid_t  file_id, dset_id, filespace;  
  hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL);
  file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
  H5Pclose(plist_id);

  hsize_t dimsf[8];
  dimsf[0] = size0;
  dimsf[1] = size1;
  dimsf[2] = size2;
  dimsf[3] = size3;
  for(int i =0; i < 4; i++) dimsf[4+i] = lattsize[i];
  int RANK = 8;

  filespace = H5Screate_simple(RANK,dimsf, NULL);
  dset_id = H5Dcreate(file_id, dset_name, H5T_NATIVE_DOUBLE, filespace,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  H5Sclose(filespace);

  hsize_t *start = new hsize_t[RANK];
  hsize_t *spacesize = new hsize_t[RANK];
  hsize_t *count = new hsize_t[RANK];
  for (int i = 0; i < 4; i++) start[i] = dims_offset[i];
  for (int i = 4; i < 8; i++) start[i] = subgrid[i-4] * coord[i-4];

  for (int i = 0; i < 4; i++) count[i] = 1;
  for (int i = 4; i < 8; i++) count[i] = subgrid[i-4];
  
  memspace = H5Screate_simple(RANK, count, NULL);
  filespace = H5Dget_space(dset_id);

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  double data_d[subgrid[0]][subgrid[1]][subgrid[2]][subgrid[3]];

  int mynode = QDP::Layout::nodeNumber();
  int numsites = QDP::Layout::sitesOnNode();

  QDP::StopWatch timer;
  timer.reset();
  timer.start();
  for(int a = 0; a < size0; a ++){
    for(int b = 0; b < size1; b ++){
      for(int c = 0; c < size2; c ++){
        for(int d = 0; d < size3; d ++){
          QDP::QDPIO::cout << "writing hdf5 data " <<a << b << c<< d<< endl;
            for(int x0 =0; x0 < subgrid[0]; x0 ++ ){
              for(int x1 = 0; x1< subgrid[1]; x1 ++){
                for(int x2 =0; x2< subgrid[2]; x2 ++){
                  for(int x3= 0; x3< subgrid[3]; x3 ++){
                    data_d[x0][x1][x2][x3] = ker_loc[a][b][c][d][x0][x1][x2][x3];
                  }
                }
              }
            }

            start[0] = a;
            start[1] = b;
            start[2] = c;
            start[3] = d;
            H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, dxpl, &(data_d[0][0][0][0]));
        }
      }
    }
  }
  timer.stop();
  QDP::QDPIO::cout<<"hdf5 writing finished ( "<<timer.getTimeInSeconds() << " s)"<<endl;
  
  H5Pclose(dxpl);
  H5Sclose(memspace);
  H5Dclose(dset_id);

  H5Sclose(filespace);
  H5Fclose(file_id);

  delete [] start;
  delete [] count;
}

/* uses the last four dimensions for the lattice */
/*
void hdf5_writer::write_lattice_generic(const char *name, const int *dims_offset, hid_t type, const void *data)
{
  int i;
  hid_t dataset = H5Dopen(file, name, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);

  int ndim = H5Sget_simple_extent_ndims(dataspace);
  assert(ndim >= 4);

  QDP::multi1d<int> subgrid = QDP::Layout::subgridLattSize();
  QDP::multi1d<int> coord   = QDP::Layout::nodeCoord();

  hsize_t *start = new hsize_t[ndim];
  for (i = 0; i < ndim-4; i++) start[i] = dims_offset[i];
  for (int j = 0; i < ndim; i++, j++) start[i] = subgrid[j] * coord[j];

  hsize_t *count = new hsize_t[ndim];
  H5Sget_simple_extent_dims(dataspace, count, NULL);
  for (i = 0; i < ndim-4; i++) count[i] = 1;
  for (int j = 0; i < ndim; i++, j++) count[i] = subgrid[j];
  H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start,NULL,count,NULL);

  hsize_t dims_mem = QDP::Layout::sitesOnNode();
  hid_t memspace = H5Screate_simple(1, &dims_mem, NULL);
  H5Sselect_elements(memspace, H5S_SELECT_SET, dims_mem, mapping);

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  H5Dwrite(dataset, type, memspace, dataspace, H5P_DEFAULT, data);

  H5Pclose(dxpl);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  delete [] start;
  delete [] count;
}
*/

void hdf5_writer::write_lattice_generic(const char *name, const int *dims_offset, hid_t type, const void *data)
{
  int i;
  hid_t dataset = H5Dopen(file, name, H5P_DEFAULT);
  hid_t dataspace = H5Dget_space(dataset);

  int ndim = H5Sget_simple_extent_ndims(dataspace);
  assert(ndim >= 4);

  QDP::multi1d<int> subgrid = QDP::Layout::subgridLattSize();
  QDP::multi1d<int> coord   = QDP::Layout::nodeCoord();

  hsize_t *start = new hsize_t[ndim];
  for (i = 0; i < ndim-4; i++) start[i] = dims_offset[i];
  for (int j = 0; i < ndim; i++, j++) start[i] = subgrid[j] * coord[j];

  hsize_t *count = new hsize_t[ndim];
  H5Sget_simple_extent_dims(dataspace, count, NULL);
  for (i = 0; i < ndim-4; i++) count[i] = 1;
  for (int j = 0; i < ndim; i++, j++) count[i] = subgrid[j];
  H5Sselect_hyperslab(dataspace,H5S_SELECT_SET,start,NULL,count,NULL);

  /*
  QDP::StopWatch swatch;
  swatch.reset();
  swatch.start();
  */
  hsize_t dims_mem = QDP::Layout::sitesOnNode();
  size_t tsize = H5Tget_size(type);
  char *buf = new char[tsize*dims_mem];
  for (int j = 0; j < dims_mem; j++)
    memcpy(buf + j*tsize, ((const char *) data) + mapping[j]*tsize, tsize);
  hid_t memspace = H5Screate_simple(1, &dims_mem, NULL);
  //H5Sselect_elements(memspace, H5S_SELECT_SET, dims_mem, mapping);
  /*
  swatch.stop();
  QDP::QDPIO::cout << "[hdf5_writer] remap data: "
		   << swatch.getTimeInSeconds() << " sec" << endl;
  */

  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);

  /*
  swatch.reset();
  swatch.start();
  */
  H5Dwrite(dataset, type, memspace, dataspace, dxpl, buf);
  /*
  swatch.stop();
  QDP::QDPIO::cout << "[hdf5_writer] write: "
		   << swatch.getTimeInSeconds() << " sec" << endl;
  */

  H5Pclose(dxpl);
  H5Sclose(memspace);
  H5Sclose(dataspace);
  H5Dclose(dataset);
  delete [] start;
  delete [] count;
  delete [] buf;
}


void hdf5_writer::write_lattice_data(const char *name, const int *dims_offset, const QDP::LatticeRealD &data)
{
  write_lattice_generic(name, dims_offset, H5T_NATIVE_DOUBLE, &data.elem(0));
}

void hdf5_writer::write_lattice_data(const char *name, int *dims_offset, int dim, const QDP::multi1d<QDP::LatticeRealD> &data)
{
  for (int i = 0; i < data.size(); i++) {
    dims_offset[dim] = i;
    write_lattice_generic(name, dims_offset, H5T_NATIVE_DOUBLE, &data[i].elem(0));
  }
}

void hdf5_writer::write_lattice_data(const char *name, int *dims_offset, int dim0, int dim1, const QDP::multi2d<QDP::LatticeRealD> &data)
{
  for (int i = 0; i < data.size2(); i++) {
    dims_offset[dim0] = i;
    for (int j = 0; j < data.size1(); j++) {
      dims_offset[dim1] = j;
      write_lattice_generic(name, dims_offset, H5T_NATIVE_DOUBLE, &data[i][j].elem(0));
    }
  }
}

void hdf5_writer::write_lattice_data(const char *name, int *dims_offset, int dim0, int dim1, int dim2, int dim3, const QDP::multi4d<QDP::LatticeRealD> &data)
{
  for (int i = 0; i < data.size4(); i++) {
    dims_offset[dim0] = i;
    for (int j = 0; j < data.size3(); j++) {
      dims_offset[dim1] = j;
      for(int k = 0; k < data.size2(); k++){
        dims_offset[dim2] = k;
        for(int l = 0 ; l < data.size1(); l++){
          dims_offset[dim3] = l;
          write_lattice_generic(name, dims_offset, H5T_NATIVE_DOUBLE, &data[i][j][k][l].elem(0));
        }
      }
    }
  }
}

void hdf5_writer::write_lattice_data(const char *name, int *dims_offset, int dim0, int dim1, int dim2, int dim3, const QDP::multi4d<QDP::LatticeComplexD> &data)
{
  for (int i = 0; i < data.size4(); i++) {
    dims_offset[dim0] = i;
    for (int j = 0; j < data.size3(); j++) {
      dims_offset[dim1] = j;
      for (int k = 0; k < data.size2(); k++) {
	dims_offset[dim2] = k;
	for (int l = 0; l < data.size1(); l++) {
	  dims_offset[dim3] = l;
	  write_lattice_generic(name, dims_offset, complex_t, &data[i][j][k][l].elem(0));
	}
      }
    }
  }
}

hdf5_writer_noflush::~hdf5_writer_noflush()
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  for (std::map<std::string,hid_t>::iterator it = datasets.begin();
       it != datasets.end(); it++)
    H5Dclose(it->second);
}

void hdf5_writer_noflush::write_generic_init(const char *name, int ndim,
					     const int *dims,
					     const char *const *axis_labels,
					     const char *const *axis_info, hid_t type)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;
  hid_t dset = write_generic_init_0(name, ndim, dims,
				    axis_labels, axis_info, type);
  datasets[std::string(name)] = dset;
}

void hdf5_writer_noflush::write_data_generic(const char *name, const int *dims_offset,
					     const int *dims_data, int ndim_data,
					     hid_t type, const void *data)
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;

  hid_t dataset = datasets[std::string(name)];
  write_data_generic_0(dataset, dims_offset, dims_data, ndim_data, type, data);
}

void hdf5_writer_noflush::flush()
{
  if (!use_mpi && QDP::Layout::nodeNumber() != 0)
    return;

  //H5Fflush(file, H5F_SCOPE_GLOBAL);
  for (std::map<std::string,hid_t>::iterator it = datasets.begin();
       it != datasets.end(); it++) {
    H5Dclose(it->second);
    it->second = H5Dopen(file, it->first.c_str(), H5P_DEFAULT);
  }
}
