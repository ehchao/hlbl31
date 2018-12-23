#ifndef HDF5_OUTPUT_H
#define HDF5_OUTPUT_H
#include <hdf5.h>
#include <qdp.h>
class hdf5_writer
{
public:
  hdf5_writer(const char *filename, bool use_mpi=false);
  hdf5_writer(const char *prefix, int cfg_no, bool use_mpi=false);
  virtual ~hdf5_writer();

  void write_axis_info(const char *name, int count, int size, const int *data);
  void write_axis_info(const char *name, int count, int size,
                       const double *data);
  void write_axis_info(const char *name, int count, const char *const *data);

  hid_t write_generic_init_0(const char *name, int ndim, const int *dims,
			     const char *const *axis_labels,
			     const char *const *axis_info, hid_t type);
  virtual void write_generic_init(const char *name, int ndim, const int *dims,
				  const char *const *axis_labels,
				  const char *const *axis_info, hid_t type);

  void write_complex_init(const char *name, int ndim, const int *dims,
                          const char *const *axis_labels = NULL,
                          const char *const *axis_info = NULL);
  void write_real_init(const char *name, int ndim, const int *dims,
                       const char *const *axis_labels = NULL,
                       const char *const *axis_info = NULL);

  void write_data(const char *name, const int *dims_offset, int axis,
                  const QDP::multi1d < QDP::ComplexD > &data);

  /* must satisfy axis0 < axis1 */
  void write_data(const char *name, const int *dims_offset,
                  int axis0, int axis1,
                  const QDP::multi2d < QDP::ComplexD > &data);

  /* must satisfy axis0 < axis1 < axis2 */
  void write_data(const char *name, const int *dims_offset,
                  int axis0, int axis1, int axis2,
                  const QDP::multi3d < QDP::ComplexD > &data);

  /* must satisfy axis0 < axis1 < axis2 < axis3 */
  void write_data(const char *name, const int *dims_offset,
		  int axis0, int axis1, int axis2, int axis3,
		  const QDP::multi4d < QDP::ComplexD > &data);

  void write_data(const char *name, const int *dims_offset, int axis,
                  const QDP::multi1d < QDP::RealD > &data);

  /* must satisfy axis0 < axis1 */
  void write_data(const char *name, const int *dims_offset,
                  int axis0, int axis1,
                  const QDP::multi2d < QDP::RealD > &data);

  /* must satisfy axis0 < axis1 < axis2 */
  void write_data(const char *name, const int *dims_offset,
                  int axis0, int axis1, int axis2,
                  const QDP::multi3d < QDP::RealD > &data);
  void write_data(const char *name, const int *dims_offset, int dim0, int dim1, int dim2, int dim3, const QDP::multi4d<QDP::RealD> &data);

  void write_data_generic_0(hid_t dataset, const int *dims_offset,
			    const int *dims_data, int ndim_data,
			    hid_t type, const void *data);
  virtual void write_data_generic(const char *name, const int *dims_offset,
				  const int *dims_data, int ndim_data,
				  hid_t type, const void *data);

  void write_lattice_generic(const char *name, const int *dims_offset,
			     hid_t type, const void *data);

  void write_lattice_data_locally(const char *name, int *dims_offset, int size0, int size1, int size2, int size3, const QDP::multi4d<QDP::LatticeRealD> &data);

  void write_lattice_data_parallel(const char *file_name, const char *dset_name, const int *dims_offset, int size0, int size1, int size2, int size3, double ********ker_loc);

  void write_lattice_data(const char *name, const int *dims_offset,
			  const QDP::LatticeRealD &data);

  void write_lattice_data(const char *name, int *dims_offset, int axis,
			  const QDP::multi1d<QDP::LatticeRealD> &data);

  void write_lattice_data(const char *name, int *dims_offset,
			  int axis0, int axis1,
			  const QDP::multi2d<QDP::LatticeRealD> &data);

  void write_lattice_data(const char *name, int *dims_offset,
			  int axis0, int axis1, int axis2, int axis3,
			  const QDP::multi4d<QDP::LatticeRealD> &data);

  void write_lattice_data(const char *name, int *dims_offset,
			  int axis0, int axis1, int axis2, int axis3,
			  const QDP::multi4d<QDP::LatticeComplexD> &data);

  virtual void flush() {}

protected:
    hid_t file;
  hid_t complex_t;
  hid_t axes_group;
  hsize_t *mapping;
  bool use_mpi;

  void init(const char *filename, bool mpi);
};

class hdf5_writer_noflush : public hdf5_writer
{
public:
  hdf5_writer_noflush(const char *filename, bool use_mpi=false)
    : hdf5_writer(filename, use_mpi) {}
  hdf5_writer_noflush(const char *prefix, int cfg_no, bool use_mpi=false)
    : hdf5_writer(prefix, cfg_no, use_mpi) {}
  virtual ~hdf5_writer_noflush();

  virtual void write_generic_init(const char *name, int ndim, const int *dims,
				  const char *const *axis_labels,
				  const char *const *axis_info, hid_t type);
  virtual void write_data_generic(const char *name, const int *dims_offset,
				  const int *dims_data, int ndim_data,
				  hid_t type, const void *data);
  virtual void flush();  
  
private:
  std::map<std::string,hid_t> datasets;
};
#endif
