#include <qdp.h>
#include <vector>

void read_hlbl_g_minus_2_infile(int argc, char *argv[]);
void get_func_zrho_sine(QDP::LatticeRealD &func, int rho, int sigma,
  const QDP::multi1d<int> &pos_zero, const QDP::multi1d<int> &y_coord);
void get_func_zrho(QDP::LatticeRealD &func, int rho,
  const QDP::multi1d<int> &pos_zero, const QDP::multi1d<int> &y_coord);
void get_func_zrho_delta(QDP::LatticeRealD &func, int rho, int sigma,
  const QDP::multi1d<int> &pos_zero, const QDP::multi1d<int> &y_coord, const QDP::multi1d<int> &delta_coord);
template <class LatticeType>
void LatticeToRegion(LatticeType &func, int region, const QDP::multi1d<int> &pos_zero);
template <class T>
void llxx_seqprop(QDP::LatticePropagatorD &prop_seq,
               const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			   const QDP::LatticePropagatorD &prop_fwd,
			   const T &func,
			   int musnk,
               const QDP::multi1d<int> &srcpt,
               const QDP::multi1d<int> &pos_zero,
			   int LatticeRegion);
void connected_4pt_llxx(
            std::vector<QDP::LatticeRealD> &corr, /* corr[kernel] */
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::LatticePropagatorD &prop_fwd_y,
			const QDP::LatticePropagatorD &prop_seq,
			const QDP::LatticePropagatorD &prop_seq_y,
			const QDP::LatticeRealD &func_zrho,
			const std::vector<QDP::multi3d<QDP::LatticeRealD> > &func_kernels_mnl,
            const QDP::multi1d<int> &pos_zero,
			const QDP::multi1d<int> &y,
			int sigma,
			int LatticeRegion);
void connected_4pt_llxx(
            std::vector<QDP::LatticeRealD> &corr, /* corr[kernel] */
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::LatticePropagatorD &prop_fwd_y,
			const QDP::LatticePropagatorD &prop_seq,
			const QDP::LatticePropagatorD &prop_seq_y,
			const QDP::LatticeRealD &func_zrho,
			const std::vector<QDP::multi3d<QDP::LatticeRealD> > &func_kernels_mnl,
            const QDP::multi1d<int> &pos_zero,
			const QDP::multi1d<int> &y,
			int sigma,
			char Vx_lc,
			int LatticeRegion);
void connected_4pt_llxx(QDP::multi2d<QDP::LatticeRealD> &corr, /* corr[svt][0xyz_0xzy_0yxz] */
                        const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
                        const QDP::LatticePropagatorD &prop_fwd,
                        const QDP::LatticePropagatorD &prop_fwd_y,
                        const QDP::multi4d<QDP::LatticeRealD> &func_kernel_mnl, /* connected to local current */
						const QDP::multi1d<int> &pos_zero,
                        const QDP::multi1d<int> &y,
                        int rho, int sigma);
void connected_4pt_llxx(QDP::multi1d<QDP::LatticeRealD> &corr, /* corr[svt] */
                        const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
                        const QDP::LatticePropagatorD &prop_fwd,
                        const QDP::LatticePropagatorD &prop_fwd_y,
                        const QDP::multi4d<QDP::LatticeRealD> &func_kernel_mnl, /* connected to local current */
						const QDP::multi1d<int> &pos_zero,
                        const QDP::multi1d<int> &y,
                        int rho, int sigma);
void check_current_conservation_V_sigma_z(
            const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
            const QDP::LatticePropagatorD &prop_fwd,
            const QDP::LatticePropagatorD &prop_fwd_y,
            const QDP::multi1d<int> &pos_zero,
            const QDP::multi1d<int> &y);
void write_hlbl_g_minus_2_init(hdf5_writer &writer, QDP::multi2d<int> &src_pos);
void write_hlbl_g_minus_2_init_noboxtype(hdf5_writer &writer, QDP::multi2d<int> &src_pos);
void write_hlbl_g_minus_2(hdf5_writer &writer, const QDP::multi2d<QDP::LatticeRealD> &corr,
  const QDP::multi1d<int> &src_pos, int src_num);
void write_hlbl_g_minus_2(hdf5_writer &writer, const QDP::multi1d<QDP::LatticeRealD> &corr,
  const QDP::multi1d<int> &src_pos, int src_num);
