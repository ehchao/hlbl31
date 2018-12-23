#include <qdp.h>

void read_light_by_light_infile(int argc, char *argv[]);

void connected_4pt_lccc(QDP::multi2d<QDP::LatticeRealD> &corr,
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::LatticeRealD &func1,
			const QDP::LatticeRealD &func2,
			int mu1, int mu2);
void connected_4pt_lccc(QDP::multi2d<QDP::LatticeComplexD> &corr,
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::LatticeComplexD &func1,
			const QDP::LatticeComplexD &func2,
			int mu1, int mu2);
void connected_4pt_lccc(QDP::multi4d<QDP::LatticeComplexD> &corr,
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::multi1d<QDP::LatticeComplexD> &func1,
			const QDP::multi1d<QDP::LatticeComplexD> &func2);
void connected_4pt_lccc(QDP::multi4d<QDP::LatticeComplexD> &corr,
			const QDP::multi1d<QDP::LatticeColorMatrixD> &u,
			const QDP::LatticePropagatorD &prop_fwd,
			const QDP::multi1d<QDP::LatticeComplexD> &func1);

void write_light_by_light_4pt(int config_no,
			      const QDP::multi2d<QDP::LatticeRealD> &corr);
void write_light_by_light_4pt(int config_no,
			      const QDP::multi4d<QDP::LatticeComplexD> &corr);
void write_light_by_light_FT_and_3dsums_init(hdf5_writer &writer,
				     const QDP::multi1d<int> &max_mom_comp,
				     const QDP::multi2d<int> &src_pos);
void write_light_by_light_FT_and_3dsums(hdf5_writer &writer,
		const QDP::multi4d<QDP::LatticeComplexD> &corr,
		const QDP::multi1d<int> &src_pos,
	        const QDP::multi1d<int> &max_mom_comp,
		int isrc);
