#include "observer.h"
#include <mpi.h>

namespace QDP
{



/* --------------------------------------------------------------------------*/
/**
* \brief  Read in cls config. This is legacy code (very slow!)
*
* \param u
* \param cfg_file
*/
/* ----------------------------------------------------------------------------*/
  void readCLS(multi1d < LatticeColorMatrixD > &u, const string & cfg_file)
  {
    u.resize(Nd);

    int i, j, k, l;
      multi1d < int >nrow(Nd);
      multi1d < int >site(Nd);
      nrow = Layout::lattSize();
    string nullstring;
      QDP::BinaryFileReader txtfile;
    ColorMatrixD matrix;
    int lat[Nd];
    bool swap_b = false;

      try
    {
      txtfile.open(cfg_file);

      // checks we are using the correct lattice size
      for (i = 0; i < Nd; i++)
      {
        txtfile.QDP::BinaryReader::read(lat[i]);
        /* QDP assumes BigEndian */
        /* CLS configs saved in little Endian */
        QDPUtil::byte_swap((void *) &lat[i], sizeof(int), 1);
        //QDPIO::cout << " Value is  " << lat[i] << endl;
      }
      if (lat[0] != nrow[3] || lat[1] != nrow[0] || lat[2] != nrow[1]
          || lat[3] != nrow[2])
      {
        QDPIO::cout << "Error: wrong lattice extensions." << endl;
        QDP_abort(1);
      }
      // iterates through all links
      double plaq;
      int old;
      if (!QDPUtil::big_endian())
      {
        swap_b = true;
      }
      // Skip Average Plaquette
      txtfile.QDP::BinaryReader::read(plaq);
      QDPUtil::byte_swap((void *) &plaq, sizeof(plaq), 1);
      for (i = 0; i < lat[0]; i++)
      {
        for (j = 0; j < lat[1]; j++)
        {
          for (k = 0; k < lat[2]; k++)
          {
            for (l = 0; l < lat[3]; l++)
            {
              if ((i + j + k + l) % 2 == 0)
                continue;
              site[0] = j;
              site[1] = k;
              site[2] = l;
              site[3] = i;
              // Read Color Matrix 
              // + 0

              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              pokeSite(u[3], matrix, site);
              // - 0
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              old = site[3];
              if (site[3] == 0)
              {
                site[3] = lat[0] - 1;
              }
              else
              {
                site[3] = old - 1;
              }
              pokeSite(u[3], matrix, site);
              site[3] = old;
              // Read Color Matrix 
              // + 1 
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              pokeSite(u[0], matrix, site);
              // - 1
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);

              old = site[0];
              if (site[0] == 0)
              {
                site[0] = lat[1] - 1;
              }
              else
              {
                site[0] = old - 1;
              }
              pokeSite(u[0], matrix, site);
              site[0] = old;

              // Read Color Matrix 
              // + 2 
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              pokeSite(u[1], matrix, site);
              // - 2
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              old = site[1];
              if (site[1] == 0)
              {
                site[1] = lat[2] - 1;
              }
              else
              {
                site[1] = old - 1;
              }
              pokeSite(u[1], matrix, site);
              site[1] = old;

              // Read Color Matrix 
              // + 3 
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              pokeSite(u[2], matrix, site);
              // - 3
              read(txtfile, matrix);
              if (swap_b)
                QDPUtil::byte_swap((void *) &matrix.elem(), sizeof(double),
                                   2 * Nc * Nc);
              old = site[2];
              if (site[2] == 0)
              {
                site[2] = lat[3] - 1;
              }
              else
              {
                site[2] = old - 1;
              }
              pokeSite(u[2], matrix, site);
              site[2] = old;

            }
          }
        }
      }

      QDPIO::cout << "Checks: ======= " << endl;
      ComplexD val = 0;
      for (int nu = 0; nu < 4; ++nu)
      {
        ComplexD tmp = sum(trace(u[nu] * adj(u[nu])));
        val += tmp;
      }
      QDPIO::cout << "Unitarity Check, should be ( 1 , 0 ):  " <<
        val / (lat[0] * lat[1] * lat[2] * lat[3]) / Nc / Nd << endl;
      QDP::Double w_plaq;
      w_plaq = 0;

      // Compute the average plaquettes
      for (int mu = 1; mu < Nd; ++mu)
      {
        for (int nu = 0; nu < mu; ++nu)
        {
          LatticeColorMatrix tmp_0 =
            shift(u[nu], FORWARD, mu) * adj(shift(u[mu], FORWARD, nu));
          LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);
          Double tmp = sum(real(trace(u[mu] * tmp_1)));
          w_plaq += tmp;
        }
      }
      QDPIO::cout << "Check Plaquette, should be 0: " << (w_plaq /
                                                          (6 *
                                                           (lat[0] * lat[1] *
                                                            lat[2] * lat[3])) -
                                                          plaq) << endl;

      txtfile.close();
    }                           // end try
    catch(const std::string & e)
    {
      QDPIO::cout << " caught exception " << e << endl;
      QDP_abort(1);
    }
  }



/* --------------------------------------------------------------------------*/
/**
* \brief  Copy single gauge link at site x form openQCD to QDP
*
* \param mat
* \param u
*/
/* ----------------------------------------------------------------------------*/
  static void copy_ufld_dble(multi1d < ColorMatrixD > &mat, su3_dble * u)
  {
    int i, k;
    ComplexD tmp;

    for (i = 0; i < 4; i++)
    {
      if (i == 0)
      {
        k = 3;
      }
      else
      {
        k = i - 1;
      }

      tmp = QDP::cmplx(Real(u[i].c11.re), Real(u[i].c11.im));
      pokeColor(mat[k], tmp, 0, 0);
      tmp = QDP::cmplx(Real(u[i].c12.re), Real(u[i].c12.im));
      pokeColor(mat[k], tmp, 0, 1);
      tmp = QDP::cmplx(Real(u[i].c13.re), Real(u[i].c13.im));
      pokeColor(mat[k], tmp, 0, 2);

      tmp = QDP::cmplx(Real(u[i].c21.re), Real(u[i].c21.im));
      pokeColor(mat[k], tmp, 1, 0);
      tmp = QDP::cmplx(Real(u[i].c22.re), Real(u[i].c22.im));
      pokeColor(mat[k], tmp, 1, 1);
      tmp = QDP::cmplx(Real(u[i].c23.re), Real(u[i].c23.im));
      pokeColor(mat[k], tmp, 1, 2);

      tmp = QDP::cmplx(Real(u[i].c31.re), Real(u[i].c31.im));
      pokeColor(mat[k], tmp, 2, 0);
      tmp = QDP::cmplx(Real(u[i].c32.re), Real(u[i].c32.im));
      pokeColor(mat[k], tmp, 2, 1);
      tmp = QDP::cmplx(Real(u[i].c33.re), Real(u[i].c33.im));
      pokeColor(mat[k], tmp, 2, 2);
    }

  }

  static void copy_color_dble(ColorMatrixD & mat, su3_dble * u)
  {
    ComplexD tmp;


    tmp = peekColor(mat, 0, 0);
    u->c11.re = toDouble(real(tmp));
    u->c11.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 0, 1);
    u->c12.re = toDouble(real(tmp));
    u->c12.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 0, 2);
    u->c13.re = toDouble(real(tmp));
    u->c13.im = toDouble(imag(tmp));

    tmp = peekColor(mat, 1, 0);
    u->c21.re = toDouble(real(tmp));
    u->c21.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 1, 1);
    u->c22.re = toDouble(real(tmp));
    u->c22.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 1, 2);
    u->c23.re = toDouble(real(tmp));
    u->c23.im = toDouble(imag(tmp));

    tmp = peekColor(mat, 2, 0);
    u->c31.re = toDouble(real(tmp));
    u->c31.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 2, 1);
    u->c32.re = toDouble(real(tmp));
    u->c32.im = toDouble(imag(tmp));
    tmp = peekColor(mat, 2, 2);
    u->c33.re = toDouble(real(tmp));
    u->c33.im = toDouble(imag(tmp));

  }

  static void copy_color_dble(multi1d < ColorMatrixD > &mat, su3_dble ** u)
  {
    int i, k;
    ComplexD tmp;

    for (i = 0; i < 4; i++)
    {
      if (i == 0)
      {
        k = 3;
      }
      else
      {
        k = i - 1;
      }

      tmp = peekColor(mat[k], 0, 0);
      u[i]->c11.re = toDouble(real(tmp));
      u[i]->c11.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 0, 1);
      u[i]->c12.re = toDouble(real(tmp));
      u[i]->c12.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 0, 2);
      u[i]->c13.re = toDouble(real(tmp));
      u[i]->c13.im = toDouble(imag(tmp));

      tmp = peekColor(mat[k], 1, 0);
      u[i]->c21.re = toDouble(real(tmp));
      u[i]->c21.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 1, 1);
      u[i]->c22.re = toDouble(real(tmp));
      u[i]->c22.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 1, 2);
      u[i]->c23.re = toDouble(real(tmp));
      u[i]->c23.im = toDouble(imag(tmp));

      tmp = peekColor(mat[k], 2, 0);
      u[i]->c31.re = toDouble(real(tmp));
      u[i]->c31.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 2, 1);
      u[i]->c32.re = toDouble(real(tmp));
      u[i]->c32.im = toDouble(imag(tmp));
      tmp = peekColor(mat[k], 2, 2);
      u[i]->c33.re = toDouble(real(tmp));
      u[i]->c33.im = toDouble(imag(tmp));
    }

  }


  static multi2d<int> qdp_to_openqcd_map; // [x,0] rank, [x,1] buffer-idx
  static multi1d<int> qdp_to_openqcd_counts;
  static multi1d<int> qdp_to_openqcd_dest;
  static int          qdp_to_openqcd_ndest;
  static multi2d<int> openqcd_to_qdp_map; // [x,0] rank, [x,1] buffer-idx
  static multi1d<int> openqcd_to_qdp_counts;
  static multi1d<int> openqcd_to_qdp_dest;
  static int          openqcd_to_qdp_ndest;
  static bool mappings_allocated = false;

  static void alloc_mappings()
  {
    if (mappings_allocated)
      return;
    /****************************************************************
     * In communications buffers, we will always store data sorted
     * according to the QDP site order, whether QDP is the source
     * or destination.
     ****************************************************************/
    int rank, nrank;
    MPI_Comm_size(MPI_COMM_WORLD, &nrank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int openqcd_sites = L0*L1*L2*L3;
    int nnode = Layout::numNodes();
    assert(nnode == nrank);
    int node  = Layout::nodeNumber();
    assert(node  == rank);
    int qdp_sites = Layout::sitesOnNode();
    qdp_to_openqcd_map.resize(qdp_sites,2);
    multi1d<bool> qdp_to_openqcd_xfer(nrank);
    qdp_to_openqcd_xfer = false;
    qdp_to_openqcd_ndest = 0;

    /* For QDP->OpenQCD communications:
     * Count how many local sites we have to send to each destination rank.
     * Also, store a mapping from the local site index into temporary
     * communications buffers.
     */
    for (int i = 0; i < qdp_sites; i++) {
      multi1d<int> qdp_coords = Layout::siteCoords(node, i);
      int x[4]; // openqcd coordinate, on the global lattice
      x[0] = qdp_coords[3];
      x[1] = qdp_coords[0];
      x[2] = qdp_coords[1];
      x[3] = qdp_coords[2];
      int openqcd_rank, openqcd_site;
      ipt_global(x,&openqcd_rank,&openqcd_site);
      if (!qdp_to_openqcd_xfer[openqcd_rank]) {
	qdp_to_openqcd_xfer[openqcd_rank] = true;
	qdp_to_openqcd_ndest++;
      }
    }
    qdp_to_openqcd_dest.resize(qdp_to_openqcd_ndest);
    multi1d<int> qdp_to_openqcd_destbuf(nrank);
    qdp_to_openqcd_destbuf = -1;
    for (int i = 0, j = 0; i < nrank; i++) {
      if (qdp_to_openqcd_xfer[i]) {
	qdp_to_openqcd_dest[j]    = i;
	qdp_to_openqcd_destbuf[i] = j;
	j++;
      }
    }
    qdp_to_openqcd_counts.resize(qdp_to_openqcd_ndest);
    qdp_to_openqcd_counts = 0;
    for (int i = 0; i < qdp_sites; i++) {
      multi1d<int> qdp_coords = Layout::siteCoords(node, i);
      int x[4]; // openqcd coordinate, on the global lattice
      x[0] = qdp_coords[3];
      x[1] = qdp_coords[0];
      x[2] = qdp_coords[1];
      x[3] = qdp_coords[2];
      int openqcd_rank, openqcd_site;
      ipt_global(x,&openqcd_rank,&openqcd_site);
      int destbuf = qdp_to_openqcd_destbuf[openqcd_rank];
      qdp_to_openqcd_map[i][0] = destbuf;
      qdp_to_openqcd_map[i][1] = qdp_to_openqcd_counts[destbuf]++;
    }


    multi1d<bool> openqcd_to_qdp_xfer(nrank);
    openqcd_to_qdp_xfer = false;
    openqcd_to_qdp_ndest = 0;
    openqcd_to_qdp_map.resize(openqcd_sites,2);
    /* Same as above, but for OpenQCD->QDP communications
     */
    for (int x0 = 0; x0 < L0; x0++)
    for (int x1 = 0; x1 < L1; x1++)
    for (int x2 = 0; x2 < L2; x2++)
    for (int x3 = 0; x3 < L3; x3++) {
      multi1d<int> qdp_coords(4);
      qdp_coords[0] = cpr[1]*L1 + x1;
      qdp_coords[1] = cpr[2]*L2 + x2;
      qdp_coords[2] = cpr[3]*L3 + x3;
      qdp_coords[3] = cpr[0]*L0 + x0;
      int qdp_node = Layout::nodeNumber(qdp_coords);
      if (!openqcd_to_qdp_xfer[qdp_node]) {
	openqcd_to_qdp_xfer[qdp_node] = true;
	openqcd_to_qdp_ndest++;
      }
    }
    openqcd_to_qdp_dest.resize(openqcd_to_qdp_ndest);
    multi1d<int> openqcd_to_qdp_destbuf(nnode);
    openqcd_to_qdp_destbuf = -1;
    for (int i = 0, j = 0; i < nnode; i++) {
      if (openqcd_to_qdp_xfer[i]) {
	openqcd_to_qdp_dest[j]    = i;
	openqcd_to_qdp_destbuf[i] = j;
	j++;
      }
    }
    openqcd_to_qdp_counts.resize(openqcd_to_qdp_ndest);
    openqcd_to_qdp_counts = 0;
    for (int x0 = 0; x0 < L0; x0++)
    for (int x1 = 0; x1 < L1; x1++)
    for (int x2 = 0; x2 < L2; x2++)
    for (int x3 = 0; x3 < L3; x3++) {
      int i = ipt[x3 + L3*x2 + L2*L3*x1 + L1*L2*L3*x0];
      multi1d<int> qdp_coords(4);
      qdp_coords[0] = cpr[1]*L1 + x1;
      qdp_coords[1] = cpr[2]*L2 + x2;
      qdp_coords[2] = cpr[3]*L3 + x3;
      qdp_coords[3] = cpr[0]*L0 + x0;
      int qdp_node = Layout::nodeNumber(qdp_coords);
      int destbuf = openqcd_to_qdp_destbuf[qdp_node];
      openqcd_to_qdp_map[i][0] = destbuf;
      /* This doesn't necessarily match the QDP site order.
       openqcd_to_qdp_map[i][1] = openqcd_to_qdp_counts[destbuf];
       */
      openqcd_to_qdp_counts[destbuf]++;
    }
    /* For each node, we need to re-sort the openqcd_to_qdp_map[i][1].
     * Rather than actually performing a sort, we can get the information
     * by looping over qdp sites again, and communicating it.
     */
    multi1d< multi1d<int> > send_tmp(qdp_to_openqcd_ndest);
    multi1d< multi1d<int> > recv_tmp(openqcd_to_qdp_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      send_tmp[i].resize(qdp_to_openqcd_counts[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recv_tmp[i].resize(openqcd_to_qdp_counts[i]);
    multi1d<int> count_tmp(qdp_to_openqcd_ndest);
    count_tmp = 0;
    for (int i = 0; i < qdp_sites; i++) {
      multi1d<int> qdp_coords = Layout::siteCoords(node, i);
      int x[4]; // openqcd coordinate, on the global lattice
      x[0] = qdp_coords[3];
      x[1] = qdp_coords[0];
      x[2] = qdp_coords[1];
      x[3] = qdp_coords[2];
      int openqcd_rank, openqcd_site;
      ipt_global(x,&openqcd_rank,&openqcd_site);
      int destbuf = qdp_to_openqcd_destbuf[openqcd_rank];
      send_tmp[destbuf][count_tmp[destbuf]++] = openqcd_site;
    }
    multi1d<MPI_Request> sends(qdp_to_openqcd_ndest);
    multi1d<MPI_Request> recvs(openqcd_to_qdp_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Isend(&send_tmp[i][0], qdp_to_openqcd_counts[i], MPI_INT,
		qdp_to_openqcd_dest[i], ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Irecv(&recv_tmp[i][0], openqcd_to_qdp_counts[i], MPI_INT,
		openqcd_to_qdp_dest[i], ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    //MPI_Waitall(openqcd_to_qdp_ndest, &recvs[0], &statuses[0]);
    //for (int i = 0; i < openqcd_to_qdp_ndest; i++) {
    //  MPI_Wait(&recvs[i], &statuses[i]);
    int i;
    while (MPI_Waitany(openqcd_to_qdp_ndest, &recvs[0], &i, &statuses[0]),
	   i != MPI_UNDEFINED) {
      for (int j = 0; j < openqcd_to_qdp_counts[i]; j++) {
	int site = recv_tmp[i][j];
	assert(openqcd_to_qdp_map[site][0] == i);
	openqcd_to_qdp_map[site][1] = j;
      }
    }

    MPI_Waitall(qdp_to_openqcd_ndest, &sends[0], &statuses[0]);
    mappings_allocated = true;
  }


/* --------------------------------------------------------------------------*/
/**
* \brief  Copies gauge links from QDP to openQCD.
*
* \param u QDO gauge links
* \param ud openQCD gauge links
*/
/* ----------------------------------------------------------------------------*/
  void copy_to_openQCD(const multi1d < LatticeColorMatrixD > &u, su3_dble * ud)
  {
    start_time(copy_to_openQCD_gauge_links);
    alloc_mappings();

    int ofs[4];
    ofs[0] = VOLUME + (FACE0 / 2);
    ofs[1] = ofs[0] + (FACE0 / 2) + (FACE1 / 2);
    ofs[2] = ofs[1] + (FACE1 / 2) + (FACE2 / 2);
    ofs[3] = ofs[2] + (FACE2 / 2) + (FACE3 / 2);

    int snu[4];
    snu[0] = 0;
    snu[1] = snu[0] + (FACE0 / 2);
    snu[2] = snu[1] + (FACE1 / 2);
    snu[3] = snu[2] + (FACE2 / 2);

    multi1d<multi1d<ColorMatrixD> > uu(4);
    for (int m = 0; m < 4; m++) {
      uu[m].resize(QDP::Layout::sitesOnNode());
      QDP_extract(uu[m], u[m], all);
    }

    multi1d<multi2d<su3_dble> > sendbuf(qdp_to_openqcd_ndest);
    multi1d<multi2d<su3_dble> > recvbuf(openqcd_to_qdp_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      sendbuf[i].resize(qdp_to_openqcd_counts[i],4);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recvbuf[i].resize(openqcd_to_qdp_counts[i],4);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      int dest_buf = qdp_to_openqcd_map[i][0];
      int buf_idx  = qdp_to_openqcd_map[i][1];
      su3_dble *udb[4];
      multi1d<ColorMatrixD> mat(4);
      for (int m = 0; m < 4; m++) {
	udb[m] = &sendbuf[dest_buf][buf_idx][m];
	mat[m] = uu[m][i];
      }
      copy_color_dble(mat, udb);
    }
    multi1d<MPI_Request> sends(qdp_to_openqcd_ndest);
    multi1d<MPI_Request> recvs(openqcd_to_qdp_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Isend(&sendbuf[i][0][0],
		qdp_to_openqcd_counts[i]*4*sizeof(su3_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Irecv(&recvbuf[i][0][0],
		openqcd_to_qdp_counts[i]*4*sizeof(su3_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(openqcd_to_qdp_ndest, &recvs[0], &statuses[0]);
    for (int x0 = 0; x0 < L0; x0++)
    for (int x1 = 0; x1 < L1; x1++)
    for (int x2 = 0; x2 < L2; x2++)
    for (int x3 = 0; x3 < L3; x3++) {
      int i = ipt[x3 + L3*x2 + L2*L3*x1 + L1*L2*L3*x0];
      int source_buf = openqcd_to_qdp_map[i][0];
      int buf_idx    = openqcd_to_qdp_map[i][1];
      if (((x0 + x1 + x2 + x3) & 0x1) == 1) {
	assert(i < VOLUME);
	for (int m = 0; m < 4; m++)
	  ud[8*(i-VOLUME/2) + 2*m] = recvbuf[source_buf][buf_idx][m];
      }	else {
	for (int m = 0; m < 4; m++) {
	  int iu = iup[i][m];
	  if (iu >= VOLUME)
	    iu = 4*VOLUME + snu[m] + iu - ofs[m] - BNDRY/2;
	  else
	    iu = 8*(iu-VOLUME/2) + 2*m+1;
	  ud[iu] = recvbuf[source_buf][buf_idx][m];
	}
      }
    }
    /* We still need to retrieve the backward links on the lower faces.
     * These are currently stored as boundary links on adjacent ranks.
     */
    multi1d<multi1d<su3_dble> > bdrybuf(4);
    bdrybuf[0].resize(FACE0/2);
    bdrybuf[1].resize(FACE1/2);
    bdrybuf[2].resize(FACE2/2);
    bdrybuf[3].resize(FACE3/2);
    multi1d<MPI_Request> bsends(4), brecvs(4);
    ctag = mpi_tag();
    MPI_Isend(ud + 4*VOLUME + snu[0], FACE0/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*0+1], ctag, MPI_COMM_WORLD, &bsends[0]);
    MPI_Irecv(&bdrybuf[0][0],         FACE0/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*0  ], ctag, MPI_COMM_WORLD, &brecvs[0]);
    MPI_Isend(ud + 4*VOLUME + snu[1], FACE1/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*1+1], ctag, MPI_COMM_WORLD, &bsends[1]);
    MPI_Irecv(&bdrybuf[1][0],         FACE1/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*1  ], ctag, MPI_COMM_WORLD, &brecvs[1]);
    MPI_Isend(ud + 4*VOLUME + snu[2], FACE2/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*2+1], ctag, MPI_COMM_WORLD, &bsends[2]);
    MPI_Irecv(&bdrybuf[2][0],         FACE2/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*2  ], ctag, MPI_COMM_WORLD, &brecvs[2]);
    MPI_Isend(ud + 4*VOLUME + snu[3], FACE3/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*3+1], ctag, MPI_COMM_WORLD, &bsends[3]);
    MPI_Irecv(&bdrybuf[3][0],         FACE3/2*sizeof(su3_dble)/sizeof(double),
	      MPI_DOUBLE, npr[2*3  ], ctag, MPI_COMM_WORLD, &brecvs[3]);
    multi1d<MPI_Status> bstatuses(4);
    MPI_Waitall(4, &brecvs[0], &bstatuses[0]);
    // we could instead loop over the boundary explicitly
    for (int x0 = 0; x0 < L0; x0++)
    for (int x1 = 0; x1 < L1; x1++)
    for (int x2 = 0; x2 < L2; x2++)
    for (int x3 = 0; x3 < L3; x3++) {
      if (((x0 + x1 + x2 + x3) & 0x1) == 1) {
	int i = ipt[x3 + L3*x2 + L2*L3*x1 + L1*L2*L3*x0];
	for (int m = 0; m < 4; m++) {
	  int id = idn[i][m];
	  if (id >= VOLUME) {
	    /* We have to fully-qualify the global "map" from OpenQCD, since
	     * 'using namespace QDP' imports the STL map.
	     */
	    int io = iup[::map[id-VOLUME]][m] - ofs[m] - BNDRY/2;
	    ud[8*(i-VOLUME/2) + 2*m+1] = bdrybuf[m][io];
	  }
	}
      }
    }      

    MPI_Waitall(qdp_to_openqcd_ndest, &sends[0], &statuses[0]);
    MPI_Waitall(4, &bsends[0], &bstatuses[0]);

    set_flags(UPDATED_UD);
    // update the remaining boundary links
    copy_bnd_ud();
    stop_time(copy_to_openQCD_gauge_links);
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copies gauge links from openQCD to QDP.
*
* \param u: Destination gauge links
* \param ud: Source gauge links
*/
/* ----------------------------------------------------------------------------*/
  void copy_from_openQCD(multi1d < LatticeColorMatrixD > &u, const su3_dble * ud)
  {
    /*******************************************************************
     * Assumes that the "boundary" links in ud have already been filled!
     *******************************************************************/
    start_time(copy_from_openQCD_gauge_links); 
    alloc_mappings();

    int ofs[4];
    ofs[0] = VOLUME + (FACE0 / 2);
    ofs[1] = ofs[0] + (FACE0 / 2) + (FACE1 / 2);
    ofs[2] = ofs[1] + (FACE1 / 2) + (FACE2 / 2);
    ofs[3] = ofs[2] + (FACE2 / 2) + (FACE3 / 2);

    int snu[4];
    snu[0] = 0;
    snu[1] = snu[0] + (FACE0 / 2);
    snu[2] = snu[1] + (FACE1 / 2);
    snu[3] = snu[2] + (FACE2 / 2);

    multi2d<ColorMatrixD> uu(4,QDP::Layout::sitesOnNode());

    multi1d<multi2d<su3_dble> > sendbuf(openqcd_to_qdp_ndest);
    multi1d<multi2d<su3_dble> > recvbuf(qdp_to_openqcd_ndest);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      sendbuf[i].resize(openqcd_to_qdp_counts[i],4);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      recvbuf[i].resize(qdp_to_openqcd_counts[i],4);
    for (int x0 = 0; x0 < L0; x0++)
    for (int x1 = 0; x1 < L1; x1++)
    for (int x2 = 0; x2 < L2; x2++)
    for (int x3 = 0; x3 < L3; x3++) {
      int i = ipt[x3 + L3*x2 + L2*L3*x1 + L1*L2*L3*x0];
      int dest_buf = openqcd_to_qdp_map[i][0];
      int buf_idx  = openqcd_to_qdp_map[i][1];
      if (((x0 + x1 + x2 + x3) & 0x1) == 1) {
	assert(i < VOLUME);
	for (int m = 0; m < 4; m++)
	  sendbuf[dest_buf][buf_idx][m] = ud[8*(i-VOLUME/2) + 2*m];
      }	else {
	for (int m = 0; m < 4; m++) {
	  int iu = iup[i][m];
	  if (iu >= VOLUME)
	    iu = 4*VOLUME + snu[m] + iu - ofs[m] - BNDRY/2;
	  else
	    iu = 8*(iu-VOLUME/2) + 2*m+1;
	  sendbuf[dest_buf][buf_idx][m] = ud[iu];
	}
      }
    }
    multi1d<MPI_Request> sends(openqcd_to_qdp_ndest);
    multi1d<MPI_Request> recvs(qdp_to_openqcd_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Isend(&sendbuf[i][0][0],
		openqcd_to_qdp_counts[i]*4*sizeof(su3_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Irecv(&recvbuf[i][0][0],
		qdp_to_openqcd_counts[i]*4*sizeof(su3_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(qdp_to_openqcd_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      int source_buf = qdp_to_openqcd_map[i][0];
      int buf_idx    = qdp_to_openqcd_map[i][1];
      multi1d<ColorMatrixD> tmpCol(4);
      copy_ufld_dble(tmpCol, &recvbuf[source_buf][buf_idx][0]);
      for (int m = 0; m < 4; m++)
	uu[m][i] = tmpCol[m];
    }
    for (int m = 0; m < 4; m++)
      QDP_insert(u[m], uu[m], all);

    MPI_Waitall(openqcd_to_qdp_ndest, &sends[0], &statuses[0]);
    stop_time(copy_from_openQCD_gauge_links);
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copy propagator at site x. Time component is mapped 0 <-> 3
*
* \param mat Destination propagator
* \param psi1 Source spinor_dble array
*/
/* ----------------------------------------------------------------------------*/
  static void copy_full_spinor_dble(PropagatorD & mat, spinor_dble * psi1)
  {
    int i, k;
    ComplexD tmp;
    ColorMatrixD cmat;
    for (i = 0; i < 4; i++)
    {
      cmat = 0;
      for (k = 0; k < 3; k++)
      {
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c1.c1.re),
                     Real(psi1[k + i * 3].c1.c1.im));
        pokeColor(cmat, tmp, 0, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c1.c2.re),
                     Real(psi1[k + i * 3].c1.c2.im));
        pokeColor(cmat, tmp, 1, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c1.c3.re),
                     Real(psi1[k + i * 3].c1.c3.im));
        pokeColor(cmat, tmp, 2, k);
      }
      pokeSpin(mat, cmat, 0, i);

      cmat = 0;
      for (k = 0; k < 3; k++)
      {
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c2.c1.re),
                     Real(psi1[k + i * 3].c2.c1.im));
        pokeColor(cmat, tmp, 0, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c2.c2.re),
                     Real(psi1[k + i * 3].c2.c2.im));
        pokeColor(cmat, tmp, 1, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c2.c3.re),
                     Real(psi1[k + i * 3].c2.c3.im));
        pokeColor(cmat, tmp, 2, k);
      }
      pokeSpin(mat, cmat, 1, i);



      cmat = 0;
      for (k = 0; k < 3; k++)
      {
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c3.c1.re),
                     Real(psi1[k + i * 3].c3.c1.im));
        pokeColor(cmat, tmp, 0, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c3.c2.re),
                     Real(psi1[k + i * 3].c3.c2.im));
        pokeColor(cmat, tmp, 1, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c3.c3.re),
                     Real(psi1[k + i * 3].c3.c3.im));
        pokeColor(cmat, tmp, 2, k);
      }
      pokeSpin(mat, cmat, 2, i);



      cmat = 0;
      for (k = 0; k < 3; k++)
      {
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c4.c1.re),
                     Real(psi1[k + i * 3].c4.c1.im));
        pokeColor(cmat, tmp, 0, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c4.c2.re),
                     Real(psi1[k + i * 3].c4.c2.im));
        pokeColor(cmat, tmp, 1, k);
        tmp =
          QDP::cmplx(Real(psi1[k + i * 3].c4.c3.re),
                     Real(psi1[k + i * 3].c4.c3.im));
        pokeColor(cmat, tmp, 2, k);
      }
      pokeSpin(mat, cmat, 3, i);
    }
  }

  static void copy_propagator(PropagatorD & mat, spinor_dble * psi)
  {
    int i, k;
    ComplexD tmp;
    ColorMatrixD cmat;
    for (i = 0; i < 4; i++)
    {
      cmat = peekSpin(mat, 0, i);
      for (k = 0; k < 3; k++)
      {
        tmp = peekColor(cmat, 0, k);
        psi[k + i * 3].c1.c1.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c1.c1.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 1, k);
        psi[k + i * 3].c1.c2.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c1.c2.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 2, k);
        psi[k + i * 3].c1.c3.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c1.c3.im = QDP::toDouble(imag(tmp));
      }

      cmat = peekSpin(mat, 1, i);
      for (k = 0; k < 3; k++)
      {
        tmp = peekColor(cmat, 0, k);
        psi[k + i * 3].c2.c1.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c2.c1.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 1, k);
        psi[k + i * 3].c2.c2.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c2.c2.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 2, k);
        psi[k + i * 3].c2.c3.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c2.c3.im = QDP::toDouble(imag(tmp));
      }

      cmat = peekSpin(mat, 2, i);
      for (k = 0; k < 3; k++)
      {
        tmp = peekColor(cmat, 0, k);
        psi[k + i * 3].c3.c1.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c3.c1.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 1, k);
        psi[k + i * 3].c3.c2.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c3.c2.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 2, k);
        psi[k + i * 3].c3.c3.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c3.c3.im = QDP::toDouble(imag(tmp));
      }

      cmat = peekSpin(mat, 3, i);
      for (k = 0; k < 3; k++)
      {
        tmp = peekColor(cmat, 0, k);
        psi[k + i * 3].c4.c1.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c4.c1.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 1, k);
        psi[k + i * 3].c4.c2.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c4.c2.im = QDP::toDouble(imag(tmp));
        tmp = peekColor(cmat, 2, k);
        psi[k + i * 3].c4.c3.re = QDP::toDouble(real(tmp));
        psi[k + i * 3].c4.c3.im = QDP::toDouble(imag(tmp));
      }
    }
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copies propagator from QDP to openQCD
*
* \param u LatticePropagator
* \param psi Destination pointer to spinor_dble array
*/
/* ----------------------------------------------------------------------------*/
  void copy_to_openQCD(const LatticePropagatorD & u, spinor_dble ** psi)
  {
    alloc_mappings();

    multi1d < PropagatorD > props(QDP::Layout::sitesOnNode());
    QDP_extract(props, u, all);

    multi1d<multi2d<spinor_dble> > sendbuf(qdp_to_openqcd_ndest);
    multi1d<multi2d<spinor_dble> > recvbuf(openqcd_to_qdp_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      sendbuf[i].resize(qdp_to_openqcd_counts[i],12);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recvbuf[i].resize(openqcd_to_qdp_counts[i],12);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      PropagatorD tmpprop = Gamma(2) * props[i] * Gamma(2);
      int dest_buf = qdp_to_openqcd_map[i][0];
      int buf_idx  = qdp_to_openqcd_map[i][1];
      copy_propagator(tmpprop, &sendbuf[dest_buf][buf_idx][0]);
    }
    multi1d<MPI_Request> sends(qdp_to_openqcd_ndest);
    multi1d<MPI_Request> recvs(openqcd_to_qdp_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Isend(&sendbuf[i][0][0],
		qdp_to_openqcd_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Irecv(&recvbuf[i][0][0],
		openqcd_to_qdp_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(openqcd_to_qdp_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < L0*L1*L2*L3; i++) {
      int source_buf = openqcd_to_qdp_map[i][0];
      int buf_idx    = openqcd_to_qdp_map[i][1];
      for (int j = 0; j < 12; j++)
	psi[j][i] = recvbuf[source_buf][buf_idx][j];
    }

    MPI_Waitall(qdp_to_openqcd_ndest, &sends[0], &statuses[0]);
      
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copies lattice propagator, i.e. an array of 12 lattice fermion fields
 from openQCD to QDP. This is used to copy the solve from openQCD.
*
* \param u: Destination LatticePropagator
* \param psi: Pointer to the solution array of spinor_dble[12]
*/
/* ----------------------------------------------------------------------------*/
  void copy_from_openQCD(LatticePropagatorD & u, const spinor_dble *const * psi)
  {
    alloc_mappings();

    multi1d < PropagatorD > props(QDP::Layout::sitesOnNode());

    multi1d<multi2d<spinor_dble> > sendbuf(openqcd_to_qdp_ndest);
    multi1d<multi2d<spinor_dble> > recvbuf(qdp_to_openqcd_ndest);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      sendbuf[i].resize(openqcd_to_qdp_counts[i],12);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      recvbuf[i].resize(qdp_to_openqcd_counts[i],12);
    for (int i = 0; i < L0*L1*L2*L3; i++) {
      int dest_buf = openqcd_to_qdp_map[i][0];
      int buf_idx  = openqcd_to_qdp_map[i][1];
      for (int j = 0; j < 12; j++)
	sendbuf[dest_buf][buf_idx][j] = psi[j][i];
    }
    multi1d<MPI_Request> sends(openqcd_to_qdp_ndest);
    multi1d<MPI_Request> recvs(qdp_to_openqcd_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Isend(&sendbuf[i][0][0],
		openqcd_to_qdp_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Irecv(&recvbuf[i][0][0],
		qdp_to_openqcd_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(qdp_to_openqcd_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      int source_buf = qdp_to_openqcd_map[i][0];
      int buf_idx    = qdp_to_openqcd_map[i][1];
      PropagatorD tmpprop;
      copy_full_spinor_dble(tmpprop, &recvbuf[source_buf][buf_idx][0]);
      props[i] = Gamma(2) * tmpprop * Gamma(2);
    }
    QDP_insert(u, props, all);

    MPI_Waitall(openqcd_to_qdp_ndest, &sends[0], &statuses[0]);
      
  }


  static void copy_spinor_dble(FermionD & mat, spinor_dble psi1)
  {
    ComplexD tmp;
    ColorVectorD cmat;



    tmp = QDP::cmplx(Real(psi1.c1.c1.re), Real(psi1.c1.c1.im));
    pokeColor(cmat, tmp, 0);
    pokeSpin(mat, cmat, 0);
    tmp = QDP::cmplx(Real(psi1.c1.c2.re), Real(psi1.c1.c2.im));
    pokeColor(cmat, tmp, 1);
    pokeSpin(mat, cmat, 0);
    tmp = QDP::cmplx(Real(psi1.c1.c3.re), Real(psi1.c1.c3.im));
    pokeColor(cmat, tmp, 2);
    pokeSpin(mat, cmat, 0);

    tmp = QDP::cmplx(Real(psi1.c2.c1.re), Real(psi1.c2.c1.im));
    pokeColor(cmat, tmp, 0);
    pokeSpin(mat, cmat, 1);
    tmp = QDP::cmplx(Real(psi1.c2.c2.re), Real(psi1.c2.c2.im));
    pokeColor(cmat, tmp, 1);
    pokeSpin(mat, cmat, 1);
    tmp = QDP::cmplx(Real(psi1.c2.c3.re), Real(psi1.c2.c3.im));
    pokeColor(cmat, tmp, 2);
    pokeSpin(mat, cmat, 1);

    tmp = QDP::cmplx(Real(psi1.c3.c1.re), Real(psi1.c3.c1.im));
    pokeColor(cmat, tmp, 0);
    pokeSpin(mat, cmat, 2);
    tmp = QDP::cmplx(Real(psi1.c3.c2.re), Real(psi1.c3.c2.im));
    pokeColor(cmat, tmp, 1);
    pokeSpin(mat, cmat, 2);
    tmp = QDP::cmplx(Real(psi1.c3.c3.re), Real(psi1.c3.c3.im));
    pokeColor(cmat, tmp, 2);
    pokeSpin(mat, cmat, 2);

    tmp = QDP::cmplx(Real(psi1.c4.c1.re), Real(psi1.c4.c1.im));
    pokeColor(cmat, tmp, 0);
    pokeSpin(mat, cmat, 3);
    tmp = QDP::cmplx(Real(psi1.c4.c2.re), Real(psi1.c4.c2.im));
    pokeColor(cmat, tmp, 1);
    pokeSpin(mat, cmat, 3);
    tmp = QDP::cmplx(Real(psi1.c4.c3.re), Real(psi1.c4.c3.im));
    pokeColor(cmat, tmp, 2);
    pokeSpin(mat, cmat, 3);
  }

  static void copy_fermion(FermionD & mat, spinor_dble * psi1)
  {
    ComplexD tmp;
    RealD res;
    ColorVectorD cmat;

    cmat = peekSpin(mat, 0);
    tmp = peekColor(cmat, 0);
    psi1->c1.c1.re = toDouble(real(tmp));
    psi1->c1.c1.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 1);
    psi1->c1.c2.re = toDouble(real(tmp));
    psi1->c1.c2.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 2);
    psi1->c1.c3.re = toDouble(real(tmp));
    psi1->c1.c3.im = toDouble(imag(tmp));

    cmat = peekSpin(mat, 1);
    tmp = peekColor(cmat, 0);
    psi1->c2.c1.re = toDouble(real(tmp));
    psi1->c2.c1.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 1);
    psi1->c2.c2.re = toDouble(real(tmp));
    psi1->c2.c2.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 2);
    psi1->c2.c3.re = toDouble(real(tmp));
    psi1->c2.c3.im = toDouble(imag(tmp));

    cmat = peekSpin(mat, 2);
    tmp = peekColor(cmat, 0);
    psi1->c3.c1.re = toDouble(real(tmp));
    psi1->c3.c1.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 1);
    psi1->c3.c2.re = toDouble(real(tmp));
    psi1->c3.c2.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 2);
    psi1->c3.c3.re = toDouble(real(tmp));
    psi1->c3.c3.im = toDouble(imag(tmp));

    cmat = peekSpin(mat, 3);
    tmp = peekColor(cmat, 0);
    psi1->c4.c1.re = toDouble(real(tmp));
    psi1->c4.c1.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 1);
    psi1->c4.c2.re = toDouble(real(tmp));
    psi1->c4.c2.im = toDouble(imag(tmp));
    tmp = peekColor(cmat, 2);
    psi1->c4.c3.re = toDouble(real(tmp));
    psi1->c4.c3.im = toDouble(imag(tmp));
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copies a lattice fermion field array of size 12 from QDP to openQCD. 
* This is  used for example in copying source fields to the solver of openQCD
*
* \param u: Fermion field array, e.g. source field
* \param psi: Destimation spinor has to be allocated using reserve_wsd(12);
*/
/* ----------------------------------------------------------------------------*/
  void copy_to_openQCD(const multi1d < LatticeFermionD > u, spinor_dble ** psi)
  {
    alloc_mappings();

    multi1d<multi1d<FermionD> > ferms(12);
    for (int i = 0; i < 12; i++) {
      ferms[i].resize(QDP::Layout::sitesOnNode());
      QDP_extract(ferms[i], u[i], all);
    }

    multi1d<multi2d<spinor_dble> > sendbuf(qdp_to_openqcd_ndest);
    multi1d<multi2d<spinor_dble> > recvbuf(openqcd_to_qdp_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      sendbuf[i].resize(qdp_to_openqcd_counts[i],12);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recvbuf[i].resize(openqcd_to_qdp_counts[i],12);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      int dest_buf = qdp_to_openqcd_map[i][0];
      int buf_idx  = qdp_to_openqcd_map[i][1];
      for (int j = 0; j < 12; j++) {
	FermionD tmpferm = Gamma(2) * ferms[j][i];
	copy_fermion(tmpferm, &sendbuf[dest_buf][buf_idx][j]);
      }
    }
    multi1d<MPI_Request> sends(qdp_to_openqcd_ndest);
    multi1d<MPI_Request> recvs(openqcd_to_qdp_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Isend(&sendbuf[i][0][0],
		qdp_to_openqcd_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Irecv(&recvbuf[i][0][0],
		openqcd_to_qdp_counts[i]*12*sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(openqcd_to_qdp_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < L0*L1*L2*L3; i++) {
      int source_buf = openqcd_to_qdp_map[i][0];
      int buf_idx    = openqcd_to_qdp_map[i][1];
      for (int j = 0; j < 12; j++)
	psi[j][i] = recvbuf[source_buf][buf_idx][j];
    }

    MPI_Waitall(qdp_to_openqcd_ndest, &sends[0], &statuses[0]);
      
  }

/* --------------------------------------------------------------------------*/
/**
* \brief  Copies a lattice fermion field from QDP++ to openQCD.
*
* \param u: Fermion lattice field to be copied
* \param psi: Destination spinor_dble, has to be allocated using reserve_wsd(1);
*/
/* ----------------------------------------------------------------------------*/
  void copy_to_openQCD(const LatticeFermionD & u, spinor_dble * psi)
  {
    alloc_mappings();

    multi1d < FermionD > ferms(QDP::Layout::sitesOnNode());
    QDP_extract(ferms, u, all);

    multi1d<multi1d<spinor_dble> > sendbuf(qdp_to_openqcd_ndest);
    multi1d<multi1d<spinor_dble> > recvbuf(openqcd_to_qdp_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      sendbuf[i].resize(qdp_to_openqcd_counts[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recvbuf[i].resize(openqcd_to_qdp_counts[i]);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      FermionD tmpferm = Gamma(2) * ferms[i];
      int dest_buf = qdp_to_openqcd_map[i][0];
      int buf_idx  = qdp_to_openqcd_map[i][1];
      copy_fermion(tmpferm, &sendbuf[dest_buf][buf_idx]);
    }
    multi1d<MPI_Request> sends(qdp_to_openqcd_ndest);
    multi1d<MPI_Request> recvs(openqcd_to_qdp_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Isend(&sendbuf[i][0],
		qdp_to_openqcd_counts[i] * sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Irecv(&recvbuf[i][0],
		openqcd_to_qdp_counts[i] * sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(openqcd_to_qdp_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < L0*L1*L2*L3; i++) {
      int source_buf = openqcd_to_qdp_map[i][0];
      int buf_idx    = openqcd_to_qdp_map[i][1];
      psi[i] = recvbuf[source_buf][buf_idx];
    }

    MPI_Waitall(qdp_to_openqcd_ndest, &sends[0], &statuses[0]);
      
  }


/* --------------------------------------------------------------------------*/
/**
* \brief  Copies a lattice fermion from openQCD to QDP.
*
* \param u: Destination field
* \param psi: Source spinor_dble filed.
*/
/* ----------------------------------------------------------------------------*/
  void copy_from_openQCD(LatticeFermionD & u, const spinor_dble * psi)
  {
    alloc_mappings();

    multi1d < FermionD > ferms(QDP::Layout::sitesOnNode());

    multi1d<multi1d<spinor_dble> > sendbuf(openqcd_to_qdp_ndest);
    multi1d<multi1d<spinor_dble> > recvbuf(qdp_to_openqcd_ndest);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      sendbuf[i].resize(openqcd_to_qdp_counts[i]);
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      recvbuf[i].resize(qdp_to_openqcd_counts[i]);
    for (int i = 0; i < L0*L1*L2*L3; i++) {
      int dest_buf = openqcd_to_qdp_map[i][0];
      int buf_idx  = openqcd_to_qdp_map[i][1];
      sendbuf[dest_buf][buf_idx] = psi[i];
    }
    multi1d<MPI_Request> sends(openqcd_to_qdp_ndest);
    multi1d<MPI_Request> recvs(qdp_to_openqcd_ndest);
    int ctag = mpi_tag();
    for (int i = 0; i < openqcd_to_qdp_ndest; i++)
      MPI_Isend(&sendbuf[i][0],
		openqcd_to_qdp_counts[i] * sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, openqcd_to_qdp_dest[i],
		ctag, MPI_COMM_WORLD, &sends[i]);
    for (int i = 0; i < qdp_to_openqcd_ndest; i++)
      MPI_Irecv(&recvbuf[i][0],
		qdp_to_openqcd_counts[i] * sizeof(spinor_dble)/sizeof(double),
		MPI_DOUBLE, qdp_to_openqcd_dest[i],
		ctag, MPI_COMM_WORLD, &recvs[i]);
    multi1d<MPI_Status> statuses(max(openqcd_to_qdp_ndest,
				     qdp_to_openqcd_ndest));
    MPI_Waitall(qdp_to_openqcd_ndest, &recvs[0], &statuses[0]);
    for (int i = 0; i < Layout::sitesOnNode(); i++) {
      int source_buf = qdp_to_openqcd_map[i][0];
      int buf_idx    = qdp_to_openqcd_map[i][1];
      FermionD tmpferm;
      copy_spinor_dble(tmpferm, recvbuf[source_buf][buf_idx]);
      ferms[i] = Gamma(2) * tmpferm;
    }
    QDP_insert(u, ferms, all);

    MPI_Waitall(openqcd_to_qdp_ndest, &sends[0], &statuses[0]);
      
  }

}
