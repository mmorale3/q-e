/*
 * Copyright (C) 2004 PWSCF group 
 * Copyright (C) 2007 QMCPACK developers
 *
 * @author Jeongnim Kim http://www.mcc.uiuc.edu/qmcpack/
 * @brief Implements generic hdf5 interfaces for plane wave codes and qmcpack
 *
 */

#define F77_FUNC_(name,NAME) name ## _

#ifdef __linux__
#include <stdio.h>
#include <sys/sysinfo.h>
#include <sys/resource.h>
#endif

void F77_FUNC_(print_freemem,PRINT_FREEMEM)()
{
#ifdef __linux__
  struct sysinfo si;
  sysinfo(&si);
  si.freeram += si.bufferram;
  unsigned long mem = si.freeram >> 20;
  printf(" available memory (GB): %g\n",(double)mem/1024.0);
#else
  printf(" memory report not available\n");
#endif
}

#if defined(__HDF5) || defined(__HDF5_C)

#include <stdio.h> 
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>
#include <assert.h>
#include "hdf5.h"
#include "hdf5_hl.h"


static H5E_auto_t err_func;
static void *client_data=0;

/** create a file and write version & application
 * @param fname name of the output file
 * @param length size of the file name
 *
 * h_file is initialized.
 */
void F77_FUNC_(esh5_posthf_open_file,ESH5_POSTHF_OPEN_FILE)
( hid_t* h_file, const char* fname, const int* length, int* old)
{
  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ; 

  //hid_t access_plist = H5Pcreate(H5P_FILE_ACCESS);
  //H5Pset_fclose_degree(access_plist, H5F_CLOSE_SEMI);
  //H5Pset_fclose_degree(access_plist, H5F_CLOSE_STRONG);

  if(*h_file>=0) H5Fclose(*h_file); 
  *h_file = H5Fopen(hfname,H5F_ACC_RDWR,H5P_DEFAULT);

  if( *old == 0 ) {
    // always delete the already existing file.
    if(*h_file>=0)
    {
      // always delete the already existing file.
      printf("esh5 destory the existing %s\n",hfname);
      remove(hfname);
      *h_file=-1;
    }
    //printf("esh5 create %s\n",hfname);
    *h_file = H5Fcreate(hfname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  }
  *old=0;

  if(*h_file>=0)
  {
    /* impelements version 1.00 hdf5 format */
    int version[]={2,1,0};
    hsize_t dim=3;
    herr_t ret=H5LTmake_dataset(*h_file,"version",1,&dim,H5T_NATIVE_INT,version);
    hsize_t ns=1;
    {
      hid_t strtype = H5Tcopy (H5T_C_S1);
      ret = H5Tset_size (strtype, 7); /* create string of length 5 */
      ret=H5LTmake_dataset(*h_file,"format",1,&ns,strtype,"ES-HDF");
    }

    hid_t h_app = H5Gcreate(*h_file,"application",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    {
      hid_t strtype = H5Tcopy (H5T_C_S1);
      ret = H5Tset_size (strtype, 8); /* create string of length 5 */
      ret=H5LTmake_dataset(h_app,"code",1,&ns,strtype,"espresso");
    }
    version[0]=4;
    version[2]=4;
    ret=H5LTmake_dataset(h_app,"version",1,&dim,H5T_NATIVE_INT,version);
    H5Gclose(h_app);
  } else {
    *old=1;  // error message
  }  

  free(hfname);
}

void F77_FUNC_(esh5_posthf_open_file_read,ESH5_POSTHF__OPEN_FILE_READ)
(hid_t* h_file, const char* fname, const int* length, int* error)
{
  *error=0;
  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);
  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  if(*h_file>=0) H5Fclose(*h_file);
  *h_file = H5Fopen(hfname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(*h_file<0)
    *error=*h_file;
  free(hfname);
}

void F77_FUNC_(esh5_posthf_close_file,ESH5_POSTHF_CLOSE_FILE)(hid_t* h_file)
{
  if(*h_file>=0) H5Fclose(*h_file);
  *h_file=-1;
  H5Eset_auto (H5E_DEFAULT,err_func, client_data);
}

/** write one body hamiltonian 
 * @param m   number of bands for kpoint ik 
 * @param ik  kpoint index
 * @param H1  One body hamiltonian 
 */
void F77_FUNC_(esh5_posthf_write_h1,ESH5_POSTHF_WRITE_H1)
  (hid_t* h_file, const int* m, const int* ik, const double* H1 )
{
  char aname[64];
  sprintf(aname,"H1_kp%i",(*ik)-1);

  hsize_t dims[3];
  dims[0] = (hsize_t)(*m);
  dims[1] = (hsize_t)(*m);
  dims[2] = 2;

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham < 0) ham = H5Gcreate(*h_file,"Hamiltonian",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  herr_t ret=H5LTmake_dataset(ham,aname,3,dims,H5T_NATIVE_DOUBLE,H1);
  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_read_h1,ESH5_POSTHF_READ_H1)
  (hid_t* h_file, const int* ik, double* H1, int* error )
{
  *error=0;  

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham<0){*error=ham;printf("Error in esh5_posthf_read_h1 1\n");return;}

  char aname[64];
  sprintf(aname,"H1_kp%i",(*ik)-1);
  herr_t ret=H5LTread_dataset(ham,aname,H5T_NATIVE_DOUBLE,H1);
  if(ret<0){*error=ret;printf("Error in esh5_posthf_read_h1 2\n");return;}

  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_read_nchol,ESH5_POSTHF_READ_NCHOL)
  (hid_t* h_file, int* ncholQ, int* error)
{
  *error=0;
  char aname[64];

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham<0){*error=ham;printf("Error in esh5_posthf_read_nchol 1\n");return;}

  sprintf(aname,"NCholPerKP");
  herr_t ret=H5LTread_dataset(ham,aname,H5T_NATIVE_INT,ncholQ);
  if(ret<0){*error=ret;printf("Error in esh5_posthf_read_nchol 2\n");return;}

  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_read_norbk,ESH5_POSTHF_READ_NORBK)
  (hid_t* h_file, int* norbK, int* error)
{
  *error=0;
  char aname[64];

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham<0){*error=ham;printf("Error in esh5_posthf_read_norbK 1\n");return;}

  sprintf(aname,"NMOPerKP");
  herr_t ret=H5LTread_dataset(ham,aname,H5T_NATIVE_INT,norbK);
  if(ret<0){*error=ret;printf("Error in esh5_posthf_read_norbK 2\n");return;}

  H5Gclose(ham);
}


void F77_FUNC_(esh5_posthf_read_cholesky,ESH5_POSTHF_READ_CHOLESKY)
  (hid_t* h_file, const int* Q, const int* k0, int* nchol, const int* nijtot, double* A, int* error)
{
  *error=0;
  char aname[64];
  herr_t ret;

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham<0){*error=ham;printf("Error in esh5_posthf_read_cholesky 1\n");return;}
  hid_t kp = H5Gopen(ham,"KPFactorized",H5P_DEFAULT);
  if(kp<0){*error=kp;printf("Error in esh5_posthf_read_cholesky 2\n");return;}

  sprintf(aname,"L%i",*Q);
  hid_t dataset = H5Dopen (kp, aname,H5P_DEFAULT);
  if(dataset<0){*error=1;printf("Error in esh5_posthf_read_cholesky 4\n");return;}

  hsize_t offsets[3];
  hsize_t dims[3];
  dims[0] = (hsize_t) 1;
  dims[1] = (hsize_t) (*nchol)*(*nijtot);
  dims[2] = (hsize_t) 2;
  offsets[0] = (hsize_t) *k0;
  offsets[1] = (hsize_t) 0;
  offsets[2] = (hsize_t) 0;
  hid_t filespace = H5Dget_space(dataset);
  if(filespace<0){*error=1;printf("Error in esh5_posthf_read_cholesky 5\n");return;}
  ret      = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, dims, NULL);
  if(ret<0){*error=ret;printf("Error in esh5_posthf_read_cholesky 6\n");return;}

  dims[0] = (hsize_t) (*nchol)*(*nijtot);
  dims[1] = (hsize_t) 2;
  hid_t memspace = H5Screate_simple(2, dims, NULL);
  if(memspace<0){*error=1;printf("Error in esh5_posthf_join_all 7\n");return;}
  ret            = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, A);
  if(ret<0){*error=ret;printf("Error in esh5_posthf_join_all 8\n");return;}

  H5Dclose(memspace);
  H5Sclose(filespace);
  H5Dclose (dataset);
  H5Gclose(kp);
  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_cholesky,ESH5_POSTHF_CHOLESKY)
  (hid_t* h_file, const int* Q, const int* k0, const int* nkloc, const int* ij0, const int* nij, int* nchol, 
    const int* nktot, const int* nijtot, int* nchol_max, const double* A, int* error)
{
  *error=0;
  if(*nkloc > 1) {
    if( *nij != *nijtot ) {
      *error=1;
      printf("Error in esh5_posthf_cholesky: Inconsistent data distribution\n");
      return;  
    }
  }

  char aname[64];
  const hsize_t dim4=4;

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham < 0) ham = H5Gcreate(*h_file,"Hamiltonian",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t kp = H5Gopen(ham,"KPFactorized",H5P_DEFAULT);
  if(kp < 0) kp = H5Gcreate(ham,"KPFactorized",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  int Ldims[4];
  Ldims[0] = (*k0);
  Ldims[1] = (*nkloc);
  Ldims[2] = ((*nchol)*(*ij0));
  Ldims[3] = ((*nchol)*(*nij));
  sprintf(aname,"L%i_dims",(*Q)-1);
  herr_t ret=H5LTmake_dataset(kp,aname,1,&dim4,H5T_NATIVE_INT,Ldims);
/*
  const hsize_t dim1=1;
  herr_t ret=H5LTmake_dataset(kp,"k0",1,&dim1,H5T_NATIVE_INT,k0);
  ret=H5LTmake_dataset(kp,"nk",1,&dim1,H5T_NATIVE_INT,nkloc);
  int c0 = (*nchol)*(*ij0);
  ret=H5LTmake_dataset(kp,"c0",1,&dim1,H5T_NATIVE_INT,&c0);
  c0 = (*nchol)*(*nij);
  ret=H5LTmake_dataset(kp,"nc",1,&dim1,H5T_NATIVE_INT,&c0);
*/

  hsize_t dims[4];
  dims[0] = (hsize_t) (*nkloc);
  dims[1] = (hsize_t) ((*nchol)*(*nij));
  dims[2] = (hsize_t) (2);
  hid_t dataspace = H5Screate_simple (3, dims, NULL); 

  sprintf(aname,"L%i",(*Q)-1);
  hid_t dataset = H5Dcreate (kp, aname, H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t filespace = H5Dget_space(dataset);

  dims[0] = (hsize_t)(*nkloc);
  dims[1] = (hsize_t)(*nij);
  dims[2] = (hsize_t)(*nchol_max);
  dims[3] = (hsize_t) (2);
  hid_t memspace = H5Screate_simple(4, dims, NULL);

  dims[0] = (hsize_t)(*nkloc);
  dims[1] = (hsize_t)(*nij);
  dims[2] = (hsize_t)(*nchol);
  dims[3] = (hsize_t) (2);
  hsize_t offsets[4];
  offsets[0] = offsets[1] = offsets[2] = offsets[3] = (hsize_t) 0;
  ret            = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets, NULL, dims, NULL);
  if(ret<0){*error=1;printf("Error in esh5_posthf_cholesky: H5Sselect_hyperslab\n");return;}


  ret            = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, A);
  if(ret<0){*error=1;printf("Error in esh5_posthf_cholesky: H5Dwrite\n");return;}

  H5Dclose(memspace);
  H5Sclose(filespace);
  H5Dclose (dataset);
  H5Sclose (dataspace);
  H5Gclose(kp);
  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_cholesky_root,ESH5_POSTHF_CHOLESKY_ROOT)
  (hid_t* h_file, const int* Q, const int* k0, const int* nkloc, const int* ij0, const int* nij, int* nchol, 
    const int* nktot, const int* nijtot, int* nchol_max, const double* A, int* error)
{
  *error=0;
  if(*nkloc > 1) {
    if( *nij != *nijtot ) {
      *error=1;
      printf("Error in esh5_posthf_cholesky_root: unexpected data distribution\n");
      return;  
    }
  }

  char aname[64];
  sprintf(aname,"L%i",(*Q)-1);

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham < 0) ham = H5Gcreate(*h_file,"Hamiltonian",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  hid_t kp = H5Gopen(ham,"KPFactorized",H5P_DEFAULT);
  if(kp < 0) kp = H5Gcreate(ham,"KPFactorized",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  hsize_t dims[4];
  dims[0] = (hsize_t) (*nktot);
  dims[1] = (hsize_t) ((*nchol)*(*nijtot));
  dims[2] = (hsize_t) (2);
  hid_t dataspace = H5Screate_simple (3, dims, NULL); 

  hid_t dataset = H5Dcreate (kp, aname, H5T_NATIVE_DOUBLE, dataspace,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(dataset<0)
    dataset = H5Dopen (kp, aname,H5P_DEFAULT);
  if(dataset<0){*error=1;printf("Error in esh5_posthf_cholesky_root: dataset root\n");return;}
  
  dims[0] = (hsize_t)(*nkloc);
  dims[1] = (hsize_t)((*nij)*(*nchol));
  dims[2] = (hsize_t) (2);
  hsize_t offsets[4];
  offsets[0] = (hsize_t)(*k0);
  offsets[1] = (hsize_t)((*nchol)*(*ij0)); 
  offsets[2] = (hsize_t) 0;
  hid_t filespace = H5Dget_space(dataset);
  herr_t ret      = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, dims, NULL);
  if(ret<0){*error=1;printf("Error in esh5_posthf_cholesky_root: H5Sselect_hyperslab file\n");return;}

  dims[0] = (hsize_t)(*nkloc);
  dims[1] = (hsize_t)(*nij);
  dims[2] = (hsize_t)(*nchol_max);
  dims[3] = (hsize_t) (2);
  hid_t memspace = H5Screate_simple(4, dims, NULL);

  dims[0] = (hsize_t)(*nkloc);
  dims[1] = (hsize_t)(*nij);
  dims[2] = (hsize_t)(*nchol);
  dims[3] = (hsize_t) (2);
  offsets[0] = offsets[1] = offsets[2] = offsets[3] = (hsize_t) 0;
  ret            = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offsets, NULL, dims, NULL);
  if(ret<0){*error=1;printf("Error in esh5_posthf_cholesky_root: H5Sselect_hyperslab memspace\n");return;}


  ret            = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, A);
  if(ret<0){*error=1;printf("Error in esh5_posthf_cholesky_root: H5Dwrite\n");return;}

  H5Dclose(memspace);
  H5Sclose(filespace);
  H5Dclose (dataset);
  H5Sclose (dataspace);
  H5Gclose(kp);
  H5Gclose(ham);
}

void F77_FUNC_(esh5_posthf_join_all,ESH5_POSTHF_JOIN_ALL)
  (hid_t* h_file, const char* fname, const int* length, const int* nQ, const int* nproc, int* error)  
{
  *error=0;
  if(*nproc == 1) return;
  char aname[64];
  hsize_t dims[4];
  char * hfname = ( char * ) malloc( (*length) + 100 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  if(*h_file<0){*error=1;printf("Closed file in esh5_posthf_join_all\n");return;} 
  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham<0){*error=1;printf("Error in esh5_posthf_join_all 1\n");return;} 
  hid_t kp = H5Gopen(ham,"KPFactorized",H5P_DEFAULT);
  if(kp<0){*error=1;printf("Error in esh5_posthf_join_all 2\n");return;} 

  int datasz=0;
  double* data=NULL;

  for(int i=1; i<(*nproc); i++) {

    sprintf(hfname+*length,"_part%i\0",i);
    //printf(" Joining %s\n",hfname);
    hid_t ifile = H5Fopen(hfname,H5F_ACC_RDONLY,H5P_DEFAULT); 
    if(ifile<0){*error=1;printf("Error in esh5_posthf_join_all 3\n");return;} 

    hid_t iham = H5Gopen(ifile,"Hamiltonian",H5P_DEFAULT);
    if(iham<0){*error=1;printf("Error in esh5_posthf_join_all 4\n");return;}
    hid_t ikp = H5Gopen(iham,"KPFactorized",H5P_DEFAULT);
    if(ikp<0){*error=1;printf("Error in esh5_posthf_join_all 5\n");return;}

    for(int Q=0; Q<(*nQ); Q++) {

      int Ldims[4];
      int k0, nk, c0, nc;
      sprintf(aname,"L%i_dims",Q);
      herr_t ret=H5LTread_dataset(ikp,aname,H5T_NATIVE_INT,Ldims);
      // MAM: it is ok if LQs are missing, this means some symmetry
      //if(ret<0){*error=1;printf("Error in esh5_posthf_join_all 6\n");return;}
      if(ret<0) continue;
      k0 = Ldims[0];
      nk = Ldims[1];
      c0 = Ldims[2];
      nc = Ldims[3];

      int sz=2*nk*nc;
      if(datasz < sz) {
        if(datasz) free(data);
        data = (double*) malloc( sizeof(double)*sz );
      }        

      sprintf(aname,"L%i",Q);
      dims[0] = (hsize_t) nk;
      dims[1] = (hsize_t) nc;
      dims[2] = (hsize_t) 2;
      ret=H5LTread_dataset(ikp,aname,H5T_NATIVE_DOUBLE,data);
      if(ret<0){*error=1;printf("Error in esh5_posthf_join_all 10\n");return;}

      hid_t dataset = H5Dopen (kp, aname,H5P_DEFAULT);

      hsize_t offsets[3];
      offsets[0] = (hsize_t) k0;
      offsets[1] = (hsize_t) c0;
      offsets[2] = (hsize_t) 0;
      hid_t filespace = H5Dget_space(dataset);
      ret      = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsets, NULL, dims, NULL);
      if(ret<0){*error=ret;printf("Error in esh5_posthf_join_all 11\n");return;}

      hid_t memspace = H5Screate_simple(3, dims, NULL);
      ret            = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data);
      if(ret<0){*error=ret;printf("Error in esh5_posthf_join_all 12\n");return;}

      H5Dclose(memspace);
      H5Sclose(filespace);
      H5Dclose (dataset);

    }

    H5Gclose(ikp);
    H5Gclose(iham);
    H5Fclose(ifile);
    //printf(" Removing %s\n",hfname);
    remove(hfname);

  }
 
  H5Gclose(kp);
  H5Gclose(ham);
  if(datasz) free(data);
  free(hfname);

}

/** write one body hamiltonian 
 * @param nup   number of up electrons
 * @param ndown   number of down electrons 
 * @param nk   number of k points 
 * @param xk   list of k points 
 * @param kminus   Q -> -Q mapping
 * @param QKtoK2   {Q,K} -> {Ka, Kb} mapping
 * @param norbK   number of orbitals per K point
 * @param ncholQ   number of cholesky vectors per Q point 
 */
void F77_FUNC_(esh5_posthf_kpoint_info,ESH5_POSTHF_KPOINT_INFO)
  (hid_t* h_file, const int* nup, const int*ndown, const double* E0, const double* Efc, 
    const int* nk, const double* xk, 
    const int* kminus, const int* QKtoK2, const int* norbK, const int* ncholQ )
{
  int* qktok2_ = (int*) malloc((*nk)*(*nk)*sizeof(int));
  char aname[64];
  hsize_t dims[2];

  hid_t ham = H5Gopen(*h_file,"Hamiltonian",H5P_DEFAULT);
  if(ham < 0) ham = H5Gcreate(*h_file,"Hamiltonian",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  // kpoints
  sprintf(aname,"KPoints");
  dims[0] = (hsize_t)(*nk);
  dims[1] = (hsize_t)(3);
  herr_t ret=H5LTmake_dataset(ham,aname,2,dims,H5T_NATIVE_DOUBLE,xk);

  // kminus
  sprintf(aname,"MinusK");
  for(int i=0; i<(*nk); i++) qktok2_[i] = kminus[i]-1;
  dims[0] = (hsize_t)(*nk);
  ret=H5LTmake_dataset(ham,aname,1,dims,H5T_NATIVE_INT,qktok2_);

  // qktok2
  for(int i=0; i<(*nk); i++) {
    for(int j=0; j<(*nk); j++)
      qktok2_[i*(*nk)+j] = QKtoK2[j*(*nk)+i]-1;
  }
  sprintf(aname,"QKTok2");
  dims[0] = (hsize_t)(*nk);
  dims[1] = (hsize_t)(*nk);
  ret=H5LTmake_dataset(ham,aname,2,dims,H5T_NATIVE_INT,qktok2_);
  free(qktok2_);

  // norb_per_kp
  sprintf(aname,"NMOPerKP");
  dims[0] = (hsize_t)(*nk);
  ret=H5LTmake_dataset(ham,aname,1,dims,H5T_NATIVE_INT,norbK);

  // nchol_per_Q
  sprintf(aname,"NCholPerKP");
  dims[0] = (hsize_t)(*nk);
  ret=H5LTmake_dataset(ham,aname,1,dims,H5T_NATIVE_INT,ncholQ);

  // energies
  sprintf(aname,"Energies");
  dims[0] = (hsize_t)(2);
  double E_[2];
  E_[0]=*E0;
  E_[1]=*Efc;
  ret=H5LTmake_dataset(ham,aname,1,dims,H5T_NATIVE_DOUBLE,E_);

  // dims 
  int data[8];
  data[0]=data[1]=data[6]=data[7]=data[3]=0;
  data[2]=(*nk);
  for(int i=0; i<(*nk); i++) data[3] += norbK[i];
  data[4]=(*nup)*(*nk);
  data[5]=(*ndown)*(*nk);
  sprintf(aname,"dims");
  dims[0] = (hsize_t)(8);
  ret=H5LTmake_dataset(ham,aname,1,dims,H5T_NATIVE_INT,data);

  H5Gclose(ham);

}

void write_sparse_matrix(hid_t hg, const int nr, const int nc, const double* A);
void generate_SM(const int nspin, const int nk, const int* norbK, const int norbmax, const int nelmax, 
  const int nfull, const int* fullocc, const int npart, const int* partocc,
  const int mixed, const int, int*, const double* M, double* SM); 

/** write one body hamiltonian 
 * @param nspin        number of spin components 
 * @param nk           number of k points
 * @param norbK         number of orbitals per kpoint 
 * @param nelmax       electron dimension of Psi 
 * @param ndets        requested number of determinants
 * @param wg           weights of Kohn-Sham states in trial state (possibly degenerate) 
 *                     wg(ispin, ik, n), n in {0,nelmax-1} 
 * @param mixed        flag indicating if M is assumed identity (=0) or not 
 * @param M            Overlap matrix M(spin, ik, n, iorb) 
 *                      = < psi(ik, iorb) | Psi_DFT( spin, ik, n ) >,  n in {0,nelmax-1}
 *                      where |psi(ik,iorb)> are the basis set elements
 *                      and | Psi_DFT( spin, ik, n ) > the Kohn-Sham states
 *                      that correspond to the weights in wg
 * @param nkocc_ref     Returns the number of occupied orbitals per kpoint per spin
 *                     in the reference slater determinant    
 */
void F77_FUNC_(esh5_posthf_write_wavefunction,ESH5_POSTHF_WRITE_WAVEFUNCTION)
  (hid_t* h_file, const int* nspin, const int* nk, const int* norbK, const int* norbmax,
    const int* nelmax, const int* ndets, const double* wg, 
    const int* mix_, const double* lowcut, const double* highcut, 
    const double* M, int* nkocc_ref, int* error)
{
  int noncol_fac=1;
  if(*nspin==4) noncol_fac=2;
  int nspinSM=1;
  if(*nspin==4) nspinSM=2;
  const int mixed = *mix_;
  *error=0;
  double accn=0;
  int nfullocc[2], npartocc[2];
  nfullocc[0] = nfullocc[1] = npartocc[0] = npartocc[1] = 0;
  for(int ik=0, in=0; ik<(*nk); ik++) 
   for(int n=0; n<(*nelmax); n++, in++) { 
     if( fabs(wg[in]) >= *highcut ) nfullocc[0]++; 
     else if( fabs(wg[in]) >= *lowcut ) npartocc[0]++; 
     accn += wg[in];
   }  
  int nup = (int) round(accn), ndown=0; 
  if( *nspin == 2 ) {
    accn=0;
    for(int i=0, in=(*nk)*(*nelmax); i<(*nk); i++) 
      for(int n=0; n<(*nelmax); n++, in++) { 
        if( fabs(wg[in]) >= *highcut ) nfullocc[1]++; 
        else if( fabs(wg[in]) >= *lowcut ) npartocc[1]++; 
        accn += wg[in];
    } 
    ndown = (int) round(accn); 
  }
  printf(" esh5_posthf_write_wavefunction: \n"); 
  if( *nspin == 1)  {
    printf("    - Closed-shell calculation:");
    printf("       - nelec, full occ, partial occ: %i %i %i\n", nup, nfullocc[0],npartocc[0]);
  } else if( *nspin == 2)  {
    printf("    - Collinear calculation:");
    printf("    - Alpha orbitals: \n");
    printf("       - nelec, full occ, partial occ: %i %i %i\n", nup, nfullocc[0],npartocc[0]);
    printf("    - Beta orbitals: \n"); 
    printf("       - nelec, full occ, partial occ: %i %i %i\n", ndown, nfullocc[1],npartocc[1]);
  } else if(*nspin == 4) {
    printf("    - Non-collinear calculation: %i\n", nup);
    printf("       - nelec, full occ, partial occ: %i %i %i\n", nup, nfullocc[0],npartocc[0]);
  } else  {
    *error=1;
    return;
  } 

  int* focc = (int*) malloc(sizeof(int)*(nfullocc[0] + nfullocc[1]));
  int* pocc = (int*) malloc(sizeof(int)*(npartocc[0] + npartocc[1]));
  double* pocc_w = (double*) malloc(sizeof(double)*(npartocc[0] + npartocc[1]));

  for(int ik=0, in=0, io=0,ip=0; ik<(*nk); ik++)
    for(int n=0; n<(*nelmax); n++, in++) {
      if( fabs(wg[in]) >= *highcut ) {
        focc[io] = ik*(*nelmax)+n;  
        io++;
      } else if( fabs(wg[in]) >= *lowcut ) {
        pocc[ip] = ik*(*nelmax)+n;  
        ip++;
      }
    }  
  if( *nspin == 2 ) {
    accn=0;
    for(int ik=0, in=(*nk)*(*nelmax), io= nfullocc[0], ip=npartocc[0]; ik<(*nk); ik++) 
      for(int n=0; n<(*nelmax); n++, in++) {
        if( fabs(wg[in]) >= *highcut ) {
          focc[io] = ik*(*nelmax)+n; 
          io++;
        } else if( fabs(wg[in]) >= *lowcut ) {
          pocc[ip] = ik*(*nelmax)+n; 
          ip++;
        }
      }
  }
  
  int norbtot = 0;
  for(int i=0; i<(*nk); i++) norbtot+=norbK[i];

  int maxnel = nup;
  if( ndown > nup ) maxnel = ndown; 
  double* SM = (double*) malloc(sizeof(double)*2*maxnel*noncol_fac*norbtot);

  int nrem[2];
  nrem[0] = nup - nfullocc[0];
  nrem[1] = ndown - nfullocc[1];

  hid_t wfn = H5Gopen(*h_file,"Wavefunction",H5P_DEFAULT);
  if(wfn < 0) wfn = H5Gcreate(*h_file,"Wavefunction",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  if( nrem[0] == 0 && npartocc[0] > 0 ) {
    printf(" Error: # full occupied (up) == nup, with partially occupied states.\n");
    *error=1;
    return;
  }
  if( nrem[1] == 0 && npartocc[1] > 0 ) {
    printf(" Error: # full occupied (down) == nup, with partially occupied states.\n");
    *error=1;
    return;
  }
  if(*nspin > 1 && *mix_==0) {
    printf(" Error: nspin>1 && mixed==0.\n");
    *error=1;
    return;
  }

  char aname[64];
  hsize_t dims[3];

  if( (nrem[0] == 0 && nrem[1] == 0) ) {

    // write single determinant in NOMSD form
    hid_t msd = H5Gopen(wfn,"NOMSD",H5P_DEFAULT);
    if(msd < 0) msd = H5Gcreate(wfn,"NOMSD",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

    double ci[2];
    ci[0]=1.0; ci[1]=0.0;
    sprintf(aname,"ci_coeffs");
    dims[0] = (hsize_t)(1);
    dims[1] = (hsize_t)(2);
    herr_t ret=H5LTmake_dataset(msd,aname,2,dims,H5T_NATIVE_DOUBLE,ci);

// NMO,NAEA,NAEB,walker_type,ndets_to_read
    int data[5]; 
    data[0]=norbtot;
    data[1]=nup;
    data[4]=1;
    if(*nspin==1) {
      data[2]=nup;
      data[3]=1;
    } else if(*nspin==2) {
      data[2]=ndown;
      data[3]=2;
    } else if(*nspin==4) {
      data[2]=0; 
      data[3]=3;
    } else { 
      printf(" Inconsistent nspin value: %i\n",*nspin);
      *error=1;
      return;
    }
    sprintf(aname,"dims");
    dims[0] = (hsize_t)(5);
    ret=H5LTmake_dataset(msd,aname,1,dims,H5T_NATIVE_INT,data);

    generate_SM(nspinSM,*nk,norbK,*norbmax,*nelmax,nfullocc[0],focc,0,pocc,mixed,
                1,nkocc_ref,M,SM); 

    sprintf(aname,"Psi0_alpha");
    dims[0] = (hsize_t)(norbtot*noncol_fac);
    dims[1] = (hsize_t)(nup);
    dims[2] = (hsize_t)(2);
    ret=H5LTmake_dataset(msd,aname,3,dims,H5T_NATIVE_DOUBLE,SM);

    // PsiT Alpha 
    hid_t psiT = H5Gopen(msd,"PsiT_0",H5P_DEFAULT);
    if(psiT < 0) psiT = H5Gcreate(msd,"PsiT_0",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
    write_sparse_matrix(psiT,norbtot*noncol_fac,nup,SM);
    H5Gclose(psiT);

    if(*nspin == 2) {

      generate_SM(1,*nk,norbK,*norbmax,*nelmax,nfullocc[1],focc+nfullocc[0],
                    0,pocc,mixed,1,nkocc_ref+(*nk),M+2*(*nk)*(*norbmax)*(*nelmax),SM);

      sprintf(aname,"Psi0_beta");
      dims[0] = (hsize_t)(norbtot);
      dims[1] = (hsize_t)(ndown);
      dims[2] = (hsize_t)(2);
      herr_t ret=H5LTmake_dataset(msd,aname,3,dims,H5T_NATIVE_DOUBLE,SM);
 
      // PsiT Beta 
      hid_t psiT = H5Gopen(msd,"PsiT_1",H5P_DEFAULT);
      if(psiT < 0) psiT = H5Gcreate(msd,"PsiT_1",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
      write_sparse_matrix(psiT,norbtot,ndown,SM);
      H5Gclose(psiT);

    }

    H5Gclose(msd);

  } else {
    // write in PHMSD format, if desired in NOMSD format, use converter script
    *error=1;
    return; 

    hid_t msd = H5Gopen(wfn,"PHMSD",H5P_DEFAULT);
    if(msd < 0) msd = H5Gcreate(wfn,"PHMSD",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

// call hpsort()

    // 1. obtain strings of occupation numbers from wg representing
    //    the DFT states occupied in each determinant in the list
    //    - list is sorted according to product of weights in wg

    // 2. Construct the appropriate Slater Matrices for each set of
    //    occupation numbers and output in sparse format

  
    H5Gclose(msd);

  }

  free(focc);
  free(pocc);
  free(pocc_w);
  free(SM);
  H5Gclose(wfn);

  return;

}

// Note: QMCPACK expects SM.conj().T, so writing hermitian conjugate here
void write_sparse_matrix(hid_t hg, const int nr, const int nc, const double* A)
{
  int info[3]; // nr, nc, nnz
  int nnz=0;
  int* pbeg = (int*) malloc(sizeof(int)*(nc+1));

  for(int i=0, ij=0; i<nr; i++)
    for(int j=0; j<nc; j++, ij++)
      if( sqrt(A[2*ij]*A[2*ij] + A[2*ij+1]*A[2*ij+1]) > 1e-8 ) nnz++; 

  int* jdata = (int*) malloc(sizeof(int)*(nnz));
  double* data = (double*) malloc(sizeof(double)*2*nnz);

  int cnt=0;
  for(int i=0; i<nc; i++) {
    pbeg[i] = cnt;
    for(int j=0; j<nr; j++) {
      int ij = j*nc+i;
      if( sqrt(A[2*ij]*A[2*ij] + A[2*ij+1]*A[2*ij+1]) > 1e-8 ) {
        data[2*cnt] = A[2*ij];
        data[2*cnt+1] = -A[2*ij+1];
        jdata[cnt] = j;
        cnt++;
      } 
    }
  }
  pbeg[nc] = cnt;

  char aname[64];
  hsize_t dims[2];
  dims[0]=nnz;
  dims[1]=2;
  sprintf(aname,"data_");
  herr_t ret=H5LTmake_dataset(hg,aname,2,dims,H5T_NATIVE_DOUBLE,data);
  sprintf(aname,"jdata_");
  ret=H5LTmake_dataset(hg,aname,1,dims,H5T_NATIVE_INT,jdata);
  dims[0]=nc;
  sprintf(aname,"pointers_begin_");
  ret=H5LTmake_dataset(hg,aname,1,dims,H5T_NATIVE_INT,pbeg);
  sprintf(aname,"pointers_end_");
  ret=H5LTmake_dataset(hg,aname,1,dims,H5T_NATIVE_INT,pbeg+1);
  dims[0]=3;
  info[0]=nc;
  info[1]=nr;
  info[2]=nnz;
  sprintf(aname,"dims");
  ret=H5LTmake_dataset(hg,aname,1,dims,H5T_NATIVE_INT,info);

  free(data);
  free(jdata);
  free(pbeg);
}

/** 
 */
void F77_FUNC_(esh5_posthf_open_read,ESH5_POSTHF_OPEN_READ)
    (hid_t* h_file, const char* fname, const int* length, int* nk, int* norbK, int* nmax_dm, 
     int* nspin, int* npol, int* npwx,
     double* xk, int* grid_type, int* nr1, int* nr2, int* nr3, double* recvec, int* error)
{
  *error=0;
  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  if(*h_file>=0) H5Fclose(*h_file);
  *h_file = H5Fopen(hfname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(*h_file<0)
  {
    printf("esh5 error opening orbital file. %s\n",hfname);
    *h_file=-1;
    *error = 1;
    return;
  }

  herr_t ret;
  hid_t g = H5Gopen(*h_file,"DM",H5P_DEFAULT);
  if(g >= 0) {
    int ndm=-1;
    ret=H5LTread_dataset(g,"ndim",H5T_NATIVE_INT,&ndm);
    if(ret >= 0) *nmax_dm = ndm;
    H5Gclose(g);
  }

  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsG",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1; 
    printf("esh5 error opening OrbsG group.\n");
    return;
  }   

  int npwx_=0;
  ret=H5LTread_dataset(h_orb_grp_read,"npwx",H5T_NATIVE_INT,&npwx_);
  if(ret>=0) { *npwx = npwx_; }

  int npol_=0;
  *npol = 1;
  ret=H5LTread_dataset(h_orb_grp_read,"npol",H5T_NATIVE_INT,&npol_);
  if(ret>=0) { *npol = npol_; }

  int nspin_=1;
  *nspin = 1;
  ret=H5LTread_dataset(h_orb_grp_read,"nspin",H5T_NATIVE_INT,&nspin_);
  if(ret>=0) { *nspin = nspin_; }

  ret=H5LTread_dataset(h_orb_grp_read,"grid_type",H5T_NATIVE_INT,grid_type);
  if(ret<0){*error=1;printf("esh5 problems reading grid_type\n");return;}   

  ret=H5LTread_dataset(h_orb_grp_read,"number_of_kpoints",H5T_NATIVE_INT,nk);
  if(ret<0){*error=1;printf("esh5 problems reading number_of_kpoints\n");return;}   

  ret=H5LTread_dataset(h_orb_grp_read,"kpoints",H5T_NATIVE_DOUBLE,xk);
  if(ret<0){*error=1;printf("esh5 problems reading kpoints\n");return;}   

  int mesh[3];  
  ret=H5LTread_dataset(h_orb_grp_read,"fft_grid",H5T_NATIVE_INT,mesh);
  if(ret < 0){*error=1;printf("esh5 error reading fft_grid\n");return;}   
  *nr1=mesh[0]; 
  *nr2=mesh[1]; 
  *nr3=mesh[2]; 

  ret=H5LTread_dataset(h_orb_grp_read,"reciprocal_vectors",H5T_NATIVE_DOUBLE,recvec);
  if(ret < 0){*error=1;printf("esh5 error reading reciprocal_vectors\n");return;}   

  ret=H5LTread_dataset(h_orb_grp_read,"number_of_orbitals",H5T_NATIVE_INT,norbK);
  if(ret < 0){*error=1;printf("esh5 error reading number_of_orbitals 1\n");return;}   

  H5Gclose(h_orb_grp_read);

  free(hfname);
}

/** 
 */
void F77_FUNC_(esh5_posthf_open_read_r,ESH5_POSTHF_OPEN_READ_R)
    (hid_t* h_file, const char* fname, const int* length, int* error) 
{
  *error=0;
  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  if(*h_file>=0) H5Fclose(*h_file);
  *h_file = H5Fopen(hfname,H5F_ACC_RDONLY,H5P_DEFAULT);
  if(*h_file<0)
  {
    printf("esh5 error opening orbital file. %s\n",hfname);
    *h_file=-1;
    *error = 1;
    return;
  }

  // check group exists
  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsR",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1; 
    printf("esh5 error opening OrbsR group.\n");
    return;
  }   
  H5Gclose(h_orb_grp_read);

  free(hfname);
}

void F77_FUNC_(esh5_posthf_close_read,ESH5_POSTHF_CLOSE_READ)(hid_t* h_file)
{
  if(*h_file>=0) H5Fclose(*h_file);
  *h_file=-1;
  H5Eset_auto (H5E_DEFAULT,err_func, client_data);
}

void F77_FUNC_(esh5_posthf_read,ESH5_POSTHF_READ)
  (hid_t* h_file,  const int* ik, const int* ib, const int* nb, double* Psir, const int* lda, int* error)
{
  *error=0;
  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsG",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1;
    printf("esh5 error opening OrbsG group.\n");
    return;
  }
  char aname[64];
  for( int i=0, iend=(*nb), inc=0; i<iend; i++, inc+=2*(*lda)) {
    int bnd = (*ib) + i;
    sprintf(aname,"kp%i_b%i\0",*ik,bnd);
    herr_t ret=H5LTread_dataset(h_orb_grp_read,aname,H5T_NATIVE_DOUBLE,Psir+inc);
    if(ret < 0){*error=ret;printf("esh5 error reading orbital\n");return;}   
  }  
  H5Gclose(h_orb_grp_read);
}

void F77_FUNC_(esh5_posthf_read_norb,ESH5_POSTHF_READ_NORB)
  (hid_t* h_file, int* nb, int* error)
{
  *error=0;
  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsG",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1;
    printf("esh5 error opening OrbsG group.\n");
    return;
  }
  herr_t ret=ret=H5LTread_dataset(h_orb_grp_read,"number_of_orbitals",H5T_NATIVE_INT,nb);
  if(ret < 0){*error=1;printf("esh5 error reading number_of_orbitals\n");return;}
  H5Gclose(h_orb_grp_read);
}

void F77_FUNC_(esh5_posthf_read_nkpts,ESH5_POSTHF_READ_NKPTS)
  (hid_t* h_file, int* nkpts, int* error)
{
  *error=0;
  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsG",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1;
    printf("esh5 error opening OrbsG group.\n");
    return;
  }
  herr_t ret=ret=H5LTread_dataset(h_orb_grp_read,"number_of_kpoints",H5T_NATIVE_INT,nkpts);
  if(ret<0){*error=1;printf("esh5 problems reading number_of_kpoints\n");return;}
  H5Gclose(h_orb_grp_read);
}

void F77_FUNC_(esh5_posthf_open_write,ESH5_POSTHF_OPEN_WRITE)
  (hid_t* h_file, const char* fname, const int* length, int* error ) 
{
  *error=0;
  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  if(*h_file>=0) H5Fclose(*h_file);
  *h_file = H5Fopen(hfname,H5F_ACC_RDWR,H5P_DEFAULT);
  if(*h_file>=0)
  {
    printf("esh5 destory the existing orbital file %s\n",hfname);
    remove(hfname);
    *h_file=-1;
  }

  *h_file = H5Fcreate(hfname,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  if(*h_file<0) {
    printf("esh5 error opening orbital write file. %s\n",hfname);
    *h_file=-1;
    *error = 1;
    free(hfname);
    return;
  }
  free(hfname);
}

/** write orbitals 
 * @param 
 *
 */
void F77_FUNC_(esh5_posthf_write_meta,ESH5_POSTHF_WRITE_META)
  ( hid_t* h_file,  const char* dname, const int* dlength,
    const int* nk, const int* nspin, const int* npol, const int* npwx, 
    const double* xk, const int* grid_type, const int* nr1, const int* nr2, 
    const int* nr3, const double* lattvec, const double* recvec,
    const double* alat,
    int* error)
{
  *error=0;
  hsize_t dims[3];

  char * hdname = ( char * ) malloc( (*dlength) + 1 ) ;
  memcpy( hdname , dname , *dlength ) ;
  hdname[*dlength] = '\0' ;

  hid_t h_orb_grp_write = H5Gopen(*h_file,hdname,H5P_DEFAULT);
  if(h_orb_grp_write < 0) h_orb_grp_write = H5Gcreate(*h_file,hdname,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT); 
  if(h_orb_grp_write < 0) {
    h_orb_grp_write=-1;
    *error=1;
    printf("esh5 error opening orbital group.\n");
    free(hdname);
    return;
  }
  free(hdname);

  const hsize_t dim1=1;
  const hsize_t dim3=3;
  if(*npwx > 0) {
    herr_t ret=H5LTmake_dataset(h_orb_grp_write,"npwx",1,&dim1,H5T_NATIVE_INT,npwx);
    if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}
  }

  herr_t ret=H5LTmake_dataset(h_orb_grp_write,"npol",1,&dim1,H5T_NATIVE_INT,npol);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  ret=H5LTmake_dataset(h_orb_grp_write,"nspin",1,&dim1,H5T_NATIVE_INT,nspin);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  ret=H5LTmake_dataset(h_orb_grp_write,"grid_type",1,&dim1,H5T_NATIVE_INT,grid_type);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  ret=H5LTmake_dataset(h_orb_grp_write,"number_of_kpoints",1,&dim1,H5T_NATIVE_INT,nk);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}
  int data[3]; data[0]=*nr1; data[1]=*nr2; data[2]=*nr3;
  ret=H5LTmake_dataset(h_orb_grp_write,"fft_grid",1,&dim3,H5T_NATIVE_INT,data);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  dims[0] = (hsize_t)((*nk));
  dims[1] = (hsize_t)(3);
  ret=H5LTmake_dataset(h_orb_grp_write,"kpoints",2,dims,H5T_NATIVE_DOUBLE,xk);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  dims[0] = (hsize_t)(3);
  dims[1] = (hsize_t)(3);
  ret=H5LTmake_dataset(h_orb_grp_write,"reciprocal_vectors",2,dims,H5T_NATIVE_DOUBLE,recvec);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  ret=H5LTmake_dataset(h_orb_grp_write,"lattice_vectors",2,dims,H5T_NATIVE_DOUBLE,lattvec);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  ret=H5LTmake_dataset(h_orb_grp_write,"alat",1,&dim1,H5T_NATIVE_DOUBLE,alat);
  if(ret<0){*error=1;printf("esh5 problems writing to orbital file\n");return;}

  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_write_norb,ESH5_POSTHF_WRITE_NORB)
  (hid_t* h_file, const char* dname, const int* dlength, const int* nk, const int* norbK, int* error)
{
  *error=0;
  char * hdname = ( char * ) malloc( (*dlength) + 1 ) ;
  memcpy( hdname , dname , *dlength ) ;
  hdname[*dlength] = '\0' ;

  hid_t h_orb_grp_write = H5Gopen(*h_file,hdname,H5P_DEFAULT);
  if(h_orb_grp_write < 0) {
    h_orb_grp_write=-1;
    *error=1;
    printf("esh5 error opening orbital group.\n");
    free(hdname);
    return;
  }
  free(hdname);

  hsize_t dims[1];
  dims[0] = (hsize_t)(*nk);
  herr_t ret=H5LTmake_dataset(h_orb_grp_write,"number_of_orbitals",1,dims,H5T_NATIVE_INT,norbK);

  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_close_write,ESH5_POSTHF_CLOSE_WRITE)(hid_t* h_file)
{
  if(*h_file>=0) H5Fclose(*h_file);
  *h_file=-1;
  H5Eset_auto (H5E_DEFAULT,err_func, client_data);
}

void F77_FUNC_(esh5_posthf_write,ESH5_POSTHF_WRITE)
  (hid_t* h_file, const char* dname, const int* dlength, const int* ik, const int* npwx, const int* nb, const double* Psir, const int* lda, int* error)
{
  *error=0;
  char * hdname = ( char * ) malloc( (*dlength) + 1 ) ;
  memcpy( hdname , dname , *dlength ) ;
  hdname[*dlength] = '\0' ;

  hid_t h_orb_grp_write = H5Gopen(*h_file,hdname,H5P_DEFAULT);
  if(h_orb_grp_write < 0) {
    h_orb_grp_write=-1;
    *error=1;
    printf("esh5 error opening orbital group.\n");
    free(hdname);
    return;
  }
  free(hdname);

  char aname[64];
  hsize_t dims[2];
  dims[0] = (hsize_t)((*npwx));
  dims[1] = 2;
  for( int i=0, iend=(*nb), inc=0; i<iend; i++, inc+=2*(*lda)) {
    sprintf(aname,"kp%i_b%i",*ik,i);
    herr_t ret=H5LTmake_dataset(h_orb_grp_write,aname,2,dims,H5T_NATIVE_DOUBLE,Psir+inc);
    if(ret < 0){*error=1;printf("esh5 error writing orbital\n");}
  }

  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_write_g, ESH5_POSTHF_WRITE_G)
  (hid_t* h_file, const double* gk, const int* npw_in, const int* ik_in,
   int* error)
{
  const int ndim=3;
  const int npw = *npw_in;
  const int ik = *ik_in;
  char aname[64];
  hsize_t dims[2];
  dims[0] = (hsize_t)npw;
  dims[1] = ndim;

  // !!!! assume OrbsG exists
  *error = 0;
  hid_t h_orb_grp_write = H5Gopen(*h_file, "OrbsG",H5P_DEFAULT);
  sprintf(aname,"kp%i_g", ik);
  herr_t ret = H5LTmake_dataset(h_orb_grp_write, aname, 2, dims, H5T_NATIVE_DOUBLE, gk);
  if (ret<0)
  {
    *error = 1;
    printf("esh5 error writting g (tpiba)\n");
  }
  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_write_et,ESH5_POSTHF_WRITE_ET)
 (hid_t* h_file, const char* dname, const int* dlength, const int* ik, const int* m, const double* et, const double* wg, int* error) 
{
  *error=0;
  char * hdname = ( char * ) malloc( (*dlength) + 1 ) ;
  memcpy( hdname , dname , *dlength ) ;
  hdname[*dlength] = '\0' ;
  
  hid_t h_orb_grp_write = H5Gopen(*h_file,hdname,H5P_DEFAULT);
  if(h_orb_grp_write < 0) {
    h_orb_grp_write=-1;
    *error=1;
    printf("esh5 error opening orbital group.\n");
    free(hdname);
    return;
  }
  free(hdname);

  char aname[64];
  hsize_t dims[1];
  dims[0] = (hsize_t)(*m);

  sprintf(aname,"EigVal_kp%i",(*ik));
  herr_t ret=H5LTmake_dataset(h_orb_grp_write,aname,1,dims,H5T_NATIVE_DOUBLE,et);

  sprintf(aname,"wg_kp%i",(*ik));
  ret=H5LTmake_dataset(h_orb_grp_write,aname,1,dims,H5T_NATIVE_DOUBLE,wg);
  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_read_et,ESH5_POSTHF_READ_ET)
 (hid_t* h_file, const int* ik, double* et, double* wg, int* error)
// ( const char* dname, const int* dlength, const int* ik, double* et, double* wg, int* error)
{
  *error=0;
  hid_t h_orb_grp_read = H5Gopen(*h_file,"OrbsG",H5P_DEFAULT);
  if(h_orb_grp_read < 0) {
    *error=1;
    printf("esh5 error opening OrbsG group.\n");
    return;
  }
  char aname[64];
  sprintf(aname,"EigVal_kp%i",(*ik));
  herr_t ret=ret=H5LTread_dataset(h_orb_grp_read,aname,H5T_NATIVE_DOUBLE,et);

  sprintf(aname,"wg_kp%i",(*ik));
  ret=ret=H5LTread_dataset(h_orb_grp_read,aname,H5T_NATIVE_DOUBLE,wg);
  H5Gclose(h_orb_grp_read);
}

void F77_FUNC_(esh5_posthf_write_band,ESH5_POSTHF_WRITE_BAND)
  (hid_t* h_file, const char* dname, const int* dlength, const int* ik, const int* npwx, const int* ib, double* Psir, int* error)
{
  *error=0;
  char * hdname = ( char * ) malloc( (*dlength) + 1 ) ;
  memcpy( hdname , dname , *dlength ) ;
  hdname[*dlength] = '\0' ;

  hid_t h_orb_grp_write = H5Gopen(*h_file,hdname,H5P_DEFAULT);
  if(h_orb_grp_write < 0) {
    h_orb_grp_write=-1;
    *error=1;
    printf("esh5 error opening orbital group.\n");
    free(hdname);
    return;
  }
  free(hdname);

  char aname[64];
  hsize_t dims[2];
  dims[0] = (hsize_t)((*npwx));
  dims[1] = 2;
  sprintf(aname,"kp%i_b%i",*ik,*ib);
  herr_t ret=H5LTmake_dataset(h_orb_grp_write,aname,2,dims,H5T_NATIVE_DOUBLE,Psir);
  if(ret < 0){*error=1;printf("esh5 error writing orbital\n");}

  H5Gclose(h_orb_grp_write);
}

void F77_FUNC_(esh5_posthf_add_dm,ESH5_POSTHF_ADD_DM)
  ( const char* fname, const int* length, const char* gname, const int* glen, const int* n, const int* nk, const int* ns, double* DM, int* error)
{
  *error=0;
  char aname[64];
  hsize_t dims[2];
  int n2 = (*n)*(*n);
  dims[0] = (hsize_t)(n2);
  dims[1] = 2;

  H5Eget_auto (H5E_DEFAULT, &err_func, &client_data);
  H5Eset_auto (H5E_DEFAULT,NULL, NULL);

  char * hfname = ( char * ) malloc( (*length) + 1 ) ;
  memcpy( hfname , fname , *length ) ;
  hfname[*length] = '\0' ;

  hid_t f = H5Fopen(hfname,H5F_ACC_RDWR,H5P_DEFAULT);
  if(f<0) {
    printf("esh5 error opening orbital write file. %s\n",hfname);
    *error = 1;
    return;
  }

  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(f,gstr,H5P_DEFAULT);
  if(g < 0) g = H5Gcreate(f,gstr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

  const hsize_t dim1=1;
  herr_t ret=H5LTmake_dataset(g,"ndim",1,&dim1,H5T_NATIVE_INT,n);

  for( int is=0, isend=(*ns), inc=0; is<isend; is++) {
    for( int k=0, kend=(*nk); k<kend; k++, inc += 2*n2) {
      sprintf(aname,"s%i_kp%i",is,k);
      ret=H5LTmake_dataset(g,aname,2,dims,H5T_NATIVE_DOUBLE,DM+inc);
      if(ret < 0){*error=1;printf("esh5 error writing DM\n");}
    }
  }
  H5Gclose(g);
  H5Fclose(f);
  free(hfname);
}

// Expansion of an orbital set in terms of an underlying basis
// OMat[a,i,ik,ispin] = < basis(a,ik) | orbital(i,ik,ispin) >, ik=kpoint
void F77_FUNC_(esh5_posthf_write_orbmat,ESH5_POSTHF_WRITE_ORBMAT)
  (hid_t* h_file, const char* gname, const int* glen, const int* na, const int* ni, const int* nk, const int* ns, double* OMat, int* error)
{
  *error=0;
  char aname[64];
  hsize_t dims[2];
  int n2 = (*na)*(*ni);
  dims[0] = (hsize_t)(n2);
  dims[1] = 2;

  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(*h_file,gstr,H5P_DEFAULT);
  if(g < 0) g = H5Gcreate(*h_file,gstr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(g < 0){*error=1;printf("esh5 error making group in write_orbmat 0\n");return;}

  const hsize_t dim4=4;
  int data[4]; data[0]=*na; data[1]=*ni; data[2]=*nk; data[3]=*ns;
  herr_t ret=H5LTmake_dataset(g,"dims",1,&dim4,H5T_NATIVE_INT,data);
  if(ret < 0){*error=1;printf("esh5 error writing orbmat 1\n");return;}

  for( int is=0, nspin=(*ns), inc=0; is<nspin; is++) {
    for( int k=0, kend=(*nk); k<kend; k++, inc += 2*n2) {
      sprintf(aname,"s%i_kp%i",is,k);
      ret=H5LTmake_dataset(g,aname,2,dims,H5T_NATIVE_DOUBLE,OMat+inc);
      if(ret < 0){*error=1;printf("esh5 error writing orbmat 2\n");return;}
    }
  }
  H5Gclose(g);
}

void F77_FUNC_(esh5_posthf_read_orbmat_info,ESH5_POSTHF_READ_ORBMAT_INFO)
  (hid_t* h_file, const char* gname, const int* glen, int* na, int* ni, int* nk, int* ns, int* error)
{
  *error=0;
  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(*h_file,gstr,H5P_DEFAULT);
  if(g < 0) {*error=1;printf("esh5 error opening OMat group in orbmat_info\n");return;}

  int data[4];
  herr_t ret=H5LTread_dataset(g,"dims",H5T_NATIVE_INT,data);
  if(ret < 0){*error=1;printf("esh5 error reading dims in orbmat_info\n");}
  *na=data[0];
  *ni=data[1];
  *nk=data[2];
  *ns=data[3];
  H5Gclose(g);
}

void F77_FUNC_(esh5_posthf_read_orbmat,ESH5_POSTHF_READ_ORBMAT)
  (hid_t* h_file, const char* gname, const int* glen, const int* ik, const int* is, double* OMat, int* error)
{
  *error=0;
  char aname[64];

  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(*h_file,gstr,H5P_DEFAULT);
  if(g < 0) {*error=1;printf("esh5 error opening OMat group\n");return;}

  sprintf(aname,"s%i_kp%i",*is,*ik);
  herr_t ret=H5LTread_dataset(g,aname,H5T_NATIVE_DOUBLE,OMat);
  if(ret < 0){*error=1;printf("esh5 error reading OMat\n");}
  H5Gclose(g);
}

void F77_FUNC_(esh5_posthf_write_dm,ESH5_POSTHF_WRITE_DM)
  (hid_t* h_file, const char* gname, const int* glen, const int* n, const int* nk, const int* ns, double* DM, int* error)
{
  *error=0;
  char aname[64];
  hsize_t dims[2];
  int n2 = (*n)*(*n);
  dims[0] = (hsize_t)(n2);
  dims[1] = 2;

  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(*h_file,gstr,H5P_DEFAULT);
  if(g < 0) g = H5Gcreate(*h_file,gstr,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
  if(g < 0){*error=1;printf("esh5 error making group in write_DM 0\n");return;}

  const hsize_t dim1=1;
  herr_t ret=H5LTmake_dataset(g,"ndim",1,&dim1,H5T_NATIVE_INT,n);
  if(ret < 0){*error=1;printf("esh5 error writing DM 1\n");return;}

  for( int is=0, isend=(*ns), inc=0; is<isend; is++) {
    for( int k=0, kend=(*nk); k<kend; k++, inc += 2*n2) {
      sprintf(aname,"s%i_kp%i",is,k);
      ret=H5LTmake_dataset(g,aname,2,dims,H5T_NATIVE_DOUBLE,DM+inc);
      if(ret < 0){*error=1;printf("esh5 error writing DM 2\n");return;}
    }
  }
  H5Gclose(g);
}

void F77_FUNC_(esh5_posthf_read_dm,ESH5_POSTHF_READ_DM)
  (hid_t* h_file, const char* gname, const int* glen, const int* ik, const int* is, double* DM, int* error)
{
  *error=0;
  char aname[64];

  char gstr[256];
  memcpy( gstr , gname , *glen ) ;
  gstr[*glen] = '\0' ;
  hid_t g = H5Gopen(*h_file,gstr,H5P_DEFAULT);
  if(g < 0) {*error=1;printf("esh5 error opening DM group\n");return;}   

  sprintf(aname,"s%i_kp%i",*is,*ik);
  herr_t ret=H5LTread_dataset(g,aname,H5T_NATIVE_DOUBLE,DM);
  if(ret < 0){*error=1;printf("esh5 error reading DM\n");}
  H5Gclose(g);
}

void generate_SM(int nspin, const int nk, const int* norbK, const int norbmax, const int nelmax, 
  const int nfull, const int* fullocc, const int npart, const int* partocc,
  const int mixed, const int get_nk_occ, int* nkocc, const double* M, double* SM) 
{
  int npol = 1;
  if(nspin>1) npol=2;
  int* cpk = (int*) malloc(sizeof(int)*nk);
  int* epk = (int*) malloc(sizeof(int)*nk);
  for(int i=0; i<nk; i++) cpk[i]=0;
  for(int i=0; i<nk; i++) epk[i]=0;

  int norbtot = 0;
  for(int i=0; i<nk; i++) norbtot+=norbK[i];

  for(int i=0, ie=2*npol*(nfull+npart)*norbtot; i<ie; i++) SM[i]=0.0;

  // count electrons in each kpoint 
  for(int i=0; i<nfull; i++) {
    int indx = fullocc[i];
    int ik = indx/nelmax; // kpoint index 
    epk[ik]++;
  }
  for(int i=0; i<npart; i++) {
    int indx = partocc[i];
    int ik = indx/nelmax; // kpoint index 
    epk[ik]++;
  }
  if(get_nk_occ!=0) for(int i=0; i<nk; i++) nkocc[i]=epk[i]; 

  // SM( M(ik,iorb), N(ik,ib) ) = M(ik, ib, iorb) or delta(iorb,ib)
  for(int i=0; i<nfull; i++) {
    int indx = fullocc[i];
    int ik = indx/nelmax; // kpoint index 
    int ib = indx%nelmax;  // band index
    int M0 = 0;     // starting row of ik block
    for(int i=0; i<ik; i++) M0+=norbK[i];
    int N = cpk[ik];
    for(int ii=0; ii<ik; ii++) N += epk[ii];  
    if( mixed!=0 ) {
      // copy M into appropriate column in SM
      for(int iorb=0; iorb<norbK[ik]; iorb++) {
        int ip = ( M0 + iorb ) * ( nfull+npart ) + N;
        //int iq = ( ik*norb + iorb )*nelmax + ib; 
        int iq = ( ik*nelmax + ib )*norbmax + iorb; 
        SM[ 2*ip ]   = M[ 2*iq ];
        SM[ 2*ip+1 ] = M[ 2*iq+1 ];
      }   
    } else {
      // single contribution at iorb = ib
      int ip = ( M0 + ib ) * ( nfull+npart ) + N;
      SM[ 2*ip ]   = 1.0;
      SM[ 2*ip+1 ] = 0.0;
    }
    if(nspin > 1) { 
      int iq0 = nk*norbmax*nelmax;
      M0 += norbtot;  
      if( mixed!=0 ) {
        // copy M into appropriate column in SM
        for(int iorb=0; iorb<norbK[ik]; iorb++) {
          int ip = ( M0 + iorb ) * ( nfull+npart ) + N;
          //int iq = ( ik*norb + iorb )*nelmax + ib; 
          int iq = iq0 + ( ik*nelmax + ib )*norbmax + iorb;
          SM[ 2*ip ]   = M[ 2*iq ];
          SM[ 2*ip+1 ] = M[ 2*iq+1 ];
        }
      } else {
        // error !!!
      }
    }
    cpk[ik]++;
  }

  free(cpk);
  free(epk);
}

#endif
