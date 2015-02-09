#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _DAM_MPI
#include <mpi.h>
#endif //_DAM_MPI
#include <stdarg.h>
#include "common.h"
#include "transfer.h"
#include "perturbations.h"

void mpi_init(int *p_argc,char ***p_argv)
{
#ifdef _DAM_MPI
  MPI_Init(p_argc,p_argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&(Mpi_this_node));
  MPI_Comm_size(MPI_COMM_WORLD,&(Mpi_n_nodes));
#else //_DAM_MPI
  Mpi_this_node=0;
  Mpi_n_nodes=1;
#endif //_DAM_MPI
  Mpi_bins_allocated=0;
  Mpi_cross_allocated=0;
}

void mpi_printf(const char *fmt, ...)
{
  if(Mpi_this_node==0) {
    va_list argp;
    
    va_start(argp,fmt);
    vfprintf(stdout,fmt,argp);
    fflush(stdout);
    va_end(argp);
  }
}

void mpi_abort(const int errtyp,const char *fmt, ...)
{
    va_list argp;
    
    va_start(argp,fmt);
    vfprintf(stderr,fmt,argp);
    va_end(argp);

#ifdef _DAM_MPI
    MPI_Abort(MPI_COMM_WORLD,errtyp);
#else //_DAM_MPI
    exit(errtyp);
#endif //_DAM_MPI
}

void mpi_share_bins(int n_bins)
{
#ifdef _DAM_OPTIMIZED_ORDERING
  int ii,nbins_here;

  Mpi_nbins_total=n_bins;

  if(n_bins<Mpi_n_nodes) {
    mpi_abort(1,"MPI: too many nodes (%d) for too few bins (%d)\n",
	      Mpi_n_nodes,n_bins);
  }

  Mpi_bin0=Mpi_this_node;
  nbins_here=1;
  for(ii=0;ii<Mpi_nbins_total;ii++) {
    if((ii>Mpi_bin0)&&((ii-Mpi_bin0)%Mpi_n_nodes==0)) {
      nbins_here++;
    }
  }

  Mpi_nbins_here=nbins_here;
  Mpi_bin_ids_here=malloc(Mpi_nbins_here*sizeof(int));
  Mpi_bin_ids_here[0]=Mpi_bin0;
  nbins_here=1;
  for(ii=0;ii<Mpi_nbins_total;ii++) {
    if((ii>Mpi_bin0)&&((ii-Mpi_bin0)%Mpi_n_nodes==0)) {
      Mpi_bin_ids_here[nbins_here]=ii;
      nbins_here++;
    }
  }

  Mpi_bin_ids_total=malloc(Mpi_nbins_total*sizeof(int));
  int index=0;
  for(ii=0;ii<Mpi_n_nodes;ii++) {
    int bin_id=ii;
    while(bin_id<Mpi_nbins_total) {
      Mpi_bin_ids_total[index]=bin_id;
      bin_id+=Mpi_n_nodes;
      index++;
    }
  }
#else //_DAM_OPTIMIZED_ORDERING
  int ii,nb_extra,nb_per;

  Mpi_nbins_total=n_bins;

  if(n_bins<Mpi_n_nodes) {
    mpi_abort(1,"MPI: too many nodes (%d) for too few bins (%d)\n",
	      Mpi_n_nodes,n_bins);
  }
  nb_extra=n_bins%Mpi_n_nodes;
  nb_per=n_bins/Mpi_n_nodes;

  if(nb_extra==0) {
    Mpi_nbins_here=nb_per;
    Mpi_bin0=Mpi_this_node*Mpi_nbins_here;
  }
  else {
    if(Mpi_this_node<nb_extra) {
      Mpi_nbins_here=nb_per+1;
      Mpi_bin0=Mpi_this_node*Mpi_nbins_here;
    }
    else {
      Mpi_nbins_here=nb_per;
      Mpi_bin0=(Mpi_this_node-nb_extra)*nb_per+nb_extra*(nb_per+1);
    }
  }

  Mpi_bin_ids_here=malloc(Mpi_nbins_here*sizeof(int));
  for(ii=0;ii<Mpi_nbins_here;ii++)
    Mpi_bin_ids_here[ii]=Mpi_bin0+ii;

  Mpi_bin_ids_total=malloc(Mpi_nbins_total*sizeof(int));
  for(ii=0;ii<Mpi_nbins_total;ii++)
    Mpi_bin_ids_total[ii]=ii;
#endif //_OPTIMIZED_ORDERING

  Mpi_bins_allocated=1;
}

void mpi_end(void)
{
  if(Mpi_bins_allocated==1) {
    free(Mpi_bin_ids_here);
    free(Mpi_bin_ids_total);
  }

  if(Mpi_cross_allocated==1) {
    free(Mpi_i1_cross);
    free(Mpi_i2_cross);
  }
#ifdef _DAM_MPI
  MPI_Finalize();
#endif //_DAM_MPI
}

void mpi_distribute_spectra(int n_bins,int non_diag)
{
  int ii;
  int n_cross_total,i_cross,n_cross_extra,n_cross_per,cross0;
  int *i1_cross_total;
  int *i2_cross_total;

  n_cross_total=0;
  for(ii=0;ii<n_bins;ii++) {
    int jj;
    for(jj=ii;jj<=MIN(ii+non_diag,n_bins-1);jj++)
      n_cross_total++;
  }

  n_cross_extra=n_cross_total%Mpi_n_nodes;
  n_cross_per=n_cross_total/Mpi_n_nodes;

  if(n_cross_extra==0) {
    Mpi_cross_here=n_cross_per;
    cross0=Mpi_this_node*Mpi_cross_here;
  }
  else {
    if(Mpi_this_node<n_cross_extra) {
      Mpi_cross_here=n_cross_per+1;
      cross0=Mpi_this_node*Mpi_cross_here;
    }
    else {
      Mpi_cross_here=n_cross_per;
      cross0=(Mpi_this_node-n_cross_extra)*n_cross_per+
	n_cross_extra*(n_cross_per+1);
    }
  }
  Mpi_i1_cross=malloc(Mpi_cross_here*sizeof(int));
  Mpi_i2_cross=malloc(Mpi_cross_here*sizeof(int));
  Mpi_cross_allocated=1;

  i1_cross_total=malloc(n_cross_total*sizeof(int));
  i2_cross_total=malloc(n_cross_total*sizeof(int));

  i_cross=0;
  for(ii=0;ii<n_bins;ii++) {
    int jj;
    for(jj=ii;jj<=MIN(ii+non_diag,n_bins-1);jj++) {
      i1_cross_total[i_cross]=ii;
      i2_cross_total[i_cross]=jj;
      i_cross++;
    }
  }

  for(ii=0;ii<Mpi_cross_here;ii++) {
    Mpi_i1_cross[ii]=i1_cross_total[cross0+ii];
    Mpi_i2_cross[ii]=i2_cross_total[cross0+ii];
  }

  free(i1_cross_total);
  free(i2_cross_total);
}

void mpi_print_bins(void)
{
  int ii;

  mpi_printf("Bin order is:\n");
  for(ii=0;ii<Mpi_nbins_total;ii++)
    mpi_printf("%d %d\n",ii,Mpi_bin_ids_total[ii]);
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Node %d will compute the transfer function for bins: ",Mpi_this_node);
  for(ii=0;ii<Mpi_nbins_here;ii++)
    printf("%d, ",Mpi_bin_ids_here[ii]);
  printf("\n");

  printf("Node %d will compute the following cross-spectra: ",Mpi_this_node);
  for(ii=0;ii<Mpi_cross_here;ii++) {
    int ii1=Mpi_bin_ids_total[Mpi_i1_cross[ii]];
    int ii2=Mpi_bin_ids_total[Mpi_i2_cross[ii]];
    int id1=MIN(ii1,ii2);
    int id2=MAX(ii1,ii2);

    printf(" %dx%d (%dx%d), ",Mpi_i1_cross[ii],Mpi_i2_cross[ii],id1,id2);
  }
  printf("\n");
}
