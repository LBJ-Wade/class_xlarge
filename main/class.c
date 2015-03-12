/** @file class.c
 * Julien Lesgourgues, 17.04.2011
 */

#include "class.h"

#ifdef _DAM_MPI
#include <mpi.h>

struct transfers mpi_join_transfers(struct transfers *ptr,struct perturbs *ppt)
{
  int ii;
  struct transfers ptr_glb;

  mpi_printf("Gathering transfer functions from all nodes\n");

  //  memcpy(&ptr_glb,ptr,sizeof(struct transfers));
  //Useless now
  ptr_glb.has_nz_file=0;
  ptr_glb.has_nz_analytic=0;
  sprintf(ptr_glb.nz_file_name,"shit");
  ptr_glb.nz_size=0;
  ptr_glb.nz_z=NULL;
  ptr_glb.nz_nz=NULL;
  ptr_glb.nz_ddnz=NULL;
  ptr_glb.has_bz_file=0;
  sprintf(ptr_glb.bz_file_name,"shit");
  ptr_glb.bz_size=0;
  ptr_glb.bz_z=NULL;
  ptr_glb.bz_bz=NULL;
  ptr_glb.bz_ddbz=NULL;
  ptr_glb.has_sz_file=0;
  sprintf(ptr_glb.sz_file_name,"shit");
  ptr_glb.sz_size=0;
  ptr_glb.sz_z=NULL;
  ptr_glb.sz_sz=NULL;
  ptr_glb.sz_ddsz=NULL;
  ptr_glb.has_ez_file=0;
  sprintf(ptr_glb.ez_file_name,"shit");
  ptr_glb.ez_size=0;
  ptr_glb.ez_z=NULL;
  ptr_glb.ez_ez=NULL;
  ptr_glb.ez_ddez=NULL;

  ptr_glb.lcmb_rescale=ptr->lcmb_rescale;
  ptr_glb.lcmb_tilt=ptr->lcmb_tilt;
  ptr_glb.lcmb_pivot=ptr->lcmb_pivot;
  ptr_glb.bias=ptr->bias;
  ptr_glb.s_bias=ptr->s_bias;
  ptr_glb.has_cls=ptr->has_cls;
  ptr_glb.md_size=ptr->md_size;
  ptr_glb.tt_size=malloc(ptr_glb.md_size*sizeof(int));
  ptr_glb.l_size=malloc(ptr_glb.md_size*sizeof(int));
  for(ii=0;ii<ptr->md_size;ii++)
    ptr_glb.l_size[ii]=ptr->l_size[ii];
  ptr_glb.l_size_max=ptr->l_size_max;
  ptr_glb.l=malloc(ptr_glb.l_size_max*sizeof(int));
  for(ii=0;ii<ptr->l_size_max;ii++)
    ptr_glb.l[ii]=ptr->l[ii];
  ptr_glb.angular_rescaling=ptr->angular_rescaling;
  ptr_glb.q_size=ptr->q_size;
  ptr_glb.q=malloc(ptr_glb.q_size*sizeof(double));
  for(ii=0;ii<ptr_glb.q_size;ii++)
    ptr_glb.q[ii]=ptr->q[ii];
  ptr_glb.k=malloc(ptr_glb.md_size*sizeof(double *));
  for(ii=0;ii<ptr_glb.md_size;ii++) {
    int jj;
    ptr_glb.k[ii]=malloc(ptr_glb.q_size*sizeof(double));
    for(jj=0;jj<ptr_glb.q_size;jj++)
      ptr_glb.k[ii][jj]=ptr->k[ii][jj];
  }
  ptr_glb.index_q_flat_approximation=ptr->index_q_flat_approximation;
  ptr_glb.initialise_HIS_cache=ptr->initialise_HIS_cache;
  ptr_glb.transfer_verbose=ptr->transfer_verbose;
  sprintf(ptr_glb.error_message,"Nothing");
  printf("Node %d done copying first bunch\n",Mpi_this_node);

  //Compute total number of selection functions
  int *selection_num_allnodes=malloc(Mpi_n_nodes*sizeof(int));
  int *transfer_nelements=malloc(Mpi_n_nodes*sizeof(int));
  int *transfer_displ=malloc(Mpi_n_nodes*sizeof(int));
  MPI_Allgather(&(ptr->selection_num),1,MPI_INT,
		selection_num_allnodes,1,MPI_INT,MPI_COMM_WORLD);
  ptr_glb.selection_num=0;
  for(ii=0;ii<Mpi_n_nodes;ii++) {
#ifdef _DAM_DEBUG
    printf("Node %d: Node %d has %d bins\n",Mpi_this_node,ii,selection_num_allnodes[ii]);
#endif //_DAM_DEBUG
    ptr_glb.selection_num+=selection_num_allnodes[ii];
  }

  //Reorganize transfer indices
#ifdef _DAM_DEBUG
  printf("Node %d: Reorganising transfer indices\n",Mpi_this_node);
#else //_DAM_DEBUG
  mpi_printf("Reorganising transfer indices\n");
#endif //_DAM_DEBUG
  int index_tt=0,index_tt_common;

  class_define_index(ptr_glb.index_tt_t2,ppt->has_cl_cmb_temperature,index_tt,1);
  class_define_index(ptr_glb.index_tt_e,ppt->has_cl_cmb_polarization,index_tt,1);

  index_tt_common=index_tt;

  if(ppt->has_scalars == _TRUE_) {
    index_tt=index_tt_common;

    class_define_index(ptr_glb.index_tt_t0,ppt->has_cl_cmb_temperature,index_tt,1);
    class_define_index(ptr_glb.index_tt_t1,ppt->has_cl_cmb_temperature,index_tt,1);
    class_define_index(ptr_glb.index_tt_lcmb,ppt->has_cl_cmb_lensing_potential,
		       index_tt,1);
    class_define_index(ptr_glb.index_tt_density,ppt->has_nc_density,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_rsd,ppt->has_nc_rsd1,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_d0,ppt->has_nc_rsd2,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_d1,ppt->has_nc_rsd3,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_lens,ppt->has_nc_lens,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_g1,ppt->has_nc_gr1,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_g2,ppt->has_nc_gr2,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_g3,ppt->has_nc_gr3,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_g4,ppt->has_nc_gr4,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_nc_g5,ppt->has_nc_gr5,index_tt,
		       ptr_glb.selection_num);
    class_define_index(ptr_glb.index_tt_lensing,ppt->has_cl_lensing_potential,index_tt,
		       ptr_glb.selection_num);

    ptr_glb.tt_size[ppt->index_md_scalars]=index_tt;
  }

  if (ppt->has_vectors == _TRUE_) {
    index_tt=index_tt_common;

    class_define_index(ptr_glb.index_tt_t1,ppt->has_cl_cmb_temperature, index_tt,1);
    class_define_index(ptr_glb.index_tt_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr_glb.tt_size[ppt->index_md_vectors]=index_tt;
  }

  if (ppt->has_tensors == _TRUE_) {
    index_tt=index_tt_common;

    class_define_index(ptr_glb.index_tt_b, ppt->has_cl_cmb_polarization,index_tt,1);

    ptr_glb.tt_size[ppt->index_md_tensors]=index_tt;
  }

  //Compute l number for each transfer
  ptr_glb.l_size_tt=malloc(ptr_glb.md_size*sizeof(int *));
  for(ii=0;ii<ptr_glb.md_size;ii++) {
    ptr_glb.l_size_tt[ii]=malloc(ptr_glb.tt_size[ii]*sizeof(int));
  }

  if(ppt->has_scalars==_TRUE_) {
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_t0]=
	ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_t0];
      ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_t1]=
	ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_t1];
      ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_t2]=
	ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_t2];
    }
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_e]=
	ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_e];
    }
    if(ppt->has_cl_cmb_lensing_potential==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_lcmb]=
	ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_lcmb];
    }
    if(ppt->has_nc_density==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_density+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_density];
      }
    }
    if(ppt->has_nc_rsd1==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_rsd+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_rsd];
      }
    }
    if(ppt->has_nc_rsd2==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_d0+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_d0];
      }
    }
    if(ppt->has_nc_rsd3==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_d1+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_d1];
      }
    }
    if(ppt->has_nc_lens==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_lens+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_lens];
      }
    }
    if(ppt->has_nc_gr1==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_g1+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_g1];
      }
    }
    if(ppt->has_nc_gr2==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_g2+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_g2];
      }
    }
    if(ppt->has_nc_gr3==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_g3+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_g3];
      }
    }
    if(ppt->has_nc_gr4==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_g4+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_g4];
      }
    }
    if(ppt->has_nc_gr5==_TRUE_) {
      for(ii=0;ii<ptr_glb.selection_num;ii++) {
	ptr_glb.l_size_tt[ppt->index_md_scalars][ptr_glb.index_tt_nc_g5+ii]=
	  ptr->l_size_tt[ppt->index_md_scalars][ptr->index_tt_nc_g5];
      }
    }
  }

  if(ppt->has_vectors==_TRUE_) {
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_vectors][ptr_glb.index_tt_t1]=
	ptr->l_size_tt[ppt->index_md_vectors][ptr->index_tt_t1];
      ptr_glb.l_size_tt[ppt->index_md_vectors][ptr_glb.index_tt_t2]=
	ptr->l_size_tt[ppt->index_md_vectors][ptr->index_tt_t2];
    }
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_vectors][ptr_glb.index_tt_e]=
	ptr->l_size_tt[ppt->index_md_vectors][ptr->index_tt_e];
      ptr_glb.l_size_tt[ppt->index_md_vectors][ptr_glb.index_tt_b]=
	ptr->l_size_tt[ppt->index_md_vectors][ptr->index_tt_b];
    }
  }

  if(ppt->has_tensors==_TRUE_) {
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_tensors][ptr_glb.index_tt_t2]=
	ptr->l_size_tt[ppt->index_md_tensors][ptr->index_tt_t2];
    }
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      ptr_glb.l_size_tt[ppt->index_md_tensors][ptr_glb.index_tt_b]=
	ptr->l_size_tt[ppt->index_md_tensors][ptr->index_tt_b];
      ptr_glb.l_size_tt[ppt->index_md_tensors][ptr_glb.index_tt_e]=
	ptr->l_size_tt[ppt->index_md_tensors][ptr->index_tt_e];
    }
  }

  //Gather transfers
  ptr_glb.transfer=malloc(ptr_glb.md_size*sizeof(double *));
  int ran_out=0;
  if(ptr_glb.transfer==NULL) {
    ran_out=1;
    mpi_abort(1,"Some processes ran out of memory\n");
  }
  for(ii=0;ii<ptr_glb.md_size;ii++) {
    int transfer_size=ptr_glb.tt_size[ii]*ptr_glb.l_size[ii]*ptr_glb.q_size;
    ptr_glb.transfer[ii]=malloc(transfer_size*sizeof(double));
    if(ptr_glb.transfer[ii]==NULL){
      ran_out=1;
      printf("Node %d: out of memory %d %d %d\n",Mpi_this_node,ii,ptr_glb.tt_size[ii],transfer_size);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(ran_out)
    mpi_abort(1,"Some processes ran out of memory\n");

#ifdef _DAM_DEBUG
  printf("Node %d: Copying CMB transfers\n",Mpi_this_node);
#else //_DAM_DEBUG
  mpi_printf("Copying CMB transfers\n");
#endif //_DAM_DEBUG
  //First just copy the CMB ones
  if(ppt->has_scalars==_TRUE_) {
    int index_md=ppt->index_md_scalars;
    int size_this_transfer=ptr_glb.l_size[index_md]*ptr->q_size;
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t0*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t0*size_this_transfer]),
	     size_this_transfer);
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t1*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t1*size_this_transfer]),
	     size_this_transfer);
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t2*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t2*size_this_transfer]),
	     size_this_transfer);
    }
    if(ppt->has_cl_cmb_lensing_potential==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_lcmb*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_lcmb*size_this_transfer]),
	     size_this_transfer);
    }     
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_e*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_e*size_this_transfer]),
	     size_this_transfer);
    }     
  }
  if(ppt->has_vectors==_TRUE_) {
    int index_md=ppt->index_md_vectors;
    int size_this_transfer=ptr_glb.l_size[index_md]*ptr->q_size;
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t1*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t1*size_this_transfer]),
	     size_this_transfer);
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t2*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t2*size_this_transfer]),
	     size_this_transfer);
    }
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_e*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_e*size_this_transfer]),
	     size_this_transfer);
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_b*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_b*size_this_transfer]),
	     size_this_transfer);
    }     
  }
  if(ppt->has_tensors==_TRUE_) {
    int index_md=ppt->index_md_tensors;
    int size_this_transfer=ptr_glb.l_size[index_md]*ptr->q_size;
    if(ppt->has_cl_cmb_temperature==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_t2*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_t2*size_this_transfer]),
	     size_this_transfer);
    }
    if(ppt->has_cl_cmb_polarization==_TRUE_) {
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_e*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_e*size_this_transfer]),
	     size_this_transfer);
      memcpy(&(ptr_glb.transfer[index_md][ptr_glb.index_tt_b*size_this_transfer]),
	     &(ptr->transfer[index_md][ptr->index_tt_b*size_this_transfer]),
	     size_this_transfer);
    }     
  }

  //Now gather the number counts
#ifdef _DAM_DEBUG
  printf("Node %d: Gathering number counts transfers\n",Mpi_this_node);
#else //_DAM_DEBUG
  mpi_printf("Gathering number counts transfers\n");
#endif //_DAM_DEBUG
  if(ppt->has_scalars==_TRUE_) {
    int index_md=ppt->index_md_scalars;
    int size_this_transfer=ptr_glb.l_size[index_md]*ptr->q_size;

    transfer_nelements[0]=selection_num_allnodes[0]*size_this_transfer;
    transfer_displ[0]=0;
    for(ii=1;ii<Mpi_n_nodes;ii++) {
      transfer_nelements[ii]=selection_num_allnodes[ii]*size_this_transfer;
      transfer_displ[ii]=transfer_displ[ii-1]+transfer_nelements[ii-1];
    }

    if(ppt->has_nc_density==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_density*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_density*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_rsd1==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_rsd*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_rsd*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_rsd2==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_d0*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_d0*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_rsd3==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_d1*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_d1*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_lens==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_lens*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_lens*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_gr1==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_g1*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_g1*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_gr2==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_g2*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_g2*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_gr3==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_g3*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_g3*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_gr4==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_g4*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_g4*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
    if(ppt->has_nc_gr5==_TRUE_) {
      MPI_Allgatherv(&(ptr->transfer[index_md][ptr->index_tt_nc_g5*size_this_transfer]),
		     ptr->selection_num*size_this_transfer,MPI_DOUBLE,
		     &(ptr_glb.transfer[index_md][ptr_glb.index_tt_nc_g5*size_this_transfer]),
		     transfer_nelements,transfer_displ,MPI_DOUBLE,MPI_COMM_WORLD);
    }
  }

  free(selection_num_allnodes);
  free(transfer_nelements);
  free(transfer_displ);
#ifdef _DAM_DEBUG
  printf("Node %d: Done gathering transfers\n",Mpi_this_node);
#endif //_DAM_DEBUG

  return ptr_glb;
}
#endif //_DAM_MPI

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */
  
#ifdef _DAM_MPI
  mpi_init(&argc,&argv);
#endif //_DAM_MPI

  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    mpi_abort(1,"\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    mpi_abort(1,"\n\nError running background_init \n=>%s\n",ba.error_message);
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    mpi_abort(1,"\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    mpi_abort(1,"\n\nError in perturb_init \n=>%s\n",pt.error_message);
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    mpi_abort(1,"\n\nError in primordial_init \n=>%s\n",pm.error_message);
  }

  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    mpi_abort(1,"\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
  }

  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    mpi_abort(1,"\n\nError in transfer_init \n=>%s\n",tr.error_message);
  }

#ifdef _DAM_MPI
#ifdef _DAM_DEBUG
  printf("Node %d finished  transfers\n",Mpi_this_node);
#endif //_DAM_DEBUG
  struct transfers tr_glb=mpi_join_transfers(&tr,&pt);
  if(transfer_free(&tr) == _FAILURE_)
    mpi_abort(1,"\n\nError in transfer_free \n=>%s\n",tr.error_message);

  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr_glb,&sp) == _FAILURE_) {
    mpi_abort(1,"\n\nError in spectra_init \n=>%s\n",sp.error_message);
  }
  
  if(Mpi_this_node==0) {
    if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
      mpi_abort(1,"\n\nError in lensing_init \n=>%s\n",le.error_message);
    }
  }

  if (output_init(&ba,&th,&pt,&pm,&tr_glb,&sp,&nl,&le,&op) == _FAILURE_) {
    mpi_abort(1,"\n\nError in output_init \n=>%s\n",op.error_message);
  }
#else //_DAM_MPI
  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
    mpi_abort(1,"\n\nError in spectra_init \n=>%s\n",sp.error_message);
  }
  
  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    mpi_abort(1,"\n\nError in lensing_init \n=>%s\n",le.error_message);
  }
  
  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
    mpi_abort(1,"\n\nError in output_init \n=>%s\n",op.error_message);
  }
#endif //_DAM_MPI
  
  /****** all calculations done, now free the structures ******/
  
#ifdef _DAM_MPI
  if(Mpi_this_node==0) {
    if (lensing_free(&le) == _FAILURE_) {
      mpi_abort(1,"\n\nError in lensing_free \n=>%s\n",le.error_message);
    }
  }
    
  if (spectra_free(&sp) == _FAILURE_) {
    mpi_abort(1,"\n\nError in spectra_free \n=>%s\n",sp.error_message);
  }

  if(transfer_free(&tr_glb) == _FAILURE_) {
    mpi_abort(1,"\n\nError in transfer_free \n=>%s\n",tr_glb.error_message);
  }
#else //_DAM_MPI
  if (lensing_free(&le) == _FAILURE_) {
    mpi_abort(1,"\n\nError in lensing_free \n=>%s\n",le.error_message);
  }
  
  if (spectra_free(&sp) == _FAILURE_) {
    mpi_abort(1,"\n\nError in spectra_free \n=>%s\n",sp.error_message);
  }

  if (transfer_free(&tr) == _FAILURE_) {
    mpi_abort(1,"\n\nError in transfer_free \n=>%s\n",tr.error_message);
  }
#endif //_DAM_MPI

  if (nonlinear_free(&nl) == _FAILURE_) {
    mpi_abort(1,"\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
  }

  if (primordial_free(&pm) == _FAILURE_) {
    mpi_abort(1,"\n\nError in primordial_free \n=>%s\n",pm.error_message);
  }

  if (perturb_free(&pt) == _FAILURE_) {
    mpi_abort(1,"\n\nError in perturb_free \n=>%s\n",pt.error_message);
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    mpi_abort(1,"\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
  }

  if (background_free(&ba) == _FAILURE_) {
    mpi_abort(1,"\n\nError in background_free \n=>%s\n",ba.error_message);
  }

#ifdef _DAM_MPI
  mpi_end();
#endif //_DAM_MPI

  return _SUCCESS_;
}
