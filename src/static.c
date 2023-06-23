#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

void read_pgm_image(unsigned char **image, const char *image_name,  unsigned int * maxval, unsigned int *rows, unsigned int *cols, unsigned int* rest, unsigned int buffer_size, int procs, int rank);
void set_parameters(unsigned long size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp);
MPI_Comm initialize_procs(int argc, char* argv[], const int rows, int* my_procs, int* my_rank, int* rank_above, int* rank_below, int coords[2]);
/*
unsigned long return_pos(offset, cols, ver_dir, hor_dir){
  return (offset%cols > 0 && hor_dir < 0) || (offset%cols +1 < cols && hor_dir > 0) ?
  offset+ver_dir*cols + hor_dir: (offset+ver_dir*cols+hor_dir)*(hor_dir<0)+*(offset-cols)*(hor_dir > 0);
}
*/

/*
unsigned char compute_neighbor(unsigned long matrix_pos, unsigned char* this_batch, unsigned char* nb_above, unsigned char* nb_below, unsigned long rows, unsigned long cols, MPI_Offset first, MPI_Offset last, unsigned int maxval){
    long offset=matrix_pos-first;
    long dist = matrix_pos-last;
    short int top_left = -1;
    short int bottom_right = +1;
    
    unsigned int cl_value = matrix_pos % cols > 0 ? this_batch[offset-1]: this_batch[offset + cols - 1];
    unsigned int cr_value = (matrix_pos + 1) % cols > 0 ? this_batch[offset + 1]: this_batch[offset - cols +1];
    unsigned int tc_value = (long)(offset - cols) > 0 ? this_batch[offset-cols] : nb_above[matrix_pos%cols];

    unsigned int tl_index=return_pos(offset, cols, top_left, top_left);
    unsigned int tl_value =(long) (offset - cols) > 1 ? this_batch[tl_index]: nb_above[tl_index];
    unsigned int tr_index=return_pos(offset, cols, top_left, bottom_right);
    unsigned int tr_value = (long)(offset - cols) > - 1 ? this_batch[tr_index]: nb_above[tr_index];
    
    unsigned int bc_value = (long)(dist + cols) < 0 ? this_batch[offset+cols] : nb_below[matrix_pos%cols];
    unsigned int bl_index=return_pos(offset, cols, bottom_right, top_left);
    unsigned int bl_value = (long)(dist + cols)< 1 ? this_batch[bl_index]: nb_below[bl_index];
    unsigned int br_index=return_pos(offset, cols, bottom_right, bottom_right);
    unsigned int br_value = (long)(dist + cols) < -1 ? this_batch[br_index]: nb_below[br_index];

    unsigned int sum = tl_value+tc_value+tr_value+cl_value+cr_value+bl_value+bc_value+br_value;
    printf("Matrix pos %lu, vicini: %u %u %u %u %u %u %u %u\n", matrix_pos,tl_value, tc_value, tr_value, cl_value, cr_value, bl_value, bc_value, br_value);
    return sum >= maxval*2 && sum <= maxval * 3? (unsigned char) maxval : (unsigned char)0;
} 
*/

void evolve_static(unsigned char* ptr, unsigned int rows, unsigned int cols, unsigned int* maxval, const char* filename, int argc, char* argv[]){
  //printf("START:\n");
  //for(int grid_value=0; grid_value<size; grid_index++)
    //printf("\tPosizione %d: indirizz  o %p valore %d\n", grid_index, (grid+grid_index)->value, *((grid+grid_index)->value));

  MPI_Offset disp;
  MPI_Request req;
  MPI_Status status;
  unsigned long i=disp;
  unsigned long size=rows*cols;
  unsigned int rest, buffer_size, color_depth, mpi_provided_thread_level;
  char MagicN[3];
  int coords[2];
  int procs, rank, proc_above, proc_below;
  MPI_Comm cart_comm=initialize_procs(argc, argv, rows, &procs, &rank, &proc_above, &proc_below, coords);
  
  read_pgm_image(&ptr, filename, maxval, &rows, &cols, &rest, buffer_size, procs, rank);
  //printf("%lu, %lu\n", *rows, *cols);
  


  //size= *rows * *cols * (1 + *maxval>255);

  //printf("Processo %d, coordinate [%d,%d]\n", cart_rank, coords[0], coords[1]);
  
  set_parameters(size, procs, rank, &rest, &buffer_size, &disp);
  /*
  unsigned char * evo = malloc(buffer_size);
  unsigned long end_buffer=disp+buffer_size-1;
  unsigned char nb_above [*cols];
  unsigned char nb_below [*cols];
  //printf("Proc. %d start:%lld end:%lu\n", cart_rank, disp, end_buffer);
  MPI_Isend(ptr, *cols, MPI_UNSIGNED_CHAR, proc_above, 1, MPI_COMM_WORLD, &req);
  MPI_Recv(nb_below, *cols, MPI_UNSIGNED_CHAR, proc_below, 1 , MPI_COMM_WORLD, &status);
  MPI_Isend(ptr + (buffer_size - *cols), *cols, MPI_UNSIGNED_CHAR, proc_below, 0, MPI_COMM_WORLD, &req);
  MPI_Recv(nb_above, *cols, MPI_UNSIGNED_CHAR, proc_above, 0 , MPI_COMM_WORLD, &status);

  unsigned long j=0;
  i=disp;
  //printf("%lu\n", sizeof(nb_above));

  while(i<=end_buffer){
    evo[j]=compute_neighbor(i, ptr, nb_above, nb_below, *rows, *cols, disp, end_buffer, *maxval);
    i++;
    j++;
  }
  ptr=evo;
  printf("\n\n");
  for (int i=0; i< buffer_size; i++){
    if(evo[i]!= ptr[i])
    printf("MAIALE\n");
  }
  */
  MPI_Finalize();
    
  return;
}
