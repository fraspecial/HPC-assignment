#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#define SEED 435423

#ifndef WRITE
#define WRITE
void write_pgm_image(const MPI_Comm* comm, const char *image_name, unsigned char * buf, const unsigned int procs, const unsigned int rank, const unsigned int rows, const unsigned int cols, const unsigned int maxval);
#endif
  
void set_parameters(unsigned int rows, unsigned int cols, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp){
  *rest = rows % procs;
  *buffer_size = cols * (rows/procs+( rank < *rest));
  unsigned int partial = *buffer_size * rank;
  *disp = rank < *rest ? partial : partial + cols*(*rest);
}

unsigned int write_header(const char* filename, const unsigned long rows, const unsigned long cols, const unsigned int maxval, const char* magic){
	FILE* f;
    unsigned long length=0;
    
    f=fopen(filename, "w");
    length=fprintf(f, "%s\n%lu %lu\n%u\n",magic, cols, rows, maxval);
    fclose(f);
    return length;
}

MPI_Comm initialize_procs(const int rows, int* my_procs, int* my_rank, int* rank_above, int* rank_below, int coords[2]){
  
  //MPI_Comm cart_coord;
  MPI_Comm my_world=MPI_COMM_WORLD;
  
  const int periodic[2]={1,0};
  
  MPI_Comm_size(my_world, my_procs);
  MPI_Comm_rank(my_world, my_rank); 
  
  if(*my_procs>=rows){
    
    MPI_Group comm_world_group;
    MPI_Group my_world_group;

    int ranks[rows];
    for(int i=0; i<rows; i++){
      ranks[i]=i;
    }

    MPI_Comm_group(my_world, &comm_world_group);
    MPI_Group_incl(comm_world_group, rows, ranks, &my_world_group);
    MPI_Comm_create_group(my_world, my_world_group, 0, &my_world);
    
    if(my_world!=MPI_COMM_NULL){
      MPI_Comm_rank(my_world, my_rank); 
      MPI_Comm_size(my_world, my_procs);
    }
    else{
      *my_rank=MPI_UNDEFINED;
      *my_procs=MPI_UNDEFINED;
    }
    MPI_Group_free(&comm_world_group);
    MPI_Group_free(&my_world_group);    
  }
  int dim[2]={*my_procs, 1};
  if(my_world!=MPI_COMM_NULL){
  MPI_Cart_create(my_world, 2, dim, periodic, 1, &my_world);
  MPI_Cart_coords(my_world, *my_rank, 2, coords);
  MPI_Cart_shift(my_world, 0, 1, rank_above, rank_below);
  }
  else{
    *rank_above=MPI_UNDEFINED;
    *rank_below=MPI_UNDEFINED;
    coords[0]=MPI_UNDEFINED;
    coords[1]=MPI_UNDEFINED;
  }
  return my_world;
}

void initialize_image(const unsigned int rows,const unsigned int cols, const unsigned int* maxval, const char* image_name, int* argc, char** argv[]){
  int coords[2];
  int mpi_provided_thread_level, my_rank, my_procs, rank_above, rank_below, rest, buffer_size;
  
  MPI_Offset disp;
  MPI_Init(argc, argv);
  MPI_Comm my_world=initialize_procs(rows, &my_procs, &my_rank, &rank_above, &rank_below, coords);
  unsigned int seed = SEED *my_rank + my_procs;
  if(my_rank!=MPI_UNDEFINED){
    srand(seed);
    rest = rows % my_procs;
    buffer_size = cols * (rows/my_procs+( my_rank < rest));
    unsigned char* local_buffer =(unsigned char*)malloc(buffer_size*sizeof(unsigned char) );
  

    for ( int i = 0; i < buffer_size; i++ )
      local_buffer[i] = *maxval*round((double)rand()/(double)RAND_MAX);
   

    write_pgm_image(&my_world, image_name, local_buffer, my_procs, my_rank, rows, cols, *maxval);
    free(local_buffer);
  }
  
  MPI_Finalize();
}
