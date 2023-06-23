#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#define SEED 100

#ifndef WRITE
#define WRITE
void write_pgm_image(MPI_Comm comm, const char *image, const unsigned char* buf, const int rank, const unsigned int buf_size, const unsigned int rest);
#endif

void set_parameters(unsigned long size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp){
  *rest = size % procs;
  *buffer_size = size/procs+( rank < *rest);
  unsigned int partial = *buffer_size * rank;
  *disp = rank < *rest ? partial : partial + *rest;
  //printf("hello from proc %d! Rest=%u, buffer_size=%u, disp=%lld.\n", rank, *rest, *buffer_size, *disp);
}

unsigned int write_header(const char* filename, const unsigned long rows, const unsigned long cols, const unsigned int maxval, const char* magic){
	FILE* f;
    unsigned long length=0;
    
    f=fopen(filename, "w");
    //printf("Size of magic: %lu. strlen(magic)=%lu\n",sizeof(magic),strlen(magic));
    //printf("%s è %lu caratteri + 1 a capo, %lu è %lu caratteri + 1 spazio, %lu è %lu caratteri+ 1 a capo,\
    %u è %d caratteri +1 a capo.\n",magic, strlen(magic),rows, (unsigned long)(ceil(log10(rows))),\
    cols, (unsigned long)(ceil(log10(cols))), maxval, (int)(ceil(log10(maxval))));
    length=fprintf(f, "%s\n%lu %lu\n%u\n",magic, rows, cols, maxval);
    //printf("Ho scritto su file %lu caratteri.\n", length);
    fclose(f);
    return length;
}

MPI_Comm initialize_procs(int argc, char* argv[], const int rows, int* my_procs, int* my_rank, int* rank_above, int* rank_below, int coords[2]){
  
  MPI_Comm cart_coord;
  MPI_Comm my_world=MPI_COMM_WORLD;
  
  const int periodic[2]={1,1};
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(my_world, my_procs);
  MPI_Comm_rank(my_world, my_rank); 
  MPI_Comm_size(my_world, my_procs); 

  if(*my_procs>=rows){
    //int color = rank < rows? 0:1;  // Determine the color based on the desired number of columns
    MPI_Group comm_world_group;
    MPI_Group my_world_group;
    int ranks[rows];
    for(int i=0; i<rows; i++){
      ranks[i]=i;
    }
    MPI_Comm_group(cart_coord, &comm_world_group);
    MPI_Group_incl(comm_world_group, rows, ranks,&my_world_group);
    MPI_Comm_create_group(MPI_COMM_WORLD, my_world_group, 0, &my_world);
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
  MPI_Cart_create(my_world, 2, dim, periodic, 1, &cart_coord);
  MPI_Cart_coords(cart_coord, *my_rank, 2, coords);
  MPI_Cart_shift(cart_coord, 0, 1, rank_above, rank_below);  
  }
  else{
    *rank_above=MPI_UNDEFINED;
    *rank_below=MPI_UNDEFINED;
    coords[0]=MPI_UNDEFINED;
    coords[1]=MPI_UNDEFINED;
  }


  
  return cart_coord;
}

void initialize_image(const unsigned int rows,const unsigned int cols, const unsigned int maxval, const char* image_name, int argc, char* argv[]){
  unsigned long size = rows*cols;
  int coords[2];
  int mpi_provided_thread_level;
  int my_rank, my_procs, rank_above, rank_below;
  

  //printf("Header \"%s\", lunghezza %lu", header, header_length);
  
  MPI_Offset disp;
  MPI_Comm my_world=initialize_procs(argc, argv, rows, &my_procs, &my_rank, &rank_above, &rank_below, coords);

  if(my_rank!=MPI_UNDEFINED){
    printf("hello from proc %d/%d. Coords[%d,%d]. Rank_above %d and rank_below %d\n", my_rank, my_procs, coords[0], coords[1], rank_above, rank_below);
    int seed=SEED*my_rank+my_procs;
  //printf("buffer_size for proc %lu: %lu\n", rank, buffer_size);
    unsigned int buffer_size, header_length, rest;
  

    if(my_rank==0){
      //printf("RAND MAX IS %d\n", RAND_MAX);
      header_length=write_header(image_name, rows, cols, maxval, "P5");
      //printf("Header ha lunghezza %u\n", header_length);
    }
    
    MPI_Bcast(&header_length,1,MPI_UNSIGNED,0,my_world);
    
    //printf("Ricevuto header length processo %d: %u\n", my_rank, header_length);
    set_parameters(size, my_procs, my_rank , &rest, &buffer_size, &disp);
    disp=disp+header_length;
    //printf("Proc %d. Header length: %u; rest: %u; buffer_size:%u; proc_offset: %lld\n", my_rank, header_length, rest, buffer_size, disp);
    
    unsigned char* local_buffer =(unsigned char*)malloc(buffer_size*sizeof(unsigned char) );
    srand(seed);
    
    for ( int i = 0; i < buffer_size; i++ ){
      local_buffer[i] = maxval*round((double)rand()/(double)RAND_MAX);
    }
    
    write_pgm_image(my_world, image_name, local_buffer, my_rank, buffer_size, disp);
    free(local_buffer);
  }
  
  MPI_Finalize();
}

/*
unsigned char * generate_image(struct Cell * grid, const  int size){
  unsigned char* cImage;

  cImage = calloc(size, sizeof(unsigned char) );
  for ( int i = 0; i < size; i++ ){
    unsigned char value = *((grid+i)->value);
    cImage[i] = value;
  }
  return cImage;
}
*/