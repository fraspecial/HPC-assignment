#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define SEED 100

unsigned char * initialise_image(const unsigned int maxval,const unsigned int xsize,const unsigned int ysize, const char* image_name, int argc, char** argv){

  unsigned int procs, rank, mpi_provided_thread_level;
  unsigned int size=xsize*ysize;

  //printf("Generale - xsize: %d\n",xsize);
  MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED,
  &mpi_provided_thread_level);

  if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
    printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
    MPI_Finalize();
    exit( 1 );
  }

  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  unsigned int rest=size%procs;
  //printf("hello from proc %d/%d\n", rank, procs);
  //printf("rest for proc %d: %d\n", rank, rest);

  //printf("xsize: %d\n",xsize);
  unsigned int buffer_size=size/procs+(rank<rest);
  int seed=SEED*rank+procs;
  srand(seed);
  printf("buffer_size for proc %d: %d\n", rank, buffer_size);

  unsigned char* local_buffer = malloc(buffer_size*sizeof(unsigned char) );

  for ( int i = 0; i < buffer_size; i++ ){
    local_buffer[i] = maxval*round((double)rand()/(double)RAND_MAX);
    //printf("Processor %d iteration %d = %d\n",rank, i, local_buffer[i]);
  }
  write_pgm_image(image_name, rank, buffer_size, rest);
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
}*/
