#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#define SEED 100

#ifndef WRITE
#define WRITE
void write_header(const char* filename, const char* header);
void write_pgm_image(const char *image, unsigned char* buf, unsigned short int header_size,unsigned int rank, unsigned int buf_size, const unsigned int rest);
#endif


void initialise_image(const unsigned int maxval,const unsigned int xsize,const unsigned int ysize, const char* image_name, int* argc, char** argv[]){

  unsigned int procs, rank, mpi_provided_thread_level;
  unsigned int size=xsize*ysize;
  unsigned int rest;
  //printf("Header \"%s\", lunghezza %u", header, header_length);

  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //printf("hello from proc %d/%d\n", rank, procs);
  int seed=SEED*rank+procs;
  //printf("buffer_size for proc %u: %u\n", rank, buffer_size);
  unsigned int buffer_size;

  MPI_Offset disp;

  set_parameters(&size, &procs, &rank , &rest, &buffer_size, &disp);
  unsigned char* local_buffer =(unsigned char*)malloc(buffer_size*sizeof(unsigned char) );

  srand(seed);
  for ( int i = 0; i < buffer_size; i++ ){
    local_buffer[i] = maxval*round((double)rand()/(double)RAND_MAX);
    //printf("Processor %d iteration %d = %d\n",rank, i, local_buffer[i]);
  }

  if(rank==0){
    char* header=(char*)malloc(30);
    sprintf(header, "P5\n%u %u\n%u\n", xsize, ysize, maxval);
    write_header(image_name, header);
  }

  write_pgm_image(image_name, local_buffer, rank, buffer_size, disp+strlen(header)+1);
  free(local_buffer);
  free(header);
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
