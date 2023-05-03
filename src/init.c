#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#define SEED 100

#ifndef WRITE
#define WRITE
void write_pgm_image(const char *image, unsigned char* buf, unsigned int rank, unsigned int buf_size, const unsigned int rest);
#endif

void set_parameters(unsigned int size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp){

 
  *buffer_size = size/procs+( rank < *rest);
  //printf("hello from proc %d/%d\n", rank, procs);
  *disp= (rank * *buffer_size) + *rest * (rank + 1 > *rest);
}

unsigned int write_header(const char* filename, const unsigned int xsize, const unsigned int ysize, const unsigned int maxval, const char* magic){
	FILE* f;
    unsigned int length=0;
    
    f=fopen(filename, "w");
    printf("Size of magic: %u. strlen(magic)=%u\n",sizeof(magic),strlen(magic));
    printf("%s è %u caratteri + 1 a capo, %u è %u caratteri + 1 spazio, %u è %u caratteri+ 1 a capo, %u è %u caratteri +1 a capo.\n",magic, strlen(magic),xsize, (int)(ceil(log10(xsize))), ysize, (int)(ceil(log10(ysize))), maxval, (int)(ceil(log10(maxval))));
    length=fprintf(f, "%s\n%u %u\n%u\n",magic, xsize, ysize, maxval);
    printf("Ho scritto su file %u caratteri.\n", length);
    fclose(f);
    return length;
}

void initialise_image(const unsigned int maxval,const unsigned int xsize,const unsigned int ysize, const char* image_name, int* argc, char** argv[]){

  unsigned int procs, rank, mpi_provided_thread_level;
  
  //printf("Header \"%s\", lunghezza %u", header, header_length);

  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Offset disp;
  printf("hello from proc %d/%d. Image name %s\n", rank, procs, image_name);

  int seed=SEED*rank+procs;
  //printf("buffer_size for proc %u: %u\n", rank, buffer_size);
  unsigned int buffer_size, header_length, rest, size;

  if(rank==0){
    header_length=write_header(image_name, xsize, ysize, maxval, "P5");
    printf("Header ha lunghezza %u\n", header_length);
  }
  
  MPI_Bcast(&header_length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
  printf("Ricevuto header length processo %u: %u\n", rank, header_length);
  size=xsize*ysize*1+(maxval>255);
  rest=size%procs;
  set_parameters(size, procs, rank , &rest, &buffer_size,  &disp);
  disp=disp+header_length;
  printf("Proc %u. Header length: %u; rest: %u; buffer_size:%u; proc_offset: %u\n", rank, header_length, rest, buffer_size, disp);
  unsigned char* local_buffer =(unsigned char*)malloc(buffer_size*sizeof(unsigned char) );

  srand(seed);
  for ( int i = 0; i < buffer_size; i++ ){
    local_buffer[i] = maxval*round((double)rand()/(double)RAND_MAX);
  }

  write_pgm_image(image_name, local_buffer, rank, buffer_size, disp);
  free(local_buffer);
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
