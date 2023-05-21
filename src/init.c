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

void set_parameters(unsigned long size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp){

 
  *buffer_size = size/procs+( rank < *rest);
  //printf("hello from proc %d/%d\n", rank, procs);
  *disp= (rank * *buffer_size) + *rest * (rank + 1 > *rest);
}

unsigned int write_header(const char* filename, const unsigned long rows, const unsigned long cols, const unsigned int maxval, const char* magic){
	FILE* f;
    unsigned long length=0;
    
    f=fopen(filename, "w");
    printf("Size of magic: %lu. strlen(magic)=%lu\n",sizeof(magic),strlen(magic));
    printf("%s è %lu caratteri + 1 a capo, %lu è %lu caratteri + 1 spazio, %lu è %lu caratteri+ 1 a capo,\
    %u è %d caratteri +1 a capo.\n",magic, strlen(magic),rows, (unsigned long)(ceil(log10(rows))),\
    cols, (unsigned long)(ceil(log10(cols))), maxval, (int)(ceil(log10(maxval))));
    length=fprintf(f, "%s\n%lu %lu\n%u\n",magic, rows, cols, maxval);
    printf("Ho scritto su file %lu caratteri.\n", length);
    fclose(f);
    return length;
}

void initialise_image(const unsigned int maxval,const unsigned long rows,const unsigned long cols, const char* image_name, int* argc, char** argv[]){

  int procs, rank, mpi_provided_thread_level;
  MPI_Comm cart_coord;
  int dim[]={rows, cols};
  int periodic[]={1,1};
  //printf("Header \"%s\", lunghezza %lu", header, header_length);

  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periodic, 1, &cart_coord);
  MPI_Comm_rank(cart_coord, &rank);
  
  MPI_Offset disp;
  printf("hello from proc %d/%d. Image name %s\n", rank, procs, image_name);

  int seed=SEED*rank+procs;
  //printf("buffer_size for proc %lu: %lu\n", rank, buffer_size);
  unsigned int buffer_size, header_length, rest;
  unsigned long size;

  if(rank==0){
    header_length=write_header(image_name, rows, cols, maxval, "P5");
    printf("Header ha lunghezza %u\n", header_length);
  }
  
  MPI_Bcast(&header_length,1,MPI_UNSIGNED,0,MPI_COMM_WORLD);
  printf("Ricevuto header length processo %d: %u\n", rank, header_length);
  size=rows*cols*1+(maxval>255);
  rest=size%procs;
  set_parameters(size, procs, rank , &rest, &buffer_size,  &disp);
  disp=disp+header_length;
  printf("Proc %d. Header length: %u; rest: %u; buffer_size:%u; proc_offset: %lld\n", rank, header_length, rest, buffer_size, disp);
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
