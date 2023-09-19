#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>



unsigned int write_header(const char* filename, const unsigned int rows, const unsigned int cols, const unsigned int maxval, const char magic[2]);
void set_parameters(unsigned int rows, unsigned int cols, unsigned int procs, unsigned int rank, unsigned int* rest, unsigned int* buffer_size, MPI_Offset* disp);

void write_pgm_image(const MPI_Comm* comm, const char *image_name, unsigned char* local_buffer, const unsigned int procs, const unsigned int rank, const unsigned int rows, const unsigned int cols, const unsigned int maxval)
/*
* image        : a pointer to the memory region that contains the image
* maxval       : either 255 or 65536
* rows, cols : x and y dimensions of the image
* image_name   : the name of the file to be written
*
*/
{
  int err;
  MPI_File fh;
  MPI_Offset offset;
  MPI_Status status;
  MPI_Request req;
  unsigned int buffer_size, header_length, rest;

  if(rank==0){
      //printf("RAND MAX IS %d\n", RAND_MAX);
      header_length=write_header(image_name, rows, cols, maxval, "P5");
      //printf("Header ha lunghezza %u\n", header_length);
    }
    
  MPI_Bcast(&header_length,1,MPI_UNSIGNED,0,*comm);
    
    //printf("Ricevuto header length processo %d: %u\n", my_rank, header_length);
  set_parameters(rows, cols, procs, rank , &rest, &buffer_size, &offset);
  offset=offset+header_length;
    //printf("Proc %d. Header length: %u; rest: %u; buffer_size:%u; proc_offset: %lld\n", my_rank, header_length, rest, buffer_size, disp);
    
    
    
  MPI_File_open(*comm, image_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  //CHECK_ERR(MPI_File_set_view);
  MPI_File_write(fh, local_buffer, buffer_size, MPI_CHAR,&status);
 //CHECK_ERR(MPI_File_write_all);
  MPI_File_close(&fh);
}

void read_header(unsigned char *image, const char* image_name, unsigned int *rows, unsigned int *cols, unsigned int* maxval, MPI_Offset* initial_offset){
    char MagicN[3];
    FILE* image_file;
    *initial_offset=0;
    char   *line = NULL;
    size_t  k, n = 0;
    unsigned int color_depth=1+(*maxval>255);
    *rows = *cols = *maxval = 0;
    image = NULL;

    image_file = fopen(image_name, "r");
  
    if(fscanf(image_file, "%2s%*c", MagicN ) >0){
	    *initial_offset=(*initial_offset)+1+strlen(MagicN);
    }

    k = getline( &line, &n, image_file);
    *initial_offset+=k;

    while ( (k > 0) && (line[0]=='#') ){
      k = getline( &line, &n, image_file);
      *initial_offset+=k;
    }

    if (k > 0)
    {
      k = sscanf(line, "%u%*c%u%*c%d%*c", cols, rows, maxval);
      if ( k < 3 ){
        k=getline( &line, &n, image_file);
        *initial_offset+=k;
       }
      sscanf(line, "%d%*c", maxval);
    }

    else
    {
      *maxval = -1;         // this is the signal that there was an I/O error
      // while reading the image header
      free( line );
      return;
    }
    free( line );
    fclose(image_file);
    
    
    //printf("Header length: %lld\n", initial_offset);
    return;
}

unsigned char * read_pgm_image(MPI_Comm* my_world, unsigned char *image, const char* image_name, MPI_Offset initial_offset, int* maxval, unsigned int *rows, unsigned int *cols, unsigned int* rest, unsigned int* buffer_size, int procs, int rank){
  
  
  MPI_File fh;
  MPI_Status status;
  unsigned int size;
  MPI_Offset proc_offset=0;
  
  //printf("P%d RICEVE %u rows e %u columns\n", rank, *rows, *cols);
  size=(*rows) * (*cols);
  set_parameters(*rows, *cols, procs,rank, rest, buffer_size, &proc_offset);
  proc_offset+= initial_offset;
  
  if  ((image = (unsigned char*)malloc((2*(*cols) + *buffer_size)*sizeof(unsigned char))) == NULL ){
    *maxval = -2;         // this is the signal that memory was insufficient
    *rows  = 0;
    *cols  = 0;
    return NULL;
  }
  
  MPI_File_open(*my_world, image_name, MPI_MODE_RDONLY,MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, proc_offset, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_CHAR, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, &(image[*cols]), *buffer_size, MPI_UNSIGNED_CHAR, &status);
  MPI_File_close(&fh);
  return image;
}

void write_csv(int size, int nthreads, unsigned int dim, int n, double global, const char* datacsv){

    FILE *csv_file = fopen(datacsv, "a");
    if (csv_file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

  fprintf(csv_file, "%d,%d,%d,%d,%lf\n", size, nthreads,dim, n,global);

  fclose(csv_file);

  return;

}
