#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void write_header(const char* filename, const char * header){
  FILE* f;
  f = fopen(filename, "w");
  fprintf(f, header);
  fclose(f);
}

unsigned int set_parameters(unsigned int size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, unsigned int* disp){

  *rest= size % procs;
  *buffer_size = size/procs+( rank < *rest);
  //printf("hello from proc %d/%d\n", rank, procs);
  *disp=(rank * *buffer_size) + *rest * (rank + 1 > *rest);
}

void write_pgm_image(const char *image_name, const unsigned char* local_buffer, const unsigned int rank, const unsigned int buffer_size, const unsigned int offset)
/*
* image        : a pointer to the memory region that contains the image
* maxval       : either 255 or 65536
* xsize, ysize : x and y dimensions of the image
* image_name   : the name of the file to be written
*
*/
{
  int err;
  MPI_File fh;
  MPI_Status status;
  MPI_Request req;

  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  //CHECK_ERR(MPI_File_set_view);
  MPI_File_write_all(fh, local_buffer, buffer_size, MPI_CHAR,&status);
 //CHECK_ERR(MPI_File_write_all);
  printf("Processor %d wrote %d bytes\n", rank, buffer_size);

  MPI_File_close(&fh);


    // Writing header
    // The header's format is as follows, all in ASCII.
    // "whitespace" is either a blank or a TAB or a CF or a LF
    // - The Magic Number (see below the magic numbers)
    // - the image's width
    // - the height
    // - a white space
    // - the image's height
    // - a whitespace
    // - the maximum color value, which must be between 0 and 65535
    //
    // if he maximum color value is in the range [0-255], then
    // a pixel will be expressed by a single byte; if the maximum is
    // larger than 255, then 2 bytes will be needed for each pixel
    //
    // Writing file

}


void read_pgm_image(unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name, int * argc, char ** argv[])
/*
* image        : a pointer to the pointer that will contain the image
* maxval       : a pointer to the int that will store the maximum intensity in the image
* xsize, ysize : pointers to the x and y sizes
* image_name   : the name of the file to be read
*
*/
{
  unsigned int procs, rank, buffer_size, mpi_provided_thread_level;
  unsigned int size= *xsize * *ysize * color_depth;;
  MPI_Init(NULL,NULL);
  if(rank==0){
    FILE* image_file;
    char    MagicN[2];
    char   *line = NULL;
    size_t  k, n = 0;
    unsigned int color_depth=1+(*maxval>255);

    image_file = fopen(image_name, "r");

    *image = NULL;
    *xsize = *ysize = *maxval = 0;

    // get the Magic Number
    k = fscanf(image_file, "%2s%*c", MagicN );
    offset+=k;
    // skip all the comments
    k = getline( &line, &n, image_file);
    offset+=k;

    while ( (k > 0) && (line[0]=='#') ){
      k = getline( &line, &n, image_file);
      offset+=k;
    }
    if (k > 0)
    {
      k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
      if ( k < 3 ){
        k=getline( &line, &n, image_file);
        offset+=k;
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
    fclose(image_file);
    free( line );
    printf("offset:%s\n", offset);
  }

  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_File fh;
  MPI_Offset offset;
  MPI_Status status;
  set_parameters(size,procs,rank, &rest,&buffer_size, &offset);
  if ( (*image = (unsigned char*)malloc( buffer_size )) == NULL )
  {

    *maxval = -2;         // this is the signal that memory was insufficient
    *xsize  = 0;
    *ysize  = 0;
    return;
  }
  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY,MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, *image, buffer_size, MPI_CHAR, &status);
  MPI_File_close(&fh)
  MPI_File_open(MPI_COMM_WORLD, "Nuova.pgm", MPI_MODE_CREATE|MPI_MODE_WRONLY,MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  MPI_File_write_all(fh, *image, buffer_size, MPI_CHAR, &status);
  MPI_File_close(&fh);
  return;
}
