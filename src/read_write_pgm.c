#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>


unsigned int write_header(const char* filename, const unsigned int rows, const unsigned int cols, const unsigned int maxval, const char magic[2]);

void set_parameters(unsigned int size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset* disp);

void write_pgm_image(const char *image_name, const unsigned char* local_buffer, const unsigned int rank, const unsigned int buffer_size, MPI_Offset offset)
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
  MPI_Status status;
  MPI_Request req;

  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

  MPI_File_set_view(fh, offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  //CHECK_ERR(MPI_File_set_view);
  MPI_File_write_all(fh, local_buffer, buffer_size, MPI_CHAR,&status);
 //CHECK_ERR(MPI_File_write_all);
  for (int i = 0; i < buffer_size; i++){

  printf("Processor %u: byte %lld=%u\n", rank, offset+i, local_buffer[i]);
  }
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

/*
void read_pgm_image(unsigned char **image,  unsigned long* rows, unsigned long* cols, unsigned int* maxval, const char* image_name, const int procs, const int rank){

* image        : a pointer to the pointer that will contain the image
* maxval       : a pointer to the int that will store the maximum intensity in the image
* rows, cols : pointers to the x and y sizes
* image_name   : the name of the file to be read
*

  MPI_File fh;
  MPI_Status status;
  rest = size % procs;
  set_parameters(size,procs,rank, &rest,&buffer_size, &proc_offset);
  proc_offset+=initial_offset;

  printf("Proc. %u ; proc_offset: %lld\n", rank, proc_offset);

  if  ((*image = (unsigned char*)malloc(buffer_size*sizeof(unsigned char))) == NULL )
  {

    *maxval = -2;         // this is the signal that memory was insufficient
    *rows  = 0;
    *cols  = 0;
    return;
  }

  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_RDONLY,MPI_INFO_NULL, &fh);
  MPI_File_set_view(fh, proc_offset, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);
  MPI_File_read_all(fh, *image, buffer_size, MPI_CHAR, &status);
  MPI_File_close(&fh);
  for(int i=0; i<buffer_size;i++){
	  printf("Processo %u ha letto byte %lld valore %u\n", rank, proc_offset+i, (*image)[i]);
  }

  if(rank==0){
	  printf("sizeof MagicN=%u, strlen MagicN=%u\n", sizeof(MagicN),strlen(MagicN));
	  printf("MagicN[0]=%c, MagicN[1]=%c, MagicN[2]=%c, MagicN[3]=%c\n", MagicN[0],MagicN[1],MagicN[2], MagicN[3]);
	  unsigned int hed=write_header("Nuovissima.pgm", *rows, *cols, *maxval, MagicN);
	  printf("Verità %u offset %u soggetivo %u\n", hed, initial_offset, proc_offset-initial_offset);
  }
  write_pgm_image("Nuovissima.pgm", *image, rank, buffer_size, proc_offset);

  return;
}
*/