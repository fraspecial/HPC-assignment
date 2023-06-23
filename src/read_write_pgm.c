#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>


unsigned int write_header(const char* filename, const unsigned int rows, const unsigned int cols, const unsigned int maxval, const char magic[2]);

void set_parameters(unsigned long rows, unsigned long cols, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset* disp);

void write_pgm_image(const MPI_Comm comm, const char *image_name, const unsigned char* local_buffer, const int rank, const unsigned int buffer_size, const MPI_Offset offset)
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

  MPI_File_open(comm, image_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
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


void read_pgm_image(unsigned char **image, const char* image_name, unsigned int* maxval, unsigned int *rows, unsigned int *cols, unsigned int* rest, unsigned int buffer_size, int procs, int rank){
  
  MPI_File fh;
  MPI_Status status;
  unsigned int size, color_depth;
  MPI_Offset initial_offset=0;
  MPI_Offset proc_offset=0;
  char MagicN[3];
  
  if(rank==0){
    FILE* image_file;
    char   *line = NULL;
    size_t  k, n = 0;

    color_depth=1+(*maxval>255);
    image_file = fopen(image_name, "r");

    *image = NULL;
    *rows = *cols = *maxval = 0;

    // get the Magic Number
    if(fscanf(image_file, "%2s%*c", MagicN ) >0){
	    initial_offset=initial_offset+1+strlen(MagicN);
    }
    //printf("Offset dopo magic number: %lld\n",initial_offset);
    // skip all the comments
    k = getline( &line, &n, image_file);
    initial_offset+=k;

    //printf("Offset dopo seconda riga: %lld\n",initial_offset);
    while ( (k > 0) && (line[0]=='#') ){
      k = getline( &line, &n, image_file);
      initial_offset+=k;
    }

    if (k > 0)
    {
      k = sscanf(line, "%u%*c%u%*c%d%*c", rows, cols, maxval);
      if ( k < 3 ){
        k=getline( &line, &n, image_file);
        initial_offset+=k;
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
    printf("Header length: %lld\n", initial_offset);
    
  }
  MPI_Bcast(&initial_offset, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(rows, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  MPI_Bcast(cols, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  //printf("P%d RICEVE %lu rows e %lu columns\n", rank, *rows, *cols);
  set_parameters(*rows, *cols, procs,rank, rest, &buffer_size, &proc_offset);
  //printf("%lu, %lu\n", *rows, *cols);
  proc_offset+=initial_offset;
  
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
	  printf("Processo %d ha letto byte %lld valore %u\n", rank, proc_offset+i, (*image)[i]);
  }
/*
  if(rank==0){
	  printf("sizeof MagicN=%u, strlen MagicN=%u\n", sizeof(MagicN),strlen(MagicN));
	  printf("MagicN[0]=%c, MagicN[1]=%c, MagicN[2]=%c, MagicN[3]=%c\n", MagicN[0],MagicN[1],MagicN[2], MagicN[3]);
	  unsigned int hed=write_header("Nuovissima.pgm", *rows, *cols, *maxval, MagicN);
	  printf("Verità %u offset %u soggetivo %u\n", hed, initial_offset, proc_offset-initial_offset);
  }
  write_pgm_image("Nuovissima.pgm", *image, rank, buffer_size, proc_offset);
*/
  return;
}
