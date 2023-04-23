#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

void write_pgm_image(const char *image_name, unsigned char* local_buffer,char * header, unsigned short int header_size, unsigned int rank, unsigned int buffer_size, const unsigned int rest)
/*
* image        : a pointer to the memory region that contains the image
* maxval       : either 255 or 65536
* xsize, ysize : x and y dimensions of the image
* image_name   : the name of the file to be written
*
*/
{
  MPI_File fh;
 
  MPI_File_open(MPI_COMM_WORLD, image_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
  if(rank==0)
	  MPI_File_write_at(fh, 0, header, header_size, MPI_CHAR, MPI_STATUS_IGNORE);
  MPI_File_set_view(fh, (header_size+1)+(rank*buffer_size) + rest*(rank + 1 > rest), MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL);

  //  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_write_all(fh, local_buffer, buffer_size, MPI_CHAR, MPI_STATUS_IGNORE);

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


void read_pgm_image(unsigned char **image, int *maxval, int *xsize, int *ysize, const char *image_name)
/*
* image        : a pointer to the pointer that will contain the image
* maxval       : a pointer to the int that will store the maximum intensity in the image
* xsize, ysize : pointers to the x and y sizes
* image_name   : the name of the file to be read
*
*/
{
  FILE* image_file;
  image_file = fopen(image_name, "r");

  *image = NULL;
  *xsize = *ysize = *maxval = 0;

  char    MagicN[2];
  char   *line = NULL;
  size_t  k, n = 0;

  // get the Magic Number
  k = fscanf(image_file, "%2s%*c", MagicN );

  // skip all the comments
  k = getline( &line, &n, image_file);
  while ( (k > 0) && (line[0]=='#') )
  k = getline( &line, &n, image_file);

  if (k > 0)
  {
    k = sscanf(line, "%d%*c%d%*c%d%*c", xsize, ysize, maxval);
    if ( k < 3 )
    fscanf(image_file, "%d%*c", maxval);
  }
  else
  {
    *maxval = -1;         // this is the signal that there was an I/O error
    // while reading the image header
    free( line );
    return;
  }
  free( line );

  int color_depth = 1 + ( *maxval > 255 );
  unsigned int size = *xsize * *ysize * color_depth;

  if ( (*image = (unsigned char*)malloc( size )) == NULL )
  {
    fclose(image_file);
    *maxval = -2;         // this is the signal that memory was insufficient
    *xsize  = 0;
    *ysize  = 0;
    return;
  }

  if ( fread( *image, 1, size, image_file) != size )
  {
    free( image );
    image   = NULL;
    *maxval = -3;         // this is the signal that there was an i/o error
    *xsize  = 0;
    *ysize  = 0;
  }

  fclose(image_file);
  return;
}
