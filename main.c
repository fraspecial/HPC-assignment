#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

#define MAXVAL 255
#define K 200
#define EVO 1
#define N 10
#define STEPS 1


#ifndef INIT
#define INIT
void initialize_image(const unsigned long rows, const unsigned long cols, const unsigned int maxval, const char* filename, int* argc, char** argv[]);
#endif

#ifndef READ
#define READ

void read_pgm_image(unsigned char **image, unsigned long* xisize, unsigned long* cols, unsigned int* maxval, const char *image_name, int * argc, char ** argv[]);
#endif


#ifndef EVOLVE
#define EVOLVE
void evolve_static(unsigned char** ptr, const unsigned int rows, const unsigned int cols, const unsigned int maxval, const char* filename, int* argc, char** argv[]);
#endif


int main( int argc, char **argv )
{
  unsigned int rows=K, cols=K;
  unsigned int maxval=MAXVAL, n=10, steps=STEPS, evo=EVO, ini_flag=0, run_flag=0, opt;
  unsigned char* ptr;
  char* filename=NULL;
  opterr=0;

  while ((opt = getopt (argc, argv, ":irk:e:f:n:s:")) != -1)
  switch (opt){
    case 'i':
      ini_flag = 1;
      break;
    case 'r':
      run_flag = 1;
      break;
    case 'k':
      rows=strtoul(optarg, NULL, 10);
      cols=rows;
      break;
    case 'e':
      evo=strtoul(optarg, NULL, 2);
      break;
    case 'f':
      filename=optarg;
      break;
    case 'n':
      n=strtoul(optarg, NULL, 10);
      break;
    case 's':
      steps=strtoul(optarg, NULL, 10);
      break;

    case '?':
      if (isprint(optopt))
        fprintf (stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
      return 1;
    default:
      abort ();
    }


    if(ini_flag==1){
	    //printf("rows: %d\n",rows);
	    remove(filename);
	    initialize_image(rows, cols, maxval, filename,&argc, &argv);
    }
    if(run_flag==1 && n>=steps){
      evolve_static(&ptr, rows, cols, maxval, filename, &argc, &argv);
      //find_neighbors_all_cells(ptr, grid, rows, cols);
      //evolve(grid,maxval, rows, cols,n,steps);
      //free(grid);
    }
    //write_pgm_image( ptr,  maxval, rows, cols, "check.pgm", ini_flag, run_flag);


    //if(run_flag==1)


    // ---------------------------------------------------
    //

    /*
    rows = 0;
    cols = 0;
    maxval = 0;

    read_pgm_image( &ptr, &maxval, &rows, &cols, "image.pgm");
    printf("The gradient has been read in again\n");

    free( ptr );

    read_pgm_image( &ptr, &maxval, &rows, &cols, "check_me.pgm");
    printf("The imaget has been read\n");

    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
    swap_image( ptr, rows, cols, maxval);

    // do something on the image (for instance, blur it)
    // ...
    //

    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
    swap_image( ptr, rows, cols, maxval);

    write_pgm_image( ptr, maxval, rows, cols, "check_me.back.pgm");
    printf("The imaget has been written back\n");

    free(ptr);
    */

    return 0;
  }
