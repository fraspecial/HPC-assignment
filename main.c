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

int main( int argc, char **argv )
{
  int xsize=K, ysize=K, n=10, steps=STEPS, evo=EVO, ini_flag=0, run_flag=0, opt;
  char* filename=NULL;
  unsigned char* ptr=NULL;
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
      xsize=strtoul(optarg, NULL, 10);
      ysize=xsize;
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
      ptr=initialise_image( maxval, xsize, ysize, argc, argv);
      write_pgm_image( ptr,  maxval, xsize, ysize, filename);
    }
    MPI_Finalize();


    else if(run_flag==1 && n>=steps){
      read_pgm_image(&ptr, &maxval, &xsize, &ysize, filename);
      struct Cell* grid = calloc(xsize*ysize, sizeof(struct Cell));
      find_neighbors_all_cells(ptr, grid, xsize, ysize);
      evolve(grid,maxval, xsize, ysize,n,steps);
      free(grid);
    }


    //write_pgm_image( ptr,  maxval, xsize, ysize, "check.pgm", ini_flag, run_flag);


    //if(run_flag==1)

    free(ptr);

    // ---------------------------------------------------
    //

    /*
    xsize = 0;
    ysize = 0;
    maxval = 0;

    read_pgm_image( &ptr, &maxval, &xsize, &ysize, "image.pgm");
    printf("The gradient has been read in again\n");

    free( ptr );

    read_pgm_image( &ptr, &maxval, &xsize, &ysize, "check_me.pgm");
    printf("The imaget has been read\n");

    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
    swap_image( ptr, xsize, ysize, maxval);

    // do something on the image (for instance, blur it)
    // ...
    //

    // swap the endianism
    //
    if ( I_M_LITTLE_ENDIAN )
    swap_image( ptr, xsize, ysize, maxval);

    write_pgm_image( ptr, maxval, xsize, ysize, "check_me.back.pgm");
    printf("The imaget has been written back\n");

    free(ptr);
    */
    return 0;
  }
