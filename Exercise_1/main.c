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

void initialize_image(const unsigned long rows, const unsigned long cols, const unsigned int* maxval, const char* filename, int* argc, char** argv[]);
void evolve_static(unsigned int evos, unsigned int step, unsigned char* ptr, const unsigned int rows, const unsigned int cols, const unsigned int* maxval, const char* filename, const char* output, int* argc, char** argv[]);
void evolve_ordered(unsigned int evolutions, unsigned int steps, unsigned char* ptr, unsigned int rows, unsigned int cols, unsigned int* maxval, const char* filename, const char* output_file, int* argc, char** argv[]);


int main( int argc, char **argv )
{
  unsigned int rows=K, cols=K;
  unsigned int maxval=MAXVAL, n=10, steps=STEPS, evo=EVO, ini_flag=0, run_flag=0, opt;
  unsigned char* ptr;
  char* filename=NULL;
  char* output_file=NULL;
  opterr=0;

  while ((opt = getopt (argc, argv, "irk:x:y:e:f:o:f:n:s:")) != -1)
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
    case 'x':
      rows=strtoul(optarg, NULL, 10);
      break;
    case 'y':
      cols=strtoul(optarg, NULL, 10);
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
    case 'o':
      output_file=optarg;
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
	    initialize_image(rows, cols, &maxval, filename,&argc, &argv);
    }
    if(run_flag==1 && n>=steps){
      if(evo==0)
        evolve_ordered(n, steps, ptr, rows, cols, &maxval, filename, output_file, &argc, &argv);
      else
        evolve_static(n, steps, ptr, rows, cols, &maxval, filename, output_file, &argc, &argv);
    }

    return 0;
}
