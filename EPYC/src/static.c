#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>
#include <time.h>

void read_header(unsigned char *image, const char* image_name, unsigned int *rows, unsigned int *cols, unsigned int* maxval,  MPI_Offset * initial_offset);
unsigned char* read_pgm_image(MPI_Comm* my_world, unsigned char *image, const char *image_name, MPI_Offset initial_offset,  unsigned int* maxval, unsigned int *rows, unsigned int *cols, unsigned int* rest, unsigned int* buffer_size, int procs, int rank);
void set_parameters(unsigned int rows, unsigned int cols, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp);
void write_pgm_image(const MPI_Comm* comm, const char *image_name, unsigned char* local_buffer, const unsigned int procs, const unsigned int rank, const unsigned int rows, const unsigned int cols, const unsigned int maxval);
MPI_Comm initialize_procs(const int rows, int* my_procs, int* my_rank, int* rank_above, int* rank_below, int coords[2]);
void write_csv(int size, int nthreads, unsigned int N, int n, double global, const char* datacsv);

unsigned char compute_neighbor(unsigned int matrix_pos, unsigned char* this_batch, unsigned int cols, unsigned int maxval){

	unsigned int left = matrix_pos % cols > 0 ? -1 : cols -1;
	unsigned int right = matrix_pos % cols < cols-1 ? 1: -cols+1;
    
    unsigned int sum = this_batch[matrix_pos-cols+left]+this_batch[matrix_pos-cols] + this_batch[matrix_pos-cols+right]+this_batch[matrix_pos+left]+this_batch[matrix_pos+right]+this_batch[matrix_pos+cols+left]+this_batch[matrix_pos+cols]+this_batch[matrix_pos+cols+right];
    return sum >= maxval*2 && sum <= maxval * 3? (unsigned char) maxval : (unsigned char)0;
} 

void evolve_static(unsigned int evolutions, unsigned int steps, unsigned char* ptr, unsigned int rows, unsigned int cols, unsigned int* maxval, const char* filename, const char* output_file, int* argc, char** argv[]){
  
  MPI_Offset initial_offset=0;
  MPI_Request req;
  MPI_Status status;
  unsigned long i, j;
  double global;
  unsigned int rest, buffer_size, color_depth;
  char evo_title[30];
  char MagicN[3];
  int coords[2];
  int procs, cart_rank, proc_above, proc_below, nthreads, mpi_provided_thread_level;
  read_header(ptr, filename, &rows, &cols, maxval, &initial_offset);
  
  MPI_Init_thread(argc, argv, MPI_THREAD_FUNNELED, &mpi_provided_thread_level);
  if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
	  printf("a problem arise when asking for MPI_THREAD_FUNNELED level\n");
	  MPI_Finalize();
	  exit( 1 );
  }
  double t_start=MPI_Wtime();
  
  MPI_Comm cart_comm=initialize_procs(rows, &procs, &cart_rank, &proc_above, &proc_below, coords);
  

  if (cart_rank!=MPI_UNDEFINED){
    if(procs > 1){
      ptr = read_pgm_image(&cart_comm, ptr, filename, initial_offset, maxval, &rows, &cols, &rest, &buffer_size, procs, cart_rank);
      unsigned char * evo = malloc(2*cols+buffer_size);

	  for(int i=0; i < evolutions; i++){
      MPI_Isend(&(ptr[cols]), cols, MPI_UNSIGNED_CHAR, proc_above, 1, cart_comm, &req);
      MPI_Isend(&(ptr[buffer_size]), cols, MPI_UNSIGNED_CHAR, proc_below, 0, cart_comm, &req);
      MPI_Recv(&(ptr[cols + buffer_size]), cols, MPI_UNSIGNED_CHAR, proc_below, 1 , cart_comm, &status);
      MPI_Recv(ptr, cols, MPI_UNSIGNED_CHAR, proc_above, 0 , cart_comm, &status);
          
      #pragma omp parallel for schedule(static)
	    for(int j=cols; j < cols+buffer_size; j++)   
		    evo[j]=compute_neighbor(j, ptr, cols, *maxval);
    
	      if(steps>0)
          if(i%steps==0){
            sprintf(evo_title, "./evos/evo_%d_of_%u.pgm", i+1, evolutions);
            write_pgm_image(&cart_comm,evo_title, &(evo[cols]), procs, cart_rank, rows, cols, *maxval); 
          }
        unsigned char * temp = evo;
        evo = ptr;
        ptr = temp;
      }	
    free(evo);
    double duration=MPI_Wtime()-t_start;
    MPI_Reduce(&duration, &global, 1, MPI_DOUBLE, MPI_MAX, 0, cart_comm);
 	}
  
  else{
    ptr = read_pgm_image(&cart_comm, ptr, filename, initial_offset, maxval, &rows, &cols, &rest, &buffer_size, procs, cart_rank);
    unsigned char * evo = malloc(2*cols+buffer_size);

    for(int i=0; i < evolutions; i++){
    #pragma omp parallel for schedule(static)
      for (int k=0; k<cols; k++){
        ptr[k]=ptr[buffer_size+k];
        ptr[cols+buffer_size+k]=ptr[cols+k];
      }

      #pragma omp parallel for schedule(static)
      for(int j=cols; j < cols+buffer_size; j++)  
		    evo[j]=compute_neighbor(j, ptr, cols, *maxval);
        
        if(steps > 0)
          if (i%steps==0){
            sprintf(evo_title, "./evos/evo_%d_of_%u.pgm", i+1, evolutions);
            write_pgm_image(&cart_comm,evo_title, &(evo[cols]), procs, cart_rank, rows, cols, *maxval);
          }         
        unsigned char * temp = evo;
        evo = ptr;
        ptr = temp;
    }
  
    free(evo);
    global=MPI_Wtime()-t_start;
  }
    if(cart_rank==0)
	    write_csv(procs, omp_get_max_threads(), rows, evolutions, global, output_file);
  }

  MPI_Finalize();
  return;
}
 
