#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <string.h>

void read_pgm_image(unsigned char **image, unsigned long* xisize, unsigned long* cols, unsigned int* maxval, const char *image_name, const int procs, const int rank);
void set_parameters(unsigned int size, unsigned int procs, unsigned int rank ,unsigned int* rest, unsigned int* buffer_size, MPI_Offset * disp);

unsigned char compute_neighbor(unsigned int* neighbors, unsigned char* this_batch, unsigned char* nb_before, unsigned char* nb_after, unsigned long cols, MPI_Offset first, MPI_Offset last){
    int sum=0;
    printf("First=%lld; Last=%lld\n", first, last);
    printf("NB_BEFORE:\n");
    for (int i=0; i < cols ; i++){
        printf("%u\n", nb_before[i]);
    }
    printf("NB_AFTER:\n");
    for(int i=0; i < cols ; i++){
        printf("%u\n", nb_after[i]);
    }


    for(int i=0; i < 8; i++){
        unsigned int nb_pos=neighbors[i];
        int dist = nb_pos-first;
        unsigned int nb_value = (nb_pos <= last) * (dist >= 0) * this_batch[dist] +
        (dist < 0) * nb_before[cols+dist] +
        (dist > (last - first)) * nb_after[nb_pos - last - 1];
        printf("Iteration %d, nb_pos %u, dist = %d,  nb_value = %u\n", i, nb_pos, dist, nb_value);
        sum+=nb_value;
    }

    
    return sum;
} 

void evolve_static(unsigned char** ptr, const char* filename, unsigned long* rows, unsigned long* cols, unsigned int* maxval, int* argc, char** argv[]){
  //printf("START:\n");
  //for(int grid_index=0; grid_index<size; grid_index++)
    //printf("\tPosizione %d: indirizzo %p valore %d\n", grid_index, (grid+grid_index)->value, *((grid+grid_index)->value));
  int procs, rank;
  MPI_Comm cart_coord;
  int periodic[]={1,1};
  int coords[2];

  MPI_Init(argc,argv);
  MPI_Comm_size(MPI_COMM_WORLD, &procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  unsigned int size, rest, color_depth, mpi_provided_thread_level;

  MPI_Offset initial_offset=0;


  char MagicN[3];
  unsigned int buffer_size;
  MPI_Offset proc_offset=0;
 if(rank==0){
    FILE* image_file;
    char   *line = NULL;
    size_t  k, n = 0;

    color_depth=1+(*maxval>255);
    image_file = fopen(filename, "r");

    *ptr = NULL;
    *rows = *cols = *maxval = 0;

    // get the Magic Number
    if(fscanf(image_file, "%2s%*c", MagicN ) >0){
	    initial_offset=initial_offset+1+strlen(MagicN);
    }
    printf("Offset dopo magic number: %lld\n",initial_offset);
    // skip all the comments
    k = getline( &line, &n, image_file);
    initial_offset+=k;

    printf("Offset dopo seconda riga: %lld\n",initial_offset);
    while ( (k > 0) && (line[0]=='#') ){
      k = getline( &line, &n, image_file);
      initial_offset+=k;
    }

    if (k > 0)
    {
      k = sscanf(line, "%lu%*c%lu%*c%d%*c", rows, cols, maxval);
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
  MPI_Bcast(&initial_offset, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(rows, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Bcast(cols, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  size= *rows * *cols * (1 + *maxval>255);

    int dims[2] = {0, 0};
    MPI_Dims_create(size, 2, dims);
 
    // Make both dimensions periodic
    int periods[2] = {1, 1};
 
    // Let MPI assign arbitrary ranks if it deems it necessary
    int reorder = 1;
 
    // Create a communicator given the 2D torus topology.
    MPI_Comm new_communicator;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &new_communicator);
 
    // My rank in the new communicator
    int my_rank;
    MPI_Comm_rank(new_communicator, &my_rank);
 
    // Get my coordinates in the new communicator
    int my_coords[2];
    MPI_Cart_coords(new_communicator, my_rank, 2, my_coords);
 
    // Print my location in the 2D torus.
    printf("[MPI process %d] I am located at (%d, %d).\n", my_rank, my_coords[0],my_coords[1]);
 
    MPI_Finalize();
  int dim[]={7, 7};
  //printf("%d %d\n,", dim[0], dim[1]);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periodic, 1, &cart_coord);
  MPI_Comm_rank(cart_coord, &my_rank);
  int return_val = MPI_Cart_create(MPI_COMM_WORLD,2,dim, periodic,1,&cart_coord);
  /*  
  read_pgm_image(ptr, rows, cols, maxval, filename, procs, rank);
  const unsigned long size= *rows* *cols;
  unsigned char* nb_before=malloc(buffer_size);
  unsigned char* nb_after=malloc(buffer_size);
  MPI_Offset disp;
  MPI_Request req;
  MPI_Status status;
  rest=size%procs;
  set_parameters(size, procs, rank , &rest, &buffer_size, &disp);
  
  unsigned int proc_before = (rank+1) < procs ? (rank + 1) : 0;
  unsigned int proc_after = (rank-1) >= 0 ? (rank - 1) : (procs - 1);
  //printf("Proc %d: up-> %u, down->%u\n", rank, proc_up, proc_down);

  //unsigned int last_up=first_up + *rows-1;
  //unsigned int last_down=;
  //unsigned int first_down=first_up+buffer_size -1 - *rows;
  
  //for(int i=0; i < buffer_size; i++){
    int i = 3;
    unsigned int matrix_pos = disp+i;
    
    unsigned int cl_index = (matrix_pos) % *cols > 0 ? matrix_pos - 1: matrix_pos + *cols - 1;
    unsigned int cr_index = (matrix_pos + 1) % *cols > 0 ? matrix_pos + 1: matrix_pos - *cols +1;

    unsigned int tc_index = matrix_pos  > *cols ? matrix_pos - *cols: matrix_pos + *cols * (*rows-1);
    unsigned int tl_index = (tc_index) % *cols > 0 ? tc_index - 1: tc_index + *cols - 1;
    unsigned int tr_index = (tc_index + 1) % *cols > 0 ? tc_index + 1: tc_index - *cols +1;
    
    unsigned int bc_index = matrix_pos/(*cols) < (*rows-1) ? matrix_pos + *cols: matrix_pos % *cols;
    unsigned int bl_index = (bc_index) % *cols > 0 ? bc_index - 1: bc_index + *cols - 1;
    unsigned int br_index = (bc_index + 1) % *cols > 0 ? bc_index + 1: bc_index - *cols +1;

    unsigned int neighbors[]={tl_index, tc_index, tr_index, cl_index, cr_index,bl_index, bc_index, br_index};
    //}
    MPI_Isend(*ptr, *cols, MPI_UNSIGNED_CHAR, proc_before, 1, MPI_COMM_WORLD, &req);
    MPI_Recv(nb_after, *cols, MPI_UNSIGNED_CHAR, proc_after, 1 , MPI_COMM_WORLD, &status);
    MPI_Isend(*ptr + (buffer_size - *cols), *cols, MPI_UNSIGNED_CHAR, proc_after, 0, MPI_COMM_WORLD, &req);
    MPI_Recv(nb_before, *cols, MPI_UNSIGNED_CHAR, proc_before, 0 , MPI_COMM_WORLD, &status);
    if(rank==3)
    compute_neighbor(neighbors, *ptr, nb_before, nb_after, *cols, disp, disp+buffer_size-1);
    */
    MPI_Finalize();
    
  return;
}