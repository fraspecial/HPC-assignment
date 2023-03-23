#define SEED 15051997

unsigned char * initialise_image( int maxval, int xsize, int ysize, int argc, char** argv){

  srand(SEED);
  unsigned int size, rank, mpi_provided_threaD_level;
  unsigned int rest=xsize%size;
  unsigned int buffer_size=(xsize/size+(rank<rest))*ysize;

  MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED,
  &mpi_provided_thread_level);

  if ( mpi_provided_thread_level < MPI_THREAD_FUNNELED ) {
    printf("a problem arise when asking for MPI_THREAD_FUNNELED
    level\n");
    MPI_Finalize();
    exit( 1 );
  }

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  unsigned char* local_buffer = malloc(buffer_size, sizeof(unsigned char) );
  for ( int i = 0; i < buffer_size; i++ ){
    unsigned char value = maxval*round((double)rand()/(double)RAND_MAX);
    cImage[i] = value;
  }

  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, 'filename',
                MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,
                &fh);
  MPI_File_set_view(fh, (size - rest) * buffer_size , MPI_T, MPI_T,
                    "native", MPI_INFO_NULL);

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_File_write_all(fh, local_buffer, N / size, MPI_T, MPI_STATUS_IGNORE);
  if (rank == 0)
      printf("Wrote %llu B in %f s, %f GiB/s\n",
             N * sizeof(T), duration,
             N * sizeof(T) / (duration * 1024 * 1024 * 1024));

  MPI_File_close(&fh);
  MPI_Finalize();
}

unsigned char * generate_image(struct Cell * grid, const  int size){
  unsigned char* cImage;

  cImage = calloc(size, sizeof(unsigned char) );
  for ( int i = 0; i < size; i++ ){
    unsigned char value = *((grid+i)->value);
    cImage[i] = value;
  }
  return cImage;
}
