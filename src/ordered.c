/*struct Cell {
  unsigned char * value;
  unsigned char * neighbors[8];
};

unsigned char** find_neighbors_single_cell(unsigned char* img, const unsigned int pos, const unsigned int xsize, const unsigned int ysize){
  const unsigned int xx= pos/ysize;
  const unsigned int yy= pos%ysize;
  unsigned char ** nb=calloc(8, sizeof(unsigned char *));

  unsigned int x_t=xx-1;
  unsigned int x_b=xx+1;
  unsigned int y_l=yy-1;
  unsigned int y_r=yy+1;

  if(xx==0){
    x_t=xsize-1;
  }
  else if(xx==(xsize-1)){
    x_b=0;
  }

  if(yy==0){
    y_l=ysize-1;
  }
  else if(yy==(ysize-1)){
    y_r=0;
  }

  unsigned int tl=x_t*ysize+y_l;
  unsigned int tc=x_t*ysize+yy;
  unsigned int tr=x_t*ysize+y_r;
  unsigned int cl=xx*ysize+y_l;
  unsigned int cr=xx*ysize+y_r;
  unsigned int bl=x_b*ysize+y_l;
  unsigned int bc=x_b*ysize+yy;
  unsigned int br=x_b*ysize+y_r;

  unsigned int neighbors[8]={tl,tc,tr,cl,cr,bl,bc,br};
  for(int i=0; i<8; i++){
    nb[i]=&(img[neighbors[i]]);
  }
  return nb;
}


void evolve_cell(struct Cell * cell){
  int sum=0;
  unsigned short int alive = *(cell->value);
  for (int i=0;i <8; i++){
    sum+=*(cell->neighbors[i]);
  }
  //printf("\tCella %s: indirizzo %p valore %d vicini %d\n", alive ? "viva":"morta", (cell->value),*(cell->value), sum/255);
  if(alive){
    if(sum < 255*2 || sum > 255*3 ){
      *(cell->value)=0;
      alive=0;
    }
  }
  else if(sum == 255*3){
    *(cell->value)=255;
    alive=255;
  }
  //printf("\tNuova cella %s: indirizzo %p valore %d vicini %d\n\n", alive ? "viva":"morta", (cell->value),*(cell->value), sum/255);
}

void evolve(struct Cell * grid, const int maxval, const int xsize, const int ysize, const int n, const int step){
  //printf("START:\n");
  //for(int grid_index=0; grid_index<size; grid_index++)
    //printf("\tPosizione %d: indirizzo %p valore %d\n", grid_index, (grid+grid_index)->value, *((grid+grid_index)->value));
  const int size=xsize*ysize;
  for(int i=0; i<n; i++){
    //printf("EVOLUZIONE:\n");
    for(int grid_index=0; grid_index<size; grid_index++){
      evolve_cell(grid + grid_index);
      //printf("\tPosizione %d: indirizzo %p valore %d\n", grid_index, (grid+grid_index)->value, *((grid+grid_index)->value));
    }
    if((i+1) % step ==0){
      unsigned char* ptr=generate_image(grid, size);
      char image_name[17];
      sprintf(image_name, "snapshot%04d.pgm", i+1);
      write_pgm_image(ptr, maxval, xsize, ysize, image_name);
    }
  }
}

void find_neighbors_all_cells(unsigned char *image, struct Cell* grid, const int xsize, const int ysize){

  int size=xsize*ysize;

  for(int grid_index=0; grid_index < size; grid_index++){
    (grid+grid_index)->value=(unsigned char*)&image[grid_index];
    unsigned char** nb=find_neighbors_single_cell(image, grid_index, xsize, ysize);
    for(int nb_index=0; nb_index<8; nb_index++)
      (grid+grid_index)->neighbors[nb_index]=nb[nb_index];
    free(nb);
  }
}*/
