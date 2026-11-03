#include <stdio.h>  
#include<stdlib.h>
#include <mpi.h>

int main(int argc, char **argv)
{
    int me,nproc,me_grid,row,dim;
    int *period, reorder,*ndim;
    MPI_Comm comm_grid; /* definizione di tipo communicator */  
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&menum);  
    MPI_Comm_size(MPI_COMM_WORLD,&nproc);

    row=3;
    dim=2;
    int coordinate[dim];
    ndim = (int*)calloc(dim,sizeof(int));

    period = (int*)calloc(dim,sizeof(int));
    period[0] = period[1] = 0; 
    reorder=0;
    /* Definizione della griglia bidimensionale */
    MPI_Cart_create(MPI_COMM_WORLD,dim,ndim,period,reorder,  &comm_grid);
    MPI_Comm_rank(comm_grid,&menum_grid);
    /* Definizione delle coordinate di ciascun processo  nella griglia bidimensionale */
    MPI_Cart_coords(comm_grid,menum_grid,dim,coordinate);
    
    MPI_Finalize();  
    return 0; 
}