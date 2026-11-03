#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

/**
* Funzione che esegue il prodotto matrice vettore
*/
void prod_mat_vett(double w[], double *a, int ROWS, int COLS, double v[])
{
    int i, j;
    
    for(i=0;i<ROWS;i++)
    {
        w[i]=0;
        for(j=0;j<COLS;j++)
        { 
            w[i] += a[i*COLS+j]* v[j];
        } 
    }    
}

int main(int argc, char **argv) {
    
}