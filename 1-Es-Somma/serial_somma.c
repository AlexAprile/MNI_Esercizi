/*
    Esercitazione 1:
	1. Sviluppare e implementare in linguaggio C un algoritmo seriale per il
	calcolo della somma di n numeri reali.

    1.1 Somma Seriale

	compilazione:
	gcc -o a serial_somma.c -I "%MSMPI_INC%" -L "%MSMPI_LIB64%" -l msmpi
	esegui:
	mpiexec -np 1 a

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h> // uso per calcolo tempi, per confronto con parallelo

int main (int argc, char **argv){
    double sum;
	double T_inizio,T_fine;
	double *vett;
	int n,i;
	
	printf("Inserire il numero di elementi da sommare: \n");
	fflush(stdout);
	scanf("%d",&n);

	vett=(double*)calloc(n,sizeof(double));
    sum=0.0;

	MPI_Init(&argc,&argv);

	double min =-3.0, max = 3.0; // scelgo intervallo per i reali
	srand((unsigned int) time(0));  // necessaria per avviare il punto di partenza altrimenti è predefinito a 1
    for(i=0;i<n;i++){
		double r = (double)rand() / (double)RAND_MAX; // normalizzo per virgola
		*(vett+i) = min + r * (max-min); // cosi' reali in un intervallo 
    }

	/* --- Sola operazione di somma: da calcolarne il tempo --- */

	// registriamo tempo di inizio
	T_inizio=MPI_Wtime();

	// somma con algoritmo seriale
	for(i=0;i<n;i++){
		sum += vett[i];
	}

	// calcoliamo tempo della sola esecuzione
	T_fine=MPI_Wtime()-T_inizio;

	// Stampa del vettore che contiene i dati da sommare, se sono meno di 100 
	if (n<100)
	{
		for (i=0; i<n; i++)
		{
			printf("Elemento %d del vettore =%.3f \n",i,*(vett+i));
		}
	}

    printf("\nLa somma e' : %.3f\n", sum);
	
	printf("\nTempo di esecuzione *Somma Seriale*: %.9f\n",T_fine);

	MPI_Finalize();
}//fine programma