/*
    1° Strategia Prodotto Matrice-Vettore
    Implementazione del prodotto matrice vettore utilizzando la strategia 1:
    ogni processo riceve un blocco di righe della matrice e l'intero vettore e ne calcola il prodotto.

    funziona solo per matrici con righe multiple dei processi, altrimenti
    tali righe vengono perse.

    Compilazione rapida creata:
    mympicc 1strat-prod-mat-vett.c
    comando rapido per eseguire:
    mpiexec -np 2 1strat-prod-mat-vett
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

// Funzione che esegue il prodotto matrice vettore
void prod_mat_vett(double w[], double *a, int ROWS, int COLS, double v[]){
    int i, j;
    
    for(i=0;i<ROWS;i++){
        w[i]=0;
        for(j=0;j<COLS;j++){ 
            w[i] += a[i*COLS+j]* v[j];
        } 
    }    
}


int main(int argc, char **argv) {

    int nproc;              // Numero di processi totale
    int me;                 // Il mio id
    int m,n;                  // Dimensione della matrice
    int local_m;            // Dimensione dei dati locali
    int i,j;                    // Iteratori vari 
    double T_inizio,T_fine;

    // Variabili di lavoro
    double *A, *v, *localA,*local_w, *w;

    /*Attiva MPI*/
    MPI_Init(&argc, &argv);
    /*Trova il numero totale dei processi*/
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    /*Da ad ogni processo il proprio numero identificativo*/
    MPI_Comm_rank (MPI_COMM_WORLD, &me);

    // Se sono radice
    if(me == 0){
        /* --- INIZIALIZZAZIIONE --- */
        printf("inserire numero di righe m e colonne n = \n"); 
        fflush(stdout);
        scanf("%d",&m); 
        
        n=m;
        //printf("inserire numero di colonne n = \n"); 
        //fflush(stdout);
        //scanf("%d",&n);
        // Numero di righe da processare - righe/nprocessi (solo per righe multiplo dei processi, altrimenti le righe in eccesso non considerate)
        local_m = m/nproc; 
        
        // Alloco spazio di memoria
        A = malloc(m * n * sizeof(double));
        v = malloc(sizeof(double)*n);
        w =  malloc(sizeof(double)*m); 
        
        for (j=0;j<n;j++)
                v[j]=j; 

        for (i=0;i<m;i++){
            for(j=0;j<n;j++){
                if (j==0)
                    A[i*n+j]= 1.0/(i+1)-1;
                else
                    A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j); 
            }
        }
        
        
        if (n<11 && m<11){  
            printf("v = \n");   
            for (j=0;j<n;j++)
                printf("%f ", v[j]);
            printf("\n");
        
        printf("A = \n"); 
            for (i=0;i<m;i++){
                for(j=0;j<n;j++)
                        printf("%f ", A[i*n+j] );
                printf("\n");
            }
        } //endif

    } // fine me==0

    fflush(stdout);

    /* --- DISTRIBUZIONE DATI --- */
    // Spedisco m, n, local_m e v
    MPI_Bcast(&m,1,MPI_INT,0,MPI_COMM_WORLD);  
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);            
    MPI_Bcast(&local_m,1,MPI_INT,0,MPI_COMM_WORLD);            

    // Se sono un figlio alloco v 
    if(me != 0)
        v = malloc(sizeof(double)*n);
        
    MPI_Bcast(&v[0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);            

    // tutti allocano A locale e il vettore dei risultati
    localA = malloc(local_m * n * sizeof(double));
    local_w = malloc(local_m * sizeof(double));

    // Adesso 0 invia a tutti un pezzo della matrice
    // righe da inviare * numero colonne = num numero tot elem da inviare per ogni processore
    int num = local_m*n;    
    MPI_Scatter(
        &A[0], num, MPI_DOUBLE,
        &localA[0], num, MPI_DOUBLE,
        0, MPI_COMM_WORLD);

    // Scriviamo la matrice locale ricevuta
    if (n<11 && m<11){
        printf("localA %d = \n", me); 
        for(i = 0; i < local_m; i++){
            for(j = 0; j < n; j++)
                    printf("%lf\t", localA[i*n+j]);
            printf("\n");
            }
    }
        
    /* --- CALCOLO PROD. PARZIALI --- */
    // sincronizzazione dei processori del contesto MPI_COMM_WORLD
    MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
    T_inizio=MPI_Wtime(); // calcolo del tempo di inizio
    
    // Effettuiamo i calcoli
    prod_mat_vett(local_w,localA,local_m,n,v);
        
    // Calcolo del tempo pdi fine prima della Gather
    MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
    T_fine = MPI_Wtime()-T_inizio; // calcolo del tempo di fine

    // --- COMBINAZIONE DEI RISULTATI
    // 0 raccoglie i risultati parziali
    MPI_Gather(&local_w[0],local_m,MPI_DOUBLE,&w[0],local_m,MPI_DOUBLE,0,MPI_COMM_WORLD);

    // --- STAMPA RISULTATO ---
    // 0 stampa la soluzione
    if(me==0){ 
        if (n<11 && m<11){
            printf("w = \n"); 
            for(i = 0; i < m; i++)
                printf("%f ", w[i]);
            printf("\n");
        }

        printf("\nProcessori impegnati: %d\n", nproc);
        printf("\nTempo calcolo locale: %lf\n", T_fine);
    }

    MPI_Finalize (); /* Disattiva MPI */
    return 0;  
}
