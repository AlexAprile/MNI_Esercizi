/*
    Esercitazione 3 - 3° Strategia Prodotto Matrice-Vettore
    Implementare in linguaggio C-MPI un programma che effettui la distribuzione 
    di dati (A e x) per il prodotto matrice vettore, secondo la strategia a blocchi
    di righe e colonne per la matrice A. Si supponga di avere a disposizione pxq 
    processi e che la matrice A abbia dimensioni mxn, con m multiplo di p ed n 
    multiplo di q. 

    Implementazione del prodotto matrice vettore utilizzando la strategia 3: 
    ogni processo riceve un blocco di righe e di colonne della matrice 
    e una parte del vettore della stessa dimensione del numero di colonne.
    Viene infine calcola il prodotto, prima localmente e poi globalmente.

    Compilazione rapida creata:
    mympicc 3strat-prod-mat-vett.c
    comando rapido per eseguire:
    mpiexec -np 2 3strat-prod-mat-vett

    provare con 2 4 e 8 processori
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h> 

void prod_mat_trasv_vett(double *w, double *AT, int m, int n, double *v);
// Funzione per calcolare la griglia ottimale
void optimal_grid(int nproc, int rows, int cols, int dim, int *ndim);

int main(int argc, char **argv) {
    int nproc;              // Numero di processi totale
    int me;                 // Il mio id
    int m,n;                  // Dimensione della matrice
    int local_m, local_n;            // Dimensione dei dati locali                  // Iteratori vari 
    double T_inizio,T_fine, T_max; // Variabili per il calcolo del tempo

    double *A, *v, *localA, *localAT, *local_v, *local_w, *w, *finalAT;
    
    // definizione della griglia
    MPI_Comm comm_grid;
    int dim = 2, reorder, *period;
    int *ndim = malloc(dim * sizeof(int));
    int coordinate[2];
    int menum_grid; // id del processo nella griglia

    MPI_Init(&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank (MPI_COMM_WORLD, &me);

    if(me == 0){
        /* --- INIZIALIZZAZIIONE --- */
        printf("inserire numero di righe m e colonne n = \n");
        fflush(stdout);
        scanf("%d",&m); 
        n=m;

        // Alloco spazio di memoria
        A = malloc(m * n * sizeof(double));
        v = malloc(sizeof(double)*n);
        w =  malloc(sizeof(double)*m); 

        // Calcolo la griglia ottimale
        optimal_grid(nproc, m, n, dim, ndim);
        local_m = m/ndim[0]; // numero di righe locale matrice
        local_n = n/ndim[1]; // numero di colonne locale matrice e dimensione locale del vettore
        
        // inizializza v
        for (int j=0;j<n;j++)
                v[j]=j; 
        
        // Inizializza matrice A
        for (int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                if (j==0)
                    A[i*n+j]= 1.0/(i+1)-1;
                else
                    A[i*n+j]= 1.0/(i+1)-pow(1.0/2.0,j); 
            }
        }
        
        // Stampa matrice e vettore se piccoli
        if (n<11 && m<11){  
            printf("v = \n");   
            for (int j=0;j<n;j++)
                printf("%f ", v[j]);
            printf("\n");
        
        printf("A = \n"); 
            for (int i=0;i<m;i++){
                for(int j=0;j<n;j++)
                        printf("%f ", A[i*n+j] );
                printf("\n");
            }
            printf("\n");
        } //endif

    } // fine me==0

    fflush(stdout);

    /* --- DISTRIBUZIONE DATI --- */
    // Spedisco m, n, local_n, local_m e ndim da P0 a tutti i processi
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&local_m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ndim, dim, MPI_INT, 0, MPI_COMM_WORLD);    
    
    /*--- Creazione Contesto ---*/
    // vettore contenente la periodicità delle dimensioni. In questo caso non periodica
    period = malloc(dim * sizeof(int));
    period[0] = period[1] = 0;
    reorder = 0;

    // Definizione della griglia bidimensionale
    MPI_Cart_create(MPI_COMM_WORLD, dim, ndim, period, reorder, &comm_grid);

    // Ogni processore conosce la propria posizione nella griglia
    MPI_Comm_rank(comm_grid, &menum_grid); // id del processore nella griglia
    MPI_Cart_coords(comm_grid, menum_grid, dim, coordinate); // coordinate cartesiane di ciascun processo nella griglia    

    // Sottogriglie per comuinicazioni tra righe e tra colonne
    MPI_Comm row_comm, col_comm;
    int keep_rows[2] = {0,1}; // tieni solo la dimensione colonne (processori nella stessa riga)
    int keep_cols[2] = {1,0}; // tieni solo la dimensione righe (processori nella stessa colonna)
    MPI_Cart_sub(comm_grid, keep_rows, &row_comm);
    MPI_Cart_sub(comm_grid, keep_cols, &col_comm);

    /*--- Invio del vettore v a tutti i processi ---*/
    // Ogni processo alloca local_v
    local_v = malloc(local_n * sizeof(double));

    // Fase 1: processi della prima riga
    if(coordinate[0] == 0) {
        MPI_Scatter(v, local_n, MPI_DOUBLE, 
            local_v, local_n, MPI_DOUBLE,
            0, row_comm);
    }

    // Fase 2: tutti i processi della prima riga inviano il vettore ai processi delle altre righe
    // In ogni riga il processo root ha rank 0
    MPI_Bcast(local_v, local_n, MPI_DOUBLE, 0, col_comm);

    /*--- Invio della Matrice A ad ogni processo ---*/
    // Fase 1: P00 invia blocchi di righe ai processi della prima colonna
    localA = malloc(local_m * n * sizeof(double)); // ogni processo alloca localA
    if(coordinate[1] == 0) {
        MPI_Scatter(A, local_m * n, MPI_DOUBLE,
            localA, local_m * n, MPI_DOUBLE,
            0, col_comm); 
    }

    // Fase 2: ogni processo della prima colonna invia il blocco di righe ai processi delle altre colonne
    // Inizialmente è necessario trasporre localA per inviare blocchi di colonne
    localAT = malloc(n * local_m * sizeof(double));
    for(int i=0; i<local_m; i++){
        for(int j=0; j<n; j++){
            localAT[j * local_m + i] = localA[i * n + j];
        }
    }

    finalAT = malloc(local_m * local_n * sizeof(double));

    // In ogni colonna il processo root ha rank 0
    MPI_Scatter(localAT, local_m * local_n, MPI_DOUBLE,
        finalAT, local_m * local_n, MPI_DOUBLE,0, row_comm); 

    // Stampa di debug
    if(n < 11 && m < 11){
        printf("\nProcesso %d (%d,%d)\n", me, coordinate[0], coordinate[1]);
        
        printf("local_v: ");
        for(int i=0; i<local_n; i++)
            printf("%.2f ", local_v[i]);
        printf("\n");

        printf("finalA\n"); 
        for(int i=0; i<local_m; i++){
            for(int j=0; j<local_n; j++){
                printf("%.4f ", finalAT[j * local_m + i]);
            }
            printf("\n");
        }
    }

    // Ogni processo alloca il vettore per il calcolo locale temporaneo
    double *tmp_w = calloc(local_m, sizeof(double));
    
    // Solo i processi nella prima colonna della griglia raccoglieranno la somma di riga
    // quindi allochiamo local_w solo per loro per evitare memory leak
    if (coordinate[1] == 0) {
        local_w = malloc(local_m * sizeof(double));
    }

    // sincronizzazione dei processori del contesto MPI_COMM_WORLD
	MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
	T_inizio=MPI_Wtime(); // calcolo del tempo di inizio

    /* --- PRODOTTO MATRICE VETTORE --- */
    // ogni processo calcola il pezzo di w relativo alle sue local_m e local_n 
    // inserendolo in tmp_w
    prod_mat_trasv_vett(tmp_w, finalAT, local_m, local_n, local_v);

    /* --- RIDUZIONE (Somma di riga) --- */
    // Sommiamo i vettori tmp_w di tutti i processi all'interno della stessa riga della griglia.
    // Il risultato per ogni riga della griglia viene salvato nel local_w del processo root di quella riga (colonna 0).
    MPI_Reduce(tmp_w, local_w, local_m, MPI_DOUBLE, MPI_SUM, 0, row_comm);

    /* --- Gather dei processi della prima colonna --- */
    // partecipano solo i processi root di ogni riga (cioè i processi con coordinate[1] == 0), 
    // che ora possiedono le somme parziali definitive per quel blocco di righe.
    if (coordinate[1] == 0) {
        MPI_Gather(local_w, local_m, MPI_DOUBLE,
                w, local_m, MPI_DOUBLE,
                0, col_comm); // il processo P00 (root di col_comm) raccoglie il risultato finale
    }

    MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
    T_fine = MPI_Wtime() - T_inizio;

    // tempo massimo
    MPI_Reduce(&T_fine, &T_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    // --- STAMPA RISULTATO ---
    // 0 stampa la soluzione
    if(me==0){ 
        printf("\nProcessori impegnati: %d\n", nproc);
        printf("\nTempo calcolo locale: %lf\n", T_fine);
        printf("\nMPI_Reduce max time: %f\n",T_max);

        if (n<11 && m<11){
            printf("\nRisultato finale w: ");
            for (int i = 0; i < m; i++) printf("%.4f ", w[i]);
            printf("\n");
        } 
        printf("\n");
    }

    // Cleanup
    free(tmp_w);
    if (coordinate[1] == 0) free(local_w);

    MPI_Finalize (); /* Disattiva MPI */
    return 0;  
}

void prod_mat_trasv_vett(double *w, double *AT, int m, int n, double *v) {
    for (int i = 0; i < n; i++) {
        double vi = v[i]; // elemento corrispondente del blocco di vettore
        for (int j = 0; j < m; j++) {
            w[j] += AT[i * m + j] * vi;
        }
    }

    if(n < 11 && m < 11){ 
        printf("local_w: ");
        for (int i = 0; i < m; i++) 
            printf("%.4f ", w[i]);
        printf("\n");
    }
}

void optimal_grid(int nproc, int rows, int cols, int dim, int *ndim){
    int fattori[nproc];
    int nfattori = 0;
    double ratio = (double)rows / cols; // Rapporto desiderato tra righe e colonne
    double best_diff = rows; // Inizializza con un valore grande

    // Trova i fattori di nproc
    for(int i = 1; i <= (nproc / 2); i++){
        if(nproc % i == 0){
            fattori[nfattori] = i;
            nfattori++;
        }
    }
     // Aggiungi nproc stesso come fattore
    fattori[nfattori]=nproc;
    nfattori++;

    // Cerca la coppia di fattori che meglio approssima il rapporto desiderato
    for(int i=0; i < nfattori; i++){
        int f1 = fattori[i];
        int f2 = nproc / f1;

        // Il rapporto di tra le dimensioni di ogni sottogriglia deve essere vicino a quello della griglia
        double current_ratio = (double)(rows/f1) / (cols/f2);
        double diff = abs(current_ratio - ratio);

        if(diff < best_diff){
            best_diff = diff;
            ndim[0] = f1;
            ndim[1] = f2;
        }
    }  

    if(rows < 11 && cols < 11)
        printf("Griglia ottimale: %d x %d\n", ndim[0], ndim[1]);
}
