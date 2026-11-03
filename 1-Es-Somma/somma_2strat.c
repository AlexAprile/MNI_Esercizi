/*
	1.2 Somma Parallela con Seconda Strategia

    Esercitazione 1:
	2. Sviluppare e implementare in linguaggio C--MPI un algoritmo parallelo per il
	calcolo della somma di n numeri reali, che utilizzi la II strategia di
	parallelizzazione.

	compilazione:
	gcc -o somma_2strat serial_2strat.c -I "%MSMPI_INC%" -L "%MSMPI_LIB64%" -l msmpi
	oppure:
	mympicc somma_2strat.c
	
	esegui:
	mpiexec -np 8 somma_2strat

*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


int main (int argc, char **argv){
    int me,nproc,tag,sendTo,recvBy;
	int n,i,p,passi,r;
	int nloc,nlocgen,pos,resto, *potenze;
	double T_inizio,T_fine,T_max;
	double *vett, *vett_parz;
	double somma_parz, somma, tmp;
	MPI_Status status;
	
    somma=0.0;
	somma_parz=0.0;

	MPI_Init(&argc,&argv);
	// Identificativo processi: assegnazione IdProcessore a me (id processo)
	MPI_Comm_rank(MPI_COMM_WORLD, &me);
	// assegna numero processori a nproc (num processi)
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	// Setup sul processore principale 0
	if(me==0){
		printf("Inserire il numero di elementi da sommare: \n");
		fflush(stdout);
		scanf("%d",&n);

		vett=(double*)calloc(n,sizeof(double));
	}

	/* invio del valore di n(addendi) a tutti i processori appartenenti a MPI_COMM_WORLD*/
	MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

	
	/* # --- FASE DISTRIBUZIONE --- # */
	
	/* numero di addendi da assegnare a ciascun processore*/
	nlocgen=n/nproc; // divisione intera
	
    resto=n%nproc; // resto della divisione

	// nloc = numero di addendi per ogni singola locazione
	/* Se resto è non nullo, i primi resto processi ricevono un addento in più */
	if(me<resto)
	{
		nloc=nlocgen+1;
	}
	else
	{
		nloc=nlocgen;
	}

	/* allocazione di memoria del vettore per le somme parziali */
	vett_parz=(double*)calloc(nloc,sizeof(double));

	// P0 generera e assegna elementi nel vettore da sommare
	if(me==0){
		double min =-3.0, max = 3.0; // scelgo intervallo per i reali
		// necessaria per avviare il punto di partenza altrimenti è predefinito a 1
		srand((unsigned int) time(0));

		for(int i=0;i<n;i++){
			double r = (double)rand() / (double)RAND_MAX; // normalizzo per virgola
			*(vett+i) = min + r * (max-min); // cosi' reali in un intervallo 
    	}

		// assegnazione solo per P0, con nloc numero di addendi
		for(i=0;i<nloc;i++){
			*(vett_parz+i)=*(vett+i);
		}
		
		// posizione in cui ci troviamo, del numero di addendi che sono stati assegnati
		pos=nloc; // ora solo pos dopo riempimento P0

		/* P0 assegna i restanti addendi agli altri processori disponibili */
		for(i=1; i<nproc; i++){
			tag=i; // tag del messaggio uguale all'id del processo che riceve
			/* i processori che riceverenno dli addendi  in sovrannumero */
			if(i<resto){
				// MPI_Send(msg, count, datatype, dest, tag, comm)
				MPI_Send(vett+pos,nloc,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
				// agggiornamneto index pos
				pos+=nloc;
			}else{	// addendi regolari non in sovrannumero
				MPI_Send(vett+pos,nlocgen,MPI_DOUBLE,i,tag,MPI_COMM_WORLD);
				pos+=nlocgen;
			}
		}
	}
	/* Processori che non sono P0: ricevono i dati da esso */
	else{
		// tag salviamo come numero del processore che riceve
		tag=me;
		// MPI_Recv(msg, count, datatype, source, tag, comm,&status)
		MPI_Recv(vett_parz,nloc,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&status);
		// source 0 = P0
	}//end if

	/* Tutti eseguono da qui */

	/* sincronizzazione dei processori del contesto MPI_COMM_WORLD*/
	MPI_Barrier(MPI_COMM_WORLD); //barrier = ferma i processi (var. interna) tutti del common world
	

	/* --- Sola operazione di somma: da calcolarne il tempo --- */

	// registriamo tempo di inizio
	T_inizio=MPI_Wtime();

	// somma parziale con algoritmo in parallelo
	for(i=0;i<nloc;i++){
		somma_parz += vett_parz[i];
	}

	// Algo 2 strategia
	p=nproc;
	passi=0;
	// calcoliamo il numero di passi che ci serviranno in base al nproc
	while(p!=1){
		p=p>>1; //shift a destra del bit divisione per 2 facciamo prima
		passi++; // numero di passi per l'algoritmo
	}

	/* creazione del vettore potenze, che contiene le potenze di 2*/
	// contiene le distanze progressive (tra i processori)
	potenze=(int*)calloc(passi+1,sizeof(int));
	for(i=0;i<=passi;i++){
		potenze[i]=p<<i; // shift sinistra moltiplicazione per 2
	}

	// calcolo delle somme parziali e combinazione risultati
	for(i=0;i<passi;i++){
		// calcolo id processore
		// resto tra n proc(me) e la grandezza gruppo del proc in quel passo
		r=me%potenze[i+1];

		// INVIA
		// resto = 2^i, processore me invia (si trova a destra dispari)
		if(r==potenze[i]){
			// calcolo dell'id del processore a cui spedire la somma locale
			sendTo=me-potenze[i];
			tag=sendTo;
			// MPI_Send(msg, count, datatype, dest, tag, comm)
			MPI_Send(&somma_parz,1,MPI_DOUBLE,sendTo,tag,MPI_COMM_WORLD);
		}else if(r==0){ // RICEVE
			// sono i proc sempre primi nel loro gruppo
			recvBy=me+potenze[i];
			tag=me;
			MPI_Recv(&tmp,1,MPI_DOUBLE,recvBy,tag,MPI_COMM_WORLD,&status);
			somma_parz+=tmp;
		}
		// else INATTIVO: gli altri hanno già eseguito e non ci sono addendi per loro
	}// end for somme e combinazione ris parziali

	// --------------------------------------------

	MPI_Barrier(MPI_COMM_WORLD); // sincronizzazione
	// calcoliamo tempo della sola esecuzione
	T_fine=MPI_Wtime()-T_inizio;

	/* calcolo del tempo totale di esecuzione*/
	MPI_Reduce(&T_fine,&T_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);

	somma=somma_parz;

	/* Un solo processore deve dare risultati e stampe ... :) */
	if(me==0){
		// Stampa del vettore che contiene i dati da sommare, se sono meno di 100 
		if (n<100)
		{
			for (i=0; i<n; i++)
			{
				printf("Elemento %d del vettore =%f \n",i,*(vett+i));
			}
		}

		printf("\nLa somma e' : %f\n", somma);
		printf("\nProcessori impegnati: %d\n", nproc);

		printf("\nTempo di esecuzione *Somma in Parallelo (2 strat)*: %.9f\n",T_fine);
		//printf("\nMPI_Reduce max time: %f\n",T_max);

		/*
		// debug: somma seriale
		double sum=0;
		double T_fine_serial;
		T_inizio=MPI_Wtime();
		for(i=0;i<n;i++){
		sum += vett[i];
		}
		T_fine_serial=MPI_Wtime()-T_inizio;
		printf("\n# --- DEBUG --- #\n");
		printf("Somma DI CONTROLLO seriale: %f\n",sum);
		//---
		
		double SpeedUP=0.0;
		SpeedUP= T_fine_serial/T_fine;
		printf("\n# SpeedUP: %f\n",SpeedUP);
		*/
	}

	MPI_Finalize();

}//fine programma