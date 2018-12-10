#pragma once
#include <iostream>
#include <vector>
#include <BM_SDE_PMCMC.h>
#include <my_matrix_SAEM.h>
#include <random>

using namespace std;



namespace BM_SDE_SAEM_Optimized {
	static void run() {
		

		static const int M = my_data::global_data.sc_M;
		static const int M1 = my_data::global_data.sc_M1;

		// devo stare attento alla visibility di questa variabile e salvare i dati prima della fine altrimenti perdo
		// le simulazioni...
		
		parametri.ModelTheta2Tzeta();  // controllare queste due funzioni
		parametri.ModelTzeta2Phi();
	
		//double* ptr_arr; // questo pointer serve per copiare i dati in PHIALL e XALL

		int r = 100;
		pmcmc.iterations(r); // set the number of iterations inside the pmcmc object
		double **hidden_sample;
		double * ptr_arr;
		// inizializza le PHI e le X per tutti i soggetti
		for (int subj = 1; subj <= my_data::global_data.sc_nSubjects; ++subj) {
			cout << subj << "\t" << my_data::global_data.len_subj[subj-1] << "\n"; //printa i risultati per debug

			//static clock_t start = clock();  // fa partire il timer 
			hidden_sample = pmcmc.run(subj);  // le variabili tzeta vengono aggiornate in BM_SDE_ModelTheta2Tzeta
			/*static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
			static double time = (double)(end - start) / CLOCKS_PER_SEC;
			cout << "\n\n######  time execution PMCMC is: " << time << "  #######" << endl;
			*/
			// devo capire cosa fare con i vettori dinamici dichiarati dento pmcmc
			parametri.BM_SDE_ModelPhi2Tzeta(hidden_sample[1]);

			cout << "fuori dal pmcmc, soggetto numero " << subj << endl;
			
			ptr_arr = parametri.phi;
			PHI.copy_phi_to_PHIALL(ptr_arr, 0,subj-1,0 );

			ptr_arr = hidden_sample[0];  // puntatore verso i campioni X samplati nel pmcmc	
			XALL.copy_newX_to_XALL(ptr_arr, 0, subj - 1, 0);
			
			
			
			
		}

		//cout << PHI.value_in(0,3,1) << endl;   //solo per fare debug

		}
	
}

