#pragma once
#include <iostream>
#include <vector>
#include <BM_SDE_PMCMC.h>
#include <my_matrix_SAEM.h>
#include <random>
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#include <BM_SDE_Stochastic_Approximation.h>

using namespace std;



namespace BM_SDE_SAEM_Optimized {
	static void run() {
		

		static const int M = my_data::global_data.sc_M;
		static const int M1 = my_data::global_data.sc_M1;
		double lambda;
		column_vector theta(6);

		// devo stare attento alla visibility di questa variabile e salvare i dati prima della fine altrimenti perdo
		// le simulazioni...
		//cout << "Omega " << parametri.Omega[0][0] << "\t" << parametri.Omega[1][1]  << "\t" << endl;

		parametri.ModelTheta2Tzeta();  // controllare queste due funzioni
		parametri.ModelTzeta2Phi();  // inizializza le variabili
	
		//cout << "Omega " << parametri.Omega[0][0] << "\t" << parametri.Omega[1][1] << "\t" << endl;
		

		int r = 100;
		pmcmc.iterations(r); // set the number of iterations inside the pmcmc object

		double **hidden_sample;
		double * ptr_arr;
		// inizializza le PHI e le X per tutti i soggetti
		for (int subj = 1; subj <= my_data::global_data.sc_nSubjects; ++subj) {
			cout << "\n\nsoggetto : " << subj << "\t" <<" lunghezza x : " << my_data::global_data.len_subj[subj-1] << "\n"; //printa i risultati per debug

			static clock_t start = clock();  // fa partire il timer 
			hidden_sample = pmcmc.run(subj, parametri.theta);  // le variabili tzeta vengono aggiornate in BM_SDE_ModelTheta2Tzeta
			static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
			static double time = (double)(end - start) / CLOCKS_PER_SEC;
			
			// devo capire cosa fare con i vettori dinamici dichiarati dento pmcmc
			parametri.BM_SDE_ModelPhi2Tzeta(hidden_sample[1]);

			//cout << "fuori dal pmcmc, soggetto numero " << subj << endl;
			
			ptr_arr = hidden_sample[1];
			cout << "phi sampled : " << ptr_arr[0] << "\t" << ptr_arr[1] << endl;
			PHI.copy_phi_to_PHIALL(ptr_arr, 0,subj-1,0 );

			ptr_arr = hidden_sample[0];  // puntatore verso i campioni X samplati nel pmcmc	
			cout << "\nX sampled : ";
			for (int it = 0; it < my_data::global_data.len_subj[subj - 1]; ++it) {
				cout << ptr_arr[it] << "\t" ;
			}
			cout << endl;
			cout << "\n######  time execution PMCMC is: " << time << "  #######" << endl << endl << endl;

			XALL.copy_newX_to_XALL(ptr_arr, 0, subj - 1, 0);
			

		}

		for (int m = 1; m < 2; ++m) { // per ora faccio una sola iterazione 
			parametri.ModelTheta2Tzeta(parametri.theta);

			for (int subj = 1; subj <= my_data::global_data.sc_nSubjects; ++subj) {
				hidden_sample = pmcmc.run(subj, parametri.theta);  
				parametri.BM_SDE_ModelPhi2Tzeta(hidden_sample[1]);

				ptr_arr = hidden_sample[1]; // puntatore verto i phi samplati nel pmcmc
				PHI.copy_phi_to_PHIALL(ptr_arr, m, subj - 1, 0);

				ptr_arr = hidden_sample[0];  // puntatore verso i campioni X samplati nel pmcmc	
				XALL.copy_newX_to_XALL(ptr_arr, m, subj - 1, 0);
			}
			
			
			saem.add_m();
			//double logPY = BM_SDE_Stochastic_Approximation(theta , m, lambda);
			//cout << "logPY per iterazione " << m << " è :" << logPY << endl;
		}
			

		
		
		


		}
	
}

