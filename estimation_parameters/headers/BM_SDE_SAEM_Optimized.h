#pragma once
#include <iostream>
#include <vector>
#include <BM_SDE_PMCMC.h>
#include <my_matrix_SAEM.h>
#include <random>

using namespace std;



namespace BM_SDE_SAEM_Optimized {
	static void run() {
		// devo stare attento alla visibility di questa variabile e salvare i dati prima della fine altrimenti perdo
		// le simulazioni...
		
		parametri.ModelTheta2Tzeta();
		parametri.ModelTzeta2Phi();
	
		//double* ptr_arr; // questo pointer serve per copiare i dati in PHIALL e XALL

		int r = 100;
		pmcmc.iterations(r); // set the number of iterations inside the pmcmc object

		// inizializza le PHI e le X per tutti i soggetti
		for (int subj = 1; subj <= my_data::global_data.sc_nSubjects; ++subj) {
			cout << subj << "\t" << my_data::global_data.len_subj[subj-1] << "\n"; //printa i risultati per debug
			pmcmc.run(subj);

			cout << "fuori dal pmcmc, soggetto numero " << subj << endl;
			/*
			ptr_arr = parametri.phi;
			PHI.copy_phi_to_PHIALL(ptr_arr, 0,subj-1,0 );

			ptr_arr = pmcmc.X;
			XALL.copy_newX_to_XALL(ptr_arr, 0, subj - 1, 0);
			
			*/
			
			
		}

		//cout << PHI.value_in(0,3,1) << endl;   //solo per fare debug

		}
	
}

