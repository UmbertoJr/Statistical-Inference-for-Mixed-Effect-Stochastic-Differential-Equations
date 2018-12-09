#pragma once
#include <BM_SDE_SMC.h>

static class Pmcmc
{
	//questa classe si occupa di svoldere tutte le funzioni relative al PMCMC
public:
	int subj, R;
	double X[my_data::global_data.sc_NMAXOBS];
	

	Pmcmc() {
		subj = 0;
		R = 100;
		for (int it = 0; it < my_data::global_data.sc_NMAXOBS; it++) {
			X[it] = 0;
		}
	}
	void iterations(int& r) {
		R = r;
	}

	void run(int& Subj)
	{
		subj = Subj;
		parametri.ModelTheta2Tzeta(&parametri_tmp); parametri.ModelTzeta2Phi();   // da sistemare
		double* phi0 = parametri.phi;

		// ora bisogna creare BM_SDE_SMC
		BM_SDE_SMC::SMC smc(subj);
		smc.run();

		cout << "fuori dal smc, soggetto numero " << subj << endl;
	}
	
} pmcmc ;
