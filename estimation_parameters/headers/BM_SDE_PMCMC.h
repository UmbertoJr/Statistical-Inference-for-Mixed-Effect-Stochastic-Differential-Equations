#pragma once
#include <BM_SDE_SMC.h>

static class Pmcmc
{
	//questa classe si occupa di svoldere tutte le funzioni relative al PMCMC
public:
	int subj, R;
	double** X;
	double **hidden_vec;
	

	Pmcmc() {
		subj = 0;
		R = 100;
	}
	~Pmcmc() {
		// qui libero la memoria dinamica dai risultati creati durante il SMC 
		//  delete[] X[0];  questa qui non la devo eliminare perche lo uso anche nel passaggio SAEM
		delete X[1];
		delete[] X;
	}

	void iterations(int& r) {
	
		R = r;
	}

	double** run(int& Subj)
	{
		subj = Subj;
		parametri.ModelTheta2Tzeta(&parametri_tmp); parametri.ModelTzeta2Phi();   // da sistemare
		double* phi0 = parametri.phi;
		
		/*  // just for DEBUG
		for (int k = 0; k < 6; ++k) {
			cout << parametri.tzeta[k]<< endl;
		}
		*/
		
		//static clock_t start = clock();  // fa partire il timer 
										 
		// ora bisogna creare BM_SDE_SMC
		///*

		double * tzetai = parametri.tzeta; // inizializa i parametri a quelli iniziali
		BM_SDE_SMC::SMC smc(subj, tzetai);
		X = smc.run();  // controllare che i parametri vengano aggiornati
		// OCCHIO MI DEVO RICORDARE DI CANCELLARE X POI PERCHE Cè UN VETTORE DINAMICO e anche i due puntatori al suo interno
		
		
		
		/*
		for (int t = 0; t < my_data::global_data.len_subj[subj - 1]; ++t) {
			cout << X[0][t] << "\t";
		}
		cout << endl;
		*/ 
		/*
		static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
		static double time = (double)(end - start) / CLOCKS_PER_SEC;
		cout << "\n\n######  time execution SMC is: " << time << "  #######" << endl;
		*/

		//cout << "fuori dal smc, soggetto numero " << subj << endl;

		hidden_vec = new double* [2];
		hidden_vec[0] = X[0]; hidden_vec[1] = phi0;
		return  hidden_vec;
	}
	
} pmcmc ;
