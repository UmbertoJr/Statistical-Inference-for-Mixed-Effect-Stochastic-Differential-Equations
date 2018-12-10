#pragma once
#include <BM_SDE_SMC.h>
#include <random>
#include <math.h>


static class Pmcmc
{
	std::mt19937 rng;
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
		
		static clock_t start = clock();  // fa partire il timer 
										 
		// ora bisogna creare BM_SDE_SMC
		///*

		double * tzetai = parametri.tzeta; // inizializa i parametri a quelli iniziali
		BM_SDE_SMC::SMC smc(subj);
		X = smc.run(tzetai);  // controllare che i parametri vengano aggiornati
		// OCCHIO MI DEVO RICORDARE DI CANCELLARE X POI PERCHE Cè UN VETTORE DINAMICO e anche i due puntatori al suo interno
		
		
		
		/*
		for (int t = 0; t < my_data::global_data.len_subj[subj - 1]; ++t) {
			cout << X[0][t] << "\t";
		}
		cout << endl;
		*/ 
		///*
		static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
		static double time = (double)(end - start) / CLOCKS_PER_SEC;
		cout << "\n\n######  time execution SMC is: " << time << "  #######" << endl;
		//*/

		//cout << "fuori dal smc, soggetto numero " << subj << endl;

		double * phic = new double[2];
		double * tzetac = new double[4];
		double ** Xc;
		tzetac[2]= tzetai[2]; tzetac[3] = tzetai[3]; 
		for (int r = 0; r < R; ++r) {
			BM_SDE_ModelProposalDensityPhisample(phi0, phic);  // from phi0 to phic
			// cout << "phi0 : " << phi0[0] << "  phic : " << phic[0] << endl; // just DEBUG
			parametri.BM_SDE_ModelPhi2Tzeta(phic, tzetac);
			//cout << " tzetai : " << tzetai[1] << " tzetac : " << tzetac[1] << endl;
			Xc = smc.run(tzetac);  // controllare che i parametri vengano aggiornati
			//cout << " X[0][4] : " << X[0][4] << " Xc[0][4] : " << Xc[0][4] << endl;

		}





		hidden_vec = new double* [2];
		hidden_vec[0] = X[0]; hidden_vec[1] = phi0;
		return  hidden_vec;
	}
	


	void BM_SDE_ModelProposalDensityPhisample(double * phi , double * save_phi) {
		// siccome laa variabile multinormale a covarianza zero samplo due volte da una normale...
		// questo nel caso generale non si può fare
		rng.seed(std::random_device()());
		std::normal_distribution<double> distribution(0, 1.0);
		for (int ind = 0; ind < 2; ++ind) {
			save_phi[ind] = phi[ind] + 0.01* pow(phi[ind], 2) * distribution(rng);
		}
	}
} pmcmc ;
