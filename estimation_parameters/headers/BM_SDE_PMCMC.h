#pragma once
#include <BM_SDE_SMC.h>
#include <random>
#include <math.h>
#include <cmath>


static class Pmcmc
{
	std::mt19937 rng;
	//questa classe si occupa di svoldere tutte le funzioni relative al PMCMC
public:
	int subj, R;
	double** X;
	double **hidden_vec;
	double *theta, *tzeta, *phi;
	double * phic = new double[2];
	double * tzetac = new double[4];
	

	Pmcmc() {
		subj = 0;
		R = 100;
	}
	
	~Pmcmc() {
		// qui libero la memoria dinamica dai risultati creati durante il SMC 
		//  delete[] X[0];  questa qui non la devo eliminare perche lo uso anche nel passaggio SAEM
		delete X[1];
		delete[] X;
		delete[] phic;
		delete[] tzetac;
	}

	void iterations(int& r) {
	
		R = r;
	}

	double** run(int& Subj, double * thetaI)
	{
		rng.seed(std::random_device()());
		subj = Subj;
		parametri.ModelTheta2Tzeta(thetaI); parametri.ModelTzeta2Phi();   // da sistemare
		double* phi0 = parametri.phi;
		
		/*  // just for DEBUG
		for (int k = 0; k < 6; ++k) {
			cout << parametri.tzeta[k]<< endl;
		}
		*/
		
		//static clock_t start = clock();  // fa partire il timer 
										 
		// ora bisogna creare BM_SDE_SMC
		

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
		/*
		static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
		static double time = (double)(end - start) / CLOCKS_PER_SEC;
		cout << "\n\n######  time execution SMC is: " << time << "  #######" << endl;
		*/

		//cout << "fuori dal smc, soggetto numero " << subj << endl;

		
		double ** Xc;
		tzetac[2]= tzetai[2]; tzetac[3] = tzetai[3]; 
		double acpt_rate = 0;
		for (int r = 0; r < R; ++r) {
			BM_SDE_ModelProposalDensityPhisample(phi0, phic);  // from phi0 to phic
			//cout << "phi0 : " << phi0[0] << "  phic : " << phic[0] << endl; // just DEBUG
			//cout << "phi0 : " << phi0[1] << "  phic : " << phic[1] << endl; // just DEBUG
			parametri.BM_SDE_ModelPhi2Tzeta(phic, tzetac);
			//cout << " tzetai : " << tzetai[1] << " tzetac : " << tzetac[1] << endl;
			Xc = smc.run(tzetac);  // controllare che i parametri vengano aggiornati
			//cout << " X[0][4] : " << X[0][4] << " Xc[0][4] : " << Xc[0][4] << endl;

			//compute probability of acceptance
			// ho utilizzato la loglikelihood
			double num_alpha = BM_SDE_ModelProposalDensityPhievaluate(phi0, phic) * (*Xc[1]) * BM_SDE_ModelAprioriDensityPhievaluate(phic);
			//cout << log(*Xc[1]) << "\t" << num_alpha;
			// ho tolto BM_SDE_ModelProposalDensityPhievaluate(phi0, phic) tanto è simmetrico

			double den_alpha = BM_SDE_ModelProposalDensityPhievaluate(phic ,phi0) * (*X[1]) * BM_SDE_ModelAprioriDensityPhievaluate(phi0);
			//cout << " \t" << log(*X[1]) << "\t" << den_alpha << "\t" << num_alpha- den_alpha << "\t" << exp(num_alpha - den_alpha) <<  endl;
			//alpha = min(1, num_alpha / den_alpha);
			double alpha = std::min( 1., num_alpha /den_alpha);
			double U = (double)rand() / RAND_MAX;
			if (U <= alpha) {
				acpt_rate += 1;
				//cout << "\n#########passato con U : " << U << " e alpha : " << alpha << endl;
				//cout << "den alpha : " << den_alpha <<endl;
				//cout << "num alpha : " << num_alpha << endl;
				//cout << endl;
				for (int it_ = 0; it_ < my_data::global_data.len_subj[subj - 1]; ++it_) {
					X[0][it_] = Xc[0][it_];
				}
				*X[1] = *Xc[1];
				phi0[0] = phic[0]; phi0[1] = phic[1];
				
			}
		}

		//std::cout << "acceptance rate : " << acpt_rate / R << std::endl << std::endl ;
		hidden_vec = new double* [2];
		hidden_vec[0] = X[0]; hidden_vec[1] = phi0;
		return  hidden_vec;
	}
	
	double BM_SDE_ModelAprioriDensityPhievaluate(double * phi) {
		// il denominatore non lo considero tanto non varia nel PMCMC
		double * mu = parametri.mu;
		double * Omega[2] = { parametri.Omega[0],parametri.Omega[1] };
		//cout << "Omega " << Omega[0][0] << "\t" << Omega[1][1] << "\t" << endl;
		double phi_mu[2];
		phi_mu[0] = phi[0] - mu[0]; phi_mu[1] = phi[1] - mu[1];
		double foo = pow(phi_mu[0], 2)*(1 / Omega[0][0]) + pow(phi_mu[1], 2)*(1 / Omega[1][1]);
		double q = exp( -foo / 2);
		//cout << "density evaluation "<< phi_mu[0] << "\t" << foo << "\t" << q << endl;
		//cout << "Omega " << Omega[0] << "\t" << Omega[3] << "\t" << q << endl;
		return q;
	}

	double BM_SDE_ModelProposalDensityPhievaluate(double * phi1,  double * phi2) {
		double variance[2];
		variance[0] = 0.01*pow(phi2[0], 2); variance[1] = 0.01*pow(phi2[1], 2);
		for (int j = 0; j < 2; ++j) {
			if (variance[j] < 1.e-12) {
				variance[j] = 1e-12;
			}
		}
		double delta_phi[2]; delta_phi[0] = phi1[0] - phi2[0]; delta_phi[1] = phi1[1] - phi2[1];
		double foo = (1 / (sqrt(2 * M_PI*variance[0] * variance[1])));
		double exp_part = pow(delta_phi[0], 2) * (1 / variance[0]) + pow(delta_phi[1], 2) * (1 / variance[1]);
		double prob = foo * exp(-exp_part / 2);
		//cout << prob << endl;
		return prob;
	}


	void BM_SDE_ModelProposalDensityPhisample(double * phi , double * save_phi) {
		// siccome laa variabile multinormale a covarianza zero samplo due volte da una normale...
		// questo nel caso generale non si può fare
		std::normal_distribution<double> distribution(0, 1.0);
		for (int ind = 0; ind < 2; ++ind) {
			double variance = 0.01* pow(phi[ind], 2);
			if (variance < 1e-12) {
				variance = 1e-12;
			}
			save_phi[ind] = phi[ind] + variance * distribution(rng);
		}
	}
} pmcmc ;
