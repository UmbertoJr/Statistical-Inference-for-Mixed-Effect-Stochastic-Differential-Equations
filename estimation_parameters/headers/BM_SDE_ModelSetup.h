#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

static struct par_mom {
	double theta[6], tzeta[4], mu[2], Omega[4], phi[2];

} parametri_tmp;

static class parameters {
	// da fare funzione che legge e scrive i parametri del modello
	double alpha; 
	double kappa; 
	double omega_alpha; 
	double omega_kappa; 
	double beta;
	double sigma; 
		

public:

	double tzeta[4];
	double mu[2];
	static const int len_mu = (sizeof(mu) / sizeof(*mu));;
	double Omega[2][2];
	double theta[6];
	double PARMIN[6];
	double PARMAX[6];
	double phi[2];

	

	
	void initialization() {
		tzeta[0] = alpha;
		tzeta[1] = kappa;
		tzeta[2] = beta;
		tzeta[3] = sigma;
			

		mu[0] = alpha;
		mu[1] = kappa;
			

		Omega[0][0]  = pow(omega_alpha, 2);
		Omega[0][1] = 0;
		Omega[1][0] = 0;
		Omega[1][1] = pow(omega_kappa, 2);
			

		theta[0] = alpha;       PARMIN[0] = 0.01;   PARMAX[0] = 1;
		theta[1] = kappa;       PARMIN[1] = 1E-6;   PARMAX[1] = 1E-2;
		theta[2] = omega_alpha; PARMIN[2] = 0.0001; PARMAX[2] = 0.05;
		theta[3] = omega_kappa; PARMIN[3] = 1E-5;   PARMAX[3] = 1E-2;
		theta[4] = beta;        PARMIN[4] = 0;      PARMAX[4] = 1;
		theta[5] = sigma;       PARMIN[5] = 0;      PARMAX[5] = 1E-2;

	};

	void ModelTheta2Tzeta() {
		tzeta[0] = theta[0];
		tzeta[1] = theta[1];
		tzeta[2] = theta[4];
		tzeta[3] = theta[5];

		mu[0] = theta[0];
		mu[1] = theta[1];

		Omega[0][0] = pow(theta[2], 2);
		Omega[0][1] = 0;
		Omega[1][0] = 0;
		Omega[1][1] = pow(theta[3], 2);

	};

	void ModelTheta2Tzeta(double* theta) {
		tzeta[0] = theta[0];
		tzeta[1] = theta[1];
		tzeta[2] = theta[4];
		tzeta[3] = theta[5];

		mu[0] = theta[0];
		mu[1] = theta[1];

		Omega[0][0] = pow(theta[2], 2);
		Omega[0][1] = 0;
		Omega[1][0] = 0;
		Omega[1][1] = pow(theta[3], 2);

	};

	void ModelTheta2Tzeta(par_mom* par) {
		*(par->tzeta) = theta[0];
		*(par->tzeta + 1) = theta[1];
		*(par->tzeta + 2) = theta[4];
		*(par->tzeta + 3) = theta[5];

		*(par->mu) = theta[0];
		*(par->mu + 1) = theta[1];

		*(par->Omega) = pow(theta[2], 2);
		*(par->Omega + 1) = 0;
		*(par->Omega + 2) = 0;
		*(par->Omega + 3) = pow(theta[3], 2);

	};

	void ModelTzeta2Phi() {
		phi[0] = tzeta[0];
		phi[1] = tzeta[1];
	};

	void BM_SDE_ModelPhi2Tzeta(double * phii_sampled) {
		// questo passaggio deve passare anche nel PMCMC ???
		tzeta[0] = phii_sampled[0]; tzeta[1] = phii_sampled[1];
	}

	void BM_SDE_ModelPhi2Tzeta(double * phii_sampled, double * tzeta_ptr) {
		tzeta_ptr[0] = phii_sampled[0]; tzeta_ptr[1] = phii_sampled[1];
	}

	void ModelTzeta2Phi(double* tzeta) {
		phi[0] = tzeta[0];
		phi[1] = tzeta[1];
	};

	void ModelTzeta2Phi(double* tzeta, par_mom* par) {
		*(par->phi) = tzeta[0];
		*(par->phi + 1) = tzeta[1];
	};

	void read_parameters_from_file(string& namefile) {
		ifstream file(namefile);
		file.seekg(0);
		string str;
		std::getline(file, str);	//salta gli headers

		char split_char = '\t';
		std::vector<std::string> tokens;
		double parameters[6];
		int i = 0;
		for (std::string each; std::getline(file, each, split_char); tokens.push_back(each)) {
			parameters[i] = std::stod(each);
			++i;
			//std::cout << i << "\t" << parameters[i - 1] << "\n";
		};
		alpha = parameters[0];
		kappa = parameters[1];
		omega_alpha = parameters[2];
		omega_kappa = parameters[3];
		beta = parameters[4];
		sigma = parameters[5];
		// devo fare anche Parmin e Parmax

		file.close();
	};

	void write_parameters_to_file(string& filename) {
		ofstream file(filename);
		file << "alpha\tkappa\tomega_alpha\tomega_kappa\tbeta\tsigma\n";
		file << alpha << "\t" << kappa << "\t" << omega_alpha << "\t" << omega_kappa << "\t" << beta << "\t" << sigma << "\n";
		file.close();
	};



		
} parametri ;

