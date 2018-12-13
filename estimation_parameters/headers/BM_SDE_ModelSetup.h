#pragma once
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <dlib/optimization.h>
# define M_PI           3.14159265358979323846  /* pi */
typedef dlib::matrix<double, 0, 1> column_vector;

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
		// questo passaggio deve passare anche nel PMCMC ???      /////////##########     OCCHIO 
		tzeta[0] = phii_sampled[0]; tzeta[1] = phii_sampled[1];
	}

	


	void BM_SDE_ModelPhi2Tzeta(double * phii_sampled, double * tzeta_ptr) {
		tzeta_ptr[0] = phii_sampled[0]; tzeta_ptr[1] = phii_sampled[1];
	}
	void ModelPhi2Tzeta( int &subj, int & m) {
		tzeta[0] = value_in_PHI(m, subj, 0);  tzeta[1] = value_in_PHI(m, subj, 1);
	}
	

	void ModelTzeta2Phi(double* tzeta) {
		phi[0] = tzeta[0];
		phi[1] = tzeta[1];
	};

	void ModelTzeta2Phi(double* tzeta, par_mom* par) {
		*(par->phi) = tzeta[0];
		*(par->phi + 1) = tzeta[1];
	};

	//######## funzioni utili per salvare array in matric mcmc
	double PHI[my_data::global_data.sc_M][my_data::global_data.sc_nSubjects][len_mu];
	double*ptr_phi = &PHI[0][0][0];

	double XALL[my_data::global_data.sc_M][my_data::global_data.sc_nSubjects][my_data::global_data.sc_NMAXOBS];
	double*ptr_xall = &XALL[0][0][0];


	void copy_phi_to_PHIALL(double* arr, int a, int b, int c) {
		//int start = a * my_data::global_data.sc_nSubjects*parametri.len_mu + b * parametri.len_mu + c;
		for (int i = 0; i < parametri.len_mu; i++) {
			PHI[a][b][i] = arr[i];
			//std::cout << "pointer phi iteration :" << a << "\t subj : " << b << " \t is :" << ptr_phi + start + i << " \t valore inserito : "<<arr[i]<< "\n";
		}

	}

	double value_in_PHI(int a, int b, int c) {
		return PHI[a][b][c];
	}

	void copy_newX_to_XALL(double* arr, int a, int b, int c) {
		int start = a * my_data::global_data.sc_nSubjects*my_data::global_data.sc_NMAXOBS + b * my_data::global_data.sc_NMAXOBS + c;
		for (int i = 0; i < my_data::global_data.sc_NMAXOBS; i++) {
			*(ptr_xall + start + i) = arr[i];
		}
	}

	double value_in_XALL(int &a, int& b, int c) {
		return XALL[a][b][c];
	}

	// ###### queste funzioni servono per implementare il SAEM ############
	double BM_SDE_Compute_likelihood(column_vector &theta, int m  ) {
		double logpY = 0;
		double pYijgivenXij = 0;
		for (int subj = 0; subj < my_data::global_data.sc_nSubjects; ++subj) {
			ModelTheta2Tzeta(theta);
			double *tzetai = ModelPhi2Tzeta_opt(subj, m);  // aggiorna tzeta con i "phi" semplati per il soggetto "subj" durante l'iterazione "m"
			
			//double *phii = PHI[m][subj];
			//std::cout << "iterazione : " << m << "\t soggetto : " << subj << " \t phi pointer : " << phii << " \t with value : "<< *phii<< "\n";

			double pXijgivenPhii = 1, pXigivenPhii = 1, logpYigivenXi, logpXigivenPhii, logpPhii;

			double x_value = value_in_XALL(m, subj, 0);
			double pYigivenXi = ModelCondDensityYgivenXevaluate(my_data::global_data.YOBSi[subj][0], x_value);


			for (int time_ = 1; time_ < my_data::global_data.len_subj[subj]; ++time_) {
				x_value = value_in_XALL(m, subj, time_);
				pYijgivenXij = ModelCondDensityYgivenXevaluate(my_data::global_data.YOBSi[subj][time_], x_value);
				pYigivenXi = pYigivenXi * pYijgivenXij;
				//std::cout << pYijgivenXij << "\t " << pYigivenXi << std::endl;

				pXijgivenPhii = ModelTransitionDensityXevaluate(m, subj, time_ , parametri.tzeta);
				pXigivenPhii = pXigivenPhii * pXijgivenPhii;
				//std::cout << "pXigivenPhii = " << pXigivenPhii << "\t pYigivenXi = " << pYigivenXi << "\n";
			}
			if (pYigivenXi == 0) {
				pYigivenXi = 1.e-300;
			}
			if (pXigivenPhii == 0) {
				pXigivenPhii = 1.e-300;
			}
			logpYigivenXi = log(pYigivenXi);
			logpXigivenPhii = log(pXigivenPhii);
		
			logpPhii = log(ModelAprioriDensityPhievaluate( m, subj )); // occhiooooo phii deve essere assolutamente staccato

			logpY +=  (logpYigivenXi + logpXigivenPhii + logpPhii);
			//std::cout << "logpYigivenXi = " << logpYigivenXi << "\t logpXigivenPhii = " << logpXigivenPhii << "\t logpPhii = " << logpPhii << "\t logpY  = " << logpY << std::endl;
			delete[] tzetai;
		}
		return logpY;
	}

	double PreviousQm(column_vector & theta, int & m, double & lambda) {
		double logPY = BM_SDE_Compute_likelihood(theta, 0);
		std::cout << "ok 0 \n";
		double q_1 = logPY;
		double logPYm, q = q_1;
		for (int k = 1; k < m; ++k) {
			logPYm = BM_SDE_Compute_likelihood(theta,k);
			std::cout << "ok "<< k <<" \n";
			q = q_1 + lambda * (logPYm - q_1);
			q_1 = q;
			std::cout << "Q per iterazione " << k << " e' : " << q << std::endl;
		}
		return q;
	}



	double ModelCondDensityYgivenXevaluate(double &Yij, double& Xij) {
		double norm_pdf = 1 / ( sqrt(2 * M_PI) * tzeta[3] ) * exp( -pow( (Yij - Xij) / tzeta[3] , 2 ) / 2);
		if (norm_pdf == 0) {
			norm_pdf = 1.e-100;
		}
		return norm_pdf;
	}

	double ModelTransitionDensityXevaluate(int &m , int & subj, int & time, double * tzetai) {
		double timeij = my_data::global_data.TIMEi[subj][time], timeij_1 = my_data::global_data.TIMEi[subj][time - 1];
		double deltatime = timeij - timeij_1;
		double xij_1 = value_in_XALL(m, subj, time - 1), xij = value_in_XALL(m, subj, time);
		double beta = tzetai[2];
		double uj = EulerStep(xij_1, timeij_1, timeij, tzetai);
		double sj = beta * uj*sqrt(deltatime);
		double p = (1 / ( sqrt(2 * M_PI)* sj )  )  *  exp(-pow(((xij - uj) / sj), 2) / 2);
		if (p == 0) {
			p = 1.e-300;
		}
		return p;
	}

	double ModelAprioriDensityPhievaluate(int& m, int& subj) {
		double * phii = PHI[m][subj];
		//std::cout << "phii[0] = " <<phii[0] << "\t phii[1] = " << phii[1] << "\t mu[0] = "<< mu[0]<< " \n";
		double phi_mu[2]; phi_mu[0] = phii[0] - mu[0]; phi_mu[1] = phii[1] - mu[1];
		double foo = pow(phi_mu[0],2)* Omega[0][0] + pow(phi_mu[1], 2)* Omega[1][1];
		double q = 1 / (Omega[0][0] * Omega[1][1] * sqrt(2 * M_PI))* exp(-foo / 2);
		if (q == 0) {
			q = 1.e-300;
		}
		// std::cout << "prob = " << q<< "\t foo = "<< foo<< "\n";
		return q;
	}

	double * ModelPhi2Tzeta_opt(int & subj, int & m) {
		// questo passaggio serve per calcolare 
		double * tzetai = new double[4];
		tzetai[0] = value_in_PHI(m, subj, 0); tzetai[1] = value_in_PHI(m, subj, 1); tzetai[2] = tzeta[2]; tzetai[3] = tzetai[3];
		return tzetai;
	}

	void ModelTheta2Tzeta(column_vector &theta) {
		tzeta[0] = theta(0);
		tzeta[1] = theta(1);
		tzeta[2] = theta(4);
		tzeta[3] = theta(5);

		mu[0] = theta(0);
		mu[1] = theta(1);

		Omega[0][0] = std::pow(theta(2), 2);
		Omega[0][1] = 0;
		Omega[1][0] = 0;
		Omega[1][1] = std::pow(theta(3), 2);
	}

	double EulerStep(double& xij_1, double&  timeij_1, double&  timeij, double* tzeta) {
		double deltat = (timeij - timeij_1) / 10, dt;
		double Xhat1 = xij_1;
		double time = timeij_1;
		while (time < timeij) {
			if ((timeij - time) > deltat) {
				dt = deltat;
			}
			else {
				dt = (timeij - time);
			}
			Xhat1 += ModelDrift(Xhat1, tzeta) * dt;
			time += dt;
			//cout << time << endl;
		}
		//cout << Xhat1 << "#";
		return Xhat1;
	}

	double ModelDrift(double& X, double* tzeta) {   // non serve double& time,
		double alpha = tzeta[0], kappa = tzeta[1];
		double dXdt = alpha * log(kappa / X + 1)*X;
		return dXdt;
	}

	// ###### LEGGE E SCRIVERE PARAMETRI IN FILE TXT
	void read_parameters_from_file(std::string& namefile) {
		std::ifstream file(namefile);
		file.seekg(0);
		std::string str;
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

	void write_parameters_to_file(std::string& filename) {
		std::ofstream file(filename);
		file << "alpha\tkappa\tomega_alpha\tomega_kappa\tbeta\tsigma\n";
		file << alpha << "\t" << kappa << "\t" << omega_alpha << "\t" << omega_kappa << "\t" << beta << "\t" << sigma << "\n";
		file.close();
	};

		
} parametri ;

