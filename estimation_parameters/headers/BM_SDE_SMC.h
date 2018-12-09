#pragma once
# define M_PI           3.14159265358979323846  /* pi */
#include <random>


namespace BM_SDE_SMC {
	struct my_matrix {
		double* M;
		int t_r, t_c;
		my_matrix(const int& t_row, int& t_col) {
			M = new double[t_row * t_col]; 
			t_r = t_row; t_c = t_col;
		}
		
		void copy_from(double* arr){
			for (int itera = 0; itera < t_r*t_c; ++itera) {
				M[itera] = arr[itera];
			}
			
		}
		
		void insert_value_in_row_col(double& val, int& row, int& col) {
			M[row*t_c + col] = val;
		}

		void insert_value_in_row_col(double& val, int row, int col) {
			M[row*t_c + col] = val;
		}

		double get(int& row, int& col) {
			return M[row*t_c + col];
		}
		double get(int row, int col) {
			return M[row*t_c + col];
		}
	};

	
	struct SMC {
		//double Xprototype = 1E-07;
		static const int nStateVars = 1, K = 5;   //questa parte va inizializzata a seconda del modello
		int Ji, Subj_i;
		float Delta = 0.5;
		double MargDistrY, SumWeigths, pYgivenX;
		double  Xj_1, Xj;
		my_matrix* PARTICLEXS_0[K], *PARTICLEXS_1[K],  *ptr_PART[K];
		double PARTICLEWS[K];
		std::default_random_engine generator;
		
		
		// ho fissato la dimensione della matrice per essere più rapido
		
		SMC(int& subj) {
			Subj_i = subj -1;
			Ji = my_data::global_data.len_subj[Subj_i];		
			
			for (int k = 0; k < K; ++k) {
				PARTICLEXS_0[k] = new my_matrix(nStateVars, Ji );
				PARTICLEXS_1[k] = new my_matrix(nStateVars, Ji);
				ptr_PART[k] = PARTICLEXS_0[k];
			}
		}


		~SMC() {
			for (int k = 0; k < K; ++k) {
				//cout << "destructor called" << endl;
				delete[] PARTICLEXS_0[k]->M;
				delete[] PARTICLEXS_1[k]->M;
				delete PARTICLEXS_0[k];
				delete PARTICLEXS_1[k];
			}
		}
		
		void run() {
		
			// ho saltato la funzione BM_SDE_ModelInitialDensityXsample(tzetai);
			// nStateVars = sizeof(Xprototype) / sizeof(*Xprototype);
			
			
			//				Computation for time 0
			SumWeigths = 0;
			MargDistrY = 1;
			for (int k = 0; k < K; ++k) {
				Xj = 1e-7;  // praticamente non usato BM_SDE_ModelInitialDensityXsample(tzetai);   Xprototype
				
				(*ptr_PART[k]).insert_value_in_row_col(Xj, 0, 0);

				//double px0 = 1;  //funzione praticamente inutile BM_SDE_ModelInitialDensityXevaluate
				pYgivenX = BM_SDE_ModelCondDensityYgivenXevaluate(my_data::global_data.TIMEi[Subj_i][0], Xj, parametri_tmp.tzeta[3]);

				//cout << Subj_i + 1 << "\t" << k << "\t" << pYgivenX << endl;  // DEBUG
				// qXgivenYX_1 non lo considero visto che 1
				PARTICLEWS[k] = pYgivenX;
				SumWeigths += PARTICLEWS[k];
			}

			MargDistrY *= SumWeigths / K;
			cout << Subj_i + 1 << "\t" << MargDistrY << endl;

			//               Computation for following times
			for (int time = 1; time < my_data::global_data.len_subj[Subj_i]; ++time) {
				
				SumWeigths = 0;
				cout << Subj_i + 1 << "\t" << time << endl;
				for (int k = 0; k < K; ++k) {
					Xj_1 = (*ptr_PART[k]).get(0,time-1);
					//cout <<"\t" << Xj_1 << " -->";

					double tj = my_data::global_data.TIMEi[Subj_i][time], tj_1 = my_data::global_data.TIMEi[Subj_i][time - 1];
					double* tzetai = parametri.tzeta;  // questo passaggio deve essere aggiornato ad ogni iterazione del PMCMC
					Xj = BM_SDE_ModelTransitionDensityXsample(Xj_1, tj , tj_1, tzetai);
					(*ptr_PART[k]).insert_value_in_row_col(Xj, 0, time);
					

					pYgivenX = BM_SDE_ModelCondDensityYgivenXevaluate(my_data::global_data.YOBSi[Subj_i][time], Xj, tzetai[3]);
					//cout << Xj << " \twith prob --> "<< pYgivenX << "\t";		// DEBUG
					PARTICLEWS[k] = pYgivenX;
					SumWeigths += PARTICLEWS[k];
				}
				//cout << SumWeigths << endl;
				
				



				BM_mnrnd(SumWeigths, time);   // QUESTA FUNZIONE SI OCCUPA DI FARE IL SAMPLING TRA LE PARTICELLE E COPIARE IL PERCORSO SAMPLLATO
				// NELLE PARTICELLE FIGLIE, SICCOME AL PRIMO GIRO SONO TUTTE EQUIPROBABILI LO FACCIO ALLA FINE

				//cout << time << endl;
				
			}
			
		}
			

		double BM_SDE_ModelCondDensityYgivenXevaluate(double& Y, double& X, double& sigma) {
			if (sigma < 0) {
				cout << "errore sigma non può essere inferiore a 0" << endl;
				return 0.0f;
			}
			// cout << Y-X << "\t" << sigma << endl ;   // DEBUG
			double pdf_gaussian = (1 / (sigma*sqrt(2 * M_PI))) * exp(-0.5 * pow((X- Y)/sigma, 2.0));
			if (pdf_gaussian < 1e-100) {
				pdf_gaussian = 1e-100;
			}


			return pdf_gaussian;
		}

		double BM_SDE_ModelTransitionDensityXsample(double& xXj_1, double& Tj, double& Tj_1, double* tzeta) {
			double deltatime = Tj - Tj_1;
			double uj =	 BM_SDE_EulerStep(xXj_1, Tj_1, Tj, tzeta);
			//cout <<"\t" << xXj_1 << " -->" << uj << "\t";
			//cout << "fuori da euler step" << endl;
			double sj = tzeta[2] * uj * sqrt(deltatime);
			std::normal_distribution<double> distribution(0, 1.0);
			double number = distribution(generator);
			//cout << number;
			double xij = uj + sj * distribution(generator);
			if (xij < 1e-100) {
				xij = BM_SDE_ModelTransitionDensityXsample(xXj_1, Tj, Tj_1, tzeta);
			}
			return xij;
		}

		double BM_SDE_EulerStep(double& xij_1, double&  timeij_1, double&  timeij, double* tzeta) {  
			double deltat = (timeij- timeij_1)/10,  dt;
			double Xhat1 = xij_1;
			double time = timeij_1;
			while (time < timeij) {
				if ((timeij - time) > deltat) {
					dt = deltat;
				}
				else {
					dt = (timeij - time);
				}
				Xhat1 += BM_SDE_ModelDrift(Xhat1, tzeta) * dt;
				time += dt;
				//cout << time << endl;
			}
			//cout << Xhat1 << "#";
			return Xhat1;
		}

		double BM_SDE_ModelDrift(double& X, double* tzeta) {   // non serve double& time,
			double alpha = tzeta[0], kappa = tzeta[1];
			double dXdt = alpha * log(kappa / X + 1)*X;
			return dXdt;
		}

		void BM_mnrnd(double& sum,  int& time) {
			// questo passaggio è leggermente complicato devo spiegarlo bene
			if (time % 2 == 1) {
				cout << "si va dai pointers 0 a quelli 1" << endl;
				// si va dai pointers 0 a quelli 1
				for (int k = 0; k < K; ++k) {
					int j = my_random_ptr(sum);
					cout << "random pointer : " << j <<  " pointer # " << PARTICLEXS_0[j] << endl;
					// qui si copiano gli elementi dalle matrici indicate da ptr_PART alle matrici indicate da PARTICLEXS_1
					(*PARTICLEXS_1[k]).copy_from((*PARTICLEXS_0[j]).M);
					ptr_PART[k] = PARTICLEXS_1[k];
				}

				

			}
			else {
				cout << "si va dai pointers 1 a quelli 0" << endl;
				// si va dai pointers 1 a quelli 0
				for (int k = 0; k < K; ++k) {
					int j = my_random_ptr(sum);
					cout << "random pointer : " << j << " pointer # " << PARTICLEXS_1[j] << endl;
					// qui si copiano gli elementi dalle matrici indicate da ptr_PART alle matrici indicate da PARTICLEXS_0
					(*PARTICLEXS_0[k]).copy_from((*PARTICLEXS_1[j]).M);
					ptr_PART[k] = PARTICLEXS_0[k];
				}

			}
			
		}

		int my_random_ptr(double& sum) {
			double f = (double)rand()* sum / RAND_MAX;
			for (int i = 0; i < K; i++) {
				if (f < PARTICLEWS[i])
					return i;
				f -= PARTICLEWS[i];
			}
			cout << "should never get here"<< endl;
			return my_random_ptr(sum);
		}

	};
	
}