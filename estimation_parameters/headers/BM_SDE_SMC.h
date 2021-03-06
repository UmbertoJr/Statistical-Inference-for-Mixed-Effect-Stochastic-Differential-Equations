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
		static const int nStateVars = 1, K = 50;   //questa parte va inizializzata a seconda del modello
		int Ji, Subj_i;
		float Delta = 0.5;
		double SumWeigths, pYgivenX;
		double *MargDistrY = new double;
		double  Xj_1, Xj;
		my_matrix* PARTICLEXS_0[K], *PARTICLEXS_1[K],  *ptr_PART[K];
		double PARTICLEWS[K];
		std::mt19937 rng;
		double * X_vec, * Tzeta;
		
		
		// ho fissato la dimensione della matrice per essere pi� rapido
		
		SMC(int& subj) {
			srand(time(0)); // init seed for random
			rng.seed(std::random_device()());
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
		
		double**  run(double* tzeta) {
			Tzeta = tzeta;
			// ho saltato la funzione BM_SDE_ModelInitialDensityXsample(tzetai);
			// nStateVars = sizeof(Xprototype) / sizeof(*Xprototype);
			
			
			//				Computation for time 0
			SumWeigths = 0;
			*MargDistrY = 1;
			for (int k = 0; k < K; ++k) {
				Xj = 1e-7;  // praticamente non usato BM_SDE_ModelInitialDensityXsample(tzetai);   Xprototype
				
				(*ptr_PART[k]).insert_value_in_row_col(Xj, 0, 0);

				//double px0 = 1;  //funzione praticamente inutile BM_SDE_ModelInitialDensityXevaluate

				//cout << Subj_i + 1 << "\t" << k << "\t" << pYgivenX << endl;  // DEBUG
				// qXgivenYX_1 non lo considero visto che 1
				PARTICLEWS[k] = ModelCondDensityYgivenXevaluate(my_data::global_data.TIMEi[Subj_i][0], Xj, Tzeta[3]);
				SumWeigths += PARTICLEWS[k];
			}

			*MargDistrY *= SumWeigths / K;
			//cout << Subj_i + 1 << "\t" << MargDistrY << endl;

			//               Computation for following times
			for (int time = 1; time < my_data::global_data.len_subj[Subj_i]; ++time) {
				
				SumWeigths = 0;
				//cout << Subj_i + 1 << "\t" << time << endl;
				for (int k = 0; k < K; ++k) {
					Xj_1 = (*ptr_PART[k]).get(0,time-1);
					//cout <<"\t" << Xj_1 << " -->";

					double tj = my_data::global_data.TIMEi[Subj_i][time], tj_1 = my_data::global_data.TIMEi[Subj_i][time - 1];
					Xj = ModelTransitionDensityXsample(Xj_1, tj , tj_1, Tzeta);
					(*ptr_PART[k]).insert_value_in_row_col(Xj, 0, time);
					

					pYgivenX = ModelCondDensityYgivenXevaluate(my_data::global_data.YOBSi[Subj_i][time], Xj, Tzeta[3]);
					//cout << Xj << " \twith prob --> "<< pYgivenX << "\t";		// DEBUG
					PARTICLEWS[k] = pYgivenX;
					SumWeigths += PARTICLEWS[k];
				}
				//cout << SumWeigths << endl;
				*MargDistrY *= SumWeigths / K;
				



				BM_mnrnd(SumWeigths, time);   // QUESTA FUNZIONE SI OCCUPA DI FARE IL SAMPLING TRA LE PARTICELLE E COPIARE IL PERCORSO SAMPLLATO
				// NELLE PARTICELLE FIGLIE, SICCOME AL PRIMO GIRO SONO TUTTE EQUIPROBABILI LO FACCIO ALLA FINE

				//cout << time << endl;
				
			}
			//cout << Subj_i + 1 << "\t" << MargDistrY << endl;

			X_vec = sample_X(SumWeigths, Ji);
			double ** ptr_ret = new double*[2]; 
			ptr_ret[0] = X_vec; ptr_ret[1] = MargDistrY;
			// << "pointer to ret " << ptr_ret << " pointer to 0 " << ptr_ret[0] << " -> "<< (ptr_ret[0])[3] <<
				//" pointer to 1 " << ptr_ret[1] << " -> " << ptr_ret[1][3] << endl;
			return ptr_ret;

		}
			
		double* sample_X(double& sum, int& len_s) {
			// questo passaggio � leggermente complicato devo spiegarlo bene
			double* X  = new double[len_s];
			int j = my_random_ptr(sum);
			for (int t=0; t < len_s; ++t){
				//cout << "random pointer : " << j <<  " pointer # " << PARTICLEXS_0[j] << endl; // DEBUG
				// qui si copiano gli elementi dalle matrici indicate da ptr_PART alle matrici indicate da PARTICLEXS_1
				X[t] = (*ptr_PART[j]).get(0, t);
				//cout << X[t] << "\t";
			}
			//cout << endl;
			return X;
		}


		double ModelCondDensityYgivenXevaluate(double& Y, double& X, double& sigma) {
			if (sigma < 0) {
				std::cout << "errore sigma non pu� essere inferiore a 0" << std::endl;
				system("pause");
				return 0.0f;
			}
			// cout << Y-X << "\t" << sigma << endl ;   // DEBUG
			double pdf_gaussian = (1 / (sigma*std::sqrt(2 * M_PI))) * exp(-0.5 * std::pow((X- Y)/sigma, 2.0));
			if (pdf_gaussian < 1e-100) {
				pdf_gaussian = 1e-100;
			}


			return pdf_gaussian;
		}

		double ModelTransitionDensityXsample(double& xXj_1, double& Tj, double& Tj_1, double* tzeta) {
			double deltatime = Tj - Tj_1;
			double uj =	 EulerStep(xXj_1, Tj_1, Tj, tzeta);
			//cout <<"\t" << xXj_1 << " -->" << uj << "\t";
			//cout << "fuori da euler step" << endl;
			double sj = tzeta[2] * uj * sqrt(deltatime);
			std::normal_distribution<double> distribution(0, 1.0);
			double number = distribution(rng);  // ho cambiato il generator   distribution(generator)
			//cout << number;
			double xij = uj + sj * number;
			if (xij < 1e-100) {
				xij = ModelTransitionDensityXsample(xXj_1, Tj, Tj_1, tzeta);
			}
			return xij;
		}

		double EulerStep(double& xij_1, double&  timeij_1, double&  timeij, double* tzeta) {  
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

		void BM_mnrnd(double& sum,  int& time) {
			// questo passaggio � leggermente complicato devo spiegarlo bene
			if (time % 2 == 1) {
				//cout << "si va dai pointers 0 a quelli 1" << endl;
				// si va dai pointers 0 a quelli 1
				for (int k = 0; k < K; ++k) {
					int j = my_random_ptr(sum);
					//cout << "random pointer : " << j <<  " pointer # " << PARTICLEXS_0[j] << endl; // DEBUG
					// qui si copiano gli elementi dalle matrici indicate da ptr_PART alle matrici indicate da PARTICLEXS_1
					(*PARTICLEXS_1[k]).copy_from((*PARTICLEXS_0[j]).M);
					ptr_PART[k] = PARTICLEXS_1[k];
				}

				

			}
			else {
				//cout << "si va dai pointers 1 a quelli 0" << endl;
				// si va dai pointers 1 a quelli 0
				for (int k = 0; k < K; ++k) {
					int j = my_random_ptr(sum);
					//cout << "random pointer : " << j << " pointer # " << PARTICLEXS_1[j] << endl;  // DEBUG
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
			double s=0;
			/*for (int i = 0; i < K; i++) {
				s += PARTICLEWS[i];
			}*/
			//cout << "should never get here  "<< sum  <<" s : "<< s << " f : " << f << endl;
			return K - 1;
			//return my_random_ptr(sum);
		}

	};
	
}