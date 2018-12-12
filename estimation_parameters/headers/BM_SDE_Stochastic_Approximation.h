#pragma once

typedef dlib::matrix<double, 0, 1> column_vector;

class SAEM {


	int m;
	double  lambda;
	
public:
	SAEM() {
		m = 0; lambda = 1;
	}

	void add_m() {
		++m;
		// std::cout << m << std::endl;
	}

	double Stochastic_Approximation(column_vector &theta)
		/*
			Questa funzione calcola il valore dell'approssimazione stocastica della log likelihood
		*/
	{
		// la parte con l'huge momentaneamente la commento, poichè posso definire i vincoli nella funzione "find_max_gloabl"
		/*
		for (int i = 0; i < 6; ++i) {
			if (theta(i) < parametri.PARMIN[i] || theta(i) > parametri.PARMAX[i]) {
				return -1.e+150;
			}
		}
		*/
		if (m <= my_data::global_data.sc_M1) {
			lambda = 1;
		}
		else {
			lambda = 1 / pow((m - my_data::global_data.sc_M1), 0.8);
		}

		double logPY = parametri.BM_SDE_Compute_likelihood(theta, m);
		double Qm_1 = parametri.PreviousQm(theta);
		return Qm_1 + lambda * (logPY - Qm_1);

	}

} saem ;

