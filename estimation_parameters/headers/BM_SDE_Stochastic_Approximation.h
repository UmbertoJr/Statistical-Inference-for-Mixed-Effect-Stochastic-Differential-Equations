#pragma once

typedef dlib::matrix<double, 0, 1> column_vector;

class SAEM {


	static int m;
	double  lambda;
	double BM_SDE_Compute_likelihood(column_vector &theta) {
		double logpY = 0;
		for (int subj = 0; subj < my_data::global_data.sc_nSubjects; ++subj) {

		}
		return logpY;
	}

	double PreviousQm(column_vector & theta) {
		return 0;
	}
public:
	SAEM() {
		m = 0; lambda = 1;
	}

	void add_m() {
		++m;
		cout << m << endl;
	}

	double BM_SDE_Stochastic_Approximation(column_vector &theta)
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

		double logPY = BM_SDE_Compute_likelihood(theta);
		double Qm_1 = PreviousQm(theta);
		return Qm_1 + lambda * (logPY - Qm_1);

	}

} saem ;

