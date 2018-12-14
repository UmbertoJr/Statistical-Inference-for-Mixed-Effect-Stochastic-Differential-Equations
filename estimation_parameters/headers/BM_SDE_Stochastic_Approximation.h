#pragma once
#include <dlib/optimization.h>
#include <dlib/global_optimization.h>
#if PR_DEBUG==1
#define log_space(x) printer(#x, (x));std::cout << "\t"
#define log_line(x) printer(#x, (x));std::cout <<  "\n"
#endif

void printer(const char *name, double value) {
	printf("%s\t =  %f", name, value);
}

typedef dlib::matrix<double, 0, 1> column_vector;

class SAEM {


public:

	int m;
	double  lambda;
	//double lambda;

	SAEM() {
		m = 0; lambda = 1;
		
	}

	void add_m() {
		++m;
		// std::cout << m << std::endl;
	}

	void run(int number_iterations) 
	{	
		column_vector  parmin(6), parmax(6); // theta(5),
		for (int i = 0; i < 6; ++i) {
			//theta(i) = parametri.theta[i];
			parmin(i) = parametri.PARMIN[i]; parmax(i) = parametri.PARMAX[i];
			//std::cout << theta(i) << std::endl;
		}

		auto Stochastic_Approximation = [&](const column_vector &theta)
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
			double * thetai = new double[6];
			for (int i = 0; i < 6; ++i) {
				thetai[i] = theta(i);
			}
			if (saem.m <= my_data::global_data.sc_M1) {
				lambda = 1;
			}
			else {
				lambda = 1 / pow((m - my_data::global_data.sc_M1), 0.8);
			}

			double logPY = parametri.BM_SDE_Compute_likelihood(thetai, m);
			double Qm_1 = parametri.PreviousQm(thetai, m, lambda);
			//log_space(Qm_1); //log_line(seam.m);
			double Qm = Qm_1 + lambda * (logPY - Qm_1);
			return Qm;

		};

		auto result = dlib::find_max_global(Stochastic_Approximation, parmin, parmax, dlib::max_function_calls(number_iterations));
		std::cout.precision(9);
		// These cout statements will show that find_min_global() found the
		// globally optimal solution to 9 digits of precision:
		std::cout << "function solution y : " << result.y << std::endl;
		std::cout << "function solution x:\n" << result.x << std::endl;

		for (int i = 0; i < 6; ++i) {
			parametri.theta[i] = result.x(i);
		}
	}

} saem ;




