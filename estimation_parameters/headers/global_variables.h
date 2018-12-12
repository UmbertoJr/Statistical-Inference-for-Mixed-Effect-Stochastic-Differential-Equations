#pragma once
#include <windows.h>
#include <string>
#include <iostream>



void check_memory_alloc(double* point);

namespace directory {
	// variabile che definisce la directory dove si trovano le cartele di input e output 
	static std::string path = "./";

	// questa funziona cerca la directory attuale dove gira l'algoritmo 
	static std::string ExePath() {
		char buffer[MAX_PATH];
		GetModuleFileName(NULL, buffer, MAX_PATH);
		std::string::size_type pos = std::string(buffer).find_last_of("\\/");
		return std::string(buffer).substr(0, pos);
	}

}


namespace my_data {
	static struct global {
		int VRBL;
		double PARTICLEXS;
		double PARTICLEWS;
		int PARMASK;
		static const int sc_nSubjects = 6;		// numero totale di soggetti
		static const int sc_tot_num_obs = 86;
		double* TIME;
		double* YOBS;
		double* SUBJ;
		int len_subj[sc_nSubjects];
		static const int sc_NMAXOBS = 18;		// numero osservazioni massimo per soggetto
		double TIMEi[sc_nSubjects][sc_NMAXOBS];
		double YOBSi[sc_nSubjects][sc_NMAXOBS];
		std::string MYOPT;
		int PhiAll;
		int XAll;
		static const int sc_M = 10;  // numero di iterazioni da svolgere per SAEM... da leggere da input
		static const int sc_M1 = 6;  // numero di iterazioni dopo il quale la constante lambda inizia a decrescere  ##### qui c'è da cambiare
		int NSTATEVARS;
		int NFREEPAR;

		void define_length_subjs() {
			// questa funzione trova la lunghezza del vettore osservazione per ciascun soggetto 
			//  e ritorna un vettore di lunghezza sc_nSubjects.
			for (int subj_i = 1; subj_i <= sc_nSubjects; subj_i++) {
				int count = 0;
				for (int j = 0; j < sc_tot_num_obs; j++) {
					if (SUBJ[j] == subj_i) {
						count++;
					}
				}
				len_subj[subj_i - 1] = count;
			}
		}

		void init_timei_and_yobsi() {
			
			for (int subj_n = 1; subj_n <= sc_nSubjects; subj_n++) {
				int j = 0;
				//cout << subj_n << "\t" << len_subj[subj_n-1] << "\t" << endl;
				for (int i = 0; i < sc_tot_num_obs; i++) {
					if (SUBJ[i] == subj_n) {
						TIMEi[subj_n -1][j] = TIME[i];
						YOBSi[subj_n - 1][j] = YOBS[i];
						j++;
						//cout << j << "\t" << YOBSi[subj_n - 1][j-1] << "\t";   // printa i risultati per debug
					}
				}
				//cout << endl << endl;
			}
			
		}
	} global_data;
}



void check_memory_alloc(double* point ) {
	if (!(point)) {
		std::cout << "Error: out of memory." << std::endl;
		exit(1);
	}
}
