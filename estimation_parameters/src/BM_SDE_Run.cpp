#include <ctime>
#include <global_variables.h>
#include <BM_SDE_ModelData.h>
#include <BM_SDE_ModelSetup.h>
#include <BM_SDE_SAEM_Optimized.h>



int main() {
	if (directory::path.size() < 3) {
		std::cout << "the directory where it's running the algorithm is \n"  << directory::ExePath() << std::endl << std::endl;
		// se si vuole cambiare directory basta andare nel file intestazione global_variables.h
		// e cambiare la directory da li...
	}
	else {
		std::cout << "the directory where it's running the algorithm is \n" << directory::path << std::endl << std::endl;
	}

	static clock_t start = clock();  // fa partire il timer 

	BM_SDE_ModelData::load(); // legge i dati per modificare guardare [BM_SDE_ModelData]	
	std::string fname = "../data/parametri_modello.txt";
	parametri.read_parameters_from_file(fname);
	parametri.initialization();    // inizializza i parametri (theta, tzeta, mu, Omega) per modificare guardare [BM_SDE_ModelSetup]

	BM_SDE_SAEM_Optimized::run();



	std::string file_da_scrivere = "../data/parametri_modello_stimati.txt";
	parametri.write_parameters_to_file(file_da_scrivere);

	std::string file_phi = "../data/last_phi.txt";
	parametri.save_last_phi(file_phi);

	std::string file_X = "../data/last_X.txt";
	parametri.save_last_X(file_X);
	//cout << my_data::global_data.YOBSi[5][3] << endl;


	static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
	static double time = (double) (end - start)/ CLOCKS_PER_SEC ;
	std::cout << "\n\n######  time execution total algorithm is: " << time << "  #######" << std::endl;



	system("pause");
	return 0;
}


