#include <ctime>
#include <global_variables.h>
#include <BM_SDE_ModelData.h>
#include <BM_SDE_ModelSetup.h>
#include <BM_SDE_SAEM_Optimized.h>


using namespace std;

int main() {
	if (directory::path.size() < 3) {
		cout << "the directory where it's running the algorithm is \n"  << directory::ExePath() << endl << endl;
		// se si vuole cambiare directory basta andare nel file intestazione global_variables.h
		// e cambiare la directory da li...
	}
	else {
		cout << "the directory where it's running the algorithm is \n" << directory::path << endl << endl;
	}

	static clock_t start = clock();  // fa partire il timer 

	BM_SDE_ModelData::load(); // legge i dati per modificare guardare [BM_SDE_ModelData]	
	string fname = "../data/parametri_modello.txt";
	parametri.read_parameters_from_file(fname);
	parametri.initialization();    // inizializza i parametri (theta, tzeta, mu, Omega) per modificare guardare [BM_SDE_ModelSetup]

	BM_SDE_SAEM_Optimized::run();



	string file_da_scrivere = "../data/parametri_modello.txt";
	parametri.write_parameters_to_file(file_da_scrivere);





	static clock_t end = clock();  // stop per il timer e print del tempo algoritmo
	static double time = (double) (end - start)/ CLOCKS_PER_SEC ;
	cout << "\n\n######  time execution is: " << time << "  #######" << endl;

	system("pause");
	return 0;
}


