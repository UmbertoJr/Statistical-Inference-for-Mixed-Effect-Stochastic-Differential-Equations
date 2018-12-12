#pragma once
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <stdio.h>
#include <vector>

bool exists_test(const std::string& name);

int BM_getdata(std::string& name_file, double array[], int size) {
	/* HASCOVARS = BM_getdata(datafilename)
		
		loads data and covariates
		
		IN:     datafilename
		
		OUT : HASCOVARS : 0 no, 1 yes  there are covariates
		
		*/
	std::string path = "../data/";
	std::string extention = ".txt";
	int data_present = exists_test(path + name_file + extention);
	if (data_present == 0) {
		printf("\n\n BM_getdata: datafile  %s%s%s  NOT FOUND \n\n",path.c_str(), name_file.c_str(), extention.c_str());
		return 0;
	}
	else {
		// qui leggo il file e lo salvo in un vettore
		std::ifstream file;
		file.open(path + name_file + extention);
		file.seekg(0);
		double num = 0;
		for (int ite = 0; ite < size; ++ite) { // Iterates through the file till the pointer reads the last line of the file.
			file >> num;
			array[ite] = num;
			//cout << array[ite];

		}

		file.close();
	}
	
	return 1;
}


bool exists_test(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

