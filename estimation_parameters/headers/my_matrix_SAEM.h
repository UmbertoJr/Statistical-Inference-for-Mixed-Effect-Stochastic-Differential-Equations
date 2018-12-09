#pragma once

static class Phi {

	double data[my_data::global_data.sc_M][my_data::global_data.sc_nSubjects][parametri.len_mu];
	double*ptr = &data[0][0][0];

public:
	void copy_phi_to_PHIALL(double* arr, int a, int b, int c) {
		int start = a * my_data::global_data.sc_nSubjects*parametri.len_mu + b * parametri.len_mu + c;
		for (int i = 0; i < parametri.len_mu; i++) {
			*(ptr + start + i) = arr[i];
		}

	}

	double value_in(int a, int b, int c) {
		return data[a][b][c];
	}

} PHI;



static class XALL {
	double data[my_data::global_data.sc_M][my_data::global_data.sc_nSubjects][my_data::global_data.sc_NMAXOBS];
	double*ptr = &data[0][0][0];
public:
	void copy_newX_to_XALL(double* arr, int a, int b, int c) {
		int start = a * my_data::global_data.sc_nSubjects*my_data::global_data.sc_NMAXOBS + b * my_data::global_data.sc_NMAXOBS + c;
		for (int i = 0; i < my_data::global_data.sc_NMAXOBS; i++) {
			*(ptr + start + i) = arr[i];
		}

	}
}XALL;