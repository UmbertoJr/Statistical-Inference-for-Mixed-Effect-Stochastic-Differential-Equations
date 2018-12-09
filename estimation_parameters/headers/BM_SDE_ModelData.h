#pragma once

#include "BM_getdata.h"



namespace BM_SDE_ModelData {
	void load() {
		// init dynamical memory
		my_data::global_data.TIME = new double[my_data::global_data.sc_tot_num_obs];
		my_data::global_data.YOBS = new double[my_data::global_data.sc_tot_num_obs];
		my_data::global_data.SUBJ = new double[my_data::global_data.sc_tot_num_obs];

		string Datafile = "time";
		BM_getdata(Datafile, my_data::global_data.TIME, my_data::global_data.sc_tot_num_obs);
		Datafile = "YOBS";
		BM_getdata(Datafile, my_data::global_data.YOBS, my_data::global_data.sc_tot_num_obs);
		Datafile = "SUBJ";
		BM_getdata(Datafile, my_data::global_data.SUBJ, my_data::global_data.sc_tot_num_obs);
		
		// init the matrices related to len subjs and timei and yobsi
		my_data::global_data.define_length_subjs();
		my_data::global_data.init_timei_and_yobsi();

		//delete dynamic memory
		delete[] my_data::global_data.TIME, my_data::global_data.YOBS, my_data::global_data.SUBJ;
	}
}