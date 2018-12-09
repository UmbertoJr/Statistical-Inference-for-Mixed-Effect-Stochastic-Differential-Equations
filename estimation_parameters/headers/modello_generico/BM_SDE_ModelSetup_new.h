#pragma once

//#include "global_variables.h"

namespace BM_SDE_ModelSetup {
	class parameters {
		double Tzero = -90; // (bigtheta 001);
		double Tend = 360; // (bigtheta 002);
		double Tdelta = 1; // (bigtheta 003);
		double Vg = 0.19; // (bigtheta 004);
		double Ragut0 = 0; // (bigtheta 005);
		double Rkid0 = 0; // (bigtheta 006);
		double Insu0 = 100; // (bigtheta 007);
		double kxGI1 = 2.5e-05; // (bigtheta 008);
		double kxGI2 = 0.6e-05; // (bigtheta 009);
		double kxGI3 = 1.0e-06; // (bigtheta 010);
		double Weight = 120; // (bigtheta 011);
		double Subj = 4; // (bigtheta 012);
		double lambdaIG = 0.012; // (bigtheta 013);
		double gammaI = 0.58; // (bigtheta 014);
		double alpha = 0.46; // (bigtheta 015);
		double Rinf0 = 0; // (bigtheta 016);
		double kb = 0.0346; // (bigtheta 017)
		double kbH = 0.0179; // (bigtheta 018)
		double Gluc0 = 5; // (bigtheta 019);
		double kappa3max = 0.000487025; // (bigtheta 020);
		double Raliv0 = 0.000475; // (bigtheta 021);
		double Rupt0 = 0.0025; // (bigtheta 022);
		double RhinfConst = 0.00022; // (bigtheta 023);
		double DHot0 = 0.022; // (bigtheta 024);
		double Hot0 = 0.115789; // (bigtheta 025);
		double TTR0 = 0.0231579; // (bigtheta 026);
		//---- - the last elements of tzeta----//
		double betapar = 0.0178;
		double sigma = 0.009;
		//---------------------------------- - //

		double omega_Vg = 0.1*Vg;
		double omega_kxGI1 = 0.1*kxGI1;
		double omega_kxGI2 = 0.1*kxGI2;
		double omega_kxGI3 = 0.1*kxGI3;
		double omega_lambdaIG = 0.1*lambdaIG;
		double omega_gammaI = 0.1*gammaI;
		double omega_alpha = 0.1*alpha;
		double omega_kb = 0.1*kb;


	public:

		double tzeta[28][1];

		/*IsoglicemieParvals2Tzeta     //da fare questo file

		IsoglicemieInitializeParmaskPar % Determine which parameters in tzeta are active   //da fare questo file

		global_data.NFREEPAR = sum(global_data.PARMASK);    questo va messo nelle global_variables*/

		double *mu = new double[global_data.NFREEPAR];
		double *Omega = new double[global_data.NFREEPAR][NFREEPAR];
		double *theta = new double[ 2 * global_data.NFREEPAR + 2];

	};
}
