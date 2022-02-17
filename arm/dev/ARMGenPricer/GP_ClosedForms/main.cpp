/*
 *
 *
 *
 * Version initiale 9/7/2003
 *
 */
/*----------------------------------------------------------------------------*
 
     Main.cpp

     GP_ClosedFormExtension main function
 
     
 
 
 
     Copyright (c) 2003 
 
*----------------------------------------------------------------------------*/

#include <glob/firsttoinc.h>

//#include <libCCTools++\CCString.h>

 #include "swapleg.h"


// #include <stdio.h>
// #include <windows.h>
#include <iostream>
// #include <cstring>
#include <vector>
#include <string>
#include <map>

// #include "gpinfra/dealdescription.h"
// #include "gpbase/ostringstream.h"
// #include "gpinfra/typedef.h"

#include <cmath>
#include <iomanip>
//
/*
#include <vector>
#include "basic_distributions.h"
#include "blackscholes.h"
#include "add_submodel.h"
#include "change_numeraire.h"
#include "input_extender.h"
#include "product.h"
#include "spreadoption_lognormal.h"
#include "spreadoption_lognormal_interface.h"
#include "gaussian_integrals.h"
#include "smile_sabr.h"
#include "spreadoption_normal.h"
#include "spreadoption_sabr.h"
#include "spreadoption_sabr_interface.h"
#include "spreadoption_shiftedlognormal.h"
#include "spreadoption_shiftedlognormal_interface.h"
#include "merton.h"
#include "extended_sabr.h"
#include "heston.h"
#include "integrals.h"
#include "hypergeometric.h"
#include "whittaker.h"
#include "bessel.h"
#include "cev.h"
#include "tri_spreadoption_lognormal.h"

#include "long_double.h"
#include "gamma.h"
#include "erf.h"
#include "inverse.h"
  */

#include "spreadoption_lognormal_interface.h"
#include "spreadoption_sabr_interface.h"
#include "tri_spreadoption_lognormal.h"
#include "Barriere_bs.h"
#include "Barriere_bs_formula.h"

using namespace std;
using namespace ARM;
// using namespace ARM_GP;

int main()		
{
	double S=95.;
	double X=90;
	double H=100;
	double K=3.;
	double T=1.0;
	double t1=0.5;
	double sig=0.25;
	double r=0.1;
	double b=0.1;
	int     callput=K_CALL;
	int     optiontype=ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_OUT;

	cout << "\n option="<<BS_PartialTime_End_SingleBarrierCall( S, X,  H, K, t1, T,  sig, r, b, optiontype);

	cout << "\n done" << endl;

	return 0;
}

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/