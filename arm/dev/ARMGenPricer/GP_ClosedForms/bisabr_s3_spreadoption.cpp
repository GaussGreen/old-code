
/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_s3_spreadoption.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Nov 2006
 */

#include <glob/firsttoinc.h>

#include <cmath>
#include <complex>

#include <vector>
#include <iostream>
#include <iomanip>
#include <glob/expt.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/bisabr_s3_spreadoption.h"


using namespace std; 


CC_BEGIN_NAMESPACE(ARM)


double Packaged_BiSABR_S3_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
								  double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag)
{
		return 0;
}





CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/