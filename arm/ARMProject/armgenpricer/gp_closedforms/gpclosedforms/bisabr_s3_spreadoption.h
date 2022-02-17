/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_s3_spreadoption.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Nov 2006
 */
 
#ifndef _GP_CF_BISABR_S3_SPREADOPTION_H
#define _GP_CF_BISABR_S3_SPREADOPTION_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include <complex>
using namespace std; 

CC_BEGIN_NAMESPACE(ARM)

double Packaged_BiSABR_SpreadOption(double F1,double alpha1,double beta1,double rho1,double nu1,double F2,double alpha2,double beta2,double rho2,double nu2,
								  double K,double T,int CallPut,double rhos,double rhov,double rhoc12,double rhoc21,int flag);



CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


