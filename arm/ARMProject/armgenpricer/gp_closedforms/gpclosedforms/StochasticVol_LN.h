/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file stochasticvol_ln.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_STOCHASTICVOL_LN_H
#define _GP_CF_STOCHASTICVOL_LN_H


#include "firsttoinc.h"
#include "gpbase/port.h"
#include <vector>
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"


#include "expt.h"
using std::vector;

CC_BEGIN_NAMESPACE(ARM)



double StochasticVol_LN_Arithmetic_VanillaOption(double f,double K,double T,double r,double sig,double VolDrift,
												 double VolOfVol,double averaging,double callput,int LegendreNb);



double StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(double f,double K,double T,double r,double sig,double VolDrift,
												 double VolOfVol,double averaging,double reset,double callput,int LegendreNb);

double StochasticVol_LN_Geometric_VanillaOption(double f,double K,double T,double r,double sig,double VolDrift,
												double VolOfVol,double averaging,double callput,int LegendreNb);

double StochasticVol_LN_Geometric_VanillaOption_with_Reset(double f,double K,double T,double r,double sig,double VolDrift,
												 double VolOfVol,double averaging,double reset,double callput,int LegendreNb);








CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


