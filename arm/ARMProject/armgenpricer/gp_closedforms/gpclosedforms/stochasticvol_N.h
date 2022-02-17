/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file stochasticvol_n.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_STOCHASTICVOL_N_H
#define _GP_CF_STOCHASTICVOL_N_H


#include "firsttoinc.h"
#include "gpbase/port.h"
#include <vector>
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"


#include "expt.h"
using std::vector;

CC_BEGIN_NAMESPACE(ARM)


double StochasticVol_N_VanillaOption(double f,double K,double T,double drift,double sig,
								 double VolDrift,double VolVol,double callput,int nbsteps);

double Export_StochasticVol_N_VanillaOption(double f,double K,double T,double drift,double sig,
								 double VolDrift,double VolVol,double callput,int nbsteps);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


