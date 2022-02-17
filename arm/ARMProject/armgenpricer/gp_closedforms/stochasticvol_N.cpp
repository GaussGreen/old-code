/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file merton.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>

#include "gpclosedforms/StochasticVol_N.h"
#include "gpclosedforms/vanilla_normal.h"


CC_BEGIN_NAMESPACE(ARM)

double StochasticVol_N_VanillaOption(double f,double K,double T,double r,double sig,double VolDrift,double VolOfVol,double callput,int LegendreNb)
{
	if(LegendreNb>1)
		{
			
			double mub=log(sig*sig*T)+(VolDrift-VolOfVol*VolOfVol/2.0)*T;
			double sigmab=ARM_NumericConstants::ARM_TWO_DIVIDED_BY_SQRT_3*VolOfVol*sqrt(T);
			ReducedGaussHermite_Coefficients c(LegendreNb);
			double Sum=0;
			double x,Vol;
			int i;
			for(i=0;i<LegendreNb;i++){
				x=exp(sigmab*c.get_point(i)+mub);
				Vol=sqrt(x/T);
				/// change the volatility according to the integrated volatility
				Sum+= VanillaOption_N(f,Vol,K,T,callput)*c.get_weight(i);
			}
			return Sum*exp(-r*T);
		}
		else
		{	
			return  exp(-r*T)*VanillaOption_N(f,sig,K,T,callput);
		}
}



double Export_StochasticVol_N_VanillaOption(double f,double K,double T,double drift,double sig,
											double VolDrift,double VolVol,double callput,int nbsteps)
{
	return StochasticVol_N_VanillaOption( f, K, T, drift, sig,
		VolDrift, VolVol, callput, nbsteps);
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
