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

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/merton.h"
#include "gpclosedforms/gamma.h"
#include "gpclosedforms/vanilla_bs.h"


CC_BEGIN_NAMESPACE(ARM)



double merton(double F,double K,double t,double sigma, double lambda, double muJ, double sigmaJ,int nb)
{
	double lambdap=lambda*exp(muJ+sigmaJ*sigmaJ/2);
	int i;
	double sum=0;
	long_double r;
	if(lambda<0.000000001)
	{
		return BS(F,K,t,sigma);
	}
	else
	{
		for(i=0;i<nb;i++)
		{
			r=long_double(i*log(lambda*t),1)/gammaLD((long double) (i+1));
			sum+= r.todouble()*
				BS(F*exp((lambda-lambdap)*t+i*(muJ+sigmaJ*sigmaJ/2)),K,t,sqrt(sigma*sigma+i*sigmaJ*sigmaJ/t))*
				exp(-lambda*t);
		}
		return sum;
	}
}

double MertonOption(double S, double K, int CallPut, double t, double sigma, double lambda1, double U1, double lambda2, double U2, int nb)
{
	if(CallPut == -1)
	{
		return MertonOption(S, K, 1, t, sigma, lambda1, U1, lambda2, U2, nb) + K - S;
	}

	int i, j;

	double sum = 0.;
	double fac = exp(-(lambda1 + lambda2) * t), facij = exp(-U1*lambda1*t + U2*lambda2*t);
	double Lij;
	double coeffi = 1., coeffj, coeffbsi = 1., coeffbsj;
	double faci = 1., facj;

	for(i = 0; i < nb; i++)
	{
		coeffj = 1.;
		coeffbsj = 1.;
		facj = 1.;
		for(j = 0; j < nb; j++)
		{
			if(faci*facj > 1e12) break;
			Lij = S*coeffi*coeffj*facij;
			sum += BS(Lij,K,t,sigma,1) * coeffbsj * coeffbsi / faci / facj;
			coeffj *= 1-U2;
			coeffbsj *= lambda2*t;
			facj *= j+1;
		}
		coeffi *= 1.+U1;
		coeffbsi *= lambda1*t;
		faci *= i+1;
	}

	return sum*fac;
}


CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
