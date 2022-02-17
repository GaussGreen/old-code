/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bivariate_normal.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include <cmath>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

double bivariate_cdfNormal_phi(double a, double b, double rho)
{
	double x[5]={0.24840615,.39233107,.21141819,.033246660,.00082485334};
	double y[5]={0.10024215,.48281397,1.0609498,1.7797294,2.6697604};
	double a1=a/sqrt(2.*(1-rho*rho));
	double b1=b/sqrt(2.*(1-rho*rho));
	double sum=0.;
	int i,j;
	for(i=0;i<5;i++) for(j=0;i<5;i++)
	{
		sum+=x[i]*x[j]*exp(a1*(2.*y[i]-a1)+b1*(2.*y[j]-b1)+2.*rho*(y[i]-a1)*(y[j]-b1));
	}
	return sqrt(1-rho*rho)/ARM_NumericConstants::ARM_PI*sum;
}

double bivariate_cdfNormal_aux(double a, double b, double rho)
{
	if ((a<=0)&&(b<=0)&&(rho<=0)) return bivariate_cdfNormal_phi(a,b,rho);
	if ((a<=0)&&(b>=0)&&(rho>=0)) return ARM_GaussianAnalytics::cdfNormal(a)-bivariate_cdfNormal_phi(a,-b,-rho);
	if ((a>=0)&&(b<=0)&&(rho>=0)) return ARM_GaussianAnalytics::cdfNormal(b)-bivariate_cdfNormal_phi(-a,b,-rho);
	if ((a>=0)&&(b>=0)&&(rho<=0)) return ARM_GaussianAnalytics::cdfNormal(a)+ARM_GaussianAnalytics::cdfNormal(b)-1.+bivariate_cdfNormal_phi(-a,-b,-rho);
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bivariate_cdfNormal_aux: should not reach here");
	return 0;

}

double bivariate_cdfNormal(double a, double b, double rho)
{
	if(a*b*rho<=0) return bivariate_cdfNormal_aux(a,b,rho);
	else
	{
		double rho1,rho2, delta;
		if (a>=0) rho1=(rho*a-b)/sqrt(a*a-2.*rho*a*b+b*b); else rho1=-(rho*a-b)/sqrt(a*a-2.*rho*a*b+b*b);
		if (b>=0) rho2=(rho*b-a)/sqrt(a*a-2.*rho*a*b+b*b); else rho2=-(rho*b-a)/sqrt(a*a-2.*rho*a*b+b*b);
		if (((a>=0)&&(b>=0)) || ((a<0)&&(b<0)))  delta= 0; else delta=0.5;
		return bivariate_cdfNormal_aux(a,0,rho1)+bivariate_cdfNormal_aux(b,0,rho2)-delta;
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"bivariate_cdfNormal_aux: should not reach here");
	return 0;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/