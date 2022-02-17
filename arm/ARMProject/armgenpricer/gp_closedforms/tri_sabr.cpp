/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file bisabr_spreadoption.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2006
 */

#include "firsttoinc.h"

#include <cmath>
#include <complex>

#include <vector>
#include <iostream>
#include <iomanip>
#include "expt.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/erf.h"
#include "gpclosedforms/hypergeometric.h"
#include "gpclosedforms/tri_sabr.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/extended_sabr_interface.h"

using namespace std; 


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_K_SHIFT_FOR_DERIVATION 0.0000001

/*double power(double a,double b)
{	
	if(b==-1)	return 1./a;
	if(b==2)	return a*a;
	if(b==3)	return a*a*a;
	if(b==4)	return a*a*a*a;
	return pow(a,b);
}*/



#define power(a,b) (((b)==-1)?1./(a):(((b)==2)?(a)*(a):(((b)==3)?(a)*(a)*(a):(((b)==4)?(a)*(a)*(a)*(a):pow((a),(b))))))

void TriSABRSumOptionObservables(
			double F1,double alpha1, double beta1, double rho1,double nu1,
			double F2,double alpha2, double beta2, double rho2,double nu2,
			double F3,double alpha3, double beta3, double rho3,double nu3,
			double rhos12, double rhos23, double rhos13,
			double rhov12, double rhov23, double rhov13,
			double rhoc12, double rhoc13,
			double rhoc21, double rhoc23,
			double rhoc31, double rhoc32,
			double beta,
			double* Stddev,double* varSABR,double* covarSABR)
{
	double SABRdrift,F1b1=pow(F1,beta1),F2b2=pow(F2,beta2),F3b3=pow(F3,beta3),Sb,S=F1+F2-F3;
	Sb=pow(F1+F2-F3,beta);
	
	*Stddev=sqrt(2*alpha1*F1b1*(alpha2*F2b2*rhos12 - alpha3*F3b3*rhos13) - 
		2*alpha2*alpha3*F2b2*F3b3*rhos23 + power(alpha1,2)*power(F1b1,2) + 
		power(alpha2,2)*power(F2b2,2) + power(alpha3,2)*power(F3b3,2));
	
	*varSABR=power(Sb,-4)*(-2*beta*Sb*power(alpha1,5)*power(F1b1,5)*power(S,-1)*
		(Sb*(nu1*rho1 + alpha1*beta1*F1b1*power(F1,-1)) + 3*beta*(-(alpha2*F2b2*rhos12) + alpha3*F3b3*rhos13)*Sb*power(S,-1)) - 
		2*beta*Sb*power(alpha2,5)*power(F2b2,5)*power(S,-1)*
		(Sb*(nu2*rho2 + alpha2*beta2*F2b2*power(F2,-1)) + 3*alpha3*beta*F3b3*rhos23*Sb*power(S,-1)) + 
		power(alpha1,6)*power(beta,2)*power(F1b1,6)*power(S,-2)*power(Sb,2) + 
		power(alpha2,6)*power(beta,2)*power(F2b2,6)*power(S,-2)*power(Sb,2) + 
		power(alpha3,4)*power(F3b3,4)*((2*alpha3*beta3*F3b3*nu3*rho3*power(F3,-1) + 
		power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2) + power(nu3,2))*power(Sb,2) + 
        power(alpha3,2)*power(beta,2)*power(F3b3,2)*power(S,-2)*power(Sb,2) + 
        2*alpha3*beta*F3b3*(nu3*rho3 + alpha3*beta3*F3b3*power(F3,-1))*power(S,-1)*power(Sb,2)) + 
		power(alpha1,4)*power(F1b1,4)*((2*alpha1*beta1*F1b1*nu1*rho1*power(F1,-1) + 
		power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2) + power(nu1,2))*power(Sb,2) + 
        3*power(beta,2)*(-2*alpha2*alpha3*F2b2*F3b3*(4*rhos12*rhos13 + rhos23) + power(alpha2,2)*power(F2b2,2)*(1 + 4*power(rhos12,2)) + 
		power(alpha3,2)*power(F3b3,2)*(1 + 4*power(rhos13,2)))*power(S,-2)*power(Sb,2) - 
        2*beta*(alpha2*F2b2*(nu1*rhoc21 + 3*nu1*rho1*rhos12 + nu2*rhoc12*rhos12 + 4*alpha1*beta1*F1b1*rhos12*power(F1,-1) + 
		alpha2*beta2*F2b2*power(F2,-1)*power(rhos12,2)) - 
		alpha3*F3b3*(nu1*rhoc31 + 3*nu1*rho1*rhos13 + nu3*rhoc13*rhos13 + 4*alpha1*beta1*F1b1*rhos13*power(F1,-1) + 
		alpha3*beta3*F3b3*power(F3,-1)*power(rhos13,2)))*power(S,-1)*power(Sb,2)) - 
		2*alpha2*F2b2*power(alpha3,3)*power(F3b3,3)*(rhos23*
		(nu3*(nu3 + nu2*rhov23) + alpha3*beta3*F3b3*(2*nu3*rho3 + nu2*rhoc32)*power(F3,-1) + 
		alpha2*beta2*F2b2*power(F2,-1)*(nu3*rhoc23 + alpha3*beta3*F3b3*rhos23*power(F3,-1)) + 
		power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2))*power(Sb,2) + 
        3*rhos23*power(alpha3,2)*power(beta,2)*power(F3b3,2)*power(S,-2)*power(Sb,2) + 
        alpha3*beta*F3b3*(nu3*rhoc23 + 3*nu3*rho3*rhos23 + nu2*rhoc32*rhos23 + 4*alpha3*beta3*F3b3*rhos23*power(F3,-1) + 
		alpha2*beta2*F2b2*power(F2,-1)*power(rhos23,2))*power(S,-1)*power(Sb,2)) + 
		power(alpha2,4)*power(F2b2,4)*((2*alpha2*beta2*F2b2*nu2*rho2*power(F2,-1) + 
		power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2) + power(nu2,2))*power(Sb,2) + 
        3*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(1 + 4*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
        2*alpha3*beta*F3b3*(nu2*rhoc32 + 3*nu2*rho2*rhos23 + nu3*rhoc23*rhos23 + 4*alpha2*beta2*F2b2*rhos23*power(F2,-1) + 
		alpha3*beta3*F3b3*power(F3,-1)*power(rhos23,2))*power(S,-1)*power(Sb,2)) - 
		2*alpha3*F3b3*power(alpha2,3)*power(F2b2,3)*(rhos23*
		(power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2) + nu2*(nu2 + nu3*rhov23 + alpha3*beta3*F3b3*rhoc32*power(F3,-1)) + 
		alpha2*beta2*F2b2*power(F2,-1)*(2*nu2*rho2 + nu3*rhoc23 + alpha3*beta3*F3b3*rhos23*power(F3,-1)))*power(Sb,2) + 
        2*rhos23*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(3 + 2*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
        alpha3*beta*F3b3*(nu2*rho2 + nu3*rhoc23 + nu3*rho3*rhos23 + 3*nu2*rhoc32*rhos23 + 2*nu2*rho2*power(rhos23,2) + 
		2*nu3*rhoc23*power(rhos23,2) + 2*alpha3*beta3*F3b3*rhos23*power(F3,-1)*(1 + power(rhos23,2)) + 
		beta2*F2b2*power(F2,-1)*(alpha2 + 5*alpha2*power(rhos23,2)))*power(S,-1)*power(Sb,2)) + 
		power(alpha2,2)*power(alpha3,2)*power(F2b2,2)*power(F3b3,2)*
		((2*nu2*nu3*rhov23 + 2*nu2*nu3*rhov23*power(rhos23,2) + power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2)*power(rhos23,2) + 
		power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2)*power(rhos23,2) + power(nu2,2)*power(rhos23,2) + 
		power(nu3,2)*power(rhos23,2) + 2*alpha3*beta3*F3b3*power(F3,-1)*
		(nu3*rho3*power(rhos23,2) + nu2*rhoc32*(1 + power(rhos23,2))) + 
		2*alpha2*beta2*F2b2*power(F2,-1)*(nu2*rho2*power(rhos23,2) + nu3*rhoc23*(1 + power(rhos23,2)) + 
		alpha3*beta3*F3b3*rhos23*power(F3,-1)*(1 + power(rhos23,2))))*power(Sb,2) + 
        3*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(1 + 4*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
        2*alpha3*beta*F3b3*(nu3*rho3 + nu2*rhoc32 + nu2*rho2*rhos23 + 3*nu3*rhoc23*rhos23 + 2*nu3*rho3*power(rhos23,2) + 
		2*nu2*rhoc32*power(rhos23,2) + 2*alpha2*beta2*F2b2*rhos23*power(F2,-1)*(1 + power(rhos23,2)) + 
		beta3*F3b3*power(F3,-1)*(alpha3 + 5*alpha3*power(rhos23,2)))*power(S,-1)*power(Sb,2)) + 
		2*power(alpha1,3)*power(F1b1,3)*(-(beta*Sb*power(alpha2,2)*power(F2b2,2)*power(S,-1)*
		(Sb*(nu1*rho1 + nu2*rhoc12 + nu2*rho2*rhos12 + 3*nu1*rhoc21*rhos12 + 2*nu1*rho1*power(rhos12,2) + 
		2*nu2*rhoc12*power(rhos12,2) + 2*alpha2*beta2*F2b2*rhos12*power(F2,-1)*(1 + power(rhos12,2)) + 
		beta1*F1b1*power(F1,-1)*(alpha1 + 5*alpha1*power(rhos12,2))) + 
		6*alpha3*beta*F3b3*Sb*(rhos13 + 2*rhos12*rhos23 + 2*rhos13*power(rhos12,2))*power(S,-1))) + 
        2*rhos12*power(alpha2,3)*power(beta,2)*power(F2b2,3)*(3 + 2*power(rhos12,2))*power(S,-2)*power(Sb,2) + 
        alpha2*F2b2*(rhos12*(power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2) + 
		nu1*(nu1 + nu2*rhov12 + alpha2*beta2*F2b2*rhoc21*power(F2,-1)) + 
		alpha1*beta1*F1b1*power(F1,-1)*(2*nu1*rho1 + nu2*rhoc12 + alpha2*beta2*F2b2*rhos12*power(F2,-1)))*power(Sb,2) + 
		6*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(rhos12 + 2*rhos13*rhos23 + 2*rhos12*power(rhos13,2))*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(3*nu1*rhoc31*rhos12 + nu2*rhoc32*rhos12 + 3*nu1*rhoc21*rhos13 + nu3*rhoc23*rhos13 + 
		4*nu1*rho1*rhos12*rhos13 + 2*nu2*rhoc12*rhos12*rhos13 + 2*nu3*rhoc13*rhos12*rhos13 + 2*nu1*rho1*rhos23 + 
		nu2*rhoc12*rhos23 + nu3*rhoc13*rhos23 + 2*alpha1*beta1*F1b1*(5*rhos12*rhos13 + rhos23)*power(F1,-1) + 
		2*alpha2*beta2*F2b2*rhos12*(rhos12*rhos13 + rhos23)*power(F2,-1) + 2*alpha3*beta3*F3b3*rhos13*rhos23*power(F3,-1) + 
		2*alpha3*beta3*F3b3*rhos12*power(F3,-1)*power(rhos13,2))*power(S,-1)*power(Sb,2)) - 
        alpha3*F3b3*(rhos13*(power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2) + 
		nu1*(nu1 + nu3*rhov13 + alpha3*beta3*F3b3*rhoc31*power(F3,-1)) + 
		alpha1*beta1*F1b1*power(F1,-1)*(2*nu1*rho1 + nu3*rhoc13 + alpha3*beta3*F3b3*rhos13*power(F3,-1)))*power(Sb,2) + 
		2*rhos13*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(3 + 2*power(rhos13,2))*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(nu1*rho1 + nu3*rhoc13 + nu3*rho3*rhos13 + 3*nu1*rhoc31*rhos13 + 2*nu1*rho1*power(rhos13,2) + 
		2*nu3*rhoc13*power(rhos13,2) + 2*alpha3*beta3*F3b3*rhos13*power(F3,-1)*(1 + power(rhos13,2)) + 
		beta1*F1b1*power(F1,-1)*(alpha1 + 5*alpha1*power(rhos13,2)))*power(S,-1)*power(Sb,2))) + 
		power(alpha1,2)*power(F1b1,2)*(-2*beta*Sb*power(alpha2,3)*power(F2b2,3)*power(S,-1)*
		(Sb*(nu2*rho2 + nu1*rhoc21 + nu1*rho1*rhos12 + 3*nu2*rhoc12*rhos12 + 2*nu2*rho2*power(rhos12,2) + 2*nu1*rhoc21*power(rhos12,2) + 
		2*alpha1*beta1*F1b1*rhos12*power(F1,-1)*(1 + power(rhos12,2)) + beta2*F2b2*power(F2,-1)*(alpha2 + 5*alpha2*power(rhos12,2)))
		+ 6*alpha3*beta*F3b3*Sb*(2*rhos12*rhos13 + rhos23 + 2*rhos23*power(rhos12,2))*power(S,-1)) + 
        3*power(alpha2,4)*power(beta,2)*power(F2b2,4)*(1 + 4*power(rhos12,2))*power(S,-2)*power(Sb,2) + 
        power(alpha3,2)*power(F3b3,2)*((2*nu1*nu3*rhov13 + 2*nu1*nu3*rhov13*power(rhos13,2) + 
		power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2)*power(rhos13,2) + 
		power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2)*power(rhos13,2) + power(nu1,2)*power(rhos13,2) + 
		power(nu3,2)*power(rhos13,2) + 2*alpha3*beta3*F3b3*power(F3,-1)*
		(nu3*rho3*power(rhos13,2) + nu1*rhoc31*(1 + power(rhos13,2))) + 
		2*alpha1*beta1*F1b1*power(F1,-1)*(nu1*rho1*power(rhos13,2) + nu3*rhoc13*(1 + power(rhos13,2)) + 
		alpha3*beta3*F3b3*rhos13*power(F3,-1)*(1 + power(rhos13,2))))*power(Sb,2) + 
		3*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(1 + 4*power(rhos13,2))*power(S,-2)*power(Sb,2) + 
		2*alpha3*beta*F3b3*(nu3*rho3 + nu1*rhoc31 + nu1*rho1*rhos13 + 3*nu3*rhoc13*rhos13 + 2*nu3*rho3*power(rhos13,2) + 
		2*nu1*rhoc31*power(rhos13,2) + 2*alpha1*beta1*F1b1*rhos13*power(F1,-1)*(1 + power(rhos13,2)) + 
		beta3*F3b3*power(F3,-1)*(alpha3 + 5*alpha3*power(rhos13,2)))*power(S,-1)*power(Sb,2)) + 
        power(alpha2,2)*power(F2b2,2)*((2*nu1*nu2*rhov12 + 2*nu1*nu2*rhov12*power(rhos12,2) + 
		power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2)*power(rhos12,2) + 
		power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2)*power(rhos12,2) + power(nu1,2)*power(rhos12,2) + 
		power(nu2,2)*power(rhos12,2) + 2*alpha2*beta2*F2b2*power(F2,-1)*
		(nu2*rho2*power(rhos12,2) + nu1*rhoc21*(1 + power(rhos12,2))) + 
		2*alpha1*beta1*F1b1*power(F1,-1)*(nu1*rho1*power(rhos12,2) + nu2*rhoc12*(1 + power(rhos12,2)) + 
		alpha2*beta2*F2b2*rhos12*power(F2,-1)*(1 + power(rhos12,2))))*power(Sb,2) + 
		6*power(alpha3,2)*power(beta,2)*power(F3b3,2)*
		(1 + 8*rhos12*rhos13*rhos23 + 2*power(rhos12,2) + 2*power(rhos13,2) + 2*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
		2*alpha3*beta*F3b3*(nu1*rhoc31 + nu2*rhoc32 + nu1*rho1*rhos13 + 2*nu2*rhoc12*rhos13 + nu3*rhoc13*rhos13 + 
		2*nu2*rho2*rhos12*rhos13 + 4*nu1*rhoc21*rhos12*rhos13 + 2*nu3*rhoc23*rhos12*rhos13 + nu2*rho2*rhos23 + 
		2*nu1*rhoc21*rhos23 + nu3*rhoc23*rhos23 + 2*nu1*rho1*rhos12*rhos23 + 4*nu2*rhoc12*rhos12*rhos23 + 
		2*nu3*rhoc13*rhos12*rhos23 + 4*alpha3*beta3*F3b3*rhos12*rhos13*rhos23*power(F3,-1) + 2*nu1*rhoc31*power(rhos12,2) + 
		2*nu2*rhoc32*power(rhos12,2) + 2*alpha1*beta1*F1b1*power(F1,-1)*(rhos13 + 2*rhos12*rhos23 + 3*rhos13*power(rhos12,2)) + 
		2*alpha2*beta2*F2b2*power(F2,-1)*(2*rhos12*rhos13 + rhos23 + 3*rhos23*power(rhos12,2)) + 
		alpha3*beta3*F3b3*power(F3,-1)*power(rhos13,2) + alpha3*beta3*F3b3*power(F3,-1)*power(rhos23,2))*power(S,-1)*power(Sb,2)) - 
        2*alpha2*alpha3*F2b2*F3b3*((nu1*nu2*rhos12*rhos13*rhov12 + nu1*nu2*rhos23*rhov12 + nu1*nu3*rhos12*rhos13*rhov13 + 
		nu1*nu3*rhos23*rhov13 + nu2*nu3*rhos12*rhos13*rhov23 + 
		rhos12*rhos13*power(alpha1,2)*power(beta1,2)*power(F1,-2)*power(F1b1,2) + 
		alpha3*beta3*F3b3*nu1*rhoc31*rhos12*rhos13*power(F3,-1) + alpha3*beta3*F3b3*nu2*rhoc32*rhos12*rhos13*power(F3,-1) + 
		alpha3*beta3*F3b3*nu1*rhoc31*rhos23*power(F3,-1) + 
		alpha2*beta2*F2b2*power(F2,-1)*(nu3*rhoc23*rhos12*rhos13 + nu1*rhoc21*(rhos12*rhos13 + rhos23) + 
		alpha3*beta3*F3b3*rhos12*rhos13*rhos23*power(F3,-1)) + 
		alpha1*beta1*F1b1*power(F1,-1)*(2*nu1*rho1*rhos12*rhos13 + nu2*rhoc12*rhos12*rhos13 + nu3*rhoc13*rhos12*rhos13 + 
		nu2*rhoc12*rhos23 + nu3*rhoc13*rhos23 + alpha2*beta2*F2b2*rhos12*(rhos12*rhos13 + rhos23)*power(F2,-1) + 
		alpha3*beta3*F3b3*rhos13*(rhos12*rhos13 + rhos23)*power(F3,-1)) + rhos12*rhos13*power(nu1,2))*power(Sb,2) + 
		6*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(2*rhos12*rhos13 + rhos23 + 2*rhos23*power(rhos13,2))*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(nu1*rhoc21 + nu3*rhoc23 + nu1*rho1*rhos12 + nu2*rhoc12*rhos12 + 2*nu3*rhoc13*rhos12 + 
		2*nu3*rho3*rhos12*rhos13 + 4*nu1*rhoc31*rhos12*rhos13 + 2*nu2*rhoc32*rhos12*rhos13 + nu3*rho3*rhos23 + 
		2*nu1*rhoc31*rhos23 + nu2*rhoc32*rhos23 + 2*nu1*rho1*rhos13*rhos23 + 2*nu2*rhoc12*rhos13*rhos23 + 
		4*nu3*rhoc13*rhos13*rhos23 + 4*alpha3*beta3*F3b3*rhos12*rhos13*power(F3,-1) + 2*alpha3*beta3*F3b3*rhos23*power(F3,-1) + 
		2*nu1*rhoc21*power(rhos13,2) + 2*nu3*rhoc23*power(rhos13,2) + 6*alpha3*beta3*F3b3*rhos23*power(F3,-1)*power(rhos13,2) + 
		2*alpha1*beta1*F1b1*power(F1,-1)*(rhos12 + 2*rhos13*rhos23 + 3*rhos12*power(rhos13,2)) + 
		alpha2*beta2*F2b2*power(F2,-1)*(4*rhos12*rhos13*rhos23 + power(rhos12,2) + power(rhos23,2)))*power(S,-1)*power(Sb,2))) - 
		2*alpha1*F1b1*(beta*Sb*power(alpha2,4)*power(F2b2,4)*power(S,-1)*
		(Sb*(nu2*rhoc12 + 3*nu2*rho2*rhos12 + nu1*rhoc21*rhos12 + 4*alpha2*beta2*F2b2*rhos12*power(F2,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*power(rhos12,2)) + 3*alpha3*beta*F3b3*(rhos13 + 4*rhos12*rhos23)*Sb*power(S,-1)) - 
        3*rhos12*power(alpha2,5)*power(beta,2)*power(F2b2,5)*power(S,-2)*power(Sb,2) - 
        alpha2*F2b2*power(alpha3,2)*power(F3b3,2)*((nu1*nu2*rhos13*rhos23*rhov12 + nu1*nu3*rhos12*rhov13 + nu1*nu3*rhos13*rhos23*rhov13 + 
		nu2*nu3*rhos12*rhov23 + nu2*nu3*rhos13*rhos23*rhov23 + alpha3*beta3*F3b3*nu1*rhoc31*rhos12*power(F3,-1) + 
		alpha3*beta3*F3b3*nu2*rhoc32*rhos12*power(F3,-1) + 2*alpha3*beta3*F3b3*nu3*rho3*rhos13*rhos23*power(F3,-1) + 
		alpha3*beta3*F3b3*nu1*rhoc31*rhos13*rhos23*power(F3,-1) + alpha3*beta3*F3b3*nu2*rhoc32*rhos13*rhos23*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*(nu3*rhoc13*rhos12 + nu2*rhoc12*rhos13*rhos23 + nu3*rhoc13*rhos13*rhos23 + 
		alpha2*beta2*F2b2*rhos12*rhos13*rhos23*power(F2,-1) + alpha3*beta3*F3b3*rhos13*(rhos12 + rhos13*rhos23)*power(F3,-1)) + 
		alpha2*beta2*F2b2*power(F2,-1)*(nu1*rhoc21*rhos13*rhos23 + nu3*rhoc23*(rhos12 + rhos13*rhos23) + 
		alpha3*beta3*F3b3*rhos23*(rhos12 + rhos13*rhos23)*power(F3,-1)) + 
		rhos13*rhos23*power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2) + rhos13*rhos23*power(nu3,2))*power(Sb,2) + 
		3*(rhos12 + 4*rhos13*rhos23)*power(alpha3,2)*power(beta,2)*power(F3b3,2)*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(2*nu3*rho3*rhos12 + nu1*rhoc31*rhos12 + nu2*rhoc32*rhos12 + nu1*rhoc21*rhos13 + 3*nu3*rhoc23*rhos13 + 
		nu2*rhoc12*rhos23 + 3*nu3*rhoc13*rhos23 + 4*nu3*rho3*rhos13*rhos23 + 2*nu1*rhoc31*rhos13*rhos23 + 
		2*nu2*rhoc32*rhos13*rhos23 + 2*alpha1*beta1*F1b1*rhos13*(rhos12 + rhos13*rhos23)*power(F1,-1) + 
		2*alpha2*beta2*F2b2*rhos23*(rhos12 + rhos13*rhos23)*power(F2,-1) + 2*alpha3*beta3*F3b3*rhos12*power(F3,-1) + 
		10*alpha3*beta3*F3b3*rhos13*rhos23*power(F3,-1))*power(S,-1)*power(Sb,2)) + 
        power(alpha3,3)*power(F3b3,3)*(rhos13*(nu3*(nu3 + nu1*rhov13) + alpha3*beta3*F3b3*(2*nu3*rho3 + nu1*rhoc31)*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*(nu3*rhoc13 + alpha3*beta3*F3b3*rhos13*power(F3,-1)) + 
		power(alpha3,2)*power(beta3,2)*power(F3,-2)*power(F3b3,2))*power(Sb,2) + 
		3*rhos13*power(alpha3,2)*power(beta,2)*power(F3b3,2)*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(nu3*rhoc13 + 3*nu3*rho3*rhos13 + nu1*rhoc31*rhos13 + 4*alpha3*beta3*F3b3*rhos13*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*power(rhos13,2))*power(S,-1)*power(Sb,2)) - 
        power(alpha2,3)*power(F2b2,3)*(rhos12*(nu2*(nu2 + nu1*rhov12) + alpha2*beta2*F2b2*(2*nu2*rho2 + nu1*rhoc21)*power(F2,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*(nu2*rhoc12 + alpha2*beta2*F2b2*rhos12*power(F2,-1)) + 
		power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2))*power(Sb,2) + 
		6*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(rhos12 + 2*rhos13*rhos23 + 2*rhos12*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(nu1*rhoc31*rhos12 + 3*nu2*rhoc32*rhos12 + 2*nu2*rho2*rhos13 + nu1*rhoc21*rhos13 + nu3*rhoc23*rhos13 + 
		3*nu2*rhoc12*rhos23 + nu3*rhoc13*rhos23 + 4*nu2*rho2*rhos12*rhos23 + 2*nu1*rhoc21*rhos12*rhos23 + 
		2*nu3*rhoc23*rhos12*rhos23 + 2*alpha1*beta1*F1b1*rhos12*(rhos13 + rhos12*rhos23)*power(F1,-1) + 
		2*alpha2*beta2*F2b2*(rhos13 + 5*rhos12*rhos23)*power(F2,-1) + 2*alpha3*beta3*F3b3*rhos13*rhos23*power(F3,-1) + 
		2*alpha3*beta3*F3b3*rhos12*power(F3,-1)*power(rhos23,2))*power(S,-1)*power(Sb,2)) + 
        alpha3*F3b3*power(alpha2,2)*power(F2b2,2)*((nu1*nu2*rhos13*rhov12 + nu1*nu2*rhos12*rhos23*rhov12 + nu1*nu3*rhos12*rhos23*rhov13 + 
		nu2*nu3*rhos13*rhov23 + nu2*nu3*rhos12*rhos23*rhov23 + 
		rhos12*rhos23*power(alpha2,2)*power(beta2,2)*power(F2,-2)*power(F2b2,2) + 
		alpha3*beta3*F3b3*nu2*rhoc32*rhos13*power(F3,-1) + alpha3*beta3*F3b3*nu1*rhoc31*rhos12*rhos23*power(F3,-1) + 
		alpha3*beta3*F3b3*nu2*rhoc32*rhos12*rhos23*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*(nu2*rhoc12*rhos13 + nu2*rhoc12*rhos12*rhos23 + nu3*rhoc13*rhos12*rhos23 + 
		alpha2*beta2*F2b2*rhos12*(rhos13 + rhos12*rhos23)*power(F2,-1) + alpha3*beta3*F3b3*rhos12*rhos13*rhos23*power(F3,-1)) + 
		alpha2*beta2*F2b2*power(F2,-1)*(nu3*rhoc23*rhos13 + 2*nu2*rho2*rhos12*rhos23 + nu3*rhoc23*rhos12*rhos23 + 
		nu1*rhoc21*(rhos13 + rhos12*rhos23) + alpha3*beta3*F3b3*rhos23*(rhos13 + rhos12*rhos23)*power(F3,-1)) + 
		rhos12*rhos23*power(nu2,2))*power(Sb,2) + 
		6*power(alpha3,2)*power(beta,2)*power(F3b3,2)*(rhos13 + 2*rhos12*rhos23 + 2*rhos13*power(rhos23,2))*power(S,-2)*power(Sb,2) + 
		alpha3*beta*F3b3*(nu2*rhoc12 + nu3*rhoc13 + nu2*rho2*rhos12 + nu1*rhoc21*rhos12 + 2*nu3*rhoc23*rhos12 + nu3*rho3*rhos13 + 
		nu1*rhoc31*rhos13 + 2*nu2*rhoc32*rhos13 + 2*nu3*rho3*rhos12*rhos23 + 2*nu1*rhoc31*rhos12*rhos23 + 
		4*nu2*rhoc32*rhos12*rhos23 + 2*nu2*rho2*rhos13*rhos23 + 2*nu1*rhoc21*rhos13*rhos23 + 4*nu3*rhoc23*rhos13*rhos23 + 
		2*alpha3*beta3*F3b3*rhos13*power(F3,-1) + 4*alpha3*beta3*F3b3*rhos12*rhos23*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*(4*rhos12*rhos13*rhos23 + power(rhos12,2) + power(rhos13,2)) + 
		2*nu2*rhoc12*power(rhos23,2) + 2*nu3*rhoc13*power(rhos23,2) + 6*alpha3*beta3*F3b3*rhos13*power(F3,-1)*power(rhos23,2) + 
		2*alpha2*beta2*F2b2*power(F2,-1)*(rhos12 + 2*rhos13*rhos23 + 3*rhos12*power(rhos23,2)))*power(S,-1)*power(Sb,2))))*
		power(2*alpha1*F1b1*(alpha2*F2b2*rhos12 - alpha3*F3b3*rhos13) - 2*alpha2*alpha3*F2b2*F3b3*rhos23 + power(alpha1,2)*power(F1b1,2) + 
		power(alpha2,2)*power(F2b2,2) + power(alpha3,2)*power(F3b3,2),-1);

	*covarSABR=(-(beta*Sb*power(alpha1,4)*power(F1b1,4)*power(S,-1)) - beta*Sb*power(alpha2,4)*power(F2b2,4)*power(S,-1) - 
		power(alpha3,3)*power(F3b3,3)*(Sb*(nu3*rho3 + alpha3*beta3*F3b3*power(F3,-1)) + alpha3*beta*F3b3*Sb*power(S,-1)) + 
		power(alpha1,3)*power(F1b1,3)*(Sb*(nu1*rho1 + alpha1*beta1*F1b1*power(F1,-1)) + 
		4*beta*(-(alpha2*F2b2*rhos12) + alpha3*F3b3*rhos13)*Sb*power(S,-1)) + 
		power(alpha2,3)*power(F2b2,3)*(Sb*(nu2*rho2 + alpha2*beta2*F2b2*power(F2,-1)) + 4*alpha3*beta*F3b3*rhos23*Sb*power(S,-1)) + 
		alpha2*F2b2*power(alpha3,2)*power(F3b3,2)*(Sb*(nu3*rhoc23 + nu3*rho3*rhos23 + nu2*rhoc32*rhos23 + 
		2*alpha3*beta3*F3b3*rhos23*power(F3,-1) + alpha2*beta2*F2b2*power(F2,-1)*power(rhos23,2)) + 
		4*alpha3*beta*F3b3*rhos23*Sb*power(S,-1)) - alpha3*F3b3*power(alpha2,2)*power(F2b2,2)*
		(Sb*(nu2*rhoc32 + nu2*rho2*rhos23 + nu3*rhoc23*rhos23 + 2*alpha2*beta2*F2b2*rhos23*power(F2,-1) + 
		alpha3*beta3*F3b3*power(F3,-1)*power(rhos23,2)) + 2*alpha3*beta*F3b3*Sb*(1 + 2*power(rhos23,2))*power(S,-1)) + 
		alpha1*F1b1*(-4*beta*rhos12*Sb*power(alpha2,3)*power(F2b2,3)*power(S,-1) + 
		power(alpha3,2)*power(F3b3,2)*(Sb*(nu3*rhoc13 + nu3*rho3*rhos13 + nu1*rhoc31*rhos13 + 2*alpha3*beta3*F3b3*rhos13*power(F3,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*power(rhos13,2)) + 4*alpha3*beta*F3b3*rhos13*Sb*power(S,-1)) + 
		power(alpha2,2)*power(F2b2,2)*(Sb*(nu2*rhoc12 + nu2*rho2*rhos12 + nu1*rhoc21*rhos12 + 2*alpha2*beta2*F2b2*rhos12*power(F2,-1) + 
		alpha1*beta1*F1b1*power(F1,-1)*power(rhos12,2)) + 4*alpha3*beta*F3b3*(rhos13 + 2*rhos12*rhos23)*Sb*power(S,-1)) - 
		alpha2*alpha3*F2b2*F3b3*(Sb*(nu1*rhoc31*rhos12 + nu2*rhoc32*rhos12 + nu1*rhoc21*rhos13 + nu3*rhoc23*rhos13 + nu2*rhoc12*rhos23 + 
		nu3*rhoc13*rhos23 + 2*alpha1*beta1*F1b1*rhos12*rhos13*power(F1,-1) + 2*alpha2*beta2*F2b2*rhos12*rhos23*power(F2,-1) + 
		2*alpha3*beta3*F3b3*rhos13*rhos23*power(F3,-1)) + 4*alpha3*beta*F3b3*(rhos12 + 2*rhos13*rhos23)*Sb*power(S,-1))) + 
		power(alpha1,2)*power(F1b1,2)*(-2*beta*Sb*power(alpha2,2)*power(F2b2,2)*(1 + 2*power(rhos12,2))*power(S,-1) + 
		alpha2*F2b2*(Sb*(nu1*rhoc21 + nu1*rho1*rhos12 + nu2*rhoc12*rhos12 + 2*alpha1*beta1*F1b1*rhos12*power(F1,-1) + 
		alpha2*beta2*F2b2*power(F2,-1)*power(rhos12,2)) + 4*alpha3*beta*F3b3*(2*rhos12*rhos13 + rhos23)*Sb*power(S,-1)) - 
		alpha3*F3b3*(Sb*(nu1*rhoc31 + nu1*rho1*rhos13 + nu3*rhoc13*rhos13 + 2*alpha1*beta1*F1b1*rhos13*power(F1,-1) + 
		alpha3*beta3*F3b3*power(F3,-1)*power(rhos13,2)) + 2*alpha3*beta*F3b3*Sb*(1 + 2*power(rhos13,2))*power(S,-1))))*power(Sb,-2)*
		power(2*alpha1*F1b1*(alpha2*F2b2*rhos12 - alpha3*F3b3*rhos13) - 2*alpha2*alpha3*F2b2*F3b3*rhos23 + power(alpha1,2)*power(F1b1,2) + 
		power(alpha2,2)*power(F2b2,2) + power(alpha3,2)*power(F3b3,2),-0.5);

	return;
	}
		
		
	double TriSABR_VanillaOption(
		double F1,
		double alpha1,
		double beta1,
		double rho1,
		double nu1,
		
		double F2,
		double alpha2,
		double beta2,
		double rho2,
		double nu2,
		
		double F3,
		double alpha3,
		double beta3,
		double rho3,
		double nu3,
		
		double rhos12,
		double rhos23,
		double rhos13,
		double rhov12,
		double rhov23,
		double rhov13,
		double rhoc12,
		double rhoc21,
		double rhoc23,
		double rhoc32,
		double rhoc13,
		double rhoc31,
		double K,
		double T,
		int callput,
		int flag,
		double nbsteps)
	{
		double alpha,beta,rho,nu,F,Fbeta;
		F=F1+F2-F3;
		beta=beta1;			/// Hypothese raisonable

		double Stddev, varSABR,covarSABR;
		TriSABRSumOptionObservables(
			F1, alpha1, beta1,  rho1, nu1,
			F2, alpha2, beta2,  rho2, nu2,
			F3, alpha3, beta3,  rho3, nu3,
			rhos12,  rhos23,  rhos13,
			rhov12,  rhov23,  rhov13,
			rhoc12,  rhoc13,
			rhoc21,  rhoc23,
			rhoc31,  rhoc32,
			beta,
			&Stddev, &varSABR, &covarSABR);

		Fbeta=pow(F,beta);
		alpha=Stddev/Fbeta;
		nu=sqrt(varSABR)/alpha;
		rho=covarSABR/(alpha*alpha*nu*Fbeta);
		return Export_SABR_VanillaOption(F, K,T,alpha,beta,rho,nu, callput, flag, nbsteps);
	}
		

		





CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/