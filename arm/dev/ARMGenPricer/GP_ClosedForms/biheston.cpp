/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file biheston.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date Avril 2007
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include <cmath>
#include <complex>

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/oscillatory_integrals.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/vanilla_normal.h"

#include "gpbase/numericconstant.h"

using std::sqrt;
using std::log;
using std::real;
using std::imag;
using std::exp;
using std::norm;
using std::complex;

//using std::complex<double>;
//using std::complex<long double>;

///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S dW1
///			dV=lambda*(Vinf-V)dt +nu*V1^(1/2) dW2 
///			dW1.dW2=rho1*dt
///
///			dS2= S2 dW3
///			dV2=lambda2*(Vinf2-V2)dt +nu2*V2^(1/2) dW4 
///			dW3.dW4=rho2*dt
///
///			dW1.dW3=rhos*dt
///			dW1.dW4=rhoc12*dt
///			dW2.dW3=rhoc21*dt
///			dW2.dW4=rho2v*dt

/////////////////////////////////////////////////////////////////////////////:

CC_BEGIN_NAMESPACE(ARM)

double BiHeston_VanillaOption(
									 double		F1,
									 double		V1,
									 double		Vinfini1,
									 double		lambda1,
									 double		nu1,
									 double		rho1,
									 double		F2,
									 double		V2,
									 double		Vinfini2,
									 double		lambda2,
									 double		nu2,
									 double		rho2,
									 double		rhos,
									 double		rhov,
									 double		rhoc12,
									 double		rhoc21,
									 double		k,
									 double		T,
									 double		CallPut,
									 double		LambdaB,
									 double		Flag,
									 int nbfirst,
									 int nb,
									 int NbStage,
									 int NbOscill,
									 double prec
									 )
{
	/*
		int nbfirst=40;
		int nb=40;
		int NbStage=10;
		int NbOscill=5;
		double prec=1.0e-12;
		*/

	double S, V, VarianceVariance, CovarianceVariance, LongTermVariance1,LongTermVariance2,
		LongTermCorrelation, LongTermVarianceVariance, NuS, RhoS, VinfiniS, LambdaS, 
		SQRV,rhos2,LongTermVarianceVariance0;
	rhos2=rhos*rhos;
	S = F1 - F2;
	SQRV =sqrt(V1*V2);
	V = F1*F1*V1 + F2*F2*V2 - 2.0*F1*F2*SQRV*rhos;

	double aux1,aux2,aux30,aux31,aux32,bx1,bx2,bx3;
	aux1=2.0*F1*F1*F1*F2*V1*rhos*(4.0*V1*V1*V2+V2*nu1*nu1+4.0*V1*V2*nu1*rho1+2.0*V1*SQRV*nu2*rhoc12+
		2.0*V2*SQRV*nu1*rhoc21+4.0*V1*V2*SQRV*rhos+SQRV*nu1*nu2*rhov);
	aux2=2.0*F1*F2*F2*F2*V2*rhos*(V1*(4.0*V2*V2+nu2*(nu2+2.0*SQRV*rhoc12)+
		4.0*V2*(nu2*rho2+SQRV*rhos))+SQRV*nu1*(2.0*V2*rhoc21+nu2*rhov));
	aux30=4.0*V1*V1*V2*(nu2*rhoc12*(1.0+rhos*rhos)+rhos*(SQRV*rhos+2.0*V2*(1.0+rhos*rhos)));
	aux31=4.0*V2*V2*(SQRV*rhos*rhos+nu1*rhoc21*(1.0+rhos*rhos));

	aux32=2.0*V2*(2.0*SQRV*nu2*rho2*rhos*rhos+nu1*(2.0*SQRV*rho1*rhos*rhos+nu2*(1.0+rhos*rhos)*rhov));

	VarianceVariance=1/SQRV*(F1*F1*F1*F1*V1*SQRV*(4.0*V1*V1+nu1*nu1+4.0*V1*nu1*rho1)+
		F2*F2*F2*F2*V2*SQRV*(4.0*V2*V2+nu2*nu2+4.0*V2*nu2*rho2)-aux1-aux2+F1*F1*F2*F2*(V2*SQRV*nu1*nu1*rhos*rhos
		+aux30+V1*(SQRV*nu2*nu2*rhos*rhos+aux31+aux32)));

	bx1=F1*F1*V1*(2.0*F1*V1-2.0*F2*SQRV*rhos)-F1*F2*SQRV*rhos*(2.0*F1*V1-2.0*F2*SQRV*rhos);
	bx2=-F2*F2*V2*(2.0*F2*V2-2.0*F1*SQRV*rhos)+F1*F2*SQRV*rhos*(2.0*F2*V2-2.0*F1*SQRV*rhos);
	bx3=F1*V1*nu1*rho1*(F1*F1-F1*F2*V2*rhos/SQRV)-F2*SQRV*nu1*rhoc21*(F1*F1-F1*F2*V2*rhos/SQRV);
	CovarianceVariance=-F2*V2*nu2*rho2*(F2*F2-F1*F2*V1*rhos/SQRV)+F1*SQRV*nu2*rhoc12*(F2*F2-F1*F2*V1*rhos/SQRV)+bx1+bx2+bx3;


	
	
	NuS=sqrt(VarianceVariance/V);
	RhoS=CovarianceVariance/(V*sqrt(VarianceVariance/V));
	
	if(Flag==0)
	{
		LongTermVariance1=nu1*nu1*V1/(2.0*lambda1);
		LongTermVariance2=nu2*nu2*V2/(2.0*lambda2);
		LongTermCorrelation=rhov*sqrt(lambda1*lambda2)/(lambda1+lambda2);
		LongTermVarianceVariance0=2*F1*F1*Vinfini1 + 2*F1*F2*rhos*Vinfini2;
		LongTermVarianceVariance=LongTermVariance2*(2*F1*F2*rhos*Vinfini1 + 2*F2*F2*Vinfini2) + 
			2*LongTermCorrelation*sqrt(LongTermVariance1*LongTermVariance2)*(2*F1*F2*rhos*Vinfini1 + 2*F2*F2*Vinfini2)*
			(2*F1*F1*Vinfini1 + 2*F1*F2*rhos*Vinfini2) + LongTermVariance1*LongTermVarianceVariance0*LongTermVarianceVariance0;
		VinfiniS=  sqrt(F1*F1*F1*F1*Vinfini1*Vinfini1 - (2*F1*F1*F2*F2*sqrt(lambda1*lambda2)*rhos*Vinfini1*Vinfini2)/
			(lambda1 + lambda2) + F2*F2*F2*F2*Vinfini2*Vinfini2);
		LambdaS= VarianceVariance/(2.*LongTermVarianceVariance);
		return NormalHeston(RhoS,LambdaS,VinfiniS,NuS,V,S,k, T,CallPut,LambdaB,
			nbfirst, nb, NbStage,  NbOscill, prec);
	}
	if(Flag==1)
	{
		VinfiniS=  sqrt(F1*F1*F1*F1*Vinfini1*Vinfini1 - (4.0*F1*F1*F2*F2*sqrt(lambda1*lambda2)*rhos*Vinfini1*Vinfini2)/
			(lambda1 + lambda2) + F2*F2*F2*F2*Vinfini2*Vinfini2);
		LambdaS=0.5*(lambda1+lambda2);
		return NormalHeston(RhoS,LambdaS,VinfiniS,NuS,V,S,k, T,CallPut,LambdaB,
			nbfirst, nb, NbStage,  NbOscill, prec);
	}
	if(Flag==2)
	{
		return VanillaOption_N(S,sqrt(V),k, T, CallPut);
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BiHeston_VanillaOption : flag , UnImplemented");
		
}

double BiShiftedHeston_VanillaOption(
									 double		F1,
									 double		V1,
									 double		Vinfini1,
									 double		lambda1,
									 double		nu1,
									 double		rho1,
									 double		gamma1,
									 double		F2,
									 double		V2,
									 double		Vinfini2,
									 double		lambda2,
									 double		nu2,
									 double		rho2,
									 double		gamma2,
									 double		rhos,
									 double		rhov,
									 double		rhoc12,
									 double		rhoc21,
									 double		k,
									 double		T,
									 double		CallPut,
									 double		LambdaB,
									 double		Flag,
									 int nbfirst,
									 int nb,
									 int NbStage,
									 int NbOscill,
									 double prec
									 )
{
	return BiHeston_VanillaOption(
									 	F1+gamma1,
									 	V1,
									 	Vinfini1,
									 	lambda1,
									 	nu1,
									 	rho1,
									 	F2+gamma2,
									 	V2,
									 	Vinfini2,
									 	lambda2,
									 	nu2,
									 	rho2,
									 	rhos,
									 	rhov,
									 	rhoc12,
									 	rhoc21,
									 	k+gamma1-gamma2,
									 	T,
									 	CallPut,
										LambdaB,
										Flag,
										nbfirst,
										nb,
										NbStage,
										NbOscill,
										prec
									 );
}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
