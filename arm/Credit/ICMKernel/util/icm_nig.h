#ifndef _ICM_NIG_H_
#define _ICM_NIG_H_

/*********************************************************************************/
/*! \ "ICM_Nig.h"
 *  \author Fakher Ben Atig
 *	\version 1.0
 *	\date   October 2005
 *	\file   icm_nig.h
 *	\Normal Inverse Gaussian (NIG) functions
 *  \Creates a MC generator for the NIG Copula
/*****************************************************************************************/

#include "ICMKernel/util/icm_root_generator.h"
#include "ARMKernel\util\gaussian.h"
#include <stdio.h>
#include <nag.h>

//static double __stdcall FunctionToIntegrate (double y,Nag_User *comm);


class ICM_Nig : public ARM_Object
{

private :
	
	double	itsAlpha;	//1st parameter
	double	itsBetaNig;	//2nd parameter
	double	itsMu;		//3rd parameter
	double	itsDelta;	//4th parameter
	
	double	itsX	;	//Used for integration purpose

	double itsRho ;		//Correlation factor, itsRho * itsRho = Correlation

protected: 


public : 

	void Init()
	{
		itsAlpha = 0.;
		itsBetaNig	= 0. ;
		itsMu	= 0. ;
		itsDelta = 0.;
		itsX	=  0.;
		itsRho	= 0.;
	}

	ICM_Nig() { Init() ;}
	
	ICM_Nig (double	alpha, double BetaNig, double Mu, double Delta, double Rho) ;

	// Constructor for Nig Copula used in CDO pricing

	ICM_Nig (double	alpha, double BetaNig, double Rho) ;

	ICM_Nig(ICM_Nig & GenNig) ;

	void setParameters(double	alpha, double BetaNig, double Mu, double Delta, double Rho) ;

	void Set_NIG_Alpha(double alpha) { itsAlpha=alpha;} 
	void Set_NIG_Beta(double BetaNig) {itsBetaNig=BetaNig;} 
	void Set_NIG_Mu(double mu) {itsMu=mu;} 
	void Set_NIG_Delta(double delta) {itsDelta=delta;} 
	void SetX(double x) {itsX=x;} 

	void SetRho(double rho) {itsRho = rho;} 
	
	double Get_NIG_Alpha() {return itsAlpha;}
	double Get_NIG_Beta() {return itsBetaNig;}
	double Get_NIG_Mu() {return itsMu;}
	double Get_NIG_Delta() {return itsDelta;}
	double GetitsX() {return itsX;}

	double GetitsRho() {return itsRho;}

	double generateRandom() ;

	//Density function
	double GetDensity(double x) ;

	// Cumulative Distribution Function (CDF)
	double GetDistribution(double x) ;
	
	// Conditional Probability
	double getConditionalDefaultProbability(double level, double MktDefProb);

	//Inverse CDF of  NIG
	double NIG_Inverse(double y);

	//Implied Barrier (by the Market Default Probability of the issuer)
	double GetDefaultBarrier(double MktDefProb);

	double IntegrateByGaussLaguerre (double a, double b, double(__stdcall f)(double x,Nag_User *comm));

	
	~ICM_Nig() {};


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_Nig * GenNig = (ICM_Nig *) src ;

		Set_NIG_Alpha(GenNig->Get_NIG_Alpha()) ;
		Set_NIG_Beta(GenNig->Get_NIG_Beta()) ;
		Set_NIG_Mu(GenNig->Get_NIG_Mu()) ;
		Set_NIG_Delta(GenNig->Get_NIG_Delta());
	}


	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
 
		BitwiseCopy(src) ;
	}

	// --------------
	//	Clone Method
	// --------------

	ARM_Object * Clone(void)
	{
		ICM_Nig * theClone = new ICM_Nig() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

};


#endif
/*----------------------------------------------------------------------------*/
/*---- End of file ----*/