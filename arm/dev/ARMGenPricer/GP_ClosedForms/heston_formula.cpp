/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.cpp
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

#include "gpclosedforms/heston.h"
#include "gpclosedforms/Simple_Heston.h"

#include "gpclosedforms/heston_formula.h"

#include <glob/expt.h>

///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=(theta(longtermVar-V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///			V(0)=Var0
///		
/////////////////////////////////////////////////////////////////////////////:


CC_BEGIN_NAMESPACE(ARM)



//////////////////////////////////////////////////////////////////////////////////////////
///
///      Simple Heston 
///
//////////////////////////////////////////////////////////////////////////////////////////

double ARM_CF_Heston_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_Heston_JumpDiffusion_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_Formula::value  : bad argsize");
	 }
	 int callput=a[ARM_CF_Heston_JumpDiffusion_Formula::CALLORPUT];
	 int nb1=4;
	 int NbStage=4;
	 int NbOscill=7;
	 double prec=1e-7;

switch (callput)
	{
	case K_CALL :
		{
			return Simple_Heston(
				a[ARM_CF_Heston_Formula::INDEX],					//INDEX			
				a[ARM_CF_Heston_Formula::STRIKE],					//STRIKE
				a[ARM_CF_Heston_Formula::TIMETOMATURITY],			//TIMETOMATURITY
				a[ARM_CF_Heston_Formula::INITIALVOLSQUARE],				//INITIALVOL
				a[ARM_CF_Heston_Formula::VOLSQUAREREVERTINGSPEED],//THETA	
				a[ARM_CF_Heston_Formula::LONGTERMVOLSQUARE],		//LONGTERMVOL	
				a[ARM_CF_Heston_Formula::VOLSQUAREVOLATILITY],	//KSI	
				a[ARM_CF_Heston_Formula::CORRELATION],			//RHO
				a[ARM_CF_Heston_Formula::NBTERMS],					//NBTERMS
				nb1,
				NbStage,
				NbOscill,
				prec
				);
			break;
		}
	case K_PUT :
		{
			return Simple_Heston(
				a[ARM_CF_Heston_Formula::INDEX],					//INDEX			
				a[ARM_CF_Heston_Formula::STRIKE],					//STRIKE
				a[ARM_CF_Heston_Formula::TIMETOMATURITY],			//TIMETOMATURITY
				a[ARM_CF_Heston_Formula::INITIALVOLSQUARE],				//INITIALVOL
				a[ARM_CF_Heston_Formula::VOLSQUAREREVERTINGSPEED],//THETA	
				a[ARM_CF_Heston_Formula::LONGTERMVOLSQUARE],		//LONGTERMVOL	
				a[ARM_CF_Heston_Formula::VOLSQUAREVOLATILITY],	//KSI	
				a[ARM_CF_Heston_Formula::CORRELATION],			//RHO
				a[ARM_CF_Heston_Formula::NBTERMS],					//NBTERMS	
				nb1,												// NB terms for the first stage
				NbStage,
				NbOscill,											// NB os oscillation handled per stage
				prec												// precision to be obtained before leaving the 
				)-(a[ARM_CF_Heston_Formula::INDEX]-a[ARM_CF_Heston_Formula::STRIKE]);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_Formula : callput , bad input :");
			break;
		}
	}
}



////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return standard_first_derivative(ARM_CF_Heston_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_Heston_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_Heston_Formula::INDEX :
		return 0.00001;		// INDEX
	case ARM_CF_Heston_Formula::STRIKE :
		return 0.00001;		// STRIKE
	case ARM_CF_Heston_Formula::INITIALVOLSQUARE :
		return 0.00001;		// INITIALINITIALVOLSQUAREVOL
	case ARM_CF_Heston_Formula::TIMETOMATURITY :
		return 0.00001;		// TIMETOMATURITY
	case ARM_CF_Heston_Formula::LONGTERMVOLSQUARE :  
		return 0.00001;		// LONGTERMVOL
	case ARM_CF_Heston_Formula::VOLSQUAREREVERTINGSPEED :
		return 0.00001;		// VOLSQUAREREVERTINGSPEED
	case ARM_CF_Heston_Formula::VOLSQUAREVOLATILITY :
		return 0.00001;		// VOLSQUAREVOLATILITY
	case ARM_CF_Heston_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_Heston_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_Heston_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
	if(a[ARM_CF_Heston_Formula::INDEX]<=0) return ArgumentList_Checking_Result(false," INDEX  not positive");					//			positivity of INDEX
	if(a[ARM_CF_Heston_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE  not positive");					//			positivity of STRIKE
	if(a[ARM_CF_Heston_Formula::INITIALVOLSQUARE]<=0) return ArgumentList_Checking_Result(false," INITIALVOL   not positive");		//			positivity of INITIALVOL
	if(a[ARM_CF_Heston_Formula::TIMETOMATURITY]<=0) return ArgumentList_Checking_Result(false," TIMETOMATURITY   not positive");		//			positivity of TIMETOMATURITY
	if(a[ARM_CF_Heston_Formula::LONGTERMVOLSQUARE]<=0) return ArgumentList_Checking_Result(false," LONGTERMVOLSQUARE   not positive or zero");	//			positivity of LONGTERMVOL
	if(a[ARM_CF_Heston_Formula::VOLSQUAREREVERTINGSPEED]<0) return ArgumentList_Checking_Result(false," VOLSQUAREREVERTINGSPEED   not positive or zero");	//			positivity of VOLSQUAREREVERTINGSPEED
	if(a[ARM_CF_Heston_Formula::VOLSQUAREVOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLSQUAREVOLATILITY   not positive or zero");	//			positivity of VOLSQUAREVOLATILITY
	if (fabs(a[ARM_CF_Heston_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1
	 if ((a[ARM_CF_Heston_Formula::CALLORPUT]!=K_CALL) && (a[ARM_CF_Heston_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); //	
	
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_Heston_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>=ARM_CF_Heston_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=theta*(lgtvol^2-V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///			V(0)=sig^2
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:

//////////////////////////////////////////////////////////////////////////////////////////
///
///      Generalized Heston 
///
//////////////////////////////////////////////////////////////////////////////////////////

double ARM_CF_Heston_JumpDiffusion_Formula::value(const ArgumentList& a)
	{
	 int argsize=a.size();
	 if (argsize!=ARM_CF_Heston_JumpDiffusion_Formula::Nb_Parameters)
	 {
		 throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_JumpDiffusion_Formula::value  : bad argsize");
	 }
	 int callput=a[ARM_CF_Heston_JumpDiffusion_Formula::CALLORPUT];

switch (callput)
	{
	case K_CALL :
		{
			return GeneralizedHeston(
				a[ARM_CF_Heston_JumpDiffusion_Formula::INDEX],					//INDEX			
				a[ARM_CF_Heston_JumpDiffusion_Formula::STRIKE],					//STRIKE
				a[ARM_CF_Heston_JumpDiffusion_Formula::INITIALVOL],				//INITIALVOL
				a[ARM_CF_Heston_JumpDiffusion_Formula::TIMETOMATURITY],			//TIMETOMATURITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::LONGTERMVOL],		//LONGTERMVOL	
				a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREREVERTINGSPEED],//THETA	
				a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREVOLATILITY],	//KSI	
				a[ARM_CF_Heston_JumpDiffusion_Formula::CORRELATION],			//RHO
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPMEAN],				//JUMPMEAN		
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPVOLATILITY],			//JUMPVOLATILITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPPROBABILITY],		//JUMPPROBABILITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::NBTERMS]					//NBTERMS			
				);
			break;
		}
	case K_PUT :
		{
			return GeneralizedHeston(
				a[ARM_CF_Heston_JumpDiffusion_Formula::INDEX],					//INDEX			
				a[ARM_CF_Heston_JumpDiffusion_Formula::STRIKE],					//STRIKE
				a[ARM_CF_Heston_JumpDiffusion_Formula::INITIALVOL],				//INITIALVOL
				a[ARM_CF_Heston_JumpDiffusion_Formula::TIMETOMATURITY],			//TIMETOMATURITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::LONGTERMVOL],		//LONGTERMVOL	
				a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREREVERTINGSPEED],//THETA	
				a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREVOLATILITY],	//KSI	
				a[ARM_CF_Heston_JumpDiffusion_Formula::CORRELATION],			//RHO
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPMEAN],				//JUMPMEAN		
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPVOLATILITY],			//JUMPVOLATILITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPPROBABILITY],		//JUMPPROBABILITY
				a[ARM_CF_Heston_JumpDiffusion_Formula::NBTERMS]					//NBTERMS			
				)-(a[ARM_CF_Heston_JumpDiffusion_Formula::INDEX]-a[ARM_CF_Heston_JumpDiffusion_Formula::STRIKE]);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_JumpDiffusion_Formula : callput , bad input :");
			break;
		}
	}
}



////////////////////////////////////////////////////////////////////////
///
///      First Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_JumpDiffusion_Formula::value(int i,const ArgumentList& a, double s)
{
	 		return standard_first_derivative(ARM_CF_Heston_JumpDiffusion_Formula::value,i,a,s);
}



////////////////////////////////////////////////////////////////////////
///
///      Second Derivative
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_JumpDiffusion_Formula::value(int i,int j,const ArgumentList& a, double s1,double s2)
{
	 		return standard_second_derivative(ARM_CF_Heston_JumpDiffusion_Formula::value,i,j,a,s1,s2);
}


////////////////////////////////////////////////////////////////////////
///
///      Description of the shifting used for the derivation
///
////////////////////////////////////////////////////////////////////////


double ARM_CF_Heston_JumpDiffusion_Formula::specific_shift(int i) {
	switch(i)
	{
		///		Here the limitation of the dependence in the smile and the copula is that they 
		///	are described by alpha,beta,rho,nu   and the copula described by one number (called here correlation)
		///
		
	case ARM_CF_Heston_JumpDiffusion_Formula::INDEX :
		return 0.00001;		// INDEX
	case ARM_CF_Heston_JumpDiffusion_Formula::STRIKE :
		return 0.00001;		// STRIKE
	case ARM_CF_Heston_JumpDiffusion_Formula::INITIALVOL :
		return 0.00001;		// INITIALVOL
	case ARM_CF_Heston_JumpDiffusion_Formula::TIMETOMATURITY :
		return 0.00001;		// TIMETOMATURITY
	case ARM_CF_Heston_JumpDiffusion_Formula::LONGTERMVOL :  
		return 0.00001;		// LONGTERMVOL
	case ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREREVERTINGSPEED :
		return 0.00001;		// VOLSQUAREREVERTINGSPEED
	case ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREVOLATILITY :
		return 0.00001;		// VOLSQUAREVOLATILITY
	case ARM_CF_Heston_JumpDiffusion_Formula::CORRELATION :
		return 0.00001;		// CORRELATION
	case ARM_CF_Heston_JumpDiffusion_Formula::JUMPPROBABILITY :  
		return 0.00001;		// JUMPPROBABILITY
	case ARM_CF_Heston_JumpDiffusion_Formula::JUMPMEAN :
		return 0.00001;		// JUMPMEAN
	case ARM_CF_Heston_JumpDiffusion_Formula::JUMPVOLATILITY :
		return 0.00001;		// JUMPVOLATILITY
	default :
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"ARM_CF_Heston_JumpDiffusion_Formula::specific_shift : incorrect input");
	}
}



////////////////////////////////////////////////////////////////////////
///
///      Checking of the arguments : description of the domain of validity
///
////////////////////////////////////////////////////////////////////////


ArgumentList_Checking_Result ARM_CF_Heston_JumpDiffusion_Formula::check_argument(const ArgumentList& a)
{
	///		Here the limitation of the dependence in the smile and the copula is that they 
	///	are described by alpha,beta,nu  all positive number and rho, included in {-1,1}
	///  and the copula described by one number (called here correlation)
	///
	
	if(a.size()!=ARM_CF_Heston_JumpDiffusion_Formula::Nb_Parameters) 
	{
		return ArgumentList_Checking_Result(false,"Bad number of arguments  ");
	}
	
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::INDEX]<=0) return ArgumentList_Checking_Result(false," INDEX  not positive");					//			positivity of INDEX
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::STRIKE]<=0) return ArgumentList_Checking_Result(false," STRIKE  not positive");					//			positivity of STRIKE
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::INITIALVOL]<=0) return ArgumentList_Checking_Result(false," INITIALVOL   not positive");		//			positivity of INITIALVOL
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::TIMETOMATURITY]<=0) return ArgumentList_Checking_Result(false," TIMETOMATURITY   not positive");		//			positivity of TIMETOMATURITY
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::LONGTERMVOL]<=0) return ArgumentList_Checking_Result(false," LONGTERMVOL   not positive or zero");	//			positivity of LONGTERMVOL
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREREVERTINGSPEED]<0) return ArgumentList_Checking_Result(false," VOLSQUAREREVERTINGSPEED   not positive or zero");	//			positivity of VOLSQUAREREVERTINGSPEED
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::VOLSQUAREVOLATILITY]<=0) return ArgumentList_Checking_Result(false," VOLSQUAREVOLATILITY   not positive or zero");	//			positivity of VOLSQUAREVOLATILITY
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPPROBABILITY]<0) return ArgumentList_Checking_Result(false," JUMPPROBABILITY    not positive or 0");	//			positivity of JUMPPROBABILITY
	if(a[ARM_CF_Heston_JumpDiffusion_Formula::JUMPVOLATILITY]<0) return ArgumentList_Checking_Result(false," JUMPVOLATILITY   not positive ");		//			positivity of JUMPVOLATILITY
	if (fabs(a[ARM_CF_Heston_JumpDiffusion_Formula::CORRELATION])>1.) return ArgumentList_Checking_Result(false," Abs(CORRELATION) not <=1"); //		Abs(CORRELATION)<=1
	 if ((a[ARM_CF_Heston_JumpDiffusion_Formula::CALLORPUT]!=K_CALL) && (a[ARM_CF_Heston_JumpDiffusion_Formula::CALLORPUT]!=K_PUT))
		  return ArgumentList_Checking_Result(false," CALLORPUT should be CALL or PUT"); //	
	
	return ArgumentList_Checking_Result(true,string(""));
}


////////////////////////////////////////////////////////////////////////
///
///      Checking of the dimension rank with respect to which one want ot derive
///
////////////////////////////////////////////////////////////////////////


 ArgumentList_Checking_Result ARM_CF_Heston_JumpDiffusion_Formula::check_dimension(int rank)
{
	if ((rank<0)||(rank>=ARM_CF_Heston_JumpDiffusion_Formula::Nb_Derivable_Parameters)) return ArgumentList_Checking_Result(false," Invalide derivative dimension"); // Invalide derivative dimension	
	
	return ArgumentList_Checking_Result(true,string(""));
}


 

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
