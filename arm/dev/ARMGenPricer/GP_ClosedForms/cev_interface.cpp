/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file SABR_Analytics.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"

#include "gpclosedforms/cev_formula.h"


CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV VanillaOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_VanillaOption(double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_VanillaOption_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  CEV VanillaOption
///

double Export_CEV_VanillaOption(int i,double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps)
{
	ArgumentList a(T,f,K,sig,beta,drift,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_VanillaOption_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  CEV VanillaOption
///

double Export_CEV_VanillaOption(int i,int j,double f, double K, double T,double drift, double sig, double beta,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_VanillaOption_Formula> y;
	return y(i,j,a);
}




/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV DoubleBarrierOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_DoubleBarrierOption(double f, double K, double T,double barrdown,double barrup, double drift, double sig, double beta,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,barrdown,barrup,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_DoubleBarrierOption_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_DoubleBarrierOption(int i,double f, double K, double T,double barrdown,double barrup, double drift, double sig, double beta,int callput, int nbsteps)
{
	ArgumentList a(T,f,K,sig,beta,drift,barrdown,barrup,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_DoubleBarrierOption_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  CEV DoubleBarrierOption
///

double Export_CEV_DoubleBarrierOption(int i,int j,double f, double K, double T,double barrdown,double barrup, double drift, double sig, double beta,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,barrdown,barrup,callput,nbsteps);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_CEV_DoubleBarrierOption_Formula> y;
	return y(i,j,a);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  CEV BarrierOption
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_CEV_BarrierOption(double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,barrier,optype,callput,nbsteps);  //  of the enum "ArgumentType" of the class 

	Power_Expression<ARM_CF_CEV_BarrierOption_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  CEV BarrierOption
///

double Export_CEV_BarrierOption(int i,double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,barrier,optype,callput,nbsteps);  //  of the enum "ArgumentType" of the class 

	Power_Expression<ARM_CF_CEV_BarrierOption_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  CEV BarrierOption
///

double Export_CEV_BarrierOption(int i,int j,double f, double K, double T,double barrier, double drift, double sig, double beta, int optype,int callput, int nbsteps)

{
	ArgumentList a(T,f,K,sig,beta,drift,barrier,optype,callput,nbsteps);  //  of the enum "ArgumentType" of the class 
	Power_Expression<ARM_CF_CEV_BarrierOption_Formula> y;
	return y(i,j,a);
}






CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/