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
#include <glob/firsttoinc.h>
#include "gpbase/port.h"


#include "gpclosedforms/merton_formula.h"


CC_BEGIN_NAMESPACE(ARM)



 ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



/// callput =  1 (K_CALL) for call
/// callput = -1 (K_PUT) for put

double Export_Merton_JumpDiffusion(double F,
						double K,
						double t,
						double sigma,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb)
{
	ArgumentList a(F,K,t,sigma,lambda,muJ,sigmaJ,callorput,nb);

	Power_Expression<ARM_CF_Merton_JumpDiffusion_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Merton_JumpDiffusion(int i,
						double F,
						double K,
						double t,
						double sigma,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb)
{
	ArgumentList a(F,K,t,sigma,lambda,muJ,sigmaJ,callorput,nb);
	
	Power_Expression<ARM_CF_Merton_JumpDiffusion_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_Merton_JumpDiffusion(int i,int j,
							double F,
							double K,
							double t,
							double sigma,
							double lambda,
							double muJ, 
							double sigmaJ,
							int callorput,int nb)
{
	ArgumentList a(F,K,t,sigma,lambda,muJ,sigmaJ,callorput,nb);
	
	Power_Expression<ARM_CF_Merton_JumpDiffusion_Formula> y;
	return y(i,j,a);
}




CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
