/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file change_numeraire.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/tri_spreadoption_lognormal_formula.h"

 

CC_BEGIN_NAMESPACE(ARM)


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value Tri Spread Option
///		which pays   (A0+A1*S1+A2*S2+A3*S3)^+
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadOption(double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,K_CALL,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  Tri Spread Option
///

double Export_LogNormal_TriSpreadOption(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n)
{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,K_CALL,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  Tri Spread Option
///

double Export_LogNormal_TriSpreadOption(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
										double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
										double a0,double a1,double a2,double a3,double t,int n)
										{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,K_CALL,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption_Formula> y;
	return y(i,j,a);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  Digital Tri Spread Digital option
///		which pays  1  and A0+A1*S1+A2*S2+A3*S3>0
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadDigitalOption(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int callput,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,callput,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int callput,int n)
{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,callput,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a1,double a2,double a3,double t,int callput,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a1,a2,a3,callput,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption_Formula> y;
	return y(i,j,a);
}




 

/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value  Digital Tri Spread Digital option2 
///		which pays  1 if B0+B1*S1>0 and A0+A2*S2+A3*S3>0
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadDigitalOption2(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption2_Formula> y;
	return y(a);
}


///  
///			1st Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption2(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)
{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption2_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives  Digital Tri Spread Digital option
///

double Export_LogNormal_TriSpreadDigitalOption2(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadDigitalOption2_Formula> y;
	return y(i,j,a);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value   Tri Spread  option2 
///		which pays  A0+A2*S2+A3*S3 if B0+B1*S1>0 and A0+A2*S2+A3*S3>0
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadOption2(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption2_Formula> y;
	return y(a);
}


///  
///			1st Derivatives   Tri Spread  option 2
///

double Export_LogNormal_TriSpreadOption2(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)
{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption2_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives   Tri Spread  option 2
///

double Export_LogNormal_TriSpreadOption2(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption2_Formula> y;
	return y(i,j,a);
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
///
///   Forward value   Tri Spread  option3 
///		which pays  B0+B1*S1 if B0+B1*S1>0 and A0+A2*S2+A3*S3>0
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////


double Export_LogNormal_TriSpreadOption3(double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption3_Formula> y;
	return y(a);
}


///  
///			1st Derivatives   Tri Spread  option 3
///

double Export_LogNormal_TriSpreadOption3(int i,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)
{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption3_Formula> y;
	return y(i,a);
}


///  
///			2nd Derivatives   Tri Spread  option 3
///

double Export_LogNormal_TriSpreadOption3(int i,int j,double S1,double S2,double S3,double sig1,double sig2,double sig3,
											   double rho12,double rho13,double rho23,double mu1,double mu2,double mu3,
											   double a0,double a2,double a3,double b0,double b1,double t,int n)

{
	ArgumentList a(t,S1,S2,S3,sig1,sig2,sig3,rho12,rho13,rho23,mu1,mu2,mu3,a0,a2,a3,b0,b1,n);  //  of the enum "ArgType" of the class 

	Power_Expression<ARM_CF_TriSpreadOption3_Formula> y;
	return y(i,j,a);
}





CC_END_NAMESPACE()




/*---------------------------------------------------------------------------*/
/*---- End of file ----*/