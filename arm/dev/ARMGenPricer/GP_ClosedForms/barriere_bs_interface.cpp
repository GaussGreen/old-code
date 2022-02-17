/*---------------------------------------------------------------------------*/
/*---- End of file ----*/





/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\Barriere_bs_interface.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */

#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/Barriere_bs_formula.h"





CC_BEGIN_NAMESPACE(ARM)




///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 

///
///////////////////////////////////////////////////////////////////////
double Export_BS_EuroBarriere(double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown)
{
	int optiontype;
	if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT;
	else
			if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT;
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BS_EuroBarriere : inout/updown : bad input :");

	ArgumentList a(f,k,b,r,v,t,rate,callput,optiontype);
	Power_Expression<ARM_CF_BS_SingleBarriere_Formula> y;
	return y(a);
}


double Export_BS_EuroBarriere_ImpliedVol(double f,double k, double b, double r, double opt, double t,double rate,
										 int callput, int inout, int updown)
{		
	struct PricingFunctionToInverse : public DoubleToDoubleFunc 
	{
		double f0;
		double k0;
		double b0;
		double r0; 
		double t0;
		double rate0;
		int callput0; 
		int inout0; 
		int updown0;
		PricingFunctionToInverse(double fa,double ka, double ba, double ra, double ta,double ratea,
								   int callputa, int inouta, int updowna):
		f0(fa),k0(ka),b0(ba),r0(ra),t0(ta),rate0(ratea),callput0(callputa),inout0(inouta),updown0(updowna)
		{}

		virtual double operator() (double v)  const
		{
			int optiontype;
			if((inout0==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown0==ARM_CF_BS_SingleBarriere_Formula::UP)) 
				optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN;
			else
				if((inout0==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown0==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
					optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN;
				else
					if((inout0==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown0==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
						optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT;
					else
						if((inout0==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown0==ARM_CF_BS_SingleBarriere_Formula::UP)) 
							optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT;
						else
							throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BS_EuroBarriere : inout/updown : bad input :");

			ArgumentList a(f0,k0,b0,r0,v,t0,rate0,callput0,optiontype);
			Power_Expression<ARM_CF_BS_SingleBarriere_Formula> y;
			return y(a);
		}
	};
	
	PricingFunctionToInverse x(f,k,b,r,t,rate,callput,inout,updown);
	double typical=0.2;
	double typical_change=0.02;
	double tolerance=1e-12;
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(opt,typical,typical_change,tolerance); 
};


double Export_BS_EuroDoubleBarriere(double f,double k, double bup, double bdown, double v, double t,double r,double b,
								   int callput)
{
	ArgumentList a(f,k,bup,bdown,v,t,r,b,callput);
	Power_Expression<ARM_CF_BS_DoubleBarriere_Formula> y;
	return y(a);
}


double Export_BS_EuroDoubleBarriere_ImpliedVol(double f,double k, double bup, double bdown, double opt, double t,double r,double b,
								   int callput)
{
	struct PricingFunctionToInverse : public DoubleToDoubleFunc 
	{
		double f0;
		double k0;
		double bup0;
		double bdown0;
		double t0;
		double r0;
		double b0;

		int callput0; 
		PricingFunctionToInverse(double fa,double ka, double bupa,double bdowna,  double ta,double ra,double ba,
								   int callputa):
		f0(fa),k0(ka),bup0(bupa),bdown0(bdowna),t0(ta),r0(ra),b0(ba),callput0(callputa)
		{}

		virtual double operator() (double v)  const
		{
			ArgumentList a(f0,k0,bup0,bdown0,v,t0,r0,b0,callput0);
			Power_Expression<ARM_CF_BS_DoubleBarriere_Formula> y;
			return y(a);
		}
	};
	
	PricingFunctionToInverse x(f,k,bup,bdown,t,r,b,callput);
	double typical=0.2;
	double typical_change=0.02;
	double tolerance=1e-12;
	return Inverse(x,Inverse::ALWAYSPOSITIVE)(opt,typical,typical_change,tolerance); 
}

double Export_BS_PartialTime_Start_SingleBarrier(double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bendtime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula> y;
	return y(a);
}

double Export_BS_PartialTime_End_SingleBarrier(double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bstarttime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_End_SingleBarrier_Formula> y;
	return y(a);
}

double Export_BS_SingleBarrier_2Asset(double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype)
{
	ArgumentList a( f1, k1,  f2, k2,  v1,  v2, corr,t, callput,  optype);
	Power_Expression<ARM_CF_BS_SingleBarrier_2Asset_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_BS_EuroBarriere(int i,double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown)
{
	int optiontype;
	if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT;
	else
			if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT;
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BS_EuroBarriere : inout/updown : bad input :");

	ArgumentList a(f,k,b,r,v,t,rate,callput,optiontype);
	Power_Expression<ARM_CF_BS_SingleBarriere_Formula> y;
	return y(i,a);
}
double Export_BS_EuroDoubleBarriere(int i,double f,double k, double bup, double bdown, double v, double t,
								   int callput)
{
	ArgumentList a(f,k,bup,bdown,v,t,callput);
	Power_Expression<ARM_CF_BS_DoubleBarriere_Formula> y;
	return y(i,a);
}

double Export_BS_PartialTime_Start_SingleBarrier(int i,double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bendtime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula> y;
	return y(i,a);
}

double Export_BS_PartialTime_End_SingleBarrier(int i,double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bstarttime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_End_SingleBarrier_Formula> y;
	return y(i,a);
}

double Export_BS_SingleBarrier_2Asset(int i,double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype)
{
	ArgumentList a( f1, k1,  f2, k2,  v1,  v2, corr,t, callput,  optype);
	Power_Expression<ARM_CF_BS_SingleBarrier_2Asset_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_BS_EuroBarriere(int i,int j,double f,double k, double b, double r, double v, double t,double rate,
								   int callput, int inout, int updown){
	int optiontype;
	if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_IN)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN;
	else
		if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::DOWN)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT;
	else
			if((inout==ARM_CF_BS_SingleBarriere_Formula::KNOCK_OUT)&&(updown==ARM_CF_BS_SingleBarriere_Formula::UP)) 
		optiontype=ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT;
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_BS_EuroBarriere : inout/updown : bad input :");

	ArgumentList a(f,k,b,r,v,t,rate,callput,optiontype);
	Power_Expression<ARM_CF_BS_SingleBarriere_Formula> y;
	return y(i,j,a);
}

double Export_BS_EuroDoubleBarriere(int i,int j,double f,double k, double bup, double bdown, double v, double t,
								   int callput)
{
	ArgumentList a(f,k,bup,bdown,v,t,callput);
	Power_Expression<ARM_CF_BS_DoubleBarriere_Formula> y;
	return y(i,j,a);
}


double Export_BS_PartialTime_Start_SingleBarrier(int i,int j,double f,double k, double barrier, double rebate, double v,
												 double bendtime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bendtime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula> y;
	return y(i,j,a);
}

double Export_BS_PartialTime_End_SingleBarrier(int i,int j,double f,double k, double barrier, double rebate, double v,
												 double bstarttime, double t,int callput, int optype)
{
	ArgumentList a(f,k,barrier,rebate,v,bstarttime,t,callput,optype);
	Power_Expression<ARM_CF_BS_PartialTime_End_SingleBarrier_Formula> y;
	return y(i,j,a);
}

double Export_BS_SingleBarrier_2Asset(int i,int j,double f1,double k1, double f2,double k2, double v1, double v2,double corr,
												double t,int callput, int optype)
{
	ArgumentList a( f1, k1,  f2, k2,  v1,  v2, corr,t, callput,  optype);
	Power_Expression<ARM_CF_BS_SingleBarrier_2Asset_Formula> y;
	return y(i,j,a);
}

CC_END_NAMESPACE()
 
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

