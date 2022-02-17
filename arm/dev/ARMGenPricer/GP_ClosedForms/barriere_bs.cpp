/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\Barriere_bs.cpp
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
#include "gpclosedforms/bivariate_normal.h"


#include"gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/Barriere_bs.h"
#include "gpclosedforms/Barriere_bs_formula.h"					/// for the flags
#include "gpclosedforms/normal.h"


CC_BEGIN_NAMESPACE(ARM)

//////////////////////////////////////////////////////////////////////////////////////////
///
///   Pricing Functions
///
//////////////////////////////////////////////////////////////////////////////////////////
double BS_SingleBarrierOption(double S,		/// underlying
							  double X,		/// strike
							  double H,		/// barrier
							  double K,     /// rebate
							  double T,     /// expiry
							  double sig,   /// volatility
							  double r,     /// rate risk free ( discount rate)
							  double b,     /// dividend
							  int callput,
							  int optiontype)
{
	double phi;
	double eta;
	switch (optiontype)
	{
	case ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN :
		{
			if(H<=0) return 0.0;
			if(S<=H) return BlackSholes_Formula(S,sig,exp(-r*T),X,T,callput);
			eta=1.;
			if (callput==K_CALL) phi=1.; else phi=-1.;
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN :
		{
			if(H<=0) return BlackSholes_Formula(S,sig,exp(-r*T),X,T,callput);
			if(S>=H) return BlackSholes_Formula(S,sig,exp(-r*T),X,T,callput);
			eta=-1.;
			if (callput==K_CALL) phi=1.; else phi=-1.;
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT :
		{
			if(H<=0) return BlackSholes_Formula(S,sig,exp(-r*T),X,T,callput);
			if(S<=H) return 0;
			eta=1.;
			if (callput==K_CALL) phi=1.; else phi=-1.;
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT :
		{
			if(H<=0) return 0.0;
			if(S>=H) return 0;
			eta=-1.;
			if (callput==K_CALL) phi=1.; else phi=-1.;
			break;
		}
	}

	double totalvol=sig*sqrt(T);
	double Sb=S*exp((b-r)*T);
	double Xb=X*exp(-r*T);
	double mu=(b-0.5*sig*sig)/sig/sig;
	double lambda=sqrt(mu*mu+2.*r/sig/sig);
	double z=log(H/S)/totalvol+lambda*totalvol;
	double x1=log(S/X)/totalvol+(1.+mu)*totalvol;
	double x2=log(S/H)/totalvol+(1.+mu)*totalvol;
	double y1=log(H*H/X/S)/totalvol+(1.+mu)*totalvol;
	double y2=log(H/S)/totalvol+(1.+mu)*totalvol;
	double A=phi*Sb*NormalCDF(phi*x1)-phi*Xb*NormalCDF(phi*x1-phi*sig*sqrt(T));
	double B=phi*Sb*NormalCDF(phi*x2)-phi*Xb*NormalCDF(phi*x2-phi*sig*sqrt(T));
	double C=phi*Sb*pow(H/S,2.*(mu+1.))*NormalCDF(eta*y1)-phi*Xb*pow(H/S,2.*mu)*NormalCDF(eta*y1-eta*sig*sqrt(T));
	double D=phi*Sb*pow(H/S,2.*(mu+1.))*NormalCDF(eta*y2)-phi*Xb*pow(H/S,2.*mu)*NormalCDF(eta*y2-eta*sig*sqrt(T));
	double E=K*exp(-r*T)*(NormalCDF(eta*x2-eta*sig*sqrt(T))-pow(H/S,2.*mu)*NormalCDF(eta*y2-eta*sig*sqrt(T)));
	double F=K*(pow(H/S,mu+lambda)*NormalCDF(eta*z)+pow(H/S,mu-lambda)*NormalCDF(eta*z-2.*eta*lambda*sig*sqrt(T)));

	switch (optiontype)
	{
	case ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_IN :
		{
			if (callput==K_CALL)
			{
				if(X>H) return C+E; else return A-B+D+E;
			}else{
				if(X>H) return B-C+D+E; else return A+E;
			}
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::UP_AND_IN :
		{
			if (callput==K_CALL)
			{
				if(X>H) return A+E; else return B-C+D+E;
			}else{
				if(X>H) return A-B+D+E; else return C+E;
			}
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::DOWN_AND_OUT :
		{
			if (callput==K_CALL)
			{
				if(X>H) return A-C+F; else return B-D+F;
			}else{
				if(X>H) return A-B+C-D+F; else return F;
			}
			break;
		}
	case ARM_CF_BS_SingleBarriere_Formula::UP_AND_OUT :
		{
			if (callput==K_CALL)
			{
				if(X>H) return F; else return A-B+C-D+F;
			}else{
				if(X>H) return B-D+F; else return A-C+F;
			}
			
			break;
		}
	}
	return 0;
}

double BS_DoubleBarrierCall(double S,
							double X, 
							double L,
							double U,
							double T,
							double sig,
							double r,
							double b,
							double delta1,
							double delta2,
							int nbterms)
{
	double totalvol=sig*sqrt(T);
	double F=U*exp(delta1*T);
	double mu1,mu2,mu3,d1,d2, d3,d4;
	double sum1=0.,sum2=0.;
	for(int n=-nbterms;n<=nbterms;n++)
	{
		mu1=2*(b-delta2-n*(delta1-delta2))/(sig*sig)+1;
		mu2=2*n*(delta1-delta2)/(sig*sig);
		mu3=2*(b-delta2+n*(delta1-delta2))/(sig*sig)+1;
		d1=(log(S*pow(U/L,2.*n)/X)+(b*b+sig*sig/2.)*T)/totalvol;
		d2=(log(S*pow(U/L,2.*n)/F)+(b*b+sig*sig/2.)*T)/totalvol;
		d3=(log(L*L*pow(L/U,2.*n)/X/S)+(b*b+sig*sig/2.)*T)/totalvol;
		d4=(log(L*L*pow(L/U,2.*n)/F/S)+(b*b+sig*sig/2.)*T)/totalvol;
		sum1+=pow(U/L,n*mu1)*pow(L/S,mu2)*(NormalCDF(d1)-NormalCDF(d2))-pow(pow(L/S,n)*L/S,mu3)*(NormalCDF(d3)-NormalCDF(d4));
		sum2+=pow(U/L,n*(mu1-2.))*pow(L/S,mu2)*(NormalCDF(d1-totalvol)-NormalCDF(d2-totalvol))-pow(pow(L/S,n)*L/S,mu3-2.)*(NormalCDF(d3-totalvol)-NormalCDF(d4-totalvol));
	}
	return S*exp((b-r)*T)*sum1-X*exp(-r*T)*sum2;
}

double BS_DoubleBarrierPut(double S,
						   double X,
						   double L,
						   double U,
						   double T, 
						   double sig,
						   double r,
						   double b,
						   double delta1,
						   double delta2,
						   int nbterms)
{
	double totalvol=sig*sqrt(T);
	double E=L*exp(delta2*T);
	double mu1,mu2,mu3,y1,y2, y3,y4;
	double sum1=0.,sum2=0.;
	for(int n=-nbterms;n<=nbterms;n++)
	{
		mu1=2*(b-delta2-n*(delta1-delta2))/(sig*sig)+1;
		mu2=2*n*(delta1-delta2)/(sig*sig);
		mu3=2*(b-delta2+n*(delta1-delta2))/(sig*sig)+1;
		y1=(log(S*pow(U/L,2.*n)/E)+(b*b+sig*sig/2.)*T)/totalvol;
		y2=(log(S*pow(U/L,2.*n)/X)+(b*b+sig*sig/2.)*T)/totalvol;
		y3=(log(L*L*pow(L/U,2.*n)/E/S)+(b*b+sig*sig/2.)*T)/totalvol;
		y4=(log(L*L*pow(L/U,2.*n)/X/S)+(b*b+sig*sig/2.)*T)/totalvol;
		sum1+=pow(U/L,n*(mu1-2.))*pow(L/S,mu2)*(NormalCDF(y1-totalvol)-NormalCDF(y2-totalvol))-pow(pow(L/U,n)*L/S,mu3-2.)*(NormalCDF(y3-totalvol)-NormalCDF(y4-totalvol));
		sum2+=pow(U/L,n*mu1)*pow(L/S,mu2)*(NormalCDF(y1)-NormalCDF(y2))-pow(pow(L/U,n)*L/S,mu3)*(NormalCDF(y3)-NormalCDF(y4));
	}
	return X*exp(-r*T)*sum1-S*exp((b-r)*T)*sum2;
}

double BS_DoubleBarrierOption(double S,
							  double X, 
							  double L, 
							  double U,
							  double T, 
							  double sig,
							  double r,
							  double b,
							  double delta1,
							  double delta2,
							  int callput,
							  int nbterms)
{
	switch (callput)
	{
	case K_CALL :
		{
			return BS_DoubleBarrierCall(S,X,L,U,T,sig,r,b,delta1,delta2,nbterms);
			break;
		}
	case K_PUT :
		{
			return BS_DoubleBarrierPut(S,X,L,U,T,sig,r,b,delta1,delta2,nbterms);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_DoubleBarrierOption : callput , bad input :");
			break;
		}
	}
	
}


double BS_PartialTime_Start_SingleBarrierCall(double S,
											  double X,
											  double H,
											  double K,
											  double t1,
											  double T, 
											  double sig,
											  double r,
											  double b,
											  int optiontype)
{
	double eta;
	switch (optiontype)
	{
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_IN :
		{
			eta=1.;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_IN :
		{
			eta=-1.;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_OUT :
		{
			eta=1.;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_OUT :
		{
			eta=-1.;
			break;
		}
	}
	double totalvol=sig*sqrt(T);
	double totalvol1=sig*sqrt(t1);
	double mu=(b-sig*sig/2.)/sig/sig;
	double rho=sqrt(t1/T);
	double e1=(log(S/H)+(b+sig*sig/2.)*t1)/totalvol1;
	double e2=e1-totalvol1;
	double e3=e1+2.*log(H/S)/totalvol1;
	double e4=e1-totalvol1;
	double f1=(log(S/X)+2.*log(H/S)+(b+sig*sig/2.)*T)/totalvol;
	double f2=f1-totalvol;
	double d1=(log(S/X)+(b+sig*sig/2.)*T)/totalvol;
	double d2=d1-totalvol;
	double ca=S*exp((b-r)*T)*(bivariate_cdfNormal(d1,eta*e1,eta*rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(f1,eta*e3,eta*rho))-
		X*exp(-r*T)*(bivariate_cdfNormal(d2,eta*e2,eta*rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(f2,eta*e4,eta*rho));
	switch (optiontype)
	{
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_IN :
		{
			double f=S*exp((b-r)*T);
			double d1 = (log(f/X))/totalvol+0.5*totalvol ;
			double bond=exp(-r*T);
			double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol));
			return call-ca;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_IN :
		{
			double f=S*exp((b-r)*T);
			double d1 = (log(f/X))/totalvol+0.5*totalvol ;
			double bond=exp(-r*T);
			double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol));
			return call-ca;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::DOWN_AND_OUT :
		{
			return ca;
			break;
		}
	case ARM_CF_BS_PartialTime_Start_SingleBarrier_Formula::UP_AND_OUT :
		{
			return ca;
			break;
		}
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_Start_SingleBarrierCall: should not reach here");
	return 0;
}


double BS_PartialTime_Start_SingleBarrierOption(double S,
												double X,
												double H,
												double K,
												double t1,
												double T, 
												double sig,double r,double b,int callput,int optiontype)
{
	switch (callput)
	{
	case K_CALL :
		{
			return BS_PartialTime_Start_SingleBarrierCall(S,X,H,K,t1,T,sig,r,b,optiontype);
			break;
		}
	case K_PUT :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_Start_SingleBarrierOption : put not implemented :");
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_Start_SingleBarrierOption : callput , bad input :");
			break;
		}
	}
	
}



double BS_PartialTime_End_SingleBarrierCall(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int optiontype)
{
	double totalvol=sig*sqrt(T);
	double totalvol1=sig*sqrt(t1);
	double mu=(b-sig*sig/2.)/sig/sig;
	double rho=sqrt(t1/T);
	double e1=(log(S/H)+(b+sig*sig/2.)*t1)/totalvol1;
	double e2=e1-totalvol1;
	double e3=e1+2.*log(H/S)/totalvol1;
	double e4=e1-totalvol1;
	double f1=(log(S/X)+2.*log(H/S)+(b+sig*sig/2.)*T)/totalvol;
	double f2=f1-totalvol;
	double d1=(log(S/X)+(b+sig*sig/2.)*T)/totalvol;
	double d2=d1-totalvol;
	double g1=log(S/H)+(b+sig*sig/2.)*T/totalvol;
	double g2=g1-totalvol;
	double g3=g1+2.*log(H/S)/totalvol;
	double g4=g3-totalvol;
	double ca;
	if ((optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_OUT) ||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_IN))
	{
		if (X>H)
		{
			ca=S*exp((b-r)*T)*(bivariate_cdfNormal(d1,e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(f1,-e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(d2,e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(f2,-e4,-rho));
		}
		else
		{
			ca=S*exp((b-r)*T)*(bivariate_cdfNormal(-g1,-e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(-g3,e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(-g2,-e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(-g4,e4,-rho))-
				S*exp((b-r)*T)*(bivariate_cdfNormal(-d1,-e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(-f1,e3,-rho)) +
				X*exp(-r*T)*(bivariate_cdfNormal(-d2,-e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(-f2,e4,-rho))+
				S*exp((b-r)*T)*(bivariate_cdfNormal(g1,e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(g3,-e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(g2,e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(g4,-e4,-rho));
		}
	}
	if ((optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_OUT) ||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_IN))
	{
		if (X>H)
		{
			ca=S*exp((b-r)*T)*(bivariate_cdfNormal(d1,e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(f1,-e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(d2,e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(f2,-e4,-rho));
		}
		else
		{
			ca=S*exp((b-r)*T)*(bivariate_cdfNormal(g1,e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(g3,-e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(g2,e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(g4,-e4,-rho));
		}
	}
	if ((optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_OUT) ||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_IN))
	{
		if (X>H)
		{
			ca=0;	
		}
		else
		{
			ca=S*exp((b-r)*T)*(bivariate_cdfNormal(-g1,-e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(-g3,e3,-rho)) -
				X*exp(-r*T)*(bivariate_cdfNormal(-g2,-e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(-g4,e4,-rho))-
				S*exp((b-r)*T)*(bivariate_cdfNormal(-d1,-e1,rho)-pow(H/S,2.*(mu+1.))*bivariate_cdfNormal(e3,-f1,-rho)) +
				X*exp(-r*T)*(bivariate_cdfNormal(-d2,-e2,rho)-pow(H/S,2.*mu)*bivariate_cdfNormal(e4,-f2,-rho));
		}
	}
	if((optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_OUT)||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_OUT) ||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_OUT))
	{
		return ca;
	}
	if((optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_UP_AND_IN)||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::INSIDE_DOWN_AND_IN) ||
		(optiontype==ARM_CF_BS_PartialTime_End_SingleBarrier_Formula::CROSS_AND_IN))
	{
		double f=S*exp((b-r)*T);
		double d1 = (log(f/X))/totalvol+0.5*totalvol ;
		double bond=exp(-r*T);
		double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol));
		return call-ca;
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_End_SingleBarrierCall: should not reach here");
	return 0;
}

double BS_PartialTime_End_SingleBarrierOption(double S,double X, double H,double K,double t1,double T, double sig,double r,double b,int callput,int optiontype)
{
	switch (callput)
	{
	case K_CALL :
		{
			return BS_PartialTime_End_SingleBarrierCall(S,X,H,K,t1,T,sig,r,b,optiontype);
			break;
		}
	case K_PUT :
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_End_SingleBarrierOption : put not implemented :");
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"BS_PartialTime_End_SingleBarrierOption : callput , bad input :");
			break;
		}
	}
	
}


double TwoAsset_Single_Barrier_Option(double S1,double S2,double X,double H,double T,double sig1, 
									  double sig2, double rho,
									double r,double b1,double b2, int callput, int optiontype)
{
	double eta,phi;
	switch (optiontype)
	{
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_IN :
		{
			if (callput==K_CALL)
			{
				eta=1.;
				phi=-1.;
			} else {
				eta=-1.;
				phi=-1.;
			};
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_IN :
		{
			if (callput==K_CALL)
			{
				eta=1.;
				phi=1.;
			} else {
				eta=-1.;
				phi=1.;
			};
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_OUT :
		{
			if (callput==K_CALL)
			{
				eta=1.;
				phi=-1.;
			} else {
				eta=-1.;
				phi=-1.;
			};
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_OUT :
		{
			if (callput==K_CALL)
			{
				eta=1.;
				phi=1.;
			} else {
				eta=-1.;
				phi=1.;
			};
			break;
		}
	}
	double mu1=b1-sig1*sig1/2.;
	double mu2=b2-sig2*sig2/2.;
	double totalvol1=sig1*sqrt(T);
	double totalvol2=sig2*sqrt(T);
	double e1=(log(H/S2)-(mu2+rho*sig1*sig2)*T)/totalvol2;
	double e2=e1+rho*totalvol1;
	double e3=e1-2.*log(H/S2)/totalvol2;
	double e4=e2-2*log(H/S2)/totalvol2;
	double d1=(log(S1/X)+(mu1+sig1*sig1)*T)/totalvol1;
	double d2=d1-totalvol1;
	double d3=d1+2*rho*log(H/S2)/totalvol2;
	double d4=d2+2*rho*log(H/S2)/totalvol2;
	double ca=eta*S1*exp((b1-r)*T)*(bivariate_cdfNormal(eta*d1,phi*e1,-eta*phi*rho)-
		exp((2.*(mu2+rho*sig1*sig2)*log(H/S2))/(sig2*sig2))*bivariate_cdfNormal(eta*d3,phi*e3,-eta*phi*rho))-
		eta*X*exp(-r*T)*(bivariate_cdfNormal(eta*d2,phi*e2,-eta*phi*rho)-
		exp((2.*mu2*log(H/S2))/(sig2*sig2))*bivariate_cdfNormal(eta*d4,phi*e2,-eta*phi*rho));
	switch (optiontype)
	{
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_IN :
		{
			if (callput==K_CALL)
			{
				double f=S1*exp((b1-r)*T);
				double d1 = (log(f/X))/totalvol1+0.5*totalvol1 ;
				double bond=exp(-r*T);
				double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol1));
				return call-ca;
			} else {
				double f=S2*exp((b2-r)*T);
				double d1 = (log(f/X))/totalvol2+0.5*totalvol2 ;
				double bond=exp(-r*T);
				double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol2));
				double put=call-(f-X)*bond;
				return put-ca;
			};
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_IN :
		{
			if (callput==K_CALL)
			{
				double f=S1*exp((b1-r)*T);
				double d1 = (log(f/X))/totalvol1+0.5*totalvol1 ;
				double bond=exp(-r*T);
				double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol1));
				return call-ca;
			} else {
				double f=S2*exp((b2-r)*T);
				double d1 = (log(f/X))/totalvol2+0.5*totalvol2 ;
				double bond=exp(-r*T);
				double call = bond*(f*NormalCDF(d1)-X*NormalCDF(d1-totalvol2));
				double put=call-(f-X)*bond;
				return put-ca;
			};
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::DOWN_AND_OUT :
		{
			return ca;
			break;
		}
	case ARM_CF_BS_SingleBarrier_2Asset_Formula::UP_AND_OUT :
		{
			return ca;
			break;
		}
	}
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"TwoAsset_Single_Barrier_Option: should not reach here");
	return 0;
}


CC_END_NAMESPACE()
 

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

