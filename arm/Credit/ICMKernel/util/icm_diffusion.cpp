#include "icm_diffusion.h"


double TauxVasicek(double t,
				   double T,
				   double a,
				   double s,
				   double IRVol,
				   vector<double> ZCCurve)
{

	// ERROR
	if ((fabs(a)<1e-10) || (fabs(T-t)<1e-10))
	{
		ICMTHROW(ERR_FUNCTION_EVAL_OVFL,"Division by 0."); 
	}

	double TxVasicek = 0.;

	double Beta = (1. - exp(-a * (T - t))) / a;

	double D1 = 1./(2*a) * (exp(-2.*a*(T-t)) - exp(-2.*a*T));
	double D2 = (-1./(2.*a))*(1-exp(-2.*a*t));
	double D3 = (-2./a)*(exp(-a*(T-t)) - exp(-a*T));
	double D4 = (2./a)*(1-exp(-a*t));
	double Delta =  D1 + D2 + D3 + D4;

	double B = (ZC(T,ZCCurve)/ZC(t,ZCCurve)) * exp(-Beta * s) * exp(-0.5*IRVol*IRVol*Delta/(a*a));

	TxVasicek = -log(B)/(T-t);

	return TxVasicek;
}

double ZC(double time, vector<double> ZCCurve)
{
	double value = 0.;
	int indice = int(4*time);
	
	if (indice>=  ZCCurve.size())
	{
		ICMTHROW(ERR_INDEX_OUT_OF_RANGE,"INDEX OUT OF RANGE"); 
	}

	value = ZCCurve[indice] + (4*time - indice) * (ZCCurve[indice + 1] - ZCCurve[indice]);

	return value;
}


double Forward(double T,
			   double dt,
			   vector<double> ZCCurve)
{
//	double Fwd = (-1/ZC(T,ZCCurve)) * (ZC(T,ZCCurve) - ZC(T-dt,ZCCurve))/dt;
	double Fwd = -(1/ZC(T,ZCCurve)) * (ZC(T+dt,ZCCurve) - ZC(T,ZCCurve))/dt;
	return Fwd;
}

double Diff(double Vinit,
			double Mean,
			double SpeedReversion,
			double Vol,
			double dt,
			double Norm01)
{
	double result = 0.;
	result = Vinit*(1 - SpeedReversion*dt) + SpeedReversion*Mean*dt + Vol*sqrt(dt)*Norm01;
	return result;
}
