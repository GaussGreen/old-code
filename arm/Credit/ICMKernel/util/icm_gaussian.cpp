
#include "ICMKernel/glob/icm_maths.h" 
# include <math.h>
#include <vector>

// #include "ournag.h"
#include <nags.h>

#include "icm_gaussian.h"
using namespace std;





// ----------------------------------------------------------------------------------
//	class		ICM_gaussian
//  author		L. Jacquel
//	version		1.0
//	date		September 2004
//	file		ICM_gaussian.cpp

//	brief		Useful functions for Standard Gaussian handling
// ----------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------
double NCum(double	x)
{
	double k;
	if (x > 0.0)
	{
		if (x < 15)
		{
			k = 1.0 / (1.0 + NORCUMA * x);
			return 1.0 - exp(- 0.5 * x * x ) * ((NORCUMB1 + (NORCUMB2 + (NORCUMB3 + (NORCUMB4 + NORCUMB5 * k) * k) * k) * k) * k) * ONEOVERSQRT2PI;
		}
		else
			return 1.0;
	}
	else if (x < 0.0)
	{
		if (-x < 15)
		{
			k = 1.0 / (1.0 - NORCUMA * x);
			return exp(- 0.5 * x * x ) * ((NORCUMB1 + (NORCUMB2 + (NORCUMB3 + (NORCUMB4 + NORCUMB5 * k) * k) * k) * k) * k) * ONEOVERSQRT2PI;
		}
		else
			return 0.0;
	}
	else 
		return 0.5;
}
// ----------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------
double InvNCum(double InValue)
{
	double X, S, R;
                        
	X = InValue - 0.5;
	if (fabs(X) < 0.42)
	{
		R = X * X;
		S = 1.0 + (INVNORCUMB1 + (INVNORCUMB2 + (INVNORCUMB3 + INVNORCUMB4 * R) * R) * R) * R;
		if (fabs(S) < INVNORCUMEPS) 
			return INVNORCUMDEF;
		return X * (INVNORCUMA0 + (INVNORCUMA1 + (INVNORCUMA2 + INVNORCUMA3 * R) * R) * R) / S;
	}

	if (fabs(X) >= 0.5) return INVNORCUMDEF;
	R = (X > 0.0) ? 1.0 - InValue : InValue;
	R = log(-log(R));
	R = INVNORCUMC0 + (INVNORCUMC1 + (INVNORCUMC2 + (INVNORCUMC3 + (INVNORCUMC4 + (INVNORCUMC5 + (INVNORCUMC6
	  + (INVNORCUMC7 + INVNORCUMC8 * R) * R) * R) * R) * R) * R) * R) * R;
	return (X < 0.0) ? -R : R;
}
// ----------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------
double StdNorDensity(double x)
{
	if (fabs(x) > 15)
		return 0.0;
	else
		return exp(- 0.5 * x * x) * ONEOVERSQRT2PI;
}
// ----------------------------------------------------------------------------------


// ----------------------------------------------------------------------------------
double ExpMinusHalfXSquare(double x)
{
	if (fabs(x) > 15)
		return 0.0;
	else
		return exp(- 0.5 * x * x);	//	 * ONEOVERSQRT2PI;
}
// ----------------------------------------------------------------------------------
//----------------Binormale cumulée--------------------------------------------------
double IntermFunction(double a,
					  double b,
					  double rho,
					  double x,
				      double y)

{
	double res = exp(a*(2*x-a)+b*(2*y-b)+2*rho*(x-a)*(y-b));
	return res;
}

double BinormalGenericCase(double a,
						   double b,
						   double rho)

{

	vector<double> A;
	A.resize(4);
	vector<double> B;
	B.resize(4);

	
	double res = 0.;
	double tmpij = 0.;
	int i = 0;
	int j = 0;

	A[0]=0.3253030;
	A[1]=0.4211071;
	A[2]=0.1334425;
	A[3]=0.006374323;

	B[0]=0.1337764;
	B[1]=0.6243247;
	B[2]=1.3425378;
	B[3]=2.2626645;

double	_a=a/sqrt(2*(1-rho*rho));
double	_b=b/sqrt(2*(1-rho*rho));

			for(i=0;i<4;i++)
			{
				for(j=0;j<4;j++)
				{
					tmpij = IntermFunction(_a,_b,rho,B[i],B[j]);
					res  += A[i]*A[j]*tmpij;
				}
			}

			res *= sqrt(1-rho*rho)/(3.14159265358979);
	
	return res;
}

double BinormalCumuleNegativeorzerocase(double a,
										double b,
										double rho)
{

double res=0.0;
/*
double	Na = ep::MathSrv::cumNorm(a);
double	Nb = ep::MathSrv::cumNorm(b);
*/
double Na = NAG_cumul_normal(a);
double Nb = NAG_cumul_normal(b);

			if ((a<=0.0) && (b<=0.0) && (rho<=0.0))
				res= BinormalGenericCase(a,b,rho);

			if ((a<=0.0) && (b>=0.0) && (rho>=0.0))
				res=Na-BinormalGenericCase(a,-b,-rho);

			if ((a>=0.0) && (b<=0.0) && (rho>=0.0))
				res=Nb-BinormalGenericCase(-a,b,-rho);

			if ((a>=0.0) && (b>=0.0) && (rho<=0.0))
				res=Na+Nb-1+BinormalGenericCase(-a,-b,rho);
	


return res;

}


double BinormalCumule(double a,
					  double b,
					  double rho)
{
double res=0.0,test=a*b*rho;
double test1,test2;
double rho1=(rho*a-b)*SIGN(1,a)/sqrt(a*a-2*rho*a*b+b*b);
double rho2=(rho*b-a)*SIGN(1,b)/sqrt(a*a-2*rho*a*b+b*b);
double delta=(1-SIGN(1,a)*SIGN(1,b))/4.0;

		if (test<=0.0)
		{
			res=BinormalCumuleNegativeorzerocase(a,b,rho);
		}
		else
		{
			test1=BinormalCumuleNegativeorzerocase(a,0,rho1);
			test2=BinormalCumuleNegativeorzerocase(b,0,rho2);
			res=test1+test2-delta;
		}
return res;
}
	
