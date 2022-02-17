
#ifndef _GP_CF_VBPROXY_H
#define _GP_CF_VBPROXY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gpvector.h"
#include "gpbase/vectormanip.h"
#include "gpbase/typedef.h"
#include "gpbase/numericconstant.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpnumlib/solver.h"
#include "gpclosedforms/inverse.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_VBMinMaxProxy : public ARM_RootObject
{
private:
	double				AsOfDate;
	std::vector<double>		itsResetDates;
	std::vector<double>		itsFwdRates;
	std::vector<double>		itsTotalVol;
	std::vector<double>		itsLeftVol;
	std::vector<double>		itsRightVol;
	std::vector<double>		itsNu;
	std::vector<double>		itsRho;

	std::vector<double>		itsStartLev;
	std::vector<double>		itsStartAdd;

	std::vector<double>		itsCapStartLev;
	std::vector<double>		itsCapStartAdd;
	std::vector<double>		itsVBLev;

	int					itsMaxChoice;
	int					itsMinChoice;

	int					itsNbSimul;

	std::vector<double>		itsMaxRate;
	std::vector<double>		itsMinRate;
	std::vector<double>		itsStdVB1;
	std::vector<double>		itsStdVB2;
	std::vector<double>		itsStdVB3;
	std::vector<double>		itsCoupon;

	int					itsType1sens;
	int					itsType2sens;

	double				itsPrice;

	bool				itsSABRDiff;
	
	int					itsPriceOpt;	// = 1 vol bond type 1 : Max(S(T) - S(t) - K, 0)
										// = 2 vol bond type 2 : Max(S(T) - S(T vu de t) - K, 0)
										// = 3 vol bond type 3 : Max(S) sur [t,T] - Min(S) sur [t,T]
										// = 12 type 1 + type 2
										// = 13 type 1 + type 3
										// = 23 type 2 + type 3
										// = 123 type 1 + type 2 + type 3

	bool				itsPriceType1;
	bool				itsPriceType2;
	bool				itsPriceType3;

	double				itsMinMaxFreq;

public:

	ARM_VBMinMaxProxy(double AsOf, const std::vector<double>* resetDates, const std::vector<double>* fwdRates, 
		const std::vector<double>* totalvol, const std::vector<double>* leftvol, const std::vector<double>* rightvol,
		const std::vector<double>* nu, const std::vector<double>* rho, int nbSimul, bool sabrdiff,
		int priceopt, int type1sens, int type2sens,
		const std::vector<double>* rate1Leverage, const std::vector<double>* rate1Add,
		const std::vector<double>* capLeverage, const std::vector<double>* capAdd, const std::vector<double>* vbLeverage,
		int maxchoice, int minchoice, double freq)
	{
		AsOfDate			= AsOf;
		itsResetDates		= resetDates;
		itsFwdRates			= fwdRates;
		itsTotalVol			= totalvol;
		itsLeftVol			= leftvol;
		itsRightVol			= rightvol;
		itsNu				= nu;
		itsRho				= rho;
		itsNbSimul			= nbSimul;
		itsSABRDiff			= sabrdiff;

		itsStartLev			= rate1Leverage;
		itsStartAdd			= rate1Add;
		itsCapStartLev		= capLeverage;
		itsCapStartAdd		= capAdd;
		itsVBLev			= vbLeverage;

		itsMaxChoice		= maxchoice;
		itsMinChoice		= minchoice;
		itsMinMaxFreq		= freq;

		itsPriceOpt			= priceopt;
		itsType1sens		= type1sens;
		itsType2sens		= type2sens;

		itsPriceType1		= false;
		itsPriceType2		= false;
		itsPriceType3		= false;

		switch(priceopt)
		{
		case 1:
			itsPriceType1 = true;
			break;

		case 2:
			itsPriceType2 = true;
			break;

		case 3:
			itsPriceType3 = true;
			break;

		case 12:
			itsPriceType1 = itsPriceType2 = true;
			break;

		case 13:
			itsPriceType1 = itsPriceType3 = true;
			break;

		case 23:
			itsPriceType2 = itsPriceType3 = true;
			break;

		case 123:
			itsPriceType1 = itsPriceType2 = itsPriceType3 = true;
			break;
		}

		if(itsPriceType1) itsStdVB1.resize(itsResetDates.size()-1, 0.);
		if(itsPriceType2) itsStdVB2.resize(itsResetDates.size()-1, 0.);
		itsCoupon.resize(itsResetDates.size()-1,0.);

		if(itsPriceType3)
		{
			itsStdVB3.resize(itsResetDates.size()-1, 0.);
			itsMaxRate.resize(itsResetDates.size()-1, 0.);
			itsMinRate.resize(itsResetDates.size()-1, 0.);
		}

		price();
	}

	ARM_VBMinMaxProxy(const ARM_VBMinMaxProxy& rhs) : itsMaxRate(rhs.itsMaxRate), itsMinRate(rhs.itsMinRate), 
		itsStdVB1(rhs.itsStdVB1), itsStdVB2(rhs.itsStdVB2), itsStdVB3(rhs.itsStdVB3), itsCoupon(rhs.itsCoupon)
	{
	}

	~ARM_VBMinMaxProxy() {};

	virtual ARM_Object*		Clone() const { return new ARM_VBMinMaxProxy(*this); };
	virtual string			ExportShortName() const;
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	double	GetMaxRate(int i) const {return itsMaxRate[i];};
	double	GetMinRate(int i) const {return itsMinRate[i];};
	double	GetStdVB1(int i) const {return itsStdVB1[i];};
	double	GetStdVB2(int i) const {return itsStdVB2[i];};
	double	GetStdVB3(int i) const {return itsStdVB3[i];};
	double	GetCoupon(int i) const {return itsCoupon[i];};

	const std::vector<double>*	GetMaxRate() const {return itsMaxRate;};
	const std::vector<double>*	GetMinRate() const {return itsMinRate;};
	const std::vector<double>*	GetStdVB1() const {return itsStdVB1;};
	const std::vector<double>*	GetStdVB2() const {return itsStdVB2;};
	const std::vector<double>*	GetStdVB3() const {return itsStdVB3;};
	const std::vector<double>*	GetCoupon() const {return itsCoupon;};

	double	E_Max(double fwd, double var);
	double	E_Min(double fwd, double var);
	
	double	repartmax(double fwd, double stddev, double x);
	double	repartmin(double fwd, double stddev, double x);
	double	densmax(double fwd, double stddev, double x);
	double	densmin(double fwd, double stddev, double x);

	double	Integrale(double fwd, double stddev, double lbound, double ubound,int sens);

	void	price();

	void	sabrFwdRate(std::vector<double>* fwdrate, const std::vector<double>* u, double fwd, double t, double alpha, double rho, double nu, int nbPoints);
};

inline double ARM_VBMinMaxProxy::repartmax(double fwd, double stddev, double x)
{
	double d0 = log(fwd/x)/stddev - 0.5*stddev;
	double d1 = log(fwd/x)/stddev + 0.5*stddev;

	double r = (fwd/x) * NormalCDF(d1) + NormalCDF(d0);

	return 1. - r;
}

inline double ARM_VBMinMaxProxy::repartmin(double fwd, double stddev, double x)
{
	double d0 = log(fwd/x)/stddev - 0.5*stddev;
	double d1 = log(fwd/x)/stddev + 0.5*stddev;

	double r = NormalCDF(-d0) + (fwd/x) * NormalCDF(-d1);

	return r;
}

inline double ARM_VBMinMaxProxy::Integrale(double fwd, double stddev, double lbound, double ubound, int sens)
{
	/*
	double (ARM_VBMinMaxProxy::*repart)(double,double,double);
	sens == 1 ? repart = ARM_VBMinMaxProxy::repartmax : repart = ARM_VBMinMaxProxy::repartmin;

	GaussLegendre_Coefficients c(itsNbLegPts);

	double scale = ubound - lbound, sum = 0., x;

	for(int i = 0; i < itsNbLegPts; i++)
	{
		x = lbound + scale * (0.5 + 0.5 * c.get_point(i));
		sum += (1. - (this->*repart)(fwd,stddev,x)) * 0.5 * c.get_weight(i);
	}
	
	sum *= scale;

	return sum;
	*/

	return 0.;
}

inline double ARM_VBMinMaxProxy::E_Max(double fwd, double var)
{
	/*
	double ubound, lbound = fwd;
	double integ = 0.;

	for(int k = 0; k < 5; k++)
	{
		ubound = fwd*exp((k+1)*sqrt(var));
		integ += Integrale(fwd, sqrt(var), lbound, ubound, 1);
		lbound = ubound;
	}

	if(ubound < fwd*10.) integ += Integrale(fwd,sqrt(var),ubound,fwd*10.,1);

	return integ + fwd;
	*/

	double stddev = sqrt(var);
	double N = NormalCDF(0.5*stddev);
	double n = exp(-0.5*(0.5*stddev)*(0.5*stddev)) / ARM_NumericConstants::ARM_SQRT_2_PI;

	double r = 2. * fwd * N + stddev * fwd * (0.5 * stddev * N + n);

	return r;
}

inline double ARM_VBMinMaxProxy::E_Min(double fwd, double var)
{
	/*
	double ubound = fwd, lbound = fwd;
	double integ = 0.;

	for(int k = 0; k < 5; k++)
	{
		lbound = fwd*exp(-(k+1)*sqrt(var));
		integ += Integrale(fwd, sqrt(var), lbound, ubound,-1);
		ubound = lbound;
	}

	if(lbound > 1e-8) integ += Integrale(fwd, sqrt(var), 1e-8, lbound, -1);

	return integ;
	*/

	double stddev = sqrt(var);
	double N = NormalCDF(-0.5*stddev);
	double n = exp(-0.5*(0.5*stddev)*(0.5*stddev)) / ARM_NumericConstants::ARM_SQRT_2_PI;

	double r = 2. * fwd * N - stddev * fwd * (-0.5 * stddev * N + n);

	return r;
}

inline string ARM_VBMinMaxProxy::ExportShortName() const
{
	return "LPROX";
}

CC_END_NAMESPACE()

#endif
