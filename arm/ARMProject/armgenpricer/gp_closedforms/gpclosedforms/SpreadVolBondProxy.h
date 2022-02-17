
#ifndef _GP_CF_SPREADVBPROXY_H
#define _GP_CF_SPREADVBPROXY_H

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
#include "gpclosedforms/spreadoption_lognormal_interface.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_SpreadVBProxy : public ARM_RootObject
{
private:
	double				AsOfDate;
	std::vector<double>		itsResetDates;
	std::vector<double>		itsFwdRates1;
	std::vector<double>		itsFwdRatesVol1;
	std::vector<double>		itsFwdRatesFwdVol1;
	std::vector<double>		itsVols1;
	std::vector<double>		itsNu1;
	std::vector<double>		itsRho1;
	std::vector<double>		itsFwdRates2;
	std::vector<double>		itsFwdRatesVol2;
	std::vector<double>		itsFwdRatesFwdVol2;
	std::vector<double>		itsVols2;
	std::vector<double>		itsNu2;
	std::vector<double>		itsRho2;
	std::vector<double>		itsCorrel12;
	std::vector<double>		itsfwdLeverage;
	std::vector<double>		itsfwdStrikes;

	int					itsNbSimul;

	std::vector<double>		itsStdVB1;
	std::vector<double>		itsStdVB2;
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

	double				itsCorrMeanRev;
	double				itsCorrVol;

	std::vector<double>		itsLevier;
	std::vector<double>		itsFix;

public:

	ARM_SpreadVBProxy(double AsOf, const std::vector<double>& resetDates, const std::vector<double>& fwdRates1, const std::vector<double>& fwdRates2, 
		const std::vector<double>& fwdratevols1, const std::vector<double>& fwdratefwdvols1, const std::vector<double>& vols1,
		const std::vector<double>& nu1, const std::vector<double>& rho1,
		const std::vector<double>& fwdratevols2, const std::vector<double>& fwdratefwdvols2, const std::vector<double>& vols2,
		const std::vector<double>& nu2, const std::vector<double>& rho2,
		const std::vector<double>& Correl12,
		int NbSimul, const std::vector<double>& fwdLev, const std::vector<double>& fwdK, int priceopt, int type1sens, int type2sens, bool sabrdiff,
		double corrmeanrev, double corrvol,
		const std::vector<double>& levier, const std::vector<double>& fix)
	{
		AsOfDate			= AsOf;
		itsResetDates		= resetDates;
		itsFwdRates1		= fwdRates1;
		itsFwdRates2		= fwdRates2;
		
		itsFwdRatesVol1		= fwdratevols1;
		itsFwdRatesFwdVol1	= fwdratefwdvols1;
		itsVols1			= vols1;
		itsNu1				= nu1;
		itsRho1				= rho1;
		
		itsFwdRatesVol2		= fwdratevols2;
		itsFwdRatesFwdVol2	= fwdratefwdvols2;
		itsVols2			= vols2;
		itsNu2				= nu2;
		itsRho2				= rho2;

		itsCorrel12			= Correl12;

		itsNbSimul			= NbSimul;
		itsSABRDiff			= sabrdiff;
		itsfwdLeverage		= fwdLev;
		itsfwdStrikes		= fwdK;
		itsPriceOpt			= priceopt;
		itsType1sens		= type1sens;
		itsType2sens		= type2sens;

		itsCorrMeanRev		= corrmeanrev;
		itsCorrVol			= corrvol;

		itsLevier			= levier;
		itsFix				= fix;

		itsPriceType1		= false;
		itsPriceType2		= false;

		switch(priceopt)
		{
		case 1:
			itsPriceType1 = true;
			break;

		case 2:
			itsPriceType2 = true;
			break;

		case 12:
			itsPriceType1 = itsPriceType2 = true;
			break;
		}

		if(itsPriceType1) itsStdVB1.resize(itsResetDates.size()-1, 0.);
		if(itsPriceType2) itsStdVB2.resize(itsResetDates.size()-1, 0.);
		itsCoupon.resize(itsResetDates.size()-1, 0.);

		price();
	}

	ARM_SpreadVBProxy(const ARM_SpreadVBProxy& rhs) : 
		itsStdVB1(rhs.itsStdVB1), itsStdVB2(rhs.itsStdVB2), itsCoupon(rhs.itsCoupon)
	{
	}

	~ARM_SpreadVBProxy() {};

	virtual ARM_Object*		Clone() const { return new ARM_SpreadVBProxy(*this); };
	virtual string			ExportShortName() const;
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	double	GetStdVB1(int i) const {return itsStdVB1[i];};
	double	GetStdVB2(int i) const {return itsStdVB2[i];};
	double	GetCoupon(int i) const {return itsCoupon[i];};

	const std::vector<double>&	GetStdVB1() const {return itsStdVB1;};
	const std::vector<double>&	GetStdVB2() const {return itsStdVB2;};
	const std::vector<double>&	GetCoupon() const {return itsCoupon;};

	void	price();

	void	sabrFwdRate(std::vector<double>& fwdrate, const std::vector<double>& u, double fwd, double t, double alpha, double rho, double nu, int nbPoints);
};

inline string ARM_SpreadVBProxy::ExportShortName() const
{
	return "LPROX";
}

CC_END_NAMESPACE()

#endif
