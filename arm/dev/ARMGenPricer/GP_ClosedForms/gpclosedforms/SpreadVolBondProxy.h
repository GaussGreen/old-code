
#ifndef _GP_CF_SPREADVBPROXY_H
#define _GP_CF_SPREADVBPROXY_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

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
	ARM_GP_Vector		itsResetDates;
	ARM_GP_Vector		itsFwdRates1;
	ARM_GP_Vector		itsFwdRatesVol1;
	ARM_GP_Vector		itsFwdRatesFwdVol1;
	ARM_GP_Vector		itsVols1;
	ARM_GP_Vector		itsNu1;
	ARM_GP_Vector		itsRho1;
	ARM_GP_Vector		itsFwdRates2;
	ARM_GP_Vector		itsFwdRatesVol2;
	ARM_GP_Vector		itsFwdRatesFwdVol2;
	ARM_GP_Vector		itsVols2;
	ARM_GP_Vector		itsNu2;
	ARM_GP_Vector		itsRho2;
	ARM_GP_Vector		itsCorrel12;
	ARM_GP_Vector		itsfwdLeverage;
	ARM_GP_Vector		itsfwdStrikes;

	int					itsNbSimul;

	ARM_GP_Vector		itsStdVB1;
	ARM_GP_Vector		itsStdVB2;
	ARM_GP_Vector		itsCoupon;

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

	ARM_GP_Vector		itsLevier;
	ARM_GP_Vector		itsFix;

public:

	ARM_SpreadVBProxy(double AsOf, const ARM_GP_Vector& resetDates, const ARM_GP_Vector& fwdRates1, const ARM_GP_Vector& fwdRates2, 
		const ARM_GP_Vector& fwdratevols1, const ARM_GP_Vector& fwdratefwdvols1, const ARM_GP_Vector& vols1,
		const ARM_GP_Vector& nu1, const ARM_GP_Vector& rho1,
		const ARM_GP_Vector& fwdratevols2, const ARM_GP_Vector& fwdratefwdvols2, const ARM_GP_Vector& vols2,
		const ARM_GP_Vector& nu2, const ARM_GP_Vector& rho2,
		const ARM_GP_Vector& Correl12,
		int NbSimul, const ARM_GP_Vector& fwdLev, const ARM_GP_Vector& fwdK, int priceopt, int type1sens, int type2sens, bool sabrdiff,
		double corrmeanrev, double corrvol,
		const ARM_GP_Vector& levier, const ARM_GP_Vector& fix)
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

	const ARM_GP_Vector&	GetStdVB1() const {return itsStdVB1;};
	const ARM_GP_Vector&	GetStdVB2() const {return itsStdVB2;};
	const ARM_GP_Vector&	GetCoupon() const {return itsCoupon;};

	void	price();

	void	sabrFwdRate(ARM_GP_Vector& fwdrate, const ARM_GP_Vector& u, double fwd, double t, double alpha, double rho, double nu, int nbPoints);
};

inline string ARM_SpreadVBProxy::ExportShortName() const
{
	return "LPROX";
}

CC_END_NAMESPACE()

#endif
