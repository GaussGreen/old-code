
#ifndef _GP_CF_BERM2DATESPROXY_H
#define _GP_CF_BERM2DATESPROXY_H

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

CC_BEGIN_NAMESPACE( ARM )

class ARM_Berm2DatesProxy : public ARM_RootObject
{
private:
	double				AsOfDate;
	ARM_GP_Vector		itsResetDates;
	ARM_GP_Vector		itsFwdRates;
	ARM_GP_Vector		itsFwdRatesVol;
	ARM_GP_Vector		itsFwdRatesFwdVol;
	ARM_GP_Vector		itsVols;
	ARM_GP_Vector		itsNu;
	ARM_GP_Vector		itsRho;
	double		 		itsRateCorrel;
	double				itsStrike;
	ARM_GP_Vector		itsAnnuity;

	int					itsNbSimul;

	int					itsDiffType;	// 1 = log normal
										// 2 = gaussien
										// 3 = vol sto

	double				itsBerm;
	ARM_GP_Vector		itsEuropean;

public:

	ARM_Berm2DatesProxy(double AsOf, const ARM_GP_Vector& resetDates, const ARM_GP_Vector& fwdRates, double strike,
		const ARM_GP_Vector& annuity,
		const ARM_GP_Vector& fwdratevols, const ARM_GP_Vector& fwdratefwdvols, const ARM_GP_Vector& vols,
		const ARM_GP_Vector& vvol, const ARM_GP_Vector& rho, double ratecorrel, 
		int NbSimul,  int difftype)
	{
		AsOfDate			= AsOf;
		itsResetDates		= resetDates;
		itsFwdRates			= fwdRates;
		itsStrike			= strike;
		itsAnnuity			= annuity;
		itsFwdRatesVol		= fwdratevols;
		itsFwdRatesFwdVol	= fwdratefwdvols;
		itsVols				= vols;
		itsNu				= vvol;
		itsRho				= rho;
		itsRateCorrel		= ratecorrel;
		itsNbSimul			= NbSimul;
		itsDiffType			= difftype;

		itsEuropean.resize(itsResetDates.size(),0.);

		price();
	}

	ARM_Berm2DatesProxy(const ARM_Berm2DatesProxy& rhs) : 
			itsBerm(rhs.itsBerm), itsEuropean(rhs.itsEuropean)
	{
	}

	~ARM_Berm2DatesProxy() {};

	virtual ARM_Object*		Clone() const { return new ARM_Berm2DatesProxy(*this); };
	virtual string			ExportShortName() const;
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	double	GetPrice() const {return itsBerm;};
	double	GetEuro(int i) const {return itsEuropean[i];};
	void	price();

	void	sabrFwdRate(ARM_GP_Vector& fwdrate, const ARM_GP_Vector& u, double fwd, double t, double alpha, double rho, double nu, int nbPoints);
};

inline string ARM_Berm2DatesProxy::ExportShortName() const
{
	return "LPROX";
}

inline string ARM_Berm2DatesProxy::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "ARM_VBMinMaxProxy\n";

	return os.str();	
}

CC_END_NAMESPACE()

#endif
