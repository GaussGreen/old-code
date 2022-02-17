/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

#ifndef _INGPMODELS_MODELPARAMSVBGM_H
#define _INGPMODELS_MODELPARAMSVBGM_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"
#include "gpbase/vectormanip.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsSVBGM : public ARM_ModelParams 
{
protected:
	
	double		itsRecorrel;	// recorrélation taux / taux
	
	double		itsMinRatio;	// ratio d'explication pour l'ACP

	int			itsNbFactors;

public:

	ARM_ModelParamsSVBGM(const ARM_ModelParamVector& params, double recorrel = 0., int nbfactors = 0, double minratio = 1.) : 
		ARM_ModelParams(params), itsRecorrel(recorrel), itsMinRatio(minratio), itsNbFactors(nbfactors)
	{
	}

	virtual ~ARM_ModelParamsSVBGM()
	{
	}

	ARM_ModelParamsSVBGM(const ARM_ModelParamsSVBGM& rhs) : 
		ARM_ModelParams(rhs), itsRecorrel(rhs.itsRecorrel), itsMinRatio(rhs.itsMinRatio), itsNbFactors(rhs.itsNbFactors)
	{
	}

public:

	virtual void		PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void		PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const	{ return new ARM_ModelParamsSVBGM(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	virtual size_t		FactorCount() const	{ return itsNbFactors;}

	void				checkFactorCount(int nbResetDates);

	void				UpdateParamValues(ARM_GP_Vector * times, ARM_GP_Vector * shift, ARM_GP_Vector * alpha, ARM_GP_Vector * rho, ARM_GP_Vector * nu);

	// récupération des paramètres
	double				GetShift(int i) const;
	double				GetAlpha(int i) const;
	double				GetRho(int i) const;
	double				GetNu(int i) const;

	// les covariances
	double				RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const;
	double				RateVolCorrel(double t, int ithRate, int jthVol) const;
	double				VolVolCorrel(double t, int ithVol, int jthVol) const;

	double				GetMinRatio() const;

	void				SetFactorCount(int factorsNb);
};

inline double ARM_ModelParamsSVBGM::GetShift(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift) ).GetCurve()->GetOrdinates()[i];
}

inline double ARM_ModelParamsSVBGM::GetAlpha(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Alpha) ).GetCurve()->GetOrdinates()[i];
}

inline double ARM_ModelParamsSVBGM::GetRho(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter) ).GetCurve()->GetOrdinates()[i];
}

inline double ARM_ModelParamsSVBGM::GetNu(int i) const
{
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolOfVol) ).GetCurve()->GetOrdinates()[i];
}

inline double ARM_ModelParamsSVBGM::RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const
{
	if(ithRate == jthRate) return 1.;

	double beta		= ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::BetaCorrelation)).GetCurve()->Interpolate(t);
	
	double dt		= (ithResetTime - jthResetTime) / K_YEAR_LEN;
	double delta	= (ithResetTime + jthResetTime - 2. * t) / K_YEAR_LEN;

	return fabs(itsRecorrel) < K_NEW_DOUBLE_TOL || delta < K_NEW_DOUBLE_TOL ? exp(- beta * fabs(dt)) : exp( - beta * fabs(dt) / pow(delta, itsRecorrel));
}

inline double ARM_ModelParamsSVBGM::RateVolCorrel(double t, int ithRate, int jthVol) const
{
	if(ithRate == jthVol)
	{
		return GetRho(ithRate);
	}
	else
	{
		return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::CrossFactor)).GetCurve()->Interpolate(t);
	}
}

inline double ARM_ModelParamsSVBGM::VolVolCorrel(double t, int ithVol, int jthVol) const
{
	if(ithVol == jthVol)
	{
		return 1.;
	}
	else
	{
		return ((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::ReCorrelation)).GetCurve()->Interpolate(t);
	}
}

inline double ARM_ModelParamsSVBGM::GetMinRatio() const
{
	return itsMinRatio;
}

inline void ARM_ModelParamsSVBGM::checkFactorCount(int nbResetDates)
{
	itsNbFactors = itsNbFactors <= 0 || itsNbFactors > 2 * nbResetDates ? 2 * nbResetDates : itsNbFactors;
}

inline void ARM_ModelParamsSVBGM::UpdateParamValues(ARM_GP_Vector * times,
													ARM_GP_Vector * shift, ARM_GP_Vector * alpha, ARM_GP_Vector * rho, ARM_GP_Vector * nu)
{
	((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Shift)).SetValuesAndTimes(times, shift);

	((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::Alpha)).SetValuesAndTimes(times, alpha);

	((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::VolOfVol)).SetValuesAndTimes(times, nu);

	((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::QParameter)).SetValuesAndTimes(times, rho);
}

inline void ARM_ModelParamsSVBGM::SetFactorCount(int factorsNb)
{
	itsNbFactors = factorsNb;
}

CC_END_NAMESPACE()

#endif
