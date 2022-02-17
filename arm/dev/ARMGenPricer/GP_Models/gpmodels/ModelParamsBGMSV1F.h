/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

#ifndef _INGPMODELS_MODELPARAMBGMSV1F_H
#define _INGPMODELS_MODELPARAMBGMSV1F_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"
#include "gpbase/vectormanip.h"
#include "gpcalib/numerical.h"
#include "gpclosedforms/smile_calibration.h"
#include "gpclosedforms/heston_pricer.h"
#include "gpclosedforms/vanilla_bs.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsBGMSV1F : public ARM_ModelParams 
{
protected:
	
	double			itsRecorrel;	// recorrélation taux / taux
	
	double			itsMinRatio;	// ratio d'explication pour l'ACP

	int				itsNbFactors;

	bool			itsLocalRhoCalib;
	bool			itsGlobalCalib;

	ARM_GP_Vector	itsStdDev;

	bool			itsCalibrated;
	ARM_GP_Vector	itsVVar;

public:

	ARM_ModelParamsBGMSV1F(const ARM_ModelParamVector& params, const ARM_GP_Vector& stddev, double recorrel = 0., int nbfactors = 0, double minratio = 1., bool localrhocalib = true,
		bool globalCalib = true);

	virtual ~ARM_ModelParamsBGMSV1F();

	ARM_ModelParamsBGMSV1F(const ARM_ModelParamsBGMSV1F& rhs);

public:

	virtual void		PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void		PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const	{ return new ARM_ModelParamsBGMSV1F(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	virtual size_t		FactorCount() const	{ return itsNbFactors + 1;}
	int					EffFactorCount() const {return itsNbFactors;};

	void				checkFactorCount(int nbResetDates);

	// récupération des paramètres
	double				GetV0() const;
	double				GetKappa() const;
	double				GetTheta() const;
	double				GetVVol(double time) const;
	double				GetVVol(int i) const;
	double				GetRho(int i) const;
	double				GetShift(int i) const;
	double				GetLevel(int i) const;

	ARM_Curve *			GetVVolt(double evalTime = 0.) const;

	// les covariances
	double				RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const;

	double				GetMinRatio() const;

	void				SetFactorCount(int factorsNb);

	double				GetRecorrel() const {return itsRecorrel;};

	// calibration
	void				Calibrate(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_VanillaSecDensityPtrVector& CalibSecDensities);

	bool				GetIsCalibrated() const {return itsCalibrated;};

	void				SetParams(const ARM_GP_Vector& resetTimes, double v0, double kappa, double theta,
							const ARM_GP_Vector& shifts, const ARM_GP_Vector& levels, const ARM_GP_Vector& vvols,
							const ARM_GP_Vector& rhos);

protected:

	void				GetCalibrationData(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_VanillaSecDensityPtrVector& CalibSecDensities,
							ARM_GP_Vector& calibTimes, ARM_GP_Vector& calibFwd,
							ARM_VectorVector& calibStrikes, ARM_VectorVector& calibVols, ARM_GP_Vector& calibATM,
							ARM_GP_Vector& calibWeights, ARM_GP_Vector& atmvols);

	void				finishGlobalCalibration(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_GP_Vector& atmvols, ARM_GP_Vector& calibWeights,
							const ARM_GP_Vector& calibTimes, const ARM_SmileCalibration_Params_Heston& params);

	void				finishLocalCalibration(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_GP_Vector& atmvols, ARM_GP_Vector& calibWeights,
							const ARM_GP_Vector& calibTimes, const ARM_GP_Vector& shifts, 
							const ARM_GP_Vector& levels, const ARM_GP_Vector& rhos, const ARM_GP_Vector& vvols);


	class CATMLocalConstr : public DoubleToDoubleFunc
	{
	protected:
		ARM_HestonOptionPricerVVolt	*	pricer;
		ARM_ImpliedVolBS *				inverse;

	public:
		CATMLocalConstr(ARM_HestonOptionPricerVVolt * hestonpricer, ARM_ImpliedVolBS * invfunc) : 
		  pricer(hestonpricer), inverse(invfunc)
		{}

	public:
		double operator()(double x) const
		{
			pricer->SetLevel(x);
			return inverse->vol(pricer->price());
		}

	};
};



CC_END_NAMESPACE()

#endif
