#ifndef _INGPMODELS_MODELPARAMBGMSV2F_H
#define _INGPMODELS_MODELPARAMBGMSV2F_H

#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/curve.h"
#include "typedef.h"
#include "gpbase/vectormanip.h"
#include "gpcalib/numerical.h"
#include "gpclosedforms/smile_calibration.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsBGMSV2F : public ARM_ModelParams 
{
protected:
	
	double		itsRecorrel;	// recorrélation taux / taux
	
	double		itsMinRatio;	// ratio d'explication pour l'ACP

	int			itsNbFactors;


	double			itsV01;
	double			itsKappa1;
	double			itsTheta1;
	double			itsVVol1;
	ARM_GP_Vector	itsRho1;

	double			itsV02;
	double			itsKappa2;
	double			itsTheta2;
	double			itsVVol2;
	ARM_GP_Vector	itsRho2;

	ARM_GP_Vector	itsShift;
	ARM_GP_Vector	itsSigma;

	bool			itsLocalRho1Calib;
	bool			itsLocalRho2Calib;
	bool			itsCalibShift;
	bool			itsCalibRho1;
	bool			itsCalibRho2;

	ARM_GP_Vector	itsStdDev;

	bool			itsCalibrated;

public:

	ARM_ModelParamsBGMSV2F(const ARM_ModelParamVector& params, 
		double v01, double kappa1, double rho1, double v02, double kappa2, double rho2, double shift, 
		const ARM_GP_Vector& stddev, double recorrel = 0., int nbfactors = 0, double minratio = 1., 
		bool localrho1calib = true, bool localrho2calib = true);

	virtual ~ARM_ModelParamsBGMSV2F();

	ARM_ModelParamsBGMSV2F(const ARM_ModelParamsBGMSV2F& rhs);

public:

	virtual void		PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void		PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const	{ return new ARM_ModelParamsBGMSV2F(*this);};
	virtual string toString(const string& indent="",const string& nextIndent="") const;


	// les covariances
	double				RateRateCorrel(double t, double ithResetTime, int ithRate, double jthResetTime, int jthRate) const;

	double				GetMinRatio() const {return itsMinRatio;};

	void				SetFactorCount(int factorsNb);
	virtual size_t		FactorCount() const	{ return itsNbFactors + 2;}
	int					EffFactorCount() const {return itsNbFactors;};
	void				checkFactorCount(int nbResetDates);

	double				GetRecorrel() const {return itsRecorrel;};

	double				GetV01() const {return itsV01;};
	double				GetV02() const {return itsV02;};
	double				GetKappa1() const {return itsKappa1;};
	double				GetKappa2() const {return itsKappa2;};
	double				GetTheta1() const {return itsTheta1;};
	double				GetTheta2() const {return itsTheta2;};
	double				GetVVol1() const {return itsVVol1;};
	double				GetVVol2() const {return itsVVol2;};
	double				GetRho1(int i) const {return itsRho1[i];};
	double				GetRho2(int i) const {return itsRho2[i];};
	double				GetShift(int i) const {return itsShift[i];};
	double				GetSigma(int i) const {return itsSigma[i];};


	// calibration
	void				Calibrate(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_VanillaSecDensityPtrVector& CalibSecDensities);

	bool				GetIsCalibrated() const {return itsCalibrated;};

private:

	void				GetCalibrationData(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_VanillaSecDensityPtrVector& CalibSecDensities,
							ARM_GP_Vector& calibTimes, ARM_GP_Vector& calibFwd,
							ARM_VectorVector& calibStrikes, ARM_VectorVector& calibVols, ARM_GP_Vector& calibATM,
							ARM_GP_Vector& calibWeights, ARM_GP_Vector& atmvols);

	void				finishCalibration(const ARM_GP_Vector& resetTimes, const ARM_GP_Vector& fwdRates, 
							const ARM_GP_Vector& atmvols, ARM_GP_Vector& calibWeights,
							const ARM_GP_Vector& calibTimes, const ARM_SmileCalibration_Params_Heston2b& params);
};

CC_END_NAMESPACE()

#endif