/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SABR_Model.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_SABR_MODEL_H
#define _INGPMODELS_SABR_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticIRModel.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpcalib/densityfunctors.h"

CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_SABR_ModelParams;
struct ARM_VanillaArg;

class ARM_SABR_Model : public ARM_AnalyticIRModel
{
private:
	size_t itsIntegrationStep;
	int itsTypeOfModel;

    /// --------------- fast calibration caching arguments -----------------
    CC_IS_MUTABLE ARM_VanillaArg* itsCurrentArg;

public:
	ARM_SABR_Model( const ARM_ZeroCurvePtr& zc, 
		const ARM_SABR_ModelParams& params, 
		const long& TypeOfModel= ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT, 
		size_t IntegrationStep = 140 );

	ARM_SABR_Model(const ARM_SABR_Model& rhs);
    ARM_SABR_Model& operator=(const ARM_SABR_Model& rhs);
	virtual ~ARM_SABR_Model();

	//  Update the density functor at expiryTime and tenor
	virtual void UpdateDensityFunctor(double fwd, double expiryTime, double tenor);

    /// Reconstruction formula
    virtual ARM_VectorPtr VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;

    virtual double ImpliedVol(const ARM_VanillaArg& arg) const;

    /// accessors
    inline ARM_VanillaArg* GetVanillaArg() const { return itsCurrentArg; }

    /// To calculate Implied Volatility partial derivatives
    /// in ralation to each ModelParam
	virtual bool HasClosedFormsDerivatives( ARM_ModelParamType::ParamNb paramType, size_t factorNb ){ return true; }
	virtual double PartialDerivative( const ARM_ModelParam& modelParam, size_t number, 
        size_t factorNb,const ARM_VanillaArg& arg,
        ARM_MktTargetType targetFuncType = ARM_CalibrationTarget::PriceTarget);
    
    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter);
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;


	/// --------- numerical method bool flags + validation for using SABR with CFnumMethod
    virtual bool SupportBackwardInduction() const {	return true; }
    virtual bool SupportForwardInduction()  const {	return true; }

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string GetExportShortName() const { return "LSAMO"; }
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
