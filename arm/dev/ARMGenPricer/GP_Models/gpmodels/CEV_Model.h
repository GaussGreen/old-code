/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file CEV_MODEL.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_CEV_MODEL_H
#define _INGPMODELS_CEV_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticIRModel.h"


CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_CEV_ModelParams;

class ARM_CEV_Model : public ARM_AnalyticIRModel
{
private:
	size_t itsIntegrationStep;

public:
	ARM_CEV_Model( const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, size_t IntegrationStep = 140 );
	ARM_CEV_Model( const ARM_CEV_Model& rhs);
    ARM_CEV_Model& operator=(const ARM_CEV_Model& rhs);
	virtual ~ARM_CEV_Model();

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
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const;

	virtual ARM_VectorPtr VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const;
    
    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {}
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string GetExportShortName() const { return "LCEVM";}
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
