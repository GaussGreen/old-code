/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Model.h
 *
 *  \brief 
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_HESTON_MODEL_H
#define _INGPMODELS_HESTON_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "AnalyticIRModel.h"

CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_Heston_ModelParams;

class ARM_Heston_Model : public ARM_AnalyticIRModel
{
private:
	size_t itsIntegrationStep;

public:
	ARM_Heston_Model( const ARM_ZeroCurvePtr& zc, 
		const ARM_Heston_ModelParams& params, 
		size_t IntegrationStep = 140 );
	ARM_Heston_Model(const ARM_Heston_Model& rhs);
	ASSIGN_OPERATOR(ARM_Heston_Model)
	virtual ~ARM_Heston_Model() {};

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
	virtual string GetExportShortName() const { return "LHEMO";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

