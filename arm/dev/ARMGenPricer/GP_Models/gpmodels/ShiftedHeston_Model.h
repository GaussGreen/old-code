/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ShiftedHeston_Model.h
 *
 *  \brief 
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_SHIFTEDHESTON_MODEL_H
#define _INGPMODELS_SHIFTEDHESTON_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticIRModel.h"

CC_BEGIN_NAMESPACE( ARM )

//// forward declaration
class ARM_ShiftedHeston_ModelParams;

class ARM_ShiftedHeston_Model : public ARM_AnalyticIRModel
{
private:
	size_t itsIntegrationStep;
	bool itsIsMCPrice;
	size_t itsNbStepsMC;
	size_t itsNbSimulations;

public:
	ARM_ShiftedHeston_Model( const ARM_ZeroCurvePtr& zc, const ARM_ShiftedHeston_ModelParams& params, size_t IntegrationStep = 140 ,bool isMCPrice = false, size_t nbSteps = 12, size_t nbSimulations = 100000);
	ARM_ShiftedHeston_Model(const ARM_ShiftedHeston_Model& rhs);
    ARM_ShiftedHeston_Model& operator=(const ARM_ShiftedHeston_Model& rhs);
	virtual ~ARM_ShiftedHeston_Model();

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
		bool isConstantNotional = true ,
		bool isConstantSpread = true ,
		bool isConstantStrike = true ) const;
    
	ARM_VectorPtr MCPrice(double swapResetTime, double swapRate, double strike, int callPut) const;
	virtual bool SupportForwardInduction()  const {	return true;}

    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {}
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const;

	/// method to advise the break point times
    virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb=0 );
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);

    /// Standard ARM object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	virtual string ExportShortName() const { return "LSHMO";}
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/

