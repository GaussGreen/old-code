/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarketIRModel.h
 *
 *  \brief 
 *	\author  A. Chaix
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPMODELS_MARKETIR_MODEL_H
#define _INGPMODELS_MARKETIR_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "AnalyticIRModel.h"

CC_BEGIN_NAMESPACE( ARM )
///
///
/// Idea :  This class is designed to enable new methods for the pricing of market instruments.
///		    It uses the ARM_MarketData_ManagerRep* attribute of ARM_PricingModel and ignores
///			its model params.
///			It has been created for the pricing of variable notional swaptions using
///			BS vols and correlations
///
class ARM_MarketIRModel : public ARM_AnalyticIRModel
{
public:
	/// how do we choose the strikes at which normal vols are computed
	enum VnsPricingMethod
	{
		MONEYNESS = 0,
		ATM 
	};


private:
	string itsZeroCurveKey;
	string itsBsModelKey;
	
	VnsPricingMethod itsVnsPricingMethod;

	/// To weigh up the moyeness level
	double itsMoyenessLevel;

	/// mutable attributes to keep memory of last variable notional swaption pricing
	mutable double itsVnsForward;
	mutable double itsVnsStrike;
	mutable ARM_GP_Matrix itsVnsBasketCorrels;
	mutable double itsVnsBasketVol;
	mutable double itsVnsNumeraire;
	mutable vector<ARM_Date> itsVnsStartDates;
	mutable vector<ARM_Date> itsVnsEndDates;
	mutable ARM_GP_Vector itsVnsFloatNotionals;
	mutable ARM_GP_Vector itsVnsStrikes;
	mutable ARM_GP_Vector itsVnsForwards;
	mutable ARM_GP_Vector itsVnsAtmNormVols;
	mutable ARM_GP_Vector itsVnsNormVols;
	mutable ARM_GP_Vector itsVnsBasketCoefs;

public:
	ARM_MarketIRModel (const ARM_MarketData_ManagerRep& mktDataManager,
		const ARM_StringVector& mdmKeys, 
		VnsPricingMethod method = MONEYNESS,
		double moyenessLevel = 1.0);
	ARM_MarketIRModel (const ARM_MarketIRModel& rhs);
	ASSIGN_OPERATOR(ARM_MarketIRModel)
	virtual ~ARM_MarketIRModel();

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

	/// Variable Notional swaptions can be priced under Analytic IR model
	virtual bool ClosedFormulaSwaptionFlag(
		bool isConstantNominal,
		bool isConstantSpread,
		bool isConstantstrike) const { return isConstantSpread&&isConstantstrike;}

	/// to bypass default key setting
	virtual void SetZeroCurveKey (const string& zeroCurveKey) {itsZeroCurveKey = zeroCurveKey;};
	virtual void SetBsModelKey (const string& bsModelKey){itsBsModelKey = bsModelKey;};
    string GetBsModelKey() const {return itsBsModelKey;};

	double GetVnsForward() const { return itsVnsForward; }
	double GetVnsBasketVol() const { return itsVnsBasketVol; }
	double GetVnsNumeraire() const { return itsVnsNumeraire; }
	const ARM_GP_Vector& GetBasketCoefs() const {return itsVnsBasketCoefs;}
	
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
	//virtual string GetExportShortName() const { return "LMKMO";}
	virtual string ExportShortName() const { return "LMKMO";}

	//double CptCvxAdj(double cmRate, int tenor, double expiry, double IndexVol, int payFreq=1) const;
	//double CptTimeLagAdj(double expiry, double lag, double IndexVol, double Forward, double DiscountVol, double Rho=0.7, int methUsed=1) const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
