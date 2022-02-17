/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarketHybridModel.h
 *
 *  \brief 
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date March 2006
 */


#ifndef _INGPMODELS_MARKETHYBRIDMODEL_H
#define _INGPMODELS_MARKETHYBRIDMODEL_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/typedef.h"

#include "AnalyticIRModel.h"
#include "gpinfra/pricingfunctionequity.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_MarketHybridModel
// \brief
//  Class for market model for hybrid models
//-----------------------------------------------------------------------------
      
class ARM_MarketHybridModel :	public ARM_AnalyticIRModel,
								public ARM_PricingFunctionEquity
{
private:
	ARM_PricingModelPtr itsDomZcModel;
	ARM_PricingModelPtr itsForZcModel;

	int itsRefEqFxOptionIdx;
	bool itsIsMarketVNS;
	bool itsIsLogNorRates;
	mutable size_t itsIrStrikeStatus;

	mutable double itsEqFxRate;
	mutable double itsEqFxStrike;
	mutable double itsIrRate;
	mutable double itsIrStrike;

public:
    enum mdmKeysAlias
    {
        MarketIrModelKey = 0,
        MarketEqFxModelKey,
		HybridDatasKey,

		NbKeys
	};

    enum zcModelDatasAlias
    {
		Time = 0,
        DomForCor,
        DomFxCor,
		ForFxCor,
		IrFxSpreadCor,
		DomO1LogNorVol,
		ForO1LogNorVol,
		IrNorVol,

		NbDatas
	};

	ARM_MarketHybridModel(const ARM_Date& asOfDate,const ARM_ObjectVector& mktDatas,const ARM_StringVector& mktKeys=ARM_StringVector(0));
	ARM_MarketHybridModel(const ARM_MarketHybridModel& rhs );
    ARM_MarketHybridModel& operator=( const ARM_MarketHybridModel& rhs )
	{
		if (&rhs != this)
		{ 
			this->~ARM_MarketHybridModel();
			new (this) ARM_MarketHybridModel (rhs);
		}
		return *this;
	}
	virtual ~ARM_MarketHybridModel() {}

	void SetZcModels(const ARM_PricingModelPtr& domZcModel,const ARM_PricingModelPtr& forZcModel)
	{ itsDomZcModel = domZcModel; itsForZcModel = forZcModel; }

	ARM_PricingModelPtr GetDomZcModel() const { return itsDomZcModel; }
	ARM_PricingModelPtr GetForZcModel() const { return itsForZcModel; }

	void SetHybridDatas(ARM_GP_T_Vector< ARM_VectorPtr >& hybridDatas);
	int GetRefEqFxOptionIdx() const { return itsRefEqFxOptionIdx; }
	void SetRefEqFxOptionIdx(int refEqFxOptionIdx) { itsRefEqFxOptionIdx = refEqFxOptionIdx; }
	bool IsMarketVNS() const { return itsIsMarketVNS; }
	void SetIsMarketVNS(bool isMarketVNS) { itsIsMarketVNS=isMarketVNS; }
	bool IsLogNorRates() const { return itsIsLogNorRates; }
	void SetIsLogNorRates(bool isLogNorRates) { itsIsLogNorRates=isLogNorRates; }
	size_t GetIrStrikeStatus()const { return itsIrStrikeStatus; }
	void ResetIrStrikeStatus() { itsIrStrikeStatus=0; }

	double GetEqFxRate() const		{ return itsEqFxRate; }
	double GetEqFxStrike() const	{ return itsEqFxStrike; }
	double GetIrRate() const		{ return itsIrRate; }
	double GetIrStrike() const		{ return itsIrStrike; }

    /// Calibration purpose
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {}
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const { return true;}

	virtual double ImpliedVol(const ARM_VanillaArg& arg) const ;

	virtual ARM_Object* Clone() const { return new ARM_MarketHybridModel(*this); }
	virtual string toString(const string& indent="",const string& nextIndent="") const { return string("ARM_MarketHybridModel"); }

	/// If one wants to export the object under XL (further use)
	virtual string ExportShortName() const { return "LMKHY";}

	virtual ARM_VectorPtr HybridCallVectorial(
		const string& modelName,
		double evalTime,
		double expiryTime,
		int callPut,
		const std::vector<double>& strikesPerState,

		/// Strip of forwards FX (or equity)
		const std::vector<double>& fxExpiryTimes,
		const std::vector<double>& fxSettlementTimes,
		const std::vector<double>& fxPayTimes,
		const std::vector<double>& fxNotionals,

		/// IR Swap
		double swapResetTime,
		const std::vector<double>& fixNotionals,
		const std::vector<double>& floatNotionals,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,

		const ARM_PricingStatesPtr& states,
		ARM_PricingContext* context=NULL) const;

	/// Not implemented specific IR or FX pricing function
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
		const ARM_PricingStatesPtr& states) const
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : not implemented");
		}

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
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : not implemented");
		}

	virtual ARM_VectorPtr Forward(
		const string& modelName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : not implemented");
		}

	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const std::vector<double>& strikesPerState,
        int callPut,
	    double payTime,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : not implemented");
		}

	virtual string GetSettlementCalendar(const string& modelName="") const;
	virtual double GetSettlementGap(const string& modelName="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
