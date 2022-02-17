/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Local_Model.h
 *
 *  \brief base class for local analytical model
 *	\author  A. Chaix
 *	\version 1.0
 *	\date June 2005
 */


#ifndef _INGPMODELS_LOCAL_MODEL_H
#define _INGPMODELS_LOCAL_MODEL_H

#include "gpbase/env.h"
#include "gpbase/port.h"

#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingfunctionir.h"
#include "gpinfra/pricingfunctionequity.h"
#include "gpinfra/pricingfunctioninflation.h"
#include "gpinfra/pricingmodeltype.h"

#define UNIMPLEMENTED_PRICING_FUNCTION {ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function"); return ARM_VectorPtr();}
#define UNIMPLEMENTED_PRICING_MATFUNCTION {ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function"); return ARM_GP_MatrixPtr();}

class ARM_Portfolio;
class ARM_Security;
class ARM_CapFloor;
class ARM_SpreadOption;
class ARM_Option;
class ARM_CorridorLeg;

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_MultiAssetsModel;
class ARM_ModelParams ;
class ARM_DensityFunctor;
class ARM_GenCalculator;

class ARM_Local_Model : public ARM_PricingModel,
						public ARM_PricingFunctionIR ,
						public ARM_PricingFunctionEquity
						/// if needed, derivation w.r.t. ARM_PricingFuncInflation can be added
						/* public ARM_PricingFuncInflation */
{
// ---- attributes

private:
	// numerical model (HW, SFRM,...) used to provide (unadjusted) forward rates 
	CC_IS_MUTABLE ARM_PricingModel*		itsNumericalModel; 
	ARM_GP_Matrix						itsVarianceSqueezInfo; 
	bool								itsVarianceSqueezeStatus;
	bool								itsResetCalib;
	bool								itsCorrelUnSqueezer;
	bool								itsVolUnSqueezer;
	bool								itsPVAdjuster;

protected:
	ARM_VectorPtrVector					itsFunctionals;
	ARM_GP_VectorPtr					itsGrid;
	ARM_GP_Vector						itsFwds;
	ARM_GP_Vector						itsVols;
	ARM_GP_Vector						itsShifts;
	ARM_GP_Vector						itsResetTimes;
	bool								itsFunctionalFlag;

public:
	///constructors/destructors
	ARM_Local_Model (const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params);
	ARM_Local_Model (const ARM_Local_Model& rhs);
	virtual ~ARM_Local_Model() {};

	/// utilities
	virtual string toString(const string& indent="",const string& nextIndent="") const;
	
	///accessors
	ARM_PricingModel* GetNumericalModel () const;
	void SetNumericalModel (ARM_PricingModel* numericalModel);
	ARM_GP_Matrix GetVarianceSqueezInfo() const { return itsVarianceSqueezInfo;}
	void SetVarianceSqueezInfo (ARM_GP_Matrix& variancesqueezInfo) {itsVarianceSqueezInfo = variancesqueezInfo;};
	void push_backVarianceSqueezInfo( const ARM_GP_Vector& variancesqueezInfo) {itsVarianceSqueezInfo.push_backRow(variancesqueezInfo);}
	bool GetVarianceSqueezeStatus() const { return itsVarianceSqueezeStatus;}
	void SetVarianceSqueezeStatus (bool variancesqueezStatus) {itsVarianceSqueezeStatus = variancesqueezStatus;};

	bool GetResetCalib() const { return itsResetCalib; }
	void SetResetCalib(bool resetCalib) { itsResetCalib=resetCalib; }

	bool GetCorrelUnSqueezer() const { return itsCorrelUnSqueezer; }
	void SetCorrelUnSqueezer(bool correlUnSqueezer) { itsCorrelUnSqueezer=correlUnSqueezer; }

	bool GetVolUnSqueezer() const { return itsVolUnSqueezer; }
	void SetVolUnSqueezer(bool volUnSqueezer) { itsVolUnSqueezer=volUnSqueezer; }

	bool GetPVAdjuster() const { return itsPVAdjuster; }
	void SetPVAdjuster(bool PVAdjuster) { itsPVAdjuster=PVAdjuster; }

	/// services

	/// update links
	virtual void UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel );
	
	virtual void ResetModelParams () {	ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented"); }

	//// returns the DF computed from itsNumericalModel
	virtual ARM_VectorPtr DiscountFactor (	const string& curveName,
											double evalTime, 
											double maturityTime,
											const ARM_PricingStatesPtr& states) const;


	virtual void CalibrateLocalModel (	const ARM_Portfolio& portfolio,
										const ARM_GP_Vector& evalTimes);

	virtual void CalibrateLocalModel (	const ARM_Security& security, 
										double targetPrice, 
										const ARM_GP_Vector& evalTimes,
										size_t secIdx=0) {	ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented"); }

	virtual void CalibrateLocalModel (	const ARM_GenCalculator& calculator, 
										const ARM_GP_Vector& evalTimes) { ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented"); }

// NEW METHODS for local models
	
	virtual void CalibrateLocalModelFunctional (
										vector<ARM_Security*> securities,
										vector<ARM_DensityFunctor*> densities,
										int sizeGrid = 501,
										double nbStdDev = 6.,
										bool rescaling = true);

	virtual void CalibrateLocalModel (	const ARM_Security& security,ARM_DensityFunctor& density,bool rescaling) 
											{	ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented"); }

	inline void SetFunctionalFlag ( bool flag )			{ itsFunctionalFlag = flag; };
	inline void setShift ( double shift )				{ itsShifts.push_back(shift); }
	inline void setVol ( double vol )					{ itsVols.push_back(vol); }
	inline void setFwd ( double fwd )					{ itsFwds.push_back(fwd); }

	inline void setFunc(ARM_GP_VectorPtr func)			{ itsFunctionals.push_back(func); }
	inline void setTime(double time)					{ itsResetTimes.push_back(time); }

	inline const ARM_GP_VectorPtr& Grid() const			{ return itsGrid; }
	inline const ARM_GP_Vector& ResetTimes() const		{ return itsResetTimes; }
	inline size_t GetGridSize()	const					{ return itsGrid->size(); }
	inline double GetGrid(size_t i)	const				{ return (*itsGrid)[i]; }
	inline bool GetFunctionalFlag() const				{ return itsFunctionalFlag; };

	inline double GetShift(size_t i) const				{ return itsShifts[i]; }
	inline double GetVol(size_t i) const				{ return itsVols[i]; }
	inline double GetFwd(size_t i) const				{ return itsFwds[i]; }
	

	inline double FuncValue(size_t idx,size_t j) const
											{ return itsFunctionals[idx]->Elt(j); }
	virtual ARM_VectorPtr Func(double evalTime,const ARM_GP_VectorPtr& values) const
											{	ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented"); }
// END NEW METHODS


/// utilities for calibration
protected:
	/// caplet / floorlet
	virtual void GetCapFloorLetData (	ARM_CapFloor* capFloor,
										size_t periodIdx,
										// outputs
										double& resetTime,
										double& payTime,
										double& payPeriod,
										double& payNotional,
										double& liborStartTime,
										double& liborEndTime,
										double& liborPeriod,
										double& strike,
										int&	callPut) const ;
	/// spreadoptionlet
	virtual void GetSpreadOptionLetData (	ARM_SpreadOption* spreadOption,
											size_t periodIdx,
											// outputs
											int&			callPut,
											double&			resetTime,
											double&			payTime,
											double&			payPeriod,
											double&			payNotional,
											double&			coeffLong,
											double&			coeffShort,
											double&			strike,
											double&			swapLongFloatStartTime,
											double&			swapLongFloatEndTime,
											ARM_GP_Vector&	swapLongFixPayTimes,
											ARM_GP_Vector&	swapLongFixPayPeriods,
											double&			swapShortFloatStartTime,
											double&			swapShortFloatEndTime,
											ARM_GP_Vector&	swapShortFixPayTimes,
											ARM_GP_Vector&	swapShortFixPayPeriods)  const ;

	/// spreadoptionlet
	virtual void GetCorridorSpreadOptionData (	ARM_SpreadOption* spreadOption,
										// outputs
										int&			callPut,
										ARM_GP_Vector&	resetTimes,
										ARM_GP_Vector&	payTimes,
										ARM_GP_Vector&	payPeriods,
										ARM_GP_Vector&	payNotionals,
										ARM_GP_Vector&	coeffsLong,
										ARM_GP_Vector&	coeffsShort,
										ARM_GP_Vector&	strikes,
										ARM_GP_Vector&	swapLongFloatStartTimes,
										ARM_GP_Vector&	swapLongFloatEndTimes,
										vector<ARM_GP_Vector>&	swapsLongFixPayTimes,
										vector<ARM_GP_Vector>&	swapsLongFixPayPeriods,
										ARM_GP_Vector&	swapShortFloatStartTimes,
										ARM_GP_Vector&	swapShortFloatEndTimes,
										vector<ARM_GP_Vector>&	swapShortFixPayTimes,
										vector<ARM_GP_Vector>&	swapShortFixPayPeriods,
										ARM_GP_Vector& fixValues,
										double& spread1,
										double& spread2,
										int& payIndexType,
										ARM_GP_Vector& payIndexLeverages,
										ARM_GP_Vector& swapPayFloatStartTime,
										ARM_GP_Vector& swapPayFloatEndTime,
										vector<ARM_GP_Vector>& swapPayFixPayTimes,
										vector<ARM_GP_Vector>& swapPayFixPayPeriods,
										ARM_GP_Vector& payIndexResetTimes,
										ARM_IntVector& periodIndex,
										int& rateCallPut,
										ARM_GP_Vector& rateStrikes,
										ARM_GP_Vector&	rateFloatStartTimes,
										ARM_GP_Vector&	rateFloatEndTimes,
										vector<ARM_GP_Vector>&	rateFixPayTimes,
										vector<ARM_GP_Vector>&	rateFixPayPeriods,
										double& rateSpread1,double& rateSpread2)  const ;

	virtual void GetCorridorLegData (	ARM_CorridorLeg* corridor,
											// outputs
											int&			callPut,
											double&			resetTime,
											double&			payTime,
											int&			indexPaymentType,
											double&			fwdPaymentPeriod,
											ARM_GP_Vector&	refIdxResettimes ,
											ARM_GP_Vector&	refIdxStarttimes,
											ARM_GP_Vector&	refIdxEndtimes,
											ARM_GP_Vector&	refFwdPeriods,
											ARM_GP_Vector&	refIndexWeight,
											ARM_GP_Vector&	DownBarrier,
											ARM_GP_Vector&	UpBarrier,
											double&			couponMargin,
											double&			payNotional)  const ;


	/// Equity/Fx option
	virtual void GetEqFxOptionData(	ARM_Option* option,
									// outputs
									double& resetTime,
									double& settlementTime,
									double& payTime,
									double& strike,
									int&	callPut) const ;

	//------------------------------
	//--- implemented IR functions
	//------------------------------

	// libor
	virtual ARM_VectorPtr Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const = 0;

	// caplet
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
		const ARM_PricingStatesPtr& states) const = 0;

    /// Vanilla swaption function (vectorial strike version)
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
		const ARM_GP_Vector& fixTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_GP_Matrix& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		bool isConstantNotional = true,
		bool isConstantSpread = true,
		bool isConstantStrike = true) const = 0;

	// spread option
	virtual ARM_VectorPtr  VanillaSpreadOptionLet(
			const string& curveName,
			double evalTime,
			int callPut,
			double startTime,
			double endTime,
			double resetTime,
			double payTime,
			double payPeriod,
			double notional,
			double coeffLong,
			double coeffShort,
			const ARM_GP_Vector& strikes,
			double swapLongFloatStartTime,
			double swapLongFloatEndTime,
			const ARM_GP_Vector& swapLongFixPayTimes,
			const ARM_GP_Vector& swapLongFixPayPeriods,
			double swapShortFloatStartTime,
			double swapShortFloatEndTime,
			const ARM_GP_Vector& swapShortFixPayTimes,
			const ARM_GP_Vector& swapShortFixPayPeriods,
			const ARM_PricingStatesPtr& states) const = 0;



	// check params
	virtual bool ValidateModelParams(const ARM_ModelParams& params) const = 0;

	
	//------------------------------
	//--- unimplemented IR functions
	//------------------------------
	
	/// Annuity function
	virtual ARM_VectorPtr Annuity(
		const string& curveName, 
        double evalTime,
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_FUNCTION;

	/// function to avoid computing twice the fixLegAnnuity
	/// if you want to keep the annuity value, make sure you have cloned it before!
	virtual ARM_VectorPtr SwapRateInPlaceWithComputedAnnuity(
		const string& curveName, 
		double evalTime,
		double floatStartTime, 
		double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes,
        const ARM_GP_Vector& fwdEndTimes,
        const ARM_GP_Vector& fwdPayPeriods,
		const ARM_GP_Vector& floatPayTimes,
        const ARM_GP_Vector& floatPayPeriods,
        const ARM_GP_Vector& margin,
        bool isDbleNotional,
		const ARM_VectorPtr& PreComputedAnnuity,
		const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_FUNCTION;

	/// Swap function (vectorial fixed rate version)
	virtual ARM_VectorPtr NPVSwap(
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& floatPayTimes, 
		const ARM_GP_Vector& floatPayPeriods, 
		const ARM_GP_Vector& margin,
        bool isDbleNotional,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Vector& FloatNotional,
		const ARM_GP_Matrix&  strikesPerState,
        int payRec,
        const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_FUNCTION;

	virtual ARM_GP_MatrixPtr NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fwdStartTimes, 
		const ARM_GP_Vector& fwdEndTimes, 
		const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& PayTimes, 
		const ARM_GP_Vector& PayPeriods, 
		const ARM_GP_Vector& margin, 
		const ARM_GP_Vector& notional, 
		const ARM_PricingStatesPtr& states) const  UNIMPLEMENTED_PRICING_MATFUNCTION;

	virtual ARM_GP_MatrixPtr NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_MATFUNCTION;

    /// Implied volatility function to calculate 
    /// sqrt[1/resetTime * integral(evalTime-> resettime
    /// sigme*sigme) du]
    virtual double ImpliedVol(const ARM_VanillaArg& arg) const {ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function"); return 0;};

   	
	virtual ARM_VectorPtr LiborWithControlVariate( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
        const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_FUNCTION;



	//----------------------------------
	//--- unimplemented EQ/FX functions
	//----------------------------------
	
	/// Forward function
	virtual ARM_VectorPtr Forward(
		const string& curveName, 
        double evalTime,
	    double expiryTime,
	    double settlementTime,
	    double payTime,
        const ARM_PricingStatesPtr& states) const UNIMPLEMENTED_PRICING_FUNCTION;

	/// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& curveName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
        int callPut,
	    double payTime,
	    const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const UNIMPLEMENTED_PRICING_FUNCTION;


	/// convention support
	virtual string GetSettlementCalendar(const string& modelName="") const	{ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function"); return string();};
	virtual double GetSettlementGap(const string& modelName="") const		{ARM_THROW( ERR_INVALID_ARGUMENT, "unimplemented pricing function"); return 0;};;

	
	//-----------------------------------------------------------------------------
	//--- unfortunately, all services of ARM_PricingModel must be implemented...
	//-----------------------------------------------------------------------------
	virtual void PreProcessing(ARM_ModelFitter& modelFitter) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {}
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter) {};
    virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod) {};
	virtual void AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* modelParam, size_t factorNb = 0 ) {};
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos){return ARM_PricingStatesPtr();};
	virtual ARM_PricingStatesPtr FirstPricingStates( size_t bucketSize ) const{return ARM_PricingStatesPtr();};
	virtual ARM_GP_Vector* ComputeModelTimes( const ARM_TimeInfoPtrVector& timeInfos ){return NULL;};
// FIXMEFRED: mig.vc8 (25/05/2007 11:28:24):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}
	virtual bool SupportBackwardInduction() const {return true;};
	virtual bool SupportForwardInduction() const  {return true;};	
	virtual bool SupportAnalyticMarginal() const  {return true;};
	virtual void NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const {};
	virtual void ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,ARM_MatrixVector& localVariances ) const {};
	virtual bool NeedsToCholeskyDecomposeFactors( ) const {return false;}
    virtual double VarianceToTime(double var,double minTime=0.0,double maxTime=5*K_YEAR_LEN) const {return 0.0;};
	virtual void MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const {};
	virtual void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const {};
	virtual int GetType() const {return MT_NON_STOCHASTIC_MODEL;};

};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
