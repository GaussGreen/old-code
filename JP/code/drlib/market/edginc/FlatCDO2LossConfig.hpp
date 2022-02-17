//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   File name   : FlatCDO2LossConfig.hpp
//
//   Date        : 04-Oct-2006
//
//   Description : Interface for class FlatCDO2LossConfig, which defines a 
//                 LossConfig consisting of piecewise linear payoffs defined 
//                 for at different points in time and on the loss of a 
//                 CDOPortfolio 
//                 
//----------------------------------------------------------------------------

#ifndef QLIB_FLATCDO2LOSSCONFIG_HPP
#define QLIB_FLATCDO2LOSSCONFIG_HPP

#include "edginc/Atomic.hpp"
#include "edginc/Array.hpp"
#include "edginc/Theta.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/MarketWrapper.hpp"
#include "edginc/ICreditLossConfig.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/IEffectiveCurveLossModelConfig.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"
#include "edginc/CDOParallelStrikeShift.hpp"
#include "edginc/CreditLossConfigMC.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/CreditEngineParameters.hpp"
#include "edginc/CDOParallelStrikeShift.hpp"
#include "edginc/Interpolator.hpp"
#include "edginc/ITranchesCombinationPayoff.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IProtectionProvider);
FORWARD_DECLARE(GeneralisedCDO);
FORWARD_DECLARE(CDOPortfolio);
FORWARD_DECLARE(SingleCreditAsset);
FORWARD_DECLARE(FlatCDO2LossConfig);

class MARKET_DLL FlatCDO2LossConfig : 
	public CObject,
    virtual public ICreditLossConfig,
	// virtual public ICreditLossConfigSVGenMC,  // so that it supports MC
	// virtual public ICreditLossConfigIndexedSVGenMC, out for now
    virtual public IEffectiveCurveLossModelConfig::IIntoLossGen,
	virtual public ITranchesCombinationPayoff
{

public:

	FORWARD_DECLARE(LinearLossProfile);

	FORWARD_DECLARE(TrancheProfile);

    //	typedef smartPtr<LinearLossProfileArray> LinearLossProfileArraySP;

	FlatCDO2LossConfig(CClassConstSP clazz = TYPE) : 
	CObject(clazz) {};

	static bool FlatCDO2LossConfigLoad();

	static void load(CClassSP& clazz); 

	virtual ~FlatCDO2LossConfig();
	
	FlatCDO2LossConfig(
		DateTime valDate,
		CDOPortfolioSP driver,
		LinearLossProfileArraySP profiles,
		DateTimeArrayConstSP timeLine,
		CClassConstSP clazz = TYPE);

	static IObject* defaultFlatCDO2LossConfig() {	return new FlatCDO2LossConfig(TYPE);};

    /** All instances have a 'name' which can be viewed as the name of
        the portfolio or the name of the single curve. */
    virtual string getName() const;

	////////////////////////////////

    /** Returns the historic losses on the contingent leg, due to the
        settlement of this name's default. In the presence of shorts 
        the 'losses' can be negative.
        If there are no losses, it may return an empty (0 size) 
        CtgLegLossPerDefaultArraySP or a "null" SP.
        Note: the cashflow dates are the CALCULATION DATES (as opposed to 
        settlement dates) when the corresponding losses are determined. */
    virtual CtgLegLossPerDefaultArraySP historicContingentLegLosses(
         CIntConstSP triggerDelay,
         CIntConstSP defaultToCalculationDelay,
         const DateTime& lastTriggerDate,
         IBadDayAdjusterConstSP bda,
		 const IProtectionProvider* const protect) const 
	{
		return 	CtgLegLossPerDefaultArraySP(new CtgLegLossPerDefaultArray());
	};

    /** Returns the notional reductions on the fee leg (due to losses and/or
        recovered notional) and the corresponding fee rebates if any. In the
        presence of shorts the 'losses' can be negative
        The accrualPeriods are potentially required to determine the actual  
        date of the notional reductions. */
    virtual FeeLegReductionPerDefaultArraySP historicFeeLegReductions(
        CIntConstSP               triggerDelay,
        CIntConstSP               defaultToCalculationDelay,
        const double              temporaryLossAmount,
        const DateTime&           lastTriggerDate,    
        AccrualPeriodArrayConstSP accrualPeriods,
        IBadDayAdjusterConstSP    bda,
		const bool                recoverNotional) const 
	{
        return FeeLegReductionPerDefaultArraySP();
	};

    /** Compute the rebate payments caused by historic fee leg reductions */
    virtual CashFlowArraySP historicRebatePayments(
        const IRebateCalculator* const   rebateCalc,
        FeeLegReductionPerDefaultArraySP reductions,
        IForwardRatePricerSP             model,
		const bool						 recoverNotional) const 
		{
			return	CashFlowArraySP(new CashFlowArray());
		};

	/** Returns the engine parameters for this loss cfg (NB Model needs to
		have selected appropriate set during getMarket()). The parameter
		can specify what type of engine parameters are required. An
		exception is thrown if the relevant type is not found */
	virtual CreditEngineParametersConstSP getEngineParams(
		CClassConstSP engineParamsType) const
	{
		static const string method("FlatCDO2LossConfig::getEngineParams");
		try {
			if (!engineParamsType){
				throw ModelException(method, "The type of engine parameters "
									"requested is Null");
			}
			if (!engineParams){
				throw ModelException(method, "Null engine parameters found but "
									"require type: " + engineParamsType->getName());
			}

			return engineParams->getEngineParams(engineParamsType);
		}
		catch (exception& e) {
			throw ModelException(e, method, "Failed when trying to retrieve credit "
								"engine parameters of type: " +
								engineParamsType->getName());
		}
	};
	////////////////////////////////

	virtual int numInnerLossConfigs() const;

	virtual CDOPortfolioConstSP getPortfolio() const
	{
		return lossDriver;
	};

	virtual void getMarket(const IModel* model, const MarketData* market);

	double map(int slice, double x) const;

	static CClassConstSP const TYPE;

	virtual DateTime getToday() const { return valueDate;};

	virtual SingleCreditAssetConstSP nameAsset(int index) const {
		return lossDriver->nameAsset(index);
	}; // lossDriver->getPortfolio->getName(i) ?

	virtual IEffectiveCurveLossGenSP lossGenerator(
		IEffectiveCurveLossModelConfigConstSP effCurveLossModelConfig) const;

	virtual ICreditLossConfigConstSP getInnerLossConfig(const int index) const 
	{
		return lossDriver->getInnerLossConfig(index);
	}

	virtual double maxLossGivenDefault() const {return tslcNotional;};

	virtual double notional() const {return tslcNotional;};

	DateTimeArrayConstSP getTimeline() const { return timePoints;};

	LinearLossProfileArrayConstSP getLossProfiles() const 
	{
		return lossProfiles;
	};

    
	virtual double maxPossibleLoss() const {return tslcNotional;};

	virtual void currentLossAndRecovered(
		double& loss,                                        // (O)
        double& recoveredNotional,							 // (O)
        const bool recoverNotional) const;
    
	double getPrepaidNotional(const DateTime &prepayDate) const 
	{
		return 0;
	};

	int findProfileIdx(DateTime time) const;

	double trancheOutstandingNotional(const bool recoverNotional) const;

	void validatePop2Object();

private:

	// a function representation, limited to piece-linear for now

	CDOPortfolioSP lossDriver;

	DateTime valueDate;      // Today

	DateTimeArrayConstSP timePoints;

	int firstFutureTP;

	LinearLossProfileArraySP lossProfiles;

	CreditEngineParametersWrapper engineParams; // Parameters for this loss config

	double tslcNotional;

public:

	class MARKET_DLL LinearLossProfile : public CObject
	{

		double maxValue;
		
		DoubleArray midPoints;
		
		DoubleArray midValues;

		DoubleArray stdDevs;

		DoubleArray lowerBounds;

		DoubleArray upperBounds;

	public:

		friend class FlatCDO2LossConfig;

		virtual ~LinearLossProfile() {};

		static CClassConstSP const TYPE;

		LinearLossProfile(CClassConstSP clazz = TYPE);

		double map(double x) const;

		static void load(CClassSP& clazz);

		// InterpolantNonVirtual interp; 

		LinearLossProfile(
			const DoubleArray & mx,
			const DoubleArray & my,
			CClassConstSP clazz = TYPE);

		LinearLossProfile(
			const DoubleArray & mx,
			const DoubleArray & my,
			const DoubleArray & lowerBounds,
			const DoubleArray & upperBounds,
			const DoubleArray & stdDevs,
			CClassConstSP clazz = TYPE); 

		// need to improve this, it should be the max loss
		double maxLoss() const;

		static IObject* defaultLinearLossProfile();

	};

	virtual void linearDecomposition(
        const DateTime& time,
        DoubleArray& baseStrikes,        /* output */
        DoubleArray& baseStrikesWeights, /* output */
        double& expectedLossOffset       /* output */) const;

	class MARKET_DLL TrancheProfile : public CObject {

		DateTime date;
		DoubleArray baseStrikes;
		DoubleArray baseStrikesWeights;
		double expectedLossOffset;

	public:

		virtual ~TrancheProfile() {};

		static CClassConstSP const TYPE;

		TrancheProfile(CClassConstSP clazz = TYPE);

		TrancheProfile(
			DateTime date,
			DoubleArray baseStrikes,
			DoubleArray baseStrikesWeights,
			double expectedLossOffset,
			CClassConstSP clazz = TYPE);

		static IObject* defaultTrancheProfile();

		static void load(CClassSP& clazz);
	};

	TrancheProfileArraySP getTrancheDecomposition() const;

};

DRLIB_END_NAMESPACE

#endif
