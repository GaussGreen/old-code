//----------------------------------------------------------------------------
//
//   Group       : CH Quantitative Research
//
//   File name   : ConditionalLossModel.hpp
//
//   Date        : 04-Oct-2006
//
//   Description : Interface for class ConditionalLossModel, which provides a 
//                 method whereby we use MC to simulate losses from a 
//                 lossConfig and a portfolio, and use the samples to build a 
//                 regressor loss config which we price using a second model
//                 
//----------------------------------------------------------------------------

#ifndef EDR_CONDITIONALLOSSMODEL_HPP
#define EDR_CONDITIONALLOSSMODEL_HPP


#include "edginc/Control.hpp"
#include "edginc/MonteCarlo.hpp"
#include "edginc/ModelConfigMapper.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/LastSensDate.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/ITimelineSpec.hpp"
#include "edginc/CreditTrancheLossConfig.hpp"
#include "edginc/StdDevForLossMapping.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE

DECLARE(MonteCarlo);

class MCARLO_DLL ConditionalLossModel :
					public CModel,
                    virtual public LastProductSensDate,
					virtual public ITweakableWithRespectTo<StdDevForLossMapping>
{

	//CInstrumentSP instrument;

	ICreditLossConfigSP genCDO; 

	CDOPortfolioSP regressor;

	int numBuckets;

	ITimelineSpecSP timelineSpec;

    DateTimeArrayConstSP timeline; // dates to sample losses on

	double alphaTweak;        // coefficient for stddev of mapping

public:

	static CClassConstSP const TYPE;

	ConditionalLossModel(CClassConstSP clazz = TYPE) : 
		CModel(clazz),
			numBuckets(100),
			alphaTweak(0.0)
		{};

	virtual ~ConditionalLossModel() {};

    DateTimeArrayConstSP getTimeLine() const;

	int numberBuckets() const {return numBuckets;};

     /** Implementation of method in IModel */
    virtual void Price(
		CInstrument*  instrument,
        Control*      control, 
        CResults*     results);

	ConditionalLossModel(
		IModel       *mcModel,
		IModelConfigMapperConstSP mapper);

	static void load(CClassSP& clazz);

	static IObject* defaultConditionalLossModel();

	virtual IModel::WantsRiskMapping wantsRiskMapping() const {
		return IModel::riskMappingIrrelevant;
	};

	virtual DateTime endDate(
		const CInstrument *instrument, 
		const Sensitivity *sensitivity) const;

	virtual void validatePop2Object();

	void getMarket(const MarketData*  market,
                           IInstrumentCollectionSP instruments)
	{
		mcModel->getMarket(market, instruments);
	};

protected:

	MonteCarloSP   mcModel;
	IModelSP       impliedModel;  // specifies second model

public:

	MonteCarlo* getMCModel() const
		{return mcModel.get();}; 

	IModel* getImpliedModel() const
		{return impliedModel.get();}; 

	class MCARLO_DLL IProduct {

	public:

		virtual void price(Control* control, Results* results) = 0;

	};

	DECLARE(IProduct)

public:

	// this is the class instruments derive from in order to implement 
	// the conditional loss method
	class IIntoProduct :  virtual public CModel::IModelIntoProduct {

	public:

		static CClassConstSP const TYPE;

		virtual ~IIntoProduct() {};

		virtual IProduct* createProduct(const ConditionalLossModel* model) const = 0;
	};

	virtual string sensName(const StdDevForLossMapping *param) const;

	virtual TweakOutcome sensShift(const PropertyTweak<StdDevForLossMapping>& shift);

};

DRLIB_END_NAMESPACE

#endif
