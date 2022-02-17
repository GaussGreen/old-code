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
#include "edginc/config.hpp"
#include "edginc/ConditionalLossModel.hpp"
#include "edginc/FlatCDO2LossConfig.hpp"
#include "edginc/Control.hpp"
#include "edginc/Results.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/SingleCreditAsset.hpp"
#include "edginc/IProtectionProvider.hpp"
#include "edginc/CtgLegLossPerDefault.hpp"
#include "edginc/FeeLegReductionPerDefault.hpp"
#include "edginc/IRebateCalculator.hpp"
#include "edginc/IEffectiveCurveLossGen.hpp"
#include "edginc/CmOnlyParameters.hpp"


DRLIB_BEGIN_NAMESPACE

bool ConditionalLossModelLoad() {return ConditionalLossModel::TYPE != NULL;};

CClassConstSP const ConditionalLossModel::TYPE = 
	CClass::registerClassLoadMethod(
		"ConditionalLossModel", 
		typeid(ConditionalLossModel), 
		ConditionalLossModel::load);

void ConditionalLossModel::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(ConditionalLossModel, clazz);
	SUPERCLASS(CModel);

	FIELD(mcModel, "Monte Carlo Model");	

	FIELD(impliedModel, "Model to use on loss map");	

	FIELD(timelineSpec, "Timeline specifier");
	FIELD_MAKE_OPTIONAL(timelineSpec)

    FIELD(timeline, "Timeline for loss profiles");	
	FIELD_MAKE_OPTIONAL(timeline);

	FIELD(numBuckets, "Number of points in profile, "
						"default = 100");

	FIELD(alphaTweak, "Multiple to std dev to shift mapping by");
	FIELD_MAKE_TRANSIENT(alphaTweak);

	EMPTY_SHELL_METHOD(defaultConditionalLossModel);
}

IObject* ConditionalLossModel::defaultConditionalLossModel() 
{	
	return new ConditionalLossModel(TYPE);
};

void ConditionalLossModel::validatePop2Object()
{

	// DateTimeArraySP dates = DateTimeArraySP(new DateTimeArray());

	if (numBuckets < 1)
		throw ModelException ("Must have at least 1 bucket");

	int i;

	if ( !(timeline.get()) && !(timelineSpec.get()))
	{
		throw ModelException("Need to provide timeline or timeline specifier"
							 " to conditional loss method");
	};

	if ( !(timeline->size()) && !(timelineSpec.get()))
	{
		throw ModelException("Empty timeline passed" 
							 " to conditional loss method");
	};

	for (i = 1; i < timeline->size(); i++)
	{
		if (i>0)
			if ((*timeline)[i] <= (*timeline)[i-1])
				throw ModelException("Dates must be increasing");
	};

};

DateTime ConditionalLossModel::endDate(
		const CInstrument *instrument, 
		const Sensitivity *sensitivity) const
{
	return mcModel->endDate(instrument, sensitivity);
};

DateTimeArrayConstSP
ConditionalLossModel::getTimeLine() const 
{ 
	return timeline;
};

void ConditionalLossModel::Price(
	CInstrument*  instrument,
    Control*      control, 
    CResults*     results)
{

	const IIntoProduct* intoProd = 
            dynamic_cast<const IIntoProduct*>(instrument);

	if (!intoProd)
		throw ModelException(
			"ConditionalLossModel::Price",
			"Method not valid for this instrument");

	if (timelineSpec.get())
		timeline = timelineSpec->timeline(*instrument); 
	else
	{
		DateTimeArraySP dates = DateTimeArraySP(new DateTimeArray());
		
		dates->push_back(getValueDate());

		int i;
		for (i = 0; i < timeline->size(); i++)
		{
			if ((*timeline)[i] > getValueDate())
			{
				dates->push_back((*timeline)[i]);
			};
		};

		timeline = dates; 
	};

	IProduct *prod = intoProd->createProduct(this);

	prod->price(control, results);
};

string 
ConditionalLossModel::sensName(const StdDevForLossMapping *param) const 
{
	return "StdDevShiftForLossMapping";
};

TweakOutcome 
ConditionalLossModel::sensShift(const PropertyTweak<StdDevForLossMapping>& shift)
{
	alphaTweak = shift.coefficient;
	return TweakOutcome(shift.coefficient, false);
};

DRLIB_END_NAMESPACE
