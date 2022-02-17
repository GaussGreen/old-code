//----------------------------------------------------------------------------
//
//   Group       : Quatitative Research
//
//   Filename    : ImpliedLossModel.cpp
//
//   Description : Gives index prices of off-market tranches by interpolating set of CDOQuotes
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/CCMPriceUtil.hpp"
#include "edginc/CreditFeeLegWithPV.hpp"
#include "edginc/ImpliedLossModel.hpp"

DRLIB_BEGIN_NAMESPACE


void ImpliedLossModel::load (CClassSP& clazz) {
    clazz->setPublic();			// make visible to EAS/spreadsheet
    REGISTER(ImpliedLossModel, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultImpliedLossModel);

	FIELD(trancheQuoteInterpolator, "Model to interpolate tranche quote term structure");

	FIELD(expLossInterpolator, "Model to interpolate expected losses");
	FIELD_MAKE_OPTIONAL(expLossInterpolator);

}

IObject* ImpliedLossModel::defaultImpliedLossModel() {
    return new ImpliedLossModel();
}

CClassConstSP const ImpliedLossModel::TYPE = 
    CClass::registerClassLoadMethod("ImpliedLossModel", 
                                    typeid(ImpliedLossModel), 
                                    load);

/** private constructor */
ImpliedLossModel::ImpliedLossModel() : CObject(TYPE) , trancheQuoteInterpolator(0), expLossInterpolator(0)
{}



/** Destructor */
ImpliedLossModel::~ImpliedLossModel()
{}


const string ImpliedLossModel::SPREAD ="SPREAD";		/** par spread */
const string ImpliedLossModel::ANNUITY = "ANNUITY";		/** anuuity / duration */
const string ImpliedLossModel::CONT_LEG = "CONT_LEG";	/** contingent leg */
const string ImpliedLossModel::UPFRONT = "UPFRONT";		/** upfront */
const string ImpliedLossModel::LOSS = "LOSS";			/** loss */
const string ImpliedLossModel::RISKY_ZERO = "RISKY_ZERO"; /** riksy zero */

// default coupon at 500 bps
const double ImpliedLossModel::DEFAULT_COUPON = 0.05;

 /** return array of outputType. No coupons*/
DoubleArraySP ImpliedLossModel::getValues(
		const CDOQuotes & marketQuotes,		/** set of market tranche qutoes */
		const ICDSParSpreads & index,		/** index swap spreads			*/
		const ExpiryArray & expiries,		/** set of expiries for output  */
		const DoubleArray & lowStrikes,		/** set of low strikes for output */
		const DoubleArray & highStrikes,	/** set of high sitrkes for output  */
		const StringArray &	outputType,     /** output type */
        IForwardRatePricerSP model          /** for calculating fees */
		) const
{
	

	DoubleArray coupons(outputType.size());
	for(int i = 0;i< outputType.size();i++)
	{
		coupons[i] = DEFAULT_COUPON;
	}

	return getValues(
		marketQuotes,
		index,
		expiries,
		lowStrikes,
		highStrikes,
		outputType,
		coupons,
        model
		);
}


/** return array of outputType */
DoubleArraySP ImpliedLossModel::getValues(
		const CDOQuotes & marketQuotes,		/** set of market tranche qutoes */
		const ICDSParSpreads & index,		/** index swap spreads			*/
		const ExpiryArray & expiries,		/** set of expiries for output  */
		const DoubleArray & lowStrikes,		/** set of low strikes for output */
		const DoubleArray & highStrikes,	/** set of high sitrkes for output  */
		const StringArray &	outputType,		/** output type */
		const DoubleArray &	coupons,		/** coupons */
        IForwardRatePricerSP model          /** for calculating fees */
		) const
{
	static const string method = "ImpliedLossModel::getValues";
	try
	{
		int i;
		/** size of input array */
		int numQuotes = expiries.size();

		// check that lowStrikes and highStrikes have same size ----------------------------------------
		if(lowStrikes.size() != numQuotes)
		{
			throw ModelException("lowStrikes size ("+Format::toString(lowStrikes.size())+
				") is not equal to numQuotes ("+Format::toString(numQuotes)+")");
		}
		if(highStrikes.size() != numQuotes)
		{
			throw ModelException("highStrikes size ("+Format::toString(highStrikes.size())+
				") is not equal to numQuotes ("+Format::toString(numQuotes)+")");
		}


		// get losses from stage A (time interpolation) ---------------------------------------------------
		ExpectedLossSurfaceSP lossSurface  = trancheQuoteInterpolator->getELSurface(marketQuotes, index);

		// do strike interpolation if LossInterp model is available
		if(expLossInterpolator.get() != NULL)
		{
			lossSurface = expLossInterpolator->getELSurface(*lossSurface);
		}

		
			
		// varaibles needed later ----------------------------------------------------------------
		/**  ouptut array */
		DoubleArraySP out = DoubleArraySP(new DoubleArray(numQuotes));

		/** value date */
		DateTime valueDate = marketQuotes.getValueDate();

		/** dates for effective curve */
		DateTimeArrayConstSP dates = lossSurface->getDates();

		/** yearFracs (Act/365) corresponding to dates */
		DoubleArray effCurveTimes(dates->size());
		for(i = 0; i<dates->size();i++)
		{
			effCurveTimes[i] =  valueDate.yearFrac((*dates)[i]);
		}
		
		/** yield curve */
		YieldCurveConstSP yieldCurveSP = marketQuotes.getDiscount();

		/** times (Act.365 year fracs) at which have zero rates */
		DoubleArray yearFracsDf(0);
		/** discount factors with associated with yearFracsDf */
		DoubleArray dfCurve(0);

		// populate
		for(i = 0;i<yieldCurveSP->zeroDates().size();i++) 
		{
			dfCurve.push_back(yieldCurveSP->pv(yieldCurveSP->zeroDates()[i]));
			yearFracsDf.push_back(valueDate.yearFrac(yieldCurveSP->zeroDates()[i]));
		}


		// loop through quotes --------------------------------------------------------------------
		for(i = 0;i<numQuotes;i++)
		{
			if(outputType[i] == LOSS)
			{
				double highLoss = lossSurface->getBaseEL(highStrikes[i], expiries[i]->toDate(valueDate));
				double lowLoss = lossSurface->getBaseEL(lowStrikes[i], expiries[i]->toDate(valueDate));
				(*out)[i] = highLoss - lowLoss;

			}
			else if(outputType[i] == RISKY_ZERO)
			{
				if(highStrikes[i] == lowStrikes[i])
				{
					throw ModelException("highStrikes = lowStrikes at index: "+Format::toString(i));
				}
				double highLoss = lossSurface->getBaseEL(highStrikes[i], expiries[i]->toDate(valueDate));
				double lowLoss = lossSurface->getBaseEL(lowStrikes[i], expiries[i]->toDate(valueDate));
				double expNotional = 1 - (highLoss - lowLoss)/(highStrikes[i]-lowStrikes[i]);

				(*out)[i] = yieldCurveSP->pv(expiries[i]->toDate(valueDate))*expNotional; 

			}
			else
			{
				/** effective curve for fee leg. Assumes loss of recovered notional */
				DoubleArraySP feeEffCurve;
				feeEffCurve  = lossSurface->getFeeLegEffectiveCurve(lowStrikes[i],highStrikes[i]);

				/** effective curve for contingent leg. Assumes standard notional write down */
				DoubleArraySP contEffCurve;
				contEffCurve  = lossSurface->getEffectiveCurve(lowStrikes[i],highStrikes[i]);
				
				/** maturity */		
				DateTime endDate = expiries[i]->toDate(valueDate);
				
				/** generic fee leg with unit spread and 0 upfront */
				CreditFeeLegWithPVSP indexFeeLeg = 
					CreditFeeLegWithPVSP::dynamicCast(marketQuotes.generateFeeLegOverride(1,0,endDate));

				double feeLeg =  indexFeeLeg->pvRisky((*dates), 
								  (*feeEffCurve), 
								  valueDate,
								  yieldCurveSP,
                                  model);
				
				
				// TODO build contingent leg from quotes with correct conventions
				double protectionStart = valueDate.yearFrac(valueDate.rollDate(1));

				double contLeg = CCMPriceUtil::contingentPrice(
												protectionStart, // protection start
												valueDate.yearFrac(endDate), //protection end
												0, //delay
											   yearFracsDf, dfCurve,
											   effCurveTimes, (*contEffCurve),
											   CCMPriceUtil::LINEAR);
			
				// add result to out array
				if(outputType[i] == SPREAD)
				{
					if(feeLeg == 0) throw ModelException("feeLeg[" + Format::toString(i) +"] = 0. Par spread is not defined.");
					(*out)[i] = contLeg/feeLeg;
				}
				else if (outputType[i] == ANNUITY) (*out)[i] = feeLeg;
				else if (outputType[i] == CONT_LEG) (*out)[i] = contLeg;
				else if (outputType[i] == UPFRONT)	(*out)[i] = contLeg - coupons[i]*feeLeg;
				else
				{
					throw ModelException("outputType "+outputType[i]+" is not supported");
				}
			}
		} // end quote loop ---------------------------------------------------------------

		return out;
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}




/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool ImpliedLossModelLoad() {
    return (ImpliedLossModel::TYPE != 0);
}


DRLIB_END_NAMESPACE
