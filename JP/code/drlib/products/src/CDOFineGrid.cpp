//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOFineGrid.cpp
//
//   Description :  concrete CDO fine grid class responsible for
//                 "extending" market quotes
//
//   Author      : Matthias Arnsdorf
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/CDOFineGrid.hpp"
#include "edginc/ClosedFormForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE
 
// default coupon value
const double CDOFineGrid::DEFAULT_COUPON = 0.05;


void CDOFineGrid::load (CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDOFineGrid, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ICDOFineGrid);
    EMPTY_SHELL_METHOD(defaultCDOFineGrid);

	FIELD(expiries, "Expiries for fine grid quotes");
	
	FIELD(lowStrikes, "Low strikes for fine grid quotes");

	FIELD(highStrikes, "High strikes for fine grid quotes");

	FIELD(outputTypes, "Output types for fine grid quote ("+
		ImpliedLossModel::SPREAD+" or "+
		ImpliedLossModel::UPFRONT+"). [default = "+
		ImpliedLossModel::SPREAD+"]");
	FIELD_MAKE_OPTIONAL(outputTypes);

	FIELD(coupons, "Running coupon (absolute fraction) if outputType is UPFRONT otherwise ignored. [default = "+
		Format::toString(DEFAULT_COUPON)+"]");
	FIELD_MAKE_OPTIONAL(coupons);

	FIELD(indexSpreads, "Index swap spread wrapper");

	FIELD(model,"Interpolation model");

}

IObject* CDOFineGrid::defaultCDOFineGrid() {
    return new CDOFineGrid();
}

CClassConstSP const CDOFineGrid::TYPE = 
    CClass::registerClassLoadMethod("CDOFineGrid", 
                                    typeid(CDOFineGrid), 
                                    load);
/** private constructor */
CDOFineGrid::CDOFineGrid() : CObject(TYPE) , 
expiries(0),
lowStrikes(0),
highStrikes(0),
outputTypes(0),
coupons(0),
model(0)
{}


/** Destructor */
CDOFineGrid::~CDOFineGrid()
{}



/** validate */
void CDOFineGrid::validatePop2Object()
{
	static const string method = "CDOFineGrid::validatePop2Object";
	try
	{
		int i;
		/** number of quotes */
		int numQuotes = expiries->size();

		// check that all arrays have same size
		if(lowStrikes->size() != numQuotes)
		{
			throw ModelException("lowStrikes size ("+Format::toString(lowStrikes->size())+
				") is not equal to number of expiries ("+Format::toString(numQuotes)+")");
		}
		if(highStrikes->size() != numQuotes)
		{
			throw ModelException("highStrikes size ("+Format::toString(highStrikes->size())+
				") is not equal to number of expiries ("+Format::toString(numQuotes)+")");
		}

		// check if have outputTypes otherwise set default values
		if(outputTypes.get() != NULL)
		{
			if(outputTypes->size() != numQuotes)
			{
				throw ModelException("outputTypes size ("+Format::toString(outputTypes->size())+
					") is not equal to number of expiries ("+Format::toString(numQuotes)+")");
			}
		}
		else
		{
			outputTypes = StringArraySP(new StringArray(numQuotes));
			for(i = 0;i<numQuotes;i++)
			{
				(*outputTypes)[i] = ImpliedLossModel::SPREAD;
			}
		}

		// check if have coupons otherwise set default values
		if(coupons.get() != NULL)
		{
			if(coupons->size() != numQuotes)
			{
				throw ModelException("coupons size ("+Format::toString(coupons->size())+
					") is not equal to number of expiries ("+Format::toString(numQuotes)+")");
			}
		}
		else
		{
			coupons = DoubleArraySP(new DoubleArray(numQuotes));
			for(i = 0;i<numQuotes;i++)
			{
				(*coupons)[i] = DEFAULT_COUPON;
			}
		}


	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
/** get market */
////////////////////////////////////////////////////////////////////////////////////////
void CDOFineGrid::getMarket(MarketData *market)
{
		// Dummy model (note that ImpliedLossModel does not derive off CModel yet)
		CModelSP  dummy(new NonPricingModel());

		// get market data -----------------------------------------
		indexSpreads.getData(dummy.get(), market);

}

////////////////////////////////////////////////////////////////////////////////////////////////
/**
Extend quotes
*/
///////////////////////////////////////////////////////////////////////////////////////////////
CDOQuotesSP CDOFineGrid::extendQuotes(CDOQuotesConstSP marketQuotes)
{
	static const string method = "CDOFineGrid::extendQuotes";
	try
	{
		int i;
		int N = highStrikes->size();
		
        // Require a model for calculating fees
        // so just provide the default closed form
        ClosedFormForwardRatePricerSP cfPricer =
            ClosedFormForwardRatePricerSP(
                new ClosedFormForwardRatePricer());

		/** model outputs */
		DoubleArraySP outputs = model->getValues(
			*marketQuotes,
			*indexSpreads.getSP(),
			*expiries,
			*lowStrikes,
			*highStrikes,
			*outputTypes,
			*coupons,
			cfPricer);

		DoubleArraySP fineGridSpreads = DoubleArraySP(new DoubleArray(N));
		DoubleArraySP fineGridUpfronts = DoubleArraySP(new DoubleArray(N));

		for(i = 0;i<N;i++)
		{
			if((*outputTypes)[i] == ImpliedLossModel::SPREAD)
			{
				(*fineGridSpreads)[i] = (*outputs)[i];
				(*fineGridUpfronts)[i] = 0.0;
			}
			else if((*outputTypes)[i] == ImpliedLossModel::UPFRONT)
			{
				(*fineGridSpreads)[i] = (*coupons)[i];
				(*fineGridUpfronts)[i] = (*outputs)[i];
			}
			else
			{
				throw ModelException("OutputType "+(*outputTypes)[i]+ " is not supported. Allowed values are 'SPREAD' or 'UPFRONT'");
			}
		}
		
		return marketQuotes->createWithNewQuotes(
			expiries,
			lowStrikes,
			highStrikes,
			fineGridSpreads,
			fineGridUpfronts);
	
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}



/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
 * script to force the linker to include this file */
bool CDOFineGridLoad() {
    return (CDOFineGrid::TYPE != 0);
}


DRLIB_END_NAMESPACE
