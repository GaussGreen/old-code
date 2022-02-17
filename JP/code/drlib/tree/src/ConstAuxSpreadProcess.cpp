//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ConstAuxSpreadProcess.cpp
//
//   Description : Constant (trivial) auxiliary spread process. Instance of IAuxSpreadProcess
//
//   Author      : Matthias Arnsdorf
//
//   Date        : July 2006
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Format.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/ConstAuxSpreadProcess.hpp"

DRLIB_BEGIN_NAMESPACE


/** public constructor */
ConstAuxSpreadProcess::ConstAuxSpreadProcess(
	double intensity
	) : CObject(TYPE), intensity(intensity), range()
{
	validatePop2Object();
}

/** private constructor */
ConstAuxSpreadProcess::ConstAuxSpreadProcess() : CObject(TYPE), intensity(0.01)
{
}

/** Destructor */
ConstAuxSpreadProcess::~ConstAuxSpreadProcess()
{}

/** validate */
void ConstAuxSpreadProcess::validatePop2Object()
{
	static const string method = "ConstAuxSpreadProcess::validatePop2Object";
	try
	{

		transProbs = DoubleArrayArrayArraySP(new DoubleArrayArrayArray(1));
		(*transProbs)[0].resize(1);
		(*transProbs)[0][0].resize(1);
		(*transProbs)[0][0][0] = 1.0;

		outerBot = IntArrayArraySP(new IntArrayArray(1));
		(*outerBot)[0] = IntArraySP(new IntArray(1));
		(*(*outerBot)[0])[0] = 0;
	}
	catch(exception & e)
	{
		throw ModelException(e, method);
	}
}

/** chance to pass market data to the spread process */
void ConstAuxSpreadProcess::getMarketData(const MarketData* market, int stepsPerYear) 
{
    // populate value date
    valueDate = market->GetReferenceDate(); 
   
    if(indexCurve.isEmpty())
    {
        fwdDates.resize(1);
        fwdDates[0] = valueDate;
        cleanSpreads = DoubleArraySP(new DoubleArray(1));
        (*cleanSpreads)[0] = intensity;
    }
    else
    {
         // Dummy model 
        CModelSP  dummy(new NonPricingModel());
        indexCurve.getData(dummy.get(), market);


        // intialise fwds -------------------------------------------------------------
        DefaultRatesSP defRates = indexCurve.getSP()->defaultRates();
        CashFlowArraySP cleanSpreadFlows = defRates->getCleanSpreadCurve();

        fwdDates = CashFlow::dates(*cleanSpreadFlows);
        cleanSpreads = CashFlow::amounts(*cleanSpreadFlows);

    }
    
} 



/** initialise the model. Called prior to any requests. */
void ConstAuxSpreadProcess::setupModel(
                        TimeLineSP timeLine		// timeline to use
                        ) 
{
    numDates = timeLine->NumOfStep+1;
    // set up range
    range = TreeSliceGeneral::Range::create( 1, 0, 0 );

    // check that getMarketData has been called
    if (!cleanSpreads)
    {
        throw ModelException("indexCurve has not been initialised. Has getMarketData() been called?");
    }


    // assume that the clean spread are FLAT_FORWARD (this is correct currently: 19.9.06)

    /** fwdIdx[t] gives index of first date int spreadDates >= stepDates[t] */ 
    vector<int> fwdIdx = DateTime::getCeilingProjection(timeLine->StepDates, fwdDates);

    fwds.resize(fwdIdx.size());
    int t = 0;
    for(; t < (int)fwdIdx.size(); t++)
    {
        fwds[t] = (*cleanSpreads)[fwdIdx[t]];
    }
}


/** get dates at which spread values change 
* can be called after getMarketData
*/
DateTimeArray ConstAuxSpreadProcess::getSpreadDates() const 
{
    return fwdDates;
}

/** return the spread value at current time step and for given index.
Can be called after update */
double ConstAuxSpreadProcess::spread(int index) const 
{
    return intensity;
}


void ConstAuxSpreadProcess::load (CClassSP& clazz) {
	clazz->setPublic(); // make visible to EAS/spreadsheet
	REGISTER(ConstAuxSpreadProcess, clazz);
	SUPERCLASS(CObject);
	IMPLEMENTS(IAuxSpreadProcess);
	EMPTY_SHELL_METHOD(defaultConstAuxSpreadProcess);

	FIELD(intensity, "intial intensity [default = 100bps]");
	FIELD_MAKE_OPTIONAL(intensity);

    FIELD(indexCurve, "Mean reversion level as term structure of index spreads");
    FIELD_MAKE_OPTIONAL(indexCurve);

	FIELD(valueDate,"");
	FIELD_MAKE_TRANSIENT(valueDate);

    FIELD(fwdDates,"");
    FIELD_MAKE_TRANSIENT(fwdDates);

    FIELD(cleanSpreads,"");
    FIELD_MAKE_TRANSIENT(cleanSpreads);

	FIELD(transProbs,"");
	FIELD_MAKE_TRANSIENT(transProbs);

	FIELD(outerBot,"");
	FIELD_MAKE_TRANSIENT(outerBot);


	


}

IObject* ConstAuxSpreadProcess::defaultConstAuxSpreadProcess() {
	return new ConstAuxSpreadProcess();
}

CClassConstSP const ConstAuxSpreadProcess::TYPE = 
CClass::registerClassLoadMethod("ConstAuxSpreadProcess", 
								typeid(ConstAuxSpreadProcess), 
								load);


/** Included in ProductsLib-modified::linkInClasses() via the productSrcsMap
* script to force the linker to include this file */
bool ConstAuxSpreadProcessLoad() {
	return (ConstAuxSpreadProcess::TYPE != 0);
}








DRLIB_END_NAMESPACE
