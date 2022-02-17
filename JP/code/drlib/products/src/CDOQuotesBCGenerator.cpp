//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesBCGenerator.hpp
//
//   Description : CDOQuotesBCGenerator is a generator of pseudo CDO quotes 
//                  from base correlation model inputs
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDOQuotesBCGenerator.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/CDOQuotesBootstrapper.hpp"
#include "edginc/BaseCorrelationOnlyParameters.hpp"
#include "edginc/CDO.hpp"
#include "edginc/CreditMetricsBaseCorrelation.hpp"
#include "edginc/SkewSurfaceBootstrapper.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE


/** Destructor */
CDOQuotesBCGenerator::~CDOQuotesBCGenerator() {}


/** Get the market data needed 
  * This function actually does all the job of building the instruments from the
  * base correlation skew surface, price the instruments using BC and store the
  * quotes in the object as a field retrieved later.
  */
void CDOQuotesBCGenerator::getMarket(const IModel* model, const MarketData* market)
{
    try
    {
        market->GetReferenceDate(valueDate);

        // need to build a base Correlation model if it is not present
        if (!bcModel.get())
        {
            CreditMetricsBaseCorrelationSP bc = CreditMetricsBaseCorrelationSP(new CreditMetricsBaseCorrelation());
            bc->validatePop2Object();
            IObjectSP tmpModel = bc;
            bc->convert(tmpModel, ConvolutionEngine::TYPE);

            bcModel = IModelSP(dynamic_cast<ConvolutionEngine*>(tmpModel.get()));
        }

        // call model->getInstrumentAndModelMarket to initiate market data selection
        bcModel->getMarket(market,IInstrumentCollectionSP(0));
        

        indexSkew.getData(bcModel.get(), market);
        indexSkew->getMarket(bcModel.get(), market);

        portfolio->getMarket(bcModel.get(), market);
        trancheConventions->getMarket(bcModel.get(), market);    


        // Creates field to calibrate so that it can be passed below
        // It does not matter at this point whether we are using skews or skewFast
        // since we are not actually doing a calibration.
        // The output will be used only to create the bootstrapper and the
        // corresponding instruments
        CFieldConstSP fieldToCalibrate = SkewSurface::TYPE->getDeclaredField("skews");
     
#if 0 // DISABLED (MA) IInstanceIDBootstrapperSP has been removed
        // First get a idxSkewBootstrapper from the idxSkewSurface
        IInstanceIDBootstrapperSP idxSkewBootstrapper = 
            indexSkew->createBootstrapperFromSkewSurface(
            fieldToCalibrate, 
            SkewSurfaceBootstrapper::BOOTSTRAP_TIME_STRIKE);
#endif
        ExpiryArraySP  expiries(new ExpiryArray());
        CDoubleArraySP lowStrikes(new CDoubleArray());
        CDoubleArraySP highStrikes(new CDoubleArray());
        CDoubleArraySP tempSpreads(new CDoubleArray());
        CDoubleArraySP tempUpfronts(new CDoubleArray());

#if 0
        //Loop through the surface points to extract the points (maturity, K1, K2)
        for (idxSkewBootstrapper->init(); !idxSkewBootstrapper->end(); idxSkewBootstrapper->next())
        {
            TimePoint2D *tp = DYNAMIC_CAST(TimePoint2D, ((idxSkewBootstrapper->getCurrentState()).get()));
            (*expiries).push_back(ExpirySP(new BenchmarkDate(tp->getDate())));
            (*lowStrikes).push_back(tp->getCoord1());
            (*highStrikes).push_back(tp->getCoord2());
            (*tempSpreads).push_back(runningEquityCoupon);
            (*tempUpfronts).push_back(0.0);
        }
#endif

        // need to build the base correlation engine parameters from the index skew
        IndexSkewWrapperArraySP idxSkewWrapperArray(new IndexSkewWrapperArray(1));
        DoubleArraySP           weightsArray(new DoubleArray(1));

        (*idxSkewWrapperArray)[0] = indexSkew;
        (*weightsArray)[0]        = 1.0;
        IndexWeightsSP indexWeights = IndexWeightsSP(new IndexWeights(idxSkewWrapperArray, weightsArray));

        CreditEngineParametersWrapper engineParameters = CreditEngineParametersWrapper(
                            new BaseCorrelationOnlyParameters(indexSkew->getName(), indexWeights));

        // build the generic quotes generator
        quotesGenerator = CDOQuotesGeneratorSP( new CDOQuotesGenerator(
                                                                name,
                                                                bcModel,
                                                                lowStrikes,
                                                                highStrikes,
                                                                expiries,
                                                                portfolio,
                                                                engineParameters,
                                                                trancheConventions,
                                                                runningEquityCoupon));

        // call the getMarket of the generic Quotes Generaor
        quotesGenerator->getMarket(model, market);

        // populates the field with the corresponding field in the generic quotes Generator
        cdoQuotes = quotesGenerator->getCDOQuotes();
    }
    catch (exception& e)
    {
        throw ModelException(e, "CDOQuotes::getMarket");
    }
}

/** Get the name of the market object */
string CDOQuotesBCGenerator::getName() const
{
    return name;
}


/** Called immediately after object constructed */
void CDOQuotesBCGenerator::validatePop2Object() {
    try {
    	//TODO
    } catch (exception& e){
        throw ModelException(e, "CDOQuotes::validatePop2Object");
    }
}

/** Access to portfolio */
CDOPortfolioConstSP CDOQuotesBCGenerator::getPortfolio() const {
    return portfolio;
} 

/** Access to value date */
const DateTime& CDOQuotesBCGenerator::getValueDate() const {
    return valueDate;
}

/** Main method : only retrives the CDO quotes built in GetMarket */
CDOQuotesConstSP CDOQuotesBCGenerator::buildCDOQuotes() const
{
    if (cdoQuotes.get())
    {
        return cdoQuotes;
    }
    else
    {
        throw ModelException("CDOQuotesBCGenerator::buildCDOQuotes",
            "internal error: CDO quotes have not been built");
    }
}

/** Invoked when Class is 'loaded' */
void CDOQuotesBCGenerator::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDOQuotesBCGenerator, clazz);
    IMPLEMENTS(ICDOQuotesGenerator);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultConstructor);	
    FIELD(portfolio,
        "Portfolio (including single names betas and par spreads)");
    FIELD(name, "Name");
    FIELD(indexSkew, "base correlation index skew");
    FIELD(bcModel, "base correlation model");
    FIELD_MAKE_OPTIONAL(bcModel);
    FIELD(runningEquityCoupon, "running coupon used for the equity pricing");
    FIELD(trancheConventions, "Tranche conventions");
    FIELD(valueDate, "Value date");
    FIELD_MAKE_TRANSIENT(valueDate);
    FIELD(cdoQuotes, "internally built quotes");
    FIELD_MAKE_TRANSIENT(cdoQuotes);
    FIELD(quotesGenerator, "generic quotes generator");
    FIELD_MAKE_TRANSIENT(quotesGenerator);
}

/** Private constructor (only build instances of that class using reflection) */
CDOQuotesBCGenerator::CDOQuotesBCGenerator() :
    MarketObject(TYPE),
    name(""),
    portfolio(0),
    trancheConventions(0),
    bcModel(0),
    indexSkew(0),
    runningEquityCoupon(0.05)
{}

/** Default constructor */
IObject* CDOQuotesBCGenerator::defaultConstructor() {
    return new CDOQuotesBCGenerator();
}

/** TYPE for CDOQuotesBCGenerator */
CClassConstSP const CDOQuotesBCGenerator::TYPE = CClass::registerClassLoadMethod(
    "CDOQuotesBCGenerator", typeid(CDOQuotesBCGenerator), CDOQuotesBCGenerator::load);
    
/* external symbol to allow class to be forced to be linked in */
bool CDOQuotesBCGeneratorLoad(){
    return (CDOQuotesBCGenerator::TYPE != 0);
}

DRLIB_END_NAMESPACE
