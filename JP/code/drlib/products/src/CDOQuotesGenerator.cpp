//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesGenerator.hpp
//
//   Description : CDOQuotesGenerator is a generator of pseudo CDO quotes 
//                  from base correlation model inputs
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CDOQuotesGenerator.hpp"
#include "edginc/BenchmarkDate.hpp"
#include "edginc/CDOQuotesBootstrapper.hpp"
#include "edginc/CDO.hpp"
#include "edginc/CreditMetricsBaseCorrelation.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE


/** Destructor */
CDOQuotesGenerator::~CDOQuotesGenerator() {}


/** Get the market data needed 
  * This function actually does all the job of building the instruments from the
  * base correlation skew surface, price the instruments using  and store the
  * quotes in the object as a field retrieved later.
  */
void CDOQuotesGenerator::getMarket(const IModel* model, const MarketData* market)
{
    try
    {
        market->GetReferenceDate(valueDate);

        // call model->getInstrumentAndModelMarket to initiate market data selection
        this->model->getMarket(market,IInstrumentCollectionSP(0));

        portfolio->getMarket(this->model.get(), market);
        trancheConventions->getMarket(this->model.get(), market);    

        CDoubleArraySP tempSpreads(new CDoubleArray(lowStrikes->size()));
        CDoubleArraySP tempUpfronts(new CDoubleArray(lowStrikes->size()));


        //Loop through the surface points to extract the points (maturity, K1, K2)
        int i = 0;
        for (i=0; i<lowStrikes->size(); ++i)
        {
            (*tempSpreads)[i]  = runningEquityCoupon;
            (*tempUpfronts)[i] = 0.0;
        }

        // build the cdo quotes 
        CDOQuotesSP dummyCdoQuotes = CDOQuotesSP(new CDOQuotes(
                                                    valueDate,
                                                    portfolio,
                                                    trancheConventions,
                                                    expiries,
                                                    lowStrikes,
                                                    highStrikes,
                                                    tempSpreads,
                                                    tempUpfronts));

        // we need to put quotes in the cdoQuotes field because a clone on it
        // is used below in "createWithNewQuotes"
        cdoQuotes = dummyCdoQuotes;

        // get a cdo quotes bootstrapper
        CDOQuotesBootstrapperSP quotesBootstrapper = CDOQuotesBootstrapperSP(new CDOQuotesBootstrapper(dummyCdoQuotes));

        // spreads and upfronts will contain the actual quotes computed below
        CDoubleArraySP spreads(new CDoubleArray());
        CDoubleArraySP upfronts(new CDoubleArray());

        // build the controls neededs to compute the actual quotes
        OutputRequestArraySP outputRequest(new OutputRequestArray(0));
        outputRequest->push_back(OutputRequestSP(new OutputRequest(OutputRequest::TRANCHE_RISKY_DURATION)));
        outputRequest->push_back(OutputRequestSP(new OutputRequest(OutputRequest::TRANCHE_IMPLIED_SPREAD)));
        
        CControlSP control = CControlSP(new Control(
            SensitivityArraySP(new SensitivityArray(0)),
            outputRequest,
            false,
            ""));

        //Loop through the quotes
        i = 0;
        for (quotesBootstrapper->init(); !quotesBootstrapper->end(); quotesBootstrapper->next())
        {
            // get instrument corresponding to the quote
            // the resulting instrument array is of size 1 always since we are
            // using a BOOTSTRAP_TIME_STRIKE type of bootstrapper
            CInstrumentArraySP tempInst = quotesBootstrapper->buildCurrentInstruments();
            CDO *currentCDOInst = DYNAMIC_CAST(CDO, ((*tempInst)[0]).get() );
            
            CMarketDataSP mkt(const_cast<MarketData*>(market));

            currentCDOInst->GetMarket(this->model.get(), mkt);

            CResultsSP results = CResultsSP(
                        this->model->Run(
			            currentCDOInst,
			            control.get()));

            double impliedSpreadVal =
                            CDoubleConstSP::dynamicCast(
                            results->retrieveRequestResult("TRANCHE_IMPLIED_SPREAD")
                            )->doubleValue();
            double riskyDurationVal =
                            CDoubleConstSP::dynamicCast(
                            results->retrieveRequestResult("TRANCHE_RISKY_DURATION")
                            )->doubleValue();

            // now build the actual CDO quotes
            if ((*lowStrikes)[i] == 0.0) // if this is an equity tranche
            {
                (*spreads).push_back(runningEquityCoupon);
                (*upfronts).push_back(riskyDurationVal * (impliedSpreadVal - runningEquityCoupon));
            }
            else
            {
                (*spreads).push_back(impliedSpreadVal);
                (*upfronts).push_back(0.0); 
            }
        }

        // build the actual quotes now
        cdoQuotes = cdoQuotes->createWithNewQuotes(
            expiries,
            lowStrikes,
            highStrikes,
            spreads,
            upfronts);
    }
    catch (exception& e)
    {
        throw ModelException(e, "CDOQuotes::getMarket");
    }
}

/** Get the name of the market object */
string CDOQuotesGenerator::getName() const
{
    return name;
}

/** Get the underlying built CDOQuotes */
CDOQuotesConstSP CDOQuotesGenerator::getCDOQuotes() const
{
    return cdoQuotes;
}


/** Called immediately after object constructed */
void CDOQuotesGenerator::validatePop2Object() {
    try {
    	//TODO
    } catch (exception& e){
        throw ModelException(e, "CDOQuotes::validatePop2Object");
    }
}

/** Access to portfolio */
CDOPortfolioConstSP CDOQuotesGenerator::getPortfolio() const {
    return portfolio;
} 

/** Access to value date */
const DateTime& CDOQuotesGenerator::getValueDate() const {
    return valueDate;
}

/** Main method : only retrives the CDO quotes built in GetMarket */
CDOQuotesConstSP CDOQuotesGenerator::buildCDOQuotes() const
{
    if (cdoQuotes.get())
    {
        return cdoQuotes;
    }
    else
    {
        throw ModelException("CDOQuotesGenerator::buildCDOQuotes",
            "internal error: CDO quotes have not been built");
    }
}

/** Invoked when Class is 'loaded' */
void CDOQuotesGenerator::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CDOQuotesGenerator, clazz);
    IMPLEMENTS(ICDOQuotesGenerator);
    SUPERCLASS(MarketObject);
    EMPTY_SHELL_METHOD(defaultConstructor);	
    FIELD(portfolio,
        "Portfolio (including single names betas and par spreads)");
    FIELD(name, "Name");
    FIELD(model, "base correlation model");

    FIELD(lowStrikes, "list of low strikes for the quotes");
    FIELD(highStrikes, "list of high strikes for the quotes");
    FIELD(expiries,   "list of expiries for the quotes");

    FIELD(engineParameters, "engine Parameters to be used with the model");

    FIELD(runningEquityCoupon, "running coupon used for the equity pricing");
    FIELD(trancheConventions, "Tranche conventions");
    FIELD(valueDate, "Value date");
    FIELD_MAKE_TRANSIENT(valueDate);
    FIELD(cdoQuotes, "internally built quotes");
    FIELD_MAKE_TRANSIENT(cdoQuotes);
}

/** Public constructor*/
CDOQuotesGenerator::CDOQuotesGenerator( string name,
                                        IModelSP model,
                                        DoubleArraySP lowStrikes,
                                        DoubleArraySP highStrikes,
                                        ExpiryArraySP expiries,
                                        CDOPortfolioSP  portfolio,
                                        CreditEngineParametersWrapper engineParameters, 
                                        ICreditLegConventionSP trancheConventions,
                                        double                 runningEquityCoupon) :
    MarketObject(TYPE),
    name(name),
    model(model),
    lowStrikes(lowStrikes),
    highStrikes(highStrikes),
    expiries(expiries),
    portfolio(portfolio),
    engineParameters(engineParameters),
    trancheConventions(trancheConventions),
    runningEquityCoupon(runningEquityCoupon)
{}

/** Private constructor (only build instances of that class using reflection) */
CDOQuotesGenerator::CDOQuotesGenerator() :
    MarketObject(TYPE),
    name(0),
    model(0),
    lowStrikes(0),
    highStrikes(0),
    expiries(0),
    portfolio(0),
    engineParameters(0),
    trancheConventions(0),
    runningEquityCoupon(0.05)
{}

/** Default constructor */
IObject* CDOQuotesGenerator::defaultConstructor() {
    return new CDOQuotesGenerator();
}

/** TYPE for CDOQuotesGenerator */
CClassConstSP const CDOQuotesGenerator::TYPE = CClass::registerClassLoadMethod(
    "CDOQuotesGenerator", typeid(CDOQuotesGenerator), CDOQuotesGenerator::load);

/* external symbol to allow class to be forced to be linked in */
bool CDOQuotesGeneratorLoad(){
    return (CDOQuotesGenerator::TYPE != 0);
}

DRLIB_END_NAMESPACE
