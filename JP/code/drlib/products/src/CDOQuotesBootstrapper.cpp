//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesBootstrapper.cpp
//
//   Description : CDOQuotesBootstrapper is a kind of "iterator" over CDOQuotes
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/CDOQuotesBootstrapper.hpp"
#include "edginc/CDO.hpp"
#include "edginc/CounterPartyCredit.hpp"
#include "edginc/TimePoint2D.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

/** Virtual destructor */
CDOQuotesBootstrapper::~CDOQuotesBootstrapper() {}

/** Constructor */
CDOQuotesBootstrapper::CDOQuotesBootstrapper(CDOQuotesConstSP cdoQuotes, bool ignore100pc) :
    CObject(TYPE),
    cdoQuotes(cdoQuotes),
    ignore100pc(ignore100pc),
    maturities(0),
    lowStrikes(0),
    highStrikes(0),
    currentPointIdx(0),
    portfolioNotional(0.0),
    sharedPortfolio(0),
	instruments(0)
{
    int trancheIdx, maturityIdx;
    double lowStrike, highStrike;
    DateTimeArrayConstSP trancheMaturities;
    CDOTrancheQuotesConstSP trancheQuote;
    
    // initialise portfolioNotional
    CDOPortfolioConstSP portfolio = cdoQuotes->getPortfolio();
    for (int i = 0; i < portfolio->nbName(); ++i) {
        portfolioNotional += portfolio->getName(i)->getNameNotional();
    }

    // initialise expiries, lowStrikes and highStrikes in 2 steps:
    // step 1 - split the CDOQuotes structure into an ordered set
    // step 2 - flatten the set and build expiries, lowStrikes and highStrikes
        
    CDOTrancheQuotesArrayConstSP cdoTrancheQuotes = cdoQuotes->getTrancheQuotes();
    
    // step 1 - split the CDOQuotes structure into an ordered set
    set<TimePoint2D, TimePoint2D::CompareDateC1C2> trancheQuotePointSet;
    
    for (trancheIdx = 0; trancheIdx < cdoTrancheQuotes->size(); ++trancheIdx) {
        trancheQuote = (*cdoTrancheQuotes)[trancheIdx];
        trancheMaturities = trancheQuote->getMaturityDates();
        lowStrike = trancheQuote->getLowStrike();
        highStrike = trancheQuote->getHighStrike();
        
        for(maturityIdx = 0; maturityIdx < trancheMaturities->size(); ++maturityIdx) 
        {
            // ignore points with highStrike >= 1.0 if ignore100pc = true
            if (!ignore100pc || highStrike < 1.0) 
            {
                trancheQuotePointSet.insert(TimePoint2D(
                    (*trancheMaturities)[maturityIdx],
                    lowStrike,
                    highStrike));
            }
        }
    }

    // step 2 - flatten the set and build expiries, lowStrikes and highStrikes
    DateTimeArraySP maturitiesTmp(new DateTimeArray());
    lowStrikes = DoubleArraySP(new DoubleArray());
    highStrikes = DoubleArraySP(new DoubleArray());
    set<TimePoint2D, TimePoint2D::CompareDateC1C2>::iterator iter = trancheQuotePointSet.begin();
    for(;iter != trancheQuotePointSet.end(); iter++) {
        lowStrikes->push_back(iter->getCoord1());
        highStrikes->push_back(iter->getCoord2());
        maturitiesTmp->push_back(iter->getDate());
    }
    maturities = maturitiesTmp;


	// set up cache of instruments //////////////////////////
	instruments = CInstrumentArraySP(new CInstrumentArray(maturities->size()));
	for (init(); !end(); next()) 
	{
		(*instruments)[currentPointIdx] = buildCurrentInstrument();
	}
	
	// reset index
	currentPointIdx = 0;
}

/** Constructor (internal - used by reflection) */
CDOQuotesBootstrapper::CDOQuotesBootstrapper() :
    CObject(TYPE),
    cdoQuotes(0),
    ignore100pc(true),
    maturities(0),
    lowStrikes(0),
    highStrikes(0),
    currentPointIdx(0),
    portfolioNotional(0.0),
    sharedPortfolio(0) {}

/**
 * Method called before first step of the loop
 * [Implements IBootstrapper]
 * */
void CDOQuotesBootstrapper::init() {
    currentPointIdx = 0;
    sharedPortfolio = CDOPortfolioSP::dynamicCast(
        IObjectSP(cdoQuotes->getPortfolio()->clone()));
	
}

/**
 * Method called after each step of the loop (ie go to the next quote)
 * [Implements IBootstrapper]
 * */
void CDOQuotesBootstrapper::next() {
    currentPointIdx++;
}

/**
 * Method called to test the end of the loop (returns true when there is no more quote)
 * [Implements IBootstrapper]
 * */
bool CDOQuotesBootstrapper::end() const 
{
    return (currentPointIdx >= maturities->size());
}

/**
 * Returns "state" corresponding to current step of the loop
 * [Implements IBootstrapper]
 * */
IObjectSP CDOQuotesBootstrapper::getCurrentState() const 
{
#ifdef DEBUG
    QLIB_VERIFY(currentPointIdx < maturities->size(), "currentPointIdx out of bounds");
#endif
    double highStrikePercent = (*highStrikes)[currentPointIdx];
    double lowStrikePercent = (*lowStrikes)[currentPointIdx];

    return TimePoint2DSP(new TimePoint2D(
        (*maturities)[currentPointIdx], lowStrikePercent, highStrikePercent));
}

/** Returns the instrument corresponding to the current quote */
CInstrumentSP CDOQuotesBootstrapper::buildCurrentInstrument() {
    static string method = "CDOQuotesBootstrapper::buildCurrentInstrument";
    try {
        if (end()) {
            // no more quotes
            throw ModelException(method,
                "Internal error: no more quotes to build instrument.");
        } else {
            // generate CDO components from current quote
            double lowStrike = (*lowStrikes)[currentPointIdx];
            double highStrike = (*highStrikes)[currentPointIdx];
            string discountName = cdoQuotes->getDiscountName();          
            ICreditFeeLegSP fLeg = cdoQuotes->generateFeeLeg(
                lowStrike,
                highStrike,
                (*maturities)[currentPointIdx]);
            ICreditContingentLegSP cLeg = cdoQuotes->generateContingentLeg(
                lowStrike,
                highStrike,
                (*maturities)[currentPointIdx]);
            
            return CInstrumentSP(new CDO(
                lowStrike * portfolioNotional,
                highStrike * portfolioNotional,
                true,
                CBoolSP(CBool::create(cdoQuotes->getRecoverNotional())),
                smartPtr<CounterPartyCredit>(0),
                cLeg,
                fLeg,
                sharedPortfolio,
                discountName));
        }
    } catch (exception& e){
        throw ModelException(e, "CDOQuotesBootstrapper::buildCurrentInstrument");
    }
}
/** Returns the instrument corresponding to the current quote */
CInstrumentArraySP CDOQuotesBootstrapper::buildCurrentInstruments() {
	static string method = "CDOQuotesBootstrapper::buildCurrentInstruments";
	try {
			if(!instruments)
			{
				throw ModelException("instruments are not initialised");
			}
			if(currentPointIdx >= (int)instruments->size())
			{
				throw ModelException("currentPointIdx out of bounds");
			}
			CInstrumentArraySP outArray(new CInstrumentArray(1));
			(*outArray)[0] = (*instruments)[currentPointIdx];
			
			return outArray;

	} catch (exception& e){
		throw ModelException(e, "CDOQuotesBootstrapper::buildCurrentInstruments");
	}
}

/** Invoked when Class is 'loaded' */
void CDOQuotesBootstrapper::load(CClassSP& clazz) {
    clazz->setPrivate(); // don't make visible to EAS/spreadsheet
    REGISTER(CDOQuotesBootstrapper, clazz);
    SUPERCLASS(CObject);
    //IMPLEMENTS(IBootstrapper);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD_NO_DESC(cdoQuotes);
    FIELD_NO_DESC(ignore100pc);
    FIELD_NO_DESC(maturities);
    FIELD_NO_DESC(lowStrikes);
    FIELD_NO_DESC(highStrikes);
    FIELD_NO_DESC(currentPointIdx);
    FIELD_NO_DESC(portfolioNotional);
    FIELD_NO_DESC(sharedPortfolio);
	FIELD_NO_DESC(instruments);
}

/** Default constructor */
IObject* CDOQuotesBootstrapper::defaultConstructor() {
    return new CDOQuotesBootstrapper();
}

/** TYPE (for reflection) */
CClassConstSP const CDOQuotesBootstrapper::TYPE =
    CClass::registerClassLoadMethod(
        "CDOQuotesBootstrapper",
        typeid(CDOQuotesBootstrapper),
        CDOQuotesBootstrapper::load);

DRLIB_END_NAMESPACE

