//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOTrancheQuotes.hpp
//
//   Description : CDOTrancheQuotes is a market data container for CDO tranche quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------


#ifndef CDO_TRANCHE_QUOTES_HPP
#define CDO_TRANCHE_QUOTES_HPP

#include "edginc/YieldCurve.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

// Forward declaration of CDOTrancheQuotes
class CDOTrancheQuotes;

// Support for smart pointers
typedef smartPtr<CDOTrancheQuotes>      CDOTrancheQuotesSP;
typedef smartConstPtr<CDOTrancheQuotes> CDOTrancheQuotesConstSP;

// Support for array
typedef array<CDOTrancheQuotesSP, CDOTrancheQuotes> CDOTrancheQuotesArray;
typedef smartPtr<CDOTrancheQuotesArray> CDOTrancheQuotesArraySP;
typedef smartConstPtr<CDOTrancheQuotesArray> CDOTrancheQuotesArrayConstSP;

// Support for wrapper
typedef MarketWrapper<CDOTrancheQuotes>   CDOTrancheQuotesWrapper;
typedef smartPtr<CDOTrancheQuotesWrapper> CDOTrancheQuotesWrapperSP;

// Support for wrapper arrays
typedef array<CDOTrancheQuotesWrapperSP,CDOTrancheQuotesWrapper> CDOTrancheQuotesWrapperArray; //note array of SP
typedef smartPtr<CDOTrancheQuotesWrapperArray>                   CDOTrancheQuotesWrapperArraySP;

/**
 * Container for CDO tranche market quotes
 * Eg1: contains spreads, expiries and upfront payments for ITRAXX 0%-3%
 * Eg2: contains spreads and expiries for CDX 3%-6%
 * */
class MARKET_DLL CDOTrancheQuotes : public MarketObject {
public:
    /** Creates CDOTrancheQuotes with imput quotes (and empty market name) */
    static CDOTrancheQuotesSP create(
        const DateTime& valueDate,
        double lowStrike,
        double highStrike,
        ExpiryArraySP expiries,
        CDoubleArraySP spreads,
        CDoubleArraySP upfronts = CDoubleArraySP(   )); // optional

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;
    
    /** Destructor */
    virtual ~CDOTrancheQuotes();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();
    
    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);
    
    /** Returns the name of this object */
    virtual string getName() const;
    
    /** Access to low strike */
    double getLowStrike() const;

    /** Access to high strike */
    double getHighStrike() const;
    
    /** Access to spreads */
    double getSpread(const DateTime& maturity) const;

    /** Access to upfronts */
    double getUpfront(const DateTime& maturity) const;
    
    /** Returns true if bid / ask data (spreads and upfronts) is available */
    bool hasBidAsk() const;

    /** Access to ask spreads */
    double getAskSpread(const DateTime& maturity) const;

    /** Access to ask upfronts */
    double getAskUpfront(const DateTime& maturity) const;

    /** Access to bid spreads */
    double getBidSpread(const DateTime& maturity) const;

    /** Access to bid upfronts */
    double getBidUpfront(const DateTime& maturity) const;
    
    /** Access to expiries as DateTime */
    DateTimeArraySP getMaturityDates() const;
    
    /** Specific clone method to copy "maturityoExpiryIdxMap" field */
    virtual IObject* clone() const;
    
private:
    /** Initialise maturityoExpiryIdxMap */
    void initMaturityoExpiryIdxMap();

    /**
     * Retrieve array index corresponding to "maturity",
     * throw an exception if not found
     * */
    int getExpiryIdx(const DateTime& maturity) const;

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Only build instances of that class using reflection */
    CDOTrancheQuotes();
    
    /** Default constructor */
    static IObject* defaultConstructor();

	// ------
	// FIELDS
	// ------

	/** Name of the market object [Mandatory] */
	string name;

	/** Attachment point [Mandatory] */
	double lowStrike;

	/** Detachment point [Mandatory] */
	double highStrike;

	/** Array of expiries [Mandatory] */
	ExpiryArraySP expiries;

	/** Array of (mid) spreads [Mandatory] */
	CDoubleArraySP spreads;
	
	/** Array of (mid) upfront payments [Optional] */
	CDoubleArraySP upfronts; 
	
	/** Array of bid spreads [Optional] */
	CDoubleArraySP bidSpreads;
	
	/** Array of bid upfront payments [Optional] */
	CDoubleArraySP bidUpfronts; 

	/** Array of ask spreads [Optional] */
	CDoubleArraySP askSpreads;
	
	/** Array of ask upfront payments [Optional] */
	CDoubleArraySP askUpfronts; 

    /** Value date [Transient] */
    DateTime valueDate;
    
    /**
     * Internal map between maturity and expiry index;
     * initialised by "getMarket()" [non exposed]
     * */
    map<const DateTime, int> maturityoExpiryIdxMap;     // $unregistered
};

DRLIB_END_NAMESPACE

#endif /*CDO_TRANCHE_QUOTES_HPP*/

