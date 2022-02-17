//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotes.hpp
//
//   Description : CDOQuotes is a market data container for index tranche quotes
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------


#ifndef CDO_QUOTES_HPP
#define CDO_QUOTES_HPP

#include "edginc/CDOTrancheQuotes.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/ICreditLegConvention.hpp"
#include "edginc/CDOPortfolio.hpp"
#include "edginc/ICDOQuotesGenerator.hpp"

DRLIB_BEGIN_NAMESPACE

// Forward declaration of CDOQuotes
class CDOQuotes;
class ICDOQuotesGenerator;

// Support for smart pointers
typedef smartPtr<CDOQuotes> CDOQuotesSP;
typedef smartConstPtr<CDOQuotes> CDOQuotesConstSP;

// Support for wrapper
typedef MarketWrapper<CDOQuotes> CDOQuotesWrapper;

/**
 * Market data container for index tranche quotes.
 * For a given index (eg: ITRAXX), contains the par spreads for
 * all liquid tranches (eg: 0-3%, 3-6%, 6-9%, 9-12% and 12-22%)
 * */
class PRODUCTS_DLL CDOQuotes : public MarketObject,
                               public virtual ICDOQuotesGenerator
{
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CDOQuotes();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Populates the object with the market data that this object needs */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Returns the name of this object */
    virtual string getName() const;

	/** Access to tranche quotes */
	CDOTrancheQuotesArrayConstSP getTrancheQuotes() const;

    /** Access to portfolio */
    CDOPortfolioConstSP getPortfolio() const;

    /** Access to discount curve name */
    string getDiscountName() const;

    /** Access to discount curve */
    YieldCurveConstSP getDiscount() const;

    /** Access to value date */
    const DateTime& getValueDate() const;

    /**
     * Generate a fee leg corresponding to the given
     * (lowStrike, highStrike, maturity) point.
     * Will fail if that point is not quoted.
     * */
    ICreditFeeLegSP generateFeeLeg(
        double lowStrike, double highStrike, const DateTime& maturity) const;

    /**
     * Generate a fee leg using this CDOQuotes convention but with
     * overriden spread and upfront payment
     * */
    ICreditFeeLegSP generateFeeLegOverride(
        double spread, double upfront, const DateTime& maturity) const;

    /** access to recover notional flag in creditLegConventions */
    bool getRecoverNotional() const;
    /**
     * Generate a contingent leg corresponding to the given
     * (lowStrike, highStrike, maturity) point.
     * Will fail if that point is not quoted.
     * */
    ICreditContingentLegSP generateContingentLeg(
        double lowStrike, double highStrike, const DateTime& maturity) const;

    /**
     * Generate a contingent leg using this CDOQuotes convention but with
     * overriden spread and upfront payment
     * */
    ICreditContingentLegSP generateContingentLegOverride(
        double spread, double upfront, const DateTime& maturity) const;

    /** Specific clone method to copy "strikesToTrancheQuotesIdxMap" field */
    virtual IObject* clone() const;


    /**
     * public Constructor.
     * Creates CDOQuotes object explicitely from
     * a flattened trancheQuotes
     * */
    CDOQuotes(
        DateTime         valueDate,
        CDOPortfolioSP   portfolio,
        ICreditLegConventionSP trancheConventions,
        ExpiryArraySP expiries,
        CDoubleArraySP lowStrikes,
        CDoubleArraySP highStrikes,
        CDoubleArraySP spreads,
        CDoubleArraySP upfronts = CDoubleArraySP(0));

    /**
     * Creates a new CDOQuotes object using the same conventions but
     * different tranche quotes.
     * The new CDOQuotes will have no "market name".
     * */
    CDOQuotesSP createWithNewQuotes(
        ExpiryArraySP expiries,
        CDoubleArraySP lowStrikes,
        CDoubleArraySP highStrikes,
        CDoubleArraySP spreads,
        CDoubleArraySP upfronts = CDoubleArraySP(0)) const;

    /** Main method : builds the CDOquotes object */
    CDOQuotesConstSP buildCDOQuotes() const;

private:
    /** Init strikesToTrancheQuotesIdxMap */
    void initStrikesToTrancheQuotesIdxMap();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    CDOQuotes();

    /** Default constructor */
    static IObject* defaultConstructor();

	// ------
	// FIELDS
	// ------

	/** Name of the market object [Mandatory] */
	string name;

	/** Tranche quotes [Mandatory] */
	CDOTrancheQuotesWrapperArraySP trancheQuotes;

    /** Portfolio (including single names betas and par spreads) [Mandatory] */
    CDOPortfolioSP portfolio;

    /** Tranche conventions [Mandatory] */
    ICreditLegConventionSP trancheConventions;

    /** Value date [Transient] */
    DateTime valueDate;

    /**
     * Internal map between (low strike, hisgh strike) pair and
     * tranche quotes index [non exposed]
     * */
    map<pair<double, double>, int> strikesToTrancheQuotesIdxMap;     // $unregistered
};

DRLIB_END_NAMESPACE

#endif /*CDO_QUOTES_HPP*/

