//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesGenerator.hpp
//
//   Description : CDOQuotesGenerator is a generator of pseudo CDO quotes 
//                  from any model
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------


#ifndef CDO_QUOTES_GENERATOR_HPP
#define CDO_QUOTES_GENERATOR_HPP

#include "edginc/ICDOQuotesGenerator.hpp"
#include "edginc/CDOQuotes.hpp"
#include "edginc/IModel.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Quote generator is a builder of CDOquotes
 * used whenever the market quotes are unavailable.
 * The CdoQuotes generator will generate the cdo quotes from 
 * the given model and engine parameters
 * */
class PRODUCTS_DLL CDOQuotesGenerator :
    public ICDOQuotesGenerator,
    public MarketObject
{
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CDOQuotesGenerator();

    /** Called immediately after object constructed */
    virtual void validatePop2Object();

    /** Access to portfolio */
    CDOPortfolioConstSP getPortfolio() const;

    /** Access to value date */
    const DateTime& getValueDate() const;

    /** Get the market data needed */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Get the name of the market object */
    virtual string getName() const;

    /** Get the underlying built CDOQuotes */
    virtual CDOQuotesConstSP getCDOQuotes() const;

    /** Main method : builds the CDOquotes object */
    CDOQuotesConstSP buildCDOQuotes() const;

    /** public constructor */
    CDOQuotesGenerator( string name,
                        IModelSP model,
                        DoubleArraySP lowStrikes,
                        DoubleArraySP highStrikes,
                        ExpiryArraySP expiries,
                        CDOPortfolioSP  portfolio,
                        CreditEngineParametersWrapper engineParameters, 
                        ICreditLegConventionSP trancheConventions,
                        double                 runningEquityCoupon);


private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    CDOQuotesGenerator();

    /** Default constructor */
    static IObject* defaultConstructor();

	// ------
	// FIELDS
	// ------

	/** Name of the market object [Mandatory] */
	string name;

    /** wrapper for the ithe engine parameters */
    CreditEngineParametersWrapper engineParameters;

    /** fixed running coupon defined for the equity tranche */
    double runningEquityCoupon;

    /** Portfolio (including single names betas and par spreads) [Mandatory] */
    CDOPortfolioSP portfolio;

    /** Model used to generate the quotes */
    IModelSP model;

    /** Tranche conventions [Mandatory] */
    ICreditLegConventionSP trancheConventions;

	/** expiries of the built quotes */
	ExpiryArraySP expiries;

	/** low strikes of built quotes */
	DoubleArraySP lowStrikes;

	/** high strikes of built quotes */
	DoubleArraySP highStrikes;

    /** Value date [Transient] */
    DateTime valueDate;

    /** CDO quotes that have been built */
    CDOQuotesConstSP cdoQuotes;
};

DECLARE(CDOQuotesGenerator);

typedef MarketWrapper<CDOQuotesGenerator> CDOQuotesGeneratorWrapper;

DRLIB_END_NAMESPACE

#endif /*CDO_QUOTES_BC_GENERATOR_HPP*/

