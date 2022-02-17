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


#ifndef CDO_QUOTES_BC_GENERATOR_HPP
#define CDO_QUOTES_BC_GENERATOR_HPP

#include "edginc/ICDOQuotesGenerator.hpp"
#include "edginc/CDOQuotesGenerator.hpp"
#include "edginc/CDOQuotes.hpp"
#include "edginc/IModel.hpp"
#include "edginc/CreditMetricsBaseCorrelation.hpp"
#include "edginc/IndexSkew.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Quote generator is a builder of CDOquotes
 * used whenever the market quotes are unavailable.
 * The CdoQuotes generator will generate the cdo quotes from 
 * the base correlation surface in the market environment
 * */
class PRODUCTS_DLL CDOQuotesBCGenerator :
    public ICDOQuotesGenerator,
    public MarketObject
{
public:
	
    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Destructor */
    virtual ~CDOQuotesBCGenerator();

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

    /** Main method : builds the CDOquotes object */
    CDOQuotesConstSP buildCDOQuotes() const;

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Only build instances of that class using reflection */
    CDOQuotesBCGenerator();

    /** Default constructor */
    static IObject* defaultConstructor();

	// ------
	// FIELDS
	// ------

	/** Name of the market object [Mandatory] */
	string name;

    /** wrapper for the index skew used in the base correlation model [Mandatory]*/
    IndexSkewWrapper indexSkew;

    /** fixed running coupon defined for the equity tranche [Mandatory]*/
    double runningEquityCoupon;

    /** Portfolio (including single names betas and par spreads) [Mandatory] */
    CDOPortfolioSP portfolio;

    /** Model used to generate the quotes [Optional]*/
    IModelSP bcModel;

    /** Tranche conventions [Mandatory] */
    ICreditLegConventionSP trancheConventions;

    /** Value date [Transient] */
    DateTime valueDate;

    /** CDO quotes that have been built [Transient]*/
    CDOQuotesConstSP cdoQuotes;

    /** generic CDOQuotes generator [Transient]*/
    CDOQuotesGeneratorSP quotesGenerator;
};

DECLARE(CDOQuotesBCGenerator);

typedef MarketWrapper<CDOQuotesBCGenerator> CDOQuotesBCGeneratorWrapper;

DRLIB_END_NAMESPACE

#endif /*CDO_QUOTES_BC_GENERATOR_HPP*/

