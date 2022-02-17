//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesBootstrapper.hpp
//
//   Description : CDOQuotesBootstrapper is a kind of "iterator" over CDOQuotes (see below)
//
//   Author      : Antoine Gregoire
//
//   Date        : December 2005
//
//----------------------------------------------------------------------------


#ifndef CDO_QUOTES_BOOTSTRAPPER_HPP
#define CDO_QUOTES_BOOTSTRAPPER_HPP

#include "edginc/CDOQuotes.hpp"
#include "edginc/IBootstrapper.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * CDOQuotesBootstrapper is a kind of "iterator" over CDOQuotes, responsible for 
 * creating "on the fly" the instruments corresponding to each single CDO quote
 * (i.e. {low strike, high strike, maturity}).
 * It can be used for bootstrapping. The order in which bootstrapping is done is 
 * defined by the "next()" method.
 * 
 * Reasons to define a class distinct from CDOQuotes:
 * 1. Needs to access "product information" (inside products/ directory)
 * 2. Wants to keep "status" information distinct from the market
 *    object itself
 * */
class PRODUCTS_DLL CDOQuotesBootstrapper:
    public CObject
{
public:

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~CDOQuotesBootstrapper();
    
	/** Constructor (external) */
	CDOQuotesBootstrapper(CDOQuotesConstSP cdoQuotes, bool ignore100pc = true);
	
    /** Main method that runs bootstrap calibration and returns results */
    // should remove inheritance from IBootstrapper //TODO
    void bootstrap(
        const Calibrator & calibrator,
        int & nbVars,                                 // (O) number of variables calibrated
        Calibrator::InstanceIDArray & aggregIds,    // (0)
        DoubleArray & aggregVals,                   // (O) calibrated values
        DoubleArray & obFuncValues                  // (O) objective function values
        )
        {
        throw ModelException("This class can't bootstrap", "CDOQuotesBootstrapper::bootstrap");   
        };  

    /**
     * Method called before first step of the loop
     *
     * */
    virtual void init();

    /**
     * Method called after each step of the loop (ie go to the next quote)
     * 
     * */
    virtual void next(); 

    /**
     * Method called to test the end of the loop (returns true when there is no more quote)
     * 
     * */
    virtual bool end() const;

    /**
     * Returns "state" corresponding to current step of the loop
     * 
     * */
    virtual IObjectSP getCurrentState() const;
    
	

	/** Returns the instruments corresponding to the current state */
	virtual CInstrumentArraySP buildCurrentInstruments(); 

protected:
	/** Returns the (single) instrument corresponding to the current quote */
	CInstrumentSP buildCurrentInstrument(); 

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);
    
    /** Constructor (internal - used by reflection) */
    CDOQuotesBootstrapper();
    
    /** Default constructor */
    static IObject* defaultConstructor();

	

    // ------
    // Fields
    // ------
    
    /** Pointer to CDO quotes */
    CDOQuotesConstSP cdoQuotes;

    /** flag to indicate whether we ignore 100% point. default = true */
    bool ignore100pc;

    // expiries, lowStrikes and highStrikes define an ordered
    // array of (maturity, strike1, strike2) points which corresponds to
    // the bootstrapping order
protected: // CDOQuotesBootstrapperTimeOnly needs access
    DateTimeArrayConstSP maturities;
private:
    DoubleArraySP lowStrikes;
    DoubleArraySP highStrikes;
    
    // Current index in expiries / quotes
    int currentPointIdx;
    
    /** Total notional of the CDOQuotes portfolio */
    double portfolioNotional;
    
    /** Portfolio shared by all quotes */
    CDOPortfolioSP sharedPortfolio;

	/** cache for all instruments to value */
	CInstrumentArraySP instruments;
};

// support for smart pointers
typedef smartPtr<CDOQuotesBootstrapper> CDOQuotesBootstrapperSP;

DRLIB_END_NAMESPACE

#endif /*CDO_QUOTES_BOOTSTRAPPER_HPP*/

