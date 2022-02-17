//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : CDOQuotesBootstrapperTimeOnly.hpp
//
//   Description : CDOQuotesBootstrapperTimeOnly defines an 'iterator' over CDO quotes. 
//					It differs from CDOQuotesBootstrapper in that it is for boostrapping 
//					in the time dimension only (not strike).
//
//   Author      : Matthias Arnsdorf
//
//   Date        : September 2006
//
//----------------------------------------------------------------------------


#ifndef CDO_QUOTES_BOOTSTRAPPER_TIME_ONLY_HPP
#define CDO_QUOTES_BOOTSTRAPPER_TIME_ONLY_HPP

#include "edginc/CDOQuotesBootstrapper.hpp"

DRLIB_BEGIN_NAMESPACE

/**
* CDOQuotesBootstrapperTimeOnly is a kind of "iterator" over CDOQuotes, used for
* bootstrapping in the time dimension. It creates "on the fly" the instruments 
* corresponding to the CDO quotes for a given maturity
* The order in which bootstrapping is done is 
* defined by the "next()" method.
* 
* */
class PRODUCTS_DLL CDOQuotesBootstrapperTimeOnly:
	public CDOQuotesBootstrapper
{
public:

	/** TYPE (for reflection) */
	static CClassConstSP const TYPE;

	/** Virtual destructor */
	virtual ~CDOQuotesBootstrapperTimeOnly();

	/** Constructor (external) */
	CDOQuotesBootstrapperTimeOnly(CDOQuotesConstSP cdoQuotes);

    /** Main method that runs bootstrap calibration and returns results */
    void bootstrap(
        const Calibrator & calibrator,
        int & nbVars,                                 // (O) number of variables calibrated
        Calibrator::InstanceIDArray & aggregIds,    // (0)
        DoubleArray & aggregVals,                   // (O) calibrated values
        DoubleArray & obFuncValues                  // (O) objective function values
        ){}; 

	/**
	* Method called before first step of the loop
	* [Implements IBootstrapper]
	* */
	virtual void init();

	/**
	* Method called after each step of the loop (ie go to the next quote)
	* [Implements IBootstrapper]
	* */
	virtual void next(); 

	/**
	* Method called to test the end of the loop (returns true when there is no more quote)
	* [Implements IBootstrapper]
	* */
	virtual bool end() const;

	/**
	* Returns "state" corresponding to current step of the loop
	* this is just a maturity in the time bootstrap only case
	* [Implements IBootstrapper]
	* */
	virtual IObjectSP getCurrentState() const;

	/** Returns the instruments corresponding to the current state (all quotes for current maturity) */
	virtual CInstrumentArraySP buildCurrentInstruments(); 

private:
	/** Invoked when Class is 'loaded' */
	static void load(CClassSP& clazz);

	/** Constructor (internal - used by reflection) */
	CDOQuotesBootstrapperTimeOnly();

	/** Default constructor */
	static IObject* defaultConstructor();

	// ------
	// Fields
	// ------


	// maturities for bootstrapping. Note that these are different from the maturities array in 
	// CDOQuotesBootstrpper. Here the size of the array is equal to the number of distinct maturities.
	// They are in increasing order which defines the bootstrapping order
	DateTimeArraySP bootstrapMaturities;

	// instruments for each bootstrap point
	vector< CInstrumentArraySP > instruments;


	// current index for maturity (corresponds to bootstrapping 'state'
	int currentMatIdx;

};

// support for smart pointers
typedef smartPtr<CDOQuotesBootstrapperTimeOnly> CDOQuotesBootstrapperTimeOnlySP;

DRLIB_END_NAMESPACE

#endif /*CDO_QUOTES_BOOTSTRAPPER_HPP*/

