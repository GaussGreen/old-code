//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AbsCdoParameters.hpp
//
//   Description : Portfolio level parameters needed for AbsCdo model.
//
//   Author      : Jay Wang
//
//   Date        : Dec 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_ABSCDOPARAMETERS_HPP
#define QLIB_ABSCDOPARAMETERS_HPP

#include "edginc/RationalisedCreditEngineParameters.hpp"
#include "edginc/SkewSurface.hpp"
DRLIB_BEGIN_NAMESPACE

/** Per portfolio parameters needed for AbsCdo model */
class MARKET_DLL AbsCdoParameters: public RationalisedCreditEngineParameters
{
public:
    static CClassConstSP const TYPE;
    
    /** Destructor */
    ~AbsCdoParameters();

    /** Check the model parameters */
    virtual void validatePop2Object();

    /** Returns the name of this object. This is the name with which
        it is stored in the market data cache and is the name with
        which results (eg tweaks) should be reported against */
    virtual string getName() const;

	double getLossDecBeta() const
	{ return lossDecBeta; }

	SkewSurfaceConstSP getLossSkew() const 
	{ return lossSkew.getSP(); }

	SkewSurfaceConstSP getDecSkew() const 
	{ return decSkew.getSP(); }

	SkewSurfaceConstSP getLossDecSkew() const
	{ return lossDecSkew.getSP(); }

	/** Pull out the skews the market data */
	void getMarket(const IModel* model, const MarketData* market);
    
    static AbsCdoParameters* buildDefault();

private:
    friend class StopCompilerWarning;

    AbsCdoParameters();
    AbsCdoParameters(const AbsCdoParameters& rhs); // don't use
    AbsCdoParameters& operator=(const AbsCdoParameters& rhs); // don't use

    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    string   name;         

    double   lossDecBeta;			// beta between loss and decretion
	SkewSurfaceWrapper lossSkew;	// loss BC
	SkewSurfaceWrapper decSkew;		// dec BC
	SkewSurfaceWrapper lossDecSkew;	// lossDec BC
};

// Support for smart pointers
typedef smartPtr<AbsCdoParameters> AbsCdoParametersSP;
typedef smartConstPtr<AbsCdoParameters> AbsCdoParametersConstSP;

// Support for wrapper
typedef MarketWrapper<AbsCdoParameters> AbsCdoParametersWrapper;
#ifndef QLIB_ABSCDOPARAMETERS_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<AbsCdoParameters>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<AbsCdoParameters>);
#endif

DRLIB_END_NAMESPACE

#endif
