//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierProcessSVJ.hpp
//
//   Date        : 09 March 2004
//
//
//----------------------------------------------------------------------------
#ifndef EDR_FOURIERPROCESSSVJ_HPP
#define EDR_FOURIERPROCESSSVJ_HPP

#include "edginc/VolSVJ.hpp"
#include "edginc/FourierProcessParametric.hpp"

DRLIB_BEGIN_NAMESPACE

class FOURIER_DLL FourierProcessSVJ: public FourierProcessParametric,
                         virtual public StFourierProcessLogRtn,      // implements started log return
                         virtual public FwdStFourierProcessLogRtn,   // implements fwd starting log return
                         virtual public StFourierProcessIntVar,      // implements started integrated variance
                         virtual public FwdStFourierProcessIntVar,   // implements fwd starting integrated variance
                         virtual public StFourierProcessQuadVar,     // implements started quad variance
                         virtual public FwdStFourierProcessQuadVar,  // implements fwd starting quad variance
						 virtual public FwdStFourierProcessExpQuadVar{  
public:
    static CClassConstSP const TYPE; 

    /* Started Log Return */
    virtual Complex scalelessCumulant(const StFourierProductLogRtn& product, 
                                      const Complex& z, 
                                      const DateTime& matDate) const;
    virtual double lowerRealBound(const StFourierProductLogRtn& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const StFourierProductLogRtn& product, 
                                  const DateTime& matDate) const;
    /* Fwd starting Log Return  */
    virtual Complex scalelessCumulant(const FwdStFourierProductLogRtn& product, 
                                      const Complex& z, 
                                      const DateTime& matDate) const;
    virtual double lowerRealBound(const FwdStFourierProductLogRtn& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const FwdStFourierProductLogRtn& product, 
                                  const DateTime& matDate) const;


/**ARNAUD**////////////////////////////////////////////////////////////////////////////
	/* Fwd Starting Expected Quadratic Variation */
    virtual Complex cumulant(const FwdStFourierProductExpQuadVar& product, 
                                      const Complex& z, 
                                      const DateTime& matDate) const;
    virtual double lowerRealBound(const FwdStFourierProductExpQuadVar& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const FwdStFourierProductExpQuadVar& product, 
                                  const DateTime& matDate) const;

/////////////////////////////////////////////////////////////////////////////////////////





    /* Started integrated variance */
    virtual Complex cumulant(const StFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const;   
    virtual double lowerRealBound(const StFourierProductIntVar& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const StFourierProductIntVar& product, 
                                  const DateTime& matDate) const;
    
    /* Forward Starting integrated variance */
    virtual Complex cumulant(const FwdStFourierProductIntVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const;    
    virtual double lowerRealBound(const FwdStFourierProductIntVar& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const FwdStFourierProductIntVar& product, 
                                  const DateTime& matDate) const;
    
    /* Started quad variation */
    virtual Complex cumulant(const StFourierProductQuadVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const; 
    virtual double lowerRealBound(const StFourierProductQuadVar& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const StFourierProductQuadVar& product, 
                                  const DateTime& matDate) const;
    virtual double expectation(const StFourierProductQuadVar& product,
                               const DateTime& matDate) const;

    /* Forward starting quad variation */
    virtual Complex cumulant(const FwdStFourierProductQuadVar& product, 
                             const Complex z, 
                             const DateTime& matDate) const; 
    virtual double lowerRealBound(const FwdStFourierProductQuadVar& product, 
                                  const DateTime& matDate) const;
    virtual double upperRealBound(const FwdStFourierProductQuadVar& product, 
                                  const DateTime& matDate) const;
    virtual double expectation(const FwdStFourierProductQuadVar& product,
                               const DateTime& matDate) const;

    /** Create a MarketDataFetcher which will be used by the [Fourier] model
     *  for retrieving market data etc */
    virtual MarketDataFetcherSP marketDataFetcher() const;

    virtual const TimeMetric& getTimeMetric() const;

    /** Returns the parameters of the process */
    virtual string getParameters() const;

    /** Give the process the chance to do some extra product specific 
        validation and to initialize its transient fields */
    virtual void validate(const FourierProduct* product);
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    FourierProcessSVJ();

private:

    static IObject* defaultCtor();

    /** Optional frequency range field. Is optional because
        process will not always need to provide it if
        the payoff does not require it */
    FrequencyRangeConstSP stFrequencyRange;
    FrequencyRangeConstSP fwdStFrequencyRange;
    
    FrequencyRangeConstSP stFrequencyRangeIntVar;
    FrequencyRangeConstSP fwdStFrequencyRangeIntVar;

    // transient
    VolSVJSP theVol;
    TimeMetricSP timeMetric;
};

typedef smartPtr<FourierProcessSVJ> FourierProcessSVJSP;

DRLIB_END_NAMESPACE
#endif
