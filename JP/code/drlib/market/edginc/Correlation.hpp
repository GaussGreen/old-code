
#ifndef EDR_CORRELATION_HPP
#define EDR_CORRELATION_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/CorrelationCommon.hpp"
#include "edginc/Correl.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/CorrSwapBasisAdj.hpp"		// for CorrSwapBasis
#include "edginc/CorrSwapSamplingAdj.hpp"	// for CorrSwapBasis

DRLIB_BEGIN_NAMESPACE
class Correlation;
typedef smartPtr<Correlation> CorrelationSP;
typedef smartConstPtr<Correlation> CorrelationConstSP;
typedef array<CorrelationSP, Correlation> CorrelationArray;
#ifndef QLIB_CORRELATION_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<Correlation>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<Correlation>);
EXTERN_TEMPLATE(class MARKET_DLL array<CorrelationSP _COMMA_ Correlation>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<Correlation>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<Correlation>);
INSTANTIATE_TEMPLATE(class MARKET_DLL array<CorrelationSP _COMMA_ Correlation>);
#endif

/* 
 * This is just a simple draft implementation of something which 
 * will be great and wonderful one day 
 */
/** This class is capable of being sensitive to any of the different 
    correlation tweaks. Use the configureForSensitivities to specify which 
    one */
class MARKET_DLL Correlation: public CorrelationCommon,
                              public virtual ITweakableWithRespectTo<Correl> {
public:
    static CClassConstSP const TYPE;

    static const string BENCHMARK_EXPIRY;
    static const int BENCHMARK_TIME;

    /** Validation */
    void validatePop2Object();

	virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the correlation as a simple double */
    double getCorrelation() const;

    /** returns the correlation expiry (used if there is term structure of correlation) */
    ExpirySP getCorrExpiry() const;
    
    /** returns the correlation's name */
    virtual string getName() const;

    /** Configure this correlation object so that under tweaking it behaves
        properly when it is a correlation between object of type clazz1 and
        an object of type clazz2 */
    virtual void configureForSensitivities(CClassConstSP clazz1,
                                           CClassConstSP clazz2);
    /** Returns true if this correlation is [really] sensitive to the
        supplied sensitivity */
    virtual bool isSensitiveTo(const IPerNameSensitivity* sens) const;
    
    /** Returns the name of the stock/asset - used to determine
        whether to tweak the object */
    virtual string sensName(const Correl*) const;

    /** Shifts the object using given shift */    
    virtual TweakOutcome sensShift(const PropertyTweak<Correl>& tweak);
    
    /** constructor */
    Correlation(const string& name,        // optional
                const string& nameAsset1,
                const string& nameAsset2,
                double        correlation);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultCorrelation();
    Correlation();
    enum SensType{
        NONE = 0,
        PHI,
        FX_PHI};
    ////// fields ////
    double          correlation;
    ExpirySP        corrExpiry;
    int             sensType; // transient
	// transient, but tweakable
    CorrSwapSamplingAdjSP   samplingAdj;
    CorrSwapBasisAdjSP      basisAdj1;
    CorrSwapBasisAdjSP      basisAdj2;	
    double                  timeDiff;
};

DRLIB_END_NAMESPACE
#endif // CORRELATION_HPP
