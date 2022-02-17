
#ifndef EDR_FACTOR_CORRELATION_HPP
#define EDR_FACTOR_CORRELATION_HPP

#include "edginc/MarketObject.hpp"
#include "edginc/CorrelationCommon.hpp"
/*
#include "edginc/Phi.hpp"
#include "edginc/FXPhi.hpp"
#include "edginc/CorrelationSqueeze.hpp"
#include "edginc/CorrSwapBasisAdj.hpp"		// for CorrSwapBasis
#include "edginc/CorrSwapSamplingAdj.hpp"	// for CorrSwapBasis
#include "edginc/Expiry.hpp"
*/
#include "edginc/DoubleMatrix.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class FactorCorrelation;
DECLARE(FactorCorrelation)
/*
#ifndef QLIB_CORRELATION_CPP
EXTERN_TEMPLATE(class smartConstPtr<FactorCorrelation>);
EXTERN_TEMPLATE(class smartPtr<FactorCorrelation>);
EXTERN_TEMPLATE(class array<FactorCorrelationSP _COMMA_ FactorCorrelation>);
#else
INSTANTIATE_TEMPLATE(class smartConstPtr<FactorCorrelation>);
INSTANTIATE_TEMPLATE(class smartPtr<FactorCorrelation>);
INSTANTIATE_TEMPLATE(class array<FactorCorrelationSP _COMMA_ FactorCorrelation>);
#endif
*/
/* 
 * This is just a simple draft implementation of something which 
 * will be great and wonderful one day 
 */
/** This class is capable of being sensitive to any of the different 
    correlation tweaks. Use the configureForSensitivities to specify which 
    one */
class MARKET_DLL FactorCorrelation: public CorrelationCommon
{
public:
    static CClassConstSP const TYPE;

    static const string BENCHMARK_EXPIRY;
    static const int BENCHMARK_TIME;

    /** Validation */
    void validatePop2Object();

	//virtual void getMarket(const IModel* model, const MarketData* market);

    /** returns the correlation as a simple double */
    const DoubleMatrix& getCorrelation() const;

    virtual void configureForSensitivities(CClassConstSP clazz1,
        CClassConstSP clazz2)
    {}

    /** Returns true if this correlation is [really] sensitive to the
        supplied sensitivity */
    virtual bool isSensitiveTo(const IPerNameSensitivity* sens) const
    {
        return false;
    }

    /** constructor */
    FactorCorrelation(const string& name,        // optional
                const string& nameAsset1,
                const string& nameAsset2,
                int nbFactor1,
                int nbFactor2,
                const DoubleArray& values);

    FactorCorrelation(const string& name,        // optional
        const string& nameAsset1,
        const string& nameAsset2,
        const DoubleMatrix& values);

private:
//    bool sensShiftCommon(ScalarShift* shift, bool overrideShiftSizeSign);
    static void load(CClassSP& clazz);
    static IObject* defaultCorrelation();
    FactorCorrelation();
    DoubleMatrix    correlation;
};

DRLIB_END_NAMESPACE
#endif // FACTORCORRELATION_HPP
