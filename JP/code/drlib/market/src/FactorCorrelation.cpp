
#include "edginc/config.hpp"
#define QLIB_CORRELATION_CPP
#include "edginc/FactorCorrelation.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Addin.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/Addin.hpp"
#include "edginc/MultiMarketFactors.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/MarketDataFetcherLNSpline.hpp"

DRLIB_BEGIN_NAMESPACE

/** Validation */
void FactorCorrelation::validatePop2Object(){
    static const string method("FactorCorrelation::validatePop2Object");
    if (name.empty()){
        name = asset1 < asset2? (asset1+ "_" + asset2): (asset2 + "_" + asset1);
    } else if (name == asset1 || name == asset2){
        throw ModelException(method, "FactorCorrelation's name must be different to"
                             " either asset's name");
    }
    for (int i=0; i<correlation.numCols(); ++i)
        for (int j=0; j<correlation.numRows(); ++j)
        {
            if (fabs(correlation[i][j]) > 1.0){
                throw ModelException( method, "Correlation "+name+
                                    " must be between -1.0 and 1.0.");
            }
        }
}

const DoubleMatrix& FactorCorrelation::getCorrelation() const
{
    static const string method("FactorCorrelation::getCorrelation");
    return correlation;
}

FactorCorrelation::FactorCorrelation(const string& _name,
                         const string& _nameAsset1,
                         const string& _nameAsset2,
                         int _nbFactor1,
                         int _nbFactor2,
                         const DoubleArray& _correlation):
    CorrelationCommon(FactorCorrelation::TYPE, _name, _nameAsset1, _nameAsset2), 
    correlation(_nbFactor2, _nbFactor1)
{
    int n=0;
    for (int col=0; col<_nbFactor2; ++col)
        for (int row=0; row<_nbFactor1; ++row)
        {
            correlation[col][row] = _correlation[n++];
        }
    validatePop2Object();
}

FactorCorrelation::FactorCorrelation(const string& name,        // optional
                  const string& nameAsset1,
                  const string& nameAsset2,
                  const DoubleMatrix& values)
                  : CorrelationCommon(FactorCorrelation::TYPE, name, nameAsset1, nameAsset2), correlation(values)
{
    validatePop2Object();
}

FactorCorrelation::FactorCorrelation():CorrelationCommon(TYPE), correlation() 
{}

/** Invoked when Class is 'loaded' */
void FactorCorrelation::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FactorCorrelation, clazz);
    SUPERCLASS(CorrelationCommon);
    EMPTY_SHELL_METHOD(defaultCorrelation);
    clazz->enableCloneOptimisations();
    FIELD(name, "name for correlation");
    FIELD_MAKE_OPTIONAL(name);
    FIELD(asset1,      "Name of Asset 1");
    FIELD(asset2,      "Name of Asset 2");
    FIELD(correlation, "Factor Correlation as DoubleMatrix, "
        "with different columns for different factors of Asset 2. ");
}

IObject* FactorCorrelation::defaultCorrelation(){
    return new FactorCorrelation();
}

CClassConstSP const FactorCorrelation::TYPE = CClass::registerClassLoadMethod(
    "FactorCorrelation", typeid(FactorCorrelation), load);

// initialise type for array of Correlations
template<> CClassConstSP const FactorCorrelationArray::TYPE = CClass::registerClassLoadMethod(
    "FactorCorrelationArray", typeid(FactorCorrelationArray), load);

DRLIB_END_NAMESPACE
