//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiMarketFactors.hpp
//
//   Description : MultiMarketFactors interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MULTIMARKETFACTOR_HPP
#define EDR_MULTIMARKETFACTOR_HPP
#include "edginc/MarketFactor.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
class IModel;
class MarketData;
class DateTime;
class YieldCurve;
class CInstrument;
class Sensitivity;
class CDoubleMatrix;
class OutputRequest;
class Results;
typedef smartConstPtr<CDoubleMatrix> CDoubleMatrixConstSP;

class IMultiMarketFactors;
typedef smartPtr<IMultiMarketFactors> IMultiMarketFactorsSP;
#ifndef QLIB_MULTIFACTORS_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IMultiMarketFactors>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IMultiMarketFactors>);
#endif

typedef smartConstPtr<IMultiMarketFactors> IMultiMarketFactorsConstSP;
/** A MultiMarketFactors is a view onto a collection of MarketFactors
    objects and how they interact (ie correlations) */
class MARKET_DLL IMultiMarketFactors: public virtual IObject{
public:
    static CClassConstSP const TYPE; // in MultiFactors.cpp

    virtual ~IMultiMarketFactors();

    /** Pull out the necessary (eg assets, correlations) from the market data */
    virtual void getMarket(const IModel* model, const MarketData* market) = 0;
    
    /** Returns the number of factors in this collection (this does not include
        MarketFactors which are contained within top level MarketFactors) */
    virtual int numFactors() const = 0;

    /** For backwards compatibility: same as numFactors */
    virtual int NbAssets() const = 0;

    //// Returns the correlations between the immediate factors.
    //// This is numFactors() x numFactors() 
    virtual CDoubleMatrixConstSP factorsCorrelationMatrix() const = 0;

    /** Validate that the market data and instrument are consistent with
        each other etc */
    virtual void crossValidate(const DateTime&       startDate,
                               const DateTime&       valueDate,
                               const YieldCurve*     discCcy,
                               const CInstrument*    instrument) const = 0;

    /** Returns the name of the specified 'market factor' */
    virtual string getName(int index) const = 0;
    
    /** Returns the MarketFactor with the specifed index */
    virtual IMarketFactorConstSP getFactor(int index) const = 0;

    /** Returns a list of asset indexes indicating which asset is
        sensitive to supplied SensControl. If crossAssetSensitivities
        is true then this only covers sensitivities such as phi which
        can live at this level otherwise these are excluded. The
        relevant name of what is being tweaked must be stored in the
        SensControl */
    virtual IntArray getSensitiveFactors(
        const Sensitivity* sens,
        bool               crossAssetSensitivities) const = 0;

    /** A utility to record 'forwards' at maturity (for each factor) as
        appropriate for each factor (if applicable for that factor) */
    virtual void recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const = 0;

    /**
       Others as per what goes into MarketFactors class
       Plus something for correlations - to do
     */

    /** Make a single market factor look like a IMarketFactor. The supplied
        marketFactor must already been populated with its market data */
    static IMultiMarketFactorsSP asMulti(
        IMarketFactorConstSP marketFactor);  // in MultiFactors.cpp
private:
    static void load(CClassSP& clazz); // in MultiFactors.cpp
    class AsMulti;
};
DRLIB_END_NAMESPACE
#endif
