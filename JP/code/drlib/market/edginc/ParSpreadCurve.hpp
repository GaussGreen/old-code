//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ParSpreadCurve.hpp
//
//   Description : a par credit spread curve
//
//   Author      : Andre Segger
//
//
//----------------------------------------------------------------------------

#ifndef PAR_SPREAD_CURVE_HPP
#define PAR_SPREAD_CURVE_HPP

#include "edginc/Class.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/ParSpreadParallel.hpp"
#include "edginc/ParSpreadParallelRelative.hpp"
#include "edginc/PointwiseTweakableWith.hpp"
#include "edginc/ParSpreadPointwise.hpp"
#include "edginc/ParSpreadUpfronts.hpp"
#include "edginc/ParSpreadUpfrontParallel.hpp"
#include "edginc/ParSpreadLevel.hpp"
#include "edginc/ParSpreadLevelAtGivenDate.hpp"
#include "edginc/ParSpreadParallelShift.hpp"
#include "edginc/ParSpreadPropShift.hpp"
#include "edginc/ParSpreadWeightedShift.hpp"
#include "edginc/QuasiContractualBaseCorrelation.hpp"
#include "edginc/Theta.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL ParSpreadCurve : public MarketObject,
                       virtual public ITweakableWithRespectTo<ParSpreadParallel>,
                       virtual public ITweakableWithRespectTo<ParSpreadParallelRelative>,
                       virtual public ITweakableWithRespectTo<ParSpreadPointwise>,
                       virtual public ITweakableWithRespectTo<ParSpreadUpfrontParallelTP>,
                       virtual public ITweakableWithRespectTo<ParSpreadUpfronts>,
                       virtual public QuasiContractualBaseCorrelation::IShift,
                       virtual public ParSpreadLevel::IShift,
                       virtual public ParSpreadLevelAtGivenDate::IShift,
                       virtual public ParSpreadParallelShift::IShift,
                       virtual public ParSpreadPropShift::IShift,
                       virtual public ParSpreadWeightedShift::IShift,
                       virtual public Theta::IShift,
                       public virtual Calibrator::IAdjustable
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& );

    ParSpreadCurve();

    ParSpreadCurve(const string&    name,
                   const double&    parSpread,
                   const string&    maturity);
                   
    ParSpreadCurve(const string&    name,
                   const double&    parSpread,
                   const DateTime&  maturity);
                   
    ParSpreadCurve(
        const string& name,
        ExpiryArraySP expiries,
        const CDoubleArray& spreads);                 

    // constructors with upfront fees
    ParSpreadCurve(const string&    name,
                   const double&    parSpread,
                   const double&    upfront,
                   const string&    maturity);
                   
    ParSpreadCurve(const string&       name,
                   ExpiryArraySP       expiries,
                   const CDoubleArray& spreads,
                   CDoubleArraySP      upfronts);                 

    virtual void validatePop2Object();

    virtual string getName() const;

    // Get today from the market data cache
    void getMarket(const IModel* model, const MarketData *market);

    double getCurrentSpread(const DateTime& today,
                            const DateTime& maturityDate,
                            const BadDayConvention* bdc,
                            const Holiday* hols) const;

    // Bucketwise shift
    // Outside the usual sensitivity framework
    void bucketShift(const DoubleArray& shifts);

    // Parallel tweak (sensitivity) support
    virtual string sensName(const ParSpreadParallel*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadParallel>&);

    virtual string sensName(const ParSpreadUpfrontParallelTP*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadUpfrontParallelTP>&);

    // Parallel tweak relative (sensitivity) support
    virtual string sensName(const ParSpreadParallelRelative*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadParallelRelative>&);

    // Pointwise tweak (sensitivity) support
    virtual string sensName(const ParSpreadPointwise*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadPointwise>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const ParSpreadPointwise*) const;

    virtual string sensName(const ParSpreadUpfronts*) const;
    virtual TweakOutcome sensShift(const PropertyTweak<ParSpreadUpfronts>& shift);
    virtual ExpiryWindowArrayConstSP sensQualifiers(const ParSpreadUpfronts*) const;
    

    // Proportional shift support
    virtual string sensName(ParSpreadPropShift* shift) const;
    virtual bool sensShift(ParSpreadPropShift* shift);

    // Parallel shift (scenario) support
    virtual string sensName(ParSpreadParallelShift* shift) const;
    virtual bool sensShift(ParSpreadParallelShift* shift);

    // Flat override support
    virtual string sensName(ParSpreadLevel* shift) const;
    virtual bool sensShift(ParSpreadLevel* shift);

    // Flat override support
    virtual bool sensShift(QuasiContractualBaseCorrelation* shift);
    
    // ParSpreadLevelAtGivenDate::IShift implementation
    virtual string sensName(ParSpreadLevelAtGivenDate* shift) const;
    virtual bool sensShift(ParSpreadLevelAtGivenDate* shift);

    // Additive or multiplicative shift
    virtual string sensName(ParSpreadWeightedShift* shift) const;
    virtual bool sensShift(ParSpreadWeightedShift* shift);

    // Theta shift
    virtual string sensName(Theta* shift) const;
    virtual bool sensShift(Theta* shift);

    // Get Expiries
    ExpiryArrayConstSP getExpiries() const;

    // Get Par Spreads
    DoubleArrayConstSP getParSpreads() const;

    // Get upfront fees
    DoubleArrayConstSP getUpfronts() const;

protected:
    // no additional data
    static void acceptWrapperNameCollector(const ParSpreadCurve* parSpreadCurve,
                                           WrapperNameCollector* collector);

    DateTime        today;    // Today's date - Transient
    string          name;
    ExpiryArraySP   expiries;
    CDoubleArray    spreads;
    CDoubleArraySP  upfronts;

};

typedef smartConstPtr<ParSpreadCurve> ParSpreadCurveConstSP;
typedef smartPtr<ParSpreadCurve>      ParSpreadCurveSP;
typedef MarketWrapper<ParSpreadCurve> ParSpreadCurveWrapper;
#ifndef QLIB_PARSPREADCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<ParSpreadCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<ParSpreadCurve>);
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<ParSpreadCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<ParSpreadCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<ParSpreadCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<ParSpreadCurve>);
#endif

DRLIB_END_NAMESPACE
#endif
