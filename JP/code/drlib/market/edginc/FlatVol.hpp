//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FlatVol.hpp
//
//   Description : Trivial Implementation of Vol Interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef FLATVOL_HPP
#define FLATVOL_HPP
#include "edginc/VolBase.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/RootTimeVega.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/Theta.hpp"
#include "edginc/PowerVega.hpp"
#include "edginc/PDFCalculator.hpp"
#include "edginc/VolRelativeShift.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/** Simple implementation of volatility where the volatility is captured by
    a single scalar. */
class MARKET_DLL FlatVol: public CVolBase,
               virtual public IPDFCalculator,
               virtual public IVolatilityBS,
               virtual public IRestorableWithRespectTo<VolParallel>,
               virtual public RootTimeVega::IRestorableShift,
               virtual public VolLevel::Shift,
               virtual public VolParallelShift::Shift,
               virtual public PowerVega::Shift,
               virtual public Theta::IShift,
               virtual public VolRelativeShift::IShift,
               virtual public VolAbsoluteShift::IShift {
public:
    static CClassConstSP const TYPE;
    friend class FlatVolHelper;
    friend class FlatVolAddin;
    friend class FlatInterpVol;

    /** Returns name of vol */
    virtual string getName() const;

    /** constructor needed as we no addin to begin with */
    FlatVol(const string&     name, 
            const DateTime&   baseDate,
            const TimeMetric* timeMetric,
            double            flatVol);

    virtual void validatePop2Object();

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    /** Combines market and instrument data together to give a
        Processed Vol. Here the processed volatility is a processed
        struck volatility ie it reflects the combination of this
        CVolBase together with the supplied FX asset and the
        correlation between this CVolBase and the vol of the
        FX. */
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;

    /** Returns name identifying vol for vega parallel */
    virtual string sensName(const VolParallel*) const;
    /** Shifts the object using given shift */
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);
    /** Restores the object to its original form */
    virtual void sensRestore(const PropertyTweak<VolParallel>& tweak);
    /** Returns name identifying vol for RootTimeVega */
    virtual string sensName(RootTimeVega* shift) const;
    /** Shifts the object using given shift for RootTimeVega. Just does a
        vega parallel tweak */
    virtual bool sensShift(RootTimeVega* shift);
    /** Restores the object to its original form after a RootTimeVega tweak */
    virtual void sensRestore(RootTimeVega* shift);

    /** Implements VolLevel scenario */
    /** Returns name identifying this object for VolLevel */
    virtual string sensName(VolLevel* shift) const;
    /** Shifts the object using given shift (see VolLevel::Shift)*/
    virtual bool sensShift(VolLevel* shift);

    /** Implements VolParallelShift scenario */
    /** Returns name identifying this object for VolParallelShift */
    virtual string sensName(VolParallelShift* shift) const;
    /** Shifts the object using given shift (see VolParallelShift::Shift)*/
    virtual bool sensShift(VolParallelShift* shift);

    /** Implements PowerVega scenario */
    /** Returns name identifying this object for PowerVega */
    virtual string sensName(PowerVega* shift) const;
    /** Shifts the object using given shift (see PowerVega::Shift)*/
    virtual bool sensShift(PowerVega* shift);
   
    /** Shifts the object using given shift. */
    virtual bool sensShift(Theta* shift);

    /** Implements VolRelativeShift scenario */
    /** Returns name identifying this object for VolRelativeShift */
    virtual string sensName(VolRelativeShift* shift) const;
    /** Shifts the object using given shift (see VolRelativeShift::IShift)*/
    virtual bool sensShift(VolRelativeShift* shift);
 
    /** Implements VolAbsoluteShift scenario */
    /** Returns name identifying this object for VolAbsoluteShift */
    virtual string sensName(VolAbsoluteShift* shift) const;
    /** Shifts the object using given shift (see VolAbsoluteShift::IShift)*/
    virtual bool sensShift(VolAbsoluteShift* shift);

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    virtual PDFCalculator* getPDFCalculator(
        const PDFRequest* request,
        const CAsset*     asset) const;

private:
    FlatVol();
    FlatVol(const FlatVol &rhs);
    FlatVol& operator=(const FlatVol& rhs);
    bool scalarShift(double);
    void scalarRestore(double);
    static void acceptValueDateCollector(const FlatVol* flatVol, 
                                         CValueDateCollector* collector);

    string       name;    // name of the vol
    double       flatVol;
    DateTime     baseDate;
    TimeMetricSP timeMetric;
};

typedef smartConstPtr<FlatVol> FlatVolConstSP;
typedef smartPtr<FlatVol> FlatVolSP;

DRLIB_END_NAMESPACE
#endif
