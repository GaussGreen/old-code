//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CEVJ.hpp
//
//   Description : CEV+jump volatility
//
//   Date        : 30 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef CEVJ_HPP
#define CEVJ_HPP
#include "edginc/VolBase.hpp"
#include "edginc/VolLevel.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/TimeMetric.hpp"
#include "edginc/Addin.hpp"
#include "edginc/ValueDateCollector.hpp"
#include "edginc/VolParallel.hpp"
#include "edginc/VolPointwise.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"
#include "edginc/JumpRateParallel.hpp"
#include "edginc/JumpRatePointwise.hpp"
#include "edginc/CEVPowerParallel.hpp"
#include "edginc/CEVPowerPointwise.hpp"
#include "edginc/VolParallelShift.hpp"
#include "edginc/Theta.hpp"
#include "edginc/VolAbsoluteShift.hpp"
#include "edginc/VolBenchmarkShift.hpp"
#include "edginc/PowerVega.hpp"

DRLIB_BEGIN_NAMESPACE

class CEVJProcessed;

/** CEV+jump volatility data class */
class MARKET_DLL CEVJ: public CVolBase,
            virtual public ITweakableWithRespectTo<VolParallel>,
            virtual public ITweakableWithRespectTo<VolPointwise>,
            virtual public JumpRateParallel::IShift,
            virtual public JumpRatePointwise::IShift,
            virtual public CEVPowerParallel::IShift,
            virtual public CEVPowerPointwise::IShift,
            virtual public VolLevel::Shift,
            virtual public VolParallelShift::Shift,
            virtual public Theta::IShift,
            virtual public VolAbsoluteShift::IShift,
            virtual public VolBenchmarkShift::Shift,
            virtual public PowerVega::Shift
{
public:
    static CClassConstSP const TYPE;
    friend class CEVJHelper;
    friend class CEVJAddin;
    friend class VolProcessedCEVJ;
    friend class CTree1fCEVJCalib;

    virtual void validatePop2Object();

    /** Returns name of vol */
    virtual string getName() const;

    /** populate from market cache */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** Combines market and instrument data together to give a
        Processed Vol */
    virtual CVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                           const CAsset*      asset) const;

    // not implemented
    virtual CVolProcessed* getProcessedVol(
        const CVolRequest* volRequest,
        const CAsset*      eqAsset,
        const FXAsset*     fxAsset,
        const Correlation* eqFXCorr) const;
        
    bool sensShift(Theta* shift);
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

    /** Implements VegaParallel  */
    virtual string sensName(const VolParallel*) const;
    /** Shifts the object using given shift*/
    virtual TweakOutcome sensShift(const PropertyTweak<VolParallel>& tweak);

    /** Implements JumpRateParallel  */
    virtual string sensName(JumpRateParallel* shift) const;
    /** Shifts the object using given shift*/
    virtual bool sensShift(JumpRateParallel* shift);

    /** Implements CEVPowerParallel  */
    virtual string sensName(CEVPowerParallel* shift) const;
    /** Shifts the object using given shift*/
    virtual bool sensShift(CEVPowerParallel* shift);

    /** Implements VegaPointwise */
    virtual string sensName(const VolPointwise*) const;
    /** Shifts VegaPointwise */
    virtual TweakOutcome sensShift(const PropertyTweak<VolPointwise>&);
    /** get expiries for vol */
    ExpiryWindowArrayConstSP sensQualifiers(const VolPointwise*) const;

    /** Implements JumpRatePointwise */
    virtual string sensName(JumpRatePointwise* shift) const;
    /** Shifts JumpRatePointwise */
    virtual bool sensShift(JumpRatePointwise* shift);
    /** get expiries for jump rate */
    ExpiryArrayConstSP sensExpiries(JumpRatePointwise* shift) const;

    /** Implements CEVPowerPointwise */
    virtual string sensName(CEVPowerPointwise* shift) const;
    /** Shifts CEVPowerPointwise */
    virtual bool sensShift(CEVPowerPointwise* shift);
    /** get expiries for CEVPower */
    ExpiryArrayConstSP sensExpiries(CEVPowerPointwise* shift) const;

    /** Implements VolAbsoluteShift scenario */
    /** Returns name identifying this object for VolAbsoluteShift */
    virtual string sensName(VolAbsoluteShift* shift) const;
    /** Shifts the object using given shift (see VolAbsoluteShift::IShift)*/
    virtual bool sensShift(VolAbsoluteShift* shift);

    /** Implements VolBenchmarkShift scenario */
    /** Returns name identifying this object for VolBenchmarkShift */
    virtual string sensName(VolBenchmarkShift* shift) const;
    /** Shifts the object using given shift (see VolBenchmarkShift::Shift)*/
    virtual bool sensShift(VolBenchmarkShift* shift);

    /** Implements PowerVega scenario */
    /** Returns name identifying this object for PowerVega */
    virtual string sensName(PowerVega* shift) const;
    /** Shifts the object using given shift (see PowerVega::Shift)*/
    virtual bool sensShift(PowerVega* shift);

private:
    static void acceptValueDateCollector(const CEVJ* cevj,
                                         CValueDateCollector* collector);

    friend class CEVJProcessed;
    CEVJ();
    CEVJ(const CEVJ &rhs);
    CEVJ& operator=(const CEVJ& rhs);

    string          name;    // name of the vol
    DateTime        baseDate;
    TimeMetricSP    timeMetric;

    int		    VolMode; // 0=log-normal vol, 1=CEV+jump, 2=Lambda+jump
    // CEV and jump data or Lambda diffusion and jump data
    DoubleArray     ATMVolArr; // term structure of ATM vol values
    ExpiryArraySP   ATMVolBM; // term structure diffusion vol bench mark labels
    StringArray     BenchMarkStrg; // String Array for BenchMarks.
    DoubleArray     CEVPowerArr; // CEV or Lambda, this must be either 1 element or the same number as ATMVolBM
    DoubleArray     JumpWidthArr; // jump width
    DoubleArray     JumpMeanArr; // jump mean, this must be either 1 element or the same number as ATMVolBM
    DoubleArray     JumpRateArr;	// jump rate, this must be either 1 element or the same number as ATMVolBM
    int		    JumpMode; // 0=proportional(default), else=absolute
    double	    SpotRef; // spot reference, for scaling variance so that the unit of diff vol is the usual %

    int             DEBUG_fix2ndMoment; // 0 = use aditional 2nd Var term.  1= simple varriance. others are same as Akasaka.

};

typedef smartConstPtr<CEVJ> CEVJConstSP;
typedef smartPtr<CEVJ> CEVJSP;

DRLIB_END_NAMESPACE
#endif
