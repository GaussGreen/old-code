
#ifndef EDG_SENS_CONTROL_H
#define EDG_SENS_CONTROL_H
#include "edginc/Results_forward.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/TweakOptID.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class CInstrument; // to avoid circular header dependencies
class IModel;
class Control;
FORWARD_DECLARE(OutputName);
FORWARD_DECLARE(OutputRequest);

/** An abstract base class used for explicitly driving the tweaking
    process (cf Delta Next Day). It identifies a particular type of
    tweak to make. It is also responsible for driving how that tweak
    is calculated as well as how the result is returned.  */
class RISKMGR_DLL SensControl: public Sensitivity,
                   public virtual ITweakOptID{
public:
    friend class SensControlHelper;
    static CClassConstSP const TYPE;

    virtual ~SensControl();

    /** identifies the name used storing associated results in the output */
    virtual const string& getSensOutputName() const;

    /** returns the name identifying the market data to be shifted. Returns
        null if not set or if this tweak shifts all instances rather than
        specific named types */
    virtual OutputNameConstSP getMarketDataName() const = 0;

    /** stores an the intitial value of the quantity being shifted. This is
        very useful in calculating greeks like delta where the original value
        of what was being shifted is needed. It avoid the need for a 
        separate suite of interfaces allowing the retrieval of spot prices,
        correlations, fx spots etc. */
    void setInitialValue(double val);

    /** gets the intitial value of the quantify being shifted. See the 
        comments under {@link #initialValueSet(double)} */
    double getInitialValue() const;

    /** Used by objects implementing restorable tweaks to store data needed
        to restore themselves after the tweak. Stores given data on
        the SensControl */
    void setRestoreData(const IObjectSP& restoreData);

    /** resets internally stored values associated with tweaking. This is 
        because this type of
        object (or more precisely derived types) can store information about
        what has been tweaked (eg initial values) */
    void reset();

    /** Used by objects implementing restorable tweaks to store data needed
        to restore themselves after the tweak. Retrieves the data previously
        stored on the SensControl. The data stored is then cleared */
    IObjectSP getRestoreData();

    /** Once used to make a shift, this reports the appropriate divisor
        for this sensitivity */
    virtual double divisor() const = 0;

    /** Implemented by restorableShiftInterface()->isInstance(obj) */
    virtual bool restorableShift(IObjectConstSP obj) const;

    /** returns the interface identifying what an object has to do in order
        to be support the tweak that this object represents which is also
        restorable */
    virtual CClassConstSP restorableShiftInterface() const = 0;

    /** Like calculateSens in terms of inputs but calculates a single
        first order one-sided scalar derivative. The name of what to
        be tweaked should be set in 'this' */
    double calculateOneSidedFirstDeriv(IModel*          algorithm,
                                       CInstrument*     instrument,
                                       Control*         control,
                                       Results*         results);
protected:
    /** Note SensControl is abstract. Create a sens control of type clazz
        and which uses outputName (eg VEGA_PARALLEL) to identify results */
    SensControl(CClassConstSP clazz, const string&        outputName);

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*     tweakGroup,
                           CResults*       results) = 0;

    /** shifts tweakGroup/market data by given shift, reprices, and
        restores data to original. Returns shifted price */
    double shiftAndPrice(TweakGroup*      tweakGroup,
                         double           origPrice);

    /** shifts tweakGroup/market data by given shift, reprices, and
        restores data to original. Returns shifted price. If forceReprice
        is true then price will be recalculated regardless of whether
        tweakGroup is actually shifted during the tweak */
    double shiftAndPrice(TweakGroup*      tweakGroup,
                         double           origPrice,
                         bool             forceReprice);


    /** shifts tweakGroup by this shift and then by given
        shift, reprices, and restores data to original. Returns shifted
        price */
    double shiftAndPrice(SensControl*     secondShift,
                         TweakGroup*      tweakGroup,
                         double           origPrice);

    /** Applies this shift and then does control->calculate() where the
        control is constructed using sens and outReqs, either or both 
        of which can be null */
    ResultsSP shiftAndCalculate(TweakGroup*                      tweakGroup,
                                const SensitivityArrayConstSP&   sens,
                                const OutputRequestArrayConstSP& outReqs);

    /** Calculates 'price' for this sensitivity */
    double calcSensPrice(TweakGroup*  tweakGroup);

    /** calculates a first order one-sided scalar derivative */
    double calcOneSidedFirstDeriv(TweakGroup*    tweakGroup,
                                  CResults*      results);
private:
    string     outputName; // $unregistered
    double     initialValue;
    bool       initialValSet;
    // data needed to restore after a shift
    IObjectSP  restoreData;
    // not implemented
    SensControl();
    SensControl(const SensControl &rhs);
    SensControl& operator=(const SensControl& rhs);
};

typedef SensControl                CSensControl;
typedef smartPtr<SensControl>      SensControlSP;
typedef smartConstPtr<SensControl> SensControlConstSP;
typedef array<SensControlSP, SensControl> SensControlArray;
typedef smartPtr<SensControlArray> SensControlArraySP;
typedef smartConstPtr<SensControlArray> SensControlArrayConstSP;
#ifndef QLIB_SENSCONTROL_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<SensControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<SensControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<SensControlSP _COMMA_ SensControl>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<SensControlArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<SensControlArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<SensControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<SensControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<SensControlSP _COMMA_ SensControl>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<SensControlArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<SensControlArray>);
#endif

DRLIB_END_NAMESPACE

#endif

