//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ScenarioShift.hpp
//
//   Description : Defines a scenario shift to market data
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------


#ifndef SCENARIOSHIFT_HPP
#define SCENARIOSHIFT_HPP
#include "edginc/Object.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE
/** Interface defining a scenario shift to market data */
class RISKMGR_DLL IScenarioShift: public virtual IObject {
public:
    static CClassConstSP const TYPE;

    virtual ~IScenarioShift();

    /** 'preapply' this scenario shift to the supplied object. This
        method is invoked before the market data is retrieved for the
        instrument and model. It can be used, for example, to do a
        name substitution. Most implementations will just do "return false;"
        This is, for clarity, a non-restorable
        shift i.e there is no mechanism for undoing the shift. The
        return value indicates if anything was actually shifted (true
        => yes) */
    virtual bool preapplyScenario(IObjectSP object) = 0;

    /** apply this scenario shift to the supplied object. This is, for
        clarity, a non-restorable shift ie there is no mechanism for
        undoing the shift. The return value indicates if anything was
        actually shifted (true => yes) */
    virtual bool applyScenario(IObjectSP object) = 0;
  
private:
    static void load(CClassSP& clazz);
};

/** The obvious implementation for IScenarioShift - an object that can do
    the shift + the name of what to shift */
class RISKMGR_DLL ScenarioShift: public CObject,
                     public virtual IScenarioShift {
public:
    static CClassConstSP const TYPE;

    virtual ~ScenarioShift();

    ScenarioShift(IPerturbationSP sensCtrl, const string& marketDataName);

    /** Uses the contained IPerturbation in conjunction with the marketDataName
        to find and shift the required piece of data. The return
        value indicates if anything was actually shifted (true => yes) */
    virtual bool applyScenario(IObjectSP object);

    // Nothing to do before market data is retrieved.
    virtual bool preapplyScenario(IObjectSP object);

private:
    friend class ScenarioShiftHelper;
    ScenarioShift(); // for reflection
    ScenarioShift(const ScenarioShift &rhs);
    ScenarioShift& operator=(const ScenarioShift& rhs);

    IPerturbationSP sensCtrl;
    string          marketDataName;

};

// typedef for smart pointers to ScenarioShifts
typedef smartPtr<ScenarioShift> ScenarioShiftSP;
typedef array<ScenarioShiftSP, ScenarioShift> ScenarioShiftArray;

typedef smartPtr<IScenarioShift> IScenarioShiftSP;
typedef array<IScenarioShiftSP, IScenarioShift> IScenarioShiftArray;
typedef smartPtr<IScenarioShiftArray> IScenarioShiftArraySP;
typedef smartConstPtr<IScenarioShiftArray> IScenarioShiftArrayConstSP;


DRLIB_END_NAMESPACE
#endif
