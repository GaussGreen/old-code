//////////////////////////////////////////////////////////////////////////

#ifndef EDR_SENS_CONTROL_ALL_NAMES_HPP
#define EDR_SENS_CONTROL_ALL_NAMES_HPP
#include "edginc/SensControl.hpp"
#include "edginc/ScenarioShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Slightly more specialised SensControl where the shifting is done on a
    per name basis */
class RISKMGR_DLL SensControlAllNames: public SensControl,
                           public virtual IScenarioShift,
                           public virtual IPerturbation /* IPerturbation for 
                                                           backward
                                                           compatibility only */
{
public:    
    friend class ScalarShiftHelper;
    static CClassConstSP const TYPE;

    virtual ~SensControlAllNames();

    /** Returns null (=> all names match)  */
    virtual ITweakNameResolver* nameResolver();

    /** apply this scenario shift to the supplied object */
    virtual bool applyScenario(IObjectSP object);
    virtual bool preapplyScenario(IObjectSP object);

    /** Returns null as this tweak shifts all instances rather than
        specific named types */
    virtual OutputNameConstSP getMarketDataName() const;

    /** IPerturbation implementation - for backwards compatibility only.
        Equivalent here to applyScenario(objectToShift) */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name);

protected:
    /** Note SensControl is abstract. Create a sens control of type clazz
        and which uses outputName (eg VEGA_PARALLEL) to identify results */
    SensControlAllNames(const CClassConstSP& clazz,
                        const string&        outputName);

private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif
