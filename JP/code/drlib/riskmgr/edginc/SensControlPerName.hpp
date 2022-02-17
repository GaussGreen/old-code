//////////////////////////////////////////////////////////////////////////

#ifndef EDR_SENS_CONTROL_PER_NAME
#define EDR_SENS_CONTROL_PER_NAME

#include "edginc/SensControl.hpp"
#include "edginc/IPerNameSensitivity.hpp"
#include "edginc/TweakNameListID.hpp"
#include "edginc/TweakNameResolver.hpp"
#include "edginc/Perturbation.hpp"

DRLIB_BEGIN_NAMESPACE

/** Slightly more specialised SensControl where the shifting is done on a
    per name basis */
class RISKMGR_DLL SensControlPerName: public SensControl,
                          public virtual IPerturbation,
                          public virtual IPerNameSensitivity,
                          public virtual ITweakNameListID,
                          public virtual ITweakNameResolver{
public:    
    friend class ScalarShiftHelper;
    static CClassConstSP const TYPE;

    ~SensControlPerName();

    /** Returns array of output names which need to be tweaked for this
        sensitivity. In particular, if hasOverrideNames is true then returns
        overrideNames() else generates a list based on object supplied (which
        is typically, but needn't be, a TweakGroup */
    virtual OutputNameArrayConstSP names(const IObject* tweakGroup) const;

    /** The same as the above method except that for the case where
        this has its own names then for each of those names which is
        not present NotApplicable is stored in the Results. */
    virtual OutputNameArrayConstSP names(const IObject* tweakGroup,
                                         Results*       results);

    virtual OutputNameArrayConstSP allNames(const IObject* object) const;

    /** apply this Perturbation to the object with the specified name contained
        within the supplied object. This is, for clarity, a non-restorable
        shift ie there is no mechanism for undoing the shift. Returns true
        if something actually was shifted */
    virtual bool findAndShift(IObjectSP         objectToShift, 
                              OutputNameConstSP name);

    /** Returns this => tweaks defined by classes derived from this type are
        associated with individual names */
    virtual ITweakNameResolver* nameResolver();

    /** returns the name identifying the market data to be shifted. Returns
        null if not set */
    virtual OutputNameConstSP getMarketDataName() const;

    /** sets the name identifying the market data to be shifted */
    virtual void setMarketDataName(OutputNameConstSP name);

    /** does the given name match the name identifying the market data
        to be shifted */
    bool marketDataNameMatches(OutputNameConstSP name) const;


protected:
    /** Note SensControl is abstract. Create a sens control of type clazz
        and which uses outputName (eg VEGA_PARALLEL) to identify results */
    SensControlPerName(const CClassConstSP& clazz,
                       const string&        outputName);

private:
    OutputNameConstSP marketDataName;
    static void load(CClassSP& clazz);
};

typedef smartPtr<SensControlPerName> SensControlPerNameSP;

DRLIB_END_NAMESPACE

#endif
