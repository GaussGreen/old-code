
#ifndef QRD_SCALAR_PER_NAME_SHIFT_HPP
#define QRD_SCALAR_PER_NAME_SHIFT_HPP
#include "edginc/IHypothesis.hpp"
#include "edginc/ScalarTweak.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(OutputName)

/** Interface defining 'scalar' shifts (can indeed by 'set' type
    shifts too) which act on a specific name and the list of relevant
    names can be retured. Needed by ImpliedScalarShift which wants to use
    scenarios as well as sensitivities to root finding */
class RISKMGR_DLL IScalarPerNameShift: public virtual IObject,
                           public virtual IScalarTweak{
public:    
    static CClassConstSP const TYPE;

    IScalarPerNameShift();

    virtual ~IScalarPerNameShift();

    virtual IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                                    double shiftSize,
                                                    IObjectSP world) = 0;

    virtual OutputNameArrayConstSP allNames(const IObject* object) const = 0;

private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IScalarPerNameShift>      IScalarPerNameShiftSP;
typedef smartConstPtr<IScalarPerNameShift> IScalarPerNameShiftConstSP;
DRLIB_END_NAMESPACE

#endif
