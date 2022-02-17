/**
 * @file Restorable.hpp
 */

#ifndef EDG_RISKMGR_IRESTORABLEWITHRESPECTTO_H
#define EDG_RISKMGR_IRESTORABLEWITHRESPECTTO_H

#include "edginc/ITweakableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface designating market objects which have a particular "risk property"
 * and can undo changes to it.
 *
 * See ITweakableWithRespectTo<TAG> for the rest of the methods,
 * IRiskQuantityFactory for an overview of the "declarative" sensitivities
 * framework.
 */

template <class TAG>
class IRestorableWithRespectTo:
        public virtual ITweakableWithRespectTo<TAG> {

    static void load(CClassSP &clazz) {
        REGISTER_INTERFACE(IRestorableWithRespectTo, clazz);
        EXTENDS(ITweakableWithRespectTo<TAG>);
    }

public:

    static CClassConstSP const TYPE;

    /**
     * Undo the effect of the last sensShift().
     *
     * @a property is the identical with the argument given to the last
     * sensShift().
     */

    virtual void sensRestore(const PropertyTweak<TAG>& property) = 0;
private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each array template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class TAG> CClassConstSP const 
IRestorableWithRespectTo<TAG>::TYPE =
CClass::templateRegisterClass(typeid(IRestorableWithRespectTo<TAG>));
#endif

DRLIB_END_NAMESPACE

#endif
