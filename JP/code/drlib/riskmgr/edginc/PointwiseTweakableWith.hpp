/**
 * \file PointwiseTweakableWith.hpp
 *
 *
 * $History$
 *
 */

#ifndef EDG_RISKMGR_POINTWISETWEAKABLEWITHRESPECTTO_H
#define EDG_RISKMGR_POINTWISETWEAKABLEWITHRESPECTTO_H

#include "edginc/Class.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Interface for (market data) objects which can be subjected to a tweak of a
 * certain type
 *
 * See GenericVectorShift and GenericScalarShift for info about the purpose of
 * this interface.
 *
 * For example: ParSpreadCurve implements
 * PointwiseTweakableWith<ParSpreadRhoPointwiseTweak>.
 */

template <class TWEAK>
struct PointwiseTweakableWith: public virtual IObject {
    static CClassConstSP const TYPE;

    virtual bool sensShift(TWEAK *shift) = 0;
    virtual string sensName(TWEAK *shift) const = 0;
    virtual ExpiryArrayConstSP sensExpiries(TWEAK* shift) const = 0;

    virtual ~PointwiseTweakableWith() {}

    static void load(CClassSP &clazz) {
        REGISTER_INTERFACE(PointwiseTweakableWith, clazz);
        EXTENDS(IObject);
    }
};

DRLIB_END_NAMESPACE

#endif
