/**
 * @file SimpleTweakNameListID.hpp
 */

#ifndef QLIB_SimpleTweakNameListID_H
#define QLIB_SimpleTweakNameListID_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/TweakNameListID.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(SimpleTweakNameListID)

/**
 * An ITweakNameListID for subclasses of MarketObject which simply report their
 * own name
 */

class RISKMGR_DLL SimpleTweakNameListID: public CObject,
                             public virtual ITweakNameListID {

    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    CClassConstSP subjectType; // $unregistered

public:

    SimpleTweakNameListID(CClassConstSP subjectType);
    ~SimpleTweakNameListID();

    /**
     * Returns @a subjectType supplied to SimpleTweakNameListID() constructor
     */

    CClassConstSP shiftInterface() const;

    /**
     * Casts @a obj to MarketObject and adds its
     * MarketObject::getName() to @a namesList
     */

    void appendName(OutputNameArray& namesList, IObjectConstSP obj);
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
