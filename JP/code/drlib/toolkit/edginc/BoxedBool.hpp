/**
 * @file BoxedBool.hpp
 */

#ifndef QLIB_BoxedBool_H
#define QLIB_BoxedBool_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(BoxedBool)

/**
 * A bool wrapped up as a CObject (boxed-array version)
 *
 * CBool is QLib's primary representation of an integer object; and its array
 * type is array<bool, bool>, i.e. an array of actual (unboxed) bools, because
 * that's convenient and efficient for most purposes.  However, for abstract
 * code which wants to treat integers uniformly with genuine objects, it's
 * handy to have a class whose array type is the more consistent
 * array<BoxedBoolSP, BoxedBool> of smart pointers.
 */

class TOOLKIT_DLL BoxedBool: public CObject {

    static IObject* newOne();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    const bool b;
    BoxedBool(bool b);
    BoxedBool(const BoxedBool& rhs);
    BoxedBool& operator=(const BoxedBool& rhs);

public:

    static BoxedBool* create(bool b);
    static BoxedBoolSP SP(bool b);
    ~BoxedBool();

    virtual IObject* clone() const;

    bool boolValue() const;

    string toString() const;
};

DRLIB_END_NAMESPACE

#endif
