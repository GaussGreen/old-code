/**
 * @file BoxedInt.hpp
 */

#ifndef QLIB_BoxedInt_H
#define QLIB_BoxedInt_H

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(BoxedInt)

/**
 * An int wrapped up as a CObject (boxed-array version)
 *
 * CInt is QLib's primary representation of an integer object; and its array
 * type is array<int, int>, i.e. an array of actual (unboxed) ints, because
 * that's convenient and efficient for most purposes.  However, for abstract
 * code which wants to treat integers uniformly with genuine objects, it's
 * handy to have a class whose array type is the more consistent
 * array<BoxedIntSP, BoxedInt> of smart pointers.
 */

class TOOLKIT_DLL BoxedInt: public CObject {

    static IObject* newOne();
    static void load(CClassSP& clazz);

public:

    static CClassConstSP const TYPE;

private:

    const int i;
    BoxedInt(int i);
    BoxedInt(const BoxedInt& rhs);
    BoxedInt& operator=(const BoxedInt& rhs);

public:

    static BoxedInt* create(int i);
    static BoxedIntSP SP(int i);
    ~BoxedInt();

    virtual IObject* clone() const;

    int intValue() const;

    string toString() const;
};

DRLIB_END_NAMESPACE

#endif
