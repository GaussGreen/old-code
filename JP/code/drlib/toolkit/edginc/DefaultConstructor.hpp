/**
 * @file DefaultConstructor.hpp
 */

#ifndef QLIB_DefaultConstructor_H
#define QLIB_DefaultConstructor_H

DRLIB_BEGIN_NAMESPACE

class IObject;

template <class T>
struct DefaultConstructor {
    static IObject* iObject() { return new T(); }
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
