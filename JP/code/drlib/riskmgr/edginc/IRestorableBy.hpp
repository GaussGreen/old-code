/**
 * @file IRestorableBy.hpp
 */

#ifndef DRLIB_IRestorableBy_H
#define DRLIB_IRestorableBy_H

#include "edginc/ITweakableBy.hpp"

DRLIB_BEGIN_NAMESPACE

template <class TWEAK>
class IRestorableBy: public ITweakableBy<TWEAK> {

public:

    static CClassConstSP const TYPE;

    virtual void restore(TWEAK *shift) = 0;
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
