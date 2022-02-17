//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : ParSpreadUpfrontParallelTP.hpp
//
//   Description : ParSpread curve upfront parallel shift Tweak Property
//
//   Date        : Nov 2007
//
//----------------------------------------------------------------------------

#ifndef QLIB_PARSPREADUPFRONTPARALLELTP_H
#define QLIB_PARSPREADUPFRONTPARALLELTP_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "ParSpread curve upfront parallel shift" as a property of market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 * See ParSpreadParallel for a fully documented analogous example.
 */

struct RISKMGR_DLL ParSpreadUpfrontParallelTP: CObject {
    static CClassConstSP const TYPE;
    ParSpreadUpfrontParallelTP(); ~ParSpreadUpfrontParallelTP();
    typedef Void Qualifier;
    enum {discrete = 0};
};

DRLIB_END_NAMESPACE

#endif
