//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : PrepayParallelTP.hpp
//
//   Description : ABCDS prepay curve horizontal parallel shift Tweak Property
//
//   Date        : May 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_PREPAYPARALLELTP_H
#define QLIB_PREPAYPARALLELTP_H

#include "edginc/Void.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * Tag class denoting "prepay curve horizontal shift in datetime" as a property of market names.
 *
 * This class doesn't have any methods or members --- it just acts as a "tag"
 * for instantiating a family of templates,
 *
 * See ParSpreadParallel for a fully documented analogous example.
 */

struct RISKMGR_DLL PrepayParallelTP: CObject {
    static CClassConstSP const TYPE;
    PrepayParallelTP(); ~PrepayParallelTP();
    typedef Void Qualifier;
    enum {discrete = 0};
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
