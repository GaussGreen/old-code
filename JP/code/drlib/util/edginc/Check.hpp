//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Check.hpp
//
//   Description : Various useful checking functions etc
//
//   Date        : Oct 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDG_CHECK_HPP
#define EDG_CHECK_HPP

#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE

/** Various useful checking functions */

class UTIL_DLL Check {

public:

    /** Ensures that weights are positive, that weights.size() == numWeights
        and that the weights sum to 1.0 to within reasonable
        tolerance. The difference between 1.0 and the sum of the weights
        is then distributed back among the weights. */
    static void percWeights(DoubleArray&       weights, // note not const
                            int                numToMatch, // how many weights are required
                            string             desc); // description of what the weights are aligned with e.g. "assets"
        
};

DRLIB_END_NAMESPACE

#endif

