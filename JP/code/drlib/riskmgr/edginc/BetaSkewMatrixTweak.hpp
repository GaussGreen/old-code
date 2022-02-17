//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BetaSkewMatrixTweak.hpp
//
//   Description : Tweak of each relevant point inside a skew surface
//                 This is essentially the same as BetaSkewPointwiseTweak but
//                 changed to store the array of results inside an object
//                 rather than storing the actual array.
//                 BetaSkewPointwiseTweak is now deprecated.
//
//   Author      : Jose Hilera
//
//   Cautions    : The sensitivity implemented here is derived from 
//                 BetaSkewPointwiseTweak, which is to be removed. To avoid 
//                 duplication, BetaSkewPointwiseTweak::IShift and 
//                 BetaSkewPointwiseTweak::IRestorableShift are still used (and,
//                 indirectly, also BetaSkewPointwiseTweak::ISensitivePoints 
//                 and BetaSkewPointwiseTweak::ISensitivePointsImt) - they 
//                 should be moved into this class when BetaSkewPointwiseTweak 
//                 is finally removed, along the methods in BetaSkewPointwiseTweak
//                 which this class is currently inheriting (NOT calculate, which
//                 BetaSkewMatrixTweak overrides!)
//
//   Date        : July 2005
//
//----------------------------------------------------------------------------

#ifndef QLIB_BETA_SKEW_MATRIX_TWEAK_H
#define QLIB_BETA_SKEW_MATRIX_TWEAK_H

#include "edginc/ScalarShift.hpp"
#include "edginc/BetaSkewGridPoint.hpp"
#include "edginc/Additive.hpp"
#include "edginc/BetaSkewPointwiseTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/** Tweak of each relevant point inside a skew surface */
class RISKMGR_DLL BetaSkewMatrixTweak: public BetaSkewPointwiseTweak,
                           virtual public Additive {
public:    
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** implements a one sided vector derivative for each instance of the
        market data which is sensitive to this SensControl */
    void calculate(TweakGroup* tweakGroup,
                   CResults*   results);

    /** constructor with explicit shift size */
    BetaSkewMatrixTweak(double shiftSize);

    /** constructor with explicit shift size and tweakAll*/
    BetaSkewMatrixTweak(const double shiftSize,
                        const bool   tweakAll,
                        const int    numberOfNeighbours);

private:
    // for reflection
    BetaSkewMatrixTweak();
    static IObject* defaultBetaSkewMatrixTweak();
    static void load(CClassSP& clazz);

    BetaSkewMatrixTweak(const BetaSkewMatrixTweak& rhs);
    BetaSkewMatrixTweak& operator=(const BetaSkewMatrixTweak& rhs);
};

DRLIB_END_NAMESPACE

#endif //QLIB_BETA_SKEW_MATRIX_TWEAK_H
