//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRProductsLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//----------------------------------------------------------------------------
//
//
// This file is not itself compiled and linked into QLib.  In order to support
// selective builds*, the system generates a modified version
// IRProductsLib-filtered.cpp, in which references to "FooLoad()" are removed
// unless Foo.cpp is selected for inclusion.  It's IRProductsLib-filtered.cpp
// which is actually built into the library.
//
// (NB don't edit IRProductsLib-filtered.cpp, since any changes you make to it
// will be overwritten --- edit this file.)
//
// [Technical note: the filtering is performed by
// ../../../makerules/scripts/filterProductsLib.pl, invoked from
// ../../../makerules/gnu/selective-srcs.mkh and from
// ../../QLibSolution.vsmproj:QLibPartial.checkDeps().]
//
// ---
//   *see selectivebuild.example for more info

#include "edginc/config.hpp"
#include "edginc/IRProductsLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE
extern bool CallableTurboLoad();
extern bool CallableTurbo2Load();
extern bool TurboSwapLoad();
extern bool KOptionLoad();
extern bool KComponentLoad();
extern bool KSimpleCouponLoad();
extern bool KRibLoad();
extern bool CallableRibLoad();
extern bool CallableRib2Load();
extern bool KSumLoad();
extern bool KTurboLoad();
extern bool KFloatLeg2Load();
extern bool KKnockOutLoad();
extern bool KFloatLegLoad();
extern bool KRibFloatLegLoad();
extern bool CouponSchedLoad();
extern bool OptionSchedLoad();
extern bool KRiskyComponentLoad();
extern bool ContingentIRCouponLoad();
extern bool KProtecionLegLoad();
extern bool KTarnLoad();
extern bool KKOSwapLoad();
extern bool CallableKOSwapLoad();
extern bool KAccumulatedCFLoad();
extern bool KKOOptionLoad();
extern bool KRibTurboLoad();
extern bool KCashflowLoad();
extern bool OptionGridLoad();
extern bool KForwardOptionLoad();
extern bool QuasiVanillaLoad();
extern bool SimpleRibLoad();
extern bool VnfmTestLoad();
extern bool VnfmTestTimeLoad();
extern bool VolToolsLoad();
extern bool SharedEnumsLoad();
extern bool KBarrierLoad();
extern bool KKnockIOLoad();
extern bool KBoolExprLoad();

void IRProductsLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker 
       includes all symbols out of the market directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing a product. For example Vanilla might be dropped since
       the only references might be throught the abstract parent class. */


    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        SharedEnumsLoad() &&
        CallableTurboLoad() &&
        CallableTurbo2Load() &&
        TurboSwapLoad() &&
        KOptionLoad() &&
        KComponentLoad() &&
        KSimpleCouponLoad() &&
        KRibLoad() &&
        CallableRibLoad() &&
        CallableRib2Load() &&
        KSumLoad() &&
        KTurboLoad() &&
        KFloatLeg2Load() &&
        KKnockOutLoad() &&
        KFloatLegLoad() &&
        KRibFloatLegLoad() &&
        CouponSchedLoad() &&
        OptionSchedLoad() &&
        KTarnLoad()       &&
        KKOSwapLoad()     &&
        ContingentIRCouponLoad() &&
        CallableKOSwapLoad() &&
        KAccumulatedCFLoad() &&
        KKOOptionLoad() &&
        KRibTurboLoad() &&
        KCashflowLoad() &&
        OptionGridLoad()  &&
        KForwardOptionLoad()  &&
        SimpleRibLoad() &&
		QuasiVanillaLoad() &&
		VnfmTestLoad() &&
		VnfmTestTimeLoad() &&
		VolToolsLoad() &&
		KBarrierLoad() &&
		KKnockIOLoad() &&
		KBoolExprLoad() &&
        true;

    if (!success){
        throw ModelException("IRProductsLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
