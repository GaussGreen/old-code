//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RadarLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//----------------------------------------------------------------------------
//
//
// This file is not itself compiled and linked into QLib.  In order to support
// selective builds*, the system generates a modified version
// RadarLib-filtered.cpp, in which references to "FooLoad()" are removed
// unless Foo.cpp is selected for inclusion.  It's RadarLib-filtered.cpp
// which is actually built into the library.
//
// (NB don't edit RadarLib-filtered.cpp, since any changes you make to it
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
#include "edginc/RadarLib.hpp"
#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE
extern bool RadarRepWrapperLoad();
extern bool RadarRepAddinLoad();
extern bool PolyBasisWrapperLoad();
extern bool IFuncBasisWrapperLoad();
extern bool IFittingVarTransWrapperLoad();
extern bool IdentityTransWrapperLoad();
extern bool KRadarRepGeneratorLoad();

void RadarLib::linkInClasses(){
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
		PolyBasisWrapperLoad() &&
		IFuncBasisWrapperLoad() &&
        IFittingVarTransWrapperLoad() &&
        IdentityTransWrapperLoad() &&
		RadarRepWrapperLoad() &&
		RadarRepAddinLoad() &&
        KRadarRepGeneratorLoad() &&
        true;

    if (!success){
        throw ModelException("RadarLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
