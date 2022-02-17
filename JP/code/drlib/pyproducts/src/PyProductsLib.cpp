//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PyProductsLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//----------------------------------------------------------------------------
//
//
// This file is not itself compiled and linked into QLib.  In order to support
// selective builds*, the system generates a modified version
// PYProductsLib-filtered.cpp, in which references to "FooLoad()" are removed
// unless Foo.cpp is selected for inclusion.  It's PYProductsLib-filtered.cpp
// which is actually built into the library.
//
// (NB don't edit PYProductsLib-filtered.cpp, since any changes you make to it
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
#include "edginc/PYProductsLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE

extern bool KPythonLoad();
extern bool MCPythonLoad();

void PyProductsLib::linkInClasses()
{
    bool success =
        KPythonLoad() &&
        MCPythonLoad() &&
        true ;

    if (!success){
        throw ModelException("PYProductsLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
