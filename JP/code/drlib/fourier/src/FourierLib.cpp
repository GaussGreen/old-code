//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Regis Guichard
//
//   Date        : 09 March 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FourierLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE

// to avoid header file
extern bool FourierProcessSVLoad();
extern bool FourierProcessSVJLoad();
extern bool FourierProcessSVJJLoad();
extern bool FourierProcessGammaOULoad();
extern bool FourierProcessIGOULoad();
extern bool FourierProcessCGMYHestonLoad();
extern bool FourierProcessAJDSuperLoad();
extern bool FourierProcessMertonLoad();
extern bool FourierProcessSVCJLoad();
extern bool FourierProcessVSCurveLoad();
extern bool PDFFourierLoad();

void FourierLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker 
       includes all symbols out of the market directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing volatility. For example FlatVol might be dropped since
       the only references might be throught the abstract parent class. */


    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        FourierProcessSVLoad() &&
        FourierProcessSVJLoad() &&
        FourierProcessSVJJLoad() &&
        FourierProcessGammaOULoad() &&
        FourierProcessIGOULoad() &&
        FourierProcessCGMYHestonLoad() &&
        FourierProcessAJDSuperLoad() &&
        FourierProcessSVCJLoad() &&
        FourierProcessMertonLoad() &&
        FourierProcessVSCurveLoad() &&
        PDFFourierLoad() &&
        true;

    if (!success){
        throw ModelException("FourierLib::registerClasses",
                             "Registration error");
    }

}

DRLIB_END_NAMESPACE
