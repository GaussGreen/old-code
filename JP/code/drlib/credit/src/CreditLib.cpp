//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Jay Blumenstein
//
//   Date        : 22 Aug 2002
//
//

#include "edginc/config.hpp"
#include "edginc/CreditLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/CreditEngine.hpp" // FILTER
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE

extern bool CreditAddinLoad();

void CCreditLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker 
       includes all symbols out of the credit directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing volatility. For example FlatVol might be dropped since
       the only references might be throught the abstract parent class. */


    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        CreditPathValuesIn::TYPE &&
		CreditEngine::TYPE &&
		CreditPathValuesOut::TYPE &&
        CreditAddinLoad() &&
        true;

    if (!success){
        throw ModelException("CCreditLib::registerClasses",
                             "Registration error");
    }
}

DRLIB_END_NAMESPACE
