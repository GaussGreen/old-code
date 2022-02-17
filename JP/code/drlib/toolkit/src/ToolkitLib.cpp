
#include "edginc/config.hpp"
#include "edginc/ToolkitLib.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Interval.hpp"
#include "edginc/EnergyTermPeriodDates.hpp"
#include "edginc/EnergyTermPeriodLabels.hpp"
#include "edginc/EnergyTermPeriodBenchmark.hpp"
#include "edginc/EnergyDWMYPeriod.hpp"
#include "edginc/ToolkitDebug.hpp"

DRLIB_BEGIN_NAMESPACE
extern bool RevertTypeConvertLoad();
extern bool ErrorHandlerAddinLoad();
extern bool EnumLoad();

void CToolkitLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker 
       includes all symbols out of the toolkit directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included.

       An example of symbols that could be dropped would be an entire class
       representing a product which was referenced by no other classes */

    /* Since the classes in the toolkit are naturally referenced by many other 
       classes, the list of symbols to include is very short */

    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        CashFlow::TYPE          &&
        ClientRunnable::TYPE    &&
        Interval::TYPE          &&
        IReadableMap::TYPE      &&
        IWriteableMap::TYPE     &&
        IMap::TYPE              &&
        CommandLineParams::TYPE &&
	EnergyTermPeriodDates::TYPE &&
	EnergyTermPeriodLabels::TYPE &&
	EnergyTermPeriodBenchmark::TYPE &&
        EnergyDWMYPeriod::TYPE &&
        RevertTypeConvertLoad() &&
        ErrorHandlerAddinLoad() &&
        EnumLoad() &&
        true;
           
    if (!success){
        throw ModelException("CToolkit::registerClasses",
                             "Registration error");
    }

}
DRLIB_END_NAMESPACE
