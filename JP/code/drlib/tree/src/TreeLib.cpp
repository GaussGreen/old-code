//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TreeLib.cpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : André Segger
//
//   Date        : 22 May 2001
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TreeLib.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Tree1fLN.hpp" // FILTER (directive to makerules/scripts/filterProductsLib.pl)
#include "edginc/Tree1fCEVJCalib.hpp" // FILTER
#include "edginc/Tree1fLV.hpp" // FILTER
#include "edginc/FD1FLV.hpp" // FILTER
#include "edginc/FD1FLNResettable.hpp" // FILTER
#include "edginc/FD1FLNGeneric.hpp" // FILTER
#include "edginc/FD1FE2C.hpp" // FILTER
#include "edginc/FD1FDDE.hpp" // FILTER

#include "edginc/FD2DSVCJ.hpp" // FILTER
#include "edginc/FD2DSV.hpp" // FILTER
#include "edginc/FD2DEQCR.hpp" // FILTER
#include "edginc/FD2DLN.hpp" // FILTER

#include "edginc/FD1DLN.hpp" // FILTER
#include "edginc/FD1DLV.hpp" // FILTER

#include "edginc/FD1DRetLN.hpp" // FILTER
#include "edginc/FD1DRetLV.hpp" // FILTER
#include "edginc/FD1DRetStochGarf.hpp" // FILTER

#include "edginc/SpreadLossTree.hpp" // FILTER
#include "edginc/ToolkitDebug.hpp"

#include "edginc/OUAuxSpreadProcess.hpp" // FILTER

#include "edginc/FD1DMulti.hpp" // FILTER
#include "edginc/FD1DRegimes.hpp" // FILTER




DRLIB_BEGIN_NAMESPACE

extern bool Tree1fLVGDLoad();

extern bool RateTreeLoad();
extern bool CreditTreeLoad();
extern bool Fix3Load();
extern bool Fix3CBLoad();
extern bool Hyb3Load();
extern bool Hyb3CBLoad();
extern bool Hyb3EQLoad();
extern bool Hyb3CBEQLoad();
extern bool Hyb3CBFXLoad();
extern bool Hyb4tempLoad();
extern bool Hyb4CBLoad();
extern bool SMDLoad();
extern bool Fix32QLoad();
extern bool Fix3TDLoad();
extern bool TMXLoad();
extern bool RatesUtilsLoad();
extern bool QuasiMQLoad();
extern bool OUAuxSpreadProcessLoad();

void CTreeLib::linkInClasses(){
    /* The sole purpose of this method is to ensure that the linker 
       includes all symbols out of the market directory. Many symbols
       are automatically linked because they are used by other classes
       which are already included. */

    /* there is no order dependency - this could be automated */

    // Please don't add any new ...::TYPE entries here ---
    // use FFFLoad() where FFF.cpp is the source filename.

    bool success =
        CTree1fLN::TYPE          &&
        CTree1fCEVJ::TYPE        &&
        CTree1fCEVJCalib::TYPE   &&
        CTree1fLV::TYPE          &&
        FD1F::TYPE               &&
        FD1FDDE::TYPE            &&
        FD1FLV::TYPE             &&
        FD1FLN::TYPE             &&
        FD1FLNResettable::TYPE   &&
        FD1FLNGeneric::TYPE      &&
        FD1FE2C::TYPE            &&
        FD2D::TYPE               &&
        FD2DSV::TYPE             &&
        FD2DSVCJ::TYPE           &&
        FD2DEQCR::TYPE           &&
        FD2DLN::TYPE             &&        
        FD1D::TYPE               &&
        FD1DLV::TYPE             &&
        FD1DLN::TYPE             &&

        FD1DRet::TYPE            &&
        FD1DRetLN::TYPE          &&
        FD1DRetLV::TYPE          &&
        FD1DRetStochGarf::TYPE   &&

        Tree1fLVGDLoad()         &&
        RatesUtilsLoad()         &&
        SpreadLossTree::TYPE     &&

        RateTreeLoad()           &&
        Fix3Load()               &&
        Fix3CBLoad()             &&
        Hyb3Load()               &&
        Hyb3CBLoad()             &&
        Hyb3EQLoad()             &&
        Hyb3CBEQLoad()           &&
        Hyb3CBFXLoad()           &&
        Hyb4tempLoad()           &&
        Hyb4CBLoad()             &&
        SMDLoad()                &&
        Fix3TDLoad()             &&
        Fix32QLoad()             &&
        TMXLoad()                &&
        CreditTreeLoad()         &&
        
        QuasiMQLoad()            &&

        OUAuxSpreadProcessLoad() &&

        FD1DMulti::TYPE          &&
        FD1DRegimes::TYPE        &&
        
        true;
           
    if (!success){
        throw ModelException("CTreeLib::registerClasses",
                             "Registration error");
    }

}
DRLIB_END_NAMESPACE
