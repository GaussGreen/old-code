//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix32Q.cpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Fix32Q.hpp"
#include "edginc/IrConverter.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************** Fix32Q *************************************/

/**  collect model initialisation data, set up timeline  */

void Fix32Q::retrieveFactor() {
    try{
        SingleCcyTree::populateTreeCurves();

        IRVolSelector volSelector(getIRVolRaw(IRParams), treeCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);

        volSelector.getMktVolData(mktVolData);

        // FIX3_CLASSIC is the multi asset wrapper 2.0 style fix3.
        mktVolData.ModelChoice = FIX3_ORIGINAL;
        mktVolData.IsNmrModel = false;
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

//------------------------------------------------------------------------------

void Fix32Q::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(Fix32Q, clazz);
    SUPERCLASS(SingleCcyTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

// interface reflection/class loading mechanism
CClassConstSP const Fix32Q::TYPE = CClass::registerClassLoadMethod(
    "Fix32Q", typeid(Fix32Q), Fix32Q::load);

bool Fix32QLoad(void) {
    return (Fix32Q::TYPE != 0);
}

DRLIB_END_NAMESPACE
