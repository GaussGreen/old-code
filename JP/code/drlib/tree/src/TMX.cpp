//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : TMX.cpp
//
//   Description : temporary version of fix3 to work on the claim/bank and
//                 tree starting today without cluttering up the existing Fix3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TMX.hpp"
#include "edginc/IrConverter.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************** TMX *************************************/

/**  collect model initialisation data, set up timeline  */

void TMX::retrieveFactor()
{
    try
    {
        SingleCcyTree::populateTreeCurves();

        IRVolSelector volSelector(getIRVolRaw(IRParams), treeCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);

        if (nbFactors != 1)
            throw ModelException("TMX currently only supports 1 factor  - supplied " +
                Format::toString(nbFactors) + " factors");

        volSelector.getMktVolData(mktVolData);
        IrConverter::to_MKTVOL_DATA_TMX_SMILE(mktVolData, IRParams);
        
        mktVolData.ModelChoice = FIX3_TMX;
        mktVolData.IsNmrModel = true;

    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

//------------------------------------------------------------------------------

void TMX::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(TMX, clazz);
    SUPERCLASS(SingleCcyTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
}

// interface reflection/class loading mechanism
CClassConstSP const TMX::TYPE = CClass::registerClassLoadMethod(
    "TMX", typeid(TMX), TMX::load);

bool TMXLoad(void) {
    return (TMX::TYPE != 0);
}

DRLIB_END_NAMESPACE
