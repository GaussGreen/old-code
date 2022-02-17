//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : SMD.cpp
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/SMD.hpp"
#include "edginc/IrConverter.hpp"

DRLIB_BEGIN_NAMESPACE

/************************************** SMD *************************************/

/**  collect model initialisation data, set up timeline  */

void SMD::retrieveFactor() {
    try{
        SingleCcyTree::populateTreeCurves();

        IRVolSelector volSelector(getIRVolRaw(IRParams), treeCurves[0],
                                  IRParams->volCalibIndex,
                                  cetSmoothing, IRParams);
        volSelector.getMktVolData(mktVolData);

        mktVolData.ModelChoice = FIX3_SMD;
        mktVolData.IsNmrModel = false;

        // SMD is always 2 factor - 
        // ??? make nbFactors optional in model interface, but have to change to int*
        if (nbFactors != 2)
            throw ModelException("SMD model must always be 2 factor - nbFactor "
                "received = " + Format::toString(nbFactors));
        treeData.NbFactor = 2;

        mktVolData.Afac = correlationSkew;
        mktVolData.Bfac = correlationCurvature;
        mktVolData.Cfac = correlationLevel;
        mktVolData.Dfac = correlationTermStructure;
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

//------------------------------------------------------------------------------

void SMD::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(SMD, clazz);
    SUPERCLASS(SingleCcyTree);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(correlationSkew, "");
    FIELD(correlationCurvature, "");
    FIELD(correlationLevel, "");
    FIELD(correlationTermStructure, "");
}

// interface reflection/class loading mechanism
CClassConstSP const SMD::TYPE = CClass::registerClassLoadMethod(
    "SMD", typeid(SMD), SMD::load);

bool SMDLoad(void) {
    return (SMD::TYPE != 0);
}

DRLIB_END_NAMESPACE
