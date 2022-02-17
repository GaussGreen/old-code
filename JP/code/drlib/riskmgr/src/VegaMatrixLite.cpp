//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : VegaMatrixLite.hpp
//
//   Description : Class used to drive computation of 'Lite' Vega Matrix 
//                 sensitivity
//
//   Author      : Jon Dee
//
//   Date        : 6 October 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/VegaMatrixLite.hpp"
#include "edginc/VegaMatrix.hpp"
#include "edginc/Results.hpp"
#include "edginc/SensMgr.hpp"
#include "edginc/TweakGroup.hpp"
#include "edginc/SensitiveStrikes.hpp"

DRLIB_BEGIN_NAMESPACE


// Integrator1D
void ISupportVegaMatrixLite::load(CClassSP &clazz) {
    REGISTER_INTERFACE(ISupportVegaMatrixLite, clazz);
    EXTENDS(IObject);
    clazz->setPublic();
}


CClassConstSP const ISupportVegaMatrixLite::TYPE = 
CClass::registerInterfaceLoadMethod("ISupportVegaMatrixLite",
                                    typeid(ISupportVegaMatrixLite), load);

////////////////////////////////////////////////////////////////////////////////////////////


const string VegaMatrixLite::NAME = "VEGA_MATRIX_LITE";

VegaMatrixLite::VegaMatrixLite() : Sensitivity(TYPE) {}

void VegaMatrixLite::load(CClassSP& clazz) {
    clazz->setPublic();
    REGISTER(VegaMatrixLite, clazz);
    SUPERCLASS(Sensitivity);
    EMPTY_SHELL_METHOD(&DefaultConstructor<VegaMatrixLite>::iObject);
    FIELD(shiftSize, "How big to make the tweak");
}

CClassConstSP const VegaMatrixLite::TYPE = CClass::registerClassLoadMethod(
    "VegaMatrixLite", typeid(VegaMatrixLite), load);

/** identifies the name used for storing associated results in the output*/
const string& VegaMatrixLite::getSensOutputName() const{
    return NAME;
}

bool VegaMatrixLite::discreteShift() const
{
    return false;
}

void VegaMatrixLite::calculate(TweakGroup* tweakGroup, Results* results) {
    static const string method = "VegaMatrixLite::calculate";

    try {
        Instrument* inst = tweakGroup->getInstrument();
        IModel* model = tweakGroup->getModel();

        // Check for support of interfaces on the instrument
        ISensitiveStrikes* pVM = dynamic_cast<ISensitiveStrikes*>(inst);
        ISupportVegaMatrixLite*  pVML = dynamic_cast<ISupportVegaMatrixLite *>(inst);
        
        if(!pVM){
            throw ModelException(method, "Instrument does not support the ISensitiveStrikes interface");
        } else if (!pVML) {
            throw ModelException(method, "Instrument does not support the ISupportVegaMatrixLite interface");
        } else {
            // Give the instrument a chance to throw an error if at run time VEGA_MATRIX_LITE
            // is not supported
            pVML->avoidVegaMatrixLite(model);
            // Invoke call to price method. Note the control is now saying we are doing a VEGA_MATRIX_LITE
            ResultsSP resultsTemp(new Results());
            model->Price(inst, control, resultsTemp.get());

            results->merge(getSensOutputName(), resultsTemp.get());
        }
    } catch(exception& e) {
        results->storeGreek(IObjectSP(new Untweakable(e)), 
                            getSensOutputName(),
                            OutputNameSP(new OutputName("")));
    }
}

double VegaMatrixLite::getShiftSize() const
{
    return shiftSize;
}


DRLIB_END_NAMESPACE
