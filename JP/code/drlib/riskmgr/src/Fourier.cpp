//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Fourier.cpp
//
//   Description : Fourier Process
//
//   Date        : 22 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/Fourier.hpp"

DRLIB_BEGIN_NAMESPACE

// FourierProcess
FourierProcess::FourierProcess(const CClassConstSP& clazz):
    CObject(clazz)
{}

void FourierProcess::load(CClassSP &clazz){
    REGISTER(FourierProcess, clazz);
    SUPERCLASS(CObject);
    clazz->setPublic(); // make visible to EAS/spreadsheet
}

CClassConstSP const FourierProcess::TYPE = CClass::registerClassLoadMethod(
"FourierProcess", typeid(FourierProcess), FourierProcess::load);

// FourierProcess::FrequencyRange
FourierProcess::FrequencyRange::FrequencyRange(): CObject(TYPE){}

/** Returns the parameters of the process */
string FourierProcess::getParameters() const {
    return "FourierProcess parameters not available";
}

// Gets the parameters of the process
string FourierProcess::extractParameters(const Calibrator::IAdjustable* adjustable) {
    static const string routine = "FourierProcess::extractParameters";
        
    try {
        string parameters = "Failed with parameters";

        // get the relevant fields for this object
        CFieldArray fields(Calibrator::IAdjustable::getFields(adjustable->getClass()));
        // turn pointer into SP
        IObjectConstSP obj(IObjectConstSP::attachToRef(adjustable));
        // then loop over them
        for (unsigned int i = 0; i < fields.size(); i++){
            if (fields[i]->getType() == CDouble::TYPE){
                // Double
                string val = Format::toString(fields[i]->getDouble(obj));
                string name = fields[i]->getName();
                parameters += "\n\t" + name + ": " + val;
            } else if (fields[i]->getType() == CDoubleArray::TYPE){
                // DoubleArray
                IObjectConstSP objArray = fields[i]->constGet(obj);
                CDoubleArrayConstSP array = CDoubleArrayConstSP::dynamicCast(objArray);
                if(array.get()) {
                    string val;
                    for(int idx = 0; idx < array->size(); idx++) {
                        val += Format::toString((*array)[idx]) + ", ";
                    }
                    string name = fields[i]->getName();
                    parameters += "\n\t" + name + ": " + val;
                }
            } else {
                // Fail
                throw ModelException("Only Doubles and DoubleArrays are supported");
            }
        }

        return parameters;
    } catch(exception& e) {
        throw ModelException(e, routine);
    }
}

class FourierProcessFrequencyRangeHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FourierProcess::FrequencyRange, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(lowerRealBound, "lowerRealBound");
        FIELD(upperRealBound, "upperRealBound");
    }

    static IObject* defaultCtor(){
        return new FourierProcess::FrequencyRange();
    }
};

CClassConstSP const FourierProcess::FrequencyRange::TYPE = CClass::registerClassLoadMethod(
"FourierProcess::FrequencyRange", typeid(FourierProcess::FrequencyRange), FourierProcessFrequencyRangeHelper::load);

// FwdStFourierProcessLogRtn
void FwdStFourierProcessLogRtn::load(CClassSP &clazz){
    REGISTER_INTERFACE(FwdStFourierProcessLogRtn, clazz);
    EXTENDS(IObject);
}

CClassConstSP const FwdStFourierProcessLogRtn::TYPE = CClass::registerInterfaceLoadMethod(
    "FwdStFourierProcessLogRtn", typeid(FwdStFourierProcessLogRtn), FwdStFourierProcessLogRtn::load);


// StFourierProcessLogRtn
void StFourierProcessLogRtn::load(CClassSP &clazz){
    REGISTER_INTERFACE(StFourierProcessLogRtn, clazz);
    EXTENDS(IObject);
}

CClassConstSP const StFourierProcessLogRtn::TYPE = CClass::registerInterfaceLoadMethod(
    "StFourierProcessLogRtn", typeid(StFourierProcessLogRtn), StFourierProcessLogRtn::load);

// StIntVarFourierProcess
void StFourierProcessIntVar::load(CClassSP &clazz){
    REGISTER_INTERFACE(StFourierProcessIntVar, clazz);
    EXTENDS(IObject);
}

CClassConstSP const StFourierProcessIntVar::TYPE = CClass::registerInterfaceLoadMethod(
    "StFourierProcessIntVar", typeid(StFourierProcessIntVar), StFourierProcessIntVar::load);

// FwdStIntVarFourierProcess
void FwdStFourierProcessIntVar::load(CClassSP &clazz){
    REGISTER_INTERFACE(FwdStFourierProcessIntVar, clazz);
    EXTENDS(IObject);
}

CClassConstSP const FwdStFourierProcessIntVar::TYPE = CClass::registerInterfaceLoadMethod(
    "FwdStFourierProcessIntVar", typeid(FwdStFourierProcessIntVar), FwdStFourierProcessIntVar::load);


void StFourierProcessQuadVar::load(CClassSP &clazz){
    REGISTER_INTERFACE(StFourierProcessQuadVar, clazz);
    EXTENDS(IObject);
}

CClassConstSP const StFourierProcessQuadVar::TYPE = CClass::registerInterfaceLoadMethod(
    "StFourierProcessQuadVar", typeid(StFourierProcessQuadVar), StFourierProcessQuadVar::load);

// FwdStFourierProcessQuadVar
void FwdStFourierProcessQuadVar::load(CClassSP &clazz){
    REGISTER_INTERFACE(FwdStFourierProcessQuadVar, clazz);
    EXTENDS(IObject);
}

CClassConstSP const FwdStFourierProcessQuadVar::TYPE = CClass::registerInterfaceLoadMethod(
    "FwdStFourierProcessQuadVar", typeid(FwdStFourierProcessQuadVar), FwdStFourierProcessQuadVar::load);

// FwdStFourierProcessExpQuadVar
void FwdStFourierProcessExpQuadVar::load(CClassSP &clazz){
    REGISTER_INTERFACE(FwdStFourierProcessExpQuadVar, clazz);
    EXTENDS(IObject);
}

CClassConstSP const FwdStFourierProcessExpQuadVar::TYPE = CClass::registerInterfaceLoadMethod(
    "FwdStFourierProcessExpQuadVar", typeid(FwdStFourierProcessExpQuadVar), FwdStFourierProcessExpQuadVar::load);



DRLIB_END_NAMESPACE
