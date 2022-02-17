//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Untweakable.cpp
//
//   Description : Result object to show that requested result failed
//
//   Author      : Andrew J Swain
//
//   Date        : 2 May 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_UNTWEAKABLE_CPP
#include "edginc/Untweakable.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/SpreadSheetMode.hpp"

DRLIB_BEGIN_NAMESPACE

Untweakable::Untweakable() : CObject(TYPE), message("No error message") {}

Untweakable::Untweakable(const ModelException& e) : CObject(TYPE) 
{
    try {
        char* stack = e.stackTrace();
        message = stack;
        free(stack);

        if (SpreadSheetMode::isOn()) {
            ErrorHandler::writeMsg(message);
        }
    }
    catch (exception& x) {
        throw ModelException(x, "Untweakable::Untweakable");
    }
}

Untweakable::Untweakable(const string& message) : CObject(TYPE), 
    message(message) 
{
    if (SpreadSheetMode::isOn()) {
        ErrorHandler::writeMsg(message);
    }
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void Untweakable::outputWrite(const string& linePrefix,
                              const string& prefix,
                              ostream&      stream) const{
    stream << linePrefix << prefix << ": " << "Untweakable" << endl;
}
/** scale by factor x (implementation of CombinableResult) (does nothing)*/
void Untweakable::scale(double x){
    // empty
}

/** add Untweakable object to this result (Implementation of
    CombinableResult) (does nothing) */
void Untweakable::add(const CombinableResult& x, double scaleFactor){
    // empty - could combine error message but doesn't seem worth the effort
}

/** add an object (of general type) to this result. (implementation of 
    CombinableMixedResult) (does nothing) */
IObject* Untweakable::addResult(const  IObject& x,
                                double scaleFactor) const{
    return new Untweakable(string("Untweakable::addResult: "
                                  "Component was untweakable"));
}

/** the message */
string Untweakable::getMessage() const {
    return message;
}

class UntweakableHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Untweakable, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableMixedResult);
        EMPTY_SHELL_METHOD(defaultUntweakable);
        FIELD(message, "error message");
    }
    
    static IObject* defaultUntweakable(){
        return new Untweakable();
    }
};

CClassConstSP const Untweakable::TYPE = CClass::registerClassLoadMethod(
    "Untweakable", typeid(Untweakable), UntweakableHelper::load);

DRLIB_END_NAMESPACE
