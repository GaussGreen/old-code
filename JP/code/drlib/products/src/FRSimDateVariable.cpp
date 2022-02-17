//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRSimDateVariable.cpp
//
//   Description : Variable providing access to simulation dates
//
//   Author      : Mark A Robson
//
//   Date        : 9 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRController.hpp"
#include "edginc/FRFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Variable providing access to simulation dates */
class FRSimDateVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::dateType;
    }

    /** returns a string that can be used to identify this object */
    virtual const string& getID() const{
        return name;
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRSimDateVariable::getLValue", "Variable "+
                             name+" is read-only");
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        result[0] = "SimDate";
        result[1] = name;
        result[2] = "-";

        return result;
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        
        return new FR::RConstDate(frCtrl->getDate(index).getDate());
    }

private:
    // fields
    const string      name;

    FRSimDateVariable(): LValueBase(TYPE) {}
    FRSimDateVariable(const string& name): LValueBase(TYPE), name(name) {}

    static IObject* defaultFRSimDateVariable(){
        return new FRSimDateVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRSimDateVariable::create",
                                 "Value should be empty for sim date variable");
        }
        return new FRSimDateVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRSimDateVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        EMPTY_SHELL_METHOD(defaultFRSimDateVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_SIMDATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI Sim Date Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "SimDate");
    }
    
};

CClassConstSP const FRSimDateVariable::TYPE =
CClass::registerClassLoadMethod("FRSimDateVariable", typeid(FRSimDateVariable),
                                load);

bool loadFRSimDateVariable() {
    return (FRSimDateVariable::TYPE != 0);
}

DRLIB_END_NAMESPACE
