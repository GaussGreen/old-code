//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRDateVariable.cpp
//
//   Description : Variable of date type
//
//   Date        : May 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRController.hpp"
#include "edginc/Format.hpp"
#include "edginc/FRFactory.hpp"

DRLIB_BEGIN_NAMESPACE

class FRDateVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::dateType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    FRDateVariable(): LValueBase(TYPE) {}

     /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Date type is not implemented in IMS";
        result[1] = name;
        result[2] = "Date type is not implemented in IMS";

        return result;
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        return new FR::LValueDate(name.c_str());
    }

    FRDateVariable(CClassConstSP clazz): LValueBase(clazz) {}
    FRDateVariable(CClassConstSP clazz, const string& name):
        LValueBase(clazz), name(name) {}

    // field
    const string name;

private:
    static IObject* defaultFRDateVariable(){
        return new FRDateVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDateVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRDateVariable);
        FIELD(name, "name of variable");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_DATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI Date Variable",
                                   TYPE);
    }
    
};

CClassConstSP const FRDateVariable::TYPE =
CClass::registerClassLoadMethod("FRDateVariable", typeid(FRDateVariable),
                                load);

class FRDateConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::dateType;
     }
    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "ConstDate type is not implemented in IMS";
        result[1] = id;
        result[2] = "ConstDate type is not implemented in IMS";

        return result;
    }

    /** constant date value ie same value at all timepoints.  */
    class LValueConstDate: public FR::LValueConst,
                           public virtual FRIfaces::ILValueDate{
    private:
        struct MyRT{
            TGetValue*     func;
            DateTime::Date value;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
            explicit MyRT(const DateTime::Date&   value): 
                func(&getValue), value(value){}

            static const DateTime::Date& getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                return (rt->value);
            }
        };
        MyRT* rt;
    public:
        ~LValueConstDate(){
            delete rt;
        }

        // simple constructor
        LValueConstDate(const DateTime::Date& theValue): 
            rt(new MyRT(theValue)){}
            

        virtual void setValue(const DateTime::Date& value) {
            set(IObjectSP()); //set implemented in FR::LValueConst
        }
        /** Just wrapper around getValue method */
        IObjectConstSP get(){
            return IObjectConstSP::attachToRef(&rt->value);
        }
        virtual const DateTime::Date& getValue(){
            return MyRT::getValue(rt);
        }
        //// get the run-time object to use ie cut down version of whole class
        FRIfaces::IRValueDate::RT* getRT(){
            return (FRIfaces::IRValueDate::RT*)rt;
        }
    };

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return id;
    }

    FRDateConstVariable(): LValueBase(TYPE) {}

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRDateConstVariable::getLValue", "Variable "+
                             id+" is read-only");
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        DateTime::Date* value;
        int  size = constValues->size();
        int  numDates = frCtrl->numDates();
        if (size == 1){
            value = constValues->front().get();
            FRIfaces::IRValueSP rValue(new LValueConstDate(*value));
            // store all values now - avoid creating duplicate objects
            for (int i = 0; i < numDates; i++){
                frCtrl->setRValue(i, this, rValue);
            }
            return 0; // indicates that we've saved the value ourselves
        } else if (size != numDates){
            string m("Simulation has "+Format::toString(numDates)+
                     " dates but "+Format::toString(size)+" values "+
                     "supplied");
            throw ModelException("FRDateConstVariable", m);
        } else {
            value = (*constValues)[index].get();
        }
        return new LValueConstDate(*value);
    }

private:
    // fields
    const string                  id;
    const DateTime::DateArraySP   constValues;

    static IObject* defaultFRDateVariable(){
        return new FRDateConstVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDateConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(id, "name of variable");
        FIELD(constValues, "Constant value(s) for variable");
        EMPTY_SHELL_METHOD(defaultFRDateVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_DATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI Const Date Variable",
                                   TYPE);
    }
    
};

CClassConstSP const FRDateConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRDateConstVariable", typeid(FRDateConstVariable), load);

class FRValueDateVariable: public FRDateVariable{
public:
    static CClassConstSP const TYPE;

    FRValueDateVariable(): FRDateVariable(TYPE) {}

    FRValueDateVariable(const string& name): FRDateVariable(TYPE, name) {}

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRValueDateVariable::getLValue", "Variable "+
                             name+" is read-only");
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "ValueDate";
        result[1] = name;
        result[2] = "-";

        return result;
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        // get hold of sim date
        const DateTime& valDate = frCtrl->productView()->getValueDate();
        return new FRDateConstVariable::LValueConstDate(
            DateTime::Date(valDate.getDate()));
    }

private:
    static IObject* defaultFRValueDateVariable(){
        return new FRValueDateVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRValueDateVariable::create",
                                 "Value should be empty for value "
                                 "date variable");
        }
        return new FRValueDateVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRValueDateVariable, clazz);
        SUPERCLASS(FRDateVariable);
        EMPTY_SHELL_METHOD(defaultFRValueDateVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_VALUE_DATE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI "
                                   "Value Date Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "ValueDate");
    }
    
};

CClassConstSP const FRValueDateVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRValueDateVariable", typeid(FRValueDateVariable), load);

bool loadFRDateVariable() {
    return (FRDateVariable::TYPE != 0);
}

DRLIB_END_NAMESPACE

