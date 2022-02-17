//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRSimpleVariable.cpp
//
//   Description : Variables of Simple Types for Flex Rules
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Format.hpp"
#include "edginc/FRController.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/TabulatedFunc.hpp"
#include "edginc/FRFactory.hpp"

DRLIB_BEGIN_NAMESPACE

class FRDoubleVariable: public FR::LValueBase {
public:
    static CClassConstSP const TYPE;

    ~FRDoubleVariable(){}
    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Double";
        result[1] = name;
        result[2] = "-";

        return result;
    }

    FRDoubleVariable(): LValueBase(TYPE) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        return new FR::LValueDouble(name.c_str());
    }

private:
    // field
    const string name;

    FRDoubleVariable(const string& name): LValueBase(TYPE), name(name) {}

    static IObject* defaultFRDoubleVariable(){
        return new FRDoubleVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRDoubleVariable::create",
                                 "Value should be empty for non-const double!");
        }
        return new FRDoubleVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDoubleVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRDoubleVariable);
        FIELD(name, "name of variable");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_DOUBLE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Double Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRDoubleVariable::create, "Double");
    }
    
};

CClassConstSP const FRDoubleVariable::TYPE = CClass::registerClassLoadMethod(
    "FRDoubleVariable", typeid(FRDoubleVariable), FRDoubleVariable::load);
    
class FRDoubleConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRDoubleConstVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result;
        
        result.push_back("ConstDouble");
        result.push_back(name);
        if (constValues->size()==0) {
            result.push_back("ConstDouble Variable not initialized");
        } else if (constValues->size()==1) {
            result.push_back(Format::toString((*constValues)[0]));
        } else {
            result.push_back("PerSimDate");
            for(int i=0; i<constValues->size(); i++) {
                result.push_back(Format::toString((*constValues)[i]));
            }
        }
        return result;
    }

    /** constant double value ie same value at all timepoints. */
    class LValueConstDouble: public FR::LValueConst,
                             public virtual FRIfaces::ILValueDouble{
    private:
        struct MyRT{
            TGetValue*  func;
            TSetValue*  setFunc;
            double      value;

            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }

            explicit MyRT(double   value): func(&getValue), 
                setFunc(&setValue), value(value){}

            static double getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                return (rt->value);
            }
            static void setValue(void* structure, double value){
                throw ModelException("FRDoubleConstVariable", "Cannot assign a"
                                     " value to a const variable");
            }
        };
        MyRT* rt;

    public:
        ~LValueConstDouble(){
            delete rt;
        }

        // simple constructor
        LValueConstDouble(double theValue): rt(new MyRT(theValue)){}

        virtual void setValue(double value) {
            MyRT::setValue(rt, value);
        }

        virtual IObjectConstSP get(){
            return IObjectConstSP(CDouble::create(getValue()));
        }
        virtual double getValue(){
            return MyRT::getValue(rt);
        }
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDouble::RT* getRT(){
            return (FRIfaces::IRValueDouble::RT*)rt;
        }
    };

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    FRDoubleConstVariable(): LValueBase(TYPE) {}
protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        double value;
        int    size = constValues->size();
        int    numDates = frCtrl->numDates();
        if (size == 1){
            value = constValues->front();
            FRIfaces::IRValueSP rValue(new LValueConstDouble(value));
            // store all values now - avoid creating duplicate objects
            for (int i = 0; i < numDates; i++){
                frCtrl->setRValue(i, this, rValue);
            }
            return 0; // indicates that we've saved the value ourselves
        } else if (size != numDates){
            string m("Simulation has "+Format::toString(numDates)+
                     " dates but "+Format::toString(size)+" values "+
                     "supplied");
            throw ModelException("FRDoubleConstVariable", m);
        } else {
            value = (*constValues)[index];
        }
        return new LValueConstDouble(value);
    }

private:
    // fields
    const string               name;
    const DoubleArrayConstSP   constValues;

    FRDoubleConstVariable(const string&              name,
                          const DoubleArrayConstSP   constValues): 
        LValueBase(TYPE), name(name), constValues(constValues) {}

    static IObject* defaultFRDoubleVariable(){
        return new FRDoubleConstVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        DoubleArraySP dVals(new DoubleArray(0));
        int iInit = 0;
        if (value.size()>1) {
            // i.e. data from lookup so skip the first which is the lookup key
            iInit = 1;
        }
        for(int i=iInit; i<value.size(); i++) { 
            // only supporting the simplest case here
            char*   endPos;
            double  dVal = strtod(value[i].c_str(), &endPos);
            if (*endPos != '\0') {
                // extraneous chars - error
                throw ModelException("FRDoubleConstVariable::create", 
                                     "Badly formed number (" + value[i] +
                                     ") for value of variable named : " + 
                                     varName + " [" + Format::toString(i) + "]");
            }
            dVals->push_back(dVal);
        }
        return new FRDoubleConstVariable(varName, dVals);
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDoubleConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant value(s) for variable");
        EMPTY_SHELL_METHOD(defaultFRDoubleVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_DOUBLE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Double Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRDoubleConstVariable::create, 
                                           "ConstDouble");
    }
    
};

CClassConstSP const FRDoubleConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRDoubleConstVariable", typeid(FRDoubleConstVariable), load);

/** DoubleArray variables where components are either
    double variables themselves or just live within the array */
class FRDoubleArrayVariable: public FR::LValueBase,
                             virtual public FRIfaces::IVarArrayBarrierLevelAssist{
public:
    static CClassConstSP const TYPE;

    ~FRDoubleArrayVariable(){}
    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleArrayType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
      
        if (!variables) {
            StringArray result(3);
        
            result[0] = "DoubleArray";
            result[1] = name;
            result[2] = "-";

            return result;
        }
        else {
            string variablesIMS = variables->size()>0 ? (*variables)[0] : "-";

            for(int i=1;i<variables->size();i++){
                variablesIMS = variablesIMS+","+(*variables)[i];
            }
            StringArray result(3);
        
            result[0] = "DoubleArray";
            result[1] = name;
            result[2] = variablesIMS;

            return result;
        }
    }

    FRDoubleArrayVariable(): LValueBase(TYPE) {}

    // For IVarArrayBarrierLevelAssist
    vector<FRIfaces::IVarBarrierLevelAssist*> getComponentAssists(FRController* frCtrl) const {
        static const string method("FRDoubleArrayVariable::getComponentAssists");
        vector<FRIfaces::IVarBarrierLevelAssist*> result(0);
        if (!!variables){
            for (int i = 0; i < (*variables).size(); i++){
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression((*variables)[i]);
                if (!lValueExp){
                    // really expect this to be tripped earlier in the code below
                    string m("No variable with name "+(*variables)[i]+"\n"
                             "When defining '"+name+"' you must supply the "
                             "names of the variables it represents. ");
                    throw ModelException(method, m);
                }
                result.push_back(dynamic_cast<FRIfaces::IVarBarrierLevelAssist*>(lValueExp));
            }
        }
        return result;
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method("FRDoubleArrayVariable::createRValue");
        if (!variables || variables->empty()){
            return new FR::LValueDoubleArray(name.c_str());
        } else {
            vector<FRIfaces::ILValueDouble*> lValues(variables->size());
            for (unsigned int i = 0; i < lValues.size(); i++){
                // Need to extract FRIfaces::ILValue from controller
                // So start with ILValueExpression
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression((*variables)[i]);
                if (!lValueExp){
                    string m("No variable with name "+(*variables)[i]+"\n"
                             "When defining '"+name+"' you must supply the "
                             "names of the variables it represents. ");
                    throw ModelException(method, m);
                }
                FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                                 frCtrl);
                // cast to double
                FRIfaces::ILValueDouble* lValueDb = 
                    dynamic_cast<FRIfaces::ILValueDouble*>(lValue);
                if (!lValueDb){
                    throw ModelException(method, "Variable with name "+
                                         (*variables)[i]+" must be a "
                                         "double variable");
                }
                // put into our vector
                lValues[i] = lValueDb;
            }
            // finally, build relevant object
            return new FR::LValueDoubleVarArray(lValues);
        }
    }

private:
    // field
    const string name;
    const StringArrayConstSP variables;

    FRDoubleArrayVariable(const string&            name,
                          const StringArrayConstSP variables): 
        LValueBase(TYPE), name(name), variables(variables) {}

    static IObject* defaultFRDoubleArrayVariable(){
        return new FRDoubleArrayVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRDoubleArrayVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // ither empty "value" or special "-" means no value
        StringArraySP variables;
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           StringArray::TYPE, value[0]));
            variables = StringArraySP(StringArraySP::dynamicCast(vars));
        }
        return new FRDoubleArrayVariable(varName, variables);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDoubleArrayVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRDoubleArrayVariable);
        FIELD(name, "name of variable");
        FIELD(variables, "names of component variables");
        FIELD_MAKE_OPTIONAL(variables);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_DOUBLE_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Double Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "DoubleArray");
    }
    
};

CClassConstSP const FRDoubleArrayVariable::TYPE = 
CClass::registerClassLoadMethod(
    "FRDoubleArrayVariable", typeid(FRDoubleArrayVariable), 
    FRDoubleArrayVariable::load);

class FRDoubleArrayConstVariable: public FR::LValueBase {
public:
    static CClassConstSP const TYPE;

    ~FRDoubleArrayConstVariable(){}

    /** is the variable simulation date independent ie same value at
        each simulation date.  */
    virtual bool isSDI() const{
        return true; // we are constant across sim dates
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::doubleArrayType;
    }
    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        string constValuesIMS;
        if (constValues->size()==0) {
            constValuesIMS = "ConstDoubleArray variable not initialized";
        }
        else{
            constValuesIMS = Format::toString((*constValues)[0]);
            for(int i=1;i<constValues->size();i++){
                constValuesIMS = constValuesIMS+","+Format::toString((*constValues)[i]);
            }
        }
        
        StringArray result(3);
        
        result[0] = "ConstDoubleArray";
        result[1] = name;
        result[2] = constValuesIMS;
        
        return result;
    }

    FRDoubleArrayConstVariable(): LValueBase(TYPE) {}
protected:
    class LValueConstDoubleArray:public FR::LValueConst,
                                 public virtual FRIfaces::ILValueDoubleArray{

    private:
        struct MyRT: public FRIfaces::ILValueDoubleArray::RT{
            vector<double>    values;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->values.size();
            }
            static double getValue(void* structure, int index){
                MyRT* rt = (MyRT*) structure;
                return rt->values[index];
            }
            static void setValue(void*                             structure,
                                 int                               index,
                                 FRIfaces::IRValueDoubleArray::RT* rValue){
                throw ModelException("LValueConstDoubleArray",
                                     "Cannot assign a"
                                     " value to a const variable");
            }
            explicit MyRT(vector<double>::const_iterator begin,
                          vector<double>::const_iterator end): 
                values(begin, end){
                size = &getSize;
                func = &getValue;
                setFunc = &setValue;
            }
        };
        MyRT* rt;
    public:
        ~LValueConstDoubleArray(){
            delete rt;
        }

        // simple constructor
        LValueConstDoubleArray(const DoubleArray& values):
            rt(new MyRT(values.begin(), values.end())){}

        virtual IObjectConstSP get(){
            return IObjectConstSP(new DoubleArray(rt->values.begin(), 
                                                  rt->values.end()));
        }
        ////get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueDoubleArray::RT* getRT(){
            return (FRIfaces::IRValueDoubleArray::RT*) rt;
        }
    };
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        FRIfaces::IRValueSP rValue(new LValueConstDoubleArray(*constValues));
        // store all values now - avoid creating duplicate objects
        int    numDates = frCtrl->numDates();
        for (int i = 0; i < numDates; i++){
            frCtrl->setRValue(i, this, rValue);
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string               name;
    const DoubleArrayConstSP   constValues;

    static IObject* defaultConstructor(){
        return new FRDoubleArrayConstVariable();
    }

    // 'factory' approach to building this type of variable
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRDoubleArrayConstVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // either empty "value" or special "-" means no value
        DoubleArraySP values(new DoubleArray(0));
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           DoubleArray::TYPE, value[0]));
            values = DoubleArraySP(DoubleArraySP::dynamicCast(vars));
        }
        FRDoubleArrayConstVariable* lVal = new FRDoubleArrayConstVariable();
        const_cast<string&>(lVal->name) = varName; // avoid writing constructor
        const_cast<DoubleArrayConstSP&>(lVal->constValues) = values;
        return lVal;
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRDoubleArrayConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant values for variable");
        EMPTY_SHELL_METHOD(defaultConstructor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_DOUBLE_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Double "
                                   "Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "ConstDoubleArray");
    }
};
CClassConstSP const FRDoubleArrayConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRDoubleArrayConstVariable", typeid(FRDoubleArrayConstVariable), load);

class FRIntVariable: public FR::LValueBase {
public:
    static CClassConstSP const TYPE;

    ~FRIntVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::intType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Int";
        result[1] = name;
        result[2] = "-";

        return result;
    }

    FRIntVariable(): LValueBase(TYPE) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        return new FR::LValueInt(name.c_str());
    }

private:
    // field
    const string name;

    FRIntVariable(const string& name): LValueBase(TYPE), name(name) {}

    static IObject* defaultFRIntVariable(){
        return new FRIntVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRIntVariable::create",
                                 "Value should be empty for non-const int!");
        }
        return new FRIntVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRIntVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRIntVariable);
        FIELD(name, "name of variable");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_INT",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Int Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRIntVariable::create, "Int");
    }
    
};

CClassConstSP const FRIntVariable::TYPE = CClass::registerClassLoadMethod(
    "FRIntVariable", typeid(FRIntVariable), FRIntVariable::load);
    
class FRIntConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRIntConstVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::intType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result;
        
        result.push_back("ConstInt");
        result.push_back(name);
        if (constValues->size()==0) {
            result.push_back("ConstInt Variable not initialized");
        } else if (constValues->size()==1) {
            result.push_back(Format::toString((*constValues)[0]));
        } else {
            result.push_back("PerSimDate");
            for(int i=0; i<constValues->size(); i++) {
                result.push_back(Format::toString((*constValues)[i]));
            }
        }
        return result;
    }

    /** constant int value ie same value at all timepoints. */
    class LValueConstInt: public FR::LValueConst,
                          public virtual FRIfaces::ILValueInt{
    private:
        struct MyRT{
            TGetValue*  func;
            TSetValue*  setFunc;
            int         value;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }

            explicit MyRT(int  value): func(&getValue), 
                setFunc(&setValue), value(value){}

            static void setValue(void* structure, int value){
                throw ModelException("FRIntConstVariable", "Cannot assign a"
                                     " value to a const variable");
            }
            static int getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                return (rt->value);
            }
        };
        MyRT* rt;

    public:
        ~LValueConstInt(){
            delete rt;
        }
        // simple constructor
        LValueConstInt(int theValue): rt(new MyRT(theValue)){}

        virtual void setValue(int value) {
            MyRT::setValue(rt, value);
        }

        virtual IObjectConstSP get(){
            return IObjectConstSP(CInt::create(getValue()));
        }
        virtual int getValue(){
            return MyRT::getValue(rt);
        }
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueInt::RT* getRT(){
            return (FRIfaces::IRValueInt::RT*)rt;
        }
    };

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    FRIntConstVariable(): LValueBase(TYPE) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        int value;
        int    size = constValues->size();
        int    numDates = frCtrl->numDates();
        if (size == 1){
            value = constValues->front();
            FRIfaces::IRValueSP rValue(new LValueConstInt(value));
            // store all values now - avoid creating duplicate objects
            for (int i = 0; i < numDates; i++){
                frCtrl->setRValue(i, this, rValue);
            }
            return 0; // indicates that we've saved the value ourselves

        } else if (size != numDates){
            string m("Simulation has "+Format::toString(numDates)+
                     " dates but "+Format::toString(size)+" values "+
                     "supplied");
            throw ModelException("FRIntConstVariable", m);
        } else {
            value = (*constValues)[index];
        }
        return new LValueConstInt(value);
    }

private:
    // fields
    const string                   name;
    const smartConstPtr<IntArray>  constValues;

    FRIntConstVariable(const string&              name,
                       const IntArrayConstSP   constValues): 
        LValueBase(TYPE), name(name), constValues(constValues) {}

    static IObject* defaultFRIntVariable(){
        return new FRIntConstVariable();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        IntArraySP lVals(new IntArray(0));
        int iInit = 0;
        if (value.size()>1) {
            // i.e. data from lookup so skip the first which is the lookup key
            iInit = 1;
        }
        for(int i=iInit; i<value.size(); i++) { 
            // only supporting the simplest case here
            char*   endPos;
            long  lVal = strtol(value[i].c_str(), &endPos, 10 /* base 10 */);
            if (*endPos != '\0') {
                // extraneous chars - error
                throw ModelException("FRIntConstVariable::create", 
                                     "Badly formed number (" + value[i] +
                                     ") for value of variable named : " + 
                                     varName + " [" + Format::toString(i) + "]");
            }
            lVals->push_back(lVal);
        }
        return new FRIntConstVariable(varName, lVals);
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRIntConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant value(s) for variable");
        EMPTY_SHELL_METHOD(defaultFRIntVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_INT",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Int Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRIntConstVariable::create, "ConstInt");
    }
    
};

CClassConstSP const FRIntConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRIntConstVariable", typeid(FRIntConstVariable), load);


/** IntArray variables where components are either
    int variables themselves or just live within the array */
class FRIntArrayVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRIntArrayVariable(){}
    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::intArrayType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        if (!variables) {
            StringArray result(3);
        
            result[0] = "IntArray";
            result[1] = name;
            result[2] = "-";

            return result;
        }
        else {
            string variablesIMS = variables->size()>0 ? (*variables)[0] : "-";
        
            for(int i=1;i<variables->size();i++){
                variablesIMS = variablesIMS+","+(*variables)[i];
            }
            StringArray result(3);
        
            result[0] = "IntArray";
            result[1] = name;
            result[2] = variablesIMS;

            return result;
        }
    }

    FRIntArrayVariable(): LValueBase(TYPE) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method("FRIntArrayVariable::createRValue");
        if (!variables || variables->empty()){
            return new FR::LValueIntArray(name.c_str());
        } else {
            vector<FRIfaces::ILValueInt*> lValues(variables->size());
            for (unsigned int i = 0; i < lValues.size(); i++){
                // Need to extract FRIfaces::ILValue from controller
                // So start with ILValueExpression
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression((*variables)[i]);
                if (!lValueExp){
                    throw ModelException(method,
                                         "No variable with name "+
                                         (*variables)[i]);
                }
                FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                                 frCtrl);
                // cast to int
                FRIfaces::ILValueInt* lValueDb = 
                    dynamic_cast<FRIfaces::ILValueInt*>(lValue);
                if (!lValueDb){
                    throw ModelException(method, "Variable with name "+
                                         (*variables)[i]+" must be a "
                                         "int variable");
                }
                // put into our vector
                lValues[i] = lValueDb;
            }
            // finally, build relevant object
            return new FR::LValueIntVarArray(lValues);
        }
    }

private:
    // field
    const string name;
    const StringArrayConstSP variables;

    FRIntArrayVariable(const string&            name,
                       const StringArrayConstSP variables): 
        LValueBase(TYPE), name(name), variables(variables) {}

    static IObject* defaultFRIntArrayVariable(){
        return new FRIntArrayVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRIntArrayVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // either empty "value" or special "-" means no value
        StringArraySP variables;
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           StringArray::TYPE, value[0]));
            variables = StringArraySP(StringArraySP::dynamicCast(vars));
        }
        return new FRIntArrayVariable(varName, variables);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRIntArrayVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRIntArrayVariable);
        FIELD(name, "name of variable");
        FIELD(variables, "names of component variables");
        FIELD_MAKE_OPTIONAL(variables);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_INT_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Int Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRIntArrayVariable::create, "IntArray");
    }
    
};

CClassConstSP const FRIntArrayVariable::TYPE = 
CClass::registerClassLoadMethod(
    "FRIntArrayVariable", typeid(FRIntArrayVariable), 
    FRIntArrayVariable::load);

class FRIntArrayConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRIntArrayConstVariable(){}

    /** is the variable simulation date independent ie same value at
        each simulation date.  */
    virtual bool isSDI() const{
        return true; // we are constant across sim dates
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::intArrayType;
    }
    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        string constValuesIMS;
        if (constValues->size()==0) {
            constValuesIMS = "ConstIntArray variable not initialized";
        }
        else{
            constValuesIMS = Format::toString((*constValues)[0]);
            for(int i=1;i<constValues->size();i++){
                constValuesIMS = constValuesIMS+","+Format::toString((*constValues)[i]);
            }
        }

        StringArray result(3);
        
            result[0] = "ConstIntArray";
            result[1] = name;
            result[2] = constValuesIMS;

            return result;
    }

    FRIntArrayConstVariable(): LValueBase(TYPE) {}
protected:
    class LValueConstIntArray:public FR::LValueConst,
                              public virtual FRIfaces::ILValueIntArray{
    private:
        struct MyRT: public FRIfaces::ILValueIntArray::RT{
            vector<int>    values;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->values.size();
            }
            static int getValue(void* structure, int index){
                MyRT* rt = (MyRT*) structure;
                return rt->values[index];
            }
            static void setValue(void*                          structure,
                                 int                            index,
                                 FRIfaces::IRValueIntArray::RT* rValue){
                throw ModelException("LValueConstIntArray",
                                     "Cannot assign a"
                                     " value to a const variable");
            }
            explicit MyRT(vector<int>::const_iterator begin,
                          vector<int>::const_iterator end): 
                values(begin, end){
                size = &getSize;
                func = &getValue;
                setFunc = &setValue;
            }
        };
        MyRT* rt;
    public:
        ~LValueConstIntArray(){
            delete rt;
        }

        // simple constructor
        LValueConstIntArray(const IntArray& values):
            rt(new MyRT(values.begin(), values.end())){}

        virtual IObjectConstSP get(){
            return IObjectConstSP(new IntArray(rt->values.begin(), 
                                                  rt->values.end()));
        }
        ////get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueIntArray::RT* getRT(){
            return (FRIfaces::IRValueIntArray::RT*) rt;
        }
    };
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        FRIfaces::IRValueSP rValue(new LValueConstIntArray(*constValues));
        // store all values now - avoid creating duplicate objects
        int    numDates = frCtrl->numDates();
        for (int i = 0; i < numDates; i++){
            frCtrl->setRValue(i, this, rValue);
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string            name;
    const IntArrayConstSP   constValues;

    // 'factory' approach to building this type of variable
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRIntArrayConstVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // either empty "value" or special "-" means no value
        IntArraySP values(new IntArray(0));
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           IntArray::TYPE, value[0]));
            values = IntArraySP(IntArraySP::dynamicCast(vars));
        }
        FRIntArrayConstVariable* lVal = new FRIntArrayConstVariable();
        const_cast<string&>(lVal->name) = varName; // avoid writing constructor
        const_cast<IntArrayConstSP&>(lVal->constValues) = values;
        return lVal;
    }

    static IObject* defaultConstructor(){
        return new FRIntArrayConstVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRIntArrayConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant values for variable");
        EMPTY_SHELL_METHOD(defaultConstructor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_INT_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Int "
                                   "Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "ConstIntArray");
    }
};
CClassConstSP const FRIntArrayConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRIntArrayConstVariable", typeid(FRIntArrayConstVariable), load);

class FRBoolVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRBoolVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Bool";
        result[1] = name;
        result[2] = "-";
        
        return result;
    }

    FRBoolVariable(): LValueBase(TYPE){}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        return new FR::LValueBool(name.c_str());
    }

private:
    // field
    const string name;

    FRBoolVariable(const string& name): LValueBase(TYPE), name(name) {}

    static IObject* defaultFRBoolVariable(){
        return new FRBoolVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRBoolVariable::create",
                                 "Value should be empty for non-const bool!");
        }
        return new FRBoolVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBoolVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRBoolVariable);
        FIELD(name, "name of variable");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_BOOL",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Bool Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRBoolVariable::create, "Bool");
    }
};

CClassConstSP const FRBoolVariable::TYPE =
CClass::registerClassLoadMethod("FRBoolVariable", typeid(FRBoolVariable), 
                                load);

class FRBoolConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRBoolConstVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result;
        
        result.push_back("ConstBool");
        result.push_back(name);
        if (constValues->size()==0) {
            result.push_back("ConstBool Variable not initialized");
        } else if (constValues->size()==1) {
            result.push_back(Format::toString((*constValues)[0]));
        } else {
            result.push_back("PerSimDate");
            for(int i=0; i<constValues->size(); i++) {
                result.push_back(Format::toString((*constValues)[i]));
            }
        }
        return result;
    }

    /** constant bool value ie same value at all timepobools. */
    class LValueConstBool: public FR::LValueConst,
                           public virtual FRIfaces::ILValueBool{
    private:
        struct MyRT{
            TGetValue*  func;
            TSetValue*  setFunc;
            bool        value;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }

            explicit MyRT(bool   value): func(&getValue), 
                setFunc(&setValue), value(value){}

            static bool getValue(void* structure){
                MyRT* rt = (MyRT*)structure;
                return (rt->value);
            }
            static void setValue(void* structure, bool value){
                throw ModelException("FRBoolConstVariable", "Cannot assign a"
                                     " value to a const variable");
            }
        };
        MyRT* rt;

    public:
        ~LValueConstBool(){
            delete rt;
        }

        // simple constructor
        LValueConstBool(bool theValue): rt(new MyRT(theValue)){}

        virtual void setValue(bool value) {
            MyRT::setValue(rt, value);
        }

        virtual IObjectConstSP get(){
            return IObjectConstSP(CBool::create(getValue()));
        }
        virtual bool getValue(){
            return MyRT::getValue(rt);
        }
        //// get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBool::RT* getRT(){
            return (FRIfaces::IRValueBool::RT*)rt;
        }
    };

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    FRBoolConstVariable(): LValueBase(TYPE), constValues(0) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        bool value;
        int  size = constValues->size();
        int  numDates = frCtrl->numDates();
        if (size == 1){
            value = constValues->front();
            FRIfaces::IRValueSP rValue(new LValueConstBool(value));
            // store all values now - avoid creating duplicate objects
            for (int i = 0; i < numDates; i++){
                frCtrl->setRValue(i, this, rValue);
            }
            return 0; // indicates that we've saved the value ourselves
        } else if (size != numDates){
            string m("Simulation has "+Format::toString(numDates)+
                     " dates but "+Format::toString(size)+" values "+
                     "supplied");
            throw ModelException("FRBoolConstVariable", m);
        } else {
            value = (*constValues)[index];
        }
        return new LValueConstBool(value);
    }

private:
    // fields
    const string            name; 
    const BoolArrayConstSP  constValues;

    FRBoolConstVariable(const string&              name,
                        const BoolArrayConstSP     constValues): 
        LValueBase(TYPE), name(name), constValues(constValues) {}
    
    static IObject* defaultFRBoolVariable(){
        return new FRBoolConstVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        CBoolArraySP bVals(new BoolArray(0));
        int iInit = 0;
        if (value.size()>1) {
            // i.e. data from lookup so skip the first which is the lookup key
            iInit = 1;
        }
        for(int i=iInit; i<value.size(); i++) {
            // only supporting the simplest case here
            bool    bVal;
            if (value[i] == "true") {
                bVal = true;
            } else if (value[i] == "false") {
                bVal = false;
            } else {
                throw ModelException("FRBoolConstVariable::create", 
                                     "Badly formed (" + value[i] +
                                     ") for value of variable named : " + 
                                     varName + "[" + Format::toString(i) +
                                     "]");
            }
            bVals->push_back(bVal);
        }
        return new FRBoolConstVariable(varName, bVals);
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBoolConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant value(s) for variable");
        EMPTY_SHELL_METHOD(defaultFRBoolVariable);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_BOOL",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Bool Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRBoolConstVariable::create, "ConstBool");
    }
    
};

CClassConstSP const FRBoolConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRBoolConstVariable", typeid(FRBoolConstVariable), load);

/** BoolArray variables where components are either
    bool variables themselves or just live within the array */
class FRBoolArrayVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRBoolArrayVariable(){}
    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolArrayType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        if (!variables) {
            StringArray result(3);
            
            result[0] = "BoolArray";
            result[1] = name;
            result[2] = "-";
            
            return result;
        }
        else {
            string variablesIMS = variables->size()>0 ? (*variables)[0] : "-";
            
            for(int i=1;i<variables->size();i++){
                variablesIMS = variablesIMS+","+(*variables)[i];
            }
            StringArray result(3);
            
            result[0] = "BoolArray";
            result[1] = name;
            result[2] = variablesIMS;
            
            return result;
        }
    }


    FRBoolArrayVariable(): LValueBase(TYPE) {}

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method("FRBoolArrayVariable::createRValue");
        if (!variables || variables->empty()){
            return new FR::LValueBoolArray(name.c_str());
        } else {
            vector<FRIfaces::ILValueBool*> lValues(variables->size());
            for (unsigned int i = 0; i < lValues.size(); i++){
                // Need to extract FRIfaces::ILValue from controller
                // So start with ILValueExpression
                FRIfaces::ILValueExpression* lValueExp = 
                    frCtrl->getLValueExpression((*variables)[i]);
                if (!lValueExp){
                    throw ModelException(method,
                                         "No variable with name "+
                                         (*variables)[i]);
                }
                FRIfaces::ILValue* lValue = lValueExp->getLValue(index,
                                                                 frCtrl);
                // cast to bool
                FRIfaces::ILValueBool* lValueDb = 
                    dynamic_cast<FRIfaces::ILValueBool*>(lValue);
                if (!lValueDb){
                    throw ModelException(method, "Variable with name "+
                                         (*variables)[i]+" must be a "
                                         "bool variable");
                }
                // put into our vector
                lValues[i] = lValueDb;
            }
            // finally, build relevant object
            return new FR::LValueBoolVarArray(lValues);
        }
    }

private:
    // field
    const string name;
    const StringArrayConstSP variables;

    FRBoolArrayVariable(const string&            name,
                        const StringArrayConstSP variables): 
        LValueBase(TYPE), name(name), variables(variables) {}

    static IObject* defaultFRBoolArrayVariable(){
        return new FRBoolArrayVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRBoolArrayVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // ither empty "value" or special "-" means no value
        StringArraySP variables;
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           StringArray::TYPE, value[0]));
            variables = StringArraySP(StringArraySP::dynamicCast(vars));
        }
        return new FRBoolArrayVariable(varName, variables);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBoolArrayVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRBoolArrayVariable);
        FIELD(name, "name of variable");
        FIELD(variables, "names of component variables");
        FIELD_MAKE_OPTIONAL(variables);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_BOOL_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Bool Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRBoolArrayVariable::create, "BoolArray");
    }
    
};

CClassConstSP const FRBoolArrayVariable::TYPE = 
CClass::registerClassLoadMethod(
    "FRBoolArrayVariable", typeid(FRBoolArrayVariable), 
    FRBoolArrayVariable::load);

class FRBoolArrayConstVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRBoolArrayConstVariable(){}

    /** is the variable simulation date independent ie same value at
        each simulation date.  */
    virtual bool isSDI() const{
        return true; // we are constant across sim dates
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::boolArrayType;
    }
    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        string constValuesIMS;
        if (constValues->size()==0) {
            constValuesIMS = "ConstBoolArray variable not initialized";
        }
        else{
            constValuesIMS = Format::toString((*constValues)[0]);
            for(int i=1;i<constValues->size();i++){
                constValuesIMS = constValuesIMS+","+Format::toString((*constValues)[i]);
            }
        }
        
        StringArray result(3);
        
        result[0] = "ConstBoolArray";
        result[1] = name;
        result[2] = constValuesIMS;
        
        return result;
    }

    FRBoolArrayConstVariable(): LValueBase(TYPE) {}
protected:
    class LValueConstBoolArray:public FR::LValueConst,
                               public virtual FRIfaces::ILValueBoolArray{
    private:
        struct MyRT: public FRIfaces::ILValueBoolArray::RT{
            vector<bool>    values;
            //// uses FR::MemMgr
            void* operator new(size_t size){
                return FR::MemMgr::alloc(size);
            }
            void operator delete(void *ptr){
                FR::MemMgr::dealloc(ptr);
            }
            static int getSize(void* structure){
                MyRT* rt = (MyRT*) structure;
                return rt->values.size();
            }
            static bool getValue(void* structure, int index){
                MyRT* rt = (MyRT*) structure;
                return rt->values[index];
            }
            static void setValue(void*                           structure,
                                 int                             index,
                                 FRIfaces::IRValueBoolArray::RT* rValue){
                throw ModelException("LValueConstBoolArray",
                                     "Cannot assign a"
                                     " value to a const variable");
            }
            explicit MyRT(const BoolArray& bools): 
                values(bools.size()){
                size = &getSize;
                func = &getValue;
                setFunc = &setValue;
                for (unsigned int i = 0; i < values.size(); i++){
                    values[i] = bools[i];
                }
            }
        };
        MyRT* rt;
    public:
        ~LValueConstBoolArray(){
            delete rt;
        }

        // simple constructor
        LValueConstBoolArray(const BoolArray& values):
            rt(new MyRT(values)){}

        virtual IObjectConstSP get(){
            CBoolArraySP bools(new BoolArray(rt->values.size()));
            for (unsigned int i = 0; i < rt->values.size(); i++){
                (*bools)[i] = rt->values[i];
            }
            return IObjectSP(bools);
        }
        ////get the run-time object to use ie cut down version of whole class
        virtual FRIfaces::IRValueBoolArray::RT* getRT(){
            return (FRIfaces::IRValueBoolArray::RT*) rt;
        }
    };
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        FRIfaces::IRValueSP rValue(new LValueConstBoolArray(*constValues));
        // store all values now - avoid creating duplicate objects
        int    numDates = frCtrl->numDates();
        for (int i = 0; i < numDates; i++){
            frCtrl->setRValue(i, this, rValue);
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string              name;
    const CBoolArrayConstSP   constValues;

    // 'factory' approach to building this type of variable
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine("FRBoolArrayConstVariable::create");
        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not yet supported here");
        }
        // either empty "value" or special "-" means no value
        CBoolArraySP values(new BoolArray(0));
        if (!value[0].empty() && value[0] != "-" ) {
            IArraySP vars(FRVarFactory::parseForArrayNames(routine, varName, 
                                                           BoolArray::TYPE, value[0]));
            values = CBoolArraySP(CBoolArraySP::dynamicCast(vars));
        }
        FRBoolArrayConstVariable* lVal = new FRBoolArrayConstVariable();
        const_cast<string&>(lVal->name) = varName; // avoid writing constructor
        const_cast<CBoolArrayConstSP&>(lVal->constValues) = values;
        return lVal;
    }

    static IObject* defaultConstructor(){
        return new FRBoolArrayConstVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRBoolArrayConstVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        FIELD(name, "name of variable");
        FIELD(constValues, "Constant values for variable");
        EMPTY_SHELL_METHOD(defaultConstructor);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_CONST_BOOL_ARRAY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Const Bool "
                                   "Array Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(create, "ConstBoolArray");
    }
};
CClassConstSP const FRBoolArrayConstVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRBoolArrayConstVariable", typeid(FRBoolArrayConstVariable), load);

/** A const array of ints which represent the current index. ie length =
    number of sim dates. Values are {0, 1, 2, 3, ...} */
class FRIndexVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRIndexVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::intType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    virtual const string& getID() const{
        return name;
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRValueDateVariable::getLValue", "Variable "+
                             name+" is read-only");
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Index";
        result[1] = name;
        result[2] = "-";

        return result;
    }

protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        return new FRIntConstVariable::LValueConstInt(index);
    }

private:
    // fields
    const string                   name;

    FRIndexVariable(): LValueBase(TYPE) {}

    FRIndexVariable(const string& name): LValueBase(TYPE), name(name) {}

    static IObject* defaultConstructor(){
        return new FRIndexVariable();
    }
    
    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        // require either empty "value" or special "-" meaning no value
        if (!value.empty() && !value[0].empty() && value[0] != "-" ) {
            throw ModelException("FRIndexVariable::create",
                                 "Value should be empty for Index variable!");
        }
        return new FRIndexVariable(varName);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRIndexVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "name of variable");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_INDEX",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Index Variable",
                                   TYPE);
        FRVarFactory::registerCreateMethod(FRIndexVariable::create, "Index");
    }
    
};

CClassConstSP const FRIndexVariable::TYPE =
CClass::registerClassLoadMethod(
    "FRIndexVariable", typeid(FRIndexVariable), load);

/** holds a schedule which can be interpolated at a given date */
class FRScheduleVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRScheduleVariable(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::scheduleType;
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "Schedule type is not implemented in IMS";
        result[1] = name;
        result[2] = "Schedule type is not implemented in IMS";
        
        return result;
    }

    class FRSchedule: public FRIfaces::IRValueSchedule{
    private:
        ScheduleConstSP schedule;
    public:
        ~FRSchedule(){}

        FRSchedule(const ScheduleConstSP& schedule): schedule(schedule){}

        /** Returns the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double) */
        virtual IObjectConstSP get() {
            return (IObjectConstSP)schedule;
        }

        /** Is the value of this object known before the simulation starts
            eg a constant */
        virtual bool isKnown() const{
            return true; // only support constant schedules for now
        }

        // get the variable expressed as a double
        virtual const Schedule* getValue() {
            return schedule.get();
        }
    };
protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        // save memory by using the same schedule at each point
        const FRIfaces::IProductView* productView = frCtrl->productView();
        int numDates = productView->getSimDates()->size();
        FRIfaces::IRValueSP rValue(new FRSchedule(schedule));
        for (int i = 0; i < numDates; i++){
            frCtrl->setRValue(i, this, rValue);
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // field
    const string          name;
    const ScheduleConstSP schedule;

    FRScheduleVariable(): LValueBase(TYPE) {}

    static IObject* defaultFRScheduleVariable(){
        return new FRScheduleVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRScheduleVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRScheduleVariable);
        FIELD(name, "name of variable");
        FIELD(schedule, "schedule of doubles");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_SCHEDULE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Schedule Variable",
                                   TYPE);
    }
    
};
CClassConstSP const FRScheduleVariable::TYPE = CClass::registerClassLoadMethod(
    "FRScheduleVariable", typeid(FRScheduleVariable), load);
    
/** holds a TabulatedFunc which can be interpolated at a given date */
class FRTabulatedFuncVariable: public FR::LValueBase{
public:
    static CClassConstSP const TYPE;

    ~FRTabulatedFuncVariable(){}

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const{
        return FRIfaces::tabulatedFuncType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
            result[0] = "TabulatedFunction type is not implemented in IMS";
            result[1] = name;
            result[2] = "TabulatedFunction type is not implemented in IMS";

            return result;
    }

    class FRTabulatedFunc: public FRIfaces::IRValueTabulatedFunc{
    private:
        TabulatedFuncConstSP tabulatedFunc;
    public:
        ~FRTabulatedFunc(){}

        FRTabulatedFunc(const TabulatedFuncConstSP& tabulatedFunc): 
            tabulatedFunc(tabulatedFunc){}

        /** Returns the value of the variable at this time point. More of
            indicative method as derived types should supply method allowing
            easier extraction of variable (eg as double) */
        virtual IObjectConstSP get() {
            return (IObjectConstSP)tabulatedFunc;
        }

        /** Is the value of this object known before the simulation starts
            eg a constant. */
        virtual bool isKnown() const{
            return true; // only support const tabulated funcs for now
        }

        // get the variable expressed as a double
        virtual const TabulatedFunc* getValue() {
            return tabulatedFunc.get();
        }        
    };
protected:
    /** creates IRValue which represents the value
        of this expression at a specific timepoint */
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        // save memory by using the same TabulatedFunc at each point
        const FRIfaces::IProductView* productView = frCtrl->productView();
        int numDates = productView->getSimDates()->size();
        FRIfaces::IRValueSP rValue(new FRTabulatedFunc(tabulatedFunc));
        for (int i = 0; i < numDates; i++){
            frCtrl->setRValue(i, this, rValue);
        }
        return 0; // indicates that we've saved the value ourselves
    }

private:
    // fields
    const string               name;
    const TabulatedFuncConstSP tabulatedFunc;

    FRTabulatedFuncVariable(): LValueBase(TYPE) {}

    static IObject* defaultFRTabulatedFuncVariable(){
        return new FRTabulatedFuncVariable();
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRTabulatedFuncVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        EMPTY_SHELL_METHOD(defaultFRTabulatedFuncVariable);
        FIELD(name, "name of variable");
        FIELD(tabulatedFunc, "TabulatedFunc of doubles");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_VAR_TABULATED_FUNC",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX TabulatedFunc Variable",
                                   TYPE);
    }
    
};
CClassConstSP const FRTabulatedFuncVariable::TYPE = 
CClass::registerClassLoadMethod(
    "FRTabulatedFuncVariable", typeid(FRTabulatedFuncVariable), load);

bool loadFRSimpleVariables() {
    return (FRDoubleVariable::TYPE != 0);
}

DRLIB_END_NAMESPACE
