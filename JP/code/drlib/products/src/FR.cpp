//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FR.cpp
//
//   Description : Defines concrete FlexRules objects available at the 
//                 spreadsheet side of things
//
//   Author      : Mark A Robson
//
//   Date        : 25 March 2002
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/FRFunction.hpp"
#include "edginc/GenericNFBase.hpp"
#include "edginc/FRController.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Format.hpp"
#include "edginc/FRParseException.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Null.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/CliquetVolRequest.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/GetMarket.hpp"
#include "edginc/OutputRequestUtil.hpp"
#include <fstream>

DRLIB_BEGIN_NAMESPACE

class FlexInstrument;
class Flex;

CClassConstSP const FR::RValueBase::TYPE = CClass::registerClassLoadMethod(
    "FR::RValueBase", typeid(RValueBase), load);

/** returns shallow copy */
IObject* FR::RValueBase::clone() const{
    ensurePosRefCount();
    return const_cast<FR::RValueBase*>(this);
}

FR::RValueBase::RValueBase(const CClassConstSP& clazz):
    CObject(clazz) {}

/** Used in conjunction with equals() by FlexRule to keep a
    hash of RValueExpressions. Here, objects with the pointer value
    are deemed to refer to the same IRValue. */
int FR::RValueBase::hashCode() const{
    return (size_t) this;
}
        
/** Returns true if the this and the given IRValueExpression refer
    to the same IRValue. Done via comparison of pointers */
bool FR::RValueBase::equals(const FRIfaces::IRValueExpression* name) const{
    return (name == this);
}

/** gets FRIfaces::IRValue from FRController which represents the value
    of this expression at a specific timepoint */
FRIfaces::IRValue* FR::RValueBase::getRValue(
    FRController*  frCtrl) const{
    return getRValue(frCtrl->getIndex(), frCtrl);
}

/** gets FRIfaces::IRValue from FRController which represents the value
    of this expression at the current timepoint */
FRIfaces::IRValue* FR::RValueBase::getRValue(
    int             index,
    FRController*   frCtrl) const{
    FRIfaces::IRValue* rValue = frCtrl->getRValue(this, index);
    if (!rValue){
        rValue = createRValue(index, frCtrl); // create it
        if (rValue){
            // non null indicates value needs to be saved
            frCtrl->setRValue(index, this, 
                              FRIfaces::IRValueSP(rValue)); // save it
        } else {
            // try and get it again
            rValue = frCtrl->getRValue(this, index);
            if (!rValue){
                throw ModelException("FR::RValueBase::getRValue",
                                     "Internal error");
            }
        }
    }
    return rValue;
}

void FR::RValueBase::load(CClassSP& clazz){
    REGISTER(RValueBase, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(FRIfaces::IRValueExpression);
}

CClassConstSP const FR::LValueBase::TYPE = CClass::registerClassLoadMethod(
    "FR::LValueBase", typeid(LValueBase), load);

FR::LValueBase::LValueBase(const CClassConstSP& clazz):
    RValueBase(clazz){}

/** is the variable simulation date independent ie same value at
    each simulation date. This implementation returns false.
    Inheriting classes may need to override this */
bool FR::LValueBase::isSDI() const{
    return false;
}

/** Overrides implementation in RValueBase - hashes string
    returned by getID */
int FR::LValueBase::hashCode() const{
    return hash_string(getID());
}
        
/** Returns true if both objects have the same id string */
bool FR::LValueBase::equals(const FRIfaces::IRValueExpression* name) const{
    const LValueBase* nameBase = dynamic_cast<const LValueBase*>(name);
    return (nameBase && nameBase->getID() == getID());
}

/** gets FlexRule::ILValue from FlexRule which represents the value
    of this expression at the current timepoint. Implemented via
    cast of return object from getRValue. DO NOT FREE. */
FRIfaces::ILValue* FR::LValueBase::getLValue(FRController* frCtrl) const{
    return getLValue(frCtrl->getIndex(), frCtrl);
}
    
/** gets FRController::ILValue from FRController which represents the value
    of this expression at a specific timepoint. Implemented via
    cast of return object from getRValue */
FRIfaces::ILValue* FR::LValueBase::getLValue(int           index,
                                             FRController* frCtrl) const{
    FRIfaces::IRValue* rValue = getRValue(index, frCtrl);
    FRIfaces::ILValue* lValue = dynamic_cast<FRIfaces::ILValue*>(rValue);
    if (!lValue){
        throw ModelException("FR::LValueBase::getLValue", "Expression is "
                             "not an lvalue");
    }
    return lValue;
}

void FR::LValueBase::load(CClassSP& clazz){
    REGISTER(LValueBase, clazz);
    SUPERCLASS(RValueBase);
    IMPLEMENTS(FRIfaces::ILValueExpression);
}

/** Holds the set of variables for which a 'no op' should be 
    performed - useful when setting up 'sparse' algorithms */
class FR::TPNoOp: public CObject,
          public virtual FRIfaces::ITimePointNoOp{
private:
    const StringArrayConstSP  variables;
public:
    static CClassConstSP const TYPE;

    //// this object is immutable
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<TPNoOp*>(this);
    }

    virtual void validatePop2Object(){
        // validatePop2Object has privileged access
        StringArray& vars = const_cast<StringArray&>(*this->variables);
        // remove 'empty' variables
        for (vector<string>::iterator iter = vars.begin(); 
             iter != vars.end(); /* inc in loop body */){
            if (iter->empty() || *iter == ";"){
                iter = vars.erase(iter);
            } else{
                ++iter;
            }
        }
    }

    /** Returns the number of rules for this given timepoint ie the
        number of assignments */
    virtual int numVars() const{
        return variables->size();
    }
    
    /** Returns the lValue for the specified variable ie the
        variables name */
    virtual const string& lValueID(int ruleNum) const{
        if (ruleNum < 0 || ruleNum >= variables->size()){
            throw ModelException("FR::TPNoOp::lValueID",
                                 "index out of bounds");
        }
        return (*variables)[ruleNum];
    }
    
    virtual ~TPNoOp(){}

private:
    TPNoOp(): CObject(TYPE){}

    static IObject* defaultConstructor(){
        return new TPNoOp();
    }

    static void load(CClassSP& clazz){
        REGISTER(TPNoOp, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::ITimePointNoOp);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(variables, "variables not to compute");
        Addin::registerConstructor("FR_NO_OP",
                                   Addin::FLEX_PAYOFF,
                                   "Creates object listing variables not "
                                   "to be calculated at a given timepoint",
                                   TYPE);
    }    
};
CClassConstSP const FR::TPNoOp::TYPE = CClass::registerClassLoadMethod(
    "FR::TPNoOp", typeid(TPNoOp), load);


///////////////// TimePointRules class /////////
/** Holds set of "rules" or "assignments" for a given time point. 
    Basically holds set of lvalues which are to be set and a 
    corresponding set of rvalues with which to set them */
class FR::TimePointRules: public CObject, 
          virtual public FRIfaces::ITimePointRules{ 
public:
    static CClassConstSP const TYPE;

    virtual FRIfaces::IAssignmentArray* createAssignments (
        FRController* frCtrl) const{
        static const string routine("FR::TimePointRules::createAssignments");
        auto_ptr<FRIfaces::IAssignmentArray> assignments(
            new FRIfaces::IAssignmentArray());
        assignments->reserve(variables->size());
        FRIfaces::ILValueExpression* lValueExp;
        FRIfaces::IRValueExpression* rValueExp;
        string  parseError; // to build up collection of parse errors
        int     numErrors = 0;
        for (int i = 0; i < variables->size(); i++){
            lValueExp = (*variables)[i].get();
            rValueExp = (*expressions)[i].get();
            if (!lValueExp || !rValueExp){
                throw ModelException(routine,
                                     string((!lValueExp? 
                                             "variable": "expression"))+
                                     " number "+Format::toString(i+1) +
                                     "is null");
            }
            try{
                try {
                    FRIfaces::ILValue* lValue = lValueExp->getLValue(frCtrl);
                    FRIfaces::IRValue* rValue = rValueExp->getRValue(frCtrl);
                    FRIfaces::IAssignment* assign =
                        FRController::createAssignment(lValue, rValue);
                    assignments->push_back(FRIfaces::IAssignmentSP(assign));
                } catch (ModelException& e){
                    ModelException* cause = e.getCause();
                    FRParseException* parseEx = 
                        dynamic_cast<FRParseException*>(cause? cause:&e);
                    if (parseEx){
                        // record error and continue
                        parseError += string(e.what())+"\n"+
                            "On assignment to '"+lValueExp->getID()+ "' [" + 
                            FRFunction::argTypeToString(lValueExp->getType())+
                            "]\n";
                        numErrors += parseEx->getNumErrors();
                        if (numErrors > FRParseException::MAX_NUM_ERRORS){
                            throw FRParseException(numErrors,
                                                   parseError+"Max number of "
                                                   "errors exceeded");
                        }
                    } else {
                        throw;
                    }
                }
            } catch (exception& e) {
                // The type information here is not quite what we
                // want, but we don't have instances yet, just the
                // expressions.
                string m("Failed on assignment of expression: " + 
                         rValueExp->getID() + " [" + 
                         rValueExp->getClass()->getName() + 
                         "] to \nvariable: '" + 
                         lValueExp->getID() + "' [" + 
                         lValueExp->getClass()->getName() + "]");
                throw ModelException(e, routine, m);
            }
        }
        if (!parseError.empty()){
            // throw one big one
            throw FRParseException("Encountered following parse errors:\n"+
                                   parseError);
        }
        return assignments.release();
    }
    
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<FR::TimePointRules*>(this);
    }
    
    virtual void validatePop2Object(){
        if (variables->size() != expressions->size()){
            throw ModelException("TimePointRules::validatePop2Object",
                                 Format::toString(variables->size())+
                                 " variables but "+
                                 Format::toString(expressions->size())+
                                 " expressions");
        }
    }

    /** Returns the number of rules for this given timepoint ie the
        number of assignments */
    virtual int numRules() const{
        return variables->size();
    }
    
    /** Returns the lValue for the specified rule as a string ie the
        variable name */
    virtual string lValueID(int ruleNum) const{
        checkIndex(ruleNum);
        return ((*variables)[ruleNum]->getID());
    }
    
    /** Returns the rValue for the specified rule as a string eg X+Y */
    virtual string rValueID(int ruleNum) const{
        checkIndex(ruleNum);
        return ((*expressions)[ruleNum]->getID());
    };

private:
    void checkIndex(int ruleNum) const{
        if (ruleNum < 0 || ruleNum >= variables->size()){
            throw ModelException("TimePointRules::checkIndex", "Index out "
                                 "of bounds");
        }
    }

    // fields 
    const FRIfaces::ILValueExpressionArraySP   variables;
    const FRIfaces::IRValueExpressionArraySP   expressions;


    TimePointRules(): CObject(TYPE) {};
    
    
    static IObject* defaultTimePointRules(){
        return new TimePointRules();
    }
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(TimePointRules, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::ITimePointRules);
        EMPTY_SHELL_METHOD(defaultTimePointRules);
        FIELD(variables,             "variables");
        FIELD(expressions,           "expressions");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_TIMEPOINT_RULES",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX TimePointRules",
                                   TYPE);
    }
};
 
CClassConstSP const FR::TimePointRules::TYPE = CClass::registerClassLoadMethod(
    "FR::TimePointRules", typeid(TimePointRules), load);

// work around for msvc 7 bug
typedef FR::TimePointRulesArray FRTimePointRulesArray;
DEFINE_TEMPLATE_TYPE_WITH_NAME("FR::TimePointRulesArray", FRTimePointRulesArray);

/** Simple implementation of FRIfaces::IAlgorithm. Just an array of
    time point rules - for one each simulation date */
class FR::TPRSimpleArray: public CObject,
          virtual public FRIfaces::IAlgorithm {
    
private:
    // fields
    const FRIfaces::ITimePointRulesArraySP              algorithm;
    // registration
    static void load(CClassSP& clazz){
        REGISTER(TPRSimpleArray, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::IAlgorithm);
        EMPTY_SHELL_METHOD(defaultTPRSimpleArray);
        FIELD(algorithm,   "algorithm");
        CObject::registerObjectFromArrayMethod(
            TimePointRulesArray::TYPE, TYPE, fromArray);
        clazz->setPublic(); 
        Addin::registerConstructor("FR_SIMPLE_ALGORITHM",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a simple FLEX Algorithm",
                                   TYPE);
    }

    TPRSimpleArray(): CObject(TYPE){}

    static IObject* defaultTPRSimpleArray(){
        return new TPRSimpleArray();
    }

    /** this object is immutable */
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<TPRSimpleArray*>(this);
    }

    /** for type conversion - copes with array of time point rules and builds
        a TPRSimpleArray (for backwards compatibility) */
    static IObjectSP fromArray(const IObjectSP&  object,
                               CClassConstSP     requiredType){
        // convert from object
        TimePointRulesArray& rules =
            dynamic_cast<TimePointRulesArray&>(*object);
        // create space for new array
        FRIfaces::ITimePointRulesArraySP newRules(
            new FRIfaces::ITimePointRulesArray(rules.size()));
        for (int i = 0; i < rules.size(); i++){
            (*newRules)[i] = 
                FRIfaces::ITimePointRulesSP::dynamicCast(rules[i]);
        }
        return IObjectSP(new TPRSimpleArray(newRules));
    }

public:
    static CClassConstSP const TYPE;

    // made public to allow FlexInstrument to operate
    TPRSimpleArray(const FRIfaces::ITimePointRulesArraySP& rules): 
        CObject(TYPE), algorithm(rules){}

    virtual ~TPRSimpleArray(){}

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        const FRIfaces::ITimePointRules* rules = 
            (*algorithm)[simDateIndex].get();
        return (rules->lValueID(assignmentIndex) +" = "+
                rules->rValueID(assignmentIndex));
    }

    /** create assignment array for current time point */
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const{
        if (frCtrl->numDates() != algorithm->size()){
            throw ModelException("FR::TPRSimpleArray", 
                                 "Number of simDates and number of"
                                 " TimePointRules must be the same");
        }
        const FRIfaces::ITimePointRulesArray& algo = *algorithm;
        int  i = frCtrl->getIndex();
        if (!algo[i]){
            throw ModelException("FR::TPRSimpleArray", "null set of rules at "
                                 +frCtrl->getDate(i).toString());
        }
        return (algo[i]->createAssignments(frCtrl));
    }

    // May be better done on ProductMC ?
    virtual void writeRules(const string& fileName,
                            const string& payoutVariable) const {
        // Loop over time points
        // At each time point loop over variables
        // For each variable, if a different rule print it
        // Identify the payout variable
#ifdef sun
        /** There is a bug in the assember when dealing with very long
            symbols which shows up with the map<string,string> */
        throw ModelException("writeRules not supported on solaris");
#else
        ofstream                            ruleStream(fileName.c_str());
        /* currentRules: keyed by variable ID, storing expression ID */      
        map<string, string>                 currentRules; 
        map<string, string>::const_iterator r; 
        bool                                printMe;
        bool                                printTPHeader;

        for (int t = 0; t < algorithm->size(); t++){
            // By keeping a local record of variables with rules at
            // this time point we can report those variables which
            // have been left undefined.
            map<string, string>   undefinedRules(currentRules);
            printTPHeader = true;
            int numRules = (*algorithm)[t]->numRules();
            for (int i = 0; i < numRules; i++){
                string lValueExp = (*algorithm)[t]->lValueID(i);
                string rValueExp = (*algorithm)[t]->rValueID(i);
                printMe = true;
                r = currentRules.find(lValueExp);
                if (r == currentRules.end()){
                    // variable not found so it's got a new rule this
                    // time point: add and print it
                    currentRules[lValueExp] = rValueExp;
                } else {
                    // is this time point's rule the same?
                    if (r->second == rValueExp) {
                        // is the same so don't write it
                        printMe = false;
                    }
                    // it's defined at this time point so
                    // remove from the list of undefined
                    undefinedRules.erase(lValueExp);
                    currentRules[lValueExp] = rValueExp;
                }
                if (printMe) {
                    if (printTPHeader) {
                        printTPHeader = false;
                        ruleStream << "\nTime point " << t << "\n" <<
                            "-------------------\n";
                    }
                    ruleStream << lValueExp + " = " + rValueExp + "\n";
                }
            }
            // The convention is that rules are printed when they change. 
            // Report any vars undefined at t.
            for (r = undefinedRules.begin(); r != undefinedRules.end(); ++r) {
                if (printTPHeader) {
                    printTPHeader = false;
                    ruleStream << "\nTime point " << t << "\n" << 
                        "-------------------\n";
                }
                ruleStream << r->first + " is undefined.\n";
                currentRules.erase(r->first);
            }
        }
        ruleStream << "\n===============\nPayout variable is " +
            payoutVariable + "\n";
#endif
    }

    /** Report all dates known to the algorithm */
    DateTimeArraySP getAllDates() const {
        return DateTimeArraySP(0);
    }

    StringArrayArray getIMSInput() const {
        return StringArrayArray(0);
    }

    // do nothing
    virtual void getIMSInput(
        IntArraySP    simDateRuleIds,
        IntArraySP    ruleId,
        StringArraySP ruleDefn) const {}
};

CClassConstSP const FR::TPRSimpleArray::TYPE = CClass::registerClassLoadMethod(
    "FR::TPRSimpleArray", typeid(TPRSimpleArray), load);

// Not currently functional because FRAladdinAlgorithm is declared in FRParserTPRules.cpp
//smartPtr<FRAladdinAlgorithm> FlexAlgorithm::getAlgorithm(const DateTimeArraySP simDates) {
//    static const string routine = "FlexAlgorithm::getAlgorithm";

//    if (simDates->size() != simDateRuleIds.size()) {
//        throw ModelException(routine,
//                             "#simDates (" + Format::toString(simDates->size()) +
//                             ") does not equal #sim date rules (" +
//                             Format::toString(simDateRuleIds.size()) + ")");
//    }

    // Need to decrement the simDateRuleIds
//    IntArray rulesIndex(simDateRuleIds.size());
//    for (int i=0; i<simDateRuleIds.size(); i++) {
//        rulesIndex[i] = simDateRuleIds[i] - 1;
//    }

//    smartPtr<FRAladdinAlgorithmSP> algorithm(new FRAladdinAlgorithm(simDates,
//                                                                    rulesIndex,
//                                                                    justTheRules));
//    return algorithm;    

//}


FRIfaces::IAlgorithmSP FlexAlgorithm::getAlgorithm(const DateTimeArraySP simDates) {
    static const string routine = "FlexAlgorithm::getAlgorithm";
    
    if (simDates->size() != simDateRuleIds.size()) {
        throw ModelException(routine,
                             "#simDates (" + Format::toString(simDates->size()) +
                             ") does not equal #sim date rules (" +
                             Format::toString(simDateRuleIds.size()) + ")");
    }
    // Arrange rules according to simulation dates
    FRIfaces::ITimePointRulesArraySP rulesPerSimDate(
        new FRIfaces::ITimePointRulesArray(simDateRuleIds.size()));
    for(int i=0;i<simDateRuleIds.size();i++) {
        // "-1" since simDateRuleIds base from 1, not 0
        (*rulesPerSimDate)[i] = (*justTheRules)[simDateRuleIds[i]-1];
    }
    smartPtr<FR::TPRSimpleArray> algorithm(new FR::TPRSimpleArray(rulesPerSimDate));
    return algorithm;
};

/** Just wrapper around getValue method */
IObjectConstSP FR::RValueDouble::get(){
    return IObjectConstSP(CDouble::create(getValue()));
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueDouble::isKnown() const{ 
    return false;
}

FR::RValueDouble::~RValueDouble(){}

/** Just wrapper around getValue method */
IObjectConstSP FR::RValueInt::get(){
    return IObjectConstSP(CInt::create(getValue()));
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueInt::isKnown() const{ 
    return false;
}

FR::RValueInt::~RValueInt(){}

/** Just wrapper around getValue method */
IObjectConstSP FR::RValueBool::get(){
    return IObjectConstSP(CBool::create(getValue()));
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueBool::isKnown() const{ 
    return false;
}

FR::RValueBool::~RValueBool(){}

/** Just wrapper around getValue method */
IObjectConstSP FR::RValueDate::get(){
    return IObjectConstSP::attachToRef(&getValue());
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueDate::isKnown() const{ 
    return false;
}

FR::RValueDate::~RValueDate(){}
/** Just wrapper around getValue method */
IObjectConstSP FR::RValueDoubleArray::get(){
    IRValueDoubleArray::RT* rt = getRT(); // get hold of run time object
    int size = rt->size(rt);
    DoubleArraySP dbArray(new DoubleArray(size));
    for (int i = 0; i < size; i++){
        (*dbArray)[i] = rt->func(rt, i);
    }
    return IObjectSP(dbArray);
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueDoubleArray::isKnown() const{ 
    return false;
}

FR::RValueDoubleArray::~RValueDoubleArray(){}

/** Just wrapper around getValue method */
IObjectConstSP FR::RValueIntArray::get(){
    IRValueIntArray::RT* rt = getRT(); // get hold of run time object
    int size = rt->size(rt);
    IntArraySP dbArray(new IntArray(size));
    for (int i = 0; i < size; i++){
        (*dbArray)[i] = rt->func(rt, i);
    }
    return IObjectSP(dbArray);
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueIntArray::isKnown() const{ 
    return false;
}

FR::RValueIntArray::~RValueIntArray(){}


/** Just wrapper around getValue method */
IObjectConstSP FR::RValueBoolArray::get(){
    IRValueBoolArray::RT* rt = getRT(); // get hold of run time object
    int size = rt->size(rt);
    CBoolArraySP dbArray(new BoolArray(size));
    for (int i = 0; i < size; i++){
        (*dbArray)[i] = rt->func(rt, i);
    }
    return IObjectSP(dbArray);
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RValueBoolArray::isKnown() const{ 
    return false;
}

FR::RValueBoolArray::~RValueBoolArray(){}

/** Utility Class for a constant double value */
FR::RConstDouble::MyRT::MyRT(double  value): func(&getValue), value(value){}

double FR::RConstDouble::MyRT::getValue(void* structure){
    MyRT* rt = (MyRT*)structure;
    return (rt->value);
}
// simple constructor
FR::RConstDouble::RConstDouble(double value): 
    rt(new MyRT(value)){}

// get the variable expressed as a double
double FR::RConstDouble::getValue() {
    return rt->value;
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RConstDouble::isKnown() const{ 
    return true;
}
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDouble::RT* FR::RConstDouble::getRT(){
    return (FRIfaces::IRValueDouble::RT*) rt;
}

FR::RConstDouble::~RConstDouble(){
    delete rt;
}

/** Utility Class for a constant int value */
/** Utility Class for a constant double value */
FR::RConstInt::MyRT::MyRT(int  value): func(&getValue), value(value){}

int FR::RConstInt::MyRT::getValue(void* structure){
    MyRT* rt = (MyRT*)structure;
    return (rt->value);
}
// simple constructor
FR::RConstInt::RConstInt(int value): 
    rt(new MyRT(value)){}

// get the variable expressed as a int
int FR::RConstInt::getValue() {
    return rt->value;
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RConstInt::isKnown() const{ 
    return true;
}
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueInt::RT* FR::RConstInt::getRT(){
    return (FRIfaces::IRValueInt::RT*) rt;
}

FR::RConstInt::~RConstInt(){
    delete rt;
}

/** Utility Class for a constant bool value */
FR::RConstBool::MyRT::MyRT(bool  value): func(&getValue), value(value){}

bool FR::RConstBool::MyRT::getValue(void* structure){
    MyRT* rt = (MyRT*)structure;
    return (rt->value);
}
// simple constructor
FR::RConstBool::RConstBool(bool value): 
    rt(new MyRT(value)){}

// get the variable expressed as a bool
bool FR::RConstBool::getValue() {
    return rt->value;
}
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueBool::RT* FR::RConstBool::getRT(){
    return (FRIfaces::IRValueBool::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RConstBool::isKnown() const{ 
    return true;
}

FR::RConstBool::~RConstBool(){
    delete rt;
}

/** Utility Class for a constant date value */
// simple constructor
FR::RConstDate::RConstDate(const DateTime::Date& value): 
    value(value) {}

// get the variable expressed as a bool
const DateTime::Date& FR::RConstDate::getValue() {
    return value;
}
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDate::RT* FR::RConstDate::getRT(){
    return 0; // cheat for now
}

/** Is the value of this object known before the simulation starts
    eg a constant. */
bool FR::RConstDate::isKnown() const{ 
    return true;
}

FR::RConstDate::~RConstDate(){}

//// uses FR::MemMgr
void* FR::LValueDouble::RT::operator new(size_t size){
    return FR::MemMgr::alloc(size);
}
void FR::LValueDouble::RT::operator delete(void *ptr){
    FR::MemMgr::dealloc(ptr);
}

/** generates exception saying either value not set yet or 
    already set */
ModelException FR::LValueDouble::RT::makeException(){
    if (*isSet){
        string m("The value of the variable "+string(name)+
                 " has already been set");
        return ModelException("LValueDouble::setValue", m);
    }
    string m("The value of the variable "+string(name)+
             " has not been set");
    return ModelException("LValueDouble::getValue", m);
}
// constructor for no initial value
FR::LValueDouble::RT::RT(char* isSet, const char* name):
    func(&getValue), setFunc(&setValue), isSet(isSet), name(name){
    *isSet = 0;
}

FR::LValueDouble::~LValueDouble(){
    delete rt;
}

// constructor for initial value
FR::LValueDouble::LValueDouble(const char* name): rt(new RT(&isSet, name)){}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDouble::RT* FR::LValueDouble::getRT(){
    return (FRIfaces::IRValueDouble::RT*) rt;
}
void FR::LValueDouble::setReset(char* reset){
    rt->isSet = reset;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueDouble::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as double)  */
void FR::LValueDouble::set(const IObjectConstSP& object){
    if (CDouble::TYPE->isInstance(object)){
        const CDouble* theDb = STATIC_CONST_CAST(CDouble, object.get());
        setValue(theDb->doubleValue());
    } else if (CInt::TYPE->isInstance(object)){
        const CInt* theDb = STATIC_CONST_CAST(CInt, object.get());
        setValue(theDb->intValue());
    } else {
        throw FRParseException("Cannot set double from"
                               " type "+object->getClass()->getName());
    }
}

// get the variable expressed as a double
double FR::LValueDouble::getValue() {
    if (!*rt->isSet){
        throw rt->makeException();
    }
    return rt->value;
}

// set the variable using a double
void FR::LValueDouble::setValue(double value) {
    if (*rt->isSet){
        throw rt->makeException();
    }
    rt->value = value;
    *rt->isSet = 1;
}

//// 'derived' from FRIfaces::ILValueDoubleArray::RT
struct FR::LValueDoubleArray::RT{
    TGetSize*         size;
    TGetValue*        func;
    TSetValue*        setFunc;
    char*             isSet; /* currently applies to all elements of
                                the array */
    const char*       name;
    vector<double>    values;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }

    /** generates exception saying either value not set yet or 
        already set */
    ModelException makeException(){
        if (*isSet){
            string m("The value of the variable "+string(name)+
                     " has already been set");
            return ModelException("LValueDoubleArray::setValue", m);
        }
        string m("The value of the variable "+string(name)+
                 " has not been set");
        return ModelException("LValueDoubleArray::getValue", m);
    }
    static int getSize(void* structure){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values.size());
    }

    static double getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values[index]);
    }

    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueDoubleArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (*rt->isSet){
            throw rt->makeException();
        }
        if (index < 0){
            // we copy the array over (we support dynamic length arrays)
            int size = rValue->size(rValue);
            rt->values.resize(size);
            for (int i = 0; i < size; i++){
                rt->values[i] = rValue->func(rValue, i);
            }
        } else {
            throw ModelException("FR::LValueDoubleArray::set", "Setting "
                                 "individual elements not supported yet");
        }
        *rt->isSet = 1;
    }

    // constructor for no initial value
    explicit RT(char* isSet, const char* name): 
        size(&getSize), func(&getValue), 
        setFunc(&setValue), isSet(isSet), name(name){
        *isSet = 0;
    }
};

FR::LValueDoubleArray::~LValueDoubleArray(){
    delete rt;
}
// constructor for no initial value
FR::LValueDoubleArray::LValueDoubleArray(const char* name):
    rt(new RT(&isSet, name)){}

//// allows central pooling of which variables are known
void FR::LValueDoubleArray::setReset(char* reset){
    rt->isSet = reset;
}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDoubleArray::RT* FR::LValueDoubleArray::getRT(){
    return (FRIfaces::IRValueDoubleArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueDoubleArray::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as double)  */
void FR::LValueDoubleArray::set(const IObjectConstSP& object){
    const DoubleArray* dbArray = DYNAMIC_CONST_CAST(DoubleArray, object.get());
    rt->values = vector<double>(dbArray->begin(), dbArray->end());
    *rt->isSet = 1;
}

//// 'derived' from FRIfaces::ILValueDoubleArray::RT
struct FR::LValueDoubleVarArray::RT{
    TGetSize*         size;
    TGetValue*        func;
    TSetValue*        setFunc;
    const int         length; // fixed
    FRIfaces::ILValueDouble::RT**  lValuesRT;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        if (ptr){
            FR::MemMgr::dealloc(((RT*)ptr)->lValuesRT);
        }
        FR::MemMgr::dealloc(ptr);
    }
    //// Returns length of array
    static int getSize(void* structure){
        RT* rt = (RT*)structure;
        return (rt->length);
    }

    //// Returns value of specified element
    static double getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        FRIfaces::ILValueDouble::RT* dbRT = rt->lValuesRT[index];
        return (dbRT->func(dbRT));
    }

    //// Sets value of specified element/all elements
    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueDoubleArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (index < 0){
            int size = rValue->size(rValue);
            if (size != rt->length){
                throw ModelException("FR::LValueDoubleVarArray::RT::setValue",
                                     "Trying to assign an array of length "+
                                     Format::toString(size)+" to an array of "
                                     "variables of length "+
                                     Format::toString(rt->length));
            }
            for (int i = 0; i < size; i++){
                double value = rValue->func(rValue, i);
                FRIfaces::ILValueDouble::RT* dbRT = rt->lValuesRT[i];
                dbRT->setFunc(dbRT, value);
            }
        } else {
            double value = rValue->func(rValue, index);
            FRIfaces::ILValueDouble::RT* dbRT = rt->lValuesRT[index];
            dbRT->setFunc(dbRT, value);
        }
    }

    explicit RT(const vector<FRIfaces::ILValueDouble*>& vectorLValues): 
        size(&getSize), func(&getValue), setFunc(&setValue),
        length(vectorLValues.size()),
        lValuesRT((FRIfaces::ILValueDouble::RT**)
                  FR::MemMgr::alloc(length * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < length; i++){
            FRIfaces::ILValueDouble* lValue = vectorLValues[i];
            lValuesRT[i] = (FRIfaces::ILValueDouble::RT*)lValue->getRT();
        }
    }
};

FR::LValueDoubleVarArray::~LValueDoubleVarArray(){
    delete rt;
}
//// constructor (pass in array of variables to use for array)
FR::LValueDoubleVarArray::LValueDoubleVarArray(
    const vector<FRIfaces::ILValueDouble*>& vectorLValues):
    rt(new RT(vectorLValues)), vectorLValues(vectorLValues){}

//// allows central pooling of which variables are known. Does nothing
//// because we delegate that to component variables
void FR::LValueDoubleVarArray::setReset(char* reset){}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDoubleArray::RT* FR::LValueDoubleVarArray::getRT(){
    return (FRIfaces::IRValueDoubleArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueDoubleVarArray::isKnown() const{
    for (unsigned int i = 0; i < vectorLValues.size(); i++){
        if (!vectorLValues[i]->isKnown()){
            return false;
        }
    }
    return true;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as double)  */
void FR::LValueDoubleVarArray::set(const IObjectConstSP& object){
    const DoubleArray* dbArray = DYNAMIC_CONST_CAST(DoubleArray, object.get());
    if (rt->length != dbArray->size()){
        throw ModelException("FR::LValueDoubleVarArray::set", "Trying to "
                             "assign an array of length "+
                             Format::toString(dbArray->size()) + " to an "
                             "array of double variables of length "+
                             Format::toString(rt->length));
    }
    for (int i = 0; i < rt->length; i++){
        FRIfaces::ILValueDouble::RT* dbRT = rt->lValuesRT[i];
        dbRT->setFunc(dbRT, (*dbArray)[i]);
    }
}

//// 'derived' from FRIfaces::ILValueIntArray::RT
struct FR::LValueIntArray::RT{
    TGetSize*         size;
    TGetValue*        func;
    TSetValue*        setFunc;
    char*             isSet; /* currently applies to all elements of
                                the array */
    const char*       name;
    vector<int>       values;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }

    /** generates exception saying either value not set yet or 
        already set */
    ModelException makeException(){
        if (*isSet){
            string m("The value of the variable "+string(name)+
                     " has already been set");
            return ModelException("LValueIntArray::setValue", m);
        }
        string m("The value of the variable "+string(name)+
                 " has not been set");
        return ModelException("LValueIntArray::getValue", m);
    }
    static int getSize(void* structure){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values.size());
    }

    static int getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values[index]);
    }

    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueIntArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (*rt->isSet){
            throw rt->makeException();
        }
        if (index < 0){
            // we copy the array over (we support dynamic length arrays)
            int size = rValue->size(rValue);
            rt->values.resize(size);
            for (int i = 0; i < size; i++){
                rt->values[i] = rValue->func(rValue, i);
            }
        } else {
            throw ModelException("FR::LValueIntArray::set", "Setting "
                                 "individual elements not supported yet");
        }
        *rt->isSet = 1;
    }

    // constructor for no initial value
    explicit RT(char* isSet, const char* name): 
        size(&getSize), func(&getValue), 
        setFunc(&setValue), isSet(isSet), name(name){
        *isSet = 0;
    }
};

FR::LValueIntArray::~LValueIntArray(){
    delete rt;
}
// constructor for no initial value
FR::LValueIntArray::LValueIntArray(const char* name): 
    rt(new RT(&isSet, name)){}

//// allows central pooling of which variables are known
void FR::LValueIntArray::setReset(char* reset){
    rt->isSet = reset;
}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueIntArray::RT* FR::LValueIntArray::getRT(){
    return (FRIfaces::IRValueIntArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueIntArray::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as int)  */
void FR::LValueIntArray::set(const IObjectConstSP& object){
    const IntArray* dbArray = DYNAMIC_CONST_CAST(IntArray, object.get());
    rt->values = vector<int>(dbArray->begin(), dbArray->end());
    *rt->isSet = 1;
}

//// 'derived' from FRIfaces::ILValueIntArray::RT
struct FR::LValueIntVarArray::RT{
    TGetSize*                   size;
    TGetValue*                  func;
    TSetValue*                  setFunc;
    const int                   length; // fixed
    FRIfaces::ILValueInt::RT**  lValuesRT;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        if (ptr){
            FR::MemMgr::dealloc(((RT*)ptr)->lValuesRT);
        }
        FR::MemMgr::dealloc(ptr);
    }
    //// Returns length of array
    static int getSize(void* structure){
        RT* rt = (RT*)structure;
        return (rt->length);
    }

    //// Returns value of specified element
    static int getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        FRIfaces::ILValueInt::RT* dbRT = rt->lValuesRT[index];
        return (dbRT->func(dbRT));
    }

    //// Sets value of specified element/all elements
    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueIntArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (index < 0){
            int size = rValue->size(rValue);
            if (size != rt->length){
                throw ModelException("FR::LValueIntVarArray::RT::setValue",
                                     "Trying to assign an array of length "+
                                     Format::toString(size)+" to an array of "
                                     "variables of length "+
                                     Format::toString(rt->length));
            }
            for (int i = 0; i < size; i++){
                int value = rValue->func(rValue, i);
                FRIfaces::ILValueInt::RT* dbRT = rt->lValuesRT[i];
                dbRT->setFunc(dbRT, value);
            }
        } else {
            int value = rValue->func(rValue, index);
            FRIfaces::ILValueInt::RT* dbRT = rt->lValuesRT[index];
            dbRT->setFunc(dbRT, value);
        }
    }

    explicit RT(const vector<FRIfaces::ILValueInt*>& vectorLValues): 
        size(&getSize), func(&getValue), setFunc(&setValue),
        length(vectorLValues.size()),
        lValuesRT((FRIfaces::ILValueInt::RT**)
                  FR::MemMgr::alloc(length * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < length; i++){
            FRIfaces::ILValueInt* lValue = vectorLValues[i];
            lValuesRT[i] = (FRIfaces::ILValueInt::RT*)lValue->getRT();
        }
    }
};

FR::LValueIntVarArray::~LValueIntVarArray(){
    delete rt;
}
//// constructor (pass in array of variables to use for array)
FR::LValueIntVarArray::LValueIntVarArray(
    const vector<FRIfaces::ILValueInt*>& vectorLValues):
    rt(new RT(vectorLValues)), vectorLValues(vectorLValues){}

//// allows central pooling of which variables are known. Does nothing
//// because we delegate that to component variables
void FR::LValueIntVarArray::setReset(char* reset){}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueIntArray::RT* FR::LValueIntVarArray::getRT(){
    return (FRIfaces::IRValueIntArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueIntVarArray::isKnown() const{
    for (unsigned int i = 0; i < vectorLValues.size(); i++){
        if (!vectorLValues[i]->isKnown()){
            return false;
        }
    }
    return true;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as int)  */
void FR::LValueIntVarArray::set(const IObjectConstSP& object){
    const IntArray* dbArray = DYNAMIC_CONST_CAST(IntArray, object.get());
    if (rt->length != dbArray->size()){
        throw ModelException("FR::LValueIntVarArray::set", "Trying to "
                             "assign an array of length "+
                             Format::toString(dbArray->size()) + " to an "
                             "array of int variables of length "+
                             Format::toString(rt->length));
    }
    for (int i = 0; i < rt->length; i++){
        FRIfaces::ILValueInt::RT* dbRT = rt->lValuesRT[i];
        dbRT->setFunc(dbRT, (*dbArray)[i]);
    }
}

//// 'derived' from FRIfaces::ILValueBoolArray::RT
struct FR::LValueBoolArray::RT{
    TGetSize*         size;
    TGetValue*        func;
    TSetValue*        setFunc;
    char*             isSet; /* currently applies to all elements of
                                the array */
    const char*       name;
    vector<bool>      values;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }

    /** generates exception saying either value not set yet or 
        already set */
    ModelException makeException(){
        if (*isSet){
            string m("The value of the variable "+string(name)+
                     " has already been set");
            return ModelException("LValueBoolArray::setValue", m);
        }
        string m("The value of the variable "+string(name)+
                 " has not been set");
        return ModelException("LValueBoolArray::getValue", m);
    }
    static int getSize(void* structure){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values.size());
    }

    static bool getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        if (!*rt->isSet){
            throw rt->makeException();
        }
        return (rt->values[index]);
    }

    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueBoolArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (*rt->isSet){
            throw rt->makeException();
        }
        if (index < 0){
            // we copy the array over (we support dynamic length arrays)
            int size = rValue->size(rValue);
            rt->values.resize(size);
            for (int i = 0; i < size; i++){
                rt->values[i] = rValue->func(rValue, i);
            }
        } else {
            throw ModelException("FR::LValueBoolArray::set", "Setting "
                                 "individual elements not supported yet");
        }
        *rt->isSet = 1;
    }

    // constructor for no initial value
    explicit RT(char* isSet, const char* name): 
        size(&getSize), func(&getValue), 
        setFunc(&setValue), isSet(isSet), name(name){
        *isSet = 0;
    }
};

FR::LValueBoolArray::~LValueBoolArray(){
    delete rt;
}
// constructor for no initial value
FR::LValueBoolArray::LValueBoolArray(const char* name):
    rt(new RT(&isSet, name)){}

//// allows central pooling of which variables are known
void FR::LValueBoolArray::setReset(char* reset){
    rt->isSet = reset;
}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueBoolArray::RT* FR::LValueBoolArray::getRT(){
    return (FRIfaces::IRValueBoolArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueBoolArray::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as bool)  */
void FR::LValueBoolArray::set(const IObjectConstSP& object){
    const BoolArray* dbArray = DYNAMIC_CONST_CAST(BoolArray, object.get());
    rt->values = vector<bool>(dbArray->begin(), dbArray->end());
    *rt->isSet = 1;
}

//// 'derived' from FRIfaces::ILValueBoolArray::RT
struct FR::LValueBoolVarArray::RT{
    TGetSize*                    size;
    TGetValue*                   func;
    TSetValue*                   setFunc;
    const int                    length; // fixed
    FRIfaces::ILValueBool::RT**  lValuesRT;

    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        if (ptr){
            FR::MemMgr::dealloc(((RT*)ptr)->lValuesRT);
        }
        FR::MemMgr::dealloc(ptr);
    }
    //// Returns length of array
    static int getSize(void* structure){
        RT* rt = (RT*)structure;
        return (rt->length);
    }

    //// Returns value of specified element
    static bool getValue(void* structure, int index){
        RT* rt = (RT*) structure;
        FRIfaces::ILValueBool::RT* dbRT = rt->lValuesRT[index];
        return (dbRT->func(dbRT));
    }

    //// Sets value of specified element/all elements
    static void setValue(void*                             structure, 
                         int                               index,
                         FRIfaces::IRValueBoolArray::RT* rValue){
        RT* rt = (RT*) structure;
        if (index < 0){
            int size = rValue->size(rValue);
            if (size != rt->length){
                throw ModelException("FR::LValueBoolVarArray::RT::setValue",
                                     "Trying to assign an array of length "+
                                     Format::toString(size)+" to an array of "
                                     "variables of length "+
                                     Format::toString(rt->length));
            }
            for (int i = 0; i < size; i++){
                bool value = rValue->func(rValue, i);
                FRIfaces::ILValueBool::RT* dbRT = rt->lValuesRT[i];
                dbRT->setFunc(dbRT, value);
            }
        } else {
            bool value = rValue->func(rValue, index);
            FRIfaces::ILValueBool::RT* dbRT = rt->lValuesRT[index];
            dbRT->setFunc(dbRT, value);
        }
    }

    explicit RT(const vector<FRIfaces::ILValueBool*>& vectorLValues): 
        size(&getSize), func(&getValue), setFunc(&setValue),
        length(vectorLValues.size()),
        lValuesRT((FRIfaces::ILValueBool::RT**)
                  FR::MemMgr::alloc(length * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < length; i++){
            FRIfaces::ILValueBool* lValue = vectorLValues[i];
            lValuesRT[i] = (FRIfaces::ILValueBool::RT*)lValue->getRT();
        }
    }
};

FR::LValueBoolVarArray::~LValueBoolVarArray(){
    delete rt;
}
//// constructor (pass in array of variables to use for array)
FR::LValueBoolVarArray::LValueBoolVarArray(
    const vector<FRIfaces::ILValueBool*>& vectorLValues):
    rt(new RT(vectorLValues)), vectorLValues(vectorLValues){}

//// allows central pooling of which variables are known. Does nothing
//// because we delegate that to component variables
void FR::LValueBoolVarArray::setReset(char* reset){}
        
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueBoolArray::RT* FR::LValueBoolVarArray::getRT(){
    return (FRIfaces::IRValueBoolArray::RT*) rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueBoolVarArray::isKnown() const{
    for (unsigned int i = 0; i < vectorLValues.size(); i++){
        if (!vectorLValues[i]->isKnown()){
            return false;
        }
    }
    return true;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as bool)  */
void FR::LValueBoolVarArray::set(const IObjectConstSP& object){
    const BoolArray* dbArray = DYNAMIC_CONST_CAST(BoolArray, object.get());
    if (rt->length != dbArray->size()){
        throw ModelException("FR::LValueBoolVarArray::set", "Trying to "
                             "assign an array of length "+
                             Format::toString(dbArray->size()) + " to an "
                             "array of bool variables of length "+
                             Format::toString(rt->length));
    }
    for (int i = 0; i < rt->length; i++){
        FRIfaces::ILValueBool::RT* dbRT = rt->lValuesRT[i];
        dbRT->setFunc(dbRT, (*dbArray)[i]);
    }
}

//// uses FR::MemMgr
void* FR::LValueInt::RT::operator new(size_t size){
    return FR::MemMgr::alloc(size);
}
void FR::LValueInt::RT::operator delete(void *ptr){
    FR::MemMgr::dealloc(ptr);
}

/** generates exception saying either value not set yet or 
    already set */
ModelException FR::LValueInt::RT::makeException(){
    if (*isSet){
        string m("The value of the variable "+string(name)+
                 " has already been set");
        return ModelException("LValueInt::setValue", m);
    }
    string m("The value of the variable "+string(name)+
             " has not been set");
    return ModelException("LValueInt::getValue", m);
}
// constructor for no initial value
FR::LValueInt::RT::RT(char* isSet, const char* name):
    func(&getValue), setFunc(&setValue), isSet(isSet), name(name){
    *isSet = 0;
}

FR::LValueInt::~LValueInt(){
    delete rt;
}

// constructor for initial value
FR::LValueInt::LValueInt(const char* name): rt(new RT(&isSet, name)){}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueInt::RT* FR::LValueInt::getRT(){
    return (FRIfaces::IRValueInt::RT*) rt;
}
void FR::LValueInt::setReset(char* reset){
    rt->isSet = reset;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueInt::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as int)  */
void FR::LValueInt::set(const IObjectConstSP& object){
    if (CInt::TYPE->isInstance(object)){
        const CInt* theDb = STATIC_CONST_CAST(CInt, object.get());
        setValue(theDb->intValue());
    } else {
        throw FRParseException("Cannot set int from"
                               " type "+object->getClass()->getName());
    }
}

// get the variable expressed as a int
int FR::LValueInt::getValue() {
    if (!*rt->isSet){
        throw rt->makeException();
    }
    return rt->value;
}

// set the variable using a int
void FR::LValueInt::setValue(int value) {
    if (*rt->isSet){
        throw rt->makeException();
    }
    rt->value = value;
    *rt->isSet = 1;
}

//// uses FR::MemMgr
void* FR::LValueBool::RT::operator new(size_t size){
    return FR::MemMgr::alloc(size);
}
void FR::LValueBool::RT::operator delete(void *ptr){
    FR::MemMgr::dealloc(ptr);
}

/** generates exception saying either value not set yet or 
    already set */
ModelException FR::LValueBool::RT::makeException(){
    if (*isSet){
        string m("The value of the variable "+string(name)+
                 " has already been set");
        return ModelException("LValueBool::setValue", m);
    }
    string m("The value of the variable "+string(name)+
             " has not been set");
    return ModelException("LValueBool::getValue", m);
}
// constructor for no initial value
FR::LValueBool::RT::RT(char* isSet, const char* name):
    func(&getValue), setFunc(&setValue), isSet(isSet), name(name){
    *isSet = 0;
}

FR::LValueBool::~LValueBool(){
    delete rt;
}

// constructor for initial value
FR::LValueBool::LValueBool(const char* name): rt(new RT(&isSet, name)){}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueBool::RT* FR::LValueBool::getRT(){
    return (FRIfaces::IRValueBool::RT*) rt;
}
void FR::LValueBool::setReset(char* reset){
    rt->isSet = reset;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made). Here this is driven off whether the variable is set
    or not. This allows other variables to be precalculated if this
    variable can be */
bool FR::LValueBool::isKnown() const{
    return (*rt->isSet)? true: false;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as bool)  */
void FR::LValueBool::set(const IObjectConstSP& object){
    if (CBool::TYPE->isInstance(object)){
        const CBool* theDb = STATIC_CONST_CAST(CBool, object.get());
        setValue(theDb->boolValue());
    } else {
        throw FRParseException("Cannot set bool from"
                               " type "+object->getClass()->getName());
    }
}

// get the variable expressed as a bool
bool FR::LValueBool::getValue() {
    if (!*rt->isSet){
        throw rt->makeException();
    }
    return rt->value;
}

// set the variable using a bool
void FR::LValueBool::setValue(bool value) {
    if (*rt->isSet){
        throw rt->makeException();
    }
    rt->value = value;
    *rt->isSet = 1;
}

//// uses FR::MemMgr
void* FR::LValueDate::RT::operator new(size_t size){
    return FR::MemMgr::alloc(size);
}
void FR::LValueDate::RT::operator delete(void *ptr){
    FR::MemMgr::dealloc(ptr);
}

FR::LValueDate::RT::RT(char* isSet, const char* name):
    func(&getValue), isSet(isSet), name(name){
    *isSet = 0;
}

/** generates exception saying either value not set yet or 
    already set */
ModelException FR::LValueDate::RT::makeException(){
    if (*isSet){
        string m("The value of the variable "+string(name)+
                 " has already been set");
        return ModelException("LValueDate::setValue", m);
    }
    string m("The value of the variable "+string(name)+
             " has not been set");
    return ModelException("LValueDate::getValue", m);
}

FR::LValueDate::LValueDate(const char* name): rt(new RT(&isSet, name)){}
        
FR::LValueDate::~LValueDate(){
    delete rt;
}
void FR::LValueDate::setReset(char* reset){
    rt->isSet = reset;
}

/** sets the value of the variable at this time point. More of
    indicative method as derived types should supply method allowing
    easier extraction of variable (eg as double)  */
void FR::LValueDate::set(const IObjectConstSP& object){
    const DateTime::Date& theVal =
        dynamic_cast<const DateTime::Date&>(*object);
    setValue(theVal);
}

// get the variable expressed as a date
const DateTime::Date& FR::LValueDate::getValue() {
    return RT::getValue(rt);
}
// set the variable using a double
void FR::LValueDate::setValue(const DateTime::Date& value) {
    RT::setValue(rt, value);
}
//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDate::RT* FR::LValueDate::getRT(){
    return (FRIfaces::IRValueDate::RT*) rt;
}


//// designed to ensure all memory accessed during sim run is from
//// MemMgr - so avoid vector etc here
struct FR::RDoubleArray::MyRT{
    TGetSize*                        sizeFunc;
    TGetValue*                       func;
    int                              size;
    FRIfaces::IRValueDouble::RT**    rValuesRT;
    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }
    ~MyRT(){
        FR::MemMgr::dealloc(rValuesRT);
    }

    //// simple constructor
    explicit MyRT(vector<FRIfaces::IRValueDouble*>& vectorRValues):
        sizeFunc(&getSize), func(&getValue), size(vectorRValues.size()),
        // allocate memory
        rValuesRT((FRIfaces::IRValueDouble::RT**)FR::MemMgr::
                  alloc(size * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < size; i++){
            FRIfaces::IRValueDouble* rValue = vectorRValues[i];
            rValuesRT[i] = rValue->getRT();
        }
    }

    /** For size field in RT structure */
    static int getSize(void* structure){
        MyRT* rt = (MyRT*)structure;
        return rt->size;
    }
    /** For func field in RT structure */
    static double getValue(void* structure, int index){
        MyRT* rt = (MyRT*)structure;
        FRIfaces::IRValueDouble::RT* dbRT = rt->rValuesRT[index];
        return (dbRT->func(dbRT));
    }
};

/** constructor for initially empty array */
FR::RDoubleArray::RDoubleArray(): rt(0), known(true){}

/** Add an FRIfaces::IRValueDouble to the object */
void FR::RDoubleArray::addArg(FRIfaces::IRValueDouble* rValue){
    vectorRValues.push_back(rValue);
    if (known && !rValue->isKnown()){
        known = false;
    }
}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueDoubleArray::RT* FR::RDoubleArray::getRT(){
    if (!rt){
        rt = new MyRT(vectorRValues);
    }
    return (FRIfaces::IRValueDoubleArray::RT*)rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant.  */
bool FR::RDoubleArray::isKnown() const{
    return known;
}
         
FR::RDoubleArray::~RDoubleArray(){
    delete rt;
}

//// designed to ensure all memory accessed during sim run is from
//// MemMgr - so avoid vector etc here
struct FR::RIntArray::MyRT{
    TGetSize*                        sizeFunc;
    TGetValue*                       func;
    int                              size;
    FRIfaces::IRValueInt::RT**       rValuesRT;
    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }
    ~MyRT(){
        FR::MemMgr::dealloc(rValuesRT);
    }
    //// simple constructor
    explicit MyRT(vector<FRIfaces::IRValueInt*>& vectorRValues):
        sizeFunc(&getSize), func(&getValue), size(vectorRValues.size()),
        // allocate memory
        rValuesRT((FRIfaces::IRValueInt::RT**)FR::MemMgr::
                  alloc(size * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < size; i++){
            FRIfaces::IRValueInt* rValue = vectorRValues[i];
            rValuesRT[i] = rValue->getRT();
        }
    }

    /** For size field in RT structure */
    static int getSize(void* structure){
        MyRT* rt = (MyRT*)structure;
        return rt->size;
    }
    /** For func field in RT structure */
    static int getValue(void* structure, int index){
        MyRT* rt = (MyRT*)structure;
        FRIfaces::IRValueInt::RT* dbRT = rt->rValuesRT[index];
        return (dbRT->func(dbRT));
    }
};

/** constructor for initially empty array */
FR::RIntArray::RIntArray(): rt(0), known(true){}

/** Add an FRIfaces::IRValueInt to the object */
void FR::RIntArray::addArg(FRIfaces::IRValueInt* rValue){
    vectorRValues.push_back(rValue);
    if (known && !rValue->isKnown()){
        known = false;
    }
}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueIntArray::RT* FR::RIntArray::getRT(){
    if (!rt){
        rt = new MyRT(vectorRValues);
    }
    return (FRIfaces::IRValueIntArray::RT*)rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant.  */
bool FR::RIntArray::isKnown() const{
    return known;
}
         
FR::RIntArray::~RIntArray(){
    delete rt;
}

//// designed to ensure all memory accessed during sim run is from
//// MemMgr - so avoid vector etc here
struct FR::RBoolArray::MyRT{
    TGetSize*                        sizeFunc;
    TGetValue*                       func;
    int                              size;
    FRIfaces::IRValueBool::RT**      rValuesRT;
    //// uses FR::MemMgr
    void* operator new(size_t size){
        return FR::MemMgr::alloc(size);
    }
    void operator delete(void *ptr){
        FR::MemMgr::dealloc(ptr);
    }
    ~MyRT(){
        FR::MemMgr::dealloc(rValuesRT);
    }
    //// simple constructor
    explicit MyRT(vector<FRIfaces::IRValueBool*>& vectorRValues):
        sizeFunc(&getSize), func(&getValue), size(vectorRValues.size()),
        // allocate memory
        rValuesRT((FRIfaces::IRValueBool::RT**)FR::MemMgr::
                  alloc(size * sizeof(void*))){
        // populate rt values
        for (int i = 0; i < size; i++){
            FRIfaces::IRValueBool* rValue = vectorRValues[i];
            rValuesRT[i] = rValue->getRT();
        }
    }

    /** For size field in RT structure */
    static int getSize(void* structure){
        MyRT* rt = (MyRT*)structure;
        return rt->size;
    }
    /** For func field in RT structure */
    static bool getValue(void* structure, int index){
        MyRT* rt = (MyRT*)structure;
        FRIfaces::IRValueBool::RT* dbRT = rt->rValuesRT[index];
        return (dbRT->func(dbRT));
    }
};

/** constructor for initially empty array */
FR::RBoolArray::RBoolArray(): rt(0), known(true){}

/** Add an FRIfaces::IRValueBool to the object */
void FR::RBoolArray::addArg(FRIfaces::IRValueBool* rValue){
    vectorRValues.push_back(rValue);
    if (known && !rValue->isKnown()){
        known = false;
    }
}

//// get the run-time object to use ie cut down version of whole class
FRIfaces::IRValueBoolArray::RT* FR::RBoolArray::getRT(){
    if (!rt){
        rt = new MyRT(vectorRValues);
    }
    return (FRIfaces::IRValueBoolArray::RT*)rt;
}

/** Is the value of this object known before the simulation starts
    eg a constant.  */
bool FR::RBoolArray::isKnown() const{
    return known;
}
         
FR::RBoolArray::~RBoolArray(){
    delete rt;
}

FR::LValueConst::~LValueConst(){}

/** allows default value of isKnown() to be overridden for testing
    purposes */
void FR::LValueConst::setIsValueKnown(bool isKnown){
    known = isKnown;
}

/** Is the value of this object known before the simulation starts
    eg a constant (This is useful as it allows some optimisations
    to be made) */
bool FR::LValueConst::isKnown() const{
    return known;
}
void FR::LValueConst::setReset(char* reset){}

//// throws exception
void FR::LValueConst::set(const IObjectConstSP& object){
    throw ModelException("FR::LValueConst::set", "Variable is read only");
}

/** Sets up intitial cache */
void FR::MemMgr::intialise(int numRules, int numDates){
    if (!heaps.empty()){
        throw ModelException("FR::MemMgr::intialise", "Already initialised");
    }
    heapSize = numRules * numDates * bytesPerRule;
    heapSize |= 8; // ensure it's a multiple of 8 (=size of double)
    heaps.push_back(Malloc::allocate(heapSize));
    pos = 0;
}

/** frees any memory */
void FR::MemMgr::shutdown(){
    for (unsigned int i = 0; i < heaps.size(); i++){
        Malloc::deallocate(heaps[i]);
    }
    heaps.resize(0);
    heapSize = 0;
}

/** allocates memory of requested size */
void* FR::MemMgr::alloc(size_t size){
    if (heapSize == 0){
        throw ModelException("FR::MemMgr::alloc", "Not intialised");
    }
    if (size > heapSize){
        throw ModelException("FR::MemMgr::alloc", "Request for "+
                             Format::toString((int)size)+". Too big");
    }
    if ((heapSize - pos) < size){ 
        heaps.push_back(Malloc::allocate(heapSize));
        pos = 0;
    }
    void* returnPtr = (char*)heaps.back() + pos;
#ifdef WIN32
    // only give out values on 4 byte boundaries
    if (size & (2+1)){
        size = size ^ (size & (2+1)); // remove last two bits
        size += 4;
    }
#else
    // only give out values on 8 byte boundaries
    if (size & (4+2+1)){
        size = size ^ (size & (4+2+1)); // remove last three bits
        size += 8;
    }
#endif
    pos += size;
    return returnPtr;
}

/** deallocates memory - currently does nothing */
void FR::MemMgr::dealloc(void* ptr){}

vector<void*> FR::MemMgr::heaps;    // the memory
size_t        FR::MemMgr::pos;      // where we are in heaps.back()
size_t        FR::MemMgr::heapSize; // size of heaps
const size_t  FR::MemMgr::bytesPerRule = 100;

class FR::Instrument: public GenericNFBase,
          virtual public IMCIntoProduct, 
          virtual public FlexBarrierBreach::IEventHandler {
private:
    friend class FlexInstrument;
    friend class Flex;
    friend class InputIMS;
    DateTimeArraySP                    simDates;
    FRIfaces::IAlgorithmSP             algorithm;
    FRIfaces::ILValueExpressionArraySP variables;
    FRIfaces::ILValueExpressionSP      payout;
    FRIfaces::DebugRequestSP           debugRequest;
    // cliquetDates is here as a work around for IMS where we build the 
    // vol request ourselves
    DateTimeArraySP                    cliquetDates; // transient
    bool                               debugOn;
    string                             flexCategory; // for audit purposes - not pricing
    bool                               triggerEvents; // $unregistered
public:
    static CClassConstSP const TYPE;

    virtual void Validate(){
        static const string method("FR::Instrument::Validate");
        GenericNFBase::Validate();
        if (simDates->empty()){
            throw ModelException(method, "Simulation dates empty");
        }
        if (variables->empty()){
            throw ModelException(method, "Variables empty");
        }
    }

    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        GenericNFBase::GetMarket(model, market);
        for(int i=0; i<variables->size(); i++) {
            IGetMarket* varGM = dynamic_cast<IGetMarket*>((*variables)[i].get());
            if (varGM) {
                varGM->getMarket(model, market.get());
            }
        }
    };

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        return *simDates;
    }

    /** Implementation of MonteCarlo::IntoProduct interface - the
        implementation of this is below */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const; // see below


    void writeRules(const string& fileName) const{
        algorithm->writeRules(fileName, payout->getID());
    }

    // Builds a FlexInstrument from this
    smartPtr<FlexInstrument> convert(DateTimeArraySP cliqDates);
       
    // implementation of FlexBarrierBreach::IEventHandler interface
    virtual void getEvents(const FlexBarrierBreach* breach,
                           IModel*   model, 
                           const DateTime& eventDate,
                           EventResults*   events) const;

 private:
    friend class ProductMC;
    friend class ProductMCSV;
    Instrument(): GenericNFBase(TYPE), debugOn(false), triggerEvents(false) {}
    Instrument(const Instrument& rhs); // not implemented
    Instrument& operator=(const Instrument& rhs); // not implemented

    static IObject* defaultInstrument(){
        return new Instrument();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Instrument, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        IMPLEMENTS(FlexBarrierBreach::IEventHandler);
        EMPTY_SHELL_METHOD(defaultInstrument);
        FIELD(simDates,             "Simulation Dates");
        FIELD(algorithm,            "How to compute payout");
        FIELD(variables,            "Variables used by algorithm");
        FIELD(payout,               "Variable detailing payments");
        FIELD_NO_DESC(cliquetDates);
        FIELD_MAKE_TRANSIENT(cliquetDates);
        FIELD(debugOn,       "true to switch on");
        FIELD_MAKE_OPTIONAL(debugOn);
        FIELD(debugRequest,         "Control of debug reporting");
        FIELD_MAKE_OPTIONAL(debugRequest);
        FIELD(flexCategory,  "Auditing category - not for pricing");
        FIELD_MAKE_OPTIONAL(flexCategory);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FLEX_SPI",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX instrument",
                                   TYPE);
        Addin::registerConstructor("FLEX",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX instrument",
                                   TYPE);
    }
};


CClassConstSP const FR::Instrument::TYPE = CClass::registerClassLoadMethod(
    "FR::Instrument", typeid(Instrument), Instrument::load);

//-----------------------------------------------------------------------

class FlexInstrument: public GenericNFBase,
                      virtual public IMCIntoProduct{
private:
    DateTimeArraySP      simDates;
    FlexAlgorithmSP      algorithm;
    FRVarFactorySP       variables;
    string               payoutName;
    bool                 debugOn; // eases debug though this is not seen in prod environment
    string               flexCategory; // for audit purposes - not pricing
    StringArray          assetVarNames; // will be matched up versus MultiAsset

    // Builds a FR::Instrument from this
    FR::InstrumentSP convert() {
        static const string method("FlexInstrument::convert");
        try {
            // Extend vars with ones for assets...
            FRIfaces::ILValueExpressionArraySP vars = variables->getAllVariables();
            // Note must clone vars because we are about to modify it
            vars = FRIfaces::ILValueExpressionArraySP(copy(vars.get()));
            for(int iAsset=0; iAsset<assetVarNames.size(); iAsset++) {
                // to do: add getTrueName to IMarketFactor/IMultiMarketFactors
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                if (!IGeneralAsset::TYPE->isInstance(factor)){
                    throw ModelException(method, "Only IGeneralAssets supported");
                }
                const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                          factor.get());
                FRIfaces::ILValueExpressionSP aVar =
                    FRVarFactory::createVariable("Asset",
                                                 assetVarNames[iAsset],
                                                 StringArray(1, asset->getTrueName()));
                vars->push_back(aVar);
            }

            // building a new FR::Instrument with same GenericNFBase ...
            CDataDictionarySP newInstDD(CDataDictionary::
                                        create(FR::Instrument::TYPE));
            try {
                // using 'this' to get GenericNFBase and lower elements
                IObjectSP myInst(IObjectSP(this));

                // populating them in the new instrument
                CClassConstSP c = GenericNFBase::TYPE;
                do {
                    const CFieldArray& fields = c->getDeclaredFields();
                    for (unsigned int i = 0; i < fields.size(); i++) {
                        if (!Modifier::isTransient(fields[i]->getModifiers())) {
                            IObjectSP obj = fields[i]->get(myInst);
                            if (!obj) {
                                obj = CNull::create();
                            }
                            newInstDD->put(fields[i]->getName(), obj);
                        }
                    }            
                } while ((c = c->getSuperClass()) != 0);
            }
            catch (exception &e) {
                throw ModelException (e, method, "Couldn't extract GenericNFBase "
                                      "data from FlexInstrument");
            }
            // add the FR::Instrument elements
            FRIfaces::IAlgorithmSP alg = algorithm->getAlgorithm(simDates);
            FRIfaces::ILValueExpressionSP payVar(variables->
                                                 getVariable(payoutName));

            newInstDD->put("simDates", simDates);
            newInstDD->put("algorithm", alg);
            newInstDD->put("variables", vars);
            newInstDD->put("payout", payVar);
            // and we're done
            FR::InstrumentSP FRinst(FR::InstrumentSP::dynamicCast(newInstDD->pop2Object()));
            if (!FRinst.get()) {
                throw ModelException(method, "Failed to build internal instrument!");
            }
            // well almost - copy over some optional params
            FRinst->flexCategory = flexCategory;
            FRinst->cliquetDates = cliquetDates;
            FRinst->debugOn = debugOn;

            return FRinst;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        if (inst.get()) {
            return inst->samplingDates();
        } 
        throw ModelException("FlexInstrument::samplingDates", "Internal error");
    }

protected:
    // protected in order to allow 
    // FlexCliquetInstrument::samplingDates() to see it
    FR::InstrumentSP     inst;

    DateTimeArraySP      cliquetDates; // optional
public:
    static CClassConstSP const TYPE;

    // Delegate
    virtual void Validate(){
        static const string routine("FlexInstrument::Validate");

        // At this point all manipulations of instrument should be
        // done (e.g. getMarket) so is a good point at which to build
        // our surrogate..  I originally had this in
        // validatePop2Object() but that was too "early" Translation
        // into the working instrument
        
        // XXX Could provide default - if assetVarNames empty or
        // strings empty could provide "a1", "a2", ...
        if (assets->NbAssets() != assetVarNames.size()) {
            throw ModelException(routine, 
                                 "#assetVarNames (" + 
                                 Format::toString(assetVarNames.size()) + 
                                 ") should be equal to number of assets (" +
                                 Format::toString(assets->NbAssets()) + ")");
        }

        // The inst used for pricing is built here since earlier (e.g. GetMarket)
        // introduces problems (which I've yet to fully understand).
        // The issue is that the variables need to be the ones which have had GetMarket
        // invoked. So they need to be preserved from the GetMarket.
        // Note this is an issue only for pricing - communication with IMS does
        // not care about market data.
        inst = convert();
        inst->Validate();

        GenericNFBase::Validate(); // fill out obsMap for parent
    }

    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        GenericNFBase::GetMarket(model, market);
        // This will be modifying the vars held inside the FRVarFactory.
        FRIfaces::ILValueExpressionArraySP vars = variables->getAllVariables();
        for(int i=0; i<vars->size(); i++) {
            IGetMarket* varGM = dynamic_cast<IGetMarket*>((*vars)[i].get());
            if (varGM) {
                varGM->getMarket(model, market.get());
            }
        }
    };


    /** Simply delegate */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const {
        return inst->createProduct(model);
    }
protected:
    FlexInstrument(CClassConstSP clazz): GenericNFBase(clazz), debugOn(false) {}
private:
    FlexInstrument(): GenericNFBase(TYPE), debugOn(false) {}
    FlexInstrument(const FlexInstrument& rhs); // not implemented
    FlexInstrument& operator=(const FlexInstrument& rhs); // not implemented

    friend class FR::IMSConvert;

    static IObject* defaultFlexInstrument(){
        return new FlexInstrument();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FlexInstrument, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultFlexInstrument);
        FIELD(simDates,             "Simulation Dates");
        FIELD(algorithm,            "How to compute payout");
        FIELD(variables,            "Variables used by algorithm");
        FIELD(payoutName,    "Variable name detailing payments");
        FIELD(debugOn,       "true to switch on");
        FIELD_MAKE_OPTIONAL(debugOn);
        FIELD(flexCategory,  "Auditing category - not for pricing");
        FIELD_MAKE_OPTIONAL(flexCategory);
        FIELD(assetVarNames, "Names of asset variables");
        FIELD(inst,                 "inst");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(inst);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        FIELD(cliquetDates, "Cliquet start dates");
        FIELD_MAKE_OPTIONAL(cliquetDates);
        Addin::registerConstructor("FLEX_INSTRUMENT",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Rules instrument "
                                   "(IMS style)",
                                   TYPE);
    }
};

typedef smartPtr<FlexInstrument> FlexInstrumentSP;

CClassConstSP const FlexInstrument::TYPE = CClass::registerClassLoadMethod(
    "FlexInstrument", typeid(FlexInstrument), FlexInstrument::load);

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Another class which is the same as FlexInstrument but has a different
// name. Allows us to book Flex cliquet under a different template in IMS
class FlexCliquetInstrument: public FlexInstrument,
                             virtual public IMCIntoProduct{
public:
    static CClassConstSP const TYPE;
private:
    FlexCliquetInstrument(): FlexInstrument(TYPE) {}

    static IObject* defaultConstructor(){
        return new FlexCliquetInstrument();
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        if (inst.get()) {
            return inst->samplingDates();
        } 
        throw ModelException("FlexCliquetInstrument::samplingDates", "Internal error");
    }

    virtual void validatePop2Object(){
        GenericNFBase::validatePop2Object();
        if (!cliquetDates || cliquetDates->empty()){
            throw ModelException("FlexCliquetInstrument::validatePop2Object",
                                 "No cliquet dates provided");
        }
    }
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FlexCliquetInstrument, clazz);
        SUPERCLASS(FlexInstrument);
        EMPTY_SHELL_METHOD(defaultConstructor);
    }
};

CClassConstSP const FlexCliquetInstrument::TYPE = 
CClass::registerClassLoadMethod(
    "FlexCliquetInstrument", typeid(FlexCliquetInstrument), load);

//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// Attempt at an interface that will fit nicely into IMS which supports abstraction etc

class Flex: public GenericNFBase,
            virtual public IMCIntoProduct{
private:
    FRIfaces::ILValueExpressionArraySP variables;
    string                             payoutName;
    FRIfaces::IAlgorithmSP             algorithm; // this one required to contain simDates
    FRIfaces::DebugRequestSP           debugRequest;
    bool                               debugOn; // eases debug though this is not seen in prod environment
    string                             flexCategory; // for audit purposes - not pricing
    StringArray                        assetVarNames; // will be matched up versus MultiAsset
    BoolArray                          assetIsAbsolute;   // matches 'assetVarNames' - indicates whether perf or not.

    // Builds a FR::Instrument from this
    FR::InstrumentSP convert() {
        static const string method("Flex::convert");
        try {
            // Extend vars with ones for assets...
            FRIfaces::ILValueExpressionArraySP vars(copy(variables.get()));
            for(int iAsset=0; iAsset<assetVarNames.size(); iAsset++) {
                // to do: add getTrueName to IMarketFactor/IMultiMarketFactors
                IMarketFactorConstSP factor(assets->getFactor(iAsset));
                if (!IGeneralAsset::TYPE->isInstance(factor)){
                    throw ModelException(method, "Only IGeneralAssets supported");
                }
                const IGeneralAsset* asset = DYNAMIC_CAST(IGeneralAsset,
                                                          factor.get());
                // bizarre, isn't it, that we need to go through such twists
                // just to provide a cleaner interface. Note that we validate
                // against the assetVarNames[] themselves containing any commas.
                // We could as an alternative publish a coding interface 
                // to FRAssetVariable, and call that directly.
                string varName = assetVarNames[iAsset];
                varName += assetIsAbsolute[iAsset]?",N":",Y";
                FRIfaces::ILValueExpressionSP aVar =
                    FRVarFactory::createVariable("Asset",
                                                 varName,
                                                 StringArray(1, asset->getTrueName()));
                vars->push_back(aVar);
            }

            // building a new FR::Instrument with same GenericNFBase ...
            CDataDictionarySP newInstDD(CDataDictionary::
                                        create(FR::Instrument::TYPE));
            try {
                // using 'this' to get GenericNFBase and lower elements
                IObjectSP myInst(IObjectSP(this));

                // populating them in the new instrument
                CClassConstSP c = GenericNFBase::TYPE;
                do {
                    const CFieldArray& fields = c->getDeclaredFields();
                    for (unsigned int i = 0; i < fields.size(); i++) {
                        if (!Modifier::isTransient(fields[i]->getModifiers())) {
                            IObjectSP obj = fields[i]->get(myInst);
                            if (!obj) {
                                obj = CNull::create();
                            }
                            newInstDD->put(fields[i]->getName(), obj);
                        }
                    }            
                } while ((c = c->getSuperClass()) != 0);
            }
            catch (exception &e) {
                throw ModelException (e, method, "Couldn't extract GenericNFBase "
                                      "data from Flex");
            }
            // add the FR::Instrument elements
            
            const DateTimeArraySP simDates = algorithm->getAllDates();
            if (simDates->size() < 1) {
                throw ModelException (method, "Internal error : Wrong algorithm type, "
                                      "so no sim dates identified.");
            }
            newInstDD->put("simDates", simDates);
            newInstDD->put("algorithm", algorithm);
            newInstDD->put("variables", vars);
            FRIfaces::ILValueExpressionSP payoutVariable;
            for(int j=0; j<vars->size(); j++) {
                if ((*vars)[j]->getID() == payoutName) {
                    payoutVariable = (*vars)[j];
                }
            }
            if (!payoutVariable) {
                throw ModelException(method, 
                                     "Supplied variables must contain payout variable named '" + 
                                     payoutName + "'");
            }
            newInstDD->put("payout", payoutVariable);
            // and we're done
            FR::InstrumentSP FRinst(FR::InstrumentSP::dynamicCast(newInstDD->pop2Object()));
            if (!FRinst.get()) {
                throw ModelException(method, "Failed to build internal instrument!");
            }
            // well almost - copy over some optional params
            FRinst->flexCategory = flexCategory;
            FRinst->cliquetDates = cliquetDates;
            FRinst->debugRequest = debugRequest;
            FRinst->debugOn = debugOn;

            return FRinst;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // returns the array of all sampling dates the instrument will ever need
    // excluding (possibly) the ref level dates (handled in GenericNFBase)
    const DateTimeArray samplingDates() const {
        if (inst.get()) {
            return inst->samplingDates();
        } 
        throw ModelException("Flex::samplingDates", "Internal error");
    }

protected:
    // protected in order to allow 
    // FlexCliquetInstrument::samplingDates() to see it
    FR::InstrumentSP     inst;

    DateTimeArraySP      cliquetDates; // optional
public:
    static CClassConstSP const TYPE;

    void validatePop2Object() {
        static const string routine("Flex::validatePop2Object");
        // forbid the "isPerf" being specified via the assetNames
        // instead user should use "assetIsAbsolute"
        for(int iAsset=0; iAsset<assetVarNames.size(); iAsset++) {
            size_t n = strcspn(assetVarNames[iAsset].c_str(), ",");
            size_t len = assetVarNames[iAsset].length();
            if (n < len) {
                // comma syntax forbidden via this interface
                throw ModelException(routine,
                                     "Comma syntax not supported here for asset[" +
                                     Format::toString(iAsset) +
                                     "] = " + assetVarNames[iAsset] +
                                     ". Instead use parameter assetIsAbsolute");
            }
        }
        if (assetVarNames.size() != assetIsAbsolute.size()) {
            throw ModelException(routine,
                                 "assetVarNames has " + Format::toString(assetVarNames.size()) +
                                 " entries, but assetIsAbsolute has " + Format::toString(assetIsAbsolute.size()) +
                                 ". These should be equal.");
        }
    }

    // Delegate
    virtual void Validate(){
        static const string routine("Flex::Validate");

        // At this point all manipulations of instrument should be
        // done (e.g. getMarket) so is a good point at which to build
        // our surrogate..  I originally had this in
        // validatePop2Object() but that was too "early" Translation
        // into the working instrument
        
        // XXX Could provide default - if assetVarNames empty or
        // strings empty could provide "a1", "a2", ...
        if (assets->NbAssets() != assetVarNames.size()) {
            throw ModelException(routine, 
                                 "#assetVarNames (" + 
                                 Format::toString(assetVarNames.size()) + 
                                 ") should be equal to number of assets (" +
                                 Format::toString(assets->NbAssets()) + ")");
        }

        // The inst used for pricing is built here since earlier (e.g. GetMarket)
        // introduces problems (which I've yet to fully understand).
        // The issue is that the variables need to be the ones which have had GetMarket
        // invoked. So they need to be preserved from the GetMarket.
        // Note this is an issue only for pricing - communication with IMS does
        // not care about market data.
        inst = convert();
        inst->Validate();

        GenericNFBase::Validate(); // fill out obsMap for parent
    }

    void GetMarket(const IModel*          model, 
                   const CMarketDataSP    market) {
        GenericNFBase::GetMarket(model, market);
        // This will be modifying the variables
        for(int i=0; i<variables->size(); i++) {
            IGetMarket* varGM = dynamic_cast<IGetMarket*>((*variables)[i].get());
            if (varGM) {
                varGM->getMarket(model, market.get());
            }
        }
    };


    /** Simply delegate */
    virtual IMCProduct* createProduct(
        const MonteCarlo* model) const {
        return inst->createProduct(model);
    }
protected:
    Flex(CClassConstSP clazz): GenericNFBase(clazz), debugOn(false) {}
private:
    Flex(): GenericNFBase(TYPE), debugOn(false) {}
    Flex(const Flex& rhs); // not implemented
    Flex& operator=(const Flex& rhs); // not implemented

    friend class FR::IMSConvert;

    static IObject* defaultFlex(){
        return new Flex();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Flex, clazz);
        SUPERCLASS(GenericNFBase);
        IMPLEMENTS(IMCIntoProduct);
        EMPTY_SHELL_METHOD(defaultFlex);
        FIELD(algorithm,            "How to compute payout");
        FIELD(variables,            "Variables used by algorithm");
        FIELD(payoutName,           "Name of variable detailing payments");
        FIELD(debugRequest,         "Control of debug reporting");
        FIELD_MAKE_OPTIONAL(debugRequest);
        FIELD(debugOn,       "true to switch on");
        FIELD_MAKE_OPTIONAL(debugOn);
        FIELD(flexCategory,  "Auditing category - not for pricing");
        FIELD_MAKE_OPTIONAL(flexCategory);
        FIELD(assetVarNames, "Names of asset variables");
        FIELD(assetIsAbsolute, "per asset - measure is absolute(true) or performance(false)");
        FIELD(inst,                 "inst");
        FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(inst);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        FIELD(cliquetDates, "Cliquet start dates");
        FIELD_MAKE_OPTIONAL(cliquetDates);
        Addin::registerConstructor("PRETTY_FLEX",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Rules instrument "
                                   "(IMS style)",
                                   TYPE);
    }
};

typedef smartPtr<Flex> FlexSP;

CClassConstSP const Flex::TYPE = CClass::registerClassLoadMethod(
    "Flex", typeid(Flex), Flex::load);

//-----------------------------------------------------------------------


class FR::WriteRules: public CObject{
public:
    static CClassConstSP const TYPE;

    // 2 addin parameters
    IObjectSP      instOrAlgorithm;
    string         fileName;

    /** set an object in a data dictionary */
    static IObjectSP writeRules(WriteRules* params){
        if (FR::Instrument::TYPE->isInstance(params->instOrAlgorithm)){
            FR::Instrument& inst =
                dynamic_cast<FR::Instrument&>(*params->instOrAlgorithm);
            inst.writeRules(params->fileName);
        } else if (FRIfaces::IAlgorithm::TYPE->
                   isInstance(params->instOrAlgorithm)){
            FRIfaces::IAlgorithm& algo = 
                dynamic_cast<FRIfaces::IAlgorithm&>(*params->instOrAlgorithm);
            algo.writeRules(params->fileName, "unspecified");
        }
        return IObjectSP(CString::create("OK"));
    }

    WriteRules(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(WriteRules, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultWriteRules);
        FIELD(instOrAlgorithm,   "FR instrument or algorithm");
        FIELD(fileName,    "To write to");
        Addin::registerInstanceObjectMethod(
            "FR_WRITE_RULES",
            Addin::FLEX_PAYOFF,
            "Prints the variable rules",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)writeRules);
    }

    static IObject* defaultWriteRules(){
        return new WriteRules();
    }
};

CClassConstSP const FR::WriteRules::TYPE = CClass::registerClassLoadMethod(
    "FR::WriteRules", typeid(FR::WriteRules), FR::WriteRules::load);

// * for class loading */
extern bool loadFRSimpleVariables();
extern bool loadFRMarketVariables();
extern bool loadFRBondVariable();
extern bool loadFRSimDateVariable();
extern bool loadFRDateVariable();
extern bool loadFRParser();
extern bool loadFRUtils();
extern bool loadFRSPI();
extern bool loadFRParserTPRules();
extern bool loadFRBarrierVariable();
extern bool loadFRLegalTermsVariables();
bool FRLoad() {
    return (Instrument::TYPE != 0 &&
            loadFRSimpleVariables() &&
            loadFRMarketVariables() &&
            loadFRBondVariable() &&
            loadFRSimDateVariable() &&
            loadFRDateVariable() &&
            loadFRParser() &&
            loadFRUtils() &&
            loadFRParserTPRules() &&
            loadFRSPI() &&
            loadFRBarrierVariable() &&
            loadFRLegalTermsVariables());
}

/* MC product class for Flex Instrument */
class FR::ProductMC: public IMCProduct,
          virtual public IMCProductLN,
          virtual public FRIfaces::IProductView{
private:
    const Instrument*           inst; // reference to original instrument
    mutable vector<CVolRequestLNArray>  volReqsPerAsset; // bit of a cheat
    refCountPtr<FRController>   frCtrl;
    SVGenSpotSP                    spotGen;     //!< Generator for spot
    IRefLevel::IStateVarGenSP   refLevelGen; //!< Generator for ref level
protected:
    //// MC calls this when path generator changes
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        frCtrl->pathGenUpdated(newPathGen); // delegate
    }
public:
    static const string FR_LAST_PATH;
    /** Provide IMCProduct with view of our product */
    ProductMC(const Instrument*         inst,
              bool                      useStateVars, // new style framework?
              const SimSeriesSP&        simSeries,
              bool                      triggerEvents):
        IMCProduct(inst->assets.get(),
                  inst->valueDate,
                  inst->discount.get(),
                  inst->refLevel,
                  simSeries,
                  inst->pastValues,
                  inst->instSettle.get(),
                  simSeries->getLastDate()),
        inst(inst), volReqsPerAsset(NbAssets),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())){
        frCtrl = refCountPtr<FRController>(new FRController(this, useStateVars,
                                                            true, triggerEvents));
    }

    virtual ~ProductMC() {}

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        frCtrl->collectStateVars(svCollector); // delegate
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form
        barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices) {
        if (inst->debugOn) {
            frCtrl->payoffWithDebug(pathGen, prices);
        } else {
            frCtrl->payoff(pathGen, prices);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string routine("FR::ProductMC::getVolInterp");
        if (iAsset >= (int) volReqsPerAsset.size()){
            throw ModelException(routine, "Asset index out of bounds");
        }
        // can put this block of code back in when we get rid of FRMaths
        if (volReqsPerAsset[iAsset].empty()){
            const string& name = getMultiFactors()->assetGetTrueName(iAsset);
            // helpful error message - may need to alter way vol interps are
            // collated since if an asset is not used you'll get this
            throw ModelException(routine, "No vol request info available for " 
                                 "asset '"+name+
                                 "'\nThis is because this asset is not being "
                                 "used in the rules");
        }
        return volReqsPerAsset[iAsset];
    }

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const{
        static const string method("FR::ProductMC:recordExtraOutput");
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
                frCtrl->recordPaymentDates(control, results);
        }
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished()) {
            frCtrl->recordKnownCashFlows(control, 
                                         results,
                                         discount->getCcy());
        }
        
        request = control->requestsOutput(OutputRequest::BARRIER_LEVEL);
        if (request && !request->getHasFinished()) {
            const DateTime& today = getToday();
            const DateTime& lastRefDate = getRefLevel()->getAllDates().back();
            // Only try to satisfy this request when the refLevel is 
            // complete. This makes safe an assumption to be found in
            // FRAssetVariable::getRefLevel(). Besides it's only sensible.
            if (today >= lastRefDate) {
                // sometimes can satisfy this one
                vector<string>              assetNames;
                vector<BarrierLevelArraySP> levels;
                frCtrl->getBarrierLevelReports(assetNames,
                                               levels);
                for(unsigned int iAsset=0; iAsset<assetNames.size(); iAsset++) {
                    OutputRequestUtil::recordBarrierLevels(
                        control, results,
                        assetNames[iAsset],
                        levels[iAsset].get());
                }
            }
        }            

        if (inst->debugOn) {
            // store final values of variables in debug packet
            IObjectSP  finalVals(frCtrl->getCurrentValues());
            OutputNameSP  outputName(new OutputName(FR_LAST_PATH));
            results->storeGreek(finalVals, Results::DEBUG_PACKET, outputName);

            if (inst->debugRequest.get()) {
                 // This is somewhat late but (as they say) better than never
                 // The results MAY be dodgy if greeks have been requested
                 // and the simulation split into blocks, so to be safe
                 // this simply forbids debug info with greeks. It would of course
                 // be better to check this up-front, or better still to quietly 
                 // switch off the blocking when doing debug, but that is not so
                 // simple to do.
                 SensitivityArrayConstSP sens(control->getSens());
                 if (sens->empty()){
                     if (inst->debugRequest->getFlagVar().get()) {
                         IObjectSP condVals(frCtrl->getDebugValues("Cond"));
                         if (condVals.get()) {
                             OutputNameSP  outputName(new OutputName("FR_COND_PATH"));
                             results->storeGreek(condVals, Results::DEBUG_PACKET, outputName);
                         }
                     }
                     if (inst->debugRequest->doMinValues()) {
                         IObjectSP minVals(frCtrl->getDebugValues("Min"));
                         if (minVals.get()) {
                             OutputNameSP  outputName(new OutputName("FR_MIN_PATH"));
                             results->storeGreek(minVals, Results::DEBUG_PACKET, outputName);
                         }
                     }
                     if (inst->debugRequest->doMaxValues()) {
                         IObjectSP maxVals(frCtrl->getDebugValues("Max"));
                         if (maxVals.get()) {
                             OutputNameSP  outputName(new OutputName("FR_MAX_PATH"));
                             results->storeGreek(maxVals, Results::DEBUG_PACKET, outputName);
                         }
                     }
                     if (inst->debugRequest->doAvgValues()) {
                         IObjectSP avgVals(frCtrl->getDebugValues("Avg"));
                         if (avgVals.get()) {
                             OutputNameSP  outputName(new OutputName("FR_AVG_PATH"));
                             results->storeGreek(avgVals, Results::DEBUG_PACKET, outputName);
                         }
                     }
                     if (inst->debugRequest->doDistns() &&
                         inst->debugRequest->getDistnVars().get()) {
                         IObjectSP distnVals(frCtrl->getDebugValues("Distn"));
                         if (distnVals.get()) {
                             OutputNameSP  outputName(new OutputName("FR_VAR_DISTN"));
                             results->storeGreek(distnVals, Results::DEBUG_PACKET, outputName);
                         }
                     }
                 } else {
                     // we're silently failing, so leave a clue as to why
                     // If we fail it causes hassle to client apps, and this is 
                     // hardly critical
                     IObjectSP warningMessage(CString::create("Flex Debug is only available when Price ONLY is requested.\n"
                                                              "Please turn off any sensitivity requests."));
                     if (inst->debugRequest->getFlagVar().get()) {
                         OutputNameSP  outputName(new OutputName("FR_COND_PATH"));
                         results->storeGreek(warningMessage, Results::DEBUG_PACKET, outputName);
                     }
                     if (inst->debugRequest->doMinValues()) {
                         OutputNameSP  outputName(new OutputName("FR_MIN_PATH"));
                         results->storeGreek(warningMessage, Results::DEBUG_PACKET, outputName);
                     }
                     if (inst->debugRequest->doMaxValues()) {
                         OutputNameSP  outputName(new OutputName("FR_MAX_PATH"));
                         results->storeGreek(warningMessage, Results::DEBUG_PACKET, outputName);
                     }
                     if (inst->debugRequest->doAvgValues()) {
                         OutputNameSP  outputName(new OutputName("FR_AVG_PATH"));
                         results->storeGreek(warningMessage, Results::DEBUG_PACKET, outputName);
                     }
                     if (inst->debugRequest->doDistns() &&
                         inst->debugRequest->getDistnVars().get()) {
                         OutputNameSP  outputName(new OutputName("FR_VAR_DISTN"));
                         results->storeGreek(warningMessage, Results::DEBUG_PACKET, outputName);
                     }
                 }


            }
        }
    }
    
    virtual void retrieveEvents(EventResults* events) const {
        frCtrl->retrieveEvents(events);
    }

    // IProductView implementation
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const {
        return inst->algorithm->createAssignments(frCtrl);
    }

    virtual const FRIfaces::ILValueExpression* getPayVariable(
        FRController* frCtrl) const {
        return inst->payout.get();
    }

    virtual const FRIfaces::DebugRequest* getDebugRequest(
        const FRController* frCtrl) const {
        if (inst->debugOn) {
            return inst->debugRequest.get();
        }
        return 0;
    }

    /** Returns the value date ie today */
    virtual const DateTime& getValueDate() const{
        return inst->valueDate;
    }

    /** returns the notional - used to scale values of each path by */
    virtual double getNotional() const{
        return inst->notional;
    }

    /** returns the payment date for each simulation date */
    virtual DateTimeArrayConstSP getPayDates() const{
        DateTimeArraySP payDates(new DateTimeArray(inst->simDates->size()));
        for (int i = 0; i < payDates->size(); i++){
            (*payDates)[i] = inst->instSettle->settles((*inst->simDates)[i],
                                                       0 /* asset */);
        }
        return payDates;
    }

    virtual double pvFromPaymentDate() const{
        return 1.0; // FR handles discounting
    }

    // no need to override endDate(const Sensitivity*  sensControl) method
    // as paymentDate date in IMCProduct is correct for this purpose

    /** populates array of discount factors. Each value being the 
        discount factor for a payment at the corresponding date. The return
        value indicates the value of the index which is the first to be 
        in the future */
    virtual int getDiscountFactors(
        const DateTimeArrayConstSP& payDates,
        DoubleArray&                discountFactors) const{
        int firstDateInFuture = 0;
        bool isMargin = inst->instSettle->isMargin();
        for (int i = 0; i < payDates->size(); i++){
            double df = 0.0;
            const DateTime& date = (*payDates)[i];
            if (date.isGreater(inst->valueDate)){
                df = isMargin? 1.0: inst->discount->pv(inst->valueDate, date);
            } else {
                firstDateInFuture++;
            }
            discountFactors[i] = df;
        }
        return firstDateInFuture;
    }
    
    virtual DateTimeArrayConstSP getSimDates() const {
        return inst->simDates;
    }

    /** Given asset's name and (optional) ccy treatment, returns the
        index for that asset */
    virtual int findAssetIndex(const string& assetName,
                               const string& ccyTreatment) const{
        static const string routine("FR::ProductMC::findAssetIndex");
        const IMultiFactors* mAsset = getMultiFactors();
        string               syntheticAssetName = assetName;
        if (!ccyTreatment.empty()){
            syntheticAssetName = 
                CAsset::getSyntheticName(assetName, ccyTreatment,
                                         inst->discount->getName());
        } 
        bool     found = false;
        int      assetIdx = 0; // avoid compiler warning
        for (int iAsset = 0; iAsset < NbAssets; ++iAsset) {
            const string& compareWith = ccyTreatment.empty()?
                assetName: syntheticAssetName;
            const string& compareTo =  ccyTreatment.empty()?
                mAsset->assetGetTrueName(iAsset): mAsset->assetGetName(iAsset);
            if (compareWith == compareTo){
                if (found) {
                    throw ModelException(
                        routine, 
                        "Asset "+ syntheticAssetName+ " requires ccyTreatment "
                        "qualifier since default is not well-defined.");
                }
                found = true;
                assetIdx = iAsset;
            }
        }
        if (!found) {
            throw ModelException(routine, "Asset "+ syntheticAssetName + 
                                 " not found in MultiAsset");
        }
        return assetIdx; // to do
    }
    
    /** Given asset index and (optional) vol request, returns the path
        index for this asset/vol request combination */
    virtual int findPathIndex(int                   assetIndex,
                              const CVolRequestSP&  volRequest) const{
        CVolRequestLNSP lnRequest;
        if (!volRequest){
            // In order to get Implied working (!!) need to have a LNStrike 
            // vol request : ATMVolRequest is not sufficient 
            // lnRequest = CVolRequestLNSP(new ATMVolRequest());
            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            // find something partway reasonable...
            const IMultiFactors* mf = getMultiFactors();
            double interpLevel = fwdStarting? 1.0: mf->assetGetSpot(assetIndex);
            if (!inst->cliquetDates){
                lnRequest = CVolRequestLNSP(
                    new LinearStrikeTSVolRequest(interpLevel,
                                                 startDate,
                                                 lastSimDate,
                                                 fwdStarting));
            } else {
                // need to build cliquet vol request
                if (inst->cliquetDates->empty()){
                    throw ModelException("FR::Instrument::findPathIndex",
                                         "No cliquet dates provided");
                }
                int cliqIdx = 0;
                while (cliqIdx < inst->cliquetDates->size() &&
                       !(*inst->cliquetDates)[cliqIdx].isGreater(today)){
                    cliqIdx++;
                }
                if (cliqIdx > 0){
                    // want live cliquet dates
                    cliqIdx--;
                }
                DateTimeArray liveCliqDates(
                    inst->cliquetDates->begin()+cliqIdx, 
                    inst->cliquetDates->end());
                DoubleArray interpLevels(liveCliqDates.size(), 1.0);
                interpLevels[0] = interpLevel;
                lnRequest = CVolRequestLNSP(
                    new CliquetVolRequest(liveCliqDates[0].isGreater(today),
                                          liveCliqDates,
                                          lastSimDate,
                                          interpLevels));
            }
        } else {
            lnRequest = CVolRequestLNSP::dynamicCast((IObjectSP)volRequest);
        }
        CVolRequestLNArray& requests = volReqsPerAsset[assetIndex];
        int pathIndex = requests.size();
        requests.push_back(lnRequest);
        return pathIndex;
    }

    virtual YieldCurveConstSP findYieldCurve(const string* ccyName) const {
        if (ccyName && ccyName->compare(discount->getName())) {
            throw ModelException("Supplied ccyName does not match "
                                 "discount curve of instrument.");
        }
        return inst->discount.getSP();
    }

    /* returns variables used in payoff */
    virtual const FRIfaces::ILValueExpressionArray* getVariables() const{
        return inst->variables.get();
    }

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        return inst->algorithm->getExpression(simDateIndex, date,
                                              assignmentIndex);
    }

    /** Returns a SV generator for instrument's reference level */
    virtual IRefLevel::IStateVarGenSP getRefLevelGen() const{
        return refLevelGen;
    }

    /** Returns a SV generator for instrument's assets' spots */
    virtual SVGenSpotSP getSpotGen() const{
        return spotGen;
    }

    /** Returns a SV generator for computing discount factors on the
        settlement adjusted dates supplied*/
    virtual SVGenDiscFactorSP getDiscFactorGen(
        const DateTimeArray& unadjustedDates) const{
        return SVGenDiscFactorSP(
            new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                             inst->instSettle, unadjustedDates));
    }
};
/* MC product class for Flex Instrument for SV. This is a cut and paste from
   above - which is really bad but currently you have to derive from different
   classes (to be fixed) */
class FR::ProductMCSV: public MCProductClient,
          virtual public IMCProductLN,
          virtual public FRIfaces::IProductView{
private:
    const Instrument*           inst; // reference to original instrument
    mutable vector<CVolRequestLNArray>  volReqsPerAsset; // bit of a cheat
    refCountPtr<FRController>   frCtrl;
    SVGenSpotSP                    spotGen;     //!< Generator for spot
    IRefLevel::IStateVarGenSP   refLevelGen; //!< Generator for ref level
protected:
    //// MC calls this when path generator changes
    virtual void pathGenUpdated(IStateVariableGen::IStateGen* newPathGen){
        frCtrl->pathGenUpdated(newPathGen); // delegate
    }
public:
    static const string FR_LAST_PATH;
    /** Provide IMCProduct with view of our product */
    ProductMCSV(const Instrument*         inst,
                bool                      useStateVars, // new style framework?
                const SimSeriesSP&        simSeries,
                bool                      triggerEvents):
        MCProductClient(inst->assets.get(),
                        inst->valueDate,
                        inst->discount.get(),
                        inst->refLevel,
                        simSeries,
                        inst->pastValues,
                        inst->instSettle.get(),
                        simSeries->getLastDate()),
        inst(inst), volReqsPerAsset(NbAssets),
        spotGen(new SVGenSpot(simSeries)),
        refLevelGen(inst->refLevel->createStateVarGen(getMultiFactors(), 
                                                      getToday())){
        frCtrl = refCountPtr<FRController>(new FRController(this, useStateVars,
                                                            true, triggerEvents));
    }

    virtual ~ProductMCSV() {}

    /** Appends 'true' (ie non derived) state variable generators
        required to the supplied collector.*/
    virtual void collectStateVars(IStateVariableCollectorSP svCollector) const{
        frCtrl->collectStateVars(svCollector); // delegate
    }

    /** Use this opportunity to do any LogNormal driven initialisation
        of the instrument before the main MC loop. e.g closed form
        barrier adjustment */
    void initialiseLN(const  IMCPathGenerator*  pathGen)const{
        // empty
    }

    /** Called within the simulation loop */
    void payoff(const MCPathGenerator*  pathGen,
                IMCPrices&                prices) {
        if (inst->debugOn) {
            frCtrl->payoffWithDebug(pathGen, prices);
        } else {
            frCtrl->payoff(pathGen, prices);
        }
    }

    // for the LogNormal path generator
    CVolRequestLNArray getVolInterp(const IMCPathGenerator* pathGen,
                                    int                     iAsset) const {
        static const string routine("FR::ProductMCSV::getVolInterp");
        if (iAsset >= (int) volReqsPerAsset.size()){
            throw ModelException(routine, "Asset index out of bounds");
        }
        // can put this block of code back in when we get rid of FRMaths
        if (volReqsPerAsset[iAsset].empty()){
            const string& name = getMultiFactors()->assetGetTrueName(iAsset);
            // helpful error message - may need to alter way vol interps are
            // collated since if an asset is not used you'll get this
            throw ModelException(routine, "No vol request info available for " 
                                 "asset '"+name+
                                 "'\nThis is because this asset is not being "
                                 "used in the rules");
        }
        return volReqsPerAsset[iAsset];
    }

    /** invoked after final simulated path is run. */
    virtual void recordExtraOutput(CControl*     control,
                                   Results*      results,
                                   const IMCPrices& prices) const{
        OutputRequest* request =
            control->requestsOutput(OutputRequest::PAYMENT_DATES);
        if (request && !request->getHasFinished()) {
            frCtrl->recordPaymentDates(control, results);
        }
        request =control->requestsOutput(OutputRequest::KNOWN_CASHFLOWS);
        if (request && !request->getHasFinished()) {
            frCtrl->recordKnownCashFlows(control, 
                                         results,
                                         discount->getCcy());
        }

        if (inst->debugOn) {
            // store final values of variables in debug packet
            IObjectSP  finalVals(frCtrl->getCurrentValues());
            OutputNameSP  outputName(new OutputName(FR_LAST_PATH));
            results->storeGreek(finalVals, Results::DEBUG_PACKET, outputName);

            IObjectSP condVals(frCtrl->getDebugValues("Cond"));
            if (condVals.get()) {
                OutputNameSP  outputName(new OutputName("FR_COND_PATH"));
                results->storeGreek(condVals, Results::DEBUG_PACKET, outputName);
            }
            IObjectSP minVals(frCtrl->getDebugValues("Min"));
            if (minVals.get()) {
                OutputNameSP  outputName(new OutputName("FR_MIN_PATH"));
                results->storeGreek(minVals, Results::DEBUG_PACKET, outputName);
            }
            IObjectSP maxVals(frCtrl->getDebugValues("Max"));
            if (maxVals.get()) {
                OutputNameSP  outputName(new OutputName("FR_MAX_PATH"));
                results->storeGreek(maxVals, Results::DEBUG_PACKET, outputName);
            }
            IObjectSP avgVals(frCtrl->getDebugValues("Avg"));
            if (avgVals.get()) {
                OutputNameSP  outputName(new OutputName("FR_AVG_PATH"));
                results->storeGreek(avgVals, Results::DEBUG_PACKET, outputName);
            }
            IObjectSP distnVals(frCtrl->getDebugValues("Distn"));
            if (distnVals.get()) {
                OutputNameSP  outputName(new OutputName("FR_VAR_DISTN"));
                results->storeGreek(distnVals, Results::DEBUG_PACKET, outputName);
            }
        }
    }

    virtual void retrieveEvents(EventResults* events) const {
        frCtrl->retrieveEvents(events);
    }

    // IProductView implementation
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const {
        return inst->algorithm->createAssignments(frCtrl);
    }

    virtual const FRIfaces::ILValueExpression* getPayVariable(
        FRController* frCtrl) const {
        return inst->payout.get();
    }

    virtual const FRIfaces::DebugRequest* getDebugRequest(
        const FRController* frCtrl) const {
        if (inst->debugOn) {
            return inst->debugRequest.get();
        }
        return 0;
    }

    /** Returns the value date ie today */
    virtual const DateTime& getValueDate() const{
        return inst->valueDate;
    }

    /** returns the notional - used to scale values of each path by */
    virtual double getNotional() const{
        return inst->notional;
    }

    /** returns the payment date for each simulation date */
    virtual DateTimeArrayConstSP getPayDates() const{
        DateTimeArraySP payDates(new DateTimeArray(inst->simDates->size()));
        for (int i = 0; i < payDates->size(); i++){
            (*payDates)[i] = inst->instSettle->settles((*inst->simDates)[i],
                                                       0 /* asset */);
        }
        return payDates;
    }

    virtual double pvFromPaymentDate() const{
        return 1.0; // FR handles discounting
    }

    // no need to override endDate(const Sensitivity*  sensControl) method
    // as paymentDate date in IMCProduct is correct for this purpose

    /** populates array of discount factors. Each value being the 
        discount factor for a payment at the corresponding date. The return
        value indicates the value of the index which is the first to be 
        in the future */
    virtual int getDiscountFactors(
        const DateTimeArrayConstSP& payDates,
        DoubleArray&                discountFactors) const{
        int firstDateInFuture = 0;
        bool isMargin = inst->instSettle->isMargin();
        for (int i = 0; i < payDates->size(); i++){
            double df = 0.0;
            const DateTime& date = (*payDates)[i];
            if (date.isGreater(inst->valueDate)){
                df = isMargin? 1.0: inst->discount->pv(inst->valueDate, date);
            } else {
                firstDateInFuture++;
            }
            discountFactors[i] = df;
        }
        return firstDateInFuture;
    }
    
    virtual DateTimeArrayConstSP getSimDates() const {
        return inst->simDates;
    }

    /** Given asset's name and (optional) ccy treatment, returns the
        index for that asset */
    virtual int findAssetIndex(const string& assetName,
                               const string& ccyTreatment) const{
        static const string routine("FR::ProductMCSV::findAssetIndex");
        const IMultiFactors* mAsset = getMultiFactors();
        string               syntheticAssetName = assetName;
        if (!ccyTreatment.empty()){
            syntheticAssetName = 
                CAsset::getSyntheticName(assetName, ccyTreatment,
                                         inst->discount->getName());
        } 
        bool     found = false;
        int      assetIdx = 0; // avoid compiler warning
        for (int iAsset = 0; iAsset < NbAssets; ++iAsset) {
            const string& compareWith = ccyTreatment.empty()?
                assetName: syntheticAssetName;
            const string& compareTo =  ccyTreatment.empty()?
                mAsset->assetGetTrueName(iAsset): mAsset->assetGetName(iAsset);
            if (compareWith == compareTo){
                if (found) {
                    throw ModelException(
                        routine, 
                        "Asset "+ syntheticAssetName+ " requires ccyTreatment "
                        "qualifier since default is not well-defined.");
                }
                found = true;
                assetIdx = iAsset;
            }
        }
        if (!found) {
            throw ModelException(routine, "Asset "+ syntheticAssetName + 
                                 " not found in MultiAsset");
        }
        return assetIdx; // to do
    }
    
    /** Given asset index and (optional) vol request, returns the path
        index for this asset/vol request combination */
    virtual int findPathIndex(int                   assetIndex,
                              const CVolRequestSP&  volRequest) const{
        CVolRequestLNSP lnRequest;
        if (!volRequest){
            // In order to get Implied working (!!) need to have a LNStrike 
            // vol request : ATMVolRequest is not sufficient 
            // lnRequest = CVolRequestLNSP(new ATMVolRequest());
            const DateTime&    today = getToday();
            const DateTime&    startDate = getRefLevel()->getAllDates().front();
            bool               fwdStarting = startDate.isGreater(today);
            const DateTime&    lastSimDate = getSimSeries()->getLastDate();
            // find something partway reasonable...
            const IMultiFactors* mf = getMultiFactors();
            double interpLevel = fwdStarting? 1.0: mf->assetGetSpot(assetIndex);
            if (!inst->cliquetDates){
                lnRequest = CVolRequestLNSP(
                    new LinearStrikeTSVolRequest(interpLevel,
                                                 startDate,
                                                 lastSimDate,
                                                 fwdStarting));
            } else {
                // need to build cliquet vol request
                if (inst->cliquetDates->empty()){
                    throw ModelException("FR::Instrument::findPathIndex",
                                         "No cliquet dates provided");
                }
                int cliqIdx = 0;
                while (cliqIdx < inst->cliquetDates->size() &&
                       !(*inst->cliquetDates)[cliqIdx].isGreater(today)){
                    cliqIdx++;
                }
                if (cliqIdx > 0){
                    // want live cliquet dates
                    cliqIdx--;
                }
                DateTimeArray liveCliqDates(
                    inst->cliquetDates->begin()+cliqIdx, 
                    inst->cliquetDates->end());
                DoubleArray interpLevels(liveCliqDates.size(), 1.0);
                interpLevels[0] = interpLevel;
                lnRequest = CVolRequestLNSP(
                    new CliquetVolRequest(liveCliqDates[0].isGreater(today),
                                          liveCliqDates,
                                          lastSimDate,
                                          interpLevels));
            }
        } else {
            lnRequest = CVolRequestLNSP::dynamicCast((IObjectSP)volRequest);
        }
        CVolRequestLNArray& requests = volReqsPerAsset[assetIndex];
        int pathIndex = requests.size();
        requests.push_back(lnRequest);
        return pathIndex;
    }

    virtual YieldCurveConstSP findYieldCurve(const string* ccyName) const {
        if (ccyName && ccyName->compare(discount->getName())) {
            throw ModelException("Supplied ccyName does not match "
                                 "discount curve of instrument.");
        }
        return inst->discount.getSP();
    }

    /* returns variables used in payoff */
    virtual const FRIfaces::ILValueExpressionArray* getVariables() const{
        return inst->variables.get();
    }

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        return inst->algorithm->getExpression(simDateIndex, date,
                                              assignmentIndex);
    }

    /** Returns a SV generator for instrument's reference level */
    virtual IRefLevel::IStateVarGenSP getRefLevelGen() const{
        return refLevelGen;
    }

    /** Returns a SV generator for instrument's assets' spots */
    virtual SVGenSpotSP getSpotGen() const{
        return spotGen;
    }

    /** Returns a SV generator for computing discount factors on the
        settlement adjusted dates supplied*/
    virtual SVGenDiscFactorSP getDiscFactorGen(
        const DateTimeArray& unadjustedDates) const{
        return SVGenDiscFactorSP(
            new SVGenDiscFactor(inst->valueDate,inst->discount.getSP(),
                             inst->instSettle, unadjustedDates));
    }
};

/** Implementation of MonteCarlo::IntoProduct interface */
IMCProduct* FR::Instrument::createProduct(const MonteCarlo* model) const {
    // we need to create a SimSeries object which says which assets need
    // which dates to be simulated
    SimSeriesSP simSeries(new SimSeries(assets->NbAssets())); /* create empty
                                                                 one */
    simSeries->addDates(*simDates);

    if (debugRequest.get()) {
        // This is convenient
        debugRequest->setDistnSize(model->getNbSubSamples());
    }

    if (model->stateVarUsed()){
        return new ProductMCSV(this, model->stateVarUsed(), simSeries, triggerEvents);
    } else {
        return new ProductMC(this, model->stateVarUsed(), simSeries, triggerEvents);
    }
}

void FR::Instrument::getEvents(const FlexBarrierBreach* breach,
                               IModel* model, 
                               const DateTime& eventDate,
                               EventResults* events) const {
    static const string method = "FR::Instrument::getEvents";
    try {
        // this is a bit weak ...
        MonteCarlo* mc = dynamic_cast<MonteCarlo*>(model);
        if (mc) {
            // should consider cloning model too
            InstrumentSP myCopy(copy(this));
            // set the flag so the FRController will handle events
            myCopy->triggerEvents = true;
            auto_ptr<IMCProduct> prod(myCopy->createProduct(mc));
            MCPathGeneratorSP past = prod->runPast(mc);
            prod->retrieveEvents(events);
        } else {
            throw ModelException("FR::Instrument::getEvents", 
                    "Internal error - expected Monte Carlo model for Flex pricing");
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

const string FR::ProductMC::FR_LAST_PATH = "FR_LAST_PATH";
const string FR::ProductMCSV::FR_LAST_PATH = "FR_LAST_PATH";

/** Class to allow non asset based operations to be tested easily without
    all the baggage needed by an instrument */
class FR::Tester: public CObject,
          virtual public FRIfaces::IProductView{
private:
    // fields
    DateTime                           valueDate;
    DateTimeArraySP                    simDates;
    FRIfaces::IAlgorithmSP             algorithm;
    FRIfaces::ILValueExpressionArraySP variables;
    FRIfaces::ILValueExpressionArraySP payout;
    /** following field is optional with default value false. It is used to
        in the setIsValueKnown() method of FRIfaces::ILConstValue */
    bool                               allowResolutionOfConstVars;
public:
    static CClassConstSP const TYPE;

    /** Returns the value date ie today */
    virtual const DateTime& getValueDate() const{
        return valueDate;
    }

    // create assignment array for current time point
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const{
        return algorithm->createAssignments(frCtrl);
    }

    virtual const FRIfaces::ILValueExpression* getPayVariable(
        FRController* frCtrl) const{
        return 0; // allows us to use any type of variable for payoff
    }

    virtual const FRIfaces::DebugRequest* getDebugRequest(
        const FRController* frCtrl) const {
        return 0;
    }

    /** returns the simulation dates */
    virtual DateTimeArrayConstSP getSimDates() const{
        return simDates;
    }
        
    /** returns the notional - used to scale values of each path by */
    virtual double getNotional() const{
        return 1.0; // hard coded 
    }

    /** returns the payment date for each simulation date */
    virtual DateTimeArrayConstSP getPayDates() const{
        return simDates;
    }

    /** populates array of discount factors. Each value being the 
        discount factor for a payment at the corresponding date. The return
        value indicates the value of the index which is the first to be 
        in the future */
    virtual int getDiscountFactors(
        const DateTimeArrayConstSP& payDates,
        DoubleArray&                discountFactors) const{
        for (int i = 0; i < discountFactors.size(); i++){
            discountFactors[i] = 1.0;
        }
        return 0;
    }
            
    /** Given asset's name and (optional) ccy treatment, returns the
        index for that asset */
    virtual int findAssetIndex(const string& assetName,
                               const string& ccyTreatment) const{
        throw ModelException("FR::Tester", "Asset information not available "
                             " in this simple tester");
    }

    /** Given asset index and (optional) vol request, returns the path
        index for this asset/vol request combination. To do: review this.
        Perhaps store vol requests in controller? */
    virtual int findPathIndex(int                   assetIndex,
                              const CVolRequestSP&  volRequest) const{
        throw ModelException("FR::Tester", "Asset information not available "
                             " in this simple tester");
    }

    // Given ccy name find the actual YieldCurve from instrument 
    virtual YieldCurveConstSP findYieldCurve(
        const string* ccyName) const {
        throw ModelException("FR::Tester", "Yield Curve information not "
                             "available in this simple tester");
    }

    /* returns variables used in payoff */
    virtual const FRIfaces::ILValueExpressionArray* getVariables() const{
        return variables.get();
    }

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        return algorithm->getExpression(simDateIndex, date, assignmentIndex);
    }

    /** Returns null - we haven't got one */
    virtual IRefLevel::IStateVarGenSP getRefLevelGen() const{
        return IRefLevel::IStateVarGenSP();
    }

    /**  Returns null as there are no assets */
    virtual SVGenSpotSP getSpotGen() const{
        return SVGenSpotSP();
    }

    /** Fails */
    virtual SVGenDiscFactorSP getDiscFactorGen(
        const DateTimeArray& unadjustedDates) const{
        throw ModelException("FR::Tester::getDiscFactorGen",
                             "Not supported");
    }

    virtual IObjectSP run(){
        FRController frCtrl(this, false, // should be irrelevant
                            allowResolutionOfConstVars, false);
        ObjectArraySP resultsArray(frCtrl.tester(payout.get()));
        if (resultsArray->size() == 1){
            return (*resultsArray)[0]; // for user's ease
        }
        return resultsArray;
    }

    static IObjectSP runAddin(Tester* params){
        return params->run();
    }

    Tester(): CObject(TYPE), allowResolutionOfConstVars(false) {}
    
    
    static IObject* defaultCreate(){
        return new Tester();
    }
        
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Tester, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreate);
        FIELD(valueDate,     "value Date");
        FIELD(simDates,             "simulation Dates");
        FIELD(algorithm,            "How to compute payout");
        FIELD(variables,            "Variables used by algorithm");
        FIELD(payout,               "Variable to return values for");
        FIELD(allowResolutionOfConstVars,
                     "Resolve const vars immediately");
        FIELD_MAKE_OPTIONAL(allowResolutionOfConstVars);
        Addin::registerClassObjectMethod("FR_TESTER",
                                         Addin::FLEX_PAYOFF,
                                         "Simplified FLEX tester",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)runAddin);
    }
};

CClassConstSP const FR::Tester::TYPE = CClass::registerClassLoadMethod(
    "FR::Tester", typeid(FR::Tester), load);


//-----------------------------------------------------------------------

class FR::InputIMS: public CObject{
public:
    static CClassConstSP const TYPE;

    // 2 addin parameters
    FR::InstrumentSP     instrument;

    /** set an object in a data dictionary */
    static IObjectSP createInput(InputIMS* params){
        int i;

        // the output object
        ObjectArraySP output = ObjectArraySP(new ObjectArray(0));
        StringArraySP titles(new StringArray(0));

        // simulation dates: we only keep the date removing the time of the day (e.g. EOD, SOD)
        DateTimeArraySP simDates = params->instrument->simDates;
        //ObjectArraySP simDatesIMS(new ObjectArray(simDates->size()));
        //for(i=0;i<simDates->size();i++) {
        //    (*simDatesIMS)[i] = IObjectSP(DateTimeSP((*simDates)[i]));
        //}        
        output->push_back(IObjectSP(simDates));
        titles->push_back("Simulation Dates");
        titles->push_back("Simulation Date Times");
                 
        // the variables : type, name and initial value
        FRIfaces::ILValueExpressionArraySP variables = params->instrument->variables;
        StringArraySP type(new StringArray(0));
        StringArraySP name(new StringArray(0));
        StringArraySP initialValue(new StringArray(0));
        StringArray      lookUpKeys;
        StringArrayArray lookUps;

        for(i=0;i<variables->size();i++) {
            StringArray IMSInput;
            IMSInput  = (*variables)[i]->getIMSInput();
            if (IMSInput.size() < 3) {
                throw ModelException("FR::InputIMS::createInput",
                                     "Internal error IMSInput.size()=" +
                                     Format::toString(IMSInput.size()) +
                                     " should be at least 3");
            }
            if (IMSInput.size() == 3) {
                // standard case
                if (IMSInput[0] != "Asset") {
                    type->push_back(IMSInput[0]);
                    name->push_back(IMSInput[1]);
                    initialValue->push_back(IMSInput[2]);
                }
            } else {
                // strings after the first 3 are a lookup
                type->push_back(IMSInput[0]);
                name->push_back(IMSInput[1]);

                // Delegate to FRFactory for interpretation of lookup
                StringArray IMSLookupInput(IMSInput.begin()+2, IMSInput.end());
                FRFactoryVarLookups::interpretIMSLookup(IMSLookupInput,
                                                        lookUpKeys,
                                                        lookUps); // modified with extra data
                initialValue->push_back(lookUpKeys.back());
            }
        }
        
        output->push_back(IObjectSP(type));
        titles->push_back("Variable Types");
        output->push_back(IObjectSP(name));
        titles->push_back("Variable Names");
        output->push_back(IObjectSP(initialValue));
        titles->push_back("Variable Values");

        // payoff variable
        StringArraySP payout(new StringArray(1));
        (*payout)[0]= params->instrument->payout->getID();
        output->push_back(IObjectSP(payout));
        titles->push_back("Payoff Variable Name");

        // the algorithm 
        FRIfaces::IAlgorithmSP algorithm = params->instrument->algorithm;
        //string algo;
        //algorithm.writesrules(algo);
        StringArrayArraySP result(new StringArrayArray(1));
        *result = algorithm->getIMSInput();
        output->push_back(IObjectSP((*result)[0]));
        titles->push_back("Simulation Date Rule ID");
        output->push_back(IObjectSP((*result)[1]));
        titles->push_back("Rule Definition ID");
        output->push_back(IObjectSP((*result)[2]));
        titles->push_back("Rule Definition Text");

        // cliquet dates
        DateTimeArraySP cliquetDates = params->instrument->cliquetDates;
        if (!cliquetDates) {
        }
        else {
            StringArraySP cliquetDatesIMS(new StringArray(cliquetDates->size()));
            for(i=0;i<cliquetDates->size();i++) {
                (*cliquetDatesIMS)[i+1] = DateTime::dateFormat((*cliquetDates)[i].getDate());
            }  
            output->push_back(IObjectSP(cliquetDatesIMS));
            titles->push_back("Cliquet Dates");
        } 

        // possible lookups - somewhat contrived
        // Issue being that it's not purely a function of the data; it also depends on
        // the current state of implementation of create() methods. 
        // For now, if we have a scalar with different values per sim date then 
        // implicitly IMS will need a "loopup". 
        // An empty lookup is ok.
        for(int j=0; j<lookUps.size(); j++) {
            output->push_back(IObjectSP(lookUps[j]));
            titles->push_back(lookUpKeys[j]);
        }
        // Helps with the spreadsheet. Better at the front of course...
        output->push_back(IObjectSP(titles));

        return output;
    }
    
    InputIMS(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(InputIMS, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreate);
        FIELD(instrument,   "FR instrument");
        Addin::registerInstanceObjectMethod(
            "FR_INPUT_IMS",
            Addin::FLEX_PAYOFF,
            "Create the input to create the corresponding instrument in IMS (Pyramid)",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)createInput);
    }

    static IObject* defaultCreate(){
        return new InputIMS();
    }
};

CClassConstSP const FR::InputIMS::TYPE = CClass::registerClassLoadMethod(
    "FR::InputIMS", typeid(FR::InputIMS), FR::InputIMS::load);

//-----------------------------------------------------------------------
// Maps FR::Instrument object to FlexInstrument/FlexCliquetInstrument object for interpretation by IMS
// Maps FlexInstrument object to FR::Instrument object for interpretation by Aladdin
class FR::IMSConvert: public CObject{
public:
    static CClassConstSP const TYPE;

    // FR::Instrument or FlexInstrument
    CInstrumentSP     inst;         // FR::Instrument or FlexInstrument
    CMarketDataSP     market;       // Associated market data
    DateTimeArraySP   cliqDates;    // Optional cliquet dates to build a FlexCliquetInstrument. A cumbersome alternative would be
                                    // to dig these out of CliquetVolRequests in the FRAssetVariables in the FR::Instrument
                                    
    /** set an object in a data dictionary */
    static IObjectSP convert(IMSConvert* params){
        static const string method = "FR::IMSConvert::convert";
        try {

            IObjectSP result;

            if (FlexInstrument::TYPE->isInstance(params->inst)) {
                // IMS (FlexInstrument or FlexCliquetInstrument) to Aladdin (FR::Instrument)
                // Need to make assets in FR::Instrument IGeneralAssets to get their names needed to build
                // AssetVariables
                MarketDataFetcherSP mdf(new MarketDataFetcherLN("VolPreferred"));
                NonPricingModel dummyModel(mdf);
                params->inst->GetMarket(&dummyModel, params->market);

                FlexInstrumentSP IMSimnt = FlexInstrumentSP(dynamic_cast<FlexInstrument *>(params->inst.get()));
                FR::InstrumentSP FRimnt = IMSimnt->convert();

                // But need to override algorithm with a FRAladdinAlgorithm
                //FRAladdinAlgorithmSP alg = IMSimnt->algorithm->getAlgorithm(IMSimnt->simDates);
                //alg->validatePop2Object();
                //FRimnt->algorithm = alg;

                result = IObjectSP(FRimnt.get());

            } else if (FR::Instrument::TYPE->isInstance(params->inst)) {
                // Aladdin (FR::Instrument) to IMS (FlexInstrument or FlexCliquetInstrument)
                FR::InstrumentSP FRimnt = FR::InstrumentSP(dynamic_cast<FR::Instrument *>(params->inst.get()));
                result = IObjectSP(FRimnt->convert(params->cliqDates).get());

            } else {

                throw ModelException(method, "Passed inst parameter must be of type FR::Instrument or FlexInstrument");
            }

            return result;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    IMSConvert(): CObject(TYPE), cliqDates(0) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IMSConvert, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCreate);
        FIELD(inst,   "Flex instrument");
        FIELD(market,   "Market data");
        FIELD(cliqDates, "Cliquet dates needed to build a FlexCliquetInstrument");
        FIELD_MAKE_OPTIONAL(cliqDates);
        Addin::registerInstanceObjectMethod(
            "FR_CONVERT",
            Addin::FLEX_PAYOFF,
            "Map between Aladdin and IMS Flex format",
            TYPE,
            false,
            Addin::returnHandle,
            (Addin::ObjMethod*)convert);
    }

    static IObject* defaultCreate(){
        return new IMSConvert();
    }
};

CClassConstSP const FR::IMSConvert::TYPE = CClass::registerClassLoadMethod(
    "FR::IMSConvert", typeid(FR::IMSConvert), FR::IMSConvert::load);

// Builds a FlexInstrument from this
FlexInstrumentSP FR::Instrument::convert(DateTimeArraySP cliqDates) {
    static const string method("FR::Instrument::convert");
    try {

        CDataDictionarySP newInstDD;
        CDataDictionarySP cliquetInstDD(CDataDictionary::create(FlexCliquetInstrument::TYPE));
        CDataDictionarySP InstDD(CDataDictionary::create(FlexInstrument::TYPE));

        if (!!cliqDates && cliqDates->size()>0) {
            newInstDD = cliquetInstDD;
        } else {
            // build a new FlexInstrument with same GenericNFBase ...
            newInstDD = InstDD;
        }

        // build a new FlexCliquetInstrument/FlexInstrument with same GenericNFBase ...
        try {
            // using 'this' to get GenericNFBase and lower elements
            IObjectSP myInst(IObjectSP(this));

            // populating them in the new instrument
            CClassConstSP c = GenericNFBase::TYPE;
            do {
                const CFieldArray& fields = c->getDeclaredFields();
                for (unsigned int i = 0; i < fields.size(); i++) {
                    if (!Modifier::isTransient(fields[i]->getModifiers())) {
                        IObjectSP obj = fields[i]->get(myInst);
                        if (!obj) {
                            obj = CNull::create();
                        }
                        newInstDD->put(fields[i]->getName(), obj);
                    }
                }            
            } while ((c = c->getSuperClass()) != 0);
        }
        catch (exception &e) {
            throw ModelException (e, "Couldn't extract GenericNFBase "
                                  "data from FR::Instrument");
        }

        // simDates
        newInstDD->put("simDates", simDates);

        // algorithm
        IntArraySP simDateRuleIds(new IntArray(0));
        IntArraySP ruleId(new IntArray(0));
        StringArraySP ruleDefn(new StringArray(0));
        algorithm->getIMSInput(simDateRuleIds, ruleId, ruleDefn);
        FlexAlgorithmSP imsAlg(new FlexAlgorithm(*simDateRuleIds, 
                                                 *ruleId, 
                                                 *ruleDefn));
        imsAlg->validatePop2Object();

        newInstDD->put("algorithm", imsAlg);

        // variables
        StringArraySP type(new StringArray(0));
        StringArraySP name(new StringArray(0));
        StringArraySP initialValue(new StringArray(0));
        // this pair of arrays maintain the "lookup" fudge getting round IMS limitations
        StringArray      lookUpKeys;
        StringArrayArray lookUps;

        StringArraySP assetVarNames(new StringArray(0));
        for(int i=0;i<variables->size();i++) {
            StringArray IMSInput = (*variables)[i]->getIMSInput();
            if (IMSInput.size() < 3) {
                throw ModelException(method,
                                     "Internal error IMSInput.size()=" +
                                     Format::toString(IMSInput.size()) +
                                     " but should be at least 3");
            }
            
            if (IMSInput.size() == 3) {
                // standard case
                if (IMSInput[0] != "Asset") {
                    type->push_back(IMSInput[0]);
                    name->push_back(IMSInput[1]);
                    initialValue->push_back(IMSInput[2]);
                } else {
                    assetVarNames->push_back(IMSInput[1]);
                }
            } else {
                // strings after the first 3 are a lookup
                type->push_back(IMSInput[0]);
                name->push_back(IMSInput[1]);

                // Delegate to FRFactory for interpretation of lookup
                StringArray IMSLookupInput(IMSInput.begin()+2, IMSInput.end());
                FRFactoryVarLookups::interpretIMSLookup(IMSLookupInput,
                                                        lookUpKeys,
                                                        lookUps); // modified with extra data
                initialValue->push_back(lookUpKeys.back());
            }
        }

        FRFactoryVarLookupsSP valueLookups(0);
        if (lookUps.size()>0) {
            valueLookups = FRFactoryVarLookupsSP(new FRFactoryVarLookups(lookUpKeys, lookUps));
        }
        // Builds internal variables field
        FRVarFactorySP vars(new FRVarFactory(*type,*name, *initialValue, valueLookups));

        newInstDD->put("variables", vars);

        // payoutName
        newInstDD->put("payoutName", payout->getID());

        // assetVarNames
        newInstDD->put("assetVarNames", assetVarNames);

        // copy over optional cliquet dates
        if (!!cliqDates && cliqDates->size()>0) {
            newInstDD->put("cliquetDates", cliqDates);
        }
        newInstDD->put("flexCategory", flexCategory);
        newInstDD->put("debugOn", debugOn);

        // we're done
        FlexInstrumentSP IMSinst(FlexInstrumentSP::dynamicCast(newInstDD->pop2Object()));

        if (!IMSinst.get()) {
            throw ModelException(method,
                                 "Failed to build internal instrument!");
        }

        return IMSinst;
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

// -----------------------------------------
// Debug features
void FRIfaces::DebugRequest::validatePop2Object() {
    static const string routine("FRIfaces::DebugRequest::validatePop2Object");

    int lenV = (distnVars.get() ? distnVars->size() : 0);
    int lenD = distnDates.size();
    if (lenV != lenD) {
        throw ModelException(routine,
                             "Distribution vars (#=" + Format::toString(lenV) +
                             ") and dates (#=" + Format::toString(lenD) + 
                             ") arrays must be of equal length");
    }
}

FRIfaces::ILValueExpressionSP FRIfaces::DebugRequest::getFlagVar() const {
    return captureDebug;
}

bool FRIfaces::DebugRequest::doMinValues() const {
    return captureMinValues;
}
bool FRIfaces::DebugRequest::doMaxValues() const {
    return captureMaxValues;
}
bool FRIfaces::DebugRequest::doAvgValues() const {
    return captureAvgValues;
}
bool FRIfaces::DebugRequest::doDistns() const {
    return captureDistn;
}
int FRIfaces::DebugRequest::getDistnSize() const {
    return distnSize;
}
void FRIfaces::DebugRequest::setDistnSize(int size) {
    distnSize = size;
    if (size<1) {
        throw ModelException("FRIfaces::DebugRequest::setDistnSize",
                             "Too small! Must be > 1");
    }
}
FRIfaces::ILValueExpressionArraySP FRIfaces::DebugRequest::getDistnVars() const {
    if (!captureDistn) {
        return FRIfaces::ILValueExpressionArraySP(0);
    }
    return distnVars;
}
const DateTimeArray& FRIfaces::DebugRequest::getDistnDates() const {
    return distnDates;
}

class DebugRequestHelper {
public:
    static IObject* defaultConstructor(){
        return new FRIfaces::DebugRequest();
    }

    static void load(CClassSP& clazz){
        REGISTER(FRIfaces::DebugRequest, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(captureDebug, "Variable flagging capture of debug info");
        FIELD_MAKE_OPTIONAL(captureDebug);
        FIELD(captureMinValues, "Debug to see min of vars?");
        FIELD_MAKE_OPTIONAL(captureMinValues);
        FIELD(captureMaxValues, "Debug to see max of vars?");
        FIELD_MAKE_OPTIONAL(captureMaxValues);
        FIELD(captureAvgValues, "Debug to see avg of vars?");
        FIELD_MAKE_OPTIONAL(captureAvgValues);
        FIELD(captureDistn, "Debug to see distn of selected vars?");
        FIELD_MAKE_OPTIONAL(captureDistn);
        FIELD(distnSize, "Size of distribution output array");
        FIELD_MAKE_TRANSIENT(distnSize);
        FIELD(distnVars, "Array of vars for which distribution desired");
        FIELD_MAKE_OPTIONAL(distnVars);
        FIELD(distnDates, "Date at which distribution is recorded (one per var)");
        FIELD_MAKE_OPTIONAL(distnDates);
        FIELD(captureAllValues, "Debug to see all vars at all iterations?");
        FIELD_MAKE_OPTIONAL(captureAllValues);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_DEBUG_REQUEST",
                                   Addin::FLEX_PAYOFF,
                                   "Defines debug features to be used",
                                   FRIfaces::DebugRequest::TYPE);
    }    
};
CClassConstSP const FRIfaces::DebugRequest::TYPE = CClass::registerClassLoadMethod(
    "FRIfaces::DebugRequest", typeid(DebugRequest), DebugRequestHelper::load);


DRLIB_END_NAMESPACE
