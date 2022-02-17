//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRParserTPRules.cpp
//
//   Description : Flex Rules: Set of rules for a given timepoint where rule
//                 is supplied as a single string
//
//   Author      : Mark A Robson
//
//   Date        : 1st August 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRParser.hpp"
#include "edginc/FRParseException.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/Hashtable.hpp" // for hash_string
#include "edginc/DECLARE.hpp"
#include <fstream>
#include <algorithm>
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE

struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};

/** Holds the set of rules expressing how each variable is to be
    calculated for a single time point */
class FRParserTPRules: public CObject, 
                       public virtual FRIfaces::ITimePointRules{
private:
    const StringArrayConstSP rules; // just array of string of form "a=b+2" etc
public:
    static CClassConstSP const TYPE;

    //// this object is immutable
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<FRParserTPRules*>(this);
    }

    virtual void validatePop2Object(){
        // validatePop2Object has privileged access
        StringArray& rules = const_cast<StringArray&>(*this->rules);
        // remove 'empty' rules
        for (vector<string>::iterator iter = rules.begin(); 
             iter != rules.end(); /* inc in loop body */){
            if (iter->empty() || *iter == ";"){
                iter = rules.erase(iter);
            } else{
                ++iter;
            }
        }
    }

    //// constructor - takes copy of supplied strings
    FRParserTPRules(const StringArray& expressions): 
        CObject(TYPE), rules(new StringArray(expressions)){}

    /** returns part after first '=' sign */
    static const char* getRValueExpression(const string& expression){
        // split up into lhs and rhs - just look for first equals sign
        const char* equalsSign = strchr(expression.c_str(), '=');
        if (!equalsSign){
            throw FRParseException::make("Parse error, no '=' sign found",
                                         expression, expression.size());
        }
        equalsSign++; // move onto next char
        if (*equalsSign == '='){
            throw FRParseException::make('=', expression, 
                                         equalsSign-expression.c_str());
        }
        return equalsSign;
    }
     
    /** returns part before first '=' sign - checks variable has no spaces */
    static string getLValueExpression(const string& expression,
                                      const char*   rValueExpression){
        if (rValueExpression == 0){
            rValueExpression = getRValueExpression(expression);
        }
        string lValueExpression(expression, 0, 
                                rValueExpression-1-expression.c_str());
        // trim white space - why isn't this in the stl?
        // first turn tabs into space
        replace(lValueExpression.begin(), lValueExpression.end(), '\t', ' ');
        size_t firstChar = lValueExpression.find_first_not_of(' ');
        if (firstChar == string::npos){
            // nothing before equals sign
            throw FRParseException::make('=', expression, 
                                         lValueExpression.size());
        }
        int lastChar = lValueExpression.find_last_not_of(' ');
        lValueExpression.assign(lValueExpression, firstChar,
                                lastChar-firstChar+1);
        size_t aSpace = lValueExpression.find(' ');
        if (aSpace != string::npos){
            // space in middle
            throw FRParseException::make(' ', expression, firstChar+aSpace);
        }
        return lValueExpression;
    }
        
    /** Creates an array of assignments - one for each rule at this
        timepoint. We should consider upgrading the parser to handle complete
        expressions */
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const{
        static const string routine("FRParser::TimePointRules::"
                                    "createAssignments");
        // allocate storage space
        auto_ptr<FRIfaces::IAssignmentArray> assignments(
            new FRIfaces::IAssignmentArray());
        assignments->reserve(rules->size());
        string  parseError; // to build up collection of parse errors
        int numErrors = 0;
        // loop through each rule
        try{
            for (int i = 0; i < rules->size(); i++){
                const string& expression = (*rules)[i];
                const char* rValueExpression; 
                string lValueExpression; 
                FRIfaces::ILValueExpression* lValueExp;
                FRIfaces::IRValue* rValue;
                try {
                    // split up string into two parts
                    rValueExpression = getRValueExpression(expression);
                    lValueExpression = 
                        getLValueExpression(expression, rValueExpression);
                    // reference variable name back to its declaration
                    lValueExp = frCtrl->getLValueExpression(lValueExpression);
                    if (!lValueExp){
                        throw FRParseException::make(
                            "Undefined variable '"+lValueExpression+"'",
                            expression, 0);
                    }
                } catch (FRParseException& e){
                    parseError += string(e.what())+"\n";
                    numErrors += e.getNumErrors();
                    if (numErrors > FRParseException::MAX_NUM_ERRORS){
                        throw FRParseException(numErrors,
                                               parseError+"Max number of "
                                               "errors exceeded");
                    }
                    continue; // seems easiest way
                }
                // can then get hold of lValue
                FRIfaces::ILValue* lValue = lValueExp->getLValue(frCtrl);
                // invoke parser to get rValue
                try{
                    rValue =
                        FRParser::parseRValue(rValueExpression, frCtrl);
                    // avoid a leak
                    frCtrl->store(rValue);
                    FRIfaces::IAssignment* assign =
                        FRController::createAssignment(lValue, rValue);
                    assignments-> push_back(FRIfaces::IAssignmentSP(assign));
                } catch (ModelException& e){
                    ModelException* cause = e.getCause();
                    FRParseException* parseEx = 
                        dynamic_cast<FRParseException*>(cause? cause:&e);
                    if (parseEx){
                        // record error and continue
                        parseError += string(e.what())+ "\n"+
                            "on assignment to '"+lValueExpression+ "' [" + 
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
            }
        } catch (exception& e){
            throw ModelException(e, routine);
        }
        if (!parseError.empty()){
            // throw one big one
            throw FRParseException(numErrors,
                                   "Encountered following parse errors:\n"+
                                   parseError);
        }
        return assignments.release();
    }
        
    /** Returns the number of rules for this given timepoint ie the
        number of assignments */
    virtual int numRules() const{
        return rules->size();
    }

    /** Returns the lValue for the specified rule as a string ie the
        variable name */
    virtual string lValueID(int ruleNum) const{
        checkIndex(ruleNum);
        return getLValueExpression((*rules)[ruleNum], 0);
    }

    /** Returns the rValue for the specified rule as a string eg X+Y */
    virtual string rValueID(int ruleNum) const{
        checkIndex(ruleNum);
        const char* rValueExpression = getRValueExpression((*rules)[ruleNum]);
        return string(rValueExpression);
    }

    //// Returns array of string representing expressions at each time point
    const StringArray& getExpressions() const{
        return *rules;
    }

    /** Finds index of rule for specified variable var. Returns rules.size()
        if not found */
    static int findRuleForVar(const StringArray& rules,
                              const string&      var){
        for (int i = 0; i < rules.size(); i++){
            const string& lValue = getLValueExpression(rules[i], 0);
            if (lValue == var){
                return i;
            }
        }
        return rules.size();
    }

    virtual ~FRParserTPRules(){}
private:
    void checkIndex(int ruleNum) const{
        if (ruleNum < 0 || ruleNum >= rules->size()){
            throw ModelException("FRParserTPRules::checkIndex", "Index out "
                                 "of bounds");
        }
    }

    FRParserTPRules(): CObject(TYPE) {};

    static IObject* defaultConstructor(){
        return new FRParserTPRules();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRParserTPRules, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::ITimePointRules);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(rules, "rules");
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FR_PARSER_TP_RULES",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI Parser TimePointRules",
                                   TYPE);
    }
};

CClassConstSP const FRParserTPRules::TYPE = CClass::registerClassLoadMethod(
    "FRParserTPRules", typeid(FRParserTPRules), load);

typedef smartPtr<FRParserTPRules> FRParserTPRulesSP;
typedef array<FRParserTPRulesSP, FRParserTPRules> FRParserTPRulesArray;
typedef smartPtr<FRParserTPRulesArray> FRParserTPRulesArraySP;
DEFINE_TEMPLATE_TYPE(FRParserTPRulesArray);

// Defined here for visibility of FRParserTPRules
void FlexAlgorithm::validatePop2Object() {
    static const string routine = "FlexAlgorithm::validatePop2Object";
    int i;
    int numRuleIDs = ruleId.size();
    if (numRuleIDs<1) {
        throw ModelException(routine, "No rule IDs!");
    }
    if (ruleDefn.size()!=numRuleIDs) {
        throw ModelException(routine, "Number of rule IDs (" +
                             Format::toString(numRuleIDs) + 
                             ") should equal number of rule definition lines ("+
                             Format::toString(ruleDefn.size()) + ")" );
    }
    // Check consistent IDs
    // So ruleId's should be an increasing set of integers, without jumps
    // And simDateRuleIds should not have values outside that range
    int minRuleId = 1;
    int maxRuleId = 1;
    for(i=0; i<numRuleIDs; i++) {
        if (ruleId[i]<minRuleId) {
            throw ModelException(routine,
                                 "Rule ID #" + Format::toString(i+1) +
                                 " must be at least 1");
        }
        if (ruleId[i]>maxRuleId) {
            if (ruleId[i]>maxRuleId+1) {
                throw ModelException(routine,
                                     "Rule ID #" + Format::toString(i+1) +
                                     " must be no bigger than " +
                                     Format::toString(maxRuleId+1));
            }
            maxRuleId++;
        }
    }
    for(i=0; i<simDateRuleIds.size(); i++) {
        if (simDateRuleIds[i]<minRuleId ||
            simDateRuleIds[i]>maxRuleId) {
            throw ModelException(routine,
                                 "Rule ID for sim date #"+Format::toString(i+1)+
                                 " is " + Format::toString(simDateRuleIds[i]) + 
                                 " but must be in range [" +
                                 Format::toString(minRuleId) +
                                 "," + Format::toString(maxRuleId) + "]");
        }
    }
    int numRules = maxRuleId;
    justTheRules = FRIfaces::ITimePointRulesArraySP(
        new FRIfaces::ITimePointRulesArray(numRules)); // convenient
    int r=0;
    for(i=0;i<numRules;i++) {
        // Collect rules by ID into a single array
        StringArray tpRulesText(0);
        for(;r<ruleId.size() && ruleId[r]==i+1;r++) {
            tpRulesText.push_back(ruleDefn[r]);
        }
        // Turn them into a TimePointRule
        (*justTheRules)[i] = FRIfaces::ITimePointRulesSP(
            new FRParserTPRules(tpRulesText));
    }
};



/** Implementation of FRIfaces::IAlgorithm where we specify an initial
    list of rules at the first timepoint and then make adjustments by
    adding to or removing from the rules at various dates */
class FRParserSparseAlgorithm: public CObject,
                               virtual public FRIfaces::IAlgorithm {
private:
    // fields
    const DateTimeArrayConstSP            dates; // when rules apply from
    const FRParserTPRulesArraySP          rules;
    const FRIfaces::ITimePointNoOpArraySP noOps; // what to skip
    const StringArrayConstSP              orderOfVars; // defines order
    const FRParserTPRulesArraySP          actRules; // not part of i/face
    // registration
    static void load(CClassSP& clazz){
        REGISTER(FRParserSparseAlgorithm, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::IAlgorithm);
        EMPTY_SHELL_METHOD(defaultFRParserSparseAlgorithm);
        FIELD(dates, "Defines which rules apply when");
        FIELD(rules,  "Incremental rules to apply at each "
              "corresponding date");
        FIELD(noOps,  "Incremental set of variables to skip");
        FIELD_MAKE_OPTIONAL(noOps);
        FIELD(orderOfVars, "Defines order in which to calculate variables");
        FIELD_MAKE_OPTIONAL(orderOfVars);
        FIELD(actRules, "Processed or actual set of rules");
        FIELD_MAKE_TRANSIENT(actRules);
        clazz->setPublic(); 
        Addin::registerConstructor("FR_SPARSE_ALGORITHM",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a sparse FLEX_SPI Algorithm",
                                   TYPE);
    }

    FRParserSparseAlgorithm(): CObject(TYPE){}

    static IObject* defaultFRParserSparseAlgorithm(){
        return new FRParserSparseAlgorithm();
    }

    /** this object is immutable */
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<FRParserSparseAlgorithm*>(this);
    }

    /** Given a date, find the index of the rule. Always returns a valid
        index */
    int findDateIndex(const DateTime& date) const{
        // need to find appropriate element in array - do binary search
        vector<DateTime>::const_iterator iter(
            upper_bound(dates->begin(), dates->end(), 
                        date, ComparatorLT()));
        if (iter == dates->begin()){
            throw ModelException("FRParserSparseAlgorithm::findDateIndex",
                                 "Simulation date "+date.toString()+ " is "
                                 "before first rule date "+ 
                                 (*dates)[0].toString());
        }
        int i = iter - dates->begin() -1; // calculate offset
        return i;
    }


public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method("FRParserSparseAlgorithm::"
                                   "validatePop2Object");
        int numNoOps = noOps.get()? noOps->size(): 0;
        if (dates->size() != rules->size() || 
            dates->size() != numNoOps){
            throw ModelException(method,
                                 Format::toString(dates->size())+" dates but "+
                                 Format::toString(rules->size())+
                                 " rules and "+
                                 Format::toString(numNoOps)+" no ops");
        }
        // ensure dates are not empty and are increasing
        DateTime::ensureIncreasing(*dates, "Dates defining when rules apply",
                                   true);
     
        // reserve space - validatePop2Object has special access
        // privileges
        const_cast<FRParserTPRulesArraySP&>(actRules) = 
            FRParserTPRulesArraySP(new FRParserTPRulesArray(dates->size()));

        // build up actual set of rules - note: have to preserve order which
        // prevents us using hash maps etc
        StringArray currentRules; // goes from timepoint to timepoint
        for (int i = 0; i < rules->size(); i++){
            /* add in list of rules from this timepoint (overwriting any
               existing ones). Then remove those in noops */
            if ((*rules)[i].get()){
                const FRParserTPRulesArray& theRules = *rules; // for ease
                const StringArray& expressions = theRules[i]->getExpressions();
                bool  orderChanged = false;
                if (currentRules.empty()){
                    currentRules = expressions;
                    orderChanged = true;
                } else {
                    // loop over expressions
                    for (int j = 0; j < expressions.size(); j++){
                        const string& lValue = theRules[i]->lValueID(j);
                        /* if it appears in the list, replace, else
                           append to end */
                        int index = FRParserTPRules::findRuleForVar(
                            currentRules, lValue);
                        if (index == currentRules.size()){
                            currentRules.push_back(expressions[j]);
                            orderChanged = true;
                        } else {
                            currentRules[index] = expressions[j];
                        }
                    }
                }
                if (orderChanged && orderOfVars.get()){
                    // ensure the order is as specified
                    StringArray tmpRules;
                    // loop over supplied variable names
                    for (int j = 0; j < orderOfVars->size(); j++){
                        // and locate corresponding rule
                        int index = FRParserTPRules::findRuleForVar(
                            currentRules, (*orderOfVars)[j]);
                        if (index != currentRules.size()){
                            // there is a rule for this variable
                            tmpRules.push_back(currentRules[index]);
                            currentRules.erase(currentRules.begin()+index); 
                        }
                    }
                    if (!currentRules.empty()){
                        const string& lValue = 
                            FRParserTPRules::getLValueExpression(
                                currentRules.front(), 0);
                        throw ModelException(method, "Variable "+ lValue+
                                             " is not in list of variables"
                                             " specifying rule order");
                    }
                    currentRules = tmpRules;
                }
            }
            // then remove no ops
            if (noOps.get() && (*noOps)[i].get()){
                const FRIfaces::ITimePointNoOpArray& theNoOps = *noOps;
                for (int j = 0; j < theNoOps[i]->numVars(); j++){
                    const string& lValue = theNoOps[i]->lValueID(j);
                    // find it
                    int index = FRParserTPRules::findRuleForVar(currentRules,
                                                                lValue);
                    if (index != currentRules.size()){
                        // remove if it exists
                        currentRules.erase(currentRules.begin()+index); 
                    }
                }
            }
            // finally record - validatePop2Object has special access
            // privileges
            (const_cast<FRParserTPRulesArray&>(*actRules))[i] = 
                FRParserTPRulesSP(new FRParserTPRules(currentRules));
        }
    }
    
    virtual ~FRParserSparseAlgorithm(){}
    
    struct ComparatorLT{
        bool operator()(const DateTime& d1, const DateTime& d2){
            return (d1.isLess(d2));
        }
    };
    
    /** create assignment array for current time point */
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const{
        const DateTime& currDate = frCtrl->getDate(frCtrl->getIndex());
        int i = findDateIndex(currDate);
        return ((*actRules)[i]->createAssignments(frCtrl));
    }
    
    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex,
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        int ourIndex = findDateIndex(date);
        const FRIfaces::ITimePointRules* ourRules = 
            (*actRules)[ourIndex].get();
        return (ourRules->lValueID(assignmentIndex) +" = "+
                ourRules->rValueID(assignmentIndex));
    }

    // May be better done on ProductMC ?
    virtual void writeRules(const string& fileName,
                            const string& payoutVariable) const {
        ofstream    ruleStream(fileName.c_str());
        // Loop over time points
        for (int t = 0; t < dates->size(); t++){
            ruleStream << "Time points >= " << (*dates)[t].toString() << endl;
            ruleStream << "-------------------" << endl;
            // loop over rules for timepoint
            if ((*rules)[t].get()){
                int numRules = (*rules)[t]->numRules();
                for (int i = 0; i < numRules; i++){
                    string lValueExp = (*rules)[t]->lValueID(i);
                    string rValueExp = (*rules)[t]->rValueID(i);
                    ruleStream << lValueExp << " = " << rValueExp << endl;
                }
            }
            if (noOps.get() && (*noOps)[t].get()){
                for (int j = 0; j < (*noOps)[t]->numVars(); j++){
                    const string& lValueID = (*noOps)[t]->lValueID(j);
                    ruleStream << lValueID << " is undefined" << endl;
                }
            }
            ruleStream << endl;
        }
        ruleStream << "\n===============\nPayout variable is " +
            payoutVariable + "\n";
#if 1
        // Loop over time points
        ruleStream << "\n------- DEBUG ---------------" << endl;
        for (int t1 = 0; t1 < dates->size(); t1++){
            ruleStream << "Time points >= " << (*dates)[t1].toString() << endl;
            ruleStream << "-------------------" << endl;
            // loop over rules for timepoint
            int numRules = (*actRules)[t1]->numRules();
            for (int i = 0; i < numRules; i++){
                string lValueExp = (*actRules)[t1]->lValueID(i);
                string rValueExp = (*actRules)[t1]->rValueID(i);
                ruleStream << lValueExp << " = " << rValueExp << endl;
            }
            ruleStream << endl;
        }
        ruleStream << "\n------- END OF DEBUG --------" << endl;
#endif
    }

    /** Report all dates known to the algorithm */
    DateTimeArraySP getAllDates() const {
        return DateTimeArraySP(0);
    }

    /** returns strings representing the input for IMS */
    virtual StringArrayArray getIMSInput() const {
        StringArraySP datesID(new StringArray(dates->size()+1));
        StringArraySP rulesID(new StringArray(1));
        StringArraySP rulesIMS(new StringArray(1));
        StringArrayArray result(3,StringArraySP(new StringArray(1,"Not Supported")));

        (*datesID)[0] = "Not Supported";
        (*rulesID)[0] = "Not Supported";
        (*rulesIMS)[0]   = "Not Supported";

        //result[0] = datesID;
        //result[1] = rulesID;
        //result[2] = rulesIMS;

        return result;
    }

    // do nothing
    virtual void getIMSInput(
        IntArraySP    simDateRuleIds,
        IntArraySP    ruleId,
        StringArraySP ruleDefn) const {}
};

CClassConstSP const FRParserSparseAlgorithm::TYPE = 
CClass::registerClassLoadMethod("FRParserSparseAlgorithm", 
                                typeid(FRParserSparseAlgorithm), load);

/** Implementation of FRIfaces::IAlgorithm. This class is designed to implement
    the IMSIput interface to generate the input for Pyramid
    Its purposed is to be used by Aladdin*/
class FRAladdinAlgorithm: public CObject,
          virtual public FRIfaces::IAlgorithm {
    
private:
    // fields
    const DateTimeArrayConstSP            dates; // when rules apply from
    const IntArrayConstSP                 rulesIndex; // index of the rules for each dates
    const FRParserTPRulesArraySP          rules; // set of rules
    const FRParserTPRulesArraySP          actRules; // not part of i/face
    
    // registration
    static void load(CClassSP& clazz){
        REGISTER(FRAladdinAlgorithm, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::IAlgorithm);
        EMPTY_SHELL_METHOD(defaultFRAladdinAlgorithm);
        FIELD(dates,      "Simulation dates");
        FIELD(rulesIndex, "Index of the rules corresponding to the simulation dates");
        FIELD(rules,      "Rules");
        FIELD(actRules, "Processed or actual set of rules");
        FIELD_MAKE_TRANSIENT(actRules);
        clazz->setPublic();
        Addin::registerConstructor("FR_ALADDIN_ALGORITHM",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX_SPI Algorithm, specific for Aladdin",
                                   TYPE);
    }
    
    FRAladdinAlgorithm(): CObject(TYPE){}

    static IObject* defaultFRAladdinAlgorithm(){
        return new FRAladdinAlgorithm();
    }

    /** this object is immutable */
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<FRAladdinAlgorithm*>(this);
    }
    
public:
    static CClassConstSP const TYPE;

    // Constructor needed for doing mapping back to Aladdin
//    FRAladdinAlgorithm(DateTimeArrayConstSP dates, 
//                       IntArrayConstSP rulesIndex,
//                       FRParserTPRulesArraySP rules) : CObject(TYPE), 
//                                                       dates(dates), 
//                                                       rulesIndex(rulesIndex),
//                                                       rules(rules) {}

    virtual void validatePop2Object(){
        static const string method("FRAladdinAlgorithm::"
                                   "validatePop2Object");
       
        // ensure there is one rule index for each date
        if (dates->size() != rulesIndex->size()) {
            throw ModelException(method,
                                 Format::toString(dates->size())+" dates but "+
                                 Format::toString(rulesIndex->size())+
                                 " rules index");
        }
        // ensure dates are not empty and are increasing
        DateTime::ensureIncreasing(*dates, "Dates defining when rules apply",
                                   true);
        
        // reserve space - validatePop2Object has special access
        // privileges
        const_cast<FRParserTPRulesArraySP&>(actRules) = 
            FRParserTPRulesArraySP(new FRParserTPRulesArray(dates->size()));
        
        // associate to each date the corresponding rule
        for (int i = 0; i < dates->size(); i++) {
            if ((*rulesIndex)[i] >= rules->size()) { //rulesIndex is 0-based...
                throw ModelException(method, "date "+
                                     Format::toString(i)+" has the set of rules of index "+
                                     Format::toString((*rulesIndex)[i])+
                                     " but the total number of rules is "+ 
                                     Format::toString(rules->size())+" (Rule indices begin at 0)");
            }
            else {
                // finally record - validatePop2Object has special access
                // privileges
                (const_cast<FRParserTPRulesArray&>(*actRules))[i] = 
                    FRParserTPRulesSP(new FRParserTPRules((*rules)[(*rulesIndex)[i]]->getExpressions()));
            }
        }
    }
    
    virtual ~FRAladdinAlgorithm(){}

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        const FRIfaces::ITimePointRules* ourRules = 
            (*actRules)[simDateIndex].get();
        return (ourRules->lValueID(assignmentIndex) +" = "+
                ourRules->rValueID(assignmentIndex));
    }

    /** create assignment array for current time point */
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const
    {
        int  i = frCtrl->getIndex();
        // detect whether there are more simulation dates than algorithm dates
        // specific to AladdinALgorithm
        if (i<dates->size()) {
            return ((*actRules)[i]->createAssignments(frCtrl));
        }
        else
            throw ModelException("FRAladdinAlgorithm::createAssignments",
                                 "all simulation dates must have been assigned a rule index");
    }

    // May be better done on ProductMC ?
    virtual void writeRules(const string& fileName,
                            const string& payoutVariable) const {
        ofstream    ruleStream(fileName.c_str());
        // Loop over time points
        for (int t = 0; t < dates->size(); t++){
            ruleStream << "Time points >= " << (*dates)[t].toString() << endl;
            ruleStream << "-------------------" << endl;
            // loop over rules for timepoint
            if ((*rules)[t].get()){
                int numRules = (*rules)[t]->numRules();
                for (int i = 0; i < numRules; i++){
                    string lValueExp = (*rules)[t]->lValueID(i);
                    string rValueExp = (*rules)[t]->rValueID(i);
                    ruleStream << lValueExp << " = " << rValueExp << endl;
                }
            }
            
            ruleStream << endl;
        }
        ruleStream << "\n===============\nPayout variable is " +
            payoutVariable + "\n";
#if 1
        // Loop over time points
        ruleStream << "\n------- DEBUG ---------------" << endl;
        for (int t1 = 0; t1 < dates->size(); t1++){
            ruleStream << "Time points >= " << (*dates)[t1].toString() << endl;
            ruleStream << "-------------------" << endl;
            // loop over rules for timepoint
            int numRules = (*actRules)[t1]->numRules();
            for (int i = 0; i < numRules; i++){
                string lValueExp = (*actRules)[t1]->lValueID(i);
                string rValueExp = (*actRules)[t1]->rValueID(i);
                ruleStream << lValueExp << " = " << rValueExp << endl;
            }
            ruleStream << endl;
        }
        ruleStream << "\n------- END OF DEBUG --------" << endl;
#endif
    }

    /** Report all dates known to the algorithm */
    DateTimeArraySP getAllDates() const {
        return DateTimeArraySP(0);
    }

    void getIMSInput(IntArraySP simDateRuleIds,        // O
                     IntArraySP ruleId,                // O
                     StringArraySP ruleDefn) const {   // O
        static const string method = "FRAladdinAlgorithm::getIMSInput";
        try {

            // Resize dates
            simDateRuleIds->resize(dates->size());

            // indices of rules must begin by 1 in IMS
            // create the array of simulation dates
            for(int iDates=0; iDates<dates->size(); iDates++) {
                (*simDateRuleIds)[iDates] = (*rulesIndex)[iDates] + 1;
            }

            // create the array of rules and rule indices
            int nbRules(0); // total number of rules
            for(int ID=0; ID<rules->size(); ID++) {
                ruleId->resize(ruleId->size()+(*rules)[ID]->numRules());
                ruleDefn->resize(ruleDefn->size()+(*rules)[ID]->numRules());
                for(int i=0; i<(*rules)[ID]->numRules(); i++) {
                    (*ruleId)[nbRules] = ID+1;
                    string lValueExp = (*rules)[ID]->lValueID(i);
                    string rValueExp = (*rules)[ID]->rValueID(i);
                    (*ruleDefn)[nbRules] = lValueExp+" = "+rValueExp;
                    nbRules++;
                }
            }

        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // create the equivalent input for IMS
    StringArrayArray getIMSInput() const {
        
        StringArraySP datesID(new StringArray(dates->size()));
        StringArraySP rulesID(new StringArray(0));
        StringArraySP rulesIMS(new StringArray(0));
        
        //(*datesID)[0] = "Simulation Date Rule ID";
        //(*rulesID)[0] = "Rule Definition ID";
        //(*rulesIMS)[0] = "Rule Definition Text";
        
        // indeces of rules must begin by 1 in IMS
        // create the array of simulation dates
        for(int iDates=0; iDates<dates->size(); iDates++) {
            (*datesID)[iDates] = Format::toString((*rulesIndex)[iDates]+1);
        }

        // create the array of rules and rule indices
        int nbRules(0); // total number of rules
        for(int ID=0; ID<rules->size(); ID++) {
            rulesID->resize(rulesID->size()+(*rules)[ID]->numRules());
            rulesIMS->resize(rulesIMS->size()+(*rules)[ID]->numRules());
            for(int i=0; i<(*rules)[ID]->numRules(); i++) {
                (*rulesID)[nbRules] = Format::toString(ID+1);
                string lValueExp = (*rules)[ID]->lValueID(i);
                string rValueExp = (*rules)[ID]->rValueID(i);
                (*rulesIMS)[nbRules] = lValueExp+" = "+rValueExp;
                nbRules++;
            }
        }

        StringArrayArray result(3);
        result[0] = datesID;
        result[1] = rulesID;
        result[2] = rulesIMS;

        return result;
    }
};

CClassConstSP const FRAladdinAlgorithm::TYPE = 
CClass::registerClassLoadMethod(
   "FRAladdinAlgorithm", typeid(FRAladdinAlgorithm), load);

/************************************************************************/

class FRRuleSet : public CObject {
private:
    // fields
    const string               ruleID;
    const FRParserTPRulesSP    rules; // set of rules

    // registration
    static void load(CClassSP& clazz){
        REGISTER(FRRuleSet, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFRRuleSet);
        FIELD(ruleID,       "ID for this rule set");
        FIELD(rules,  "Set of rules");
        clazz->setPublic();
        Addin::registerConstructor("FR_RULE_SET",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX RuleSet for use in PrettyAlgorithm",
                                   TYPE);
    }

    static IObject* defaultFRRuleSet(){
        return new FRRuleSet();
    }

public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method("FRRuleSet::validatePop2Object");
        if (ruleID.empty()) {
            throw ModelException(method,
                                 "Must supply a ruleID");
        }
    }

    // for reflection
    FRRuleSet(): CObject(TYPE){}

    const string& getID() const {
        return ruleID;
    }
    const FRParserTPRulesSP  getRules() {
        return rules;
    }

    virtual ~FRRuleSet(){}

}; 
DECLARE(FRRuleSet);

CClassConstSP const FRRuleSet::TYPE = 
CClass::registerClassLoadMethod(
   "FRRuleSet", typeid(FRRuleSet), load);

template <> CClassConstSP const FRRuleSetArray::TYPE = CClass::registerClassLoadMethod(
    "FRRuleSetArray", typeid(FRRuleSetArray), load);

/** Implementation of FRIfaces::IAlgorithm. This class is designed to 
    look pretty in IMS etc when driven by meta-data. */
class FRParserPrettyAlgorithm: public CObject,
                               virtual public FRIfaces::IAlgorithm {
private:
    // fields
    const DateTimeArraySP         simDates;   // when rule set applies
    const StringArray             ruleSetIDs; // [#simDates]
    const FRRuleSetArray          ruleSets; 

    // transient
    const FRParserTPRulesArraySP  actRules; 
    
    // registration
    static void load(CClassSP& clazz){
        REGISTER(FRParserPrettyAlgorithm, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(FRIfaces::IAlgorithm);
        EMPTY_SHELL_METHOD(defaultFRParserPrettyAlgorithm);
        FIELD(simDates,   "Simulation dates");
        FIELD(ruleSetIDs, "Identifiers of the rules corresponding to the simulation dates");
        FIELD(ruleSets,   "Definition of the rule sets");
        FIELD(actRules,   "Processed or actual set of rules");
        FIELD_MAKE_TRANSIENT(actRules);
        clazz->setPublic();
        Addin::registerConstructor("FR_PRETTY_ALGORITHM",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a FLEX Algorithm",
                                   TYPE);
    }
    
    FRParserPrettyAlgorithm(): CObject(TYPE){}

    static IObject* defaultFRParserPrettyAlgorithm(){
        return new FRParserPrettyAlgorithm();
    }

    /** this object is immutable */
    virtual IObject* clone() const{
        ensurePosRefCount();
        return const_cast<FRParserPrettyAlgorithm*>(this);
    }
    
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object(){
        static const string method("FRParserPrettyAlgorithm::"
                                   "validatePop2Object");
        // ensure there is one rule id for each date
        if (simDates->size() != ruleSetIDs.size()) {
            throw ModelException(method,
                                 Format::toString(simDates->size())+" dates but "+
                                 Format::toString(ruleSetIDs.size())+
                                 " rule IDs");
        }
        // ensure dates are not empty and are increasing
        DateTime::ensureIncreasing(*simDates, "SimDates",
                                   true);
        
        // reserve space - validatePop2Object has special access
        // privileges
        const_cast<FRParserTPRulesArraySP&>(actRules) = 
            FRParserTPRulesArraySP(new FRParserTPRulesArray(simDates->size()));
        
        // some helpful structures
        typedef hash_map<string, FRRuleSetSP,
            MyStringHash>        RuleSetHash;
        RuleSetHash ruleSetHash;
        for(int j=0; j<ruleSets.size(); j++) {
            const string& id = ruleSets[j]->getID();
            RuleSetHash::iterator iter = ruleSetHash.find(id);
            if (iter == ruleSetHash.end()) {
                // ok - add it
                ruleSetHash[id] = ruleSets[j];
            } else {
                throw ModelException(method, "Duplicate rule id :" +
                                     id + " in rule set #" + 
                                     Format::toString(j+1));
            }
        }

        // associate to each date the corresponding rule
        for (int i = 0; i < simDates->size(); i++) {
            RuleSetHash::iterator iter = ruleSetHash.find(ruleSetIDs[i]);
            if (iter == ruleSetHash.end()) {
                // not there - bad ID
                throw ModelException(method, "Unrecognised rule id [" +
                                     Format::toString(i) + "] = " +
                                     ruleSetIDs[i]);
            } else {
                // extract and use the rules
                // validatePop2Object has special access privileges
                (const_cast<FRParserTPRulesArray&>(*actRules))[i] = 
                    iter->second->getRules();
            }
        }
    }
    
    virtual ~FRParserPrettyAlgorithm(){}

    /** Return a string representing the rule for given simulation date
        and assignment eg a = MAX(b,c)+5. Used for error message */
    virtual string getExpression(int             simDateIndex, 
                                 const DateTime& date, // for simDateIndex
                                 int             assignmentIndex) const{
        const FRIfaces::ITimePointRules* ourRules = 
            (*actRules)[simDateIndex].get();
        return (ourRules->lValueID(assignmentIndex) +" = "+
                ourRules->rValueID(assignmentIndex));
    }

    /** create assignment array for current time point */
    virtual FRIfaces::IAssignmentArray* createAssignments(
        FRController* frCtrl) const
    {
        int  i = frCtrl->getIndex();
        // detect whether there are more simulation dates than algorithm dates
        // specific to ParserPrettyAlgorithm
        if (i<simDates->size()) {
            return ((*actRules)[i]->createAssignments(frCtrl));
        }
        else
            throw ModelException("FRParserPrettyAlgorithm::createAssignments",
                                 "all simulation dates must have been assigned a rule id");
    }

    // May be better done on ProductMC ?
    virtual void writeRules(const string& fileName,
                            const string& payoutVariable) const {
        ofstream    ruleStream(fileName.c_str());
        // Loop over time points
        for (int t = 0; t < simDates->size(); t++){
            ruleStream << "Time point " << (*simDates)[t].toString() << 
                " uses rule with ID " << ruleSetIDs[t] << endl;
            ruleStream << "-------------------" << endl;
        }
        for (int r=0; r<ruleSets.size(); r++) {
            ruleStream << "Rules for ID : " << ruleSets[r]->getID() << endl;
            ruleStream << "-------------------" << endl;
            FRParserTPRulesSP rules = ruleSets[r]->getRules();
            if (rules.get()){
                int numRules = rules->numRules();
                for (int i = 0; i < numRules; i++){
                    string lValueExp = rules->lValueID(i);
                    string rValueExp = rules->rValueID(i);
                    ruleStream << lValueExp << " = " << rValueExp << endl;
                }
            }
            ruleStream << endl;
        }
        ruleStream << "\n===============\nPayout variable is " +
            payoutVariable + "\n";
    }

    /** Report all dates known to the algorithm */
    DateTimeArraySP getAllDates() const {
        return simDates;
    }

    void getIMSInput(IntArraySP simDateRuleIds,        // O
                     IntArraySP ruleId,                // O
                     StringArraySP ruleDefn) const {   // O
        static const string method = "FRParserPrettyAlgorithm::getIMSInput";
        try {
            simDateRuleIds->resize(simDates->size());

            StringArray idMap;
            for(int i=0; i<ruleSets.size(); i++) {
                idMap.push_back(ruleSets[i]->getID());
            }
            for(int iDates=0; iDates<simDates->size(); iDates++) {
                // inefficient but doubt this will ever be used
                bool found = false;
                for(int j=0; j<idMap.size(); j++) {
                    if (ruleSetIDs[iDates]==idMap[j]) {
                        // based from 1, hence "+1"
                        (*simDateRuleIds)[iDates] = j + 1;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw ModelException(method, 
                                         "Unable to find index for simDate " +
                                         (*simDates)[iDates].toString());
                }
            }
            
            // create the array of rules and rule indices
            int nbRules(0); // total number of rules
            for (int s=0; s<ruleSets.size(); s++) {
                FRParserTPRulesSP rules = ruleSets[s]->getRules();
                ruleId->resize(ruleId->size()+rules->numRules());
                ruleDefn->resize(ruleDefn->size()+rules->numRules());
                // ruleSets[s] has index s+1 from above
                for(int j=0; j<rules->numRules(); j++) {
                    (*ruleId)[nbRules] = s+1;
                    string lValueExp = rules->lValueID(j);
                    string rValueExp = rules->rValueID(j);
                    (*ruleDefn)[nbRules] = lValueExp+" = "+rValueExp;
                    nbRules++;
                }
            }
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }

    // create the equivalent input for IMS
    StringArrayArray getIMSInput() const {
        IntArraySP simDateRuleIds(new IntArray(0));
        IntArraySP ruleId(new IntArray(0));
        StringArraySP ruleDefn(new StringArray(0));

        // use above method
        getIMSInput(simDateRuleIds,
                    ruleId,
                    ruleDefn);

        // translate into strings
        StringArraySP datesID(new StringArray(simDateRuleIds->size()));
        StringArraySP rulesID(new StringArray(ruleId->size()));

        for(int iDates=0; iDates<simDateRuleIds->size(); iDates++) {
            (*datesID)[iDates] = Format::toString((*simDateRuleIds)[iDates]);
        }
        for(int i=0; i<ruleId->size(); i++) {
            (*rulesID)[i] = Format::toString((*ruleId)[i]);
        }

        StringArrayArray result(3);
        result[0] = datesID;
        result[1] = rulesID;
        result[2] = ruleDefn;

        return result;
    }
};

CClassConstSP const FRParserPrettyAlgorithm::TYPE = 
CClass::registerClassLoadMethod(
   "FRParserPrettyAlgorithm", typeid(FRParserPrettyAlgorithm), load);

/************************************************************************/

bool loadFRParserTPRules() {
    return FRParserTPRules::TYPE != 0;
}

DRLIB_END_NAMESPACE
