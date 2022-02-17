//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRLegalTermsVariable.cpp
//
//   Description : Allows a contractual or risk variable to be used in rules
//                 according to whether pricing under LegalTerms
//
//   Date        : Aug 2005
//

//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRController.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/LegalTerms.hpp"

DRLIB_BEGIN_NAMESPACE

// Flex legal terms variable - common (type-independent) methods
class FRLegalTermsVariable: public FR::LValueBase,
                            public virtual LegalTerms::Shift  {
public:
    static CClassConstSP const TYPE;

    ~FRLegalTermsVariable(){}

    virtual void validatePop2Object() {
        static const string method = "FRLegalTermsVariable::validatePop2Object";
        // Check that types of risk/legal variables are consistent 

        // Variable cannot be self-referencing.
        // As for other variable types (e.g. DoubleArray) there is nothing preventing the user creating more 
        // complex circular dependencies involving intermediate variables. 
        if (name==legalVarName || name==riskVarName) {
            throw ModelException(method, "LegalTerms variable " + name +
                " cannot be defined in terms of itself");
        }
    }

    // Returns the type that this variable represents */
    // Ideally want to delegate this to the legal/risk variable as appropriate
    // but we can't get the variable without the FRController and getType is called
    // before CreateRValue()
    FRIfaces::VarType getType() const {
        throw ModelException("FRLegalTermsVariable::getType", 
                "Type of FRLegalTermsVariable is unknown");
    }

    /** Satisfy LegalTerms::Shift interface */
    virtual bool sensShift(LegalTerms* shift) {
        useLegalTerms = true;
        // May have other legal terms variables
        return true;
    }

    /** override RValueBase shallow copy default implementation */
    IObject* clone() const{
        return CObject::clone();
    }

    /** returns the id for this RValueBase - used by hashCode and
        equals() */
    const string& getID() const{
        return name;
    }

    virtual FRIfaces::ILValue* getLValue(FRController* frCtrl) const{
        throw ModelException("FRLegalTermsVariable::getLValue", "Variable "+
                             name +" is read-only");
    }

    FRLegalTermsVariable(CClassConstSP clazz): LValueBase(clazz), 
                                               useLegalTerms(false) 
    {
    // empty
    }

    FRLegalTermsVariable(CClassConstSP clazz, 
                         string varName, 
                         string legalVarName, 
                         string riskVarName): 
            LValueBase(clazz), name(varName), 
            legalVarName(legalVarName), riskVarName(riskVarName), 
            useLegalTerms(false) 
    {}


protected:

    // creates IRValue which represents the value
    // of this expression at a specific timepoint.
    virtual FRIfaces::IRValue* createRValue(int           index,
                                            FRController* frCtrl) const{
        static const string method = "FRLegalTermsVariable::createRValue";

        // Perform switch and return appropriate variable
        FRIfaces::IRValue* priceLevel = NULL;

        if (!frCtrl->getRValue(this, index)){
            // avoid putting values in more than once - this can happen if
            // this routine is called again after it has failed (eg missing
            // sample). A calling method might ignore the error in the [vain]
            // hope that the particular expression will not be evaluated at
            // run time (eg consider 'a = i>3? value[i-3]: 0.0')
            
            // Check both vars exist and are correct types
            FRIfaces::ILValueExpression* legalExp = frCtrl->getLValueExpression(legalVarName);
            FRIfaces::IRValueSP legalLevel;
            if (legalExp) { // a variable
                validateType(legalExp, index, frCtrl);
            } else {
                // not a variable - see if it can be parsed as a value
                legalLevel = FRIfaces::IRValueSP(parseValue(legalVarName)); // this is performing validation also
            }
            FRIfaces::ILValueExpression* riskExp = frCtrl->getLValueExpression(riskVarName);
            FRIfaces::IRValueSP riskLevel;
            if (riskExp) { // a variable
                validateType(riskExp, index, frCtrl);
            } else {
                // not a variable - see if it can be parsed as a value
                riskLevel = FRIfaces::IRValueSP(parseValue(riskVarName)); // this is performing validation also
            }

            // In the past we should always use the legal form
            const FRIfaces::IProductView* productView = frCtrl->productView();
            const DateTime&    today = productView->getValueDate();
            DateTimeArrayConstSP simDates = productView->getSimDates();
            // 'useLegalTerms' dominates. If !useLegalTerms but in the
            // past strictly before today we use legal terms anyway. 
            // Otherwise it's a bit special : if we're doing the "events" request then
            // today we use the 'risk', otherwise today counts as past and we use 'legal'
            bool myUseLegalTerms = useLegalTerms || ((*simDates)[index] < today) ||
                (!frCtrl->getTriggerEvents() && (*simDates)[index] == today);

            FRIfaces::ILValueExpression* priceExp = myUseLegalTerms?legalExp:riskExp;
            if (priceExp) {
                priceLevel = priceExp->getRValue(index, frCtrl);
            } else {
                priceLevel = myUseLegalTerms?legalLevel.release():riskLevel.release();
            }
        }

        return priceLevel;
    }

    // Validates type of var - default implementation does nothing
    // Type mismatches (between getType() and the value returned in CreateRValue()) may
    // be picked up in subsequent dynamic cast - prudent to validate explicitly
    virtual void validateType(FRIfaces::ILValueExpression* var,
                              int           index,
                              FRController* frCtrl) const
    {
    }

    // Attempts to parse a value of the appropriate type from the passed string.
    // Default implementation is to fail since the type of the value is not known.
    virtual FRIfaces::IRValue* parseValue(const string value) const {
        static const string method = "FRLegalTermsVariable::parseValue";

        string m("No variable with name " + value + "\n"
                    "When defining '" + name + "' you must supply the "
                 "names of risk and legal variables.");
        throw ModelException(method, m);
    }
        
    // Extract risk and legal var names from passed value
    static void getVarNames(const string& varName,
                             const string& value, 
                             string& legalVar, 
                             string& riskVar) 
    {
        static const string method = "FRLegalTermsVariable::getVarNames";
        const char* pos = value.c_str();
        int paramId = 0;
        while (*pos) {
            char   c;
            // skip whitespace
            while ((c = *pos) == ' ' || c == '\t'){
                pos++;
            }
            if (c == ','){
                // skip over our chosen separator and continue
                pos++;
                // count params
                paramId++;
                if (paramId>1) {
                    throw ModelException(method, 
                        "Too many commas (must have 1) with "
                        "legalterms variable " + varName + 
                        " and value of " + value);
                }
            } else {
                // next contiguous block is the text we want
                const char* varBegin = pos;
                do{
                    pos++; /* Get another character. */
                    c = *pos;
                } while (c != '\0' && c != ',' && c != ' ' && c != '\t');
                const char* varEnd = pos;
                // done with this one
                switch (paramId) {
                case 0:
                    legalVar = string(varBegin, varEnd - varBegin);
                    break;
                case 1:
                    riskVar = string(varBegin, varEnd - varBegin);
                    break;
                }
            } 
        }

        if (paramId != 1) {
            throw ModelException(method, "For legal terms variable " + varName + 
                                 " require name of legal variable (or value)"
                                 ", and name of risk variable (or value)");
        }
    }
    
    const string name;
    const string legalVarName;   // Legal variable name or value
    const string riskVarName;    // Risk variable name or value
    bool useLegalTerms;          // transient - set by LegalTerms scenario
   // mutable FRIfaces::ILValueExpressionSP relevantExp;  // for getting variable type - transient - set in createRValue()
    
private:


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRLegalTermsVariable, clazz);
        SUPERCLASS(FR::LValueBase);
        IMPLEMENTS(LegalTerms::Shift);
        FIELD(name, "name of variable");
        FIELD(legalVarName, "Legal variable name (or value)");
        FIELD(riskVarName, "Risk variable name (or value)");
        FIELD_NO_DESC(useLegalTerms);
        FIELD_MAKE_TRANSIENT(useLegalTerms);
        //FIELD_NO_DESC(relevantExp);
        //FIELD_MAKE_TRANSIENT(relevantExp);
        clazz->setPublic(); // make visible to EAS/spreadsheet
    }
    
};

// Flex legal terms double variable
// Any constituent variables must return a double
CClassConstSP const FRLegalTermsVariable::TYPE = CClass::registerClassLoadMethod(
    "FRLegalTermsVariable", typeid(FRLegalTermsVariable), FRLegalTermsVariable::load);

class FRLegalTermsDouble: public FRLegalTermsVariable {

public:
    static CClassConstSP const TYPE;

    ~FRLegalTermsDouble(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const {
        return FRIfaces::doubleType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "LegalTermsDouble";
        result[1] = name;
        result[2] = legalVarName +","+ riskVarName;

        return result;
    }

private:
    FRLegalTermsDouble(const string varName,
                       const string legalVarName,
                       const string riskVarName): 
            FRLegalTermsVariable(TYPE, varName, legalVarName, riskVarName) {}

    FRLegalTermsDouble(): FRLegalTermsVariable(TYPE)
    {}

    static IObject* defaultFRLegalTermsDouble(){
        return new FRLegalTermsDouble();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRLegalTermsDouble::create";

        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not appropriate here");
        }
        if (value[0].empty()) {
            throw ModelException(routine, "Value of " + varName + 
                        " should contain legalVarName, riskVarName "
                        "(variables or values)");
        }

        // Format of 'value' is these params in this order, comma separated
        string legalVarName; // Legal variable name or value
        string riskVarName;  // Risk variable name or value

        // Populate legalVarName and riskVarName from value
        getVarNames(varName, value[0], legalVarName, riskVarName);

        // A simple "return new ..." gives gcc compilation warning 
        FR::LValueBase* lterms = new FRLegalTermsDouble(varName,
                                                    legalVarName,
                                                    riskVarName);
        lterms->validatePop2Object();
        return lterms;
    }

    // Validates type of var is double
    virtual void validateType(FRIfaces::ILValueExpression* var,
                              int           index,
                              FRController* frCtrl) const
    {
        if (var->getType() != FRIfaces::doubleType) {
            throw ModelException("FRLegalTermsDouble::validateType", 
                            "Variables referenced in LegalTermsDouble " + name + " must be doubles");
        }
    }

    // Attempts to parse a double value from the passed string.
    virtual FRIfaces::IRValue* parseValue(const string value) const {
        static const string method = "FRLegalTermsDouble::parseValue";

        char*   endPos;
        double  dVal = strtod(value.c_str(), &endPos);
        if (*endPos == '\0') {

            return new FR::RConstDouble(dVal);

        } else {
            string m("No variable with name " + value + "\n"
                     "When defining '"+name+"' you must supply the "
                     "names of the risk/legal variables or explicit double values. ");
            throw ModelException(method, m);
        }
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FRLegalTermsDouble, clazz);
        SUPERCLASS(FRLegalTermsVariable);
        EMPTY_SHELL_METHOD(defaultFRLegalTermsDouble);

        Addin::registerConstructor("FR_VAR_LEGALTERMS_DOUBLE",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Legal Terms Double variable",
                                   TYPE);

        FRVarFactory::registerCreateMethod(FRLegalTermsDouble::create, "LegalTermsDouble");
    }

};

CClassConstSP const FRLegalTermsDouble::TYPE = CClass::registerClassLoadMethod(
    "FRLegalTermsDouble", typeid(FRLegalTermsDouble), FRLegalTermsDouble::load);

// Flex legal terms boolean variable
// Any constituent variables must return a bool (e.g. Barrier variables)
class FRLegalTermsBool: public FRLegalTermsVariable {

public:
    static CClassConstSP const TYPE;

    ~FRLegalTermsBool(){}

    /** Returns the type that this variable represents */
    FRIfaces::VarType getType() const {
        return FRIfaces::boolType;
    }

    /** returns the input to define the variable in IMS */
    virtual StringArray getIMSInput() const {
        StringArray result(3);
        
        result[0] = "LegalTermsBool";
        result[1] = name;
        result[2] = legalVarName +","+ riskVarName;

        return result;
    }

private:
    FRLegalTermsBool(const string varName,
                       const string legalVarName,
                       const string riskVarName): 
            FRLegalTermsVariable(TYPE, varName, legalVarName, riskVarName) {}

    FRLegalTermsBool(): FRLegalTermsVariable(TYPE)
    {}

    static IObject* defaultFRLegalTermsBool(){
        return new FRLegalTermsBool();
    }

    static FR::LValueBase* create(const string& varName, const StringArray& value) {
        static const string routine = "FRLegalTermsBool::create";

        if (value.size()!=1) {
            throw ModelException(routine, 
                                 "Per sim date initialisation not appropriate here");
        }
        if (value[0].empty()) {
            throw ModelException(routine, "Value of " + varName + 
                        " should contain legalVarName, riskVarName");
        }

        // Format of 'value' is these params in this order, comma separated
        string legalVarName; // Legal variable name or value
        string riskVarName;  // Risk variable name or value

        // Populate legalVarName and riskVarName from value
        getVarNames(varName, value[0], legalVarName, riskVarName);

        // A simple "return new ..." gives gcc compilation warning 
        FR::LValueBase* lterms = new FRLegalTermsBool(varName,
                                                    legalVarName,
                                                    riskVarName);
        lterms->validatePop2Object();
        return lterms;
    }

    // Validates type of var is bool
    virtual void validateType(FRIfaces::ILValueExpression* var,
                              int           index,
                              FRController* frCtrl) const
    {
        if (var->getType() != FRIfaces::boolType) {
            throw ModelException("FRLegalTermsBool::validateType", 
                "Variables referenced in LegalTermsBool " + name + " must be bools");
        }
    }

    // Attempts to parse a bool value from the passed string.
    virtual FRIfaces::IRValue* parseValue(const string value) const {
        static const string method = "FRLegalTermsBool::parseValue";

        bool    bVal;
        if (value == "true") {
            bVal = true;
        } else if (value == "false") {
            bVal = false;
        } else {
            string m("No variable with name " + value + "\n"
                     "When defining '"+name+"' you must supply the "
                     "names of the risk/legal variables or true/false.");
            throw ModelException(method, m);
        }

        return new FR::RConstBool(bVal);
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FRLegalTermsBool, clazz);
        SUPERCLASS(FRLegalTermsVariable);
        EMPTY_SHELL_METHOD(defaultFRLegalTermsBool);

        Addin::registerConstructor("FR_VAR_LEGALTERMS_BOOL",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Legal Terms boolean variable",
                                   TYPE);

        FRVarFactory::registerCreateMethod(FRLegalTermsBool::create, "LegalTermsBool");
    }

};

CClassConstSP const FRLegalTermsBool::TYPE = CClass::registerClassLoadMethod(
    "FRLegalTermsBool", typeid(FRLegalTermsBool), FRLegalTermsBool::load);


// this operates on a per-module basis
bool loadFRLegalTermsVariables() {
    return (FRLegalTermsDouble::TYPE != 0) &&
        (FRLegalTermsBool::TYPE !=0); 
}

DRLIB_END_NAMESPACE
