//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRFactory.cpp
//
//   Description : Central location to build FR variables & algorithms. Aimed at IMS interface
//
//   Date        : July 2003
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/FRFactory.hpp"
#include "edginc/Format.hpp"
#include "edginc/Hashtable.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE

struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};

// Store create methods keyed by type name
// XXX What about const?
typedef hash_map<string, void*, MyStringHash> VarCreateMethodLookup;

//**********************************************//
FRFactoryVarLookups::FRFactoryVarLookups() : CObject(TYPE) {
    lookUps.clear();
    lookUps.push_back(&perSimDate1);
    lookUps.push_back(&perSimDate2);
    lookUps.push_back(&perSimDate3);
    lookUps.push_back(&perSimDate4);
    lookUps.push_back(&perSimDate5);
    lookUps.push_back(&perSimDate6);
    lookUps.push_back(&perSimDate7);
    lookUps.push_back(&perSimDate8);
    lookUps.push_back(&perSimDate9);

    lookUps.push_back(&fixings1);
    lookUps.push_back(&fixings2);
    lookUps.push_back(&fixings3);
    lookUps.push_back(&fixings4);
    lookUps.push_back(&fixings5);
    lookUps.push_back(&fixings6);
    lookUps.push_back(&fixings7);
    lookUps.push_back(&fixings8);
    lookUps.push_back(&fixings9);
}

FRFactoryVarLookups::FRFactoryVarLookups(const StringArray&        lookUpKeys,
                                         const StringArrayArray&   lookUpData) : CObject(TYPE) {
    static const string method("FRFactoryVarLookups::FRFactoryVarLookups");

    lookUps.clear();
    lookUps.push_back(&perSimDate1);
    lookUps.push_back(&perSimDate2);
    lookUps.push_back(&perSimDate3);
    lookUps.push_back(&perSimDate4);
    lookUps.push_back(&perSimDate5);
    lookUps.push_back(&perSimDate6);
    lookUps.push_back(&perSimDate7);
    lookUps.push_back(&perSimDate8);
    lookUps.push_back(&perSimDate9);

    lookUps.push_back(&fixings1);
    lookUps.push_back(&fixings2);
    lookUps.push_back(&fixings3);
    lookUps.push_back(&fixings4);
    lookUps.push_back(&fixings5);
    lookUps.push_back(&fixings6);
    lookUps.push_back(&fixings7);
    lookUps.push_back(&fixings8);
    lookUps.push_back(&fixings9);

    if (lookUpKeys.size() != lookUpData.size()) {
        throw ModelException(method,
                             "Internal error. #keys=" + Format::toString(lookUpKeys.size()) +
                             " should equal #data arrays=" + Format::toString(lookUpData.size()));
    }
    for (int i=0; i<lookUpKeys.size() ; i++) {
        int lookUpIdx = parseLookUp(lookUpKeys[i]);
        if (lookUpIdx==0) {
            throw ModelException(method,
                                 "Unrecognised look-up key : " + lookUpKeys[i]);
        }
        if (lookUpIdx<0) {
            int idx = MAX_LOOKUPS - lookUpIdx - 1;
            // indicates fixing
            *(lookUps[idx]) = *(lookUpData[i].get());
            // maintain the CFA in step...
            switch (-lookUpIdx) {
            case 1:
                if (!fixings1.empty()) {
                    flexFixings1 = CashFlow::fromStringArray(fixings1);
                }
                break;
            case 2:
                if (!fixings2.empty()) {
                    flexFixings2 = CashFlow::fromStringArray(fixings2);
                }
                break;
            case 3:
                if (!fixings3.empty()) {
                    flexFixings3 = CashFlow::fromStringArray(fixings3);
                }
                break;
            case 4:
                if (!fixings4.empty()) {
                    flexFixings4 = CashFlow::fromStringArray(fixings4);
                }
                break;
            case 5:
                if (!fixings5.empty()) {
                    flexFixings5 = CashFlow::fromStringArray(fixings5);
                }
                break;
            case 6:
                if (!fixings6.empty()) {
                    flexFixings6 = CashFlow::fromStringArray(fixings6);
                }
                break;
            case 7:
                if (!fixings7.empty()) {
                    flexFixings7 = CashFlow::fromStringArray(fixings7);
                }
                break;
            case 8:
                if (!fixings8.empty()) {
                    flexFixings8 = CashFlow::fromStringArray(fixings8);
                }
                break;
            case 9:
                if (!fixings9.empty()) {
                    flexFixings9 = CashFlow::fromStringArray(fixings9);
                }
                break;
            default:
                throw ModelException(method, "Unsupported index : " + Format::toString(idx));
            }
        } else {
            // indicates per-sim-date
            *(lookUps[lookUpIdx-1]) = *(lookUpData[i].get());
        }
    }
}

void FRFactoryVarLookups::validatePop2Object() {
    static const string routine = "FRFactoryVarLookups::validatePop2Object";

    if (!flexFixings1.empty()) {
        fixings1 = CashFlow::toStringArray(flexFixings1);
    }
    if (!flexFixings2.empty()) {
        fixings2 = CashFlow::toStringArray(flexFixings2);
    }
    if (!flexFixings3.empty()) {
        fixings3 = CashFlow::toStringArray(flexFixings3);
    }
    if (!flexFixings4.empty()) {
        fixings4 = CashFlow::toStringArray(flexFixings4);
    }
    if (!flexFixings5.empty()) {
        fixings5 = CashFlow::toStringArray(flexFixings5);
    }
    if (!flexFixings6.empty()) {
        fixings6 = CashFlow::toStringArray(flexFixings6);
    }
    if (!flexFixings7.empty()) {
        fixings7 = CashFlow::toStringArray(flexFixings7);
    }
    if (!flexFixings8.empty()) {
        fixings8 = CashFlow::toStringArray(flexFixings8);
    }
    if (!flexFixings9.empty()) {
        fixings9 = CashFlow::toStringArray(flexFixings9);
    }
}

// Use negative numbers to indicate fixings :)
// XXX Extended to allow containment rather than just equality... ideally
// XXX though for the moment I'll just use the end of the string
int FRFactoryVarLookups::parseLookUp(const string& id) {
    char digit = '1';
    int i;
    unsigned int minSize = LOOKUP_PREFIX.size()+1;
    if (id.size()>=minSize) {
        for (i=1; i<=MAX_LOOKUPS ; i++) {
            string lastPart(id, id.size()-LOOKUP_PREFIX.size()-1, string::npos);
            if (CString::equalsIgnoreCase(LOOKUP_PREFIX+digit, lastPart)) {
                return i;
            }
            digit++;
        }
    }
    minSize = LOOKUP_FIXINGS_PREFIX.size()+1;
    if (id.size()>=minSize) {
        digit = '1';
        for (i=1; i<=MAX_LOOKUPS ; i++) {
            string lastPart(id, id.size()-LOOKUP_FIXINGS_PREFIX.size()-1, string::npos);
            if (CString::equalsIgnoreCase(LOOKUP_FIXINGS_PREFIX+digit, lastPart)) {
                return -i;
            }
            digit++;
        }
    }
    return 0; // special. Somehow I prefer this to throwing an exception
};

bool FRFactoryVarLookups::isPerSimDate(const string& id) {
    // for this purpose we allow no number
    if (id.size() == LOOKUP_PREFIX.size()) {
        if (CString::equalsIgnoreCase(LOOKUP_PREFIX, id)) {
            return true;
        }
    }
    char digit = '1';
    int i;
    unsigned int minSize = LOOKUP_PREFIX.size()+1;
    if (id.size()>=minSize) {
        for (i=1; i<=MAX_LOOKUPS ; i++) {
            string lastPart(id, id.size()-LOOKUP_PREFIX.size()-1, string::npos);
            if (CString::equalsIgnoreCase(LOOKUP_PREFIX+digit, lastPart)) {
                return true;
            }
            digit++;
        }
    }
    return false;
};

bool FRFactoryVarLookups::isLookUp(const string& id) {
    bool islookup = (parseLookUp(id) != 0);
    return islookup;
};

const StringArray& FRFactoryVarLookups::getLookUp(const string& id) const {
    int lookUpIdx = parseLookUp(id);
    if (lookUpIdx==0) {
        throw ModelException("FRFactoryVarLookups::getLookUp",
                             "This value : " + id +
                             " does not identify any data");
    }
    if (lookUpIdx<0) {
        // a fixings lookup, same array but with an offset
        int idx = MAX_LOOKUPS - lookUpIdx - 1;
        return *(lookUps[idx]);
    }

#ifdef DEBUG
    StringArray d0(*(lookUps[0]));
    StringArray d1(*(lookUps[1]));
    StringArray d2(*(lookUps[2]));
    StringArray d3(*(lookUps[3]));
    StringArray d4(*(lookUps[4]));
    StringArray d5(*(lookUps[5]));
    StringArray d6(*(lookUps[6]));
    StringArray d7(*(lookUps[7]));
    StringArray d8(*(lookUps[8]));
#endif
    return *(lookUps[lookUpIdx-1]);
};

    // Helper for IMS-friendly forms
void FRFactoryVarLookups::interpretIMSLookup(const StringArray&  IMSInput, // from the lookup key + remaining stuff
                                             StringArray&        lookUpKeys, // with specific 'id's attached
                                             StringArrayArray&   lookUpData) { // all the data so far
    static const string method("FRFactoryVarLookups::interpretIMSLookup");
    
    if (lookUpKeys.size() != lookUpData.size()) {
        throw ModelException(method,
                             "Internal error. #keys=" + Format::toString(lookUpKeys.size()) +
                             " should equal #data arrays=" + Format::toString(lookUpData.size()));
    }
    
    // Initialise local state from arrays collected so far (held by caller)
    char perSimDateId = '1';
    char fixingsId = '1';
    for(int i=0; i<lookUpKeys.size(); i++) {
        if (FRFactoryVarLookups::isPerSimDate(lookUpKeys[i])) {
            perSimDateId++;
        } else {
            fixingsId++;
        }
    }
    // If we've breached a limit we have to fail
    if (perSimDateId > '0'+ MAX_LOOKUPS) {
        throw ModelException(method,
                             "Too many perSimDate entries - maximum currently supported is " +
                             Format::toString(MAX_LOOKUPS));
    }
    if (fixingsId > '0'+ MAX_LOOKUPS) {
        throw ModelException(method,
                             "Too many fixings entries - maximum currently supported is " +
                             Format::toString(MAX_LOOKUPS));
    }
    
    // Identify IMSInput - "perSimDate" or "fixings?"
    string thisKey = IMSInput[0];
    if (FRFactoryVarLookups::isPerSimDate(thisKey)) {
        lookUpKeys.push_back(thisKey + perSimDateId);
    } else {
        lookUpKeys.push_back(thisKey + fixingsId);
    }
    StringArraySP thisLookUp(new StringArray(IMSInput.begin()+1, IMSInput.end()));
    lookUpData.push_back(thisLookUp);
};

IObject* FRFactoryVarLookups::defaultFRFactoryVarLookups(){
    return new FRFactoryVarLookups();
};

void FRFactoryVarLookups::load(CClassSP& clazz){
    REGISTER(FRFactoryVarLookups, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultFRFactoryVarLookups);
    FIELD(perSimDate1,     "Values for variable keyed as perSimDate1");
    FIELD_MAKE_OPTIONAL(perSimDate1);
    FIELD(perSimDate2,     "Values for variable keyed as perSimDate2");
    FIELD_MAKE_OPTIONAL(perSimDate2);
    FIELD(perSimDate3,     "Values for variable keyed as perSimDate3");
    FIELD_MAKE_OPTIONAL(perSimDate3);
    FIELD(perSimDate4,     "Values for variable keyed as perSimDate4");
    FIELD_MAKE_OPTIONAL(perSimDate4);
    FIELD(perSimDate5,     "Values for variable keyed as perSimDate5");
    FIELD_MAKE_OPTIONAL(perSimDate5);
    FIELD(perSimDate6,     "Values for variable keyed as perSimDate6");
    FIELD_MAKE_OPTIONAL(perSimDate6);
    FIELD(perSimDate7,     "Values for variable keyed as perSimDate7");
    FIELD_MAKE_OPTIONAL(perSimDate7);
    FIELD(perSimDate8,     "Values for variable keyed as perSimDate8");
    FIELD_MAKE_OPTIONAL(perSimDate8);
    FIELD(perSimDate9,     "Values for variable keyed as perSimDate9");
    FIELD_MAKE_OPTIONAL(perSimDate9);
    FIELD(flexFixings1,     "Values for variable keyed as flexFixings1");
    FIELD_MAKE_OPTIONAL(flexFixings1);
    FIELD(flexFixings2,     "Values for variable keyed as flexFixings2");
    FIELD_MAKE_OPTIONAL(flexFixings2);
    FIELD(flexFixings3,     "Values for variable keyed as flexFixings3");
    FIELD_MAKE_OPTIONAL(flexFixings3);
    FIELD(flexFixings4,     "Values for variable keyed as flexFixings4");
    FIELD_MAKE_OPTIONAL(flexFixings4);
    FIELD(flexFixings5,     "Values for variable keyed as flexFixings5");
    FIELD_MAKE_OPTIONAL(flexFixings5);
    FIELD(flexFixings6,     "Values for variable keyed as flexFixings6");
    FIELD_MAKE_OPTIONAL(flexFixings6);
    FIELD(flexFixings7,     "Values for variable keyed as flexFixings7");
    FIELD_MAKE_OPTIONAL(flexFixings7);
    FIELD(flexFixings8,     "Values for variable keyed as flexFixings8");
    FIELD_MAKE_OPTIONAL(flexFixings8);
    FIELD(flexFixings9,     "Values for variable keyed as flexFixings9");
    FIELD_MAKE_OPTIONAL(flexFixings9);
    FIELD_NO_DESC(fixings1);
    FIELD_MAKE_TRANSIENT(fixings1);
    FIELD_NO_DESC(fixings2);
    FIELD_MAKE_TRANSIENT(fixings2);
    FIELD_NO_DESC(fixings3);
    FIELD_MAKE_TRANSIENT(fixings3);
    FIELD_NO_DESC(fixings4);
    FIELD_MAKE_TRANSIENT(fixings4);
    FIELD_NO_DESC(fixings5);
    FIELD_MAKE_TRANSIENT(fixings5);
    FIELD_NO_DESC(fixings6);
    FIELD_MAKE_TRANSIENT(fixings6);
    FIELD_NO_DESC(fixings7);
    FIELD_MAKE_TRANSIENT(fixings7);
    FIELD_NO_DESC(fixings8);
    FIELD_MAKE_TRANSIENT(fixings8);
    FIELD_NO_DESC(fixings9);
    FIELD_MAKE_TRANSIENT(fixings9);
    clazz->setPublic(); // make visible to EAS/spreadsheet
    Addin::registerConstructor("FLEX_FACTORY_VAR_LOOKUPS",
                               Addin::FLEX_PAYOFF,
                               "Creates a Flex Factory Variable Lookup table",
                               TYPE);
};


//**********************************************//

// Just one of these
static VarCreateMethodLookup varCreateMethods;

void FRVarFactory::registerCreateMethod(TFRVarCreateMethod*  create,
                                        const string&        typeName) {
    static const string routine = "FRVarFactory::registerCreateMethod";
    try {
        VarCreateMethodLookup::iterator iter = varCreateMethods.find(typeName);
        if (iter == varCreateMethods.end()) {
            // if it's not there, add it
            varCreateMethods[typeName] = (void *)create;
        }
        else {
            // if it is we're doubly defined, so complain
            throw ModelException(routine, "Doubly defined create method for type " + typeName);
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
};

FRVarFactory::FRVarFactory(const StringArray&    varTypes,
                           const StringArray&    varNames,
                           const StringArray&    values,
                           FRFactoryVarLookupsSP valueLookups):    
        CObject(TYPE),
        varTypes(varTypes), varNames(varNames),
        values(values), valueLookups(valueLookups)  {
    validatePop2Object();
};

void FRVarFactory::validatePop2Object() {
    static const string routine = "FRVarFactory::validatePop2Object";
    int i;
    int numVars = varTypes.size();
    if (numVars<1) {
        throw ModelException(routine, "At least one variable must be indicated, none given");
    }
    if (varNames.size()!=numVars) {
        throw ModelException(routine, "#varNames=" + Format::toString(varNames.size()) +
                             " and #varTypes=" + Format::toString(numVars) +
                             " should be equal");
    }
    if (values.size()!=numVars) {
        throw ModelException(routine, "#varTypes=" + Format::toString(numVars) +
                             " and #values=" + Format::toString(values.size()) +
                             " should be equal");
    }
    // all names must be unique
    // XXX Will the parser sort this out for me?
    for(i=0; i<numVars; i++) {
        for(int j=i+1; j<numVars; j++) {
            if (varNames[i]==varNames[j]) {
                throw ModelException(routine, "Require unique variable names but " + 
                                     varNames[i] + " is repeated in entries " + 
                                     Format::toString(i+1) + " and " + 
                                     Format::toString(j+1));
            }
        }
    }
 
    // no reason to delay the real work ->
    variables = FRIfaces::ILValueExpressionArraySP(new FRIfaces::ILValueExpressionArray(0));
    for(i=0; i<numVars; i++) {
        VarCreateMethodLookup::iterator iter = varCreateMethods.find(varTypes[i]);
        if (iter == varCreateMethods.end()) {
            throw ModelException(routine, "Failed to find FR variable create method for type " +
                                 varTypes[i]);
        }
        TFRVarCreateMethod* create = (TFRVarCreateMethod*)iter->second;
        bool isLookUp = FRFactoryVarLookups::isLookUp(values[i]);
        StringArray thisValue;
        if (isLookUp) {
            // We have a value per sim date
            if (!valueLookups) {
                throw ModelException(routine, 
                                     "Expected valueLookup but none found");
            }
            thisValue = valueLookups->getLookUp(values[i]);
            // provide the original which may have more info than just the lookup id
            // e.g. for MMRate fixings
            // Note that we subsequently rely on the length of the "values" 
            // StringArray to indicate whether there is a single value to use
            // at all sim dates, or whether a lookup has been employed.
            thisValue.insert(thisValue.begin(), values[i]);
        } else {
            // We have a single value that should apply for all sim dates
            thisValue.push_back(values[i]);
        }
        FRIfaces::ILValueExpressionSP aVar = 
            FRIfaces::ILValueExpressionSP(create(varNames[i], 
                                                 thisValue));
        (*variables).push_back(aVar);
    }

};

FRIfaces::ILValueExpressionSP FRVarFactory::getVariable(string varName) {
    static const string routine = "FRVarFactory::getVariable";
    for(int i=0; i<varNames.size(); i++) {
        if (varNames[i]==varName) {
            return (*variables)[i];
        }
    }
    throw ModelException(routine, "Failed to find variable called " + varName);
};

// Allows others to invoke the constructors
FRIfaces::ILValueExpressionSP FRVarFactory::createVariable(const string&      typeName,
                                                           const string&      varName,
                                                           const StringArray& value) {
    static const string routine = "FRVarFactory::createVariable";
    VarCreateMethodLookup::iterator iter = varCreateMethods.find(typeName);
    if (iter == varCreateMethods.end()) {
        throw ModelException(routine, "Failed to find FR variable create method for type " +
                             typeName);
    }
    TFRVarCreateMethod* create = (TFRVarCreateMethod*)iter->second;
    FRIfaces::ILValueExpressionSP aVar = 
        FRIfaces::ILValueExpressionSP(create(varName, 
                                             value));
    return aVar;
};

FRIfaces::ILValueExpressionArraySP FRVarFactory::getAllVariables() {
    return variables;
};

CClassConstSP const FRVarFactory::TYPE = CClass::registerClassLoadMethod(
    "FRVarFactory", typeid(FRVarFactory), FRVarFactory::load);

CClassConstSP const FlexAlgorithm::TYPE = CClass::registerClassLoadMethod(
    "FlexAlgorithm", typeid(FlexAlgorithm), FlexAlgorithm::load);

CClassConstSP const FRFactoryVarLookups::TYPE = CClass::registerClassLoadMethod(
    "FRFactoryVarLookups", typeid(FRFactoryVarLookups), FRFactoryVarLookups::load);
const string FRFactoryVarLookups::LOOKUP_PREFIX = "PerSimDate";
const string FRFactoryVarLookups::LOOKUP_FIXINGS_PREFIX = "Fixings";
const int FRFactoryVarLookups::MAX_LOOKUPS = 9;

/******************************************************************
 * Allow var creation via Factory in spreadsheet 
 *****************************************************************/
class FRVarCreateAddin: public CObject {
public:
    static CClassConstSP const TYPE;

    // addin parameters
    string             varName;
    string             varType;
    StringArray        value;

    FRIfaces::ILValueExpressionSP createVariable(){
        return FRVarFactory::createVariable(varType, varName, 
                                            value);
    };

    FRVarCreateAddin(): CObject(TYPE) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FRVarCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(varName, "variable name");
        FIELD(varType, "variable type");
        FIELD(value, "Value per sim date as a string array");
        Addin::registerObjectMethod("FR_VAR_CREATE",
                                    Addin::FLEX_PAYOFF,
                                    "Create a variable",
                                    true,
                                    Addin::returnHandle,
                                    &FRVarCreateAddin::createVariable);
    }
    
    static IObject* defaultCtor(){
        return new FRVarCreateAddin();
    }
 
};

CClassConstSP const FRVarCreateAddin::TYPE=CClass::registerClassLoadMethod(
    "FRVarCreateAddin", typeid(FRVarCreateAddin), load);

/* 1. Change return type to be CArraySP. Done.
   2. Add additional Class (ie type) parameter (eg DoubleArray). Done.
   3. Add check for true/false in string section. Done.
   4. Create relevant object to append to array. Done.
   5. Create zero length array to begin with. Done.
   6. Hit number switch on whether we want int or double. Done */
IArraySP FRVarFactory::parseForArrayNames(
    const string&        routine,  // for error message
    const string&        varName,
    const CClassConstSP& clazz, // Array to create eg DoubleArray
    const string&        value) {
    IArraySP values(clazz->newArrayInstance(0));
    CClassConstSP componentType = clazz->getComponentType();
    bool wantDoubles = componentType == CDouble::TYPE;
    // extract var names from "value". Ignore whitespace and recognise
    // comma separated list of names
    const char* pos = value.c_str();
    while (*pos) {
        IObjectSP elt; // element to append to array
        char   c;
        // skip whitespace
        while ((c = *pos) == ' ' || c == '\t'){
            pos++;
        }
        if (c == '-' || c == '.' || isdigit (c)) {
            const char* tmpPos = pos;
            // move over minus signs
            while (*tmpPos == '-'){
                tmpPos++;
            }
            /* see if it's an int */
            while (isdigit(*tmpPos)){
                tmpPos++;
            }
            if (*tmpPos != '.'){ 
                // must be an int
                char* endPos;
                int iVal = strtol(pos, &endPos, 10 /* base 10 */);
                pos = endPos;
                IObject* ptr = wantDoubles? // support casting ints to dbs
                    ((IObject*) CDouble::create(iVal)): 
                    ((IObject*) CInt::create(iVal));
                elt = IObjectSP(ptr);
            } else {
                char* endPos;
                double dVal = strtod(pos, &endPos);
                pos = endPos;
                elt = IObjectSP(CDouble::create(dVal));
            }
        } else if (isalpha(c) || c == '_'){
            const char* varBegin = pos;
            do{
                pos++; /* Get another character. */
                c = *pos;
            } while (c != '\0' && (isalnum(c) || c == '_'));
            const char* varEnd = pos;
            // construct string holding variable name
            string aVarName(varBegin, varEnd - varBegin);
            if (aVarName == "true"){
                elt = IObjectSP(CBool::create(true));
            } else if (aVarName == "false"){
                elt = IObjectSP(CBool::create(false));
            } else {
                elt = IObjectSP(CString::create(aVarName));
            }
        } else if (c == ','){
            // skip over our chosen separator and continue
            pos++;
        } else {
            // not recognised
            throw ModelException(routine, 
                                 "Failed with " + varName + 
                                 " and value of " + value);
        }
        // and, if we've created an element, add it to the list
        if (elt.get()){
            if (!componentType->isInstance(elt)){
                throw ModelException(routine, "For variable "+varName+
                                     " require\narray of "+
                                     componentType->getName()+"s but "
                                     "found a "+elt->getClass()->getName());
            }
            values->append(elt);
        }
    }
    return values;
}



DRLIB_END_NAMESPACE

