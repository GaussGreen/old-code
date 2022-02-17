//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRFactory.hpp
//
//   Description : To satisfy an interface to IMS we need to be able to create
//                 a variable from a "simple" description such as 
//                 "type", "name", "value". The existing constructors are varied
//                 and distributed across several modules, so we have this 
//                 factory which provides a central register of all construction
//                 methods which can create variables from the description above.
//                 Each variable class needs to register a constructor method
//                 with the factory at "load" time.
//
//   Date        : July 2003
//
//
//----------------------------------------------------------------------------
#ifndef EDR_FR_FACTORY_HPP
#define EDR_FR_FACTORY_HPP
#include "edginc/FRIfaces.hpp"
#include "edginc/FR.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

// Allows scalar variable to be initialised with a value per sim date
// via the restrictive IMS interface.
// We note that this approach stops "LookUpN" being a valid scalar input.
class FRFactoryVarLookups : public CObject {
public:
    static CClassConstSP const TYPE;
    static const string LOOKUP_PREFIX;
    static const string LOOKUP_FIXINGS_PREFIX;
    static const int MAX_LOOKUPS; // currently 9, see below

    FRFactoryVarLookups(const FRFactoryVarLookups& rhs); // not implemented
    FRFactoryVarLookups& operator=(const FRFactoryVarLookups& rhs); // not implemented

    FRFactoryVarLookups(const StringArray&        lookUpKeys,
                        const StringArrayArray&   lookUpData);

    virtual void validatePop2Object();

    ~FRFactoryVarLookups() {};

    // for reflection
    FRFactoryVarLookups();

    // does 'id' indicate a lookup?
    static bool isLookUp(const string& id);
    // accessor
    const StringArray& getLookUp(const string& id) const;

    // Helper for IMS-friendly forms
    static void interpretIMSLookup(const StringArray&  IMSInput, // from the lookup key + remaining stuff
                                   StringArray&        lookUpKeys, // with specific 'id's attached
                                   StringArrayArray&   lookUpData); // all the data so far

private:

    static IObject* defaultFRFactoryVarLookups();

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    static int parseLookUp(const string& id);

    static bool isPerSimDate(const string& id);

    // IMS forces us to have this - arrays of arrays not supported
    StringArray  perSimDate1;
    StringArray  perSimDate2;
    StringArray  perSimDate3;
    StringArray  perSimDate4;
    StringArray  perSimDate5;
    StringArray  perSimDate6;
    StringArray  perSimDate7;
    StringArray  perSimDate8;
    StringArray  perSimDate9; // corresponds to MAX_LOOKUPS

    CashFlowArray flexFixings1;
    CashFlowArray flexFixings2;
    CashFlowArray flexFixings3;
    CashFlowArray flexFixings4;
    CashFlowArray flexFixings5;
    CashFlowArray flexFixings6;
    CashFlowArray flexFixings7;
    CashFlowArray flexFixings8;
    CashFlowArray flexFixings9;

    // these correspond to the above but are in a more
    // convenient format
    StringArray fixings1;
    StringArray fixings2;
    StringArray fixings3;
    StringArray fixings4;
    StringArray fixings5;
    StringArray fixings6;
    StringArray fixings7;
    StringArray fixings8;
    StringArray fixings9;

    // transient
    vector<StringArray*> lookUps; // $unregistered
};
typedef smartPtr<FRFactoryVarLookups> FRFactoryVarLookupsSP;

typedef FR::LValueBase* (TFRVarCreateMethod)(const string&      varName, 
                                             const StringArray& value);

// Idea is that empty "value" means non-const; with "value" create a const variety
// "value" is parsed into int/double etc
// For arrays value may be comma separated list, so x = "a,b,c"
// "Lookups" offer a scalar-per-sim-date initialisation
class PRODUCTS_DLL FRVarFactory : public CObject {
private:
    StringArray                        varTypes;
    StringArray                        varNames;
    StringArray                        values;
    FRFactoryVarLookupsSP              valueLookups;      

    FRIfaces::ILValueExpressionArraySP variables;

public:
    static CClassConstSP const TYPE;

    FRVarFactory(): CObject(TYPE) {}
    FRVarFactory(const FRVarFactory& rhs); // not implemented
    FRVarFactory& operator=(const FRVarFactory& rhs); // not implemented

    FRVarFactory(const StringArray&    varTypes,
                 const StringArray&    varNames,
                 const StringArray&    values,
                 FRFactoryVarLookupsSP valueLookups = FRFactoryVarLookupsSP(   ));
    
    virtual void validatePop2Object();

    static IObject* defaultFRVarFactory(){
        return new FRVarFactory();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FRVarFactory, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFRVarFactory);
        FIELD(varTypes,     "Each variable's type");
        FIELD(varNames,     "Each variable's name");
        FIELD(values,       "Each variable's value");
        FIELD(valueLookups, "[Optional] variable's value per sim date");
        FIELD_MAKE_OPTIONAL(valueLookups);
        FIELD(variables,           "Synthesised forms");
        FIELD_MAKE_TRANSIENT(variables);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FLEX_VAR_FACTORY",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a Flex Rules Variable Factory",
                                   TYPE);
    }

    static void registerCreateMethod(TFRVarCreateMethod*  create,
                                     const string&        typeName);

    // Obvious lookup. Convenient for payout variable
    FRIfaces::ILValueExpressionSP getVariable(string varName);

    // Convenient packaging of whole array
    FRIfaces::ILValueExpressionArraySP getAllVariables();

    // Allows others to invoke the constructors
    static FRIfaces::ILValueExpressionSP createVariable(const string&      typeName,
                                                        const string&      varName,
                                                        const StringArray& valuePerDate);

    static IArraySP parseForArrayNames(
        const string&        routine,  // for error message
        const string&        varName,
        const CClassConstSP& clazz, // Array to create eg DoubleArray
        const string&        value);
    
};

typedef smartPtr<FRVarFactory> FRVarFactorySP;

// The one to be used for IMS
// The parameters are all essentially basic types, or arrays of such
// Note the definition of this class is a bit dodgy - validatePop2Object
// needs visibility of FRParserTPRules and so is defined there.
// and getAlgorithm needs visibility of TPRSimpleArray and so is defined
// in FR.cpp
class PRODUCTS_DLL FlexAlgorithm : public CObject {
private:
    IntArray                           simDateRuleIds; // [#simDates]
    IntArray                           ruleId;         // [#rule lines]
    StringArray                        ruleDefn;       // [#rule lines]
    FRIfaces::ITimePointRulesArraySP   justTheRules;
public:
    static CClassConstSP const TYPE;

    virtual void validatePop2Object();

    FRIfaces::IAlgorithmSP getAlgorithm(const DateTimeArraySP simDates);

    //FRAladdinAlgorithmSP getAlgorithm(const DateTimeArraySP simDates);

    FlexAlgorithm(): CObject(TYPE) {}

    FlexAlgorithm(const IntArray& simDateRuleIds,
                  const IntArray& ruleId,
                  const StringArray& ruleDefn):    CObject(TYPE),
                  simDateRuleIds(simDateRuleIds),
                  ruleId(ruleId),
                  ruleDefn(ruleDefn)  {}

private:
    FlexAlgorithm(const FlexAlgorithm& rhs); // not implemented
    FlexAlgorithm& operator=(const FlexAlgorithm& rhs); // not implemented

    static IObject* defaultFlexAlgorithm(){
        return new FlexAlgorithm();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(FlexAlgorithm, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultFlexAlgorithm);
        FIELD(simDateRuleIds,      "Maps rule IDs to sim dates");
        FIELD(ruleId,              "Identifies each line of a rule");
        FIELD(ruleDefn,            "Defines each line of fule identified by ruleId");
        FIELD(justTheRules,               "justTheRules");
        FIELD_MAKE_TRANSIENT(justTheRules);
        clazz->setPublic(); // make visible to EAS/spreadsheet
        Addin::registerConstructor("FLEX_ALGORITHM",
                                   Addin::FLEX_PAYOFF,
                                   "Creates a set of Flex Rules (IMS style)",
                                   TYPE);
    }
};

typedef smartPtr<FlexAlgorithm> FlexAlgorithmSP;

DRLIB_END_NAMESPACE
#endif
