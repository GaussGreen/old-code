//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Addin.cpp
//
//   Description : Class for representing addins
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Modifier.hpp"
#include ext_hash_map
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** Functor for static functions (for support for old style [static] addin
    functions */
template<class ReturnType> class AddinFunctorStatic: 
    public Functor<ReturnType> {
public:
    typedef ReturnType (*Func)(void*);
    // constructor
    AddinFunctorStatic(Func func): func(func) {}
    ~AddinFunctorStatic(){} // and destructor
    // Implementation of execute
    virtual ReturnType execute(IObject* obj) {
        // cast from IObject* to the object and call method
        return func(CClass::castFromIObject(obj));
    }
private:
    Func func;      // Pointer to static function
};

/** Class for representing addins. Each instance
    of this class corresponds to one excel addin function */

/** Categories for addin functions */
const string Addin::UTILITIES = "QLIB Utilities";
const string Addin::XL_TESTS = "QLIB Excel Testers";
const string Addin::RISK = "QLIB Risk Management";
const string Addin::MARKET = "QLIB Market Data";
const string Addin::FLEX_PAYOFF = "QLIB Flex Payoff";
const string Addin::CONV_BOND = "QLIB CVB Addins";

/** Overrides "EDR " with "newPrefix ". Needs to be called before
    loadAllClasses() */
void Addin::overrideCategoryPrefix(const string& newPrefix){
    // this is really weak and hacky. But the correct solution is
    // to move the defines into a separate class which is too painful.
    // The other workarounds are pretty weak too
    const_cast<string&>(UTILITIES).replace(0, 3, newPrefix);
    const_cast<string&>(XL_TESTS).replace(0, 3, newPrefix);
    const_cast<string&>(RISK).replace(0, 3, newPrefix);
    const_cast<string&>(MARKET).replace(0, 3, newPrefix);
    const_cast<string&>(FLEX_PAYOFF).replace(0, 3, newPrefix);
    const_cast<string&>(CONV_BOND).replace(0, 3, newPrefix);
}

/** Hash function needed for strings */
struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};


typedef hash_map<string, Addin*, MyStringHash> AddinHash;

/** Holder for hash table - avoid slowing header file */
class AddinHelper{
public:
    static AddinHash addinHash;
};

AddinHash AddinHelper::addinHash;

/** The supplied method will be called once for each addin registered.
    If kitMethod == 0, then no method is called */
void Addin::registerKit(KitRegister *kitMethod){
    if (!kitMethod){
        return;
    }
    // loop though hash 
    for (AddinHash::const_iterator iter = AddinHelper::addinHash.begin();
         !(iter == AddinHelper::addinHash.end()); ++iter)
       {
        const Addin*   addin = iter->second;
        try
        {
            // need to get decription of parameters - including those for
            // all parent fields
            CFieldArray   allFields(getDataClassFields(addin->dataClass));
            kitMethod(addin->addinName, addin->category,
                        addin->dataClass, 
                        allFields,
                        addin->handleNameParam,
                        addin->returnStyle, addin->returnIsNative,
                        addin->methodType, 
                        addin->method, 
                        addin->returnType);
        }
        catch (...)
        {
            throw ModelException("Addin::registerKit", "Problem registering addin function, "+
                                 addin->addinName +" already exists.");
        }
    }
}

Addin::~Addin(){}

Addin::Addin(
    const string&          addinName,    // excel name
    const string&          category,     // category - for XL
    const string&          description,  // what the addin does
    CClassConstSP          dataClass,   // identifies class representing addin
    bool                   handleNameParam, /* true: extra input for
                                               handle name */
    ReturnStyle            returnStyle,  // how to return the output
    MethodType             methodType,   /* constructor/instanceMethod/
                                            classMethod */
    Method                 method,       // holds function pointer
    VirtualDestructorBase* objToFree,    // delete this
    bool                   returnIsNative, /* true - parameter is returned
                                              as an atomic type */
    CClassConstSP          returnType):  /* the type of the return
                                            parameter */
    addinName(addinName), category(category), description(description),
    dataClass(dataClass),
    handleNameParam(handleNameParam), returnStyle(returnStyle),
    methodType(methodType), method(method), objToFree(objToFree),
    returnIsNative(returnIsNative), returnType(returnType){}


/** Add the given addin to hashtable containing all addins. Fails if
    addin with given name exists already. Takes ownwership of memory.
    If successfully added, then invokes kit register method */
void Addin::addToHashtable(Addin*   addin){
    const static string routine = "Addin::addToHashtable";
    string addinName = addin->addinName;
    try{
        AddinHash::const_iterator iter = 
            AddinHelper::addinHash.find(addinName);
        if (!(iter == AddinHelper::addinHash.end())){
            delete addin;
            throw ModelException(routine, "Addin with name "+
                                 addinName+" already exists");
        }
        AddinHelper::addinHash[addinName] = addin;
    } catch (exception& e){
        throw ModelException(e, routine, "Failed for "+addinName);
    }
}

/** Return an array of fields representing each parameter of the addin
    Fields are returned in a particular order */
CFieldArray Addin::getDataClassFields() const{
    return getDataClassFields(dataClass);
}

/** Returns the class identifying the parameters that the addin takes */
CClassConstSP  Addin::getDataClass() const{
    return dataClass;
}

/** Return an array of fields representing each parameter in the addin -
    this includes those for all parent fields. Skips transient fields */
CFieldArray Addin::getDataClassFields(CClassConstSP clazz){    
    CFieldArray   allFields(0);
    do {
        // loop over child and parents 
        const CFieldArray& fields = clazz->getDeclaredFields();
        allFields.reserve(fields.size());
        // go backwards through the fields (we reverse the list at the end)
        for (unsigned int i = fields.size(); i > 0 ; i--){
            // skip transient fields
            if (!Modifier::isTransient(fields[i-1]->getModifiers())){
                allFields.push_back(fields[i-1]);
            }
        }
    } while ((clazz = clazz->getSuperClass()));
    // then reverse the list - this is a style decision on the user
    // input ie ensure that the parent members become before the childs
    int numFields = allFields.size();
    for (int i = 0; i < (numFields >> 1); i++){
        CFieldConstSP field = allFields[i];
        allFields[i] = allFields[numFields -1 - i];
        allFields[numFields -1 - i] = field;
    }
    return allFields;
}
 
/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a double */
void Addin::registerClassDoubleMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    DoubleMethod*      method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerClassDoubleMethod",
                             "DoubleMethod must not be null");
    }
    Method addinMethod;
    addinMethod.doubleMethod = new AddinFunctorStatic<double>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, classMethod, 
                             addinMethod, addinMethod.doubleMethod,
                             true, CDouble::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns an int */
void Addin::registerClassIntMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    IntMethod*         method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerClassIntMethod",
                             "IntMethod must not be null");
    }
    Method addinMethod;
    addinMethod.intMethod = new AddinFunctorStatic<int>(method);;
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, classMethod, 
                             addinMethod, addinMethod.intMethod,
                             true, CInt::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a bool */
void Addin::registerClassBoolMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    BoolMethod*        method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerClassBoolMethod",
                             "BoolMethod must not be null");
    }
    Method addinMethod;
    addinMethod.boolMethod = new AddinFunctorStatic<bool>(method);;
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, classMethod, 
                             addinMethod, addinMethod.boolMethod,
                             true, CBool::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns an IObject*  */
void Addin::registerClassObjectMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    bool               handleNameParam, //true: extra input for handle name
    ReturnStyle        returnStyle, // how to return the output
    ObjMethod*         method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerClassObjectMethod", 
                             "ObjMethod must not be null");
    }
    Method addinMethod;
    addinMethod.objMethod = new AddinFunctorStatic<IObjectSP>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             handleNameParam, returnStyle, classMethod, 
                             addinMethod, addinMethod.objMethod,
                             false, IObject::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a string. */
void Addin::registerClassStringMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    StringMethod*      method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerClassStringMethod", 
                             "StringMethod must not be null");
    }
    Method addinMethod;
    addinMethod.stringMethod = new AddinFunctorStatic<string>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, classMethod, 
                             addinMethod, addinMethod.stringMethod,
                             true, CString::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a double */
void Addin::registerInstanceDoubleMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    DoubleMethod*      method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerInstanceDoubleMethod",
                             "DoubleMethod must not be null");
    }
    Method addinMethod;
    addinMethod.doubleMethod = new AddinFunctorStatic<double>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, instanceMethod, 
                             addinMethod, addinMethod.doubleMethod,
                             true, CDouble::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns an IObject*  */
void Addin::registerInstanceObjectMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    bool               handleNameParam, //true: extra input for handle name
    ReturnStyle        returnStyle, // how to return the output
    ObjMethod*         method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerInstanceObjectMethod", 
                             "ObjMethod must not be null");
    }
    Method addinMethod;
    addinMethod.objMethod = new AddinFunctorStatic<IObjectSP>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             handleNameParam, returnStyle, instanceMethod, 
                             addinMethod, addinMethod.objMethod,
                             false, IObject::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a string. */
void Addin::registerInstanceStringMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    StringMethod*      method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerInstanceStringMethod", 
                             "StringMethod must not be null");
    }
    Method addinMethod;
    addinMethod.stringMethod = new AddinFunctorStatic<string>(method);
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, instanceMethod, 
                             addinMethod, addinMethod.stringMethod,
                             true, CString::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns an int */
void Addin::registerInstanceIntMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    IntMethod*         method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerInstanceIntMethod",
                             "IntMethod must not be null");
    }
    Method addinMethod;
    addinMethod.intMethod = new AddinFunctorStatic<int>(method);;
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, instanceMethod, 
                             addinMethod, addinMethod.intMethod,
                             true, CInt::TYPE);
    addToHashtable(addin); // takes ownership of memory
}
/** Record data defining an addin which, if written in C++,
    correponds to a static class method which returns a bool */
void Addin::registerInstanceBoolMethod(
    const string&      addinName,   // excel name
    const string&      category,    // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass,  /* identifies class 
                                           representing addin data */
    BoolMethod*      method)      // holds function pointer
{
    if (!method){
        throw ModelException("Addin::registerInstanceBoolMethod",
                             "BoolMethod must not be null");
    }
    Method addinMethod;
    addinMethod.boolMethod = new AddinFunctorStatic<bool>(method);;
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             false, expandSimple, instanceMethod, 
                             addinMethod, addinMethod.boolMethod,
                             true, CBool::TYPE);
    addToHashtable(addin); // takes ownership of memory
}
    

/** Record data defining an addin which, if written in C++,
    correponds to a constructor method which returns a handle */
void Addin::registerConstructor(
    const string&      addinName,       // excel name
    const string&      category,        // category - for XL
    const string&      description, // what the addin does
    CClassConstSP      addinDataClass)  /* identifies class 
                                           representing object to build */
{
    Method addinMethod;
    addinMethod.objMethod = 0;
    Addin* addin = new Addin(addinName, category, description, addinDataClass,
                             true, returnHandle, constructor, 
                             addinMethod, 0, true, IObject::TYPE);
    addToHashtable(addin); // takes ownership of memory
}

/** Same as above but addinName is defaulted to the class's name so
    eg EDR_MyObject and the description is defaulted to "creates an "
    "object of type ..." */
void Addin::registerConstructor(
    const string&      category,        // category - for XL
    CClassConstSP      addinDataClass)  /* identifies class 
                                           representing object to build */
{
    const string& clazzName = addinDataClass->getName();
    string regName = clazzName;

    for (string::size_type i=0; i<regName.size(); ++i) {
        if (regName[i]==':') regName[i]='_';
    }
    registerConstructor(regName, category,
                        "Creates an object of type "+clazzName,
                        addinDataClass);
}

/** Look up the addin corresponding to the given string. This method
    has no knowledge of the EDR_ addin prefix */
const Addin* Addin::lookUp(const string& addinName){
    AddinHash::const_iterator iter = AddinHelper::addinHash.find(addinName);
    if (iter == AddinHelper::addinHash.end()){
        throw ModelException("Addin::lookUp", "No addin known with name "+
                             addinName);
    }
    return iter->second;
}

/** Indicates whether this addin method returns a native type or not */
bool Addin::isReturnNative() const{
    return returnIsNative;
}

/** Returns the Addin::Method union for this addin (ie gives function
    pointer for method) */
Addin::Method Addin::getMethod() const{
    return method;
}

/** Returns the type of parameter that this addin returns */
CClassConstSP Addin::getReturnType() const{
    return returnType;
}

/** Indicates whether this addin has an extra parameter for a handle
    name */
bool Addin::hasHandleNameParam() const{
    return handleNameParam;
}

/** Returns a description of what the addin does */
const string& Addin::getDescription() const{
    return description;
}

/** Returns category for the addin */
const string& Addin::getCategory() const {
    return category;
}

/** Returns an array of names sorted into alphabetical order */
CStringArraySP Addin::names(){
    // create empty array
    CStringArraySP  stringArray(new CStringArray(0));
    CStringArray&   array = *stringArray;
    // loop though hash appending to string list
    for (AddinHash::const_iterator iter = AddinHelper::addinHash.begin();
         !(iter == AddinHelper::addinHash.end()); ++iter){
        array.push_back(iter->first);
    }
    // now sort them
    sort(array.begin(), array.end());
    return stringArray;
}
    

DRLIB_END_NAMESPACE
