
#include "edginc/config.hpp"
#define QLIB_ATOMIC_CPP
#include "edginc/Format.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Null.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Hashtable.hpp"
#include ext_hash_map
#include ext_hash_set

DRLIB_BEGIN_NAMESPACE
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CDouble>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CDouble>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CInt>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CInt>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CBool>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CBool>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CString>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CString>);

//// IDRObject class declared in DRObject.hpp (no source file)
static void loadInterface(CClassSP& clazz){
    clazz->setPublic(); /* want this in EDR reflection but not in
                           DRI reflection. No way to do this currently */
    REGISTER_INTERFACE(IDRObject, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IDRObject::TYPE = 
CClass::registerInterfaceLoadMethod("IDRObject", typeid(IDRObject), 
                                    loadInterface);

/** Invoked when this class is 'loaded' */
static void loadAtomic(CClassSP& classToLoad){
    classToLoad->setPublic(); // make visible to EAS/spreadsheet
    classToLoad->setNative();
}

CDouble::~CDouble(){};

/** Are these objects equal (ie contain the same data) */
bool CDouble::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const CDouble* dble = STATIC_CAST(CDouble, drObj.get());
    return (dble->db == db);
}

//// gcc messes up the hashCode function eg hashCode(0) != 0;
void gccWorkAroundFix(int d, ...){}

/** Returns a hash code value for the supplied double */
int CDouble::hashCode(double db){
    // need to verify this is ok - implicitly assuming double is 8 bytes and
    // an int is 4 bytes.
    int part1 = (const int&)db;
    gccWorkAroundFix(part1, part1);
    int part2 = *(const int*)(((const char*)&db) + 4);
    gccWorkAroundFix(part2, part2);
    return (part1 ^ part2);
}

/** Returns a hash code value for the object */
int CDouble::hashCode() const{
    return hashCode(db);
}

/** Are these objects equal (ie contain the same data) */
bool CDouble::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (obj->getClass() != TYPE){
        return false;
    }
    const CDouble* dble = STATIC_CAST(CDouble, obj);
    return (dble->db == db);
}

IObject* CDouble::clone() const{
    /* would like to just in increase refCount by 1 and return this but
       somebody has made the object no longer immutable (scale and add) */
    return new CDouble(db);
}

// wrapper for double object
CDouble* CDouble::create(double d){
    return new CDouble(d);
}

CDoubleSP CDouble::SP(double d){
    return CDoubleSP(new CDouble(d));
}
     
CDouble::CDouble(double d): CObject(TYPE), db(d){}

double CDouble::doubleValue() const{
    return db;
}
    
void CDouble::write(const string& tag, Writer* writer) const{
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            char buffer[512];  // can cope with biggest possible double
            sprintf(buffer, "%.16f", db);
            writer->write(buffer);
        }
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(&e, "CDouble::write");
    }
}

/** populate an empty object from reader */
void CDouble::import(Reader::Node* elem, Reader* reader) {
    string s = elem->value();
    db = atof(s.c_str());
}

void CDouble::outputWrite(const string& linePrefix,
                          const string& prefix, ostream& stream) const{
    
    ios::fmtflags oldFlags = stream.flags(); // save settings
    long oldPrecision = stream.precision(); // save settings
    if (db < 1 && db > -1){
        // 8 dp
        stream.flags(ios::fixed);
    } else {
        // 8 sf
        stream.unsetf(ios::fixed);
    }        
    stream.precision(8);
    stream << linePrefix << prefix << ": " << db << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void CDouble::scale(double x) {
    db *= x;
}

/** add an object to this result */
void CDouble::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
  const CDouble& dbl = dynamic_cast<const CDouble&>(static_cast<const IObject&>(x));
    db += scaleFactor * dbl.db;           
}

string CDouble::toString() const {
    return Format::toString(doubleValue());
}

class DoubleHelper{
public:
    static IObject* defaultDouble(){
        return CDouble::create(0.0);
    }

    static void loadDouble(CClassSP& classToLoad){
        loadAtomic(classToLoad);
        REGISTER(CDouble, classToLoad);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDRObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultDouble);
    }
};


    
CClassConstSP const CDouble::TYPE = CClass::registerClassWithAliasLoadMethod(
    "Double", typeid(double), typeid(CDouble), DoubleHelper::loadDouble);


/** Class for addins that take a single optional double */
class DoubleAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // default constructor
    DoubleAddin(): CObject(TYPE), doubleValue(0){}

    // single optional field
    double* doubleValue;

    ~DoubleAddin(){
        delete doubleValue;
    }

private:
    /** addin function - either creates a CDouble or a Null if no double is
        supplied */
    static IObjectSP createDouble(DoubleAddin* params){
        if (!params->doubleValue){
            return CNull::create();
        }
        return IObjectSP(CDouble::create(*params->doubleValue));
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DoubleAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDoubleAddin);
        FIELD(doubleValue, "Double Value");
        FIELD_MAKE_OPTIONAL(doubleValue);
        Addin::registerClassObjectMethod("DOUBLE",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a double",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createDouble);
    }

    static IObject* defaultDoubleAddin(){
        return new DoubleAddin();
    }
    
};

CClassConstSP const DoubleAddin::TYPE = CClass::registerClassLoadMethod(
    "DoubleAddin", typeid(DoubleAddin), load);

CInt::~CInt(){}

/** Are these objects equal (ie contain the same data) */
bool CInt::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const CInt* theInt = STATIC_CAST(CInt, drObj.get());
    return (theInt->i == i);
}

/** Returns a hash code value for the object */
int CInt::hashCode() const{
    return i;
}

/** Are these objects equal (ie contain the same data) */
bool CInt::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (obj->getClass() != TYPE){
        return false;
    }
    const CInt* theInt = STATIC_CAST(CInt, obj);
    return (theInt->i == i);
}

IObject* CInt::clone() const{
    int count = getRefCount();
    if (count == 0){
        return new CInt(i);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

CIntSP CInt::SP(int i){
    return CIntSP(new CInt(i));
}

CInt* CInt::create(const int d){
    return new CInt(d);
}

     
CInt::CInt(int d): CObject(TYPE), i(d){}

int CInt::intValue() const{
    return i;
}
    
/** write object out to writer */
void CInt::write(const string& tag,  Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            char buffer[40];
            sprintf(buffer, "%d", i);
            writer->write(buffer);
        }
        writer->objectEnd(tag, this);
    } catch (exception& e){
        throw ModelException(&e, "CInt::write");
    }
}

/** populate an empty object from reader */
void CInt::import(Reader::Node* elem, Reader* reader){
    string s = elem->value();
    // import is privileged
    const_cast<int&>(i) = atoi(s.c_str());
}

void CInt::outputWrite(const string& linePrefix,
                       const string& prefix, ostream& stream) const{
    stream << linePrefix << prefix << ": " << i << endl;
}

string CInt::toString() const {
    return Format::toString(intValue());
}

class IntHelper{
public:
    static IObject* defaultInt(){
        return CInt::create(0);
    }
    
    static void load(CClassSP& classToLoad){
        loadAtomic(classToLoad);
        REGISTER(CInt, classToLoad);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDRObject);
        EMPTY_SHELL_METHOD(defaultInt);
    }
};

CClassConstSP const CInt::TYPE = CClass::registerClassWithAliasLoadMethod(
    "Int", typeid(int), typeid(CInt), IntHelper::load);

/** Class for addins that take a single optional int */
class IntAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // default constructor
    IntAddin(): CObject(TYPE), intValue(0){}

    // single optional field
    int* intValue;

    ~IntAddin(){
        delete intValue;
    }

private:
    /** addin function - either creates a CInt or a Null if no int is
        supplied */
    static IObjectSP createInt(IntAddin* params){
        if (!params->intValue){
            return CNull::create();
        }
        return IObjectSP(CInt::create(*params->intValue));
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IntAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIntAddin);
        FIELD(intValue, "Int Value");
        FIELD_MAKE_OPTIONAL(intValue);
        Addin::registerClassObjectMethod("INT",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a int",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createInt);
    }

    static IObject* defaultIntAddin(){
        return new IntAddin();
    }
    
};

CClassConstSP const IntAddin::TYPE = CClass::registerClassLoadMethod(
    "IntAddin", typeid(IntAddin), load);

CBool::~CBool(){}

/** Are these objects equal (ie contain the same data) */
bool CBool::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const CBool* theBool = STATIC_CAST(CBool, drObj.get());
    return (theBool->b == b);
}

/** Returns a hash code value for the object */
int CBool::hashCode() const{
    return hashCode(b);
}

/** Returns a hash code value for the supplied boolean */
int CBool::hashCode(bool b){
    return (b? 1231: 1237); // as per java!
}


/** Are these objects equal (ie contain the same data) */
bool CBool::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (obj->getClass() != TYPE){
        return false;
    }
    const CBool* theBool = STATIC_CAST(CBool, obj);
    return (theBool->b == b);
}

IObject* CBool::clone() const{
    int count = getRefCount();
    if (count == 0){
        return new CBool(b);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}

CBool* CBool::create(const bool d){
    return new CBool(d);
}
     
CBool::CBool(bool d): CObject(TYPE), b(d){}

bool CBool::boolValue() const{
    return b;
}

void CBool::write(const string& tag, Writer* writer) const {
    
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(b ? "Y" : "N");
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "CBool::write");
    }
}

/* returns true for strings "Y", "YES", "1", "TRUE" and false for strings
   "N", "NO", "0", "FALSE" else throws an exception. Comparison is
   case independent */
bool CBool::fromString(const string& value){
    bool myBool;
    if (CString::equalsIgnoreCase("Y", value) ||
        CString::equalsIgnoreCase("YES", value) ||
        CString::equalsIgnoreCase("1", value) ||
        CString::equalsIgnoreCase("TRUE", value)){
        myBool = true;
    } else if (CString::equalsIgnoreCase("N", value) ||
               CString::equalsIgnoreCase("NO", value) ||
               CString::equalsIgnoreCase("0", value) ||
               CString::equalsIgnoreCase("FALSE", value)){
        myBool = false;
    } else {
        throw ModelException("CBool::fromString", 
                             "Unrecognised value: "+value);
    }
    return myBool;
}

// so can convert strings to bools (e.g. StringArray->BoolArray) when building
// via DataDictionary
IObject* CBool::makeFromString(CClassConstSP requiredType, const string& data) {
    if (!CBool::TYPE->isAssignableFrom(requiredType)){
        throw ModelException("CBool::makeFromString", requiredType->getName()+
                             " not derived from Bool");
    }

    return CBool::create(fromString(data));
}


/** populate an empty object from reader */
void CBool::import(Reader::Node* elem, Reader* reader) {
    string s = elem->value();
    // import is privileged
    const_cast<bool&>(b) = (s == "Y") || (s == "1");
}

void CBool::outputWrite(const string& linePrefix,
                        const string& prefix, ostream& stream) const{
    stream << linePrefix << prefix << ": " << (b? "TRUE": "FALSE") << endl;
}

string CBool::toString() const {
    return Format::toString(boolValue());
}

class BoolHelper{
public:
    
    static IObject* defaultBool(){
        return CBool::create(false);
    }
    
    static void load(CClassSP& classToLoad){
        loadAtomic(classToLoad);
        REGISTER(CBool, classToLoad);
        SUPERCLASS(CObject);
        IMPLEMENTS(IDRObject);
        EMPTY_SHELL_METHOD(defaultBool);
        CString::registerObjectFromStringMethod(CBool::TYPE,
                                                CBool::makeFromString);
    }
    
};

CClassConstSP const CBool::TYPE = CClass::registerClassWithAliasLoadMethod(
    "Bool", typeid(bool), typeid(CBool), BoolHelper::load);

/** Class for addins that take a single optional bool */
class BoolAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // default constructor
    BoolAddin(): CObject(TYPE), boolValue(0){}

    // single optional field
    bool* boolValue;

    ~BoolAddin(){
        delete boolValue;
    }

private:
    /** addin function - either creates a CBool or a Null if no bool is
        supplied */
    static IObjectSP createBool(BoolAddin* params){
        if (!params->boolValue){
            return CNull::create();
        }
        return IObjectSP(CBool::create(*params->boolValue));
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(BoolAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultBoolAddin);
        FIELD(boolValue, "Bool Value");
        FIELD_MAKE_OPTIONAL(boolValue);
        Addin::registerClassObjectMethod("BOOLEAN",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a bool",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createBool);
    }

    static IObject* defaultBoolAddin(){
        return new BoolAddin();
    }
    
};

CClassConstSP const BoolAddin::TYPE = CClass::registerClassLoadMethod(
    "BoolAddin", typeid(BoolAddin), load);

CString::~CString(){}

/** Are these objects equal (ie contain the same data) */
bool CString::equals(IDRObjectConstSP drObj) const{
    if (this == drObj.get()){
        return true;
    }
    if (drObj->getClass() != TYPE){
        return false;
    }
    const CString* theStr = STATIC_CAST(CString, drObj.get());
    return (theStr->s == s);
}

/** Returns a hash code value for the object */
int CString::hashCode() const{
    return hash_string(s);
}

/** Are these objects equal (ie contain the same data) */
bool CString::equalTo(const IObject* obj) const{
    if (this == obj){
        return true;
    }
    if (obj->getClass() != TYPE){
        return false;
    }
    const CString* theStr = STATIC_CAST(CString, obj);
    return (theStr->s == s);
}

IObject* CString::clone() const{
    int count = getRefCount();
    if (count == 0){
        return new CString(s);
    } else {
        return  const_cast<IObject*>((const IObject*)this);
    }
}
CString* CString::create(const string& s){
    return new CString(s);
}
     
CString::CString(string s): CObject(TYPE), s(s){}

const string& CString::stringValue() const{
    return s;
}

static bool nocase_compare(char c1, char c2) {
    return (toupper(c1) == toupper(c2));
}

bool CString::equalsIgnoreCase(const string& s1, const string& s2) {
    if (s1.size() == s2.size() &&
        equal(s1.begin(), s1.end(),
              s2.begin(),
              nocase_compare)) {
        return true;
    }
    else {
        return false;
    }
}

/** like strncmp - does s2 match first n characters of s1 ? */
bool CString::equalsIgnoreCase(const string& s1, const string& s2, int n) {
    int strLen = s1.length();
    int lenToCompare = n < strLen? n: strLen;
    for (int i = 0; i < lenToCompare; i++) {
        if (!nocase_compare(s1[i], s2[i])) {
            return false;
        }
    }
    return true;
}

void CString::write(const string& tag, Writer* writer) const {
    try {    
        IObjectConstSP obj(writer->objectStart(tag, "", this, false));
        if (obj.get()){
            writer->write(s);
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "CString::write");
    }
}

void CString::import(Reader::Node* elem, Reader* reader) {
    const_cast<string&>(s) = elem->value();
}

void CString::outputWrite(const string& linePrefix,
                          const string& prefix, ostream& stream) const{
    stream << linePrefix << prefix << ": " << s << endl;
}

/** remove duplicates and empty strings from an array. Order is otherwise
    unchanged */
StringArraySP CString::trim(const StringArray& names) {
    hash_set<string, Hashtable::StringHash> s;
    StringArraySP trimmed(new StringArray());
    for (int i = 0; i < names.size(); i++){
        const string& name = names[i];
        if (!name.empty() && s.find(name) == s.end()){
            trimmed->push_back(name);
            s.insert(name);
        }
    }
    return trimmed;
}

string CString::toString() const {
    return stringValue();
}

typedef hash_map<CClassConstSP, CString::TObjectFromString*, 
    CClass::Hash> ObjectFromStringHash;

class StringHelper{
public:
    // storage of methods for converting strings into obects
    static ObjectFromStringHash convertMethods;

    static IObject* defaultString(){
        return CString::create("");
    }
    
    static void load(CClassSP& classToLoad){
        loadAtomic(classToLoad);
        REGISTER(CString, classToLoad);
        SUPERCLASS(CObject);
        IMPLEMENTS(ITypeConvert);
        IMPLEMENTS(IDRObject);
        EMPTY_SHELL_METHOD(defaultString);    
    }
};

ObjectFromStringHash StringHelper::convertMethods;

/** Register a conversion method for a given type. Note that the
    method will be invoked for derived classes of the supplied
    type */
void CString::registerObjectFromStringMethod(CClassConstSP      targetClass,
                                             TObjectFromString* method){
    if (!method){
        throw ModelException("CClassSP::registerObjectFromStringMethod",
                             "Null method for class "+targetClass->getName());
    }
    StringHelper::convertMethods[targetClass] = method;
}

/** Converts this object to an instance of the requiredType. Throws an
    exception if a conversion to the required Type is not supported */
void CString::convert(IObjectSP&    object,
                      CClassConstSP requiredType) const
{
    CClassConstSP c   = requiredType;
    // loop through class and superclasses to see if we can find a match
    while (c){
        ObjectFromStringHash::const_iterator iter = 
            StringHelper::convertMethods.find(c);
        if (!(iter == StringHelper::convertMethods.end())){
             object = IObjectSP(iter->second(requiredType, s));
             return;
        } else {
            c = c->getSuperClass();
        }
    }
    throw ModelException("CString::convert", "Cannot convert a string into "
                         "object of type "+requiredType->getName());
}

CClassConstSP  const CString::TYPE = CClass::registerClassWithAliasLoadMethod(
    "String", typeid(string), typeid(CString), StringHelper::load);

/** Class for addins that take a single optional string */
class StringAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // default constructor
    StringAddin(): CObject(TYPE), stringValue(0){}

    // single optional field
    string* stringValue;

    ~StringAddin(){
        delete stringValue;
    }

private:
    /** addin function - either creates a CString or a Null if no string is
        supplied */
    static IObjectSP createString(StringAddin* params){
        if (!params->stringValue){
            return CNull::create();
        }
        return IObjectSP(CString::create(*params->stringValue));
    }

    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StringAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultStringAddin);
        FIELD(stringValue, "String Value");
        FIELD_MAKE_OPTIONAL(stringValue);
        Addin::registerClassObjectMethod("STRING",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a string",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createString);
    }

    static IObject* defaultStringAddin(){
        return new StringAddin();
    }
    
};

CClassConstSP const StringAddin::TYPE = CClass::registerClassLoadMethod(
    "StringAddin", typeid(StringAddin), load);

/** Addin for joining array of string together - for overcoming excel 255
    character limitation */
class StringJoinAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // default constructor
    StringJoinAddin(): CObject(TYPE){}

    // single field
    StringArraySP strings;

private:
    /** addin function */
    static IObjectSP createString(StringJoinAddin* params){
        string newString;
        for (int i = 0; i < params->strings->size(); i++){
            newString += (*params->strings)[i];
        } 
        return IObjectSP(CString::create(newString));
    }
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(StringJoinAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultStringJoinAddin);
        FIELD(strings, "String To Join");
        Addin::registerClassObjectMethod("STRING_JOIN",
                                         Addin::UTILITIES,
                                         "Constructs a handle to a string"
                                         " concatenated from supplied strings",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createString);
    }

    static IObject* defaultStringJoinAddin(){
        return new StringJoinAddin();
    }
    
};

CClassConstSP const StringJoinAddin::TYPE = CClass::registerClassLoadMethod(
    "StringJoinAddin", typeid(StringJoinAddin), load);

DRLIB_END_NAMESPACE
