
#include "edginc/config.hpp"
#include "edginc/Field.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/PublicObject.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/GDRDate.hpp"
#include "edginc/XObject.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Format.hpp"
#include <map>
#include <cstddef>
#include ext_hash_set
#include ext_hash_map
#include <algorithm>

/** Instances of the class Class represent classes and interfaces. */

DRLIB_BEGIN_NAMESPACE
struct type_info_lt{
    bool operator()(const type_info* i1, const type_info* i2) const{
#if defined(_MSC_VER)
        // non standard implementation of before method
        return ((i1->before(*i2)) != 0);
#else
        return (i1->before(*i2));
#endif
    }
};

struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};

/** How to cast from an object to an interface - use 'cast' and then recurse
    to next IfaceCast object until 'nextClass' is null. */
struct IfaceCast{
    TIfaceCastMethod* cast;
    const IfaceCast*  nextClass; // eg interface extends another interface
    IfaceCast(TIfaceCastMethod* cast, const IfaceCast* nextClass):
        cast(cast), nextClass(nextClass){}
    IfaceCast(){} // allow use in stl containers
};

typedef hash_map<string, CClassSP, MyStringHash> ClassHashTable;
typedef map<const type_info*, CClassSP, type_info_lt> ClassByTypeInfoHT;
typedef hash_map<string, CFieldSP, MyStringHash> FieldHashTable;
typedef hash_set<CClassConstSP, CClass::Hash> ClassHashSet;
typedef hash_map<CClassConstSP, IfaceCast, CClass::Hash> IfaceCastHash;
typedef hash_map<CClassConstSP, TAcceptMethod*, CClass::Hash> ActionHash;

////////////////////////////////////////////////////////////
// Support for DRI reflection //
////////////////////////////////////////////////////////////
class DRIType: public CObject{
public:
    static CClassConstSP const TYPE;
    ~DRIType(){}
    IObject* clone() const{
        // all derived classes must be immutable
        return const_cast<DRIType*>(this);
    }
protected:
    DRIType(CClassConstSP clazz): CObject(clazz){}
private:
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIType, clazz);
        SUPERCLASS(CObject);
        // no fields
    }
};

// put all the stl stuff in a separate class to keep stl headers out
// also used to minimise changes to the header file
// Class.hpp
class CClassImp{
public:
    static bool   typesBeingLoaded; /* if true we are busy trying to
                                       invoke the load method for each type */
    // hash table containing all known classes. This has to be a pointer
    // so that we can initialise it at our choice 
    static ClassHashTable    *hashTable;
    static ClassByTypeInfoHT *classByTypeInfoHT;
    static CClassVec*         classArray;
    static CClassVec*         gccBugWorkAround;

    bool   beingLoaded;   // class is being loaded
    bool   isNonVirtualBase; /* true: castToBase() is implemented for this
                                class ie probably CObject */
    bool   isEnum; // true if Class represents an Enum
    TRegisterClassLoadMethod  *loadMethod;   // how to load a class
    // how to create an empty shell for the object ie new CLASS()
    TCreateInstanceMethod     *createMethod;

    CClassConstSP   componentType; // if an array, the type of the components
    CClassConstSP   arrayType;  /* the type of the array which has components
                                   of this type */
    // how to create an empty array with this type as its component type
    TCreateArrayInstanceMethod *createArrayMethod;
    // set containing all interfaces that this class implements
    IfaceCastHash   impIfaces;
    // set containing all interfaces that this class implements directly
    ClassHashSet    parentIfaces;
    FieldHashTable  fieldsHash;  // fields as a map of vectors
    CFieldArray     fieldsVec;   // fields as an array
    ActionHash      actionHash;  // holds accept methods
    CClassVec       constructorTypes; /* types implementing IPublicObject
                                         which generate this type */
    CClassVec       assignableTypes; /* types assignable to instances
                                        of this class */

    CClassConstSP   driProxyType; /* if non null use this object for
                                           dri field reflection information
                                           (proxy types) or for parent
                                           types (IPublicObjects) */
    smartPtr<DRIType>     driType; // DRI equivalent object
    vector<CClass::EnumValue> enumValues;  // for enum, possible values

    static bool               regFailed; // true if failed to register classes
    static ModelException*    regError;  // what the error was
    CClassImp():beingLoaded(false), isNonVirtualBase(false), isEnum(false),
                loadMethod(0), createMethod(0),
                componentType(0), arrayType(0), createArrayMethod(0), 
                fieldsVec(0), driProxyType(0){}
    void createDRIType(CClassConstSP clazz); // builds DRI equivalent object
    smartPtr<DRIType> getDRIType(CClassConstSP clazz); /* gets and if needed
                                                          builds DRI equivalent
                                                          object */
};

bool CClassImp::typesBeingLoaded = false;
ClassHashTable* CClassImp::hashTable = 0;
ClassByTypeInfoHT* CClassImp::classByTypeInfoHT = 0;
CClassVec*         CClassImp::classArray = 0;
/* we hit a really strange bug with gcc 3.1 and 3.2 which meant that
   the static variables got initialised twice (so losing existing information.
   The symptoms were that various types got lost so certain test cases (eg
   AddinTest class) would not deserialise) */
CClassVec*         CClassImp::gccBugWorkAround = 0;

class ClassShutdown{
public:
    ~ClassShutdown(){
        for (size_t i = 0; i < CClassImp::classArray->size(); i++){
            delete (*CClassImp::classArray)[i];
        }
        delete CClassImp::hashTable;
        delete CClassImp::classByTypeInfoHT;
        delete CClassImp::classArray;
    }
};

// When classShutDown is destructed it frees all the memory
static ClassShutdown classShutDown;

class Array;
// static flag - use to determine error status. Used before main() is invoked
bool CClassImp::regFailed; // DO NOT INITIALISE HERE
ModelException* CClassImp::regError; // DO NOT INITIALISE HERE

CClass::EnumValue::EnumValue(int           valueAsInt, 
                             const string& valueAsString, 
                             const string& comment):
    valueAsInt(valueAsInt), valueAsString(valueAsString), comment(comment){}

/** A dummy version of registerClassLoadMethod for use by templates
    where we want to support multiple instances of the TYPE field
    across dlls.  (With MSVC it seems to be impossible to
    automatically export the TYPE field within a template). All this
    method does is to return the CClassConstSP object for the
    specified type_info reference. The idea is that the default
    template code initialises its TYPE field with a call to this
    function. Note that a call to registerClassLoadMethod must be made
    somewhere else (and once only) in order to provide the name and
    loadMethod. The reason for not putting the call to
    registerClassLoadMethod in the default template code is to a) avoid
    code bloat since the load method causes a lot of code to be created, and
    b) to avoid issues with specifying the name of the type */
CClassConstSP CClass::templateRegisterClass(
    const type_info&            classType){  // typeid of class
    try{
        // Handy for debugging more obscure area of type registration
        // cout << "templateRegisterClass: "+string(classType.name()) << endl;
        ClassByTypeInfoHT::const_iterator iter =
            CClassImp::classByTypeInfoHT->find(&classType);
        if (iter == CClassImp::classByTypeInfoHT->end()){
            // unfortunately it could be that this method is called before
            // the genuine call to registerClassLoadMethod for this type
            // So create CClass instance with empty name
            CClassSP newClass = create("", classType);
            return newClass;
        }
        return iter->second;
    } catch (exception& e){
        // flag error (in start up code at the moment so no exceptions)
        CClassImp::regFailed = true;
        if (CClassImp::regError){
            // concatenate those errors
            CClassImp::regError->addMsg(string(e.what()));
        } else {
            CClassImp::regError = new ModelException(e);
        }
        return 0;
    }
}

/** Supplied method will be invoked when the class type information
    is required. This will be typically be at start up. This solves
    the issue of the essentially random order the classes are
    initialised */
CClassConstSP CClass::registerClassLoadMethod(
    const char*                 className,
    const type_info&            typeInfo,  // typeid of class
    TRegisterClassLoadMethod*   loadMethod){
    try{
        CClassSP newClass = 0;
        if (CClassImp::classByTypeInfoHT){
            ClassByTypeInfoHT::iterator iter =
                CClassImp::classByTypeInfoHT->find(&typeInfo);
            if (iter != CClassImp::classByTypeInfoHT->end()){
                newClass = iter->second;
                // should only happen if done via a 'weak' initilisation as part
                // of a template. If it has a name, something's wrong
                if (!newClass->className.empty()){
                    throw ModelException("CClass::registerClassLoadMethod",
                                         "Attempt to register class "+
                                         string(typeInfo.name())+" twice");
                }
                // record the actual name
                newClass->className = className;
            }
        }
        // if the class doesn't exist yet then create it
        if (!newClass){
            newClass = create(className, typeInfo);
        }
        // save the load method
        newClass->data->loadMethod = loadMethod;
        return newClass;
    } catch (exception& e){
        // flag error (in start up code at the moment so no exceptions)
        CClassImp::regFailed = true;
        if (CClassImp::regError){
            // concatenate those errors
            CClassImp::regError->addMsg(string(e.what()));
        } else {
            CClassImp::regError = new ModelException(e);
        }
        return 0;
    }
}
/** Same as above but the 'atomic' classes where we want to register two
    type_info classes against the type eg 'double' and CDouble */
CClassConstSP CClass::registerClassWithAliasLoadMethod(
    const char*                 className,
    const type_info&            typeInfo,  // typeid of class
    const type_info&            typeInfoAlias,  // typeid of alias class
    TRegisterClassLoadMethod*   loadMethod){
    // Handy for debugging more obscure area of type registration
    //cout << "registerClassWithAliasLoadMethod: "+
    //    string(typeInfo.name()) << endl;
    //cout << "registerClassWithAliasLoadMethod: "+
    //    string(typeInfoAlias.name()) << endl;

    // use typeInfoAlias here since templateRegisterClass method is passed
    // the 'full' (ie QLib boxed or wrapped type). We need to find the 
    // relevant instance of CClass if templateRegisterClass() has already 
    // been called.
    CClassConstSP clazz = registerClassLoadMethod(className, 
                                                  typeInfoAlias, loadMethod);
    if (clazz){
        (*CClassImp::classByTypeInfoHT)[&typeInfo] = 
            const_cast<CClassSP>(clazz);
    }
    return clazz;
}
    
void CClass::registerCreateMethod(
    TCreateInstanceMethod*      createMethod){
    this->data->createMethod = createMethod;
}

/** Returns the modifiers for this class or interface, encoded in
    an integer. The Modifier class should be used to decode the
    modifiers. */
int CClass::getModifiers() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return modifiers;
}

/** Flags that this class is private. Private classes are, for 
    example, hidden from Analytics Service and the spreadsheet etc */
void CClass::setPrivate(){
     modifiers |= Modifier::PRIVATE; 
}

/** Flags that this class is private. Private classes are, for 
    example, hidden from Analytics Service and the spreadsheet etc */
void CClass::setPublic(){
     modifiers |= Modifier::PUBLIC;
}

/** Flags that this class is exported. Exported classes are public
    and clients may legitimately assume order of arguments (e.g.
    spreadsheet addins) */
void CClass::setExport() {
    setPublic();
    modifiers |= Modifier::EXPORT;
}    

/** for documentation */
void CClass::setDescription(const string& desc) {
    description = desc;
}

const string& CClass::getDescription() const {
    // if no description, at least return the class name
    if (!description.empty()) {
        return description;
    }
    // handle any proxy stuff
    CClassConstSP realClass = CObject::actualClass(this);
    return (realClass ? realClass->getDescription() : getName());
}

//// sets public, protected, private, abstract flags 
void CClass::setModifiers(){
    if (array){
        if (!data->componentType){
            // the component type is supposed to be set for all arrays
            throw ModelException("CClass::setModifiers", "Internal error");
        }
        if (!data->componentType->loaded){
            loadClass(data->componentType); // probably paranoia
        }
        modifiers = data->componentType->modifiers & 
            (Modifier::PUBLIC | Modifier::PROTECTED | Modifier::PRIVATE);
    } else if (!Modifier::isPublic(modifiers) &&
               !Modifier::isPrivate(modifiers) &&
               !Modifier::isProtected(modifiers)){
        if (!IPublicObject::TYPE->loaded){ // probably paranoia
            loadClass(IPublicObject::TYPE); 
        }
        if (!IPrivateObject::TYPE->loaded){ // probably paranoia
            loadClass(IPrivateObject::TYPE); 
        }
        // make default private
        if (IPrivateObject::TYPE->isAssignableFrom(this)){
            /* Protected classes are those that 
               implement IPrivateObject (which is a reasonable analogy). */
            modifiers |= Modifier::PROTECTED;
        } else if (IPublicObject::TYPE->isAssignableFrom(this)){
            modifiers |= Modifier::PUBLIC;
        } else
            modifiers |= Modifier::PRIVATE;
    }
    // crude way of determing whether abstract or not
    if ((!array && !data->createMethod) || 
        (array && !data->createArrayMethod)){
        modifiers |= Modifier::ABSTRACT; // make interfaces abstract (?)
    } else if (modifiers & Modifier::PUBLIC){
        // have create method and public so add to list of constructor classes
        data->constructorTypes.push_back(this);
    }
}

/** same as registerClassLoadMethod but for interfaces */
CClassConstSP CClass::registerInterfaceLoadMethod(
    const char*                 className,
    const type_info&            typeInfo,  // typeid of class
    TRegisterClassLoadMethod*   loadMethod)
{
    CClassConstSP iface = 
        registerClassLoadMethod(className, typeInfo, loadMethod);
    if (iface){
        const_cast<CClassSP>(iface)->isIface = true;
        const_cast<CClassSP>(iface)->modifiers |= Modifier::INTERFACE;
    }
    return iface;
}

/** load up class of given type */
void CClass::loadClass(const CClassConstSP& constClassToLoad){
    static const string method("CClass::loadClass");
    // if class is already loaded or is being loaded, then do nothing
    if (!constClassToLoad->loaded && !constClassToLoad->data->beingLoaded){
        // fun and games with whether the class is a const or not. Easier
        // to take a const (allows load on demand of any referenced class)
        CClassSP classToLoad = const_cast<CClassSP>(constClassToLoad);
        try{
            if (!classToLoad->isIface){ 
                if (!classToLoad->data->loadMethod){
                    throw ModelException(method, "NULL class load method");
                }
                // assignableTypes for interfaces are handled separately
                classToLoad->data->assignableTypes.push_back(classToLoad);
            }
            if (classToLoad->data->loadMethod){
                classToLoad->data->beingLoaded = true;
                classToLoad->data->loadMethod(classToLoad);
                // determine modifiers after class has loaded 
                classToLoad->loaded = true;
                classToLoad->data->beingLoaded = false;
                classToLoad->setModifiers();
            } else {
                classToLoad->loaded = true;
                // sort out superclass (this code previously in getSuperClass)
                if (!classToLoad->superClass && !classToLoad->isIface &&
                    classToLoad != IObject::TYPE){
                    // not sure if this is really necessary
                    classToLoad->superClass = IObject::TYPE;
                }
                classToLoad->setModifiers();
            }
        } catch (exception &e){
            classToLoad->loaded = false;
            const string& className = classToLoad->className.empty()?
                classToLoad->typeInfo->name(): classToLoad->className;
            throw ModelException(e, method, "Failed to load class "+className);
        }
    }
}
    
//// little utility to find a a key in a map based upon supplying the value 
static const type_info* findTypoInfo(CClassSP clazz){
    // find out what the type_info was the earlier class
    for (ClassByTypeInfoHT::const_iterator iter =
             CClassImp::classByTypeInfoHT->begin();
         iter !=  CClassImp::classByTypeInfoHT->end();
         ++iter){
        if (iter->second == clazz){
            return iter->first;
        }
    } 
    return 0;
}
/** Returns true if we are busy trying to invoke the load method for each
    type */
bool CClass::classLoadingInProgress(){
    return CClassImp::typesBeingLoaded;
}

/** load up all known classes */
void CClass::loadAllClasses(){
    CClassImp::typesBeingLoaded = true;
    if (CClassImp::regFailed){
        CClassImp::regError->addMsg("Registration of class load "
                                    "methods failed");
        throw (*CClassImp::regError);
    }
    // check the all the classes have set up correctly before we try
    // and load any
    for (ClassByTypeInfoHT::const_iterator iter =
             CClassImp::classByTypeInfoHT->begin();
         !(iter == CClassImp::classByTypeInfoHT->end());
         ++iter){
        if (iter->second->className.empty()){
            string m("No registerClassLoadMethod has been defined "
                     "for class "+string(iter->first->name())+"\nNote that all "
                     "specific instances of templated classes (eg Arrays) "
                     "must have their TYPE field explicitly defined");
            throw ModelException("CClass::loadAllClasses", m);
        }
    }
    if (!CClassImp::hashTable){
        // continuation of the work around for gcc3.2
        // Initialise and populate the hashtable here
        CClassImp::hashTable = new ClassHashTable();
        for (size_t i = 0; i < CClassImp::classArray->size(); i++){
            CClassSP clazz = (CClassSP) (*CClassImp::classArray)[i];
            const string& name = clazz->getName();
            ClassHashTable::value_type newEntry(name, clazz);
            pair<ClassHashTable::iterator, bool> entry = 
                CClassImp::hashTable->insert(newEntry);
            if (!entry.second){
                //entry already exists
                const type_info* type1 = findTypoInfo(entry.first->second);
                const type_info* type2 = findTypoInfo(clazz);
                throw ModelException("CClass::loadAllClasses", 
                                     "Class "+name+
                                     " ["+string(type1? 
                                                 type1->name(): "???")+"]"
                                     " already exists (when creating class "
                                     "with id "+
                                     string(type2->name())+")");
            }
        }
    }
    // now we can 'load' each class
    for (ClassHashTable::iterator myIterator = 
             CClassImp::hashTable->begin();
         !(myIterator == CClassImp::hashTable->end());
         ++myIterator){
        //const string& typeKey = myIterator->first;
        loadClass(myIterator->second);
    }
    // leave this flag as true if registration fails since we can view this
    // as having never finished
    CClassImp::typesBeingLoaded = false;

    // now initialise DRI stuff - cleaner to do once we know all the true
    // reflection stuff has been initialised.
    for (ClassHashTable::iterator myIterator = 
             CClassImp::hashTable->begin();
         !(myIterator == CClassImp::hashTable->end());
         ++myIterator){
        //const string& typeKey = myIterator->first;
        myIterator->second->getDRIType();
    }
}
    
CClass::CClass(const string &className, const type_info* typeInfo):
    CObject(TYPE),
    className(className), typeInfo(typeInfo),
    primitive(false), loaded(false), array(false),
    assignmentOperatorClones(false), // too dangerous to switch on by default
    superClass(0), modifiers(0), baseOffset(0), isIface(false) {
    data = new CClassImp();
}
    
/** flag that a class represents a native type */
void CClass::setNative(){
    primitive = true;
    assignmentOperatorClones = true;
}

/** create a new instance of a class - if not already initialised */
CClassSP CClass::create(const string& className,
                        const type_info& typeInfo){
    if (!CClassImp::gccBugWorkAround){
        CClassImp::regFailed = false; // initialise here
        CClassImp::regError = 0; // ditto
        /** The code here looks overly complicated and redundant but was put
            in to overcome a bug in gcc and/or the linker on solaris.
            See the comments above on CClassImp::gccBugWorkAround */
        if (!CClassImp::classByTypeInfoHT){
            CClassImp::classByTypeInfoHT = new ClassByTypeInfoHT();
        }
        if (!CClassImp::gccBugWorkAround){
            CClassImp::gccBugWorkAround = new CClassVec();
#if defined(DEBUG) || !defined(sun)
            // leaving this in causes problems on solaris opt
            delete CClassImp::gccBugWorkAround;
#endif
        }
        if (!CClassImp::classArray){
            CClassImp::classArray = new CClassVec();
        }
    }
    CClassSP newClass = new CClass(className, &typeInfo);
    (*CClassImp::classByTypeInfoHT)[&typeInfo] = newClass;
    CClassImp::classArray->push_back(newClass);
    return newClass;
}

/** Returns the Class object associated with the class or
    interface with the given string name. */
CClassConstSP CClass::forName(const string &className){
    ClassHashTable::const_iterator iter =
        CClassImp::hashTable->find(className);
    if (iter == CClassImp::hashTable->end()){
        throw ModelException("CClass::forName(string)",
                             "Class not found: " + className);
    }
    CClassSP theClass = iter->second;
    if (!theClass->loaded){
        loadClass(theClass); // load on demand
    }
    return theClass;
}

/** Returns the Class object associated with the class or
    interface with the given type_inf0 */
CClassConstSP CClass::forName(const type_info& classType){
    ClassByTypeInfoHT::const_iterator iter =
        CClassImp::classByTypeInfoHT->find(&classType);
    if (iter == CClassImp::classByTypeInfoHT->end()){
        string name = classType.name();
        throw ModelException("CClass::forName(type_info)",
                             "Class not found: " + name);
    }
    CClassSP theClass = iter->second;
    if (!theClass->loaded){
        loadClass(theClass); // load on demand
    }
    return theClass;
}
        
/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool CClass::isInstance(const IObjectConstSP& obj) const throw() {
    return isInstance(obj.get());
}

/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool CClass::isInstance(const IObject& obj) const throw() {
    return isInstance(&obj);
}

/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool CClass::isInstance(const IObject* obj) const throw() {
    if (!obj){
        return false;
    }
    return isAssignableFrom(obj->getClass());
}

/** Determines if the class or interface represented by this Class
    object is either the same as, or is a superclass or superinterface
    of, the class or interface represented by the specified Class
    parameter. */
bool CClass::isAssignableFrom(CClassConstSP clazz) const throw() {
    if (this == clazz){
        return true; // performance
    }
    if (!clazz){
        return false;
    }
    if (!loaded){
        loadClass(this); // load on demand
    }
    if (!clazz->loaded){
        loadClass(clazz); // load on demand
    }
    bool noMatch;
    if (!isIface){
        // we have already checked for the case 'this == clazz'
        do {
            clazz = clazz->superClass;
        } while ((noMatch = (this != clazz)) && clazz);
    } else {
        do{
            if (this == clazz){
                noMatch = false;
            } else {
                IfaceCastHash& theHash = clazz->data->impIfaces;
                IfaceCastHash::const_iterator iterator = theHash.find(this);
                noMatch = (iterator == theHash.end());
            }
        } while (noMatch && (clazz = clazz->superClass));
    }
    return !noMatch;
}

const string& CClass::getName() const throw(){
    return className;
}

const type_info& CClass::getTypeInfo() const {
    return *typeInfo;
}

/** Returns an array containing Class objects representing all
    the public classes and interfaces that are members of the
    class represented by this Class object. */
const CClassVec* CClass::getClasses() const{
    throw ModelException("CClass::getClasses", "not yet done");
    return 0;
}
    
/** Determines the interfaces implemented by the class or
    interface represented by this object. */
CClassVec CClass::getInterfaces() const{
    //CClassVec* ifaces = new CClassVec;
    CClassVec ifaces;
    IfaceCastHash& theHash = data->impIfaces;
    for (IfaceCastHash::const_iterator myIterator = theHash.begin();
         !(myIterator == theHash.end());
         ++myIterator){
        ifaces.push_back(myIterator->first);
    }
    return ifaces;
}

/** Returns the Class representing the superclass of the
    entity. If this Class represents either the IObject class, an
    interface, a primitive type, or void, then null is
    returned. If this object represents an array class then the
    Class object representing the Object class is returned.  */
CClassConstSP CClass::getSuperClass() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return superClass;
}

/* Determines if the specified Class object represents an
   interface type. */
bool CClass::isInterface() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return isIface;
}
    
/** Determines if the specified Class object represents a
    primitive type. */
bool CClass::isPrimitive() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return primitive;
}

/** Returns true if this type represented an enum */
bool CClass::isEnum() const{
    return data->isEnum;
}

/** Allow infrastructure to deduce value of doesAssignmentOperatorClone()
    by examining reflection information. This must be called before
    any fields are registered for this class. In particular, if 
    doesAssignmentOperatorClone() is true for the parent class and all
    fields are 'inline' and doesAssignmentOperatorClone() is true for the
    type of each field, then doesAssignmentOperatorClone() will be true for
    this class. */
void CClass::enableCloneOptimisations(){
    if (!superClass){
        throw ModelException("CClass::enableCloneOptimisations",
                             "Can't call enableCloneOptimisations before "
                             "superClass is set");
    }
    if ((superClass->assignmentOperatorClones || 
         superClass == IObject::TYPE) && 
        data->fieldsVec.empty() &&
        (!array || data->componentType->assignmentOperatorClones)){
        assignmentOperatorClones = true; /* NB value is updated in
                                            setDeclaredField() */
    }
}

/* does the copy constructor of this class implicity perform a
   clone Eg no pointers anywhere inc components such as DateTime or
   DateTimeArray etc */    
bool CClass::doesAssignmentOperatorClone() const{
    return assignmentOperatorClones;
}

/** Sets a Field object that reflects the specified
    declared field of the class represented by
    this Class object. */
void CClass::setDeclaredField(const CFieldSP& field)
{
    static const string method("CClass::setDeclaredField");
    const string& fieldName = field->getName();
    if (isInterface()){
        // by definition interfaces don't have fields
        throw ModelException(method, "Trying to register a field ("+fieldName+
                             ") against an interface ("+getName()+")");
    }
    // if non transient check the field is not repeated in parent classes
    if (!Modifier::isTransient(field->getModifiers())){
        CClassConstSP c = this;
        do {
            FieldHashTable::const_iterator iter = c->
                data->fieldsHash.find(fieldName);
            if (iter != c->data->fieldsHash.end() && 
                !Modifier::isTransient(iter->second->getModifiers())){
                throw ModelException(method,
                                     "Duplicate field "+fieldName+" in class "+ 
                                     className+"\nField is already defined in "+
                                     c->className);
            }
        } while ((c = c->getSuperClass()) != 0);
    }
    // store the field in the hash
    data->fieldsHash[fieldName] = field;
    // and also store it in our list - the order of the list is the order
    // in which the parameters are registered
    data->fieldsVec.push_back(field);
    
}
    
/** Specialised version of setDeclaredField for fields accessed by 
    pointers (either SP or regular pointers). Method exists just to try
    and reduce size of code generated by all the templates */
void CClass::setDeclaredPlainPointerField(
    const char*      fieldName,   
    const type_info& typeId,  // of field
    TFieldGetMethod* getMethod,
    TFieldSetMethod* setMethod,
    ptrdiff_t              offset) /* from start of structure */
{
    setDeclaredField(fieldName, typeId, getMethod, setMethod, 
                     CField::PLAIN_POINTER, offset);
}

/** Specialised version of setDeclaredField for fields accessed by 
    pointers (either SP or regular pointers). Method exists just to try
    and reduce size of code generated by all the templates */
void CClass::setDeclaredSmartPointerField(
    const char*      fieldName,   
    const type_info& typeId,  // of field
    TFieldGetMethod* getMethod,
    TFieldSetMethod* setMethod,
    ptrdiff_t              offset) /* from start of structure */
{
    setDeclaredField(fieldName, typeId, getMethod, setMethod, 
                     CField::SMART_POINTER, offset);
}

/** Specialised version of setDeclaredField for inLine fields. Method
    exists just to try and reduce size of code generated by all the
    templates */
void CClass::setDeclaredInLineField(const char*      fieldName,   
                                    const type_info& typeId,  // of field
                                    TFieldGetMethod* getMethod,
                                    TFieldSetMethod* setMethod,
                                    ptrdiff_t              offset)  /* from start of 
                                                                 structure */
{
    setDeclaredField(fieldName, typeId, getMethod, setMethod,
                     CField::INLINE, offset);
}

/** creates a field using the data given and stores on the 
    type given by typeName */
void CClass::setDeclaredField(
    const char*              fieldName,   
    const type_info&         typeId,  // of field
    TFieldGetMethod*         getMethod,
    TFieldSetMethod*         setMethod,
    CField::PointerAttribute pointerAttribute,
    ptrdiff_t                      offset) /* from start of structure */
{
    CFieldSP field = CField::create(this,
                                    fieldName,
                                    typeId,
                                    getMethod,
                                    setMethod,
                                    pointerAttribute,
                                    offset);
    if (assignmentOperatorClones && 
        (!field->getType()->assignmentOperatorClones ||
         field->getPointerAttribute() != CField::INLINE)){
        // by default the copy constructor for any class containing a pointer 
        // does not clone the pointer
        assignmentOperatorClones = false;
    }
    setDeclaredField(field);
}
        
/** sets the description for the given field of objects of this type */
void CClass::setFieldDescription(const char* fieldName,
                                 const char* description){
    CFieldSP field = const_cast<CFieldSP>(getDeclaredField(fieldName));
    field->setDescription(description);
}

/** As above but takes string */
void CClass::setFieldDescription(const char*   fieldName,
                                 const string& description){
    setFieldDescription(fieldName, description.c_str());
}

/** sets the transient flag for the field with name fieldName */
void CClass::setFieldTransient(const char* fieldName,
                               bool        isTransient){
    try{
        CFieldSP field = const_cast<CFieldSP>(getDeclaredField(fieldName));
        field->setTransient(isTransient);
    } catch (exception& e){
        // handy tip for developers
        throw ModelException(e, "CClass::setFieldTransient", "Has the "
                             "field name changed but FIELD_MAKE_TRANSIENT "
                             "not been updated?");
    }
}
/** sets the transientForIteration flag for the field
    with name fieldName */
void CClass::setFieldTransientForIteration(
    const char*   fieldName,
    bool          isTransientForIteration){
    try{
        CFieldSP field = const_cast<CFieldSP>(getDeclaredField(fieldName));
        field->setTransientForIteration(isTransientForIteration);
    } catch (exception& e){
        // handy tip for developers
        throw ModelException(e, "CClass::setFieldTransientForIteration",
                             "Has the "
                             "field name changed but FIELD_MAKE_TRANSIENT "
                             "not been updated?");
    }
}
    
    
/** sets the description for the given field of objects of this type */
const string& CClass::getFieldDescription(const string& fieldName){
    CFieldSP field = const_cast<CFieldSP>(getDeclaredField(fieldName));
    return field->getDescription();
}

/** sets the isOptional flag for the field with name fieldName */
void CClass::setFieldOptional(const char*   fieldName,
                              bool          isOptional){
    CFieldSP field = const_cast<CFieldSP>(getDeclaredField(fieldName));
    field->setOptional(isOptional);
}
        
/** Returns a Field object that reflects the specified
    declared field of the class represented by
    this Class object */
CFieldConstSP CClass::getDeclaredField(const string& fieldName) const
{
    if (!loaded){
        loadClass(this); // load on demand
    }
    CFieldConstSP field = hasDeclaredField(fieldName);
    if (!field){
        throw ModelException("CClass::getDeclaredField",
                             "No such field "+fieldName);
    }
    return field;
}
    
/** as per getDeclaredField but returns 0 if no field exists rather
    than throwing an exception */
CFieldConstSP CClass::hasDeclaredField(const string&  fieldName) const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    FieldHashTable::const_iterator iter = data->fieldsHash.find(fieldName);
    return (iter == data->fieldsHash.end()? 0: iter->second);
}
        
/** Returns an array of Field objects reflecting all the
    fields declared by the class represented by
    this Class object and which are part of the class c
    (so c can be this or a parent type). */
const CFieldArray& CClass::getDeclaredFields() const
{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return data->fieldsVec;
}

/** Create an array object of this  type and length. Note this is the
    type of the array, and not of the components */
IArray* CClass::newArrayInstance(int length) const{
    if (!array){
        throw ModelException("CClass::newArrayInstance",
                             "Class " + className + 
                             " is not an array type");
    }
    if (!data->createArrayMethod){
        throw ModelException("CClass::newArrayInstance",
                             "Can't create arrays "
                             "of type "+className);
    }
    if (length < 0){
        throw ModelException("CClass::newArrayInstance",
                             "Array length must be >= 0");
    }
    return data->createArrayMethod(length);
}

/** Create an array object of this component type and length. Note
    this is the type of the components and not the type of the array
    itself */
IArray* CClass::newArrayInstanceByComponent(int length) const{
    if (!data->arrayType){
        throw ModelException("CClass::newArrayInstanceByComponent",
                             "Can't create arrays with components "
                             "of type "+className);
    }
    return data->arrayType->newArrayInstance(length);
}


/** Returns the Class representing the component type of an array. */
CClassConstSP CClass::getComponentType() const{
    if (!array){
        throw ModelException("CClass::getComponentType",
                             "Class " + className + 
                             " is not an array type");
    }
    if (!data->componentType){
        throw ModelException("CClass::getComponentType",
                             "Component type for type "+className+
                             " not known!");
    }
    return data->componentType;
}

/** flag that a class represents an array type */
void CClass::setArrayType(const type_info&            componentType,
                          TCreateArrayInstanceMethod* createArrayMethod){
    static const string routine = "CClass::setArrayType";
    if (isIface){
        throw ModelException(routine, "Arrays cannot be interfaces. Use "
                             "registerClassLoadMethod not "
                             "registerInterfaceLoadMethod");
    }
    array = true;
    try{
        this->data->componentType = forName(componentType);
    } catch (exception& e){
        throw ModelException(e, routine, "Could not find component"
                             " type for array "+getName());
    }
    if (createArrayMethod &&
        className != data->componentType->className + "Array"){
        // the DRI spec has been forced upon us ...
        throw ModelException("CClass::setArrayType","Array for class "+
                             data->componentType->className+
                             " is "+className+"\nIt needs to be "+
                             data->componentType->className + "Array");
    }
    this->data->createArrayMethod = createArrayMethod;
    // only bother with setting component's array type if the class
    // can be instantiated - avoids problem of
    // what the component type for CArray - it is of course IObject
    // but ObjectArray has components of IObject too 
    if (createArrayMethod){
        CClassConstSP arrayType = this->data->componentType->data->arrayType;
        if (arrayType){
            throw ModelException(routine, "Class "+getName()+" already has an "
                                 "array class ("+
                                 arrayType->getName()+") defined for it");
        }
        /* tedious const problem - would need separate forName that
           didn't return a const CClass */
        CClassSP compType = const_cast<CClassSP>(this->data->componentType);
        compType->data->arrayType = this;
    }
}

/** Determines if this Class object represents an array class. */
bool CClass::isArray() const{
    return array;
}

/** Determines if this Class object represents a matrix class. */
bool CClass::isMatrix() const {
    return (className == "DoubleMatrix"); // v weak
    // we should have some CClass representation for matrices. Eg an
    // interface type that lives in the toolkit (NB DoubleMatrix doesn't).
    // However, no pressing need for it at the moment as only one matrix
    // type
    /** Do NOT do "return DRIMatrixType::TYPE->isInstance(getDRIType());"
        because this is a circular definition */
}


/** Invoked when this class is 'loaded' */
void CClass::load(CClassSP& classToLoad){
    REGISTER(CClass, classToLoad);
    SUPERCLASS(CObject);
    // need to put all code here for class class initialisation
    // initialise this class if not already
}

CClassConstSP const CClass::TYPE = CClass::registerClassLoadMethod(
    "Class", typeid(CClass), load);

/** set the superclass for types that have explicit implementation of
    castToBase() - typically only CObject */
void CClass::setSuperClass(const type_info&  superClassTypeId){
    superClass = forName(superClassTypeId);
    data->isNonVirtualBase = true;
    baseOffset = 0;
    // update list of assignable types
    CClassConstSP clazz = superClass;
    do{
        const_cast<CClassSP>(clazz)->data->assignableTypes.push_back(this);
    } while ((clazz = clazz->getSuperClass()));
}

/** set the superclass for types that are derived non virtually from
    a class that explicitly implements castToBase(). The offset is the
    offset of the child from this base class */
void CClass::setSuperClass(const type_info&  superClassTypeId,
                           ptrdiff_t               offset){
    setSuperClass(superClassTypeId); // bit lazy as this sets isNonVirtualBase
    if (superClass->isIface){
        throw ModelException("CClass::setSuperClass", "Parent class is an "
                             "interface. Superclass cannot be an interface");
    }
    // correct the value of isNonVirtualBase
    data->isNonVirtualBase = false;
    baseOffset = offset+superClass->baseOffset;
}
    
/** Flags that the supplied type should be used for DRI reflection
    information. In particular, for proxies (see CObject) this type is
    used for field reflection information whilst for IPublicObjects this
    type is used for parent information */
void CClass::setDRIProxyType(CClassConstSP proxyType){
    data->driProxyType = proxyType;
}

/** Flag that this type represents an enum. The default is false */
void CClass::setIsEnum(){
    data->isEnum = true;
}

/** Record a possible value for a Class that represents an enum */
void CClass::setEnumValue(int         enumValue,
                          const char* enumName,
                          const char* enumComment){
    data->enumValues.push_back(EnumValue(enumValue, enumName, enumComment));
}

/** Returns the first EnumValue whose valueAsString matches the given
    string. Only valid if this class represents an Enum. */
const CClass::EnumValue& CClass::getEnumValue(const string& enumAsString) const{
    static const string method("CClass::getEnumValue");
    if (!data->isEnum){
        throw ModelException(method, "Type '"+getName()+
                             "' does not represent an enum");
    }
    for (unsigned int i = 0; i < data->enumValues.size(); i++){
        if (data->enumValues[i].valueAsString == enumAsString){
            return data->enumValues[i];
        }
    }
    throw ModelException(method, "For enum of type "+getName()+
                         ", no enum found with name '"+enumAsString+"'");
}

/** Returns the EnumValue whose valueAsInt matches the given enum value (as
    an int). Only valid if this class represents an Enum. */
const CClass::EnumValue& CClass::getEnumValue(int enumAsInt) const{
    static const string method("CClass::getEnumValue");
    if (!data->isEnum){
        throw ModelException(method, getName()+" does not represent an enum");
    }
    for (unsigned int i = 0; i < data->enumValues.size(); i++){
        if (data->enumValues[i].valueAsInt == enumAsInt){
            return data->enumValues[i];
        }
    }
    throw ModelException(method, "For enum of type "+getName()+
                         ", no enum found with value "+
                         Format::toString(enumAsInt));
}

/** Returns the list of possible EnumValues for the type that this CClass
    represents. Returns null, if this CClass does not represent an enum.
    Do not delete the returned vector. */
const vector<CClass::EnumValue>* CClass::getEnumValues() const{
    return data->isEnum? &data->enumValues: NULL;
}

/** Creates a new instance of the class represented by this
    Class object. */
IObject* CClass::newInstance() const{
    if (!data->createMethod){
        throw ModelException("CClass::newInstance",
                             "Cannot instantiate class of type " + 
                             className);
    }
    return data->createMethod();
}

IObject* CClass::clone() const{
    throw ModelException("CClass::clone",
                         "Not supported - Clonable interface not"
                         " implemented");
    return 0; // keep MSVC happy
}

/** records the fact that this class implements given interface */
void CClass::setImplementedInterface(const type_info&  ifaceType,
                                     TIfaceCastMethod* toIface){
    CClassConstSP ifaceClass = forName(ifaceType);
    IfaceCast castData(toIface, 0);
    data->impIfaces[ifaceClass] = castData;
    data->parentIfaces.insert(ifaceClass);
    // now add all other interfaces that ifaceType extends to the list
    // of implemented interfaces
    IfaceCastHash& theHash = ifaceClass->data->impIfaces;
    for (IfaceCastHash::const_iterator myIterator = theHash.begin();
         !(myIterator == theHash.end());
         ++myIterator){
        IfaceCast castData(toIface, &myIterator->second);
        data->impIfaces[myIterator->first] = castData;
    }
}


/** returns the address of the start of an object (derived
    from IObject) given an IObjectSP */
void*  CClass::castFromIObject(const IObject* object){
#if defined(_MSC_VER)
    // dynamic_cast<void*> is relatively slow with MSVC
    if (!object){
        return 0;
    }
    CClassConstSP objClass = object->getClass();
    void* ptr = const_cast<void*>(objClass->staticCast(object));
    return ptr;
#else
    // dynamic_cast<void*> is v fast with gcc
    return dynamic_cast<void*>(const_cast<IObject*>(object));
#endif
}

/* Given an IObject*, do a dynamic_cast to object of the type 
   represented by this class. Only works for types that aren't
   interfaces */
const void* CClass::dynamicCast(const IObject* object) const{
    return dynamicCast(const_cast<IObject*>(object));
}

//// Declared in Object.hpp - see explanation there as to why this exists
const void* classDynamicCast(CClassConstSP clazz, const IObject* object){
    return clazz->dynamicCast(object);
}

/* Given an IObject*, do a dynamic_cast to object of the type 
   represented by this class. Only works for types that aren't
   interfaces */
void* CClass::dynamicCast(IObject* object) const{
    if (!object){
        return 0;
    }
    // casting to interfaces is implicitly dynamic
    if (isIface){
        return ifaceCast(object);
    }
    if (!isInstance(object)){
        throw ModelException("CClass::dynamicCast", "Cannot cast object of "
                             "type "+object->getClass()->getName()+" to type "+
                             getName());
    }
    return staticCast(object);
}

/* Given an IObject*, do a static_cast to object of the type 
   represented by this class. Only works for types that aren't
   interfaces */
const void* CClass::staticCast(const IObject* object) const{
    return staticCast(const_cast<IObject*>(object));
}

/* Given an IObject*, do a static_cast to object of the type 
   represented by this class. */
void* CClass::staticCast(IObject* object) const{
    // deal with null
    if (!object){
        return 0;
    }
    if (!loaded){
        loadClass(this); // load on demand
    }
    // and check for casting to interfaces
    if (isIface){
        return ifaceCast(object);
    }
#if defined(_MSC_VER)
    // dynamic_cast<void*> is relatively slow with msvc
    // first cast to non-virtual base
    char* nonVirtualBase = (char*) object->castToBase();
#else
    // first cast to actual object
    char* nonVirtualBase = (char*) dynamic_cast<void*>(object);
    // and then back to non-virtual base
    nonVirtualBase -= object->getClass()->baseOffset;
#endif
    // then add offset for this class
    nonVirtualBase += baseOffset;
    return nonVirtualBase;
}

/** cast to an interface - implicitly a 'dynamic' rather than a 'static' cast.
    'this' must be the iface class to cast to and this->isIface must be true */
void* CClass::ifaceCast(IObject* object) const{
    static const string method("CClass::ifaceCast");
    if (this == IObject::TYPE){
        return object;
    }
    CClassConstSP clazz = object->getClass();
    // First find parent class that implements interface
    do{
        IfaceCastHash& theHash = clazz->data->impIfaces;
        IfaceCastHash::const_iterator iter = theHash.find(this);
        if (iter != theHash.end()){
            // found it. Next cast from IObject to that object
            void* iface = clazz->staticCast(object);
            // then chain through IfaceCast methods
            const IfaceCast* castData = &iter->second;
            do {
                if (!castData->cast){
                    throw ModelException(method, "Registration problem - no "
                                         "cast method registered");
                }
                iface = castData->cast(iface);
            } while ((castData = castData->nextClass));
            return iface; // finished
        }
    } while ((clazz = clazz->superClass));
    throw ModelException(method, "Cannot cast object of type "+
                         object->getClass()->getName()+" to interface "+
                         getName());
}

/** Adds the supplied class to the list of types that can be used to
    build an instance of this class using the data dictionary route
    (ie supplied type implements the IPublicObject interface and the
    toPrivateObject() method returns an instance of this class */
void CClass::addConstructorClass(CClassConstSP clazz){
    if (!loaded){
        loadClass(this); // load on demand
    }
    if (!clazz->loaded){
        loadClass(clazz); // load on demand
    }
    
    if (this != clazz && !IPublicObject::TYPE->isAssignableFrom(clazz)){
        throw ModelException("CClass::addConstructorClass", clazz->className +
                             " does not implement IPublicObject");
    }
    data->constructorTypes.push_back(clazz);
}

/** Return an array of types that can be used to construct this class
    (derived classes are excluded) via the data dictionary route. */
const CClassVec& CClass::getConstructorClasses() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    return data->constructorTypes;
}

/** Returns an array of all clazzes which satisfy
    this.isAssignableFrom(clazz).  In other words, if clazz is an
    interface a list of all classes that implement this interface is
    returned wother with all interfaces that extend this
    interface. Oterhwise a list of all classes which are derived from
    this one (which includes supplied clazz) is returned. */
const CClassVec& CClass::getAssignableClasses() const{
    if (!loaded){
        loadClass(this); // load on demand
    }
    if (isIface && data->assignableTypes.empty()){
        // easier to do this by scanning all types
        const CClassVec& classes = allClasses();
        for (size_t i = 0; i < classes.size(); i++){
            if (isAssignableFrom(classes[i])){
                data->assignableTypes.push_back(classes[i]);
            }
        }
    }
    return data->assignableTypes;
}
    

/** Returns an array of types which can be used to construct this
    class, or one derived from it, via the data dictionary route. */
CClassVec CClass::listAllConstructorClasses() const{
    CClassVec constructors;
    // loop over all classes finding those derived from it
    const CClassVec& allTheClasses = allClasses();
    for (size_t i = 0; i < allTheClasses.size(); i++){
        if (isAssignableFrom(allTheClasses[i])){
            // add getConstructorClasses() to array
            const CClassVec& classes =
                allTheClasses[i]->getConstructorClasses();
            constructors.insert(constructors.end(), classes.begin(), 
                                classes.end());
        }
    }
    return constructors;
}

/** sets accept method for specific type of collector. */
void CClass::setAcceptMethod(
    const type_info& clazzID,
    const type_info& collectorID,
    TAcceptMethod*   acceptMethod)
{
    try{
        CClassConstSP clazz = forName(clazzID);
        CClassConstSP collectorClass = forName(collectorID);
        clazz->data->actionHash[collectorClass] = acceptMethod;
    } catch (exception& e){
        throw ModelException(&e, "CClass::setAcceptMethod");
    }
}

/** invokes accept method for specific type of collector. Returns true if
    an accept method was found */
bool CClass::invokeAcceptMethod(
    const IObject* instance,
    ICollector*    collector)
{
    CClassConstSP objClass;
    CClassConstSP collClass;
    if (!instance || !collector || 
        !(objClass  = instance->getClass()) ||
        !(collClass = collector->getClass()))
    {
        throw ModelException("CClass::invokeAcceptMethod", 
                             "Null inputs");
    }

    bool noMatch = true;
    ActionHash::const_iterator iterator;
    while ( noMatch && objClass ){
        iterator = objClass->data->actionHash.find(collClass);
        noMatch  = (iterator == objClass->data->actionHash.end());
        if ( noMatch ) {
            objClass = objClass->getSuperClass();
        }
    }
    if (!noMatch){
        // need to cast both objects to correct type
        void* objPtr = castFromIObject(instance);
        void* collPtr = castFromIObject(collector);
        (iterator->second)(objPtr, collPtr);
    }
    return (!noMatch);
}
   

/** Return an array of all known classes (including interfaces) */   
const CClassVec& CClass::allClasses() {
    return *CClassImp::classArray;
}

////////////////////////////////////////////////////////////
// Support for DRI reflection //
////////////////////////////////////////////////////////////
CClassConstSP const DRIType::TYPE = 
CClass::registerClassLoadMethod("DRIType", typeid(DRIType), load);
typedef array<smartPtr<DRIType>, DRIType> DRITypeArray;
DEFINE_TEMPLATE_TYPE(DRITypeArray);

class DRIScalarType: public DRIType{
    int typeId; // DRIScalarType__typeId: DR_TYPE constant
public:
    static CClassConstSP const TYPE;
    DRIScalarType(int typeId): DRIType(TYPE), typeId(typeId){}
private:
    static IObject* defaultConstructor(){
        return new DRIScalarType(0);
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIScalarType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(typeId, "DR_TYPE constant");
    }
};

CClassConstSP const DRIScalarType::TYPE = 
CClass::registerClassLoadMethod("DRIScalarType", typeid(DRIScalarType), load);

class DRIArrayType: public DRIType{
    smartPtr<DRIType> elementType; // DRIArrayType_elementType
public:
    static CClassConstSP const TYPE;
    DRIArrayType(smartPtr<DRIType> elementType): 
        DRIType(TYPE), elementType(elementType) {}
private:
    static IObject* defaultConstructor(){
        return new DRIArrayType(smartPtr<DRIType>());
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIArrayType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(elementType, "Type of elements in array");
    }
};

CClassConstSP const DRIArrayType::TYPE = 
CClass::registerClassLoadMethod("DRIArrayType", typeid(DRIArrayType), load);

class DRIMatrixType: public DRIType{
public:
    static CClassConstSP const TYPE;
    DRIMatrixType(): DRIType(TYPE) {}
private:
    static IObject* defaultConstructor(){
        return new DRIMatrixType();
    }
    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIMatrixType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
    }
};

CClassConstSP const DRIMatrixType::TYPE = 
CClass::registerClassLoadMethod("DRIMatrixType", typeid(DRIMatrixType), load);

//// For representing Enums in DRI
class DRIEnumerationType: public DRIType{
    string                            typeName;
    string                            description;
    StringArray                       identifiers;
    StringArray                       identifierDescriptions;
public:
    static CClassConstSP const TYPE;
    DRIEnumerationType(const string&                        typeName,
                       const string&                        description,
                       const vector<CClass::EnumValue>&     enumValues):
        DRIType(TYPE), typeName(typeName), description(description),
        identifiers(enumValues.size()),
        identifierDescriptions(enumValues.size())
    {
        for (unsigned int i = 0; i < enumValues.size(); i++){
            identifiers[i] = enumValues[i].valueAsString;
            identifierDescriptions[i] = enumValues[i].comment;
        }
    }
private:
    DRIEnumerationType(): DRIType(TYPE){}
    static IObject* defaultConstructor(){
        return new DRIEnumerationType();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIEnumerationType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(typeName, "Name of type");
        FIELD(description, "Description of type");
        FIELD(identifiers, "Possible values for Enums")
        FIELD(identifierDescriptions, "Description of each Enum")
    }
};

CClassConstSP const DRIEnumerationType::TYPE = 
CClass::registerClassLoadMethod(
    "DRIEnumerationType", typeid(DRIEnumerationType), load);

class DRIAbstractType: public DRIType{
    string                            typeName;
    array<smartPtr<DRIType>, DRIType> baseTypes;
    string                            description;
public:
    struct SortFunc;
    friend struct SortFunc;
    static CClassConstSP const TYPE;
    DRIAbstractType(CClassConstSP                            clazz,
                    const string&                            typeName,
                    const array<smartPtr<DRIType>, DRIType>& baseTypes,
                    const string&                            description):
        DRIType(clazz), typeName(typeName), baseTypes(baseTypes), 
        description(description){}

    DRIAbstractType(const string&                            typeName,
                    const array<smartPtr<DRIType>, DRIType>& baseTypes,
                    const string&                            description):
        DRIType(TYPE), typeName(typeName), baseTypes(baseTypes), 
        description(description){}

    /** Provide support for sorting types based upon name */
    struct SortFunc{
        bool operator()(smartPtr<DRIType> p1, smartPtr<DRIType> p2){
            if (DRIAbstractType::TYPE->isInstance(*p1) &&
                DRIAbstractType::TYPE->isInstance(*p2)){
                DRIAbstractType& abstract1 = static_cast<DRIAbstractType&>(*p1);
                DRIAbstractType& abstract2 = static_cast<DRIAbstractType&>(*p2);
                return (abstract1.typeName < abstract2.typeName);
            }
            return false;
        }
    };
protected:                    
    DRIAbstractType(CClassConstSP clazz): DRIType(clazz){}
private:
    static IObject* defaultConstructor(){
        return new DRIAbstractType(TYPE);
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIAbstractType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(typeName, "Name of type");
        FIELD(baseTypes, "Direct base classes (one level up)");
        FIELD(description, "Description of type");
    }
};

CClassConstSP const DRIAbstractType::TYPE = 
CClass::registerClassLoadMethod("DRIAbstractType", typeid(DRIAbstractType),
                                load);

class DRIArgumentDescriptor: public CObject{
    string            name;
    smartPtr<DRIType> type;
    int               groupIndex;         // set to 0
    int               alternationIndex;   // set to 0
    int               indexInAlternation; // set to 0
    bool              isOptional;
    string            description;

public:
    static CClassConstSP const TYPE;
    DRIArgumentDescriptor(CFieldConstSP field): 
        CObject(TYPE), name(field->getName()), 
        groupIndex(0), alternationIndex(0), indexInAlternation(0), 
        isOptional(field->isOptional()), description(field->getDescription()){
        CClassConstSP clazz = field->getType();
        type = clazz->data->getDRIType(clazz); // get hold of DRIType
    }
    IObject* clone() const{
        // all derived classes must be immutable
        return const_cast<DRIArgumentDescriptor*>(this);
    }
private:
    DRIArgumentDescriptor(): CObject(TYPE),
        groupIndex(0), alternationIndex(0), indexInAlternation(0) {}

    static IObject* defaultConstructor(){
        return new DRIArgumentDescriptor();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIArgumentDescriptor, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "Name of argument");
        FIELD(type, "Type of argument");
        FIELD(groupIndex, "Always 0");
        FIELD(alternationIndex, "Always 0");
        FIELD(indexInAlternation, "Always 0");
        FIELD(isOptional, "Whether parameter is optional");
        FIELD(description, "Description of argument");
    }
};

CClassConstSP const DRIArgumentDescriptor::TYPE = 
CClass::registerClassLoadMethod("DRIArgumentDescriptor", 
                                typeid(DRIArgumentDescriptor),
                                load);

typedef array<smartPtr<DRIArgumentDescriptor>, DRIArgumentDescriptor> DRIArgumentDescriptorArray;

DEFINE_TEMPLATE_TYPE(DRIArgumentDescriptorArray);

class DRIConcreteType: public DRIAbstractType{
    array<smartPtr<DRIArgumentDescriptor>, DRIArgumentDescriptor> members;
    bool                                                          isExported;
public:
    static CClassConstSP const TYPE;
    DRIConcreteType(
        const string&                            typeName,
        const array<smartPtr<DRIType>, DRIType>& baseTypes,
        const string&                            description,
        array<smartPtr<DRIArgumentDescriptor>,DRIArgumentDescriptor> members,
        bool                                                        isExported):
        DRIAbstractType(TYPE, typeName, baseTypes, description),
        members(members), isExported(isExported){}
    
protected:
    DRIConcreteType(
        CClassConstSP                            clazz,
        const string&                            typeName,
        const array<smartPtr<DRIType>, DRIType>& baseTypes,
        const string&                            description,
        array<smartPtr<DRIArgumentDescriptor>,DRIArgumentDescriptor> members,
        bool                                                        isExported):
        DRIAbstractType(clazz, typeName, baseTypes, description),
        members(members), isExported(isExported){}
    
protected:
    DRIConcreteType(CClassConstSP clazz): DRIAbstractType(clazz){}
    static IObject* defaultConstructor(){
        return new DRIConcreteType(TYPE);
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIConcreteType, clazz);
        SUPERCLASS(DRIAbstractType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(members, "Members of class");
        FIELD(isExported, "Whether is class exported");
    }
};

CClassConstSP const DRIConcreteType::TYPE = 
CClass::registerClassLoadMethod("DRIConcreteType", typeid(DRIConcreteType),
                                load);

class DRIExecutableType: public DRIConcreteType{
    smartPtr<DRIType> resultType; 
public:
    static CClassConstSP const TYPE;
    DRIExecutableType(
        const string&                            typeName,
        const array<smartPtr<DRIType>, DRIType>& baseTypes,
        const string&                            description,
        array<smartPtr<DRIArgumentDescriptor>,DRIArgumentDescriptor> members,
        bool                                                        isExported,
        smartPtr<DRIType>                                           resultType):
        DRIConcreteType(TYPE, typeName, baseTypes, description,
                        members, isExported), resultType(resultType) {}
private:
    DRIExecutableType(): DRIConcreteType(TYPE){}
    static IObject* defaultConstructor(){
        return new DRIExecutableType();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIExecutableType, clazz);
        SUPERCLASS(DRIConcreteType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(resultType, "Type of result");
    }
};

CClassConstSP const DRIExecutableType::TYPE = 
CClass::registerClassLoadMethod("DRIExecutableType", 
                                typeid(DRIExecutableType), load);
#if 0
class DRIEnumerationType: public DRIType{
    string      typeName;
    string      description;
    StringArray identifiers;
    StringArray identifierDescriptions;
public:
    static CClassConstSP const TYPE;
    DRIEnumerationType(const string& typeName,
                       const string& description,
                       const StringArray& identifiers,
                       const StringArray& identifierDescriptions):
        DRIType(TYPE), typeName(typeName), description(description), 
        identifierDescriptions(identifierDescriptions){}
private:
    DRIEnumerationType(): DRIType(TYPE){}
    static IObject* defaultConstructor(){
        return new DRIEnumerationType();
    }

    static void load(CClassSP& clazz){
        clazz->setPublic();
        REGISTER(DRIEnumerationType, clazz);
        SUPERCLASS(DRIType);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(typeName, "Name of enumerated type");
        FIELD(description, "Description of enumerated type");
        FIELD(identifiers, "Range of allowed identifiers");
        FIELD(identifierDescriptions, 
                     "Descriptions of allowed identifiers");
    }
};

CClassConstSP const DRIEnumerationType::TYPE = 
CClass::registerClassLoadMethod("DRIEnumerationType", 
                                typeid(DRIEnumerationType), load);
#endif

CClass::~CClass(){
    delete data;
}


//// builds DRI equivalent object
void CClassImp::createDRIType(CClassConstSP clazz){
    static const string method("CClassImp::createDRIType");
    if (clazz->isPrimitive()){
        if (clazz->isEnum()){
            driType = smartPtr<DRIType>(
                new DRIEnumerationType(clazz->getName(),
                                       clazz->getDescription(),
                                       *clazz->getEnumValues()));
        } else {
            int typeId;
            if (clazz == CDouble::TYPE){
                typeId = DR_DOUBLE;
            } else if (clazz == CInt::TYPE){
                typeId = DR_INT;
            } else if (clazz == CBool::TYPE){
                typeId = DR_BOOL;
            } else if (clazz == CString::TYPE){
                typeId = DR_STRING;
            } else if (clazz == GDRDate::TYPE){
                typeId = DR_DATE;
            } else {
                throw ModelException(method, "Internal error");
            }
            driType = smartPtr<DRIType>(new DRIScalarType(typeId));
        }
    } else if (clazz->isMatrix()){
        driType = smartPtr<DRIType>(new DRIMatrixType());
    } else if (XObject::TYPE->isAssignableFrom(clazz)){
        // types derived from XObject come through the interface as DR_OBJECTs
        driType = smartPtr<DRIType>(new DRIScalarType(DR_OBJECT));
    } else if (IDRObject::TYPE->isAssignableFrom(clazz)){
        // IDRObject corresponds (in a sense) to a Variant (which is what
        // DR_UNDEFINED matches to)
        driType = smartPtr<DRIType>(new DRIScalarType(DR_UNDEFINED));
    } else if (clazz->isArray() && clazz->data->createArrayMethod != 0){
        // we check createArrayMethod above as don't want to mark abstract
        // classes as arrays
        CClassConstSP eltClazz = clazz->getComponentType();
        smartPtr<DRIType> eltType(eltClazz->data->getDRIType(eltClazz));
        driType = smartPtr<DRIType>(new DRIArrayType(eltType));
    } else {
        CClassConstSP classForFields = clazz; // default
        CClassConstSP classForParents = clazz; // default
        // the procoess of identifying proxy types etc  needs tidying up
        if (driProxyType){
            // either a IPublicObject object or a proxy
            // want to hide this from DRI clients
            classForFields = CObject::proxyClass(clazz);
            if (!classForFields){
                // not a proxy
                classForFields = clazz; // use this class for data members
                classForParents = driProxyType; // use this for parents
                if (!IPublicObject::TYPE->isAssignableFrom(clazz)){
                    throw ModelException(method,
                                         "setDRIProxyType can only be used "
                                         "for proxies or for IPublicObjects");
                }
            } else if (classForFields != driProxyType){
                throw ModelException(method, "setDRIProxyType must store the "
                                     "proxy type");
            }
        }
        // Get all immediate parents - do interfaces first
        array<smartPtr<DRIType>, DRIType> parents;
        ClassHashSet& theSet = classForParents->data->parentIfaces;
        for (ClassHashSet::const_iterator myIterator = theSet.begin();
             !(myIterator == theSet.end());
             ++myIterator){
            CClassConstSP ifaceClazz = (*myIterator);
            parents.push_back(ifaceClazz->data->getDRIType(ifaceClazz));
        }
        // for consistency across platforms, sort the interfaces
        sort(parents.begin(), parents.end(), DRIAbstractType::SortFunc());
        // and then always put superclass as the first entry
        if (!classForParents->isInterface()){
            CClassConstSP superClass = classForParents->getSuperClass();
            if (superClass){
                parents.insert(parents.begin(), 
                               superClass->data->getDRIType(superClass));
            }
        }
        int modifiers = classForFields->getModifiers();
        const string& classname = clazz->getName();
        const string& desc = classForFields->getDescription();
        if (Modifier::isInterface(modifiers) || 
            Modifier::isAbstract(modifiers)){
            driType = smartPtr<DRIType>(
                new DRIAbstractType(classname, parents, desc));
        } else {
            CFieldArray fields(Addin::getDataClassFields(classForFields));
            array<smartPtr<DRIArgumentDescriptor>, DRIArgumentDescriptor> 
                members(fields.size());
            for (int i = 0; i < members.size(); i++){
                members[i] = smartPtr<DRIArgumentDescriptor>(
                    new DRIArgumentDescriptor(fields[i]));
            }
            if (ClientRunnable::TYPE->isAssignableFrom(classForParents)){
                // all client runnable methods return an IObject
                smartPtr<DRIType> returnType(
                    IObject::TYPE->data->getDRIType(IObject::TYPE));
                driType = smartPtr<DRIType>(
                    new DRIExecutableType(classname, parents,
                                          desc,
                                          members, 
                                          Modifier::isExport(modifiers),
                                          returnType));
            } else {
                driType = smartPtr<DRIType>(
                    new DRIConcreteType(classname, parents,
                                        desc, members,
                                        Modifier::isExport(modifiers)));
            }
        }
    }
}

/* gets and if needed builds DRIType for this class */
smartPtr<DRIType> CClassImp::getDRIType(CClassConstSP clazz){
    static const string method("CClassImp::getDRIType");
    if (typesBeingLoaded){
        throw ModelException(method, "DRI type information not available until"
                             " all classes have been loaded");
    }

    if (!driType){
        createDRIType(clazz);
        if (!driType){
            throw ModelException(method, "Internal error");
        }
    }
    return driType;
}

             
/** Return the DRIType for this CClass */
IObjectSP CClass::getDRIType() const{
    return data->getDRIType(this);
}


DRLIB_END_NAMESPACE
