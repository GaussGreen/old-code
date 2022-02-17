//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XClass.cpp
//
//   Description : Describes classes from other DRI libraries
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2004
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif
#include "edginc/XClass.hpp"
#include "edginc/XField.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/DRAnalyticsInterface.h"
#include "edginc/Format.hpp"
#include "edginc/XMap.hpp"
#include "edginc/XArray.hpp"
#include "edginc/DRLibrary.hpp"
#include ext_hash_set
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
struct MyStringHash {
    size_t operator()(const string& str) const {
        return hash_string(str);
    }
};
#if 0
// not needed currently
struct MyDRServiceHash{
    size_t operator()(const DRService* svc) const {return (size_t)svc;}
};
#endif
/*
Originally these typedefs were used rather than the class definitions
below. However, the solaris assembler has problems with these long names.
By having these classes we can work around the problem.
typedef hash_map<string, XClassSP, MyStringHash> ClassHashTable;
typedef hash_map<string, smartPtr<XField>, MyStringHash> FieldHashTable;
typedef hash_set<XClassConstSP, XClass::Hash> ClassHashSet;
*/
class ClassHashTable: public hash_map<string, XClassSP, MyStringHash>{};
class FieldHashTable: public hash_map<string, smartPtr<XField>, MyStringHash>{};
class ClassHashSet: public hash_set<XClassConstSP, XClass::Hash>{};

typedef hash_map<string, ClassHashTable, MyStringHash> ClassHashTableBySvc;
typedef hash_map<string, XClassVec, MyStringHash> XClassVecBySvc;


// put all the stl stuff in a separate class to keep stl headers out
// also used to minimise changes to the header file Class.hpp
class XClassImp{
public:
    // the [artificial] 'root' object across all services (DR_OBJECT type)
    static XClassSP               rootObjectType;
    /* the [artificial] 'root' object including variants across all
     * services (DR_UNDEFINED type) */
    static XClassSP               rootVariantType;
    // hash table containing all known classes by service. 
    static ClassHashTableBySvc    allClasses;
    static XClassVecBySvc         classArray;
#if 0
    DRService*      drService; // the service in which this class lives
#endif
    string          svcName;
    bool            loaded;
    bool            beingLoaded;
    string          description;
    XClassConstSP   componentType; // if an array, the type of the components
#if 0
    XClassConstSP   arrayType;  /* the type of the array which has components
                                   of this type */
#endif
    // set containing all parents that this class implements
    ClassHashSet    allParents;
    // set containing all immediate parents of this class
    ClassHashSet    parents;
    XClassVec       assignableTypes; /* types assignable to instances
                                        of this class */

    FieldHashTable  fieldsHash;  // fields as a map of vectors
    XFieldArray     fieldsVec;   // fields as an array
    //////////// methods /////////
    XClassImp(const string& svcName): 
        svcName(svcName), loaded(false), beingLoaded(false),
        componentType(0) {}
    static const ClassHashTable& loadAllClasses(DRService*    svc,
                                                const string& svcName);
    void populateAllParents(ClassHashSet& allParents, XClassSP clazz);
};
ClassHashTableBySvc XClassImp::allClasses;
XClassVecBySvc      XClassImp::classArray;
XClassSP            XClassImp::rootVariantType;
XClassSP            XClassImp::rootObjectType;

const string& XClass::getDescription() const {
    return data->description;
}

/** Simple constructor sets primitive and array to false */
XClass::XClass(const string& svcName, const string &className): 
    CObject(TYPE), className(className), primitive(false), 
    array(false), data(new XClassImp(svcName)){}

XClassSP XClass::create(const string& svcName, const string& name,
                        bool       isArray, bool isPrimitive,
                        bool       loaded){
    // look it up
    XClassSP clazz = forNameInternal(svcName, name);
    // if it doesn't exist, create it
    if (!clazz){
        clazz = new XClass(svcName, name);
        clazz->array = isArray;
        clazz->primitive = isPrimitive;
        clazz->data->loaded = loaded;
        // and store it
        XClassImp::allClasses[svcName][name] = clazz;
        XClassImp::classArray[svcName].push_back(clazz);
        // and then update XClassImp::rootVariantType and rootObjectType
        if (!isPrimitive && XClassImp::rootObjectType){
            XClassImp::rootObjectType->data->assignableTypes.push_back(clazz);
        }
        if (XClassImp::rootVariantType){
            XClassImp::rootVariantType->data->assignableTypes.push_back(clazz);
        }
    }
    return clazz;
}

XClassSP XClass::createPrimitiveType(const string& svcName, int type){
    static const string method("XClass::createPrimitiveType");
    static const char* primNames[] = {"", "Bool", "Int", "Double", "String", 
                                      "", "GDRDate"};
    if (type > DR_MAX_TYPE){
        throw ModelException(method, "Invalid DR_TYPE ("+
                             Format::toString(type)+") for primitive type");
    }
    if (type == DR_OBJECT){
        return XClassImp::rootObjectType;
    }
    if (type == DR_UNDEFINED){
        return XClassImp::rootVariantType;
    }
    // decide upon the name
    string name(primNames[type]);
    // create it
    return create(svcName, name, false, true, true);
}

XClassSP XClass::createArrayType(const string& svcName, XClassSP eltClazz){
    // decide upon the name
    string name = eltClazz->getName() + "Array";
    // create it
    XClassSP clazz = create(svcName, name, true, false, true);
    if (!clazz->data->componentType){
        clazz->data->componentType = eltClazz;
    }
    return clazz;
}

XClassSP XClass::createMatrixType(const string& svcName){
    // decide upon the name
    string name("DoubleMatrix"); // to sort out
    // create it
    return create(svcName, name, false, true /* well sort of. to do: review */,
                  true);
}

void XClass::setDeclaredField(const string& svcName, XMapSP fieldDesc){
    string fieldName(fieldDesc->getString(DRIArgumentDescriptor_name));
    try{
        XObjectSP type(fieldDesc->getXObject(DRIArgumentDescriptor_type));
        XClassSP fieldType = load(svcName, type);
        string description(fieldDesc->getString(
                               DRIArgumentDescriptor_description));
        bool isOptional(fieldDesc->getBool(DRIArgumentDescriptor_isOptional));
        smartPtr<XField> field(XField::create(this, fieldName, fieldType,
                                              description, isOptional));
        if (data->fieldsHash.find(fieldName) != data->fieldsHash.end()){
            throw ModelException("XClass::setDeclaredField",
                                 "Duplicate entry for field with name "+
                                 fieldName+" in class "+className);
        }
        // store the field in the hash
        data->fieldsHash[fieldName] = field;
        // and also store it in our list - the order of the list is the order
        // in which the parameters are registered
        data->fieldsVec.push_back(field.get());
    } catch (exception& e){
        throw ModelException(e, "XClass::setDeclaredField", "Failed for field "+
                             fieldName);
    }
}


XClassSP XClass::load(const string& svcName, XObjectSP clazzDesc){
    static const string method("XClass::load");
    // get hold of type of type descriptor
    const string& descName = clazzDesc->getXClassName();
    string typeName(descName); // for error reporting
    try {
        // turn it into a map
        XMapSP mapClassDesc(clazzDesc->toMap());
        if (descName == DRIScalarType__name){
            // pull out the DR_TYPE
            int typeId = mapClassDesc->getInt(DRIScalarType_typeId);
            return createPrimitiveType(svcName, typeId);
        } else if (descName == DRIArrayType__name){
            // pull out the DR_TYPE
            XObjectSP eltType(mapClassDesc->
                              getXObject(DRIArrayType_elementType));
            // get the corresponding XClass
            XClassSP eltClazz = load(svcName, eltType);
            return createArrayType(svcName, eltClazz);
        } else if (descName == DRIMatrixType__name){
            return createMatrixType(svcName);
        } else if (descName == DRIEnumerationType__name){
            // to do ... pretend it's an int for now?
            return createPrimitiveType(svcName, DR_INT);
        } else if (descName == DRIAbstractType__name ||
                   descName == DRIConcreteType__name ||
                   descName == DRIExecutableType__name){
            typeName = mapClassDesc->getString(DRIAbstractType_typeName);
            XClassSP clazz = create(svcName, typeName, false, false, false);
            // if class is already loaded or is being loaded, then do nothing
            if (!clazz->data->loaded && !clazz->data->beingLoaded){
                clazz->data->beingLoaded = true;
                XArraySP parents(mapClassDesc->
                                 getXArray(DRIAbstractType_baseTypes));
                int numParents = parents->size();
                for (int i = 0; i < numParents; i++){
                    IDRObjectSP elt(parents->getElt(i));
                    XObjectSP xObj(XObjectSP::dynamicCast(elt));
                    clazz->data->parents.insert(load(svcName, xObj));
                }
                string desc(mapClassDesc->
                            getString(DRIAbstractType_description));
                clazz->data->description = desc;
                if (descName != DRIAbstractType__name){
                    XArraySP fields(mapClassDesc->
                                    getXArray(DRIConcreteType_members));
                    int numFields = fields->size();
                    for (int i = 0; i < numFields; i++){
                        IDRObjectSP theElt(fields->getElt(i));
                        XObject* elt = DYNAMIC_CAST(XObject, theElt.get());
                        XMapSP fieldDesc(elt->toMap());
                        clazz->setDeclaredField(svcName, fieldDesc);
                    }
                }
                clazz->data->beingLoaded = false;
                clazz->data->loaded = true;
            }                
            return clazz;
        } else {
            throw ModelException(method, "Unrecognised type "+descName);
        }
    } catch (exception& e){
        throw ModelException(e, method, 
                             "Failed loading class of type "+typeName);
    }
}

/** Populate allParents with clazz and all its parents */
void XClassImp::populateAllParents(ClassHashSet& allParents, XClassSP clazz){
    try{
        for (ClassHashSet::iterator iter = clazz->data->parents.begin();
             iter != clazz->data->parents.end(); ++iter){
            allParents.insert(*iter);
            populateAllParents(allParents, const_cast<XClassSP>(*iter));
        }
    } catch (exception& e){
        throw ModelException(e, "XClass::populateAllParents",
                             "Failed for class "+clazz->className);
    }
}

void XClass::populateAllParents(const string& svcName){
    XClassVec& classes = XClassImp::classArray[svcName];
    for (unsigned int i = 0; i < classes.size(); i++){
        XClassSP clazz = const_cast<XClassSP>(classes[i]);
        ClassHashSet& allParents = clazz->data->allParents;
        clazz->data->populateAllParents(allParents, clazz);
    }
}

void XClass::loadAllClasses(DRService* svc){
    // look up DRLibrary first
    DRLibrarySP lib(DRLibrary::getLibrary(svc));
    // and then use service name
    XClassImp::loadAllClasses(svc, lib->getServiceName());
}

/** Load up all the classes from the specified DRI library (if not done
    already). If svc null, use DRLibrary to retrieve a svc pointer */
const ClassHashTable& XClassImp::loadAllClasses(DRService*    svc,
                                                const string& svcName){
    static const string method("XClass::loadAllClasses");
    try {
        ClassHashTableBySvc::const_iterator iter = 
            XClassImp::allClasses.find(svcName);
        if (iter == XClassImp::allClasses.end()){
            if (!svc){
                // look up svc from DRLibrary
                svc = DRLibrary::getLibrary(svcName)->getService();
            }
            XArraySP allTypes(XObject::typeList(svc));
            // loop over types
            for (int i = 0; i < allTypes->size(); i++){
                IDRObjectSP theElt(allTypes->getElt(i));
                XObjectSP clazzDesc(XObjectSP::dynamicCast(theElt));
                XClass::load(svcName, clazzDesc);
            }
            // now sort out allParents
            XClass::populateAllParents(svcName);
            iter = XClassImp::allClasses.find(svcName);
            if (iter == XClassImp::allClasses.end()){
                throw ModelException(method, "Internal error");
            }
        }
        return iter->second;
    } catch (exception& e){
        throw ModelException(e, method);
    }
}

/** Returns null if class doesn't exist yet else the class */
XClassSP XClass::forNameInternal(const string& svcName,
                                 const string &className){
    ClassHashTableBySvc::const_iterator iter =
        XClassImp::allClasses.find(svcName);
    if (iter == XClassImp::allClasses.end()){
        return 0;
    }
    ClassHashTable::const_iterator iter2 = iter->second.find(className);
    if (iter2 == iter->second.end()){
        return 0;
    }
    XClassSP theClass = iter2->second;
    return theClass;
}
    
/** Returns the Class object associated with the class
    with the given string name and the corresponding DRService. */
XClassConstSP XClass::forName(const string& svcName, const string &className){
    const ClassHashTable& classes = XClassImp::loadAllClasses(0, svcName);
    ClassHashTable::const_iterator iter = classes.find(className);
    if (iter == classes.end()){
        throw ModelException("XClass::forName",
                             "Class '"+className+"' not found in "+svcName);
    }
    XClassSP theClass = iter->second;
    return theClass;
}


/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool XClass::isInstance(IObjectConstSP obj) const throw() {
    return isInstance(obj.get());
}

/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool XClass::isInstance(const IObject& obj) const throw() {
    return isInstance(&obj);
}

/** Determines if the specified Object is assignment-compatible
    with the object represented by this Class. */
bool XClass::isInstance(const IObject* obj) const throw() {
    if (!obj){
        return false;
    }
    if (!XObject::TYPE->isInstance(obj)){
        return false;
    }
    if (XArray::TYPE->isInstance(obj)){
        if (!array){
            // obj is array but this XClass doesn't correspond to one
            return false;
        }
        // work around for magnet - it says the type of all arrays is a 
        // variant array. So go round in circles looking at the component type
        const XArray* xArray = STATIC_CAST(XArray, obj);
        XClassConstSP eltClass = xArray->elementClass();
        if (!eltClass){
            throw ModelException("XClass::isInstance", "DRI function returned "
                                 "null for arrayElementType. Upgrade your DRI "
                                 "implementation!");
        }
        // then just do isAssignableFrom on the component type of both arrays
        return data->componentType->isAssignableFrom(eltClass);
    } else {
        const XObject* xObj = STATIC_CAST(XObject, obj);
        XClassConstSP xClass = xObj->getXClass();
        return isAssignableFrom(xClass);
    }
}

/** Determines if the class or interface represented by this Class
    object is either the same as, or is a superclass or superinterface
    of, the class or interface represented by the specified Class
    parameter. */
bool XClass::isAssignableFrom(XClassConstSP clazz) const throw() {
    if (this == clazz){
        return true;
    }
    const ClassHashSet& theHash = clazz->data->allParents;
    ClassHashSet::const_iterator iterator = theHash.find(this);
    return (iterator != theHash.end());
}

const string& XClass::getName() const throw(){
    return className;
}

/** Determines if the specified Class object represents a
    primitive type. */
bool XClass::isPrimitive() const{
    return primitive;
}

    
/** Returns a Field object that reflects the specified
    declared field of the class represented by
    this Class object */
XFieldConstSP XClass::getDeclaredField(const string&fieldName) const
{
    XFieldConstSP field = hasDeclaredField(fieldName);
    if (!field){
        throw ModelException("XClass::getDeclaredField",
                             "No such field "+fieldName);
    }
    return field;
}
    
/** as per getDeclaredField but returns 0 if no field exists rather
    than throwing an exception */
XFieldConstSP XClass::hasDeclaredField(const string&  fieldName) const{
    FieldHashTable::const_iterator iter = data->fieldsHash.find(fieldName);
    return (iter == data->fieldsHash.end()? 0: iter->second.get());
}
        
/** Returns an array of Field objects reflecting all the
    fields declared by the class represented by
    this Class object and which are part of the class c
    (so c can be this or a parent type). */
const XFieldArray& XClass::getDeclaredFields() const{
    return data->fieldsVec;
}

#if 0
/** Create an array object of this  type and length. Note this is the
    type of the array, and not of the components */
IArray* XClass::newArrayInstance(int length) const{
    if (!array){
        throw ModelException("XClass::newArrayInstance",
                             "Class " + className + 
                             " is not an array type");
    }
#if 0
    if (!data->createArrayMethod){
        throw ModelException("XClass::newArrayInstance",
                             "Can't create arrays "
                             "of type "+className);
    }
    if (length < 0){
        throw ModelException("XClass::newArrayInstance",
                             "Array length must be >= 0");
    }
    return data->createArrayMethod(length);
#endif
}

/** Create an array object of this component type and length. Note
    this is the type of the components and not the type of the array
    itself */
IArray* XClass::newArrayInstanceByComponent(int length) const{
    if (!data->arrayType){
        throw ModelException("XClass::newArrayInstanceByComponent",
                             "Can't create arrays with components "
                             "of type "+className);
    }
    return data->arrayType->newArrayInstance(length);
}
#endif

/** Returns the Class representing the component type of an array. */
XClassConstSP XClass::getComponentType() const{
    if (!array){
        throw ModelException("XClass::getComponentType",
                             "Class " + className + 
                             " is not an array type");
    }
    if (!data->componentType){
        throw ModelException("XClass::getComponentType",
                             "Component type for type "+className+
                             " not known!");
    }
    return data->componentType;
}

/** Determines if this Class object represents an array class. */
bool XClass::isArray() const{
    return array;
}

/** Invoked when this class is 'loaded' */
void XClass::myLoad(CClassSP& classToLoad){
    REGISTER(XClass, classToLoad);
    SUPERCLASS(CObject);
    // need to put all code here for class class initialisation
    // initialise this class if not already
    // Now create artificial types - leave service name blank. 
    // Names are arbitrary.
    XClassImp::rootVariantType = create("", "GDRObject", false, false, true);
    XClassImp::rootObjectType = create("", "XObject", false, false, true);
}

CClassConstSP const XClass::TYPE = CClass::registerClassLoadMethod(
    "XClass", typeid(XClass), myLoad);

#if 0
/** Creates a new instance of the class represented by this
    Class object. */
IObject* XClass::newInstance() const{
#if 0
    if (!data->createMethod){
        throw ModelException("XClass::newInstance",
                             "Cannot instantiate class of type " + 
                             className);
    }
    return data->createMethod();
#endif
}
#endif

IObject* XClass::clone() const{
    return const_cast<XClass*>(this);
}

/** Returns an array of all clazzes which satisfy
    this.isAssignableFrom(clazz).  In other words, if clazz is an
    interface a list of all classes that implement this interface is
    returned together with all interfaces that extend this
    interface. Otherwise a list of all classes which are derived from
    this one (which includes supplied clazz) is returned. */
const XClassVec& XClass::getAssignableClasses() const{
    if (data->assignableTypes.empty()){
        // easier to do this by scanning all types
        const XClassVec& classes = allClasses(data->svcName);
        for (unsigned int i = 0; i < classes.size(); i++){
            if (isAssignableFrom(classes[i])){
                data->assignableTypes.push_back(classes[i]);
            }
        }
    }
    return data->assignableTypes;
}
    

/** Return an array of all known classes for the specified DRService */   
const XClassVec& XClass::allClasses(const string& svcName) {
    XClassImp::loadAllClasses(0, svcName);
    return XClassImp::classArray[svcName];
}


XClass::~XClass(){
    delete data;
}


DRLIB_END_NAMESPACE
