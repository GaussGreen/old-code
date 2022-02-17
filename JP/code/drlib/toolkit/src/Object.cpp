#include "edginc/config.hpp"
#define QLIB_OBJECT_CPP
#if defined(_MSC_VER)
// disable warning truncated decorated names information
#pragma warning(disable : 4503)
#endif
#include "edginc/Null.hpp"
#include "edginc/PublicObject.hpp"
#include "edginc/PrivateObject.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/Format.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/XArray.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE

string IObject::stringOf(IObjectConstSP o) {
    return !o ? string("NULL") : o->toString();
}

/** Creates an exception reporting that a clone failed for the supplied
    class because the input object was null */
ModelException IObject::makeExceptionForCopy(
    CClassConstSP    clazz){        // what type was being copied
    static const string routine = "copy";
    const string& name = clazz->getName();
    return ModelException(routine, "Object of type "+name+" is NULL");
}

static void makeClassPublic(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
}
    
/* static final value in IObject class */
CClassConstSP const IObject::TYPE = 
CClass::registerInterfaceLoadMethod("IObject", typeid(IObject), 
                                    makeClassPublic);

// IObject::IObject(){} is in the header file

/** Write a very brief English name or description of an object to a stream */
ostream& operator <<(ostream& o, const IObject& x) {
    o << x.toString();
    return o;
}

/* static final value in IPrivateObject class */
CClassConstSP const IPrivateObject::TYPE = 
CClass::registerInterfaceLoadMethod("IPrivateObject", 
                                    typeid(IPrivateObject), 0);

/* static final value in IPublicObject class */
CClassConstSP const IPublicObject::TYPE = 
CClass::registerInterfaceLoadMethod("IPublicObject", typeid(IPublicObject), 
                                    makeClassPublic);

/* static final value in ICollector class */
CClassConstSP const ICollector::TYPE = 
CClass::registerInterfaceLoadMethod("ICollector", typeid(ICollector), 0);

/* static final value in ITypeConvert class */
CClassConstSP const ITypeConvert::TYPE = 
CClass::registerInterfaceLoadMethod("ITypeConvert", typeid(ITypeConvert), 0);
ITypeConvert::ITypeConvert(){}
ITypeConvert::~ITypeConvert(){}
ITypeConvert::ITypeConvert(const ITypeConvert& rhs){}
ITypeConvert& ITypeConvert::operator=(const ITypeConvert& rhs){return *this;}

const string CObject::OBJECT_TYPE_ATTRIBUTE = "TYPE";

/** Create a deep copy an object by iterating over its
    components and copying those. */
IObject* CObject::clone() const{
    IObjectConstSP obj(IObjectConstSP::attachToRef(this)); 
    try {
        CClassConstSP  c = getClass();
        IObjectSP clonedObj(c->newInstance());
        // iterate over a object and all its superclasses, extracting
        // fields in each level of inheritance and clone them
        do {
            const CFieldArray& fields = c->getDeclaredFields();
            // copy all fields including transients
            for (unsigned int i = 0; i < fields.size(); i++) {
                // clone field from from obj to clonedObj
                fields[i]->copyComponent(obj, clonedObj);
            }
        } while ((c = c->getSuperClass()) != 0);
        return clonedObj.release();
    }
    catch (exception& e) {
        throw ModelException(&e, "CObject::clone", "couldn't clone " + 
                             obj->getClass()->getName());
    }
}


// CObject::CObject(const CClassConstSP& objClass) throw() is in the header
// for performance reasons

/** default - do nothing */
void  CObject::validatePop2Object(){ /* empty */}

/** default implementation - does nothing */
void CObject::fieldsUpdated(const CFieldArray& fields){}

/** Invoked when this class is 'loaded' */
void CObject::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    clazz->setSuperClass(typeid(IObject)); // should review this
    clazz->enableCloneOptimisations(); /* this is to allow derived classes to
                                          switch this on */
}


/** write object out to writer */
void CObject::write(const string& tag, Writer* writer) const {
    static const string routine = "CObject::write";
    try {
        // start tag and see if/what object we need to write out
        // approx: write <tag TYPE='objectType'> (does convertToPublicRep)
        IObjectConstSP obj(writer->objectStart(tag, "", this, true));
        if (obj.get()){
            // object has not been written out before - so write out
            CClassConstSP c = obj->getClass();
            do {
                const CFieldArray& fields = c->getDeclaredFields();
                for (unsigned int i = 0; i < fields.size(); i++) {
                    // skip transient fields
                    if (!Modifier::isTransient(fields[i]->getModifiers())){
                        const string&  id = fields[i]->getName();
                        IObjectConstSP o(fields[i]->constGet(obj));
                        if (o.get()) {
                            o->write(id, writer);
                        } else {
                            writer->writeNull(id);
                        }
                    }
                }
            } while ((c = c->getSuperClass()) != 0);
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}
  
// read in an object from reader 
IObject* CObject::read(Reader::Node* elem, Reader* reader) {
    static const string method("CObject::read");
    CClassConstSP clazz = 0;
    try {
        clazz = elem->getClass();

        IObjectSP object;
        if (!elem->isNull()) {
            if (clazz->isArray()) {
                //// not possible to instantiate XArrays without service ptr
                if (XArray::TYPE == clazz){
                    object = IObjectSP(new XArray(0)); // work around
                } else {
                    int length = elem->arrayLength();
                    object = IObjectSP(clazz->newArrayInstance(length));
                }
            } else {
                object = IObjectSP(clazz->newInstance());
            }
            object->import(elem, reader);
            // now see if object's public interface is different to
            // its internal one
            object = convertToPrivateRep(object);
            // another hack for XObject and the fact that serialisation of
            // objects from different libraries is a bit broken in DRI
            if (object->getClass() == XObject::WRAPPER_TYPE){
                // convert to XObject (or derived type)
                checkType(object, XObject::TYPE);
            }
        }  else {
            object = CNull::create();
        }
        return object.release();
    }
    catch (exception& e) {
        throw ModelException(e, method, string("couldn't build ") + 
                             (clazz?
                              clazz->getName(): "unknown type")+ " for node "+
                             elem->name());
    }    
}

/** populate an empty object from reader */
void CObject::import(Reader::Node* elem, Reader* reader) {
    static const string method("CObject::import");
    try {
        CDataDictionarySP dd(CDataDictionary::create(getClass()));
        dd->copyOnPop2Obj(false); /* no need to clone components as they
                                     are freshly created */
        string fieldname;
        try {
            Reader::NodeListSP nl(elem->children());       
            for (unsigned int i = 0; i < nl->size(); i++) {
                Reader::Node* child = (*nl)[i].get();
                fieldname = child->name();
                IObjectSP obj(reader->read(child));
                if (dd->contains(fieldname)){
                    throw ModelException(method,
                                         "Duplicate entries for field "+fieldname);
                }
                dd->put(fieldname, obj);
            }
        }
        catch (exception& e) {
            throw ModelException(e, method,
                                 "Failed whilst processing: "+fieldname);
        }    

        dd->populateObject(this);
    }
    catch (exception& e) {
        throw ModelException(e, method,
                             "Failed for type: " + getClass()->getName());
    }    
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void CObject::outputWrite(const string& linePrefix,
                          const string& prefix, ostream& stream) const{
    IObjectConstSP obj(IObjectConstSP::attachToRef(this)); 
    try {
        // convert to public object if needed
        obj = convertToPublicRep(obj);

        CClassConstSP  c = obj->getClass();
        // iterate over a object and all its superclasses, extracting
        // fields in each level of inheritance and invoking write method
        int totalPrinted = 0;
        do {
            const CFieldArray& fields = c->getDeclaredFields();
            for (unsigned int i = 0; i < fields.size(); i++)
            {
                // skip transient fields
                if (!Modifier::isTransient(fields[i]->getModifiers())){
                    IObjectConstSP o = fields[i]->constGet(obj);
                    if (o.get()){
                        // recurse
                        const string& fieldName = fields[i]->getName();
                        o->outputWrite(linePrefix, prefix == ""? fieldName:
                                       prefix+"_"+fieldName, 
                                       stream);
                        totalPrinted++;
                    }
                }
            }
        } while ((c = c->getSuperClass()) != 0);
        if (totalPrinted == 0){
            // print out type of object instead
            stream << linePrefix << prefix << "_TYPE: " 
                   << obj->getClass()->getName() << endl;
        }
    }
    catch (exception& e) {
        throw ModelException(e, "CObject::outputWrite", "Failed for object"
            " of type "+obj->getClass()->getName());
    }
}

/** Returns a hash code value for the object by combining the hashCode for
    all fields making this object.
    Note that this method might be relatively
    slow due to the recursive nature of the implementation */
int CObject::hashCode() const{
    IObjectConstSP obj(IObjectConstSP::attachToRef(this)); 
    CClassConstSP  c = obj->getClass();
    int hCode = (size_t) c; // good starting point
    // iterate over a object and all its superclasses, extracting
    // fields in each level of inheritance and invoking hashCode method
    do {
        const CFieldArray& fields = c->getDeclaredFields();
        for (unsigned int i = 0; i < fields.size(); i++) {
            IObjectConstSP o(fields[i]->constGet(obj));
            if (o.get()){
                // recurse
                hCode ^= o->hashCode();
            }
        }
    } while ((c = c->getSuperClass()) != 0);
    return hCode;
}

/** Indicates whether some other object is "equal to" this one by comparing all
    fields making up the pair of objects (assuming they are of the same type) */
bool CObject::equalTo(const IObject* obj) const{
    if (this == obj){
        return true; // the obvious comparison
    }
    CClassConstSP  c;
    if (!obj || (c = obj->getClass()) != getClass()){
        return false; // if obj null or a different type
    }
    IObjectConstSP obj1(IObjectConstSP::attachToRef(this)); 
    IObjectConstSP obj2(IObjectConstSP::attachToRef(obj)); 
    // iterate over a object and all its superclasses, extracting
    // fields in each level of inheritance and invoking equalTo method
    do {
        const CFieldArray& fields = c->getDeclaredFields();
        for (unsigned int i = 0; i < fields.size(); i++) {
            IObjectConstSP o1(fields[i]->constGet(obj1));
            IObjectConstSP o2(fields[i]->constGet(obj2));
            if ((!o1 && o2.get() != 0) || (!o2 && o1.get())){
                return false;
            }
            // recurse
            if (o1.get() && !o1->equalTo(o2.get())){
                return false;
            }
        }
    } while ((c = c->getSuperClass()) != 0);
    return true;
}

// suppose a type wants an object but is given an array
// here we allow an "object from array" method to be registered
// with a double lookup - first find the thing we're given
// and then see if it can be converted to what we want

typedef hash_map<CClassConstSP, CObject::TObjectFromArray*, 
    CClass::Hash> ObjectFromArrayLookup;

typedef hash_map<CClassConstSP, ObjectFromArrayLookup, 
    CClass::Hash> ObjectFromArrayHash;
 
static ObjectFromArrayHash convertMethods;

/** Register a conversion method from an array to a given type */
void CObject::registerObjectFromArrayMethod(CClassConstSP     givenClass,
                                            CClassConstSP     targetClass,
                                            TObjectFromArray* method) {
    static const string routine = "CObject::registerObjectFromArrayMethod";
    try {
        // see if thisClass has an entry in convertMethods
        ObjectFromArrayHash::iterator iter = convertMethods.find(givenClass);
        if (iter == convertMethods.end()) {
            // if it's not there, create an empty ObjectFromArrayLookup and
            // stuff the method in there
            ObjectFromArrayLookup lookup;
            lookup[targetClass] = method;
            convertMethods[givenClass] = lookup;
        }
        else {
            iter->second[targetClass] = method;
        }
    }
    catch (exception& e) {
        throw ModelException(e, routine);
    }
}

/** Support "proxy" interfaces
    internally use objClass, publicly use proxyClass as if it was objClass
    i.e. publicly build class A, but with class B (the proxy) interface */
struct ProxyConversion {
    ProxyConversion(CClassConstSP clazz, CObject::ProxyConvertMethod* method) :
        clazz(clazz), method(method) {}

    ProxyConversion(): clazz(0), method(0) {}
    
    CClassConstSP                clazz;
    CObject::ProxyConvertMethod* method;
};

typedef hash_map<CClassConstSP, ProxyConversion, CClass::Hash> ProxyLookup;

// lazy - have two storage containers, one for each direction
static ProxyLookup toProxyMethods;
static ProxyLookup fromProxyMethods;

/** register a proxy for a class and how to transform in both directions */
void CObject::registerObjectProxy(CClassConstSP       objClass, 
                                  CClassConstSP       proxyClass,
                                  ProxyConvertMethod* toProxy,
                                  ProxyConvertMethod* fromProxy) {
    static const string method("CObject::registerObjectProxy");
    try {
        // see if objClass has an entry in toProxyMethods
        ProxyLookup::iterator iter = toProxyMethods.find(objClass);
        if (iter == toProxyMethods.end()) {
            // if it's not there, create an empty ProxyConversion and
            // stuff the method in there, and do the other way round too
            ProxyConversion proxy(proxyClass, toProxy);
            toProxyMethods[objClass] = proxy;

            ProxyConversion object(objClass, fromProxy);
            fromProxyMethods[proxyClass] = object;
        }
        else {
            iter->second.clazz  = proxyClass;
            iter->second.method = toProxy;   

            // and the other way
            ProxyLookup::iterator iter = fromProxyMethods.find(proxyClass);
            // this really ought to be there
            if (iter == fromProxyMethods.end()) {
                throw ModelException(method, "from proxy conversion missing "
                                     "for class " + proxyClass->getName());
            }
            iter->second.clazz  = objClass;
            iter->second.method = fromProxy;               
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** given a class, see if it has a proxy, or what it is a proxy for */
CClassConstSP CObject::proxyClass(CClassConstSP objClass) {
    ProxyLookup::iterator iter = toProxyMethods.find(objClass);
    if (iter == toProxyMethods.end()) {
        return 0;
    }
    else {
        return iter->second.clazz;
    }
}

CClassConstSP CObject::actualClass(CClassConstSP proxyClass) {
    ProxyLookup::iterator iter = fromProxyMethods.find(proxyClass);
    if (iter == fromProxyMethods.end()) {
        return 0;
    }
    else {
        return iter->second.clazz;
    }
}

/** if an object has a proxy interface, convert to the proxy */
IObjectSP CObject::toProxy(const IObjectSP& object) {
    static const string method("CObject::toProxy");
    try {
        ProxyLookup::iterator iter = toProxyMethods.find(object->getClass());
        if (iter == toProxyMethods.end()) {
            return object;
        }
        else {
            return iter->second.method(object);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** if an object is a proxy, convert to the "real" object */
IObjectSP CObject::fromProxy(const IObjectSP& object) {
    static const string method("CObject::fromProxy");
    try {
        ProxyLookup::iterator iter = fromProxyMethods.find(object->getClass());
        if (iter == fromProxyMethods.end()) {
            return object;
        }
        else {
            return iter->second.method(object);
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** little utility method - converts supplied array to array of supplied type.
    Assumes that desiredType is an array and that the components of the array
    can convert themselves or we want a weaker (ie some sort of parent) type */
static IArraySP convertArray(IArray& theArray, CClassConstSP desiredType){
    // create an array and convert each component
    int arrayLen = theArray.getLength();
    CClassConstSP cmptClazz = desiredType->getComponentType();
    IArraySP newArray(desiredType->newArrayInstance(arrayLen));
    for (int i = 0; i < arrayLen; i++){
        try{
            // pull out component from array
            IObjectSP elem(theArray.get(i));
            // convert it (first cast to relevant interface type)
            if (ITypeConvert::TYPE->isInstance(elem)){
                ITypeConvert* convertableObj = 
                    DYNAMIC_CAST(ITypeConvert, elem.get());
                convertableObj->convert(elem, cmptClazz);
            }
            // and set the object into the new array
            newArray->set(i, elem);
        } catch (exception& e){
            throw ModelException(e, "CObject::convertArray",
                                 "Failed whilst processing"
                                 " element "+Format::toString(i)+" of"
                                 " array");
        }
    }
    return newArray;
}

static IObjectSP convertArrayToObject(IObjectSP     object,
                                      CClassConstSP desiredType){
    // we have an array and we want an object
    // see if given array has an array conversion methods
    ObjectFromArrayHash::iterator iter =
        convertMethods.find(object->getClass());
    if (iter != convertMethods.end()) {
        // see if the desired type is there, if not see if there's
        // an appropriate superclass registered instead
        ObjectFromArrayLookup::iterator iter2 = iter->second.find(desiredType);
        if (iter2 != iter->second.end()) {
            return iter2->second(object, desiredType);
        } else {
            // Note we support desiredType being an interface which means that
            // we can't just iterate through the superclasses of desiredType.
            // There are only a handful of methods registered anyway
            for (iter2 = iter->second.begin(); iter2 != iter->second.end();
                 ++iter2) {
                if (desiredType->isAssignableFrom(iter2->first)) {
                    return iter2->second(object,desiredType);
                }                           
            }
        }
    }
    // give up (caller gives detailed error message)
    throw ModelException("CObject::convertArrayToObject",
                         "No conversions found");
}

/** Note static method. Takes supplied object and checks to see if it
    is of the desired type. If not it will try and convert the object
    to the desired type by using the ITypeConvert interface provided
    the supplied object implements it or by converting to a private/public
    representation. If both the supplied object and
    the desired type is an array and the components of the supplied
    array implement the ITypeConvert interface then the array is
    converted element by element. Else if the object is an array, the
    table of TObjectFromArray methods will be consulted to see if
    there is a suitable method tto convert the array into an
    object. The return value reflects whether the object was converted
    or not (true = object was converted). If it is not possible to
    convert the object to the required type an exception is thrown */
bool CObject::checkType(IObjectSP&           object,
                        CClassConstSP        desiredType)
{
    static const string routine("CObject::checkType");
    // first see if there is anything to do
    if (!desiredType->isInstance(object)){
        CClassConstSP clazz = object->getClass();
        try{
            // does the object implement the TypeConvert interface?
            if (ITypeConvert::TYPE->isInstance(object)){
                ITypeConvert& convertableObj = 
                    dynamic_cast<ITypeConvert&>(*object);
                // and then convert
                convertableObj.convert(object, desiredType);
            } else if ((desiredType->isArray() &&  // if we want an array
                        clazz->isArray()) &&        // and we have an array
                       (ITypeConvert::TYPE->isAssignableFrom(  // and components
                           clazz->getComponentType())||// can convert themselves
                        desiredType->getComponentType()->isAssignableFrom(
                            clazz->getComponentType()))) // or want weaker type
            { 
                IArray&  array = dynamic_cast<IArray&>(*object);
                object = convertArray(array, desiredType);
            } 
            else if (clazz->isArray() && !desiredType->isArray()) {
                // we have an array and we want an object
                object = convertArrayToObject(object, desiredType);
            } else if (IPublicObject::TYPE->isInstance(object)){
                // try converting to a private object
                object = convertToPrivateRep(object);
                // then see if we have to any more conversions
                checkType(object, desiredType);
            } else if (IPrivateObject::TYPE->isInstance(object)){
                // try converting to a public object
                object = convertToPublicRep(object);
                // NB Do not put a call to checkType here because you might end
                // up with an infinite loop.
            }
            // double check we've succeeded
            if (!desiredType->isInstance(object)){
                throw ModelException(routine, "No conversions found");
            }
            return true;
        } catch (exception& e){
            throw ModelException(e, routine, "Can't convert Object of "
                                 "type "+clazz->getName()+
                                 " to object of type "+desiredType->getName());
        }
    }
    return false;
}

CClassConstSP const CObject::TYPE = CClass::registerClassLoadMethod(
    "Object", typeid(CObject), load);

/** method for processing collectors. - should we make this ICollectorSP?
    Looks up method in object's class */
bool CObject::accept(ICollector* collector) const{
    return CClass::invokeAcceptMethod(this, collector);
}

/** Utility method to convert an object to its private representation
    if it has one. The supplied smart pointer is modified in place */
IObjectConstSP CObject::convertToPrivateRep(const IObjectConstSP& obj){
    if (IPublicObject::TYPE->isInstance(obj)){
        const IPublicObject& pubObj = dynamic_cast<const IPublicObject&>(*obj);
        return (IObjectConstSP(pubObj.toPrivateObject()));
    }
    return obj;
}

/** Utility method to convert an object to its private representation
    if it has one. The supplied smart pointer is modified in place */
IObjectSP CObject::convertToPrivateRep(const IObjectSP& obj){
    if (IPublicObject::TYPE->isInstance(obj)){
        const IPublicObject& pubObj = dynamic_cast<const IPublicObject&>(*obj);
        return (IObjectSP(pubObj.toPrivateObject()));
    }
    return obj;
}

/** Utility method to convert an object to its public representation
    if it has one. The supplied smart pointer is modified in place */
IObjectConstSP CObject::convertToPublicRep(const IObjectConstSP& obj){
    if (IPrivateObject::TYPE->isInstance(obj)){
        const IPrivateObject& priObj =
            dynamic_cast<const IPrivateObject&>(*obj);
        return (IObjectConstSP(priObj.toPublicObject()));
    }
    return obj;
}
/** Utility method to convert an object to its public representation
    if it has one. The supplied smart pointer is modified in place */
IObjectSP CObject::convertToPublicRep(const IObjectSP& obj){
    if (IPrivateObject::TYPE->isInstance(obj)){
        const IPrivateObject& priObj =
            dynamic_cast<const IPrivateObject&>(*obj);
        return (IObjectSP(priObj.toPublicObject()));
    }
    return obj;
}

/** utility method. Throws an exception if refCount < 1. This is
    useful for ensuring that the object in question is being
    accessed via smartPointers since it means that an extra
    reference to it can be easily made (eg via attachToRef)
    so stopping the object going out of scope */
void CObject::ensurePosRefCount() const{
    if (getRefCount() < 1){
        throw ModelException("CObject::ensurePosRefCount", "Object (type: "+
                             getClass()->getName()+") is not being accessed "
                             "via smart pointers");
    }
}

/** returns this as a void* */
void* CObject::castToBase() const{
    return const_cast<CObject*>(this);
}

/** A very brief English name or description of the object,
    for use in constructing e.g. error messages */
string CObject::toString() const {
    return Format::toString("%s @%X", getClass()->getName().c_str(), this);
}

/** Addin for building [handles to] objects */
class ObjectAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    string         typeName;
    CStringArraySP fieldNames;
    ObjectArraySP  components;

    /** the 'addin function' - construct object from components */
    static IObjectSP createObject(ObjectAddin* params){
        static const string routine = "ObjectAddin::createObject";
        IObjectSP newObj;
        CStringArray& fieldNames = *params->fieldNames;
        ObjectArray&  components = *params->components;
        try{
            CDataDictionarySP dd(CDataDictionary::create(params->typeName));
            if (params->fieldNames.get()){
                for (int i = 0; i < fieldNames.size(); i++){
                    if (!fieldNames[i].empty()){
                        if (!params->components || components.size() <= i){
                            throw ModelException(routine, "Object for field "+
                                                 fieldNames[i]+" is missing");
                        }
                        if (components[i].get()){
                            // skip over null components
                            dd->put(fieldNames[i], components[i]);
                        }
                    }
                }
            }
            //// Note that we are relying on the addin layer to have cloned
            //// the inputs to this function before being called. Therefore
            //// it we can optimise the creation of the object by avoiding
            //// an additional copy of the components.
            dd->copyOnPop2Obj(false);
            newObj = dd->pop2Object();
        } catch (exception& e){
            throw ModelException(e, routine, "Failed to build object of type "+
                                 params->typeName);
        }
             
        return newObj;
    }

    /** for reflection */
    ObjectAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectAddin);
        FIELD(typeName, "Type of object to create");
        FIELD(fieldNames, "names of the components");
        FIELD_MAKE_OPTIONAL(fieldNames);
        FIELD(components, "the components themselves");
        FIELD_MAKE_OPTIONAL(components);
        Addin::registerClassObjectMethod("OBJECT",
                                         Addin::UTILITIES,
                                         "Constructs a handle to an object of "
                                         "specified type using the components"
                                         " provided",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)createObject);
    }

    static IObject* defaultObjectAddin(){
        return new ObjectAddin();
    }
    
};

CClassConstSP const ObjectAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectAddin", typeid(ObjectAddin), load);



/** Class for addin functions that take a single object as the sole parameter */
class GenericSingleObjectAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object;

    /** the 'addin function' - construct object from components */
    string getTypeKey(){
        return object->getClass()->getName();
    }

    /** another 'addin function' - invoke an 'executable' type */
    IObjectSP actionInvoke(){
        static const string routine = "GenericSingleObjectAddin::actionInvoke";
        
        if (ClientRunnable::TYPE->isInstance(object)) {
            ClientRunnable& action =
                dynamic_cast<ClientRunnable&>(*object);
            IObjectSP result(action.run());
            return result;
        } else {
            throw ModelException(routine, "Object is not client runnable");
        }
    }

    /** another 'addin function' - is the object an array? */
    bool isArrayType(){
        return (IArray::TYPE->isInstance(object.get()));
    }

    /** another 'addin function' - is the object an Enum */
    bool isEnumType(){
        return object->getClass()->isEnum();
    }

    /** for reflection */
    GenericSingleObjectAddin():  CObject(TYPE){}

    static IObject* defaultConstructor() {
        return new GenericSingleObjectAddin();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GenericSingleObjectAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(object, "object handle");
        Addin::registerStringMethod("TYPE_KEY_GET",
                                    Addin::UTILITIES,
                                    "determines the run-time type of an "
                                    " object handle",
                                    &GenericSingleObjectAddin::getTypeKey);
        Addin::registerObjectMethod(
            "ACTION_INVOKE",
            Addin::UTILITIES,
            "invokes the action on the given object",
            true,
            Addin::returnHandle,
            &GenericSingleObjectAddin::actionInvoke);
        Addin::registerObjectMethod(
            "OBJECT_EXECUTE",
            Addin::UTILITIES,
            "invokes object's action",
            true,
            Addin::returnHandle,
            &GenericSingleObjectAddin::actionInvoke);
        Addin::registerBoolMethod("TYPE_IS_ARRAY",
                                  Addin::UTILITIES,
                                  "determines whether the handle "
                                  "represents an array object",
                                  &GenericSingleObjectAddin::isArrayType);
        Addin::registerBoolMethod("TYPE_IS_ENUM",
                                  Addin::UTILITIES,
                                  "determines whether the handle "
                                  "represents an Enum object",
                                  &GenericSingleObjectAddin::isEnumType);
    }

};

CClassConstSP const GenericSingleObjectAddin::TYPE = 
CClass::registerClassLoadMethod(
    "GenericSingleObjectAddin", typeid(GenericSingleObjectAddin), load);

/** Class for addin functions that take two objects as the sole parameters */
class GenericObjectPairAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object1;
    IObjectSP  object2; // optional

    /** another 'addin function' - is the object an array? */
    bool objectEquals(){
        return (object1->equalTo(object2.get()));
    }

    /** Returns a hash code for the supplied object. If a second object is
        provided then it returns an exclusive or of the pair of hashcodes.
        This can be useful to try and remove the runtime dependence of some
        hashcodes */
    int objectHashCode(){
        return (object1->hashCode() ^ (!object2? 0: object2->hashCode()));
    }

    /** for reflection */
    GenericObjectPairAddin():  CObject(TYPE){}

    static IObject* defaultConstructor() {
        return new GenericObjectPairAddin();
    }

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(GenericObjectPairAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(object1, "first object");
        FIELD(object2, "second object");
        FIELD_MAKE_OPTIONAL(object2);
        Addin::registerBoolMethod("EQUALS",
                                  Addin::UTILITIES,
                                  "Determines whether two objects are equal",
                                  &GenericObjectPairAddin::objectEquals);
        Addin::registerIntMethod(
            "HASHCODE",
            Addin::UTILITIES,
            "Returns a hash code for the supplied object. If a second object is"
            "provided then it returns an exclusive or of the pair of hashcodes."
            "This can be useful to try and remove the runtime dependence of "
            "some hashcodes",
            &GenericObjectPairAddin::objectHashCode);
    }
};

CClassConstSP const GenericObjectPairAddin::TYPE = 
CClass::registerClassLoadMethod(
    "GenericObjectPairAddin", typeid(GenericObjectPairAddin), load);

/** Addin to determine run-time type of a handle */
class PrivateToPublicAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP object;

    /** the 'addin function' - construct object from components */
    static IObjectConstSP convertToPublic(PrivateToPublicAddin* params){
        static const string routine = "PrivateToPublicAddin::convertToPublic";

        return CObject::convertToPublicRep(params->object);
    }

    /** for reflection */
    PrivateToPublicAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(PrivateToPublicAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultPrivateToPublicAddin);
        FIELD(object, "object handle");
        Addin::registerClassObjectMethod("PRIVATE_TO_PUBLIC",
                                         Addin::UTILITIES,
                                         "converts a private object to it's public"
                                         "representation",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)convertToPublic);
    }

    static IObject* defaultPrivateToPublicAddin() {
        return new PrivateToPublicAddin();
    }
};

CClassConstSP const PrivateToPublicAddin::TYPE = CClass::registerClassLoadMethod(
    "PrivateToPublicAddin", typeid(PrivateToPublicAddin), load);


/** Addin to determine run-time type of a handle */
class XLIsPrivateTypeAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    IObjectSP  object;

    /** the 'addin function' - construct object from components */
    static bool isPrivateType(XLIsPrivateTypeAddin* params){
        static const string routine = "XLIsPrivateTypeAddin::isPrivateType";

        return (Modifier::isPrivate(params->object->getClass()->getModifiers()));
    }

    /** for reflection */
    XLIsPrivateTypeAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLIsPrivateTypeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIsPrivateTypeAddin);
        FIELD(object, "object handle");
        Addin::registerClassBoolMethod("TYPE_IS_PRIVATE",
                                       Addin::UTILITIES,
                                       "determines whether the object is of a private type",
                                       TYPE,
                                       (Addin::BoolMethod*)isPrivateType);
    }

    static IObject* defaultIsPrivateTypeAddin() {
        return new XLIsPrivateTypeAddin();
    }
};

CClassConstSP const XLIsPrivateTypeAddin::TYPE = CClass::registerClassLoadMethod(
    "XLIsPrivateTypeAddin", typeid(XLIsPrivateTypeAddin), load);


class ObjectSaveAddin: public CObject {
    static CClassConstSP const TYPE;

    string    filename;
    IObjectSP obj;

    /** the 'addin function' */
    static bool save(ObjectSaveAddin* params){
        static const string method = "ObjectSaveAddin::save";
        try {
            XMLWriter xml(params->filename);
            params->obj->write("OBJECT", &xml);
            return true;
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ObjectSaveAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectSaveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectSaveAddin);
        FIELD(filename, "filename");
        FIELD(obj, "object");

        Addin::registerClassBoolMethod("OBJECT_SAVE",
                                       Addin::UTILITIES,
                                       "save an object to file",
                                       TYPE,
                                       (Addin::BoolMethod*)save);
    }

    static IObject* defaultObjectSaveAddin(){
        return new ObjectSaveAddin();
    }   
};

CClassConstSP const ObjectSaveAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectSaveAddin", typeid(ObjectSaveAddin), load);

class ObjectLoadAddin: public CObject {
    static CClassConstSP const TYPE;

    string    filename;

    /** the 'addin function' */
    static IObjectSP loadObject(ObjectLoadAddin* params){
        static const string method = "ObjectLoadAddin::loadObject";
        try {
            XMLReader reader(params->filename, true);
            IObjectSP obj(reader.read());
            return obj;
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ObjectLoadAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectLoadAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectLoadAddin);
        FIELD(filename, "filename");

        Addin::registerInstanceObjectMethod("OBJECT_LOAD",
                                            Addin::UTILITIES,
                                            "Load an object from file",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)loadObject);
    }

    static IObject* defaultObjectLoadAddin(){
        return new ObjectLoadAddin();
    }   
};

CClassConstSP const ObjectLoadAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectLoadAddin", typeid(ObjectLoadAddin), load);


class ObjectCopyAddin: public CObject {
    static CClassConstSP const TYPE;
    IObjectSP obj;

    /** the 'addin function' */
    static IObjectSP copyObject(ObjectCopyAddin* params){
        static const string method = "ObjectCopyAddin::copyObject";
        try {
            return IObjectSP(params->obj->clone());
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    ObjectCopyAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(ObjectCopyAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultObjectCopyAddin);
        FIELD(obj, "object");

        Addin::registerInstanceObjectMethod("OBJECT_COPY",
                                            Addin::UTILITIES,
                                            "copy an object",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)copyObject);
    }

    static IObject* defaultObjectCopyAddin(){
        return new ObjectCopyAddin();
    }   
};

CClassConstSP const ObjectCopyAddin::TYPE = CClass::registerClassLoadMethod(
    "ObjectCopyAddin", typeid(ObjectCopyAddin), load);

class IsAssignableFromAddin: public CObject {
private:
    
    static CClassConstSP const TYPE;

    string  baseClassName;
    string  derivedClassName;

    /** the 'addin function' */
    static bool isAssignableFrom(IsAssignableFromAddin* params){
        static const string method = "IsAssignableFromAddin::isAssignableFrom";
        try {
            CClassConstSP baseClazz(CClass::forName(params->baseClassName));
            CClassConstSP derivedClazz(CClass::forName(params->derivedClassName));

            bool result = baseClazz->isAssignableFrom(derivedClazz);
            return result;
        } 
        catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** for reflection */
    IsAssignableFromAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(IsAssignableFromAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIsAssignableFromAddin);
        FIELD(baseClassName,    "Base class name");
        FIELD(derivedClassName, "Derived class name");

        Addin::registerClassBoolMethod("IS_ASSIGNABLE_FROM",
                                       Addin::UTILITIES,
                                       "Determines if base derives from or is equal to derived class",
                                       TYPE,
                                       (Addin::BoolMethod*)isAssignableFrom);
    }

    static IObject* defaultIsAssignableFromAddin(){
        return new IsAssignableFromAddin();
    }   
};

CClassConstSP const IsAssignableFromAddin::TYPE = CClass::registerClassLoadMethod(
    "IsAssignableFromAddin", typeid(IsAssignableFromAddin), load);

// to write to a xml file of an object, bypassing dri provided xml writer
class XMLDumper : public CObject, public ClientRunnable
{
public:
    static CClassConstSP const TYPE;
    virtual ~XMLDumper() {};
    IObjectSP run(); 

protected:
    XMLDumper(CClassConstSP clazz) : CObject(clazz), object(0) {}
    CObject* object;
    string   pathname;

private:
    XMLDumper() : CObject(TYPE), object(0) {}
    static IObject* defaultXMLDumper() { return new XMLDumper(); }
    static void load(CClassSP& clazz);
};

void XMLDumper::load(CClassSP& clazz)
{
    REGISTER(XMLDumper, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ClientRunnable);
    EMPTY_SHELL_METHOD(defaultXMLDumper);
    FIELD(pathname, "pathname for dumping the object in XML format");
    FIELD(object,   "object to be dumped");
    clazz->setPublic(); 
}

CClassConstSP const XMLDumper::TYPE = CClass::registerClassLoadMethod(
    "XMLDumper", typeid(XMLDumper), XMLDumper::load);

// * for class loading (avoid having header file) */
bool XMLDumperLoad() {
    return (XMLDumper::TYPE != 0);
}

IObjectSP XMLDumper::run()
{
    try
    {
        XMLWriter xml(pathname);
        object->write("OBJECT", &xml);
        return IObjectSP(CBool::create(true));
    }
    catch (exception& e)
    {
        return IObjectSP(CBool::create(false));
    }
}

DRLIB_END_NAMESPACE

