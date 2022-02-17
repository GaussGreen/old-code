//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DRWrapper.cpp
//
//   Description : It's DR Wrapper
//
//   Author      : Andrew J Swain
//
//   Date        : 19 April 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DRWrapper.hpp"
#include "edginc/DataDictionary.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CommandLineParams.hpp"
#include "edginc/Writer.hpp"
#include "edginc/Instrument.hpp"

using namespace std;  // string

DRLIB_BEGIN_NAMESPACE
DEFINE_TEMPLATE_TYPE(DRWrapperArray);

// the DRWrapper iterator
class DRWrapperIterator: public CObject,
                         public virtual IMap::IIterator {
public:
    static CClassConstSP const TYPE;
    // are there any elements left to iterate over
    virtual bool hasMoreElements() const {
        return !(iter == hash.end());
    }

    // get the key for the current element
    virtual const string& getKey() const {
        return iter->first;
    }

    // get the current element (returns a reference)
    virtual IObjectSP getElement() const {
        return iter->second;
    }
    // increment the iterator
    virtual void increment() {
        iter++;
    }

    DRWrapperIterator(WrapperHashTable hash): 
        CObject(TYPE), hash(hash) {
        iter = this->hash.begin();
    }
    
    static void load(CClassSP& clazz){
        REGISTER(DRWrapperIterator, clazz);
        SUPERCLASS(CObject)
        IMPLEMENTS(IMap::IIterator);
        EMPTY_SHELL_METHOD(defaultDRWrapperIterator);
    }

    static IObject* defaultDRWrapperIterator(){
        return new DRWrapperIterator();
    }

private:
    DRWrapperIterator():CObject(TYPE) {};
    // data
    WrapperHashTable                 hash;
    WrapperHashTable::const_iterator iter;
};

CClassConstSP const DRWrapperIterator::TYPE = CClass::registerClassLoadMethod(
    "DRWrapperIterator", typeid(DRWrapperIterator), DRWrapperIterator::load);



/** XML tag for recording wrapper type */
const string DRWrapper::WRAPPER_TYPE = "wrapperType";

/** Create a DRWrapper given object's name, such as "MyNewProduct" */
DRWrapper* DRWrapper::create(const string& typekey) {
    DRWrapper* wrapper = 0;
    try {
        wrapper = new DRWrapper(typekey);
    }
    catch (exception& e) {
        throw ModelException(e, "DRWrapper::create",
                             "couldn't build a DR Wrapper for" + typekey);
    }
    return wrapper;
}

/** Add an object to a DRWrapper */
void DRWrapper::put(const string& key, const IObjectSP& value) {
    hash[key] = CObject::convertToPrivateRep(value);
}

/** Create an object from a DRWrapper */
IObjectSP DRWrapper::pop2Object() const {
    static const string method = "DRWrapper::pop2Object";
    IObjectSP obj;

    try {
        CDataDictionarySP dd=CDataDictionarySP(CDataDictionary::create(type));

        for (WrapperHashTable::const_iterator myIterator = hash.begin();
             !(myIterator == hash.end());
             ++myIterator){

            if (!myIterator->second) {
                // do nothing
            }
            // if element of wrapper is a wrapper itself, convert it first
            else if (DRWrapper::TYPE->isInstance(myIterator->second)) {
                const DRWrapper& wrapper = dynamic_cast<const DRWrapper&>(static_cast<const IObject&>(*myIterator->second));

                IObjectSP popper = wrapper.pop2Object();
                dd->put(myIterator->first, popper); 
            }
            else if (DRWrapperArray::TYPE->isInstance(myIterator->second)) {
                const DRWrapperArray& wrapper = dynamic_cast<const DRWrapperArray&>(static_cast<const IObject&>(*myIterator->second));
                // look up the type of the array 
                CClassConstSP clazz = CClass::forName(type);
                CFieldConstSP field;

                // iterate over a object and all its superclasses, looking for
                // the field
                do {
                    field = clazz->hasDeclaredField(myIterator->first);
                } while (!field && (clazz = clazz->getSuperClass()) != 0);

                if (!field) {
                    throw ModelException(method,
                                         "field " + myIterator->first + 
                                         " not found in " + type);
                }

                IArraySP array;
                // if this is really an array, build one - otherwise it may be that
                // we have an array at the IMS end which gets converted into a non-array
                // so try and build an array of the first element's type
                if (field->getType()->isArray()) {
                    array = IArraySP(field->getType()->newArrayInstance(wrapper.size()));
                }
                else if (!wrapper.empty()) {
                    CClassConstSP elemClazz = CClass::forName(wrapper[0]->type);
                    array = IArraySP(elemClazz->newArrayInstanceByComponent(wrapper.size()));
                }
                else {
                    throw ModelException(method, "field "+myIterator->first+
                                         " is not an array and can't be made "
                                         "from an empty DR wrapper array");
                }
                    
                for (int i = 0; i < wrapper.size(); i++) {
                    IObjectSP popper = wrapper[i]->pop2Object();
                    array->set(i, popper);
                }
                dd->put(myIterator->first, array);                
            }
            else {
                dd->put(myIterator->first, myIterator->second);             
            }
        }

        obj = dd->pop2Object();
    }
    catch (exception& e) {
        throw ModelException(e, method, "couldn't build " + type);
    }
    return obj;
}

/** check the executable name matches the wrapper */
void DRWrapper::validateExeName() const {
    string exename = CommandLineParams::executableName();
    if (exename != string("") && exename != type) {
        throw ModelException("DRWrapper::validateExeName",
                             "executable name (" + exename + 
                             ") is not the same as the wrapper type (" + 
                             type + ")");
    }
}

/** Helper class for xml read/write for hash table of 
    <string, IObjectSP> */
class WrapperEntry: public CObject{
public:
    static CClassConstSP const TYPE;
    string             key;
    IObjectSP          value;
    static const string XML_TAG;
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(WrapperEntry, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultWrapperEntry);
        FIELD(key, "key");
        FIELD(value, "value");
        FIELD_MAKE_OPTIONAL(value);  // allow Null
    }

    static IObject* defaultWrapperEntry(){
        return new WrapperEntry();
    }
    // for reflection
    WrapperEntry(): CObject(TYPE){}

    WrapperEntry(const string&    key, 
                 const IObjectSP& value):
        CObject(TYPE), key(key), value(value){}
};
const string WrapperEntry::XML_TAG = "WrapperEntry";

CClassConstSP const WrapperEntry::TYPE =
CClass::registerClassLoadMethod("WrapperEntry", typeid(WrapperEntry), load);

/** write object out in XML format */
void DRWrapper::write(const string& tag, Writer* writer) const {
    char buffer[256];
    try {
        sprintf(buffer, "%s='%s'",
                WRAPPER_TYPE.c_str(), type.c_str());
        IObjectConstSP obj(writer->objectStart(tag, buffer, this, true));
        if (obj.get()){
            // if we have implemented IPrivateObject then we need to act
            // on obj here which will be of type IPublicObject
            for (WrapperHashTable::const_iterator myIterator = hash.begin();
                 !(myIterator == hash.end());
                 ++myIterator){

                WrapperEntry entry(myIterator->first, myIterator->second);
                entry.write(WrapperEntry::XML_TAG, writer);
            }
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "DRWrapper::xmlWrite");
    }
}

/** populate an empty object from XML description */
void DRWrapper::import(Reader::Node* elem, Reader* reader){
    const static string method = "DRWrapper::import";
    try{
        // first get type of wrapper
        type = elem->attribute(WRAPPER_TYPE);
    
        Reader::NodeListSP entryNodes(elem->children());
        for (unsigned int i = 0; i < entryNodes->size(); i++) {
            Reader::NodeSP entry((*entryNodes)[i]);
            // read in entry
            IObjectSP entryObj = IObjectSP(reader->read(entry.get()));
            WrapperEntry* wrapperEntry = 
                dynamic_cast<WrapperEntry*>(entryObj.get());
            if (!wrapperEntry){
                throw ModelException(method,
                                     "Invalid object found");
            }
            
            put(wrapperEntry->key, wrapperEntry->value);
        }
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}
 
/** override clone method */
IObject* DRWrapper::clone() const {
    try{
        DRWrapper* copy = new DRWrapper(type);
        // iterate over all elements
        for (WrapperHashTable::const_iterator myIterator = hash.begin();
             !(myIterator == hash.end());
             ++myIterator){
            string    key = myIterator->first;
            IObjectSP component(myIterator->second.clone());
            copy->hash[key] = component;
        }
        return copy;
    } 
    catch (exception& e) {
        throw ModelException(&e, "DRWrapper::clone");
    }
}

// If object is a drwrapper convert to equivalent object, otherwise returns
// original
IObjectSP DRWrapper::drWrapperToObject(const IObjectSP& object,
                                       CClassConstSP    desiredType){
    IObjectSP newObj;
    if (DRWrapper::TYPE->isInstance(object)) {
        DRWrapper& wrapper = dynamic_cast<DRWrapper&>(*object);
        if (desiredType->isAssignableFrom(Instrument::TYPE)) {
            wrapper.validateExeName();
        }
        newObj = wrapper.pop2Object();
    } else {
        newObj = object;
    }
    CObject::checkType (newObj, desiredType); 
    return newObj;
}

/** build an iterator */
IMap::IIteratorSP DRWrapper::createIterator() {
    return IMap::IIteratorSP(new DRWrapperIterator(hash));
}

/** Is this object truly a map i.e. does toObject()/toMap() return this */
bool DRWrapper::isTrueMap() const {
    return false; // we have other data too therefore we MUST return false
}

/** for reflection */
DRWrapper::DRWrapper():CObject(TYPE), type(""){} 

DRWrapper::DRWrapper(const string& type):CObject(TYPE), type(type){}


class DRWCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // three addin parameters
    string           typeName;   // gives type of wrapper to create
    // remaining two parameters are optional
    CStringArraySP   fieldNames; // identifies objects
    ObjectArraySP    objects;    // the objects themselves

    // utility method
    static void fillWrapper(DRWrapperSP&          wrapper,
                            const CStringArray&   fieldNames,
                            const ObjectArraySP   objects) { // can be null
        
        static const string routine = "fillWrapper";
        for (int i = 0; i < fieldNames.size(); i++){
            if (!fieldNames.empty()){
                if (!objects || i >= objects->size()){
                    throw ModelException(routine, "Object for field "+
                                         fieldNames[i]+" is missing");
                }
                wrapper->put(fieldNames[i], (*objects)[i]);
            }
        }
    }

    /** create a data dictionary */
    static IObjectSP create(DRWCreateAddin *params){
        DRWrapperSP wrapper(DRWrapper::create(params->typeName));
        if (params->fieldNames.get()){
            fillWrapper(wrapper, *params->fieldNames, params->objects);
        }
        return wrapper;
    }

    DRWCreateAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DRWCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDRWCreateAddin);
        FIELD(typeName, "Type of DR Wrapper");
        FIELD(fieldNames, "Identifies the object");
        FIELD_MAKE_OPTIONAL(fieldNames);
        FIELD(objects, "The object itself");
        FIELD_MAKE_OPTIONAL(objects);
        Addin::registerClassObjectMethod("DR_WRAPPER",
                                         Addin::UTILITIES,
                                         "Creates a handle to a DR "
                                         "wrapper",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }

    static IObject* defaultDRWCreateAddin(){
        return new DRWCreateAddin();
    }
};

CClassConstSP const DRWCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "DRWCreateAddin", typeid(DRWCreateAddin), load);

class DRWPutAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // three addin parameter
    DRWrapperSP        wrapper;    // the DR Wrapper
    CStringArraySP     fieldNames; // identifies objects
    ObjectArraySP      objects;    // the objects themselves

    /** set an object in a data dictionary */
    static IObjectSP set(DRWPutAddin *params){
        DRWCreateAddin::fillWrapper(params->wrapper, *params->fieldNames,
                                    params->objects);
        return params->wrapper; // return the now modified wrapper
    }

    DRWPutAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DRWPutAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDRWPutAddin);
        FIELD(wrapper, "DR Wrapper to put object into");
        FIELD(fieldNames, "Identifies the object");
        FIELD(objects, "The object itself");
        Addin::registerInstanceObjectMethod(
            "DR_WRAPPER_PUT",
            Addin::UTILITIES,
            "Adds objects to a DR Wrapper",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)set);
    }

    static IObject* defaultDRWPutAddin(){
        return new DRWPutAddin();
    }
};

CClassConstSP const DRWPutAddin::TYPE = CClass::registerClassLoadMethod(
    "DRWPutAddin", typeid(DRWPutAddin), load);

class DRWPopAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // one addin parameter
    DRWrapperSP  wrapper;  // the DR Wrapper

    /** pop data dictionary to an object */
    static IObjectSP pop(DRWPopAddin *params){
        return params->wrapper->pop2Object();
    }

    DRWPopAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DRWPopAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDRWPopAddin);
        FIELD(wrapper, "DR Wrapper to convert");
        Addin::registerInstanceObjectMethod("DR_WRAPPER_POP",
                                            Addin::UTILITIES,
                                            "Turns a DR Wrapper into the"
                                            " corresponding object",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)pop);
    }

    static IObject* defaultDRWPopAddin(){
        return new DRWPopAddin();
    }
};

CClassConstSP const DRWPopAddin::TYPE = CClass::registerClassLoadMethod(
    "DRWPopAddin", typeid(DRWPopAddin), load);


bool DRWrapper::load() {
    return (DRWrapper::TYPE && 
            DRWCreateAddin::TYPE && 
            DRWPutAddin::TYPE && 
            DRWPopAddin::TYPE);
}

// a proxy interface for DRWrapper so fits in with Global DR Interface
// for serialisation etc.
class DRWrapperProxy: public CObject {
public:
    static CClassConstSP const TYPE;
private:
    string      type;     
    HashtableSP hash;

    typedef smartPtr<DRWrapperProxy> DRWrapperProxySP;

    static IObjectSP toDRWrapper(const IObjectConstSP& obj) {
        const DRWrapperProxy* proxy = dynamic_cast<const DRWrapperProxy*>(obj.get());
        if (!proxy) {
            throw ModelException("DRWrapperProxy::toDRWrapper",
                                 "object is not a DRWrapperProxy");
        }
    
        DRWrapperSP wrapper(new DRWrapper(proxy->type));
        if (proxy->hash.get()) {
            // now fill in the wrapper data
            IMap::IIteratorSP iter(proxy->hash->createIterator());
            while (iter->hasMoreElements()) {
                wrapper->put(iter->getKey(), iter->getElement());
                iter->increment();
            }
        }
        return wrapper;
    }

    static IObjectSP fromDRWrapper(const IObjectConstSP& obj) {
        const DRWrapper* wrapper = dynamic_cast<const DRWrapper*>(obj.get());
        if (!wrapper) {
            throw ModelException("DRWrapperProxy::fromDRWrapper",
                                 "object is not a DRWrapper");
        }
    
        DRWrapper* drw = const_cast<DRWrapper*>(wrapper);

        DRWrapperProxySP proxy(new DRWrapperProxy(drw->type));
        // now fill the proxy with the wrapper's data
        IMap::IIteratorSP iter(drw->createIterator());
        while (iter->hasMoreElements()) {
            proxy->hash->put(iter->getKey(), iter->getElement());
            iter->increment();
        }
        return proxy;
    }

    DRWrapperProxy(const string& type) : CObject(TYPE),type(type), 
        hash(new Hashtable) {}

    /** for reflection */
    DRWrapperProxy():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); 
        REGISTER(DRWrapperProxy, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDRWrapperProxy);
        registerObjectProxy(DRWrapper::TYPE,
                            DRWrapperProxy::TYPE,
                            fromDRWrapper,
                            toDRWrapper);
        FIELD(type, "type");
        FIELD(hash, "hash");
        FIELD_MAKE_OPTIONAL(hash);
    }

    static IObject* defaultDRWrapperProxy(){
        return new DRWrapperProxy();
    }
    
};

CClassConstSP const DRWrapperProxy::TYPE = CClass::registerClassLoadMethod(
    "DRWrapperProxy", typeid(DRWrapperProxy), load);

class DRWrapperHelper{
public:
/** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setDRIProxyType(DRWrapperProxy::TYPE); // use proxy for dri
        REGISTER(DRWrapper, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IWriteableMap);
        EMPTY_SHELL_METHOD(defaultDRWrapper);
        FIELD(type, "type this wrapper represents");
        clazz->setPublic();  // so can make arrays of wrappers
    }
    
    static IObject* defaultDRWrapper(){
        return new DRWrapper();
    }
};

CClassConstSP const DRWrapper::TYPE = CClass::registerClassLoadMethod(
    "DRWrapper", typeid(DRWrapper), DRWrapperHelper::load);


DRLIB_END_NAMESPACE
