

#include "edginc/config.hpp"
#define QLIB_DATADICTIONARY_CPP
#include "edginc/DataDictionary.hpp"
#include "edginc/Null.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<CDataDictionary>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<CDataDictionary>);

CDataDictionary::~CDataDictionary(){}

/** A DataDictionary is a form of hash table. It is always "typed" - 
        it contains the data that describes an object's public interface
    */
    /** set flag indicating whether copies of the objects in the
            data dictionary should be used when building an object in
            {@link #pop2Object pop2Object}. Default is true */
void CDataDictionary::copyOnPop2Obj(bool copyOnPop2Obj) throw() {
    this->copyOnPop2Obj_ = copyOnPop2Obj;
}

/** Create a DataDictionary given object's name, such as "Equity" */
CDataDictionary* CDataDictionary::create(const string& typeKey){
    static const string method("CDataDictionary::create");
    try{
        CClassConstSP type = CClass::forName(typeKey);
        int modifiers = type->getModifiers();
        if (!(Modifier::isPublic(modifiers))){
            throw ModelException(method,
                                 "Creation of data dictionary for non public"
                                 " type ("+type->getName()+") disallowed");
        }
        return new CDataDictionary(type);
    } catch (exception &e){
        throw ModelException(e, method, "Failed to create data dictionary "
                             "of type "+typeKey);
    }
}

/** Create an empty DataDictionary given a class */
CDataDictionary* CDataDictionary::create(const CClassConstSP& clazz){
    return new CDataDictionary(clazz);
}

CDataDictionary::CDataDictionary(const CClassConstSP& type): 
    CObject(TYPE), type(type), copyOnPop2Obj_(true) {
    static const string routine("CDataDictionary::CDataDictionary");
    // catch abstract types now - get a better error message
    int modifiers = type->getModifiers();
    if (Modifier::isAbstract(modifiers)){
        // This setting is currently determined by whether a 
        // EMPTY_SHELL_METHOD has been set for the type or not
        throw ModelException(routine,
                             "Cannot build a data dictionary for abstract"
                             " type: "+type->getName());
    }        
    if (type->isArray()){
        throw ModelException(routine, "This function cannot be used"
                             " for building arrays");
    }
}

/** What type of DataDictionary is this ? */
CClassConstSP CDataDictionary::getType() const throw(){
    return type;
}

/** Add an object to a DataDictionary - works its way up the
        inheritance chain, matching the first data member whose name
        matches the key * @param key Data member of clazz * @param
        value Object to add */
void CDataDictionary::put(const string& key, const IObjectSP& value){
    static const string method = "CDataDictionary::put";
    try {
        if (!value) {
            throw ModelException(method, "value for " + key + " for object "
                                 "of type "+type->getName()+" is null");
        }

        CClassConstSP  c = type;
        CFieldConstSP field;
        do {
            field = c->hasDeclaredField(key);
        } while (!field && (c = c->getSuperClass()) != 0);
        
        if (!field) {
            throw ModelException("CDataDictionary::put",
                                 key+" not found in "+type->getName());
        }
        // finally set value
        put(field, value);
    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}

/** Same as above put will clone value if refCount = 0 */
void CDataDictionary::put(const string& key, IObject* value){
    IObjectSP newValue(value->getRefCount() == 0? value->clone(): value);
    put(key, newValue);
}

//// nicer versions of put for the primitives
void CDataDictionary::put(const string& key, bool b){
    put(key, IObjectSP(CBool::create(b)));
}
void CDataDictionary::put(const string& key, double d){
    put(key, IObjectSP(CDouble::create(d)));
}
void CDataDictionary::put(const string& key, int i){
    put(key, IObjectSP(CInt::create(i)));
}
void CDataDictionary::put(const string& key, const string& s){
    put(key, IObjectSP(CString::create(s)));
}
void CDataDictionary::put(const string& key, const char* s){
    put(key, string(s));
}
    
/** Add an object to a DataDictionary given a specific Field. Not
    public as it does not check that the associated object actually
    has the field supplied */
void CDataDictionary::put(CFieldConstSP field, IObjectSP value){
    static const string method("CDataDictionary::put");
    if (Modifier::isTransient(field->getModifiers())){
        throw ModelException(method, "Field "+field->getName()+" is transient"
                             " and cannot be set");
    }
    // finally set value
    hash[field] = value;
}

/** Converts this object to an instance of the requiredType. Throws an
    exception if a conversion to the required Type is not supported */
void CDataDictionary::convert(IObjectSP&    object,
                              CClassConstSP requiredType) const{
    // first convert into corresponding object
    object = IObjectSP(pop2Object());
    // then check that type
    CObject::checkType(object, requiredType);
}    

/** Add an object to a DataDictionary - uses fully qualified
    class and member names so there is no clash of member names
    with any superclass 
    * @param clazz The fully qualified class name 
    * @param member Data member of clazz 
    * @param value Object to add */
void CDataDictionary::put(const string& clazz, 
                          const string& member, IObjectSP& value){
    static const string method = "CDataDictionary::put";
    try {
        if (!value) {
            throw ModelException(method, "value for " + member + " for object "
                                 "of type "+type->getName()+" is null");
        }

        CClassConstSP clazzType = CClass::forName(clazz);
        CFieldConstSP field = clazzType->getDeclaredField(member);
        // finally set value
        put(field, value);
    } catch (exception& e){
        throw ModelException(&e, method, "Failed to set "+clazz + "." + 
                             member + " in " + type->getName());
    }
}
    
bool CDataDictionary::contains(const string& key) const{
    const static string routine("CDataDictionary::contains");
    CClassConstSP  c = type;
    CFieldConstSP  field;
    do {
        field = c->hasDeclaredField(key);
    } while (!field && (c = c->getSuperClass()) != 0);
    if (!field) {
        throw ModelException(routine,
                             key+" not found in "+type->getName());
    }
    return !(hash.find(field) == hash.end());
}
    
/** Get an object out of a DataDictionary - works from base
    class upwards looking for a field that matches the key - be
    wary of name clashes with inherited fields. Throws an exception
    if the key is invalid */
IObjectSP CDataDictionary::get(const string& key) const{
    const static string routine("CDataDictionary::get");
    CClassConstSP  c = type;
    CFieldConstSP  field;
    do {
        field = c->hasDeclaredField(key);
    } while (!field && (c = c->getSuperClass()) != 0);
    if (!field) {
        throw ModelException(routine,
                             key+" not found in "+type->getName());
    }
    return get(field);
}

/** Get an object out of a DataDictionary - uses fully
        qualified class and member names so there is no clash of
        member names with any superclass * @param clazz The fully
        qualified class name * @param member Data member of clazz */
IObjectSP CDataDictionary::get(const string& clazz,
                               const string& member) const{
    static const string routine("CDataDictionary::get");
    try{
        CClassConstSP clazzType = CClass::forName(clazz);
        CFieldConstSP field = clazzType->getDeclaredField(member);
        return get(field);
    } catch (exception& e){
        throw ModelException(&e, routine,
                             "Failed to get "+clazz + "." + member + 
                             " in " + type->getName());
    }
}

/** Get a field out of a DataDictionary. The DataDictionary has the 'public'
    view of the object so any object which are 'private' are converted on
    the fly (ie on demand) */
IObjectSP CDataDictionary::get(CFieldConstSP& field) const{
    ObjectHashTable::iterator iter = hash.find(field);
    if (iter == hash.end()){
        // DRI 1.2 
        return CNull::create(); // must distinguish between invalid fieldName
                                // and missing field
#if 0
        throw ModelException("CDataDictionary::get", 
                             "Element not found for field "+ field->getName());
#endif
    }
    IObjectSP elt(iter->second);
    IObjectSP publicElt(CObject::convertToPublicRep(elt));
    if (elt.get() != publicElt.get()){
        // since we've converted it, might as well save it
        iter->second = publicElt;
    }
    return publicElt;
}
    
/** Create a DataDictionary containing the data that makes up
    an object */
CDataDictionary* CDataDictionary::pop2DataDict(IObjectSP obj){
    try {
        // convert to object's public representation (if it has one)
        // Generally though the object should already been converted
        IObjectSP pubObj(CObject::convertToPublicRep(obj));
        CClassConstSP c = pubObj->getClass();

        CDataDictionarySP dd(CDataDictionary::create(c));
        // iterate over a object and all its superclasses, extracting
        // fields in each level of inheritance and storing them and their
        // associated values in the hash table
        
        do {
            const CFieldArray& fields = c->getDeclaredFields();

            for (unsigned int i = 0; i < fields.size(); i++) {
                if (isPoppable(fields[i])) {
                    IObjectSP o = fields[i]->get(pubObj);
                    if (!o) {
                        o = CNull::create();
                    }

                    // Make sure we don't have problems with scope
                    // e.g. if obj gets deleted.
                    // Also,
                    // DRI 1.2: Copy all containers, i.e., array and matrix.
                    CClassConstSP clazz = o->getClass();
                    if (o->getRefCount() == 0 || 
                        clazz->isArray() ||
                        clazz->isMatrix()) {
                        dd->hash[fields[i]] = IObjectSP(o->clone());
                    }
                    else {
                        dd->hash[fields[i]] = o;
                    }
                }
            }
        } while ((c = c->getSuperClass()) != 0);
        
        return dd.release();
    }
    catch (exception &e) {
        throw ModelException (e, "CDataDictionary::pop2DataDict",
                              "couldn't pop " + obj->getClass()->getName() + 
                              " into data dict");
    }
}

void CDataDictionary::checkElementType(CFieldConstSP field, 
                                       IObjectSP& value) const{
    // check that the object is of the right type.
    // Also be flexible if data dictionaries
    // are supplied instead of the actual object
    CClassConstSP reqType = field->getType(); // required type
    try{
        CObject::checkType(value, reqType);
    } catch (exception& e){
        throw ModelException(e, "DataDictionary::checkElementType",
                             "Object for field "+ field->getName()+
                             " is of type "+value->getClass()->getName()+
                             ".\nRequire type "+reqType->getName() + 
                             " for object of type " + type->getName());
    }
}

/** Populate an existing object with the contents of a data dictionary.
    As per pop2Object, all mandatory fields must be present.
    Really only for use by CObject::xmlImport */
void CDataDictionary::populateObject(IObject* obj) const{
    const static string routine = "CDataDictionary::populateObject";
    if (obj->getClass() != type){
        throw ModelException(routine,
                             "Object is of type: "+obj->getClass()->getName()+
                             "\nDataDictionary is for type: "+type->getName());
    }
    try{
        CClassConstSP   c  = type;
        // field methods want an IObjectSP - perhaps should make them more
        // flexible
        IObjectSP       object(IObjectSP::attachToRef(obj));
   
        // iterate over the type description, pull data for each item out
        // of the hash table and store it in the shell we've just made.
        // if anything's missing from the hash table, then we fail unless
        // it's optional
        do {
            const CFieldArray& fields = c->getDeclaredFields();
                
            for (unsigned int i = 0; i < fields.size(); i++) {
                if (isPoppable(fields[i])) {
                    ObjectHashTable::const_iterator iter =
                        hash.find(fields[i]);
                    IObjectSP data((iter == hash.end())?
                                   IObjectSP(   ): iter->second);
                    // if missing or 'null' then fail if not optional
                    bool realData = data.get() && 
                        !CNull::TYPE->isInstance(data);
                    if (!realData){
                        if (!fields[i]->isOptional()){
                            throw ModelException(routine,
                                                 "Item "+
                                                 fields[i]->getName()+
                                                 " is missing");
                        }
                    } else {
                        try{
                            IObjectSP alteredData;
                            // convert any public representation to private ones
                            alteredData = CObject::convertToPrivateRep(data);
                            // report a good error message if object not of the
                            // right type
                            checkElementType(fields[i], alteredData);
                            // QLib behaviour is to CLONE depending on copyOnPop2Obj setting
                            bool makeCopy = copyOnPop2Obj_ && 
                                !data->getClass()->isPrimitive() &&
                                alteredData.get() == data.get();
                            if (makeCopy){
                                alteredData = IObjectSP(alteredData.clone());
                            }
                            fields[i]->set(object, alteredData);
                        } catch (exception& e){
                            throw ModelException(e, routine, "for field "+
                                                 fields[i]->getName());
                        }
                    }
                }
            }
        } while ((c = c->getSuperClass()) != 0);
        // now invoke object's validation routine
        object->validatePop2Object();
    } catch (exception& e){
        throw ModelException(&e, routine, "for object of type "+
                             type->getName());
    }
}

/** Create an object from a DataDictionary */
IObjectSP CDataDictionary::pop2Object() const{
    try{
        IObjectSP object(type->newInstance());
        populateObject(object.get());
        // we leave the object into its public representation (rather 
        // than any private one)
        return object;
    } catch (exception& e){
        throw ModelException(&e, "CDataDictionary::pop2Object for object"
                             " of type "+type->getName());
    }
}




bool CDataDictionary::isPoppable(CFieldConstSP f){
    return (!Modifier::isTransient(f->getModifiers()));
}

CDataDictionary::ObjectHashTable::iterator 
CDataDictionary::iterationBegin(void){
    return hash.begin();
}

CDataDictionary::ObjectHashTable::iterator 
CDataDictionary::iterationEnd(void){
    return hash.end();
}

static void loadDataDictionary(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spread sheet
    REGISTER(CDataDictionary, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(ITypeConvert);
    IMPLEMENTS(IReadableMap);
    IMPLEMENTS(IWriteableMap);
}

CClassConstSP const CDataDictionary::TYPE =
CClass::registerClassLoadMethod("DataDictionary", 
                                typeid(CDataDictionary), loadDataDictionary);

//// Iterator support
class CDataDictionary::Iterator: public CObject,
                 public virtual IMap::IIterator{
public:
    static CClassConstSP const TYPE;
    // constructor
    Iterator(const CDataDictionarySP& dd):
        CObject(TYPE), iter(dd->iterationBegin()), dd(dd){}

    //// are there any elements left to iterate over
    bool hasMoreElements() const{
        return (iter != dd->iterationEnd());
    }

    /** get the next element and increment the iterator. The
        DataDictionary has the 'public' view of the object so any object
        which are 'private' are converted on the fly (ie on demand) */
    IObjectSP getElement() const{
        IObjectSP elt(iter->second);
        IObjectSP publicElt(CObject::convertToPublicRep(elt));
        if (elt.get() != publicElt.get()){
            // since we've converted it, might as well save it
            iter->second = publicElt;
        }
        return publicElt;
    }

    //// get the key for the current element
    const string& getKey() const{
        return iter->first->getName();
    }
    //// increment the iterator
    void increment(){
        ++iter;
    }
private:
    static void load(CClassSP& clazz){
        REGISTER(Iterator, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IMap::IIterator);
        // no cloning/deserialisation - at least not for now
    }

    // fields
    ObjectHashTable::iterator    iter; // $unregistered
    CDataDictionarySP            dd; // $unregistered
};

  
/** override clone method */
IObject* CDataDictionary::clone() const{
    try{
        CDataDictionarySP copy(new CDataDictionary(type));
        // iterate over all elements
        for (ObjectHashTable::const_iterator myIterator = hash.begin();
             !(myIterator == hash.end());
             ++myIterator){
            CFieldConstSP key = myIterator->first;
            IObjectSP component(myIterator->second.clone());
            copy->hash[key] = component;
        }
        return copy.release();
    } catch (exception& e) {
        throw ModelException(&e, "CDataDictionary::clone");
    }
}

/** Hashes each object inside the data dictionary */
int CDataDictionary::hashCode() const{
    int hCode = (size_t) TYPE ^ (size_t) type;
    for (ObjectHashTable::const_iterator myIterator = hash.begin();
         !(myIterator == hash.end()); ++myIterator){
        hCode ^= myIterator->second->hashCode();
    }
    return hCode;
}

//// local definition of == method to allow comparison of hash tables
static bool operator==(IObjectSP x, IObjectSP y){
    return (x->equalTo(y.get()));
}

/** Compares each object inside the data dictionary */
bool CDataDictionary::equalTo(const IObject* obj) const{
    if (obj == this){
        return true;
    }
    if (!obj || obj->getClass() != TYPE){
        return false;
    }
    const CDataDictionary* dataDict2 = STATIC_CAST(CDataDictionary, obj);
    return (dataDict2->type == type && dataDict2->hash == hash);
}

CClassConstSP const CDataDictionary::Iterator::TYPE = 
CClass::registerClassLoadMethod(
    "DataDictionary::Iterator", typeid(Iterator), load);


//// Builds an iterator
IMap::IIteratorSP CDataDictionary::createIterator(){
    return IIteratorSP(new Iterator(CDataDictionarySP::attachToRef(this)));
}

/** Is this object truly a map ie does toObject()/toMap() return this
    Returns true */
bool CDataDictionary::isTrueMap() const{
    return true;  // as far as DR Interface is concerned
}

class DDCreateAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // three addin parameters
    string           typeName; // gives type of data dictionary to create
    // remaining two parameters are optional
    CStringArraySP   fieldNames; // identifies objects
    ObjectArraySP    objects;    // the objects themselves

    // utility method
    static void fillDataDict(CDataDictionarySP&    dataDict,
                             const CStringArray&   fieldNames,
                             const ObjectArraySP   objects) { // can be null
        
        static const string routine = "fillDataDict";
        for (int i = 0; i < fieldNames.size(); i++){
            if (!fieldNames.empty()){
                if (!objects || i >= objects->size()){
                    throw ModelException(routine, "Object for field "+
                                         fieldNames[i]+" is missing");
                }
                dataDict->put(fieldNames[i], (*objects)[i]);
            }
        }
    }

    /** create a data dictionary */
    static IObjectSP create(DDCreateAddin *params){
        CDataDictionarySP dd(CDataDictionary::create(params->typeName));
        if (params->fieldNames.get()){
            fillDataDict(dd, *params->fieldNames, params->objects);
        }
        return dd;
    }

    DDCreateAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DDCreateAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDCreateAddin);
        FIELD(typeName, "Type of Data Dictionary");
        FIELD(fieldNames, "Identifies the object");
        FIELD_MAKE_OPTIONAL(fieldNames);
        FIELD(objects, "The object itself");
        FIELD_MAKE_OPTIONAL(objects);
        Addin::registerClassObjectMethod("DATA_DICTIONARY",
                                         Addin::UTILITIES,
                                         "Creates a handle to a data "
                                         "dictionary",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)create);
    }

    static IObject* defaultDDCreateAddin(){
        return new DDCreateAddin();
    }
};

CClassConstSP const DDCreateAddin::TYPE = CClass::registerClassLoadMethod(
    "DDCreateAddin", typeid(DDCreateAddin), load);

class DDPutAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // three addin parameter
    CDataDictionarySP  dataDict;  // the data dictionary
    CStringArraySP     fieldNames; // identifies objects
    ObjectArraySP      objects;    // the objects themselves

    /** set an object in a data dictionary */
    static IObjectSP set(DDPutAddin *params){
        DDCreateAddin::fillDataDict(params->dataDict, *params->fieldNames,
                                    params->objects);
        return params->dataDict; // return the now modified data dictionary
    }

    DDPutAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DDPutAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDPutAddin);
        FIELD(dataDict, "Data Dictionary to put object into");
        FIELD(fieldNames, "Identifies the object");
        FIELD(objects, "The object itself");
        Addin::registerInstanceObjectMethod(
            "DATA_DICTIONARY_PUT",
            Addin::UTILITIES,
            "Adds objects to a data dictionary",
            TYPE,
            true,
            Addin::returnHandle,
            (Addin::ObjMethod*)set);
    }

    static IObject* defaultDDPutAddin(){
        return new DDPutAddin();
    }
};

CClassConstSP const DDPutAddin::TYPE = CClass::registerClassLoadMethod(
    "DDPutAddin", typeid(DDPutAddin), load);

class DDPopAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    // one addin parameter
    CDataDictionarySP  dataDict;  // the data dictionary

    /** pop data dictionary to an object */
    static IObjectSP pop(DDPopAddin *params){
        return params->dataDict->pop2Object();
    }

    DDPopAddin(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(DDPopAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultDDPopAddin);
        FIELD(dataDict, "Data Dictionary to convert");
        Addin::registerInstanceObjectMethod("DATA_DICTIONARY_POP",
                                            Addin::UTILITIES,
                                            "Turns a data dictionary into the"
                                            " corresponding object",
                                            TYPE,
                                            true,
                                            Addin::returnHandle,
                                            (Addin::ObjMethod*)pop);
    }

    static IObject* defaultDDPopAddin(){
        return new DDPopAddin();
    }
};

CClassConstSP const DDPopAddin::TYPE = CClass::registerClassLoadMethod(
    "DDPopAddin", typeid(DDPopAddin), load);

  
DRLIB_END_NAMESPACE
