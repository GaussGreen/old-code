
#include "edginc/config.hpp"
#define QLIB_HASHTABLE_CPP
#include "edginc/Hashtable.hpp"
#include "edginc/Class.hpp"
#include "edginc/Writer.hpp"
#include ext_hash_map
DRLIB_BEGIN_NAMESPACE

INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<Hashtable>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<Hashtable>);

Hashtable::KeyEnumeration::~KeyEnumeration(){}
Hashtable::KeyEnumeration::KeyEnumeration(){}
Hashtable::ValueEnumeration::~ValueEnumeration(){}
Hashtable::ValueEnumeration::ValueEnumeration(){}



struct MyStringHash {
    size_t operator()(const string& str) const {
        return (hash_string(str.c_str()));
    }
};


typedef hash_map<string, IObjectSP, MyStringHash> HashData;

class HashtableHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Hashtable, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IReadableMap);
        IMPLEMENTS(IWriteableMap);
        EMPTY_SHELL_METHOD(defaultHashtable);
        // no fields
    }

    static IObject* defaultHashtable(){
        return new Hashtable();
    }
    /* by making all references to hash map internal we stop all clients
       from sucking in all the stl headers */
    HashData hash;

};

/** Helper class for xml read/write for hash table of 
    <string, IObjectSP> */
class Hashtable::Entry: public CObject{
public:
    static CClassConstSP const TYPE;
    string             key;
    IObjectSP          value;
    static const string XML_TAG;
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(Entry, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultEntry);
        FIELD(key, "key");
        FIELD(value, "value");
    }

    static IObject* defaultEntry(){
        return new Entry();
    }
    // for reflection
    Entry(): CObject(TYPE){}

    Entry(const string&    key, 
          const IObjectSP& value):
        CObject(TYPE), key(key), value(value){}
};
const string Hashtable::Entry::XML_TAG = "Entry";

CClassConstSP const Hashtable::Entry::TYPE =
CClass::registerClassLoadMethod("Hashtable::Entry", typeid(Entry), load);


class HashKeyEnumeration: public Hashtable::KeyEnumeration{
public:
    HashtableHelper *data; // a reference
    HashData::const_iterator iter;
    HashKeyEnumeration(HashtableHelper *data): 
        data(data), iter(data->hash.begin()){}
    
    /** Tests if this enumeration contains more elements. */
    virtual bool hasMoreElements() const{
        return !(iter == data->hash.end());
    }
    
    /** Returns the next element of this enumeration if this
        enumeration object has at least one more element to provide. */
    virtual const string& nextElement(){
        return (iter++)->first;
    }
};

class HashValueEnumeration: public Hashtable::ValueEnumeration{
public:
    HashtableHelper *data; // a reference
    HashData::const_iterator iter;
    HashValueEnumeration(HashtableHelper *data): 
        data(data), iter(data->hash.begin()){}
    /** Tests if this enumeration contains more elements. */
    virtual bool hasMoreElements() const{
        return !(iter == data->hash.end());
    }
    
    /** Returns the next element of this enumeration if this
        enumeration object has at least one more element to provide. */
    virtual IObjectSP nextElement(){
        return (iter++)->second;
    }
};

CClassConstSP const Hashtable::TYPE = CClass::registerClassLoadMethod(
    "Hashtable", typeid(Hashtable), HashtableHelper::load);

Hashtable::Hashtable(): CObject(TYPE){
    data = new HashtableHelper();
}

Hashtable::~Hashtable(){
    delete data;
}

/** write elements out in XML format */
void Hashtable::writeElts(Writer* writer) const {
    for (HashData::const_iterator myIterator = data->hash.begin();
         myIterator != data->hash.end(); ++myIterator){
        Entry entry(myIterator->first, myIterator->second);
        entry.write(Entry::XML_TAG, writer);
    }
}

/** write object out in XML format */
void Hashtable::write(const string& tag, Writer* writer) const {
    try {
        IObjectConstSP obj(writer->objectStart(tag, "", this, true)); 
        if (obj.get()){
            writeElts(writer);
        }
        writer->objectEnd(tag, this);
    }
    catch (exception& e){
        throw ModelException(&e, "Hashtable::write");
    }
}

/** populate an empty object from XML description */
void Hashtable::import(Reader::Node* elem, Reader* reader){
    const static string method = "Hashtable::import";
    try{
        Reader::NodeListSP children(elem->children());
        for (unsigned int i = 0; i < children->size(); i++) {
            Reader::Node* node = (*children)[i].get();
            // read in entry
            IObjectSP entryObj(reader->read(node));
            Entry* wrapperEntry = dynamic_cast<Entry*>(entryObj.get());
            if (!wrapperEntry){
                throw ModelException(method, "Invalid object found");
            }
            put(wrapperEntry->key, wrapperEntry->value);
        }
    } 
    catch (exception& e){
        throw ModelException(e, method);
    }
}

/** write object out in testcmp format */
void Hashtable::outputWrite(const string& linePrefix,
                            const string& prefix, ostream& stream) const {
    for (HashData::const_iterator myIterator = data->hash.begin();
         myIterator != data->hash.end(); ++myIterator){
        myIterator->second->outputWrite(linePrefix, prefix.empty()?
                                        myIterator->first:
                                        (prefix + "_" + myIterator->first),
                                        stream);
    }
}


/** creates deep copy of each element in hash table */
IObject* Hashtable::clone() const{
    try{
        // support for derived types - build object of same type/size etc
        HashtableSP clone(&dynamic_cast<Hashtable&>(*getClass()->
                                                    newInstance()));
        clone->data->hash = data->hash; // shallow copy
        // then iterate through turning into deep copies
        for (HashData::iterator myIterator = clone->data->hash.begin();
             !(myIterator == clone->data->hash.end());
             ++myIterator){
            IObjectSP eltClone(myIterator->second->clone());
            myIterator->second = eltClone;
        }
        return clone.release();
    } catch (exception& e){
        throw ModelException(e, "Hashtable::clone");
    }
}

/** Hashes each object inside the hash table */
int Hashtable::hashCode() const{
    int hCode = (size_t) TYPE;
    for (HashData::const_iterator myIterator = data->hash.begin();
         !(myIterator == data->hash.end()); ++myIterator){
        hCode ^= myIterator->second->hashCode();
    }
    return hCode;
}

//// local definition of == method to allow comparison of hash tables
static bool operator==(IObjectSP x, IObjectSP y){
    return (x->equalTo(y.get()));
}

/** Compares each object inside the hash table */
bool Hashtable::equalTo(const IObject* obj) const{
    if (obj == this){
        return true;
    }
    if (!obj || obj->getClass() != TYPE){
        return false;
    }
    const Hashtable* hashTable2 = STATIC_CAST(Hashtable, obj);
    return (hashTable2->data->hash == data->hash);
}

/** Clears this hashtable so that it contains no keys. */
void Hashtable::clear(){
    data->hash.clear();
}

/** Tests if the specified object is a key in this hashtable. */
bool Hashtable::containsKey(const string& key) const{
    HashData::const_iterator iter = data->hash.find(key);
    return (iter != data->hash.end());
}

#if 0
/** Returns true if this Hashtable maps one or more keys to this value. */
bool Hashtable::containsValue(const IObjectSP& value) const;
#endif

/** Returns an enumeration of the values in this hashtable. */
Hashtable::ValueEnumeration* Hashtable::elements() const{
    return new HashValueEnumeration(data);
}

/** Returns the value to which the specified key is mapped in this
    hashtable. */
IObjectSP Hashtable::get(const string& key) const{
    HashData::const_iterator iter = data->hash.find(key);
    if (iter == data->hash.end()){
        return IObjectSP();
    }
    return iter->second;
}


/** Tests if this hashtable maps no keys to values. */
bool Hashtable::isEmpty() const{
    return data->hash.empty();
}    

/** Returns an enumeration of the keys in this hashtable.*/
Hashtable::KeyEnumeration* Hashtable::keys() const{
    return new HashKeyEnumeration(data);
}

/** Maps the specified key to the specified value in this
    hashtable.  */
void Hashtable::put(const string& key, const IObjectSP& value){
    if (!value){
        throw ModelException("Null value supplied", "Hashtable::put");
    }
    pair<HashData::iterator, bool> currentVal = 
        data->hash.insert(HashData::value_type(key, value));

    // store the new one
    currentVal.first->second = value; // line unnecessary?
}

/** Removes the key (and its corresponding value) from this hashtable. */
IObjectSP Hashtable::remove(const string& key){
    HashData::iterator iter = data->hash.find(key);
    if (iter == data->hash.end()){
        // key didn't exist;
        return IObjectSP();
    }
    // save the existing object
    IObjectSP returnObj(iter->second);
    // remove it
    data->hash.erase(iter);
    return returnObj;
}

/** Returns the number of keys in this hashtable. */
int Hashtable::size() const{
    return (int)(data->hash.size());
}

class Hashtable::Iterator: public CObject,
                 public virtual IMap::IIterator{
public:
    static CClassConstSP const TYPE;
    //// are there any elements left to iterate over
    virtual bool hasMoreElements() const{
        return (iter != hashtable->data->hash.end());
    }
    //// get the key for the current element
    virtual const string& getKey() const{
        return iter->first;
    }
    //// get the current element (returns a reference)
    virtual IObjectSP getElement() const{
        return iter->second;
    }
    //// increment the iterator
    virtual void increment(){
        ++iter;
    }
    Iterator(HashtableSP hashtable): 
        CObject(TYPE), hashtable(hashtable), 
        iter(hashtable->data->hash.begin()) {}
private:
    static void load(CClassSP& clazz){
        REGISTER(Iterator, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IMap::IIterator);
        // no cloning/deserialisation - at least not for now
    }
    HashtableSP              hashtable; // $unregistered
    HashData::const_iterator iter; // $unregistered
};

CClassConstSP const Hashtable::Iterator::TYPE = 
CClass::registerClassLoadMethod(
    "Hashtable::Iterator", typeid(Iterator), load);

//// Builds an iterator
IMap::IIteratorSP Hashtable::createIterator(){
    return IIteratorSP(new Iterator(HashtableSP::attachToRef(this)));
}

/** Is this object truly a map ie does toObject()/toMap() return this
    Returns true */
bool Hashtable::isTrueMap() const{
    return true;
}


DRLIB_END_NAMESPACE
