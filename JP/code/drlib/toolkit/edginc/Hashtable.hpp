
#ifndef EDG_HASH_TABLE_H
#define EDG_HASH_TABLE_H
#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/WriteableMap.hpp"
#include "edginc/ReadableMap.hpp"

DRLIB_BEGIN_NAMESPACE
class HashtableHelper;

// stick this in here for now (ideally all objects would have a hash function)
// hash function for strings
TOOLKIT_DLL size_t hash_string(const char* __s);
TOOLKIT_DLL size_t hash_string(const string& str);


#if !defined(DEBUG) || defined(QLIB_HASHTABLE_CPP)
OPT_INLINE size_t hash_string(const char* __s){
  unsigned long __h = 0; 
  for ( ; *__s; ++__s)
    __h = 5*__h + *__s;
  
  return size_t(__h);
}

OPT_INLINE size_t hash_string(const string& str) {
    return (hash_string(str.c_str()));
}
#endif
    
/** Provides hash table functionality where the hash is a string
    within typed object environment. This class requires further work to
    make it fully useful. In particular, the xml read/write and the clone
    methods need to be overridden. A useful step would be to use the idea
    of an element - see riskmgr/src/Results.cpp */
class TOOLKIT_DLL Hashtable: public CObject,
                 virtual public IReadableMap,
                 virtual public IWriteableMap{
public:
    /** hash class for use in hashmap's that use strings */
    struct StringHash {
        size_t operator()(const string& str) const {
            return (hash_string(str.c_str()));
        }
    };
    
    /** Should rename this iterator and have functions that return the 
        current key and current object */
    /* The general idea is that their is a single Enumeration interface -
        however strings aren't objects so I've created two rather than go
        for the template option */
    /** Keys in this hash table are strings. This is a enumeration class
        specfically for strings. It provides a view onto the
        keys of a hash table*/
    class TOOLKIT_DLL KeyEnumeration{
    public:
        virtual ~KeyEnumeration();
        /** Tests if this enumeration contains more elements. */
        virtual bool hasMoreElements() const = 0;
          
        /** Returns the next element of this enumeration if this
            enumeration object has at least one more element to provide. */
        virtual const string& nextElement() = 0;
    protected:
        KeyEnumeration();
    private:
        KeyEnumeration(const KeyEnumeration &rhs);
        KeyEnumeration& operator=(const KeyEnumeration& rhs);
        
    };
    typedef refCountPtr<KeyEnumeration> KeyEnumerationSP;
    /** An enumeration class for IObjects. It provides a view onto the
        elements of a hash table */
    class TOOLKIT_DLL ValueEnumeration{
    public:
        virtual ~ValueEnumeration();
        /** Tests if this enumeration contains more elements. */
        virtual bool hasMoreElements() const = 0;
          
        /** Returns the next element of this enumeration if this
            enumeration object has at least one more element to provide. */
        virtual IObjectSP nextElement() = 0;
    protected:
        ValueEnumeration();
    private:
        ValueEnumeration(const ValueEnumeration &rhs);
        ValueEnumeration& operator=(const ValueEnumeration& rhs);
    };
    typedef refCountPtr<ValueEnumeration> ValueEnumerationSP;
        
    static CClassConstSP const TYPE;
    /** create an empty hash table of default size */
    Hashtable();
    virtual ~Hashtable();
    /* need to override several default methods */
    /** creates deep copy of each element in hash table */
    virtual IObject* clone() const;

    /** Hashes each object inside the hash table */
    virtual int hashCode() const;

    /** Compares each object inside the hash table */
    virtual bool equalTo(const IObject* obj) const;

    /** write object out */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from a Reader */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** write object out in testcmp format */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    //// Builds an iterator
    virtual IIteratorSP createIterator();
    /** Is this object truly a map ie does toObject()/toMap() return this
        Returns true */
    virtual bool isTrueMap() const;

    /** Clears this hashtable so that it contains no keys. */
    void clear();

    /** Tests if the specified object is a key in this hashtable. */
    bool containsKey(const string& key) const;

    /** Returns true if this Hashtable maps one or more keys to this Key. */
    bool containsValue(const IObjectSP& value) const;

    /** Returns an enumeration of the values in this hashtable. */
    ValueEnumeration* elements() const ;

    /** Returns the value to which the specified key is mapped in this
        hashtable. Note no conversion (eg from private to public is
        performed)*/
    IObjectSP get(const string& key) const;

    /** Tests if this hashtable maps no keys to values. */
    bool isEmpty() const;

    /** Returns an enumeration of the keys in this hashtable.*/
    KeyEnumeration* keys() const;

    /** Maps the specified key to the specified value in this
        hashtable. Note no conversion (eg from public to private is
        performed) */
    virtual void put(const string& key, const IObjectSP& value);

    /** Removes the key (and its corresponding value) from this hashtable. */
    IObjectSP remove(const string& key);

    /** Returns the number of keys in this hashtable. */
    int size() const;

protected:
    /** write elements out in XML format */
    void writeElts(Writer* writer) const;

private:
    friend class HashtableHelper;
    class Entry;
    class Iterator;
    friend class Iterator;
    HashtableHelper *data; // $unregistered
};

typedef Hashtable CHashtable;
typedef smartPtr<Hashtable> HashtableSP;
typedef smartConstPtr<Hashtable> HashtableConstSP;
typedef HashtableSP CHashtableSP;
typedef HashtableConstSP CHashtableConstSP;

#ifndef QLIB_HASHTABLE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Hashtable>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Hashtable>);
#endif

DRLIB_END_NAMESPACE

#endif
