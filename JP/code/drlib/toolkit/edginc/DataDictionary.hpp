
#ifndef EDR_DATADICTIONARY_HPP
#define EDR_DATADICTIONARY_HPP

#include <map>
#include "edginc/Object.hpp"
#include "edginc/Field.hpp"
#include "edginc/TypeConvert.hpp"
#include "edginc/WriteableMap.hpp"
#include "edginc/ReadableMap.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE
class CDataDictionary;
// typedef for smart pointers to DataDictionary
typedef smartConstPtr<CDataDictionary> CDataDictionaryConstSP;
typedef smartPtr<CDataDictionary> CDataDictionarySP;

/** A DataDictionary is a form of hash table. Its primary ability is to
    hold values for fields of a class ie the value of each field of a 
    class TOOLKIT_DLL can be stored in the hash table using the field's name as the
    key in the hash table.

    The two key methods are pop2DataDict and populateObject. These provide
    the ability to turn an object (derived from IObject) into a 
    DataDictionary and vice versa.
*/
class TOOLKIT_DLL CDataDictionary: public CObject,
                       virtual public ITypeConvert,
                       virtual public IReadableMap,
                       virtual public IWriteableMap{
    friend class shutTheCompilerUp;
public:
    typedef map<CFieldConstSP, IObjectSP> ObjectHashTable;

    static CClassConstSP const TYPE;

    ~CDataDictionary();

    /** override clone method */
    virtual IObject* clone() const;

    /** Hashes each object inside the data dictionary */
    virtual int hashCode() const;

    /** Compares each object inside the data dictionary */
    virtual bool equalTo(const IObject* obj) const;

    /** set flag indicating whether copies of the objects in the
        data dictionary should be used when building an object in
        {@link #pop2Object}. Default is true */
    void copyOnPop2Obj(bool copyOnPop2Obj) throw();

    /** Create a DataDictionary given object's name, such as
        "BorrowCurve". This method allows the creation of public
        types only. */
    static CDataDictionary* create(const string& typekey);

    /** Create an empty DataDictionary given an object's class. NB This
        method allows the creation of protected and private types. Use
        createPublic(string) for use by EAS or the spreadsheet */
    static CDataDictionary* create(const CClassConstSP& clazz);

    /** What type of DataDictionary is this ? */
    CClassConstSP getType() const throw();

    //// Builds an iterator
    virtual IIteratorSP createIterator();

    /** Is this object truly a map ie does toObject()/toMap() return this
        Returns false */
    virtual bool isTrueMap() const;

    /** Add an object to a DataDictionary - works its way up the
        inheritance chain, matching the first data member whose name
        matches the key * @param key Data member of clazz * @param
        value Object to add.  */
    void put(const string& key, const IObjectSP& value);
    
    /** Same as above put will clone value if refCount = 0 */
    void put(const string& key, IObject* value);

    //// nicer versions of put for the primitives
    void put(const string& key, bool b);
    void put(const string& key, double d);
    void put(const string& key, int i);
    void put(const string& key, const string& s);
    void put(const string& key, const char* s);

    /** Add an object to a DataDictionary - uses fully qualified
        class TOOLKIT_DLL and member names so there is no clash of member names
        with any superclass * @param clazz The fully qualified class
        name * @param member Data member of clazz * @param value
        Object to add */
    void put(const string& clazz, const string& member, 
             IObjectSP& value);
    
    /** Does a DataDictionary contain an object named by 'key'?
        Works from base class upwards looking for a field that 
        matches the key - be wary of name clashes with 
        inherited fields. If this returns true, get() will work. */
    bool contains(const string& key) const;

    /** Get an object out of a DataDictionary - works from base
        class TOOLKIT_DLL upwards looking for a field that matches the key - be
        wary of name clashes with inherited fields. Returns 0 if the
        object is not found */
    IObjectSP get(const string& key) const;

    /** Get an object out of a DataDictionary - uses fully
        qualified class and member names so there is no clash of
        member names with any superclass * @param clazz The fully
        qualified class name * @param member Data member of clazz */
    IObjectSP get(const string& clazz, const string& member) const;

    /** Get a field out of a DataDictionary */
    IObjectSP get(CFieldConstSP& key) const;
    
    /** Create a DataDictionary containing the data that makes up
        an object. If the object implements the
        IPrivateObject interface, then before being converted to a 
        data dictionary, the object is
        converted to an IPublicObject using the toPublicObject() method*/
    static CDataDictionary* pop2DataDict(IObjectSP obj);

    /** Create an object from a DataDictionary.
        All mandatory fields must be present. If the resulting object
        implements the IPublicObject interface, then it is converted into an
        IPrivateObject object using the toPrivateObject() method */
    IObjectSP pop2Object() const;

    /** Populate an existing object with the contents of a data dictionary.
        As per pop2Object, all mandatory fields must be present.
        Really only for use by CObject::xmlImport */
    void populateObject(IObject* obj) const;
    
    /** Converts this object to an instance of the
        requiredType. Throws an exception if a conversion to the
        required Type is not supported. Here we turn the data
        dictionary into its object equivalent and then invoke
        CObject::checkType */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const;

    /** Returns an stl iterator for the <key, value> pairs in the hash table
        starting at the start of the hash table */
    ObjectHashTable::iterator iterationBegin(void);
    /** Returns an stl iterator for the <key, value> pairs in the hash table
        starting at the end of the hash table */
    ObjectHashTable::iterator iterationEnd(void);

private:
    class Iterator;
    friend class Iterator;
    CDataDictionary(const CClassConstSP& type);
    CDataDictionary(const CDataDictionary& rhs);
    CDataDictionary& operator=(const CDataDictionary& rhs);
    static CObject* defaultDataDictionary();
    void checkElementType(CFieldConstSP field, IObjectSP& value) const;
    void put(CFieldConstSP field, IObjectSP value);

    CClassConstSP type; // the type of the object we're representing $unregistered
    mutable ObjectHashTable hash; // contains the data for the object $unregistered

    /** if copyOnPop2Obj is true then when the data dict is
        converted to the object each component is cloned before being
        inserted into the object */
    bool   copyOnPop2Obj_; // $unregistered

    static bool isPoppable(CFieldConstSP f);

};

#ifndef QLIB_DATADICTIONARY_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<CDataDictionary>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<CDataDictionary>);
#endif

DRLIB_END_NAMESPACE
#endif

        
