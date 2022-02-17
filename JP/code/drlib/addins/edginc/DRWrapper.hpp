//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DRWrapper.hpp
//
//   Description : It's DR Wrapper
//
//   Author      : Andrew J Swain
//
//   Date        : 19 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef DRWRAPPER_HPP
#define DRWRAPPER_HPP
#include <map>
#include "edginc/WriteableMap.hpp"
#include "edginc/Array.hpp"

using namespace std;


DRLIB_BEGIN_NAMESPACE

typedef map<string, IObjectSP> WrapperHashTable;


/** DR Wrapper - captures instrument data for types not in production library */
class ADDINS_DLL DRWrapper: public CObject,
                 public virtual IWriteableMap {
public:
    static CClassConstSP const TYPE;

    /** Create a DRWrapper given object's name, such as "MyNewProduct" */
    static DRWrapper* create(const string& typekey);

    /** Add an object to a DRWrapper */
    virtual void put(const string& key, const IObjectSP& value);

    /** Create an object from a DRWrapper.
        All mandatory fields must be present. If the resulting object
        implements the IPublicObject interface, then it is converted into an
        IPrivateObject object using the toPrivateObject() method */
    IObjectSP pop2Object() const;

    /** check the executable name matches the wrapper */
    void validateExeName() const;

    /** tag for recording wrapper type */
    static const string WRAPPER_TYPE;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from reader */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** override clone method */
    virtual IObject* clone() const;

    /** so the class gets linked in properly */
    static bool load();

    /** build an iterator */
    virtual IIteratorSP createIterator();

    /** Is this object truly a map ie does toObject()/toMap() return this */
    virtual bool isTrueMap() const;

    // If object is a drwrapper convert to equivalent object, otherwise returns
    // original
    static IObjectSP drWrapperToObject(const IObjectSP& object,
                                       CClassConstSP    desiredType);
private:
    friend class DRWrapperHelper;
    friend class DRWrapperProxy;
    DRWrapper(const string& type);  
    /** for reflection */
    DRWrapper();

    DRWrapper(const DRWrapper& rhs);
    DRWrapper& operator=(const DRWrapper& rhs);

    string           type;  // what class does this represent ?
    WrapperHashTable hash;  // contains the data for the object $unregistered
};

typedef smartPtr<DRWrapper> DRWrapperSP;
typedef array<DRWrapperSP, DRWrapper> DRWrapperArray;
typedef smartPtr<DRWrapperArray> DRWrapperArraySP;

DRLIB_END_NAMESPACE
#endif
