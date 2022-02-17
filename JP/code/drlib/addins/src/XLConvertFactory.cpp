//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertFactory.cpp
//
//   Description : Provide the correct instance of an XLConvert class
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/XLConvertFactory.hpp"
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
// hash function for CClassConstSP
template <class _Key> struct myPtrHash {
    size_t operator()(const _Key p) const {return (size_t)p;}
};
typedef hash_map<CClassConstSP, refCountPtr<XLConvert>, 
    myPtrHash<CClassConstSP> > XLConvertHash;

class XLConvertFactoryHelper{
public:
    static XLConvertHash   hashTableForObjects; // ie not arrays
    static XLConvertHash   hashTableForArrays;
    static bool            registered;

    //// Store XLConvert instance in supplied Hashtable if not there already
    static bool add(XLConvertHash& hashTable,
                    CClassConstSP  clazz,
                    XLConvert*     instance){
        // since we give out XLConvert pointers, can't have them go out of scope
        if (hashTable.find(clazz) == hashTable.end()){
            hashTable[clazz] = refCountPtr<XLConvert>(instance);
        }
        return true;
    }
};

XLConvertHash XLConvertFactoryHelper::hashTableForObjects; // static field
XLConvertHash XLConvertFactoryHelper::hashTableForArrays; // static field


// force linker to include all classes
bool XLConvertFactoryHelper::registered = true;

/** Create an instance of an XLConvert which corresponds to the
    given class */
const XLConvert* XLConvertFactory::create(CClassConstSP clazz){
    /** look for exact match based upon object type, if no match found then
        repeat for parent type (and so on) */
    CClassConstSP c = clazz;
    bool isArray = c->isArray(); // work on type of components if an array
    if (isArray){
        c = c->getComponentType();
    }
    // to do - make IArray be flagged as an array type with element type IObject
    // and remove this 'else' clause
    else if (c == IArray::TYPE){
        isArray = true;
        c = IObject::TYPE;
    }
        
    const XLConvert* method;
    while (!(method = lookUp(c, isArray))){
        // Recall getSuperClass gives null for interface classes
        if (!(c = c->getSuperClass()) || 
            c == IObject::TYPE /* CObject has IObject as its superclass
                                  which is looking increasingly wrong */){
            // note that we currently don't support looking up methods for
            // objects of interface types - so we'll manually cope with
            // interfaces derived from IObject - possibly we'll never need
            // anything more (this is needed for parameters eg of type
            // IPastValues)
            if (IObject::TYPE->isAssignableFrom(clazz)){
                method = lookUp(IObject::TYPE, isArray);
                if (method){
                    return method;
                }
            }
            // if you get this error when clazz is an interface then you need
            // the interface to derive from IObject
            throw ModelException("XLConvertFactory::create", "No conversion"
                                 " method available for type "+
                                 clazz->getName());
        }
    }
    return method;
}

/** Does look up for exact type. Returns null if no match */
const XLConvert* XLConvertFactory::lookUp(
    CClassConstSP clazz,
    bool          clazzIsArrayComponentType){
    const XLConvertHash& hashTable = clazzIsArrayComponentType? 
        XLConvertFactoryHelper::hashTableForArrays:
        XLConvertFactoryHelper::hashTableForObjects;
    XLConvertHash::const_iterator iterator = hashTable.find(clazz);
    return (iterator == hashTable.end())? 0: iterator->second.get();
}

/** Supply instance of XLConvert corresponding to given class which should not
    be an array.
    NB XLConvertFactory takes ownership of memory */
bool XLConvertFactory::registerXLObjectConvert(CClassConstSP clazz,
                                               XLConvert*    instance){
    return XLConvertFactoryHelper::add(
        XLConvertFactoryHelper::hashTableForObjects, clazz, instance);
}
/** Supply instance of XLConvert corresponding to arrays whose elements are of,
    or are derived from, given class. 
    NB XLConvertFactory takes ownership of memory */
bool XLConvertFactory::registerXLArrayConvert(CClassConstSP  arrayComponentType,
                                               XLConvert*    instance){
    return XLConvertFactoryHelper::add(
        XLConvertFactoryHelper::hashTableForArrays,
        arrayComponentType,
        instance);
}

DRLIB_END_NAMESPACE
