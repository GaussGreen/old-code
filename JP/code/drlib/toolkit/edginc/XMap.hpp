//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XMap.hpp
//
//   Description : Wrapper for external [DR interface] maps
//
//   Author      : Mark A Robson
//
//   Date        : 14 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XMAP_HPP
#define EDR_XMAP_HPP
#include "edginc/XObject.hpp"

DRLIB_BEGIN_NAMESPACE
//// Wrapper for external [DR interface] objects
class TOOLKIT_DLL XMap: public XObject{
public:
    static CClassConstSP const TYPE;
    virtual ~XMap();

    /** Simple constructor. */
    explicit XMap(const MyDRObjectInterfaceSP& map);
    
    /** Simple constructor. Takes ownership of memory */
    explicit XMap(DRMap map);

    /** Creates a map of the specified type using the supplied service */
    static XMapSP create(const char* typeName, DRService* svc);

    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    /** If we own the DRObject then nothing is done (NB May need to
        review this since an XMap can be modified). Otherwise copy
        the object - by converting to object and back again (which is
        a bit weak but there is no clone method) */
    virtual IObject* clone() const;

    /** CombinableResult interface: scale by factor x */
    virtual void scale(double x);

    /** CombinableResult interface: 
        add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Create a clone of this object by recursively cloning every component
        of this object. This is used for regression testing and should not
        be used in general (unnecessarily slow) */
    virtual XObjectSP recursiveClone() const;

    //// write elements of map to writer
    void xWriteElts(Writer* writer) const;

    /**  read elements of map from reader */
    void xImportElts(Reader::Node* elem, Reader* reader);

    /** Override default CObject implementation */
    void outputWrite(const string& linePrefix,
                     const string& prefix,
                     ostream&      stream) const;

    /** Returns the type of the external object that this map corresponds to */
    string getXObjectClass() const;

    /** Gets the specified item from the external map */
    IDRObjectSP getItem(const string& fieldName) const;

    /** Gets the specified item from the external map */
    IDRObjectSP getItem(const char* fieldName) const;

    /** To make life a bit easier, get methods on the map where you know what
        type you're expecting */
    bool getBool(const string& fieldName) const; 
    int getInt(const string& fieldName) const;
    double getDouble(const string& fieldName) const;
    string getString(const string& fieldName) const;
    XObjectSP getXObject(const string& fieldName) const;
    XArraySP getXArray(const string& fieldName) const;

    /** Sets the specified item from the external map */
    void setItem(const string& fieldName, const IDRObjectSP& obj) const;

    /** Sets the specified item from the external map */
    void setItem(const char* fieldName, const IDRObjectSP& obj) const;

    /** Turns map into object */
    XObjectSP toObject() const;
    //// for iterating over maps
    class TOOLKIT_DLL Iterator{
    public:
        ~Iterator();
        //// move to first/next element. Returns false if none available
        bool nextElement();
        const string& key() const; // goes out of scope once nextElement called
        IDRObjectSP value() const; // goes out of scope once nextElement called
    private:
        friend class XMap;
        Iterator(const XMap* xMap);
        class Imp;
        Imp*   my; // hide implementation
    };
    Iterator createIterator() const;
protected:
    /** Override default XObject implementation. */
    virtual void xWrite(const string& tag, Writer* writer) const;
    /** Override default XObject implementation */
    virtual void xImport(Reader::Node* elem, Reader* reader, DRService* svc);
private:
    friend class Iterator::Imp;
    friend class Iterator;
    XMap(const XMap& rhs); // don't use
    XMap& operator=(const XMap& rhs); // don't use
    static void load(CClassSP& clazz);
    XMap();
    static IObject* defaultConstructor();
    void getItem(const char* fieldName, DRValue& drValue, int reqdType) const;
};

DRLIB_END_NAMESPACE
#endif
