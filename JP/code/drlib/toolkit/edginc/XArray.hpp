//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XArray.hpp
//
//   Description : Wrapper for external [DR interface] objects which are
//                 arrays
//
//   Author      : Mark A Robson
//
//   Date        : 2 Dec 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XARRAY_HPP
#define EDR_XARRAY_HPP
#include "edginc/XObject.hpp"

DRLIB_BEGIN_NAMESPACE
//// Wrapper for external [DR interface] objects which are arrays.
class TOOLKIT_DLL XArray: public XObject,
              public virtual IArray{
public:
    static CClassConstSP const TYPE;
    virtual ~XArray();

    /** Creates an array of specified length with elements of specified type.
        eltTypeName can be empty in which case null is used in the equivalent
        DR interace function */
    static XArraySP create(const string& eltTypeName, 
                           int           numElements,
                           DRService*    svc);

    /** Creates an array from supplied vector of XObjects using specified type.
        The vector must not be empty */
    static XArraySP create(const string&            eltTypeName, 
                           const vector<XObjectSP>& elts);

    /** Simple constructor. */
    explicit XArray(const MyDRObjectInterfaceSP& theArray);

    /** Simple constructor. Takes ownership of memory */
    explicit XArray(DRArray theArray);

    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    /** If we own the DRObject then nothing is done (NB May need to review this
        since an XArray can be modified). Otherwise
        copy the object - by building empty array and getting/setting
        (which is a bit weak but there is no clone method) */
    IObject* clone() const;

    /** Override default CObject implementation */
    virtual void xWrite(const string& tag, Writer* writer) const;

    /** Override default CObject implementation */
    virtual void xImport(Reader::Node* elem, Reader* reader, DRService* svc);

    /** Override default CObject implementation */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix,
                             ostream&      stream) const;

    /** CombinableResult interface: scale by factor x */
    virtual void scale(double x);

    /** CombinableResult interface: 
        add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Returns length of array */
    virtual int getLength() const;

    /** Returns the value of the indexed component in the
        specified array object. */
    virtual IObjectConstSP get(int index) const;

    /** Non const version of above */
    virtual IObjectSP get(int index);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value */
    virtual void set(int index, IObjectSP value);

    /** Appends the supplied object to the end of supplied array. The
        length of the array is increased by 1. (aka push_back). 
        Note this implementation is particularly slow due to lack of
        append method in DRInterface */
    virtual void append(IObjectSP value);

    /** Create a clone of this object by recursively cloning every component
        of this object. This is used for regression testing and should not
        be used in general (unnecessarily slow) */
    virtual XObjectSP recursiveClone() const;

    /** Returns length of array */
    int size() const;

    /** Get's the n'th element. */
    IDRObjectSP getElt(int index) const;

    /**  Set's the n'th element. */
    void setElt(int index, IDRObjectSP obj);

    /** Returns the type of elements of this array.  An empty string
        is returned if the corresponding DR interface method returns
        null */
    string elementType() const;

    /** Returns the XClass representing the elements of this array. Beware
        returns null if the corresponding DRI function returns null (which
        it's not supposed to but Enhanced Magnet still does) */
    XClassConstSP elementClass() const;
private:
    XArray(const XArray& rhs); // don't use
    XArray& operator=(const XArray& rhs); // don't use
    static void load(CClassSP& clazz);
    XArray();
};


DRLIB_END_NAMESPACE
#endif
