//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XObject.hpp
//
//   Description : Wrapper for external [DR interface] objects
//
//   Author      : Mark A Robson
//
//   Date        : 12 Nov 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XOBJECT_HPP
#define EDR_XOBJECT_HPP
#include "edginc/Array.hpp"
#include "edginc/DRObject.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/CombinableResult.hpp"

struct DRObjectInterface_;
typedef struct DRObjectInterface_ DRObjectInterface;
typedef DRObjectInterface* DRObject;
typedef DRObject DRMap;
typedef DRObject DRArray;
typedef DRObject DRMatrix;
typedef const char* DRString;
typedef struct DRValue_  DRValue;
typedef struct DRService_ DRService;
DRLIB_BEGIN_NAMESPACE
class XObject;
typedef smartPtr<XObject> XObjectSP;
typedef smartConstPtr<XObject> XObjectConstSP;
class XMap;
typedef smartPtr<XMap> XMapSP;
typedef smartConstPtr<XMap> XMapConstSP;
class XArray;
typedef smartPtr<XArray> XArraySP;
typedef smartConstPtr<XArray> XArrayConstSP;
class XMatrix;
typedef smartPtr<XMatrix> XMatrixSP;
typedef smartConstPtr<XMatrix> XMatrixConstSP;
class XClass;
typedef const XClass* XClassConstSP;
typedef XClass* XClassSP;
typedef vector<XClassConstSP> XClassVec;
#ifndef QLIB_XOBJECT_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<XObject>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<XMap>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<XArray>);
#else
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<XObject>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<XMap>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<XArray>);
#endif
//// same as DRObjectInterface but with destructor. Means we can put MyDRObjects
//// in refCountPtr (cast on the way in/out)
struct MyDRObjectInterface{
    DRService *svc;
    ~MyDRObjectInterface();
    void operator delete(void *ptr); // does nothing
};

typedef MyDRObjectInterface* MyDRObject;
typedef attic::oldCountPtr<MyDRObjectInterface>::type MyDRObjectInterfaceSP;
//typedef auto_ptr<MyDRObjectInterface> MyDRObjectInterfaceSP;

//// Wrapper for external [DR interface] objects
class TOOLKIT_DLL XObject: public CObject,
               public virtual IDRObject,
               public virtual IRegressionTest,
               public virtual ClientRunnable,
               public virtual CombinableResult{
public:
    static CClassConstSP const TYPE;
    // the type of the wrapper object used for serialisation (work around for
    // missing functionality in DRI)
    static CClassConstSP const WRAPPER_TYPE;

    virtual ~XObject();

    static const string XTYPE_ATTRIBUTE; // xml attribute id for getXClass()

    /** Are these objects equal (ie contain the same data) */
    virtual bool equals(IDRObjectConstSP drObj) const;

    /** If we own the DRObject then we create another XObject pointing to
        the same DRObject. Otherwise copy the object - by turning into map
        and back again */
    virtual IObject* clone() const;

    /** Override default CObject implementation */
    virtual void write(const string& tag, Writer* writer) const;

    /** Override default CObject implementation */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** Override default CObject implementation */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix,
                             ostream&      stream) const;

    /** ClientRunnable interface */
    virtual IObjectSP run();
    /** IRegressionTest interface */
    virtual IObjectSP runTest() const;

    /** Create a clone of this object by recursively cloning every component
        of this object. This is used for regression testing and should not
        be used in general (unnecessarily slow) */
    virtual XObjectSP recursiveClone() const;

    /** Returns the type of the external object */
    string getXClassName() const;

    /** Returns the type of the external object */
    XClassConstSP getXClass() const;

    /** CombinableResult interface: scale by factor x */
    virtual void scale(double x);

    /** CombinableResult interface:
        add an object (scaled by scaleFactor) to this
        result. Implementations should modify this result. If the x is
        not the same type as this then a [class cast] exception will
        be thrown */
    virtual void add(const CombinableResult& x, double scaleFactor);

    /** Turns DR Object into DR Map */
    XMapSP toMap() const;

    /** Turns object into map and retrieves relevant field */
    IDRObjectSP getField(const string& fieldName) const;
    /** Turns object into map and retrieves relevant field */
    IDRObjectSP getField(const char* fieldName) const;

    /** Invoke 'execute' on this object */
    IDRObjectSP execute() const;

   /** Turn DRValue into IDRObject. If takeOwnership is true then you
       must not free the DRValue - this is done if needed by this
       function. Alternatively if takeOwnership is false then it is up
       to the caller to free the DRValue (if necessary) */
    static IDRObjectSP toIDRObject(const DRValue& value,
                                   bool           takeOwnership,
                                   DRService*     svc);

    /** Same as above (ie toIDRObject) but for objects in DRValue */
    static XObjectSP toXObject(DRObject       drObj,
                               bool           takeOwnership);

    /** Turn IObject into DRValue - only free the DRValue for strings (in
        which case the return value is true (this is a right mess) */
    static bool fromIDRObject(
        const IDRObjectSP& obj, DRService* svc, DRValue& value);

    /** Try to convert the supplied IObject to an IDRObject. If the
        supplied object is an IDRObject already then this function
        just returns the casted object. If the object is an array of
        IDRObjects (eg array of Strings) then an XArray will be created of
        the type of the array. The service used to create the XArray is
        that in this object (ie only getService() property of this is used) */
    IDRObjectSP toIDRObject(IObjectSP obj) const;
    //// const version
    IDRObjectConstSP toIDRObject(IObjectConstSP obj) const;

    /** returns the enclosed DRObject which then must be freed. */
    DRObject getDRObject();

    //// Returns the DRService from the object
    DRService* getService() const;

    /** Return the array of types for this service */
    static XArraySP typeList(DRService* svc);

protected:
    /** Simple constructor. Takes ownership of memory */
    XObject(DRObject obj);
    static XObjectSP create(DRObject obj); /* as above but you can call this
                                              in derived class methods rather
                                              than constructors */
    //// sets pointer to null
    XObject(const CClassConstSP& clazz);

    XObject(const CClassConstSP&                    clazz,
            const MyDRObjectInterfaceSP& obj);

    XObject(const CClassConstSP& clazz,
            DRObject              obj);

    //// Derived classes should override this. Essentially X version of
    //// write()
    virtual void xWrite(const string& tag, Writer* writer) const;
    //// Derived classes should override this. Essentially X version of
    //// import()
    virtual void xImport(Reader::Node* elem, Reader* reader, DRService* svc);

    static string makeString(const char* str, DRService* svc);
    string makeString(const char* str) const;
    static void throwError(const char* method, const char *errorMsg);
    static void driErrorHandler(const char *method, const char *errMsg,
                                void *cbParam);

    //// Creates a DRString from the supplied stl string
    static DRString makeDRString(const string& str, DRService* svc);
    /** Is the DR Object a DRArray */
    static bool isArray(DRObject drObj);
    // fields
    MyDRObjectInterfaceSP object; // $unregistered
    bool owns() const { return weOwnObject;} // FIXME: eliminate using get_deleter() (see boost doc)
private:
    bool    weOwnObject; // true iff "object" is created with takeOwnership == true;

    class New;
    class Create;
    class XWriter;
    class XReader;
    class Wrapper;
    friend class Wrapper;
    /** Simple constructor. */
    explicit XObject(const MyDRObjectInterfaceSP& obj);

    XObject(const XObject& rhs); // don't use
    XObject& operator=(const XObject& rhs); // don't use
    static void load(CClassSP& clazz);
    XObject();
    static IObject* defaultConstructor();
    void outputWriteMatrix(const string& linePrefix,
                           const string& prefix,
                           ostream&      stream) const;

    /** Is the DR Object a matrix */
    static bool isMatrix(DRObject drObj, DRService* svc);
};

DRLIB_END_NAMESPACE
#endif
