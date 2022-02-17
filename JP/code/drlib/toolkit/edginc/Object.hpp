//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Object.hpp
//
//
//
//----------------------------------------------------------------------------

#ifndef EDG_OBJECT_H
#define EDG_OBJECT_H

#include "edginc/config.hpp"

DRLIB_BEGIN_NAMESPACE
// typedef's here to avoid problems with include files - avoid including
// class header file
class IObject;
class CClass;
typedef const CClass* CClassConstSP;
typedef CClass* CClassSP;
class ICollector;
// similar story with CField class
class CField;
typedef const CField* CFieldConstSP;
typedef CField* CFieldSP;
typedef vector<CFieldConstSP> CFieldArray;
// work around for smartPtr wanting to use method on CClass - although a bit
// ugly seems less work than spending hours trying to sort out include of
// headers. Only arisen when trying to move to gcc 3.4
// This function is implemented in Class.cpp
TOOLKIT_DLL const void* classDynamicCast(CClassConstSP  clazz, 
                                         const IObject* object);

DRLIB_END_NAMESPACE
#include "edginc/ModelException.hpp"
#include "edginc/smartPtr.hpp"
DRLIB_BEGIN_NAMESPACE

typedef smartPtr<IObject> IObjectSP;
typedef smartConstPtr<IObject> IObjectConstSP;

#ifndef QLIB_OBJECT_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<IObject>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IObject>);
#else
// seem to need to instantiate template before returning to cpp file - unclear
// why
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartConstPtr<IObject>);
INSTANTIATE_TEMPLATE(class TOOLKIT_DLL smartPtr<IObject>);
#endif

DRLIB_END_NAMESPACE
#include "edginc/Reader.hpp"
DRLIB_BEGIN_NAMESPACE

class Writer;

/** IObject defines the interface that all our objects must implement in
    order to be used by the infrastructure. Most classes will want to inherit
    from CObject as this provides a default implementation of all relevant
    methods. Update: getClass() and getRefCount() are no longer pure or
    virtual. This is for performance reasons - also no-one has ever overridden
    the default implementation.
 */
class TOOLKIT_DLL IObject{
public:
    /** the class representing the IObject interface */
    static CClassConstSP const TYPE;

    /**  Returns the runtime class of an object. Note inline and not virtual
         for performance */
    CClassConstSP getClass() const;

    /** Here clone means a deep copy of the given object.  Some
        classes may throw an exception. The original plan was for an
        object to support the clone method if it explicitly implements
        the Clonable marker interface - but we haven't done this */
    virtual IObject* clone() const = 0;

    /** Called after an instance of a class (which implements this
        interface) is constructed via 'pop2Object' ie reflection is
        used to fill the internal fields of the object directly */
    virtual void validatePop2Object() = 0;

    /** Called after fields within this object have been changed using the
        reflection infrastructure (chiefly through object iteration and 
        tweaking). The list of fields specifies which fields have been
        updated */
    virtual void fieldsUpdated(const CFieldArray& fields) = 0;

    virtual ~IObject();
    
    /** method for processing collectors. The general idea is to provide
        a way of collecting information from the object which is 'visited' by
        the collector. The infrastructure supports a double-dispatch mechanism
        which means that a method for each flavour of collector can be 
        specified. True is returned if the accept method was delegated to
        an appropriate handler ie if the double dispatch mechanism found a
        method for the particular flavour of collector supplied. T*/
    virtual bool accept(ICollector* collector) const = 0;

    /** write object out */
    virtual void write(const string& tag, Writer* writer) const = 0;

    /** populate an empty object from a Reader */
    virtual void import(Reader::Node* elem, Reader* reader) = 0;

    /** write object out in 'output' format - ie suitable for comparing
        regression files with. Format is essentially lines of the form
        "linePrefix <fieldName>_<fieldName>...: value". 
        Underscores are used to separate the fields. An initial underscore is
        not written if the initial prefix is "" */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix,
                             ostream&      stream) const = 0;

    /** Returns a hash code value for the object. This method is supported for
        the benefit of hashtables.
        The general contract of hashCode is: 

        (i) Whenever it is invoked on the same object more than once during an
        execution of the library, the hashCode method must consistently return
        the same integer, provided no information used in equalTo comparisons
        on the object is modified. This integer need not remain consistent 
        from one execution of an application to another execution of the same 
        application. 
        (ii) If two objects are equal according to the equalTo(IObject) method,
        then calling the hashCode method on each of the two objects must produce
        the same integer result. 
        (iii) It is not required that if two objects are unequal according to 
        the equals(IObject) method, then calling the hashCode method on each 
        of the two objects must produce distinct integer results. However, the 
        programmer should be aware that producing distinct integer results for
        unequal objects may improve the performance of hashtables. */
    virtual int hashCode() const = 0;

    /** Indicates whether some other object is "equal to" this one. 
        The equals method implements an equivalence relation on non-null objects
        (i) It is reflexive: for any non-null x, x.equals(x) should return true.
        (ii) It is symmetric: for any non-null x and y, x.equals(y) should 
        return true if and only if y.equals(x) returns true. 
        (iii) It is transitive: for any non-null reference values x, y, and z, 
        if x.equals(y) returns true and y.equals(z) returns true, 
        then x.equals(z) should return true. 
        (iv) It is consistent: for any non-null values x and y, multiple 
        invocations of x.equals(y) consistently return true or consistently
        return false, provided no information used in equalTo comparisons on 
        the objects is modified. 
        (v) For any non-null x, x.equals(null) should return false. 

        Note that it is generally necessary to override the hashCode
        method whenever this method is overridden, so as to maintain
        the general contract for the hashCode method, which states
        that equal objects must have equal hash codes.  */
    virtual bool equalTo(const IObject* obj) const = 0;

    /** Return a reference to an int which acts as a reference count for
        smart pointers - see smartPtr.hpp for meaning of return parameter.
        Ideally would be protected but issues with smartPtr template.
        Note inline and not virtual for performance */
    int& getRefCount() const;

    /** Helper function for quick casting. Function returns 'this' where
        'this' is the first class which is not derived from virtually.
        ie it should cast be equivalent to dynamic_cast<CObject*> */
    virtual void* castToBase() const = 0;

    /** A very brief English name or description of the object,
        for use in constructing e.g. error messages */
    virtual string toString() const = 0;

    /** toString() with NULL check */
    static string stringOf(IObjectConstSP);

    /** Useful wrapper around clone which gives a pointer of the right
        type as well as giving good error messges. Fails if supplied
        parameter is NULL. Implemented in Class.hpp since gcc doesn't like
        templates using classes that aren't yet defined. */
    template<class T> static T* copy(const T* orig); // in Class.hpp

    /** Creates an exception reporting that a clone failed for the supplied
        class TOOLKIT_DLL because the input object was null */
    static ModelException makeExceptionForCopy(
        CClassConstSP    clazz);        // what type was being copied

protected:
    IObject();
    IObject(const IObject& obj);
    IObject& operator=(const IObject& rhs);

    /* not done via constructor because can derive multiply from this class */
    void setObjClass(CClassConstSP objClass);
private:
    // fields
    CClassConstSP objClass; // run time type information
    mutable int   refCount;
};
//// to reduce the size of the libraries (ie .lib or .a) in debug we make
//// these methods non-inline for debug build
#if !defined(DEBUG) || defined(QLIB_OBJECT_CPP)
OPT_INLINE CClassConstSP IObject::getClass() const{ return objClass;}
OPT_INLINE IObject::~IObject(){}  // inline for performance
OPT_INLINE int& IObject::getRefCount() const{ return refCount;}
OPT_INLINE IObject::IObject(): refCount(0) {}
OPT_INLINE IObject::IObject(const IObject& obj)
    :objClass(obj.objClass), refCount(0) {}
/** Leave refCount and objClass on lhs and rhs unchanged */
OPT_INLINE IObject& IObject::operator=(const IObject& rhs){
    return *this;}
/* not done via constructor because can derive multiply from this class */
OPT_INLINE void IObject::setObjClass(CClassConstSP objClass){
    this->objClass = objClass;}
#endif

/** Write a very brief English name or description of an object to a stream */
TOOLKIT_DLL ostream& operator <<(ostream&, const IObject&);

/** Default implemententation of IObject class. All classes which wish to 
    use the facilities provided by the infrastructure must derive from
    this class. These facilities are mainly xml read/write, cloning,
    and creation from data dictionaries */
class TOOLKIT_DLL CObject: public virtual IObject{
public:
    /** attribute used to identify an object's type */
    static const string OBJECT_TYPE_ATTRIBUTE;
    static CClassConstSP const TYPE;

    virtual IObject* clone() const;

    virtual ~CObject();

    /** default implementation - does nothing */
    virtual void validatePop2Object();

    /** default implementation - does nothing */
    virtual void fieldsUpdated(const CFieldArray& fields);

    /** default implementation - delegates to CClass::invokeAcceptMethod */
    virtual bool accept(ICollector* collector) const;

    /** write object out. If the object implements the
        IPrivateObject interface, then before written out the object is
        converted to an IPublicObject using the toPublicObject() method */
    virtual void write(const string& tag, Writer* writer) const;

    /** read in an object from a Reader. If the object read in
        implements the IPublicObject object, then it is converted to an 
        IPrivateObject using  the toPrivateObject() method  */
    static IObject* read(Reader::Node* elem, Reader* reader);

     /** populate an empty object from a Reader */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix,
                             ostream&      stream) const;

    /** Returns a hash code value for the object by combining the
        hashCode for each non-transient field making this object. Note
        that this method might be relatively slow due to the recursive
        nature of the implementation */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one by
        comparing each non-transient field making up the pair of
        objects (assuming they are of the same type). Note that this
        method might be relatively slow due to the recursive nature of
        the implementation */
    virtual bool equalTo(const IObject* obj) const;

    /** Utility method to convert an object to its private representation
        if it has one. The supplied smart pointer is modified in place */
    static IObjectConstSP convertToPrivateRep(const IObjectConstSP& pubObj);

    /** Utility method to convert an object to its public representation
        if it has one. The supplied smart pointer is modified in place */
    static IObjectConstSP convertToPublicRep(const IObjectConstSP& privObj);

    /** Utility method to convert an object to its private representation
        if it has one. The supplied smart pointer is modified in place */
    static IObjectSP convertToPrivateRep(const IObjectSP& pubObj);

    /** Utility method to convert an object to its public representation
        if it has one. The supplied smart pointer is modified in place */
    static IObjectSP convertToPublicRep(const IObjectSP& privObj);

    /** Note static method. Takes supplied object and checks to see if it
        is of the desired type. If not it will try and convert the object
        to the desired type by using the ITypeConvert interface provided
        the supplied object implements it or by converting to a private/public
        representation. If both the supplied object and
        the desired type is an array and the components of the supplied
        array implement the ITypeConvert interface then the array is
        converted element by element. Else if the object is an array, the
        table of TObjectFromArray methods will be consulted to see if
        there is a suitable method tto convert the array into an
        object. The return value reflects whether the object was converted
        or not (true = object was converted). If it is not possible to
        convert the object to the required type an exception is thrown */
    static bool checkType(IObjectSP&     object,
                          CClassConstSP  desiredType);

    /** Support for conversion of objects from arrays  */
    typedef IObjectSP (TObjectFromArray)(const IObjectSP&  object,
                                         CClassConstSP     requiredType);

    /** Register a conversion method from an array to a given type. The
        given class should be the type of the array from which you want to
        convert. */
    static void registerObjectFromArrayMethod(CClassConstSP     givenClass,
                                              CClassConstSP     targetClass,
                                              TObjectFromArray* method);

    /** Support "proxy" interfaces
        internally use objClass, publicly use proxyClass as if it was objClass
        i.e. publicly build class A, but with class B (the proxy) interface */

    /** how to go from object to proxy and vice-versa */
    typedef IObjectSP (ProxyConvertMethod)(const IObjectConstSP& obj);

    /** register a proxy for a class and how to transform in both directions */
    static void registerObjectProxy(CClassConstSP       objClass, 
                                    CClassConstSP       proxyClass,
                                    ProxyConvertMethod* toProxy,
                                    ProxyConvertMethod* fromProxy);

    /** given a class, see if it has a proxy, or what it is a proxy for */
    static CClassConstSP proxyClass(CClassConstSP objClass);
    static CClassConstSP actualClass(CClassConstSP proxyClass);

    /** if an object has a proxy interface, convert to the proxy */
    static IObjectSP toProxy(const IObjectSP& object);

    /** if an object is a proxy, convert to the "real" object */
    static IObjectSP fromProxy(const IObjectSP& object);

    /** returns this as a void* */
    virtual void* castToBase() const;

    /** A very brief English name or description of the object,
        for use in constructing e.g. error messages.  The default
        implementation returns '<type> @<address>' */
    virtual string toString() const;

    /** utility method. Throws an exception if refCount < 1. This is
        useful for ensuring that the object in question is being
        accessed via smartPointers since it means that an extra
        reference to it can be easily made (eg via attachToRef)
        so stopping the object going out of scope */
    void ensurePosRefCount() const;


protected:
    //// inline for performance
    CObject(const CClassConstSP& objClass);

    //// inline for performance
    CObject(const CObject& obj);

    /** Leave objClass and refCount on lhs and rhs unchanged */
    CObject& operator=(const CObject& rhs);
private:
    static void load(CClassSP& classToLoad);
};

//// to reduce the size of the libraries (ie .lib or .a) in debug we make
//// these methods non-inline for debug build
#if !defined(DEBUG) || defined(QLIB_OBJECT_CPP)
OPT_INLINE CObject::~CObject() {} 
OPT_INLINE CObject::CObject(const CClassConstSP& objClass) { setObjClass(objClass); } 
OPT_INLINE CObject::CObject(const CObject& obj) : IObject() { setObjClass(obj.getClass()); }
/** Leave objClass and refCount on lhs and rhs unchanged */
OPT_INLINE CObject& CObject::operator=(const CObject& rhs){ return *this; }
#endif

# if defined(_MSC_VER)
// disable warning about inheriting method from base class rather than virtual
// class
#pragma warning(disable : 4250)
#endif

DRLIB_END_NAMESPACE

#endif

