//----------------------------------------------------------------------------
//
//   Group       : Equity Derivatives Research
//
//   Filename    : XClass.hpp
//
//   Description : Describes classes from other DRI libraries
//
//   Author      : Mark A Robson
//
//   Date        : 27 Feb 2004
//
//----------------------------------------------------------------------------


#ifndef EDR_XCLASS_HPP
#define EDR_XCLASS_HPP
#include "edginc/XObject.hpp"

DRLIB_BEGIN_NAMESPACE

class XField;
typedef const XField* XFieldConstSP;
typedef XField* XFieldSP;
typedef vector<XFieldConstSP> XFieldArray;
class XClass;
typedef const XClass* XClassConstSP;
typedef XClass* XClassSP;
typedef vector<XClassConstSP> XClassVec;
class XClassImp;

/** Describes classes from other DRI libraries */
class TOOLKIT_DLL XClass: public CObject{
    friend class XClassImp;
public:
    virtual ~XClass();
        
    /** will throw an exception - not supported */
    virtual IObject* clone() const;
        
    static CClassConstSP const TYPE;

    /** Load up all the classes from the specified DRI library */
    static void loadAllClasses(DRService* svc);

    /** Returns the Class object associated with the class or
        interface with the given string name. */
    static XClassConstSP forName(const string& svcName,
                                 const string& className);

    /** Return an array of all known classes (including interfaces and
        private, protected etc) */   
    static const XClassVec& allClasses(const string& svcName);

    /** Determines if the class or interface represented by this Class
        object is either the same as, or is a superclass or superinterface
        of, the class or interface represented by the specified Class
        parameter. */
    bool isAssignableFrom(XClassConstSP clazz) const throw() ;

    /** Determines if the specified Object is assignment-compatible
        with the object represented by this Class. */
    bool isInstance(IObjectConstSP obj) const throw();

    /** Determines if the specified Object is
        assignment-compatible with the object represented by this
        Class. */
    bool isInstance(const IObject* obj) const throw();

    /** Determines if the specified Object is
        assignment-compatible with the object represented by this
        Class. */
    bool isInstance(const IObject& obj) const throw();

    /** Return the name of the class */
    const string& getName() const throw();

    /** Determines if the specified Class object represents a
        primitive type. */
    bool isPrimitive() const;
    
    /** Returns a Field object that reflects the specified
        declared field of the class represented by
        this Class object */
    XFieldConstSP getDeclaredField(const string& name) const;

    /** as per getDeclaredField but returns 0 if no field exists rather
        than throwing an exception */
    XFieldConstSP hasDeclaredField(const string& name) const;

    /** Returns an array of Field objects reflecting all the
        fields declared by the class represented by
        this Class object (NB includes inherited fields) */
    const XFieldArray& getDeclaredFields() const;
                        
    /** for documentation */
    const string& getDescription() const;

    /** Returns the Class representing the component type of an array. */
    XClassConstSP getComponentType() const;

    /** Determines if this Class object represents an array class. */
    bool isArray() const;

    /** Returns an array of all clazzes which satisfy
        this.isAssignableFrom(clazz).  In other words, if clazz is an
        interface a list of all classes that implement this interface is
        returned together with all interfaces that extend this
        interface. Otherwise a list of all classes which are derived from
        this one (which includes supplied clazz) is returned. */
    const XClassVec& getAssignableClasses() const;

    /** Utility class for use in hashmap templates */
    struct Hash{
        size_t operator()(const XClassConstSP p) const {return (size_t)p;}
    };
        
private:
    // methods
    XClass(const string& svcName, const string& className);
    XClass(const XClassConstSP& rhs);
    XClass& operator=(const XClassConstSP& rhs);
    static void myLoad(CClassSP& classToLoad);
    static XClassSP create(const string& svcName, const string& className,
                           bool       isArray, bool isPrimitive,
                           bool       loaded);
    static XClassSP createPrimitiveType(const string& svcName, int type);
    static XClassSP createArrayType(const string& svcName, XClassSP eltClazz);
    static XClassSP createMatrixType(const string& svcName);
    void setDeclaredField(const string& svcName, XMapSP fieldDesc);
    static XClassSP load(const string& svcName, XObjectSP clazzDesc);
    static void populateAllParents(const string& svcName);
    static XClassSP forNameInternal(const string& svcName,
                                    const string &className);
    
    //////////////// fields
    string        className; // the name $unregistered
    bool          primitive; // is type double, long etc $unregistered
    bool          array;     // true - class represents an array $unregistered
    XClassImp*    data;      // lots of extra bits of data $unregistered
};
DRLIB_END_NAMESPACE

#endif

