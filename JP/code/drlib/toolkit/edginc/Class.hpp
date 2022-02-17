

#ifndef EDR_CLASS_HPP
#define EDR_CLASS_HPP
#include <cstddef>
#include "edginc/Object.hpp"
#include "edginc/Field.hpp"
#include "edginc/Collector.hpp"
using namespace std;
DRLIB_BEGIN_NAMESPACE

class IArray;
class CClassImp;

/** signature of method that classes must provide in order to load
    their class */
typedef void (TRegisterClassLoadMethod)(CClassSP& classToLoad);

typedef vector<CClassConstSP> CClassVec;
#ifndef QLIB_CLASS_CPP
//// commented out for now since it doesn't seem to really save anything
//// because the stl::vector methods are virtually all inline. It's been
//// left here because the syntax was a bit hard to get right
//EXTERN_TEMPLATE(class std::vector<const CClass* _COMMA_ std::allocator<const CClass *> >);
#endif

/** signature of method that classes must provide in order to create
    an empty instance of it */
typedef IObject* (TCreateInstanceMethod)();

/** signature of method that array classes must provide in order to create
    an instance of themselves of given length */
typedef IArray* (TCreateArrayInstanceMethod)(int length);

/** signature of static collect method */
typedef void (TAcceptMethod)(const void* instance, void* collector);
/** signature of object to interface casting method */
typedef void* (TIfaceCastMethod)(void* instance);

/** Instances of the class Class represent classes and interfaces in the
    library which implement the IObject interface.

    Typically instances of this class are created at application start up.
    Often classes refer to other classes and to ensure that there are no
    problems with the order which instances of CClass are created, there is
    a second step performed when the library starts up which is to invoke
    an optional method on all classes. This optional method is passed as a
    parameter to the CClass constructor. Classes should use this method
    to record information about themselves (such as what fields they have and
    what interfaces they implement). This information is used by the 
    infrastructure to allow default implementations of the IObject interface.
    */
class TOOLKIT_DLL CClass: public CObject{
    friend class CClassImp;
    friend class DRIArgumentDescriptor; // in Class.cpp
public:
    virtual ~CClass();
        
    /** will throw an exception - not supported */
    virtual IObject* clone() const;
        
    static CClassConstSP const TYPE;
    /** Supplied method will be invoked when the class type information
        is required. This will be typically be at start up. This solves
        the issue of the essentially random order the classes are
        initialised. Method can be null */
    static CClassConstSP registerClassLoadMethod(
        const char*                 className, // class name
        const type_info&            typeInfo,  // typeid of class
        TRegisterClassLoadMethod*   loadMethod);
        
    /** Same as above but the 'atomic' classes where we want to register two
        type_info classes against the type eg 'double' and CDouble */
    static CClassConstSP registerClassWithAliasLoadMethod(
        const char*                 className,
        const type_info&            typeInfo,  // typeid of class
        const type_info&            typeInfoAlias,  // typeid of alias class
        TRegisterClassLoadMethod*   loadMethod);

    /** same as registerClassLoadMethod but for interfaces */
    static CClassConstSP registerInterfaceLoadMethod(
        const char*                 className,
        const type_info&            typeInfo,  // typeid of class
        TRegisterClassLoadMethod*   loadMethod);
        
    /** A dummy version of registerClassLoadMethod for use by templates
        where we want to support multiple instances of the TYPE field
        across dlls.  (With MSVC it seems to be impossible to
        automatically export the TYPE field within a template). All this
        method does is to return the CClassConstSP object for the
        specified type_info reference. The idea is that the default
        template code initialises its TYPE field with a call to this
        function. Note that a call to registerClassLoadMethod must be made
        somewhere else (and once only) in order to provide the name and
        loadMethod. The reason for not putting the call to
        registerClassLoadMethod in the default template code is to a) avoid
        code bloat since the load method causes a lot of code to be created, and
        b) to avoid issues with specifying the name of the type */
    static CClassConstSP templateRegisterClass(
        const type_info&            classType);  // typeid of class

    /** register method which creates an empty instance of this type */
    void registerCreateMethod(TCreateInstanceMethod*  createMethod);

    /** load up all known classes ie invokes 'loadMethod' supplied by
        each class in call to registerClassLoadMethod or 
        registerInterfaceLoadMethod. This is invoked by the 
        Library::startUp() */
    static void loadAllClasses();

    /** Returns true if we are busy trying to invoke the load method for each
        type */
    static bool classLoadingInProgress();

    /** flag that a class represents a native type */
    void setNative();

    /** flag that a class represents an array type */
    void setArrayType(const type_info&            componentType,
                      TCreateArrayInstanceMethod* createArrayMethod);

    /** Flags that the supplied type should be used for DRI reflection
        information. In particular, for proxies (see CObject) this type is
        used for field reflection information whilst for IPublicObjects this
        type is used for parent information */
    void setDRIProxyType(CClassConstSP proxyType);

    /** Flag that this type represents an enum. The default is false */
    void setIsEnum();

    /** Record a possible value for a Class that represents an enum */
    void setEnumValue(int         enumValue,
                      const char* enumName,
                      const char* enumComment);

    /** Returns the Class object associated with the class or
        interface with the given string name. */
    static CClassConstSP forName(const string& className);

    /** Returns the Class object associated with the class or
        interface with the given type_inf0 */
    static CClassConstSP forName(const type_info& classType);

    /** Return an array of all known classes (including interfaces and
        private, protected etc) */   
    static const CClassVec& allClasses();

    /** Determines if the class or interface represented by this Class
        object is either the same as, or is a superclass or superinterface
        of, the class or interface represented by the specified Class
        parameter. */
    bool isAssignableFrom(CClassConstSP clazz) const throw() ;

    /** Determines if the specified Object is assignment-compatible
        with the object represented by this Class. */
    bool isInstance(const IObjectConstSP& obj) const throw();

    /** Determines if the specified Object is
        assignment-compatible with the object represented by this
        Class. */
    bool isInstance(const IObject* obj) const throw();

    /** Determines if the specified Object is
        assignment-compatible with the object represented by this
        Class. */
    bool isInstance(const IObject& obj) const throw();

    /** Return the C++ RTTI descriptor of the class */
    const type_info& getTypeInfo() const;

    /** Return the name of the class */
    const string& getName() const throw();

    /** Returns an array containing Class objects representing all
        the public classes and interfaces that are members of the
        class TOOLKIT_DLL represented by this Class object. */
    const CClassVec* getClasses() const;

    /** Determines the interfaces implemented by the class or
        interface represented by this object. */
    CClassVec getInterfaces() const;

    /** Returns the Class representing the superclass of the
        entity. If this Class represents either the Object class,
        an interface, a primitive type, or void, then null is
        returned. If this object represents an array class then
        the Class object representing the Object class is
        returned. */
    CClassConstSP getSuperClass() const;

    /** Determines if the specified Class object represents an
        interface type. */
    bool isInterface() const;

    /** Determines if the specified Class object represents a
        primitive type. Here primitive means that the C++ type actually used
        is not derived from IObject. (eg doubles, enums, strings etc) */
    bool isPrimitive() const;

    /** Returns true if this type represented an enum */
    bool isEnum() const;

    /** Allow infrastructure to deduce value of doesAssignmentOperatorClone()
        by examining reflection information. This must be called before
        any fields are registered for this class. In particular, if 
        doesAssignmentOperatorClone() is true for the parent class and all
        fields are 'inline' and doesAssignmentOperatorClone() is true for the
        type of each field, then doesAssignmentOperatorClone() will be true for
        this class. You should not call this function if you don't 
        understand what it does or your class has fields in it which are not
        registered with the reflection mechanism. */
    void enableCloneOptimisations();

    /* does the copy constructor of this class implicity perform a
       clone Eg no pointers anywhere including components. Examples include
       DateTime or DateTimeArray etc. This information is useful when cloning
       an instance of a class with a DateTimeArray (say) in it. If the 
       DateTimeArray is 'inline' then there is no need to clone the array
       before it is assigned */    
    bool doesAssignmentOperatorClone() const;

    /** Returns the modifiers for this class or interface, encoded in
        an integer. The Modifier class should be used to decode the
        modifiers. */
   int getModifiers() const;

    /** Returns a Field object that reflects the specified
        declared field of the class represented by
        this Class object */
    CFieldConstSP getDeclaredField(const string& name) const;

    /** as per getDeclaredField but returns 0 if no field exists rather
        than throwing an exception */
    CFieldConstSP hasDeclaredField(const string& name) const;

    /** Returns an array of Field objects reflecting all the
        fields declared by the class represented by
        this Class object (excludes inherited fields) */
    const CFieldArray& getDeclaredFields() const;
                        
    //// Utility class to hold information about an enum value. Class is
    //// immutable.
    class EnumValue{
    public:
        int    valueAsInt; //// the enum as an int
        string valueAsString; //// the enum as a string
        string comment; //// description for this value
        //// simple constructor - Only CClass should need this
        EnumValue(int           valueAsInt, 
                  const string& valueAsString, 
                  const string& comment);
    };

    /** Returns the first EnumValue whose valueAsString matches the given
        string. Only valid if this class represents an Enum. */
    const EnumValue& getEnumValue(const string& enumAsString) const;

    /** Returns the EnumValue whose valueAsInt matches the given enum value
        (as an int). Only valid if this class represents an Enum. */
    const EnumValue& getEnumValue(int enumAsInt) const;

    /** Returns the list of possible EnumValues for the type that this CClass
        represents. Returns null, if this CClass does not represent an enum.
        Do not delete the returned vector. */
    const vector<EnumValue>* getEnumValues() const;

    /** Sets a Field object that reflects the specified
        declared field on this class */
    void setDeclaredField(const CFieldSP& field);

    /** creates a field using the data given and stores on the 
        type given by typeName */
    void setDeclaredField(const char*              fieldName,   
                          const type_info&         typeId,  // of field
                          TFieldGetMethod*         getMethod,
                          TFieldSetMethod*         setMethod,
                          CField::PointerAttribute pointerAttribute,
                          ptrdiff_t                      offset); /* from start of 
                                                               structure */

    /** Specialised version of setDeclaredField for fields accessed by 
        plain pointers  */
    void setDeclaredPlainPointerField(const char*      fieldName,   
                                      const type_info& typeId,  // of field
                                      TFieldGetMethod* getMethod,
                                      TFieldSetMethod* setMethod,
                                      ptrdiff_t              offset); /*from start of 
                                                                  structure */
    /** Specialised version of setDeclaredField for fields accessed by 
        smart pointers  */
    void setDeclaredSmartPointerField(const char*      fieldName,   
                                      const type_info& typeId,  // of field
                                      TFieldGetMethod* getMethod,
                                      TFieldSetMethod* setMethod,
                                      ptrdiff_t              offset); /*from start of 
                                                                  structure */
    /** Specialised version of setDeclaredField for inLine fields */
    void setDeclaredInLineField(const char*      fieldName,   
                                const type_info& typeId,  // of field
                                TFieldGetMethod* getMethod,
                                TFieldSetMethod* setMethod,
                                ptrdiff_t              offset); /* from start of 
                                                             structure */

    /** Flags that this class is private. Private classes are, for 
        example, hidden from Analytics Service and the spreadsheet etc.
        Note that by default classes are private unless they implement
        the IPublicObject interface (in which case they are public) or they
        implement the IPrivateObject (in which case they are protected) */
    void setPrivate(); 

    /** Flags that this class is public. Public classes are, for 
        example, visible to Analytics Service and the spreadsheet etc */
    void setPublic();

    /** Flags that this class is exported. Exported classes are public
        and clients may legitimately assume order of arguments (e.g.
        spreadsheet addins) */
    void setExport();

    /** for documentation */
    void          setDescription(const string& desc);
    const string& getDescription() const;

    /** sets the isOptional flag for the field with name fieldName */
    void setFieldOptional(const char*   fieldName,
                          bool          isOptional);

    /** sets the transient flag for the field with name fieldName */
    void setFieldTransient(const char* fieldName,
                           bool        isTransient);

    /** sets the transientForIteration flag for the field
        with name fieldName */
    void setFieldTransientForIteration(const char*   fieldName,
                                       bool          isTransientForIteration);
    
    /** Create an array object of this  type and length. Note this is the
        type of the array, and not of the components */
    IArray* newArrayInstance(int length) const;

    /** Create an array object of this component type and length. Note
        this is the type of the components and not the type of the array
        itself */
    IArray* newArrayInstanceByComponent(int length) const;

    /** Returns the Class representing the component type of an array. */
    CClassConstSP getComponentType() const;

    /** records the fact that this class implements given interface */
    void setImplementedInterface(const type_info&  ifaceType,
                                 TIfaceCastMethod* toIface);

    /** Determines if this Class object represents an array class. */
    bool isArray() const;

    /** Determines if this Class object represents a matrix class. */
    bool isMatrix() const;

    /** set the superclass for types that have explicit implementation of
        castToBase() - typically only CObject */
    void setSuperClass(const type_info&  superClassTypeId);

    /** set the superclass for types that are derived non virtually from
        a class that explicitly implements castToBase(). The offset is the
        offset of the child from this base class */
    void setSuperClass(const type_info&  superClassTypeId,
                       ptrdiff_t               offset);

    /** Creates a new instance of the class represented by this
        Class object. */
    IObject* newInstance() const;

    /** Equivalent to object->getClass()->dynamicCast(object) except
        will return null if object null */
    static void* castFromIObject(const IObject* object);

    /* Given an IObject*, do a dynamic_cast to object of the type 
       represented by this class. Only works for types that aren't
       interfaces. Throws an exception if object of wrong type.
       If object is null, null will be returned */
    const void* dynamicCast(const IObject* object) const;

    /* Given an IObject*, do a dynamic_cast to object of the type 
       represented by this class. Only works for types that aren't
       interfaces. Throws an exception if object of wrong type.
       If object is null, null will be returned */
    void* dynamicCast(IObject* object) const;

    /* Given an IObject*, do a static_cast to object of the type 
       represented by this class. If object is null, null will be returned */
    const void* staticCast(const IObject* object) const;

    /* Given an IObject*, do a static_cast to object of the type 
       represented by this class. Only works for types that aren't
       interfaces. If object is null, null will be returned */
    void* staticCast(IObject* object) const;

    /** sets the description for the given field of objects of this type */
    void setFieldDescription(const char* fieldName,
                             const char* description);

    /** As above but takes string */
    void setFieldDescription(const char*   fieldName,
                             const string& description);

    /** sets the description for the given field of objects of this type */
    const string& getFieldDescription(const string& fieldName);

    /** Adds the supplied class to the list of types that can be used
        to build an instance of this class using the data dictionary
        route (ie supplied type implements the IPublicObject interface
        and the toPrivateObject() method returns an instance of this
        class.  It is not necessary to call this method for this class
        - it is automatically added to this list if appropriate */
    void addConstructorClass(CClassConstSP clazz);

    /** Return an array of types that can be used to construct this class
        (derived classes are excluded) via the data dictionary route. */
    const CClassVec& getConstructorClasses() const;

    /** Returns an array of types which can be used to construct this
        class, or one derived from it, via the data dictionary route. */
    CClassVec listAllConstructorClasses() const;

    /** Returns an array of all clazzes which satisfy
        this.isAssignableFrom(clazz).  In other words, if clazz is an
        interface a list of all classes that implement this interface is
        returned wother with all interfaces that extend this
        interface. Oterhwise a list of all classes which are derived from
        this one (which includes supplied clazz) is returned. */
    const CClassVec& getAssignableClasses() const;

    /** sets accept method for specific type of collector. 
        Use ClassSetAcceptMethod in preference to this */
    static void setAcceptMethod(
        const type_info& clazzID,
        const type_info& collectorID,
        TAcceptMethod*   acceptMethod);

    /** invokes accept method for specific type of collector. Returns
        true if an accept method was found. In general, use Accept on
        CObject method in preference to this */
    static bool invokeAcceptMethod(
        const IObject*    instance,
        ICollector*       collector);

    /** Return the DRIType for this CClass */
    IObjectSP getDRIType() const;
    
    /** Utility class for use in hashmap templates */
    struct Hash{
        size_t operator()(const CClassConstSP p) const {return reinterpret_cast<size_t> (p);} // VC6.0 misses uintptr_t which is more suitable than size_t here
    };
        
private:
    // methods
    /** create a new instance of a class - if not already initialised */
    static CClassSP create(const string& className,
                           const type_info& typeInfo);
    CClass(const string& className, const type_info* typeInfo);
    CClass(const CClassConstSP& rhs);
    CClass& operator=(const CClassConstSP& rhs);
    static void load(CClassSP& classToLoad);
    void setModifiers();
    /** load up class of given type */
    static void loadClass(const CClassConstSP& classToLoad);
    void* ifaceCast(IObject* object) const;
    //////////////// fields
    string           className; // the name
    const type_info* typeInfo;  // C++ RTTI descriptor
    bool             primitive; // is type double, int etc
    bool             loaded;    // has class load method been invoked
    bool             array;     // true - class represents an array
    bool             assignmentOperatorClones; /* does the copy constructor of
                                                  this class implicity perform a
                                                  clone Eg no pointers anywhere
                                                  inc components */
    CClassConstSP    superClass;
    int              modifiers; // identifies type of class
    ptrdiff_t        baseOffset;// offset from non virtual base class
    CClassImp*       data;      // lots of extra bits of data
    bool             isIface;   // true: class represents an interface
    string           description;    // what it's for
};

/** template function for registering accept method */
template<class Coll, class X> void ClassSetAcceptMethod(
    void (*method)(X*, Coll*))
{
    CClass::setAcceptMethod(typeid(X), typeid(Coll), (TAcceptMethod*) method);
}

/** class template containing static function template which will cast an
    object to a derived object (we're interested in the case where the
    derived object in an interface which is derived from virtually) */
template <class Iface> class CClassCast{
public:
    /** cast from object to interface */
    template<class Obj> static Iface* toIface(Obj* obj){
        return obj;
    }
};

/** class template containing static function template which tests that
    a class is derived non-virtually from another */
template <class Parent> class CClassUpCast{
public:
    /** cast from object to interface */
    template<class Obj> static Parent* toParent(Obj* obj){
        Parent* parent = obj; // checks Obj is derived from parent
        Obj* obj2 = static_cast<Obj*>(parent); // stop virtual inheritance
        return obj2; // stops compiler warnings
    }
};

/** For use with IMPLEMENTS macro. Stores relevent information on CClass
    object. The parameters are only used to retrieve the type information */
template <class Obj, class Iface> TIfaceCastMethod* CClassGetIfaceMethod(
    Obj* , Iface* ){
    typedef Iface* (TCastMethod)(Obj*);
    TCastMethod* method = &CClassCast<Iface>::toIface;
    return (TIfaceCastMethod*)method;
}


/** macro (for ease) for registering objects derived (directly or
    indirectly) non-virtually from CObject. Do not use this macro for
    interfaces. This will not compile if the inheritance is virtual.
    To do: change to pass pointer everywhere rather than reference to
    null pointer */
#define REGISTER(className, classVal) \
    className* instancePtr__ = 0; \
    className& instance__ = *instancePtr__; \
    instancePtr__ = &instance__; /* avoid compiler warning */ \
    /* any value for the nonNullInstancePtr__ will do, except zero */ \
    className* nonNullInstancePtr__ = (className*) 0x8000; \
    CClassSP   class__ = classVal;

/** If you get a compilation error here (or even a crash here at a
    memory address around 0x8000) then you are using the wrong macros
    or you are doing strange inheritance. You should use SUPERCLASS
    for classes that are derived (directly or indirectly) from
    CObject. You should derive non-virtually from your parent. Do not
    use this macro for interfaces. */
#define SUPERCLASS(parentName) {\
    parentName* parent__ = nonNullInstancePtr__;\
    ptrdiff_t offset = (char*)nonNullInstancePtr__ - (char*)parent__;\
    class__->setSuperClass(typeid(parentName), offset);\
    /* Next line checks we derive non-virtually from superclass */ \
    CClassUpCast<parentName >::toParent(instancePtr__); \
}

#define IMPLEMENTS(ifaceName) \
    class__->setImplementedInterface(typeid(ifaceName),\
                                     CClassGetIfaceMethod(instancePtr__, \
                                                          (ifaceName*)0))

#define EMPTY_SHELL_METHOD(method) \
    class__->registerCreateMethod(method)

#define FIELD_NO_DESC(name) \
    FieldAdd(class__, #name, instance__.name, instance__)

#define FIELD_USING_ALIAS_NO_DESC(name, alias) \
    FieldAdd(class__, #alias, instance__.name, instance__)

#define FIELD_MAKE_OPTIONAL(name) \
    class__->setFieldOptional(#name, true);

#define FIELD_MAKE_TRANSIENT(name) \
    class__->setFieldTransient(#name, true);

#define FIELD_MAKE_TWEAKABLE(name) \
    class__->setFieldTransientForIteration(#name, false);

#define FIELD_MAKE_NONTWEAKABLE(name) \
    class__->setFieldTransientForIteration(#name, true);

#define FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(name) \
    FIELD_MAKE_TRANSIENT(name) \
    FIELD_MAKE_TWEAKABLE(name) 

#define FIELD_SET_DESC(name, description) \
    class__->setFieldDescription(#name, description);

#define FIELD(name, description)\
    FIELD_NO_DESC(name); \
    FIELD_SET_DESC(name, description);

#define FIELD_USING_ALIAS(name, alias, description)\
    FIELD_USING_ALIAS_NO_DESC(name, alias); \
    FIELD_SET_DESC(alias, description);

/** For interfaces */
#define REGISTER_INTERFACE(className, classVal) \
    className* ifacePtr__ = 0; \
    CClassSP&  class__ = classVal;

#define EXTENDS(extendedIface) \
    class__->setImplementedInterface(typeid(extendedIface),\
                                     CClassGetIfaceMethod(ifacePtr__, \
                                                          (extendedIface*)0))

/** static cast macro - will convert IObject* to any object derived from 
    CObject. It is much faster than dynamic_cast<>. Note it does not
    validate the types are compatible with each other */
#define STATIC_CAST(__className, __ptr)\
    ((__className*)__className::TYPE->staticCast(__ptr))

/** as above but for const pointers */
#define STATIC_CONST_CAST(__className, __ptr)\
    ((const __className*)__className::TYPE->staticCast(__ptr))

/** dynamic cast macro - will convert IObject* to any object derived from 
    CObject. It is much faster than dynamic_cast<>.  */
#define DYNAMIC_CAST(__className, __ptr)\
    ((__className*)__className::TYPE->dynamicCast(__ptr))

/** as above but for const pointers */
#define DYNAMIC_CONST_CAST(__className, __ptr)\
    ((const __className*)__className::TYPE->dynamicCast(__ptr))


/** template function for inline field registration */
template <class T, class EnclosingClass> void FieldAdd(
    CClassSP         clazz,
    const char*      name,
    const T&         field,
    EnclosingClass&  enclosingClass)
{
    IObjectSP (*getMethod)(T*) = 
        FieldInLineMethods<T, boost::is_enum<T>::value>::get;
    void (*setMethod)(T*, IObjectSP) = 
        FieldInLineMethods<T, boost::is_enum<T>::value>::set;
    clazz->setDeclaredInLineField(
        name,
        typeid(T),
        (TFieldGetMethod*)getMethod,
        (TFieldSetMethod*)setMethod,
        reinterpret_cast<const char *>(&field) -
        reinterpret_cast<const char *>(&enclosingClass));
}        

/* template function for plain pointer field registration */
template <class T, class EnclosingClass> void FieldAdd(
    CClassSP         clazz,
    const char*      name,
    T*&              field,
    EnclosingClass&  enclosingClass)
{
    IObjectSP (*getMethod)(T**) = FieldGetPlainPtr;
    void (*setMethod)(T**, IObjectSP) = FieldSetPlainPtr;
    clazz->setDeclaredPlainPointerField(name,
                                        typeid(T),
                                        (TFieldGetMethod*)getMethod,
                                        (TFieldSetMethod*)setMethod,
                                        reinterpret_cast<const char *>(&field) - reinterpret_cast<const char *>(&enclosingClass));
}        

/* template function for smart pointer field registration */
template <class T, class EnclosingClass> void FieldAdd(
    CClassSP           clazz,
    const char*        name,
    const smartPtr<T>& field,
    EnclosingClass&    enclosingClass)
{
    IObjectSP (*getMethod)(smartPtr<T>*) = FieldGetSmartPtr;
    void (*setMethod)(smartPtr<T>*, IObjectSP) = FieldSetSmartPtr;
    clazz->setDeclaredSmartPointerField(name,
                                        typeid(T),
                                        (TFieldGetMethod*)getMethod,
                                        (TFieldSetMethod*)setMethod,
                                        reinterpret_cast<const char *>(&field) - reinterpret_cast<const char *>(&enclosingClass));
}        

/* template function for smart const pointer field registration */
template <class T, class EnclosingClass> void FieldAdd(
    CClassSP                clazz,
    const char*             name,
    const smartConstPtr<T>& field,
    EnclosingClass&         enclosingClass)
{
    IObjectSP (*getMethod)(smartPtr<T>*) = FieldGetSmartPtr;
    void (*setMethod)(smartPtr<T>*, IObjectSP) = FieldSetSmartPtr;
    clazz->setDeclaredSmartPointerField(name,
                                        typeid(T),
                                        (TFieldGetMethod*)getMethod,
                                        (TFieldSetMethod*)setMethod,
                                        reinterpret_cast<const char *>(&field) - reinterpret_cast<const char *>(&enclosingClass));
}        

/** Useful wrapper around clone which gives a pointer of the right
    type as well as giving good error messges. Fails if supplied
    parameter is NULL. The template ought to be static member of
    CObject but it fails to compile on gcc 2.8 if it is. */
template<class T> T* copy(const T* orig){
    CClassConstSP clazz = T::TYPE;
    if (!orig){
        throw IObject::makeExceptionForCopy(clazz);
    }
    IObject* theCopy = orig->clone();
    T* myCopy = ((T*) clazz->staticCast(theCopy));
    return myCopy;
}


/** Clones supplied object if the reference count is zero. Otherwise
    returns original object. Useful when you need to ensure a pointer
    won't go out of scope.  Returns NULL if supplied parameter is NULL.
    The template ought to be static member of CObject but it fails to
    compile on gcc 2.8 if it is. */
template<class T> const T* copyIfRef(const T* orig){
    CClassConstSP clazz = T::TYPE;
    if (!orig){
        return 0;
    }
    if (orig->getRefCount() == 0){
        IObject* theCopy = orig->clone();
        T* myCopy = ((T*) clazz->staticCast(theCopy));
        return myCopy;
    }
    return orig;
}

/** Useful wrapper around clone which gives a pointer of the right
    type as well as giving good error messges. Fails if supplied
    parameter is NULL. Implemented in Class.hpp since gcc doesn't like
    templates using classes that aren't yet defined. */
template<class T> T* IObject::copy(const T* orig){
    CClassConstSP clazz = T::TYPE;
    if (!orig){
        throw IObject::makeExceptionForCopy(clazz);
    }
    IObject* theCopy = orig->clone();
    T* myCopy = ((T*) clazz->staticCast(theCopy));
    return myCopy;
}

DRLIB_END_NAMESPACE

//// and then include specialisation of FieldInLineMethods for enums
#include "edginc/BoxedEnum.hpp"

#endif

