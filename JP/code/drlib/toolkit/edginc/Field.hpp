
#ifndef EDR_FIELD_HPP
#define EDR_FIELD_HPP
#include "edginc/Collector.hpp"

DRLIB_BEGIN_NAMESPACE

typedef IObjectSP (TFieldGetMethod)(void* structStart);
typedef void (TFieldSetMethod)(void* structStart, IObjectSP obj);

/** A Field provides information about, and dynamic access to, a
    single field of a class or an interface. The reflected field may
    be a class (static) field or an instance field. */
class TOOLKIT_DLL CField: public CObject{
    friend class ShutTheCompilerUp;
public:
    /** With C++ there's a range of ways a field can be declared */
    enum PointerAttribute{
        INLINE = 0, // ie not a pointer at all
        PLAIN_POINTER,
        SMART_POINTER
    };

    /** clones a field */
    virtual IObject* clone() const;

    static CClassConstSP const TYPE;

    /** creates a field using the data given */
    static CField* create(CClassSP         declaringClass,
                          const string&    name,
                          const type_info& typeId,  // of field
                          TFieldGetMethod* getMethod,
                          TFieldSetMethod* setMethod,
                          PointerAttribute pointerAttribute,
                          int              offset); /* from start
                                                       of struct */

    /** All the get/set functions take IObjectSP parameters by value.
        This allows clients to pass instance of smart pointers and the
        compiler will convert */

    /** Returns the value of the field represented by this Field,
        on the specified object. An smart pointer is returned since
        if the type is native it needs to be wrapped */
    IObjectConstSP constGet(const IObjectConstSP& obj) const;
    /** non const version - awkward if the same name as clients passing
        a non const smart pointer would have to specify which call they
        wanted */
    IObjectSP get(const IObjectSP& obj) const;

    /** Gets the value of a field as a boolean on the specified object. */
    bool getBool(const IObjectConstSP& obj) const;

    /** Gets the value of a field as a double on the specified object. */
    double getDouble(const IObjectConstSP& obj) const;

    /** Gets the value of a field as an int on the specified object. */
    int getInt(const IObjectConstSP& obj) const;

    /** Gets the value of a field as a string on the specified object. */
    const string& getString(const IObjectConstSP& obj) const;

    /** Sets the value of the field represented by this Field,
        on the specified object. */
    void set(const IObjectSP& obj, const IObjectSP& value) const;

    /** Sets the value of a field as a boolean on the specified object. */
    void setBool(const IObjectSP& obj, bool z) const;

    /** Sets the value of a field as a double on the specified object. */
    void setDouble(const IObjectSP& obj, double d) const;

    /** Sets the value of a field as an int on the specified object. */
    void setInt(const IObjectSP& obj, int i) const;

    /** Sets the value of a field as a string on the specified object. */
    void setString(const IObjectSP& obj, const string& s) const;

    /** Returns the name of the field represented by this Field object. */
    const string& getName() const throw();

    /** Returns a Class object that identifies the declared type
        for the field represented by this Field object. */
    CClassConstSP getType() const throw();

    /** Returns the Class object representing the class or
        interface that declares the field represented by this
        Field object. */
    CClassConstSP getDeclaringClass() const throw();

    /** Returns the modifiers for the field represented
        by this Field object, as an integer. (see Modifier class) */
    int getModifiers() const;

    /** indicates whether the the type of data held by the field is
        atomic */
    bool typeIsPrimitive() const throw();

    /** indicates whether the the type of data held by the field is
        an array */
    bool typeIsArray() const throw();
        
    /** flag the optionality of a field */
    void setOptional(bool isOptional);

    /** flag that the field is transient. This also sets the
        isTransientForIteration to the same value */
    void setTransient(bool isTransient);

    /** flag that the field is transient for the purposes of ObjectIteration */
    void setTransientForIteration(bool isTransient);

    /** is the field transient for the purposes of ObjectIteration */
    bool isTransientForIteration() const;

    /* is the field optional */
    bool isOptional() const;

    /** flag the optionality of a field */
    void setDescription(const string& description);
        
    /* is the field optional */
    const string& getDescription() const;

    /** Returns the pointer attribute for this field ie indicates whether the
        field is a pointer, a smart pointer or 'inline' */
    PointerAttribute getPointerAttribute() const;

    /** Copies (ie clones) component idenitified by this field from
        source to destination. Method is essntially a get followed by
        a clone and thena set but exists for performance reasons (eg
        when handling native types like doubles) */
    void copyComponent(const IObjectConstSP& source, 
                       const IObjectSP& destination) const;

private:
    // Fields must be accessed via pointers
    CField();
    //CField(const CField &rhs);
    //CField& operator=(const CField& rhs);
    void* dataAddress(const IObject*   obj,
                      const string&    routine) const;
    
    IObjectSP getObj(const IObject*    obj) const;

    // **** fields *****
    // the name of the field
    string   name;
    // to which class the field belongs
    CClassConstSP declaringClass;
    // the type of the element the field contains
    CClassConstSP type;
    int      modifiers; // holds whether transient etc
    bool     typePrimitive; //derived from type
    bool     typeArray; //derived from type
    int      offset;  // where the field goes in the structure
    TFieldGetMethod *getMethod;
    TFieldSetMethod *setMethod;
    bool     optional;
    bool     transientForIteration; /* true: field appears to be transient for
                                       ObjectIteration */
    string   description;
    PointerAttribute pointerAttribute;
};

/** template functions for getting/setting fields within objects */
template <class T> IObjectSP FieldGetInLine(T*  t){
    // don't own memory so use attach to reference
    return IObjectSP::attachToRef(t);
}
template <class T> void FieldSetInLine(T* t, IObjectSP obj){
    T* newT = (T*) T::TYPE->staticCast(obj.get());
    *t = *newT; // structure copy
}
template <class T> IObjectSP FieldGetPlainPtr(T**  t){
    // don't own memory so use attach to reference
    return IObjectSP::attachToRef(*t);
}
template <class T> void FieldSetPlainPtr(T** t, IObjectSP obj){
    IObject* newValObject;
    // must get memory from valueSP
    try{
        newValObject = obj.release();
    } catch (exception& e){
        throw ModelException(e, "setPlainPtr", 
                             "Couldn't release "
                             "memory from smart pointer");
    }
    T* newT = (T*) T::TYPE->staticCast(newValObject);
    EDR_DELETE(*t);
    *t = newT;
}

template <class T> IObjectSP FieldGetSmartPtr(T* t){
    return IObjectSP(*t);
}

template <class T> void FieldSetSmartPtr(T* t, IObjectSP obj){
    *t = T::dynamicCast(obj);
}    

/* Specialisations for native types */
template <> TOOLKIT_DLL IObjectSP FieldGetInLine(string* str);
template <> TOOLKIT_DLL void FieldSetInLine(string* str, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetInLine<double>(double* db);
template <> TOOLKIT_DLL void FieldSetInLine<double>(double* db, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetInLine<int>(int* i);
template <> TOOLKIT_DLL void FieldSetInLine<int>(int* i, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetInLine<bool>(bool* b);
template <> TOOLKIT_DLL void FieldSetInLine<bool>(bool* b, IObjectSP obj);

template <> TOOLKIT_DLL IObjectSP FieldGetPlainPtr(double** db);
template <> TOOLKIT_DLL void FieldSetPlainPtr(double** db, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetPlainPtr(int** intVal);
template <> TOOLKIT_DLL void FieldSetPlainPtr(int** intVal, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetPlainPtr(bool** boolVal);
template <> TOOLKIT_DLL void FieldSetPlainPtr(bool** boolVal, IObjectSP obj);
template <> TOOLKIT_DLL IObjectSP FieldGetPlainPtr(string** str);
template <> TOOLKIT_DLL void FieldSetPlainPtr(string** str, IObjectSP obj);

template <> TOOLKIT_DLL IObjectSP FieldGetInLine(IObjectSP* t);
template <> TOOLKIT_DLL void FieldSetInLine(IObjectSP* t, IObjectSP obj);

/** Template class to handle 'in line' methods for fields which aren't pointers
    (or smart pointers). In particular, it's a templated class to handle
    enums. We can use what's called 'type traits' to distinguish between
    enums and other types at compile time. The second parameter acts as
    a switch. If false then this is for non-enums. 
    To do: replace all FieldGetInLine and FieldSetInLine with specialisations
    of this class */
template<class T, bool = false /* not an enum */> class FieldInLineMethods{
public:
    static IObjectSP get(T* t){
        return FieldGetInLine(t);
    }
    static void set(T* t, IObjectSP obj){
        FieldSetInLine(t, obj);
    }
};
DRLIB_END_NAMESPACE

#endif
