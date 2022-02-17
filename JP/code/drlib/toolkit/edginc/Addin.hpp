//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Addin.hpp
//
//   Description : Class for representing addins
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ADDIN_HPP
#define EDG_ADDIN_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/FunctorMember.hpp"

DRLIB_BEGIN_NAMESPACE

/** Class for representing addins. Each instance
    of this class corresponds to one excel addin function. The TOOLKIT_DLL
    macro is not used because this class defines templated functions. Instead
    (because of the static data members) we must export each function/field */
class Addin{
public:
    /** Categories for addin functions */
    static TOOLKIT_DLL const string UTILITIES;
    static TOOLKIT_DLL const string RISK;
    static TOOLKIT_DLL const string XL_TESTS;
    static TOOLKIT_DLL const string MARKET;
    static TOOLKIT_DLL const string FLEX_PAYOFF;
    static TOOLKIT_DLL const string CONV_BOND;

    /* legacy typedefs for different methods */
    typedef IObjectSP (ObjMethod)(void* object);
    typedef double   (DoubleMethod)(void* object);
    typedef bool     (BoolMethod)(void* object);
    typedef string   (StringMethod)(void* object);
    typedef int      (IntMethod)(void* object);

    /** Signature of [member] methods supported by addin infrastructure */
    union Method{
        Functor<IObjectSP>* objMethod;
        Functor<double>*    doubleMethod;
        Functor<bool>*      boolMethod;
        Functor<string>*    stringMethod;
        Functor<int>*       intMethod;
    };

    /** How addin functions return data to the spreadsheet. To do: move lower
        in the class hiearchy */
    typedef enum ReturnStyle{
        returnHandle = 0,     // Return a single handle to the sheet
        expandSimple,         /* Expand output in a single
                                 'column'. For arrays this will result
                                 in each element in the array
                                 being returned on a separate row.  */
        expandMulti           /* expand output across columns and rows. There
                                 will be a column for each field within 
                                 the object. */
    } ReturnStyle;
    
    /** Identifies to what type of C++ method the addin function
        corresponds to */
    typedef enum MethodType{
        constructor = 0,
        instanceMethod,   /* acts upon an instance of a type */
        classMethod,      /* belongs with a certain type */
    } MethodType;


    /** Register an addin which returns a double */
    template <class T> static void registerDoubleMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        double (T ::*      method)()){  // holds function pointer
        Method addinMethod;
        addinMethod.doubleMethod = new FunctorMember<T, double>(method);
        addToHashtable(
            new Addin(addinName, category, description, T::TYPE,
                      false, expandSimple, classMethod, 
                      addinMethod, addinMethod.doubleMethod,
                      true, CDouble::TYPE));
    }
    
    /** Register an addin which returns an int */
    template <class T> static void registerIntMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        int (T ::*         method)()){  // holds function pointer
        Method addinMethod;
        addinMethod.intMethod = new FunctorMember<T, int>(method);
        addToHashtable(
            new Addin(addinName, category, description, T::TYPE,
                      false, expandSimple, classMethod, 
                      addinMethod, addinMethod.intMethod,
                      true, CInt::TYPE));
    }
    
    /** Register an addin which returns a bool */
    template <class T> static void registerBoolMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        bool (T ::*        method)()){  // holds function pointer
        Method addinMethod;
        addinMethod.boolMethod = new FunctorMember<T, bool>(method);
        addToHashtable(
            new Addin(addinName, category, description, T::TYPE,
                      false, expandSimple, classMethod, 
                      addinMethod, addinMethod.boolMethod,
                      true, CBool::TYPE));
    }
    
    /** Register an addin which returns a string */
    template <class T> static void registerStringMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        string (T ::*      method)()){  // holds function pointer
        Method addinMethod;
        addinMethod.stringMethod = new FunctorMember<T, string>(method);
        addToHashtable(
            new Addin(addinName, category, description, T::TYPE,
                      false, expandSimple, classMethod, 
                      addinMethod, addinMethod.stringMethod,
                      true, CString::TYPE));
    }
    
    /** Register an addin which returns a IObjectSP (or an SP of some class
        derived from IObject).
        NOTE THAT NO COPY WILL BE MADE OF THE OBJECT RETURNED - THIS 
        MEANS THAT IF YOU DO NOT WANT THE DATA TO BE MODIFIED BY 
        SUBSEQUENT OTHER ADDIN FUNCTIONS YOU MUST RETURN A COPY */
    template <class T, class ReturnType> static void registerObjectMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        bool               handleNameParam, //true: extra input for handle name
        ReturnStyle        returnStyle, // how to return the output
        smartPtr<ReturnType> (T ::* method)()){     // holds function pointer
        Method addinMethod;
        addinMethod.objMethod = 
            new FunctorMember<T, IObjectSP, smartPtr<ReturnType> >(method);
        addToHashtable(
            new Addin(addinName, category, description, T::TYPE,
                      handleNameParam, returnStyle, classMethod, 
                      addinMethod, addinMethod.objMethod,
                      false, ReturnType::TYPE));
    }
    
    /** Record data defining an addin which, if written in C++,
        correponds to a constructor method which returns a handle */
    static TOOLKIT_DLL void registerConstructor(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass);  /* identifies class 
                                                representing addin data */

    /** Deprecated - use registerDoubleMethod.
        Record data defining an addin which, if written in C++,
        correponds to a static class method which returns a double */
    static TOOLKIT_DLL void registerClassDoubleMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        DoubleMethod*      method);     // holds function pointer
    
    /** Deprecated - use registerObjectMethod.
        Record data defining an addin which, if written in C++,
        correponds to a static class method which returns an IObjectSP.
        NOTE THAT NO COPY WILL BE MADE OF THE OBJECT RETURNED - THIS 
        MEANS THAT IF YOU DO NOT WANT THE DATA TO BE MODIFIED BY 
        SUBSEQUENT OTHER ADDIN FUNCTIONS YOU MUST RETURN A COPY */
    static TOOLKIT_DLL void registerClassObjectMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        bool               handleNameParam, //true: extra input for handle name
        ReturnStyle        returnStyle, // how to return the output
        ObjMethod*         method);     // holds function pointer
    
    /** Deprecated - use registerStringMethod.
        Record data defining an addin which, if written in C++,
        correponds to a static class method which returns a string. */
    static TOOLKIT_DLL void registerClassStringMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        StringMethod*      method);     // holds function pointer
    
    /** Deprecated - use registerIntMethod.
        Record data defining an addin which, if written in C++,
        correponds to a static class method which returns an int */
    static TOOLKIT_DLL void registerClassIntMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        IntMethod*      method);     // holds function pointer
    
    /** Deprecated - use registerBoolMethod.
        Record data defining an addin which, if written in C++,
        correponds to a static class method which returns a bool. */
    static TOOLKIT_DLL void registerClassBoolMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        BoolMethod*        method);     // holds function pointer
    
    /** Deprecated - use registerDoubleMethod.
        Record data defining an addin which, if written in C++,
        correponds to an instance method (ie not static) which returns
        a double */
    static TOOLKIT_DLL void registerInstanceDoubleMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        DoubleMethod*      method);     // holds function pointer
    
    
    /** Deprecated - use registerObjectMethod.
        Record data defining an addin which, if written in C++,
        correponds to an instance method (ie not static) which returns
        an IObjectSP.  NOTE THAT NO COPY WILL BE MADE OF THE OBJECT
        RETURNED - THIS MEANS THAT IF YOU DO NOT WANT THE DATA TO BE
        MODIFIED BY SUBSEQUENT OTHER ADDIN FUNCTIONS YOU MUST RETURN A
        COPY */
    static TOOLKIT_DLL void registerInstanceObjectMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        bool               handleNameParam, //true: extra input for handle name
        ReturnStyle        returnStyle, // how to return the output
        ObjMethod*         method);     // holds function pointer
    
    /** Deprecated - use registerStringMethod.
        Record data defining an addin which, if written in C++,
        correponds to a instance method (ie not static) which returns
        a string. */
    static TOOLKIT_DLL void registerInstanceStringMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        StringMethod*      method);     // holds function pointer
    
    /** Deprecated - use registerIntMethod.
        Record data defining an addin which, if written in C++,
        correponds to a instance method (ie not static) which returns
        an int. */
    static TOOLKIT_DLL void registerInstanceIntMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        IntMethod*         method);     // holds function pointer
    
    /** Deprecated - use registerBoolMethod.
        Record data defining an addin which, if written in C++,
        correponds to a instance method (ie not static) which returns
        a bool. */
    static TOOLKIT_DLL void registerInstanceBoolMethod(
        const string&      addinName,   // excel name
        const string&      category,    // category - for XL
        const string&      description, // what the addin does
        CClassConstSP      addinDataClass,  /* identifies class 
                                               representing addin data */
        BoolMethod*         method);     // holds function pointer
    
    /** Same as above but addinName is defaulted to the class's name so
        eg EDR_MyObject and the description is defaulted to "creates an "
        "object of type ..." */
    static TOOLKIT_DLL void registerConstructor(
        const string&      category,    // category - for XL
        CClassConstSP      addinDataClass);  /* identifies class 
                                                representing addin data */
    
    /** Look up the addin corresponding to the given string. This method
        has no knowledge of the EDR_ addin prefix */
    static TOOLKIT_DLL const Addin* lookUp(const string& addinName);

    /** Returns an array of names sorted into alphabetical order */
    static TOOLKIT_DLL CStringArraySP names();

    /** Indicates whether this addin method returns a native type or not */
    TOOLKIT_DLL bool isReturnNative() const;

    /** Indicates whether this addin has an extra parameter for a handle
        name */
    TOOLKIT_DLL bool hasHandleNameParam() const;

    /** Returns the Addin::Method union for this addin (ie gives function
        pointer for method) */
    TOOLKIT_DLL Method getMethod() const;

    /** Returns a description of what the addin does */
    TOOLKIT_DLL const string& getDescription() const;

    /** Returns category for the addin */
    TOOLKIT_DLL const string& getCategory() const;

    /** Returns the type of parameter that this addin returns */
    TOOLKIT_DLL CClassConstSP getReturnType() const;

    /** Return an array of fields representing each parameter in the clazz -
        this includes those for all parent fields. Fields are returned in
        a particular order */
    TOOLKIT_DLL static CFieldArray getDataClassFields(CClassConstSP clazz);
    
    /** Return an array of fields representing each parameter of the addin
        Fields are returned in a particular order */
    TOOLKIT_DLL CFieldArray getDataClassFields() const;
    
    /** Returns the class identifying the parameters that the addin takes */
    TOOLKIT_DLL CClassConstSP  getDataClass() const;

    /** Overrides "EDR " with "newPrefix ". Needs to be called before
        loadAllClasses() */
    static TOOLKIT_DLL void overrideCategoryPrefix(const string& newPrefix);

    /** typedef for kit registration function. */
    typedef void (KitRegister)(
        const string&    addinName,
        const string&    category,
        CClassConstSP    dataClass,   // identifies class representing addin
        const CFieldArray& params,    // description for each parameter
        bool             handleNameParam, // true: extra input for handle name
        ReturnStyle      returnStyle,  // how to return the output
        bool             returnIsNative, /* true - parameter is returned
                                            as an atomic type */
        MethodType       methodType,   /* constructor/instanceMethod/
                                          classMethod */
        Method           method,       // holds function pointer
        CClassConstSP    returnType);  // the type of the return parameter

    /** The supplied method will be called once for each addin registered.
        If kitMethod == 0, then no method is called. Note that this must be
        invoked afer all the classes are loaded (as otherwise no addins are
        registered) */
    static TOOLKIT_DLL void registerKit(KitRegister *kitMethod);

private:
    string                 addinName;    // excel name
    string                 category;     // category - for XL
    string                 description;  // what the addin does
    CClassConstSP          dataClass;   // identifies class representing addin
    bool                   handleNameParam; // true: extra input for handle name
    ReturnStyle            returnStyle;  // how to return the output
    MethodType             methodType; // constructor/instanceMethod/classMethod
    Method                 method;       // holds function pointer
    VirtualDestructorBaseSP objToFree;    // holds the functor
    bool                   returnIsNative; /* true - parameter is returned
                                             as an atomic type */
    CClassConstSP          returnType;   // the type of the return parameter

    ~Addin();

    TOOLKIT_DLL Addin(
        const string&      addinName,    // excel name
        const string&      category,     // category - for XL
        const string&      description,  // what the addin does
        CClassConstSP      dataClass,   // identifies class representing addin
        bool               handleNameParam, /* true: extra input for
                                               handle name */
        ReturnStyle        returnStyle,  // how to return the output
        MethodType         methodType,   /* constructor/instanceMethod/
                                            classMethod */
        Method             method,       // holds function pointer
        VirtualDestructorBase* objToFree,    // delete this
        bool               returnIsNative, /* true - parameter is returned
                                              as an atomic type */
        CClassConstSP      returnType);  /* the type of the return
                                            parameter */
    /** Add the given addin to hashtable containing all addins. Fails if
        addin with given name exists already. Takes ownwership of memory.
    Need to export as used by templated functions above */
    static TOOLKIT_DLL void addToHashtable(Addin*   addin);    
        
};


DRLIB_END_NAMESPACE
#endif
