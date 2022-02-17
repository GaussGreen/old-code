//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : An enum wrapped up as CObject
//
//   Author      : Mark A Robson
//
//   Date        : 14 Sep 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_BOXEDENUM_HPP
#define QLIB_BOXEDENUM_HPP

#include "edginc/Enum.hpp"
#include "edginc/Class.hpp"

DRLIB_BEGIN_NAMESPACE

template <class EnumType> class BoxedEnum: public Enum{
    typedef BoxedEnum<EnumType> Self;
public:
    static CClassConstSP const TYPE;

    virtual ~BoxedEnum(){}

    /** Create a wrapped enum of the specified enumType and with specified
        enumValue */
    static smartPtr<Self> create(EnumType enumVal){
        return smartPtr<Self>(new BoxedEnum(enumVal));
    }

    //// returns the wrapped enum's value
    EnumType enumValue() const{
        return (EnumType)enumVal;
    }

    //// Returns the string corresponding to this enum for this. Useful
    //// for error messages
    static const string& toString(EnumType enumVal){
        return TYPE->getEnumValue((int)enumVal).valueAsString;
    }

private:
    DLL_FIX_FOR_TEMPLATE_TYPE; // defines DLL_TYPE on windows with DLLs
    // no implementation here, but rather each actual type produces a
    // specialised instance
    static void load(CClassSP& classToLoad);

    static void loadCommon(CClassSP& classToLoad){
        classToLoad->setNative();
        classToLoad->setIsEnum();
        REGISTER(Self, classToLoad);
        SUPERCLASS(Enum);
        EMPTY_SHELL_METHOD(defaultConstructor);
    };

    static IObject* defaultConstructor(){
        return new BoxedEnum();
    }

    BoxedEnum(): Enum(TYPE){}
    BoxedEnum(EnumType enumVal): Enum(TYPE, enumVal){}
};

#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
//// work around problem of trying to export static fields of templates
//// across dlls. For each template used there needs to be an
//// explicit specialisation of TYPE defined which specifies the name and
//// load method
template <class X> CClassConstSP const BoxedEnum<X>::TYPE =
CClass::templateRegisterClass(typeid(BoxedEnum<X>));
#endif

/** Specialisation for enums. */
template<class T> class FieldInLineMethods<T, true  /*is an enum*/>{
public:
    static IObjectSP get(T* t){
        return BoxedEnum<T>::create(*t);
    }
    static void set(T* t, IObjectSP obj){
        const BoxedEnum<T>* newEnum = 
            DYNAMIC_CONST_CAST(BoxedEnum<T>, obj.get());
        *t = newEnum->enumValue();
    }
};

///////// Magic macros for type registration                            ///////
///////// See example in Enum.cpp for use. Note that it needs to be //////
///////// used within a source file and outside of a class definition.  ///////

//// ENUM_TYPE is the C++ name of the enum. ENUM_DESC is a comment describing
//// what the enum represents. This macro does not declare the enum to be a
//// 'public' type (ie one available to clients/spreadsheets)
#define START_ENUM_DEFINITION(ENUM_TYPE, ENUM_DESC) \
  template <> void BoxedEnum<ENUM_TYPE>::load(CClassSP& classToLoad){ \
      loadCommon(classToLoad); \
      classToLoad->setDescription(ENUM_DESC);

//// Same as START_ENUM_DEFINITION but declares the type to be public
#define START_PUBLIC_ENUM_DEFINITION(ENUM_TYPE, ENUM_DESC) \
  START_ENUM_DEFINITION(ENUM_TYPE, ENUM_DESC) \
  classToLoad->setPublic(); // make visible to EAS/spreadsheet

//// This macro can only be used when the enum is not declared within a class
#define ENUM_VALUE(ENUM_VALUE, COMMENT) \
  classToLoad->setEnumValue(ENUM_VALUE, #ENUM_VALUE, COMMENT);

#define ENUM_VALUE_AND_NAME(ENUM_VALUE, ENUM_NAME, COMMENT) \
  classToLoad->setEnumValue(ENUM_VALUE, ENUM_NAME, COMMENT);


#if defined (_MSC_VER) && defined(QLIB_BUILD_DLL)
#define END_ENUM_DEFINITION(ENUM_TYPE) \
  } \
  template <> CClassConstSP const BoxedEnum<ENUM_TYPE> ::DLL_TYPE = \
  CClass::registerClassWithAliasLoadMethod(#ENUM_TYPE, typeid(ENUM_TYPE), \
                                           typeid(BoxedEnum<ENUM_TYPE>),\
                                           load);
#else
#define END_ENUM_DEFINITION(ENUM_TYPE) \
  } \
  template <> CClassConstSP const BoxedEnum<ENUM_TYPE> ::TYPE = \
  CClass::registerClassWithAliasLoadMethod(#ENUM_TYPE, typeid(ENUM_TYPE), \
                                           typeid(BoxedEnum<ENUM_TYPE>), \
                                           load);
#endif  

DRLIB_END_NAMESPACE
#endif

