/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/jni/jtypes.h
// Purpose:     Helpers for working with Java types
// Author:      Vadim Zeitlin
// Created:     06.06.03
// RCS-ID:      $Id: jtypes.h,v 1.8 2004/10/05 09:13:38 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/jni/jtypes.h
    @brief  This header defines some helpful macros for working with Java types.

    This is a semi-private header and is normally only used by ito33/jni.h. You
    may, however, use it for your own purposes as well but be warned that the
    macros here are a bit tricky.

    Basically our goal here is to allow doing something at once for all of the
    Java types and for this we need quite a lot of helpers -- defined here.

    The code here is based on the ideas from the "C++ Users Journal" Jan 2001
    issue article "JNI-C++ Integration Made Easy" by Lev Finkelstein and
    Evgeniy Gabrilovich.
 */

#ifndef _ITO33_JTYPES_H_
#define _ITO33_JTYPES_H_

#include <jni.h>

#include "ito33/longlong.h"

namespace ito33
{

namespace JNI
{

/**
    This namespace contains stuff for working with Java types.
 */
namespace Types
{

/// Make an array type name from Java type name
#define JNI_TYPE_MAKE_ARRAY(jtype)   jtype ## Array

/**
    Macro to define traits-like structs for all built-in types.

    Each struct is named after the corresponding primitive type except for the
    first letter which is capitalized and defines:
        - JavaType         Java type corresponding to this type
        - ArrayType        Java array type corresponding to this type
        - signature        type signature for this type
        - array_signature  array type signature for it
 */
#define JNI_TYPE_DEFINE(Type, jtype, sig)                                     \
  struct Type                                                               \
  {                                                                         \
    typedef jtype JavaType;                                               \
    typedef JNI_TYPE_MAKE_ARRAY(jtype) ArrayType;                         \
    static const char *signature() { return sig; }                        \
    static const char *array_signature() { return "[" sig; }              \
  }

/// Definition for jboolean
JNI_TYPE_DEFINE(Boolean, jboolean, "Z");

/// Definition for jbyte
JNI_TYPE_DEFINE(Byte, jbyte, "B");

/// Definition for jchar
JNI_TYPE_DEFINE(Char, jchar, "C");

/// Definition for jshort
JNI_TYPE_DEFINE(Short, jshort, "S");

/// Definition for jint
JNI_TYPE_DEFINE(Int, jint, "I");

/// Definition for jlong
JNI_TYPE_DEFINE(Long, jlong, "J");

/// Definition for jfloat
JNI_TYPE_DEFINE(Float, jfloat, "F");

/// Definition for jdouble
JNI_TYPE_DEFINE(Double, jdouble, "D");

/**
    Type corresponding to Java string class.

    This is not a primitive type but so commonly used that we make an
    exception for it

    Note that as jstringArray doesn't exist, we can't use JNI_TYPE_DEFINE
    here.
 */
struct String
{
  /// the Java type for this class
  typedef jstring JavaType;

  /// the signature of fields of this type
  static const char *signature() { return "Ljava/lang/String;"; }
};


/**
    Returns the Java type corresponding to this type name.

    @param Type name of the type, i.e. Boolean, Int, Long, Double, ...
 */
#define JNI_JTYPE(Type) ::ito33::JNI::Types::Type::JavaType

/**
    Returns the Java array type corresponding to this type name.

    @param Type name of the type, i.e. Boolean, Int, Long, Double, ...
 */
#define JNI_ATYPE(Type)  ::ito33::JNI::Types::Type::ArrayType

/**
    Returns the signature for the given type name.

    @param Type name of the type, i.e. Boolean, Int, Long, Double, ...
 */
#define JNI_SIGNATURE(Type)   ::ito33::JNI::Types::Type::signature()

/**
    Returns the signature for the array of values of the given type name.

    @param Type name of the type, i.e. Boolean, Int, Long, Double, ...
 */
#define JNI_ASIGNATURE(Type) ::ito33::JNI::Types::Type::array_signature()


/**
    Instantiate the given macro for all Java primtive types.

    Notice that BLOCK_MACRO parameter must take the type name (e.g. Boolean,
    Int, Long, Double, ...) as argument. Note that, unusually, this macro
    should contain any necessary semi-colons inside it and also there is no
    need to put a semi-colon after this macro itself.

    @param BLOCK_MACRO the macro to instantiate, must take type as argument
 */
#define JNI_DO_FOR_ALL_PRIMITIVE_TYPES(BLOCK_MACRO)                           \
  BLOCK_MACRO(Boolean)                                                      \
  BLOCK_MACRO(Byte)                                                         \
  BLOCK_MACRO(Char)                                                         \
  BLOCK_MACRO(Short)                                                        \
  BLOCK_MACRO(Int)                                                          \
  BLOCK_MACRO(Long)                                                         \
  BLOCK_MACRO(Float)                                                        \
  BLOCK_MACRO(Double)

} // namespace Types

/**
    This template allows to map a C++ type to a JNI type.

    It is specialized for all primitive C++ types having Java equivalents and
    has the following fields:
        - ctype C type itself (just for convenience)
        - jtype JNI type corresponding to ctype
        - atype JNI type of array with jtype elements
 */
template <typename T>
struct CTypeTraits
{
    typedef T ctype;
    typedef jobject jtype;
    typedef jobjectArray atype;
};

/**
    This template allows to map a JNI type to a C++ one.

    It is specialized for all scalar JNI types except jobject.

    @sa CTypeTraits
 */
template <typename T>
struct JTypeTraits;

#define JNI_DEFINE_TYPE_TRAITS(cpptype, jnitype)                              \
    template <>                                                               \
    struct CTypeTraits<cpptype>                                               \
    {                                                                         \
        typedef cpptype ctype;                                                \
        typedef jnitype jtype;                                                \
        typedef JNI_TYPE_MAKE_ARRAY(jnitype) atype;                           \
    };                                                                        \
                                                                              \
    template <>                                                               \
    struct JTypeTraits<jnitype>                                               \
    {                                                                         \
        typedef cpptype ctype;                                                \
        typedef jnitype jtype;                                                \
        typedef JNI_TYPE_MAKE_ARRAY(jnitype) atype;                           \
    }

JNI_DEFINE_TYPE_TRAITS(bool,            jboolean);
JNI_DEFINE_TYPE_TRAITS(unsigned char,   jbyte);
JNI_DEFINE_TYPE_TRAITS(unsigned short,  jchar);
JNI_DEFINE_TYPE_TRAITS(short,           jshort);
JNI_DEFINE_TYPE_TRAITS(int,             jint);
JNI_DEFINE_TYPE_TRAITS(ito33::LongLong, jlong);
JNI_DEFINE_TYPE_TRAITS(float,           jfloat);
JNI_DEFINE_TYPE_TRAITS(double,          jdouble);

#undef JNI_DEFINE_TYPE_TRAITS

} // namespace JNI

} // namespace ito33

#endif // _ITO33_JTYPES_H_

