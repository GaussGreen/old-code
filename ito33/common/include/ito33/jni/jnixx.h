/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/jni/jnixx.h
// Purpose:     Utility classes and functions for C++ JNI code
// Author:      Vadim Zeitlin
// Created:     06.06.03
// RCS-ID:      $Id: jnixx.h,v 1.32 2004/10/15 18:35:10 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/jni/jnixx.h
    @brief  Defines helper class and functions for C++ JNI code.

    JNI functions defined in JDK are inherently low level because they're
    C-oriented. We define here many useful classes which makes writing JNI code
    in C++ much simpler.

    Most notably, all Java exceptions are mapped to C++ exceptions provided
    that you use JNI::Env instead of the "raw" JNIEnv pointer. If you do,
    it is enough to take the entire body of all your JNI functions inside @c
    try/catch block and silently ignore any JavaException objects raised inside
    it -- if this happens it means that a Java exception is pending and the
    only thing the JNI function can do is to return as soon as possible to Java
    so that the exception could be handled there.

    There are also classes which drastically simplify working with the field
    and method ids as well as class references and much more.

    Note that all classes in this file are defined in ito33::JNI namespace.
 */
#ifndef _ITO33_JNIXX_H_
#define _ITO33_JNIXX_H_

#include "ito33/debug.h"
#include "ito33/exception.h"

#include "ito33/string.h"

#include "ito33/jni/jtypes.h"

#include <exception>

namespace ito33
{

/**
    All stuff useful for writing JNI code.
 */
namespace JNI
{

// forward declarations for JNI::Env
class JArrayAny;
class JStringChars;

template <typename T> class JavaIdAny;
template <typename T> class JavaStaticIdAny;
typedef JavaIdAny<jfieldID> FieldIdAny;
typedef JavaIdAny<jmethodID> MethodIdAny;
typedef JavaStaticIdAny<jfieldID> StaticFieldIdAny;
typedef JavaStaticIdAny<jmethodID> StaticMethodIdAny;

template <typename T> class FieldId;
template <typename T> class StaticFieldId;


/**
    An exception of this kind is thrown whenever a Java exception occurs.

    The class itself doesn't contain anything as the only thing we really
    want to do is to just return to Java code anyhow -- a Java exception
    had been already raised and will be handled in Java when the native
    method returns.
 */
class JavaException { };

// ----------------------------------------------------------------------------
// Env
// ----------------------------------------------------------------------------

/**
    JNIEnv wrapper with error checking.

    This class provides exactly the same functions as a JNIEnv pointer but
    does error checking for each call and throws a JavaException if an
    error occurs.
 */
class Env
{
protected:
  /**
      Throw a C++ exception unconditionally.

      Having @c throw only here will allow to modify it easier if we want
      to throw another exception or do something else.
   */
  void Throw() { throw JavaException(); }

  /**
      Throw a C++ exception if there is a pending Java exception.

      This function should be called after calling a JNI function which
      doesn't have an error status code. For the functions which do have
      one, it should @b not be used because it is more expensive than
      simply checking the return code.
   */
  void ThrowIfException()
  {
    if ( m_env->ExceptionCheck() )
      Throw();
  }

  /**
      Checks the given value and throws a C++ exception if it is @c NULL.

      This is just a convenient short hand for the operation of calling a
      Java function, storing its return value in a temporary variable,
      checking whether it is @c NULL and throwing an exception if it is and
      returning it to the caller otherwise.

      @param rc the return code of a JNI function
   */
  template <typename T>
  T CheckIfFailed(T rc)
  {
    if ( !rc )
      Throw();

    return rc;
  }

  /**
      Checks if there is a pending Java exception and throws a C++ exception
      if there is, simply returns the given value otherwise.

      This is similar to CheckIfFailed() but is meant to be used with the JNI
      functions without the error return code.

      @sa ThrowIfException

      @param rc the return code of a JNI function
   */
  template <typename T>
  T CheckIfException(T rc)
  {
    ThrowIfException();

    return rc;
  }

public:
  /**
      Constructor associates this object with the real JNIEnv.

      Will throw if the env pointer is @c NULL.

      @param env pointer to JNI environment, must not be @c NULL
   */
  Env(JNIEnv *env) : m_env(env)
  {
    if ( !env )
      Throw();
  }

  /// Implicit conversion to JNIEnv
  operator JNIEnv *() const { return m_env; }

  /// @name Miscellaneous JNI functions
  //@{

  jfieldID GetFieldID(jclass clazz, const char *name, const char *sig)
  {
    return CheckIfFailed(m_env->GetFieldID(clazz, name, sig));
  }

  jmethodID GetMethodID(jclass clazz, const char *name, const char *sig)
  {
    return CheckIfFailed(m_env->GetMethodID(clazz, name, sig));
  }

  jfieldID GetStaticFieldID(jclass clazz, const char *name, const char *sig)
  {
    return CheckIfFailed(m_env->GetStaticFieldID(clazz, name, sig));
  }

  jmethodID GetStaticMethodID(jclass clazz, const char *name, const char *sig)
  {
    return CheckIfFailed(m_env->GetStaticMethodID(clazz, name, sig));
  }

  jclass GetObjectClass(jobject obj)
  {
    return CheckIfFailed(m_env->GetObjectClass(obj));
  }

  jclass FindClass(const char *name)
  {
    return CheckIfFailed(m_env->FindClass(name));
  }

  jobject NewObject(jclass clazz, jmethodID methodID, ...)
  {
      va_list args;
      va_start(args, methodID);

      jobject rc = m_env->NewObjectV(clazz, methodID, args);

      va_end(args);

      return CheckIfException(rc);
  }

  jstring NewStringUTF(const char *bytes)
  {
      return CheckIfFailed(m_env->NewStringUTF(bytes));
  }

  void DeleteLocalRef(jobject obj)
  {
    // According to JNI Specification, this function never generates
    // exceptions, so we don't check for them.
    //
    // Note that even if it did, we wouldn't want to throw a C++ exception
    // from here because we may be (and often are) called from dtors of
    // local objects and dtors must not throw.
    m_env->DeleteLocalRef(obj);
  }

  void DeleteWeakGlobalRef(jobject obj)
  {
    // see above: the same considerations apply here
    m_env->DeleteWeakGlobalRef(obj);
  }

  jweak NewWeakGlobalRef(jobject obj)
  {
    // NewWeakGlobalRef() may return NULL even if no exception is pending
    // (but if obj is NULL or is a deleted weak global ref), so we can't
    // use CheckIfFailed() here
    return CheckIfException(m_env->NewWeakGlobalRef(obj));
  }


  const char *GetStringUTFChars(jstring jstr, jboolean *isCopy)
  {
    return CheckIfFailed(m_env->GetStringUTFChars(jstr, isCopy));
  }

  void ReleaseStringUTFChars(jstring jstr, const char *data)
  {
    m_env->ReleaseStringUTFChars(jstr, data);
  }

  //@}

  /// @name Calling Java methods
  //@{

#define JNI_DECLARE_CALLMETHOD(Type)                                          \
  JNI_JTYPE(Type) Call ## Type ## Method(jobject obj, jmethodID methodID, ...)\
  {                                                                           \
    va_list args;                                                             \
    va_start(args, methodID);                                                 \
    JNI_JTYPE(Type) rc = m_env->Call ## Type ## MethodV(obj, methodID, args); \
    va_end(args);                                                             \
                                                                              \
    return CheckIfException(rc);                                              \
  }

  JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DECLARE_CALLMETHOD)

#undef JNI_DECLARE_CALLMETHOD

  void CallVoidMethod(jobject obj, jmethodID methodID, ...)
  {
    va_list args;
    va_start(args, methodID);
    m_env->CallVoidMethodV(obj, methodID, args);
    va_end(args);

    ThrowIfException();
  }

  //@}

  /// @name Accessing the fields
  //@{

  /**
      This macro defines 8 set/get functions for each primitive type.

      It defines the normal Set/Get versions for the instance and class
      fields taking either jobject or jclass as the first argument and a
      jfieldID as the second one. But these methods are defined as private
      because they are type-unsafe: you may pass a jfieldID of an integer
      field to GetDoubleField(), for example.

      For type-safe access this macro also defines the same functions but
      taking the corresponding FieldId template specialization which allows
      for better type checking during compile-time.
   */
#define JNI_DECLARE_SETGET_FIELD(Type)                                        \
protected:                                                                    \
  JNI_JTYPE(Type) Get ## Type ## Field(jobject obj, jfieldID fid)           \
  {                                                                         \
    return CheckIfException(m_env->Get ## Type ## Field(obj, fid));       \
  }                                                                         \
                                       \
  JNI_JTYPE(Type) GetStatic ## Type ## Field(jclass cls, jfieldID fid)      \
  {                                                                         \
    return CheckIfException(m_env->GetStatic ## Type ## Field(cls, fid)); \
  }                                                                         \
                                       \
  void Set ## Type ## Field(jobject obj, jfieldID fid, JNI_JTYPE(Type) val) \
  {                                                                         \
    m_env->Set ## Type ## Field(obj, fid, val);                           \
                                       \
    ThrowIfException();                                                   \
  }                                                                         \
                                       \
  void                                                                      \
  SetStatic ## Type ## Field(jclass cls, jfieldID fid, JNI_JTYPE(Type) val) \
  {                                                                         \
    m_env->SetStatic ## Type ## Field(cls, fid, val);                     \
                                       \
    ThrowIfException();                                                   \
  }                                                                         \
                                       \
public:                                                                       \
  JNI_JTYPE(Type)                                                           \
  Get ## Type ## Field(jobject obj, FieldId<JNI_JTYPE(Type)>& fid);         \
                                       \
  JNI_JTYPE(Type)                                                           \
  GetStatic ## Type ## Field(jclass cls, StaticFieldId<JNI_JTYPE(Type)>& fid);    \
                                       \
  void                                                                      \
  Set ## Type ## Field(jobject obj,                                         \
            FieldId<JNI_JTYPE(Type)>& fid,                       \
            JNI_JTYPE(Type) value);                              \
                                       \
  void                                                                      \
  SetStatic ## Type ## Field(jclass cls,                                    \
               StaticFieldId<JNI_JTYPE(Type)>& fid,                 \
               JNI_JTYPE(Type) value);

  JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DECLARE_SETGET_FIELD)

#undef JNI_DECLARE_SETGET_FIELD

  /**
      Accessor for the other (than the primitive ones) types of fields.
   */
  jobject GetObjectField(jobject obj, jfieldID fid)
  {
    return CheckIfFailed(m_env->GetObjectField(obj, fid));
  }

  jobject GetStaticObjectField(jclass cls, jfieldID fid)
  {
    return CheckIfFailed(m_env->GetStaticObjectField(cls, fid));
  }

  jobject GetObjectField(jobject obj, FieldIdAny& fid);
  jobject GetStaticObjectField(jclass cls, StaticFieldIdAny& fid);

  void SetObjectField(jobject obj, jfieldID fid, jobject val)
  {
    m_env->SetObjectField(obj, fid, val);

    ThrowIfException();
  }

  void SetStaticObjectField(jclass cls, jfieldID fid, jobject val)
  {
    m_env->SetStaticObjectField(cls, fid, val);

    ThrowIfException();
  }

  JStringChars GetStringField(jobject obj, jfieldID fid);

  //@}

  /// @name Working with the arrays
  //@{

  jsize GetArrayLength(jarray array)
  {
    // JNI Specification says that no exceptions are possible here
    return m_env->GetArrayLength(array);
  }

  jobjectArray NewObjectArray(jsize length, jclass cls, jobject init)
  {
    return CheckIfFailed(m_env->NewObjectArray(length, cls, init));
  }

  jobject GetObjectArrayElement(jobjectArray arr, jsize n)
  {
    return CheckIfFailed(m_env->GetObjectArrayElement(arr, n));
  }

  void SetObjectArrayElement(jobjectArray arr, jsize n, jobject value)
  {
    m_env->SetObjectArrayElement(arr, n, value);

    ThrowIfException();
  }

  // we declare below all the standard JNI methods for working with the
  // arrays and also SetArray() function.
  //
  // NB: Get/Set<Type>ArrayElements methods are only used by JArray because
  //     they are not type safe
protected:
#define JNI_DECLARE_ARRAY_FUNCTIONS(Type)                                    \
  JNI_JTYPE(Type) *Get ## Type ## ArrayElements(JNI_ATYPE(Type) array,      \
                         jboolean *isCopy = NULL)    \
  {                                                                         \
    return                                                                \
      CheckIfFailed(m_env->Get ## Type ## ArrayElements(array, isCopy));\
  }                                                                         \
                                       \
  void Release ## Type ## ArrayElements(JNI_ATYPE(Type) array,              \
                     JNI_JTYPE(Type) *elems,             \
                     jint mode)                          \
  {                                                                         \
    /* JNI Specification says that no exceptions are possible here */     \
    m_env->Release ## Type ## ArrayElements(array, elems, mode);          \
  }                                                                         \
                                       \
public:                                                                       \
  JNI_ATYPE(Type) New ## Type ## Array(jsize size)                          \
  {                                                                         \
    return CheckIfFailed(m_env->New ## Type ## Array(size));              \
  }                                                                         \
                                       \
  void Get ## Type ## ArrayRegion(JNI_ATYPE(Type) array,                    \
                  jsize start,                              \
                  jsize len,                                \
                  JNI_JTYPE(Type) *data)                    \
  {                                                                         \
    m_env->Get ## Type ## ArrayRegion(array, start, len, data);           \
                                       \
    ThrowIfException();                                                   \
  }                                                                         \
                                       \
  void Set ## Type ## ArrayRegion(JNI_ATYPE(Type) array,                    \
                  jsize start,                              \
                  jsize len,                                \
                  JNI_JTYPE(Type) *data)                    \
  {                                                                         \
    m_env->Set ## Type ## ArrayRegion(array, start, len, data);           \
                                       \
    ThrowIfException();                                                   \
  }

  JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DECLARE_ARRAY_FUNCTIONS)

#undef JNI_DECLARE_ARRAY_FUNCTIONS

  /**
      Sets the value of an array object member field from C data.

      This function conveniently combines calls to CreateFillArray() and
      SetObjectField(). We also check for NULL pointers to avoid memory
      crashes.

      @param obj the object whose field we want to set
      @param fid the id of the field to be set
      @param nItems the number of items in pItems array
      @param pItems pointer to the array of (at least) nItems items
   */
  template <typename T>
  void SetArray(jobject obj, jfieldID fid, size_t nItems, T *pItems)
  {
   if ( pItems )
    SetObjectField(obj, fid, CreateFillArray(*this, nItems, pItems));
  }

  //@}

  /**
      @name Exceptions handling.

      Note that ConvertStdExceptionToJava() and ConvertExceptionToJava()
      functions are deprecated and kept only for backwards compatibility with
      Freebound2k Java API. New Java code uses JavaThrow().
   */
  //@{

  /**
      Throw a Java exception of given class.

      As this method is usually called from a catch handler, it never throws a
      C++ exception itself or lets one propagate through it.

      @param env pointer to JNI environment, must not be @c NULL
      @param clsname name of the Java exception class.
      @param msg the detail message for Java class ctor
   */
  static
  void JavaThrow(JNIEnv *env, const char *clsname, const std::string& msg);

  /**
      Returns the name of the Java exception class we throw.

      This is the class from com.ito33.common package which is used as the base
      class for all Java exceptions thrown by native code.
   */
  static const char *GetExceptionClassName()
  {
    return "Lcom/ito33/common/Exception;";
  }

  /**
      Returns the name of the Java exception class we used to throw before.

      This function is only used by deprecated ConvertStdExceptionToJava() and
      ConvertExceptionToJava() below.
   */
  static const char *GetOldExceptionClassName()
  {
    return "Lcom/ito33/utils/DLLErrorException;";
  }

  /**
      Throw a Java exception corresponding to the given standard one.

      The exception thrown is of the class GetExceptionClassName().

      @param env pointer to JNI environment, must not be @c NULL
      @param e a standard exception
   */
  static void JavaThrow(JNIEnv *env, const std::exception& e)
  {
    JavaThrow(env, GetExceptionClassName(), e.what());
  }

  /**
      Throw a Java exception with the associated message.

      @param env pointer to JNI environment, must not be @c NULL
      @param msg the message which will be "seen" by Java code
   */
  static void JavaThrow(JNIEnv *env, const std::string& msg)
  {
    JavaThrow(env, GetExceptionClassName(), msg);
  }

#ifdef _MSC_VER
  /**
      Overload of JavaThrow(msg) needed for some broken compilers.

      Standard library shipped with VC 7.1 adds a non-standard std::exception
      ctor taking "const char *" which results in ambiguities when writing
      JavaThrow(env, "foo") unless we provide this overload as well.
   */
  static void JavaThrow(JNIEnv *env, const char *msg)
  {
    JavaThrow(env, std::string(msg));
  }
#endif // _MSC_VER

  /**
      Throw a Java exception with a standard "unexpected" message.

      This function exists solely to avoid writing the same message all the
      time in all "catch ( ... )" handlers. It should never be called in normal
      circumstances, use it only for truely unexpected cases.

      @param env pointer to JNI environment, must not be @c NULL
   */
  static void JavaThrowUnexpected(JNIEnv *env)
  {
    JavaThrow(env, "Unexpected exception in C++ JNI code, please report as bug.");
  }

  /**
      Throw a Java exception corresponding to our C++ exception.

      The exception thrown is of the class GetExceptionClassName().

      @param env pointer to JNI environment, must not be @c NULL
      @param e an exception object
   */
  static void JavaThrow(JNIEnv *env, const ito33::Exception& e)
  {
    JavaThrow(env, GetExceptionClassName(), e.GetFullMessage());
  }

  /**
      Convert a standard C++ exception to Java one.

      This method is deprecated, use one of JavaThrow() overloads instead.

      This should be used to detect the exceptions inside the standard
      functions. It should be called from the catch handler for
      std::exception.

      @param env pointer to JNI environment, must not be @c NULL
      @param msg the message identifying the source of the exception
                 (e.what() will be appended to it)
      @param e the standard exception we caught
   */
  static void ConvertStdExceptionToJava(JNIEnv *env,
                                        const char *msg,
                                        std::exception& e)
  {
    JavaThrow(env, GetOldExceptionClassName(), std::string(msg) + e.what());
  }

  /**
      Convert our exception to Java one.

      This method is deprecated, use one of JavaThrow() overloads instead.
   */
  static void ConvertExceptionToJava(JNIEnv *env,
                                     const char *msg,
                                     ito33::Exception& e)
  {
    JavaThrow(env, GetOldExceptionClassName(),
                std::string(msg) + e.GetFullMessage());
  }

  //@}

private:
  JNIEnv *m_env;

  // give it access to our array methods
  friend class JArrayAny;
};

// ----------------------------------------------------------------------------
// JClass
// ----------------------------------------------------------------------------

/**
    This class wraps Java class reference.

    Using it allows to get a reference for the class from object with
    minimal amount of code.
 */
class JClass
{
public:
  /**
      Constructs the class reference from the class name.

      May throw JavaException if an error occurs.

      @param env the checked JNIEnv object
      @param name the name of the class to find
   */
  JClass(Env& env, const char *name)
    : m_env(env), m_cls(env.FindClass(name))
  {
  }

  /**
      Constructs the class reference from the given object.

      May throw JavaException if an error occurs.

      @param env the checked JNIEnv object
      @param obj the object whose class we're interested in
   */
  JClass(Env& env, jobject obj)
    : m_env(env), m_cls(env.GetObjectClass(obj))
  {
  }

  /**
      Destructor destroys the reference to the class.

      Although JVM will destroy the reference anyhow when the native method
      returns, doing it as soon as possible is preferred because JVM has a
      rather low limit on the total number of local references and we may
      overflow it.
   */
  ~JClass()
  {
    m_env.DeleteLocalRef(m_cls);
  }

  // default copy ctor and assignment operator are ok

  /// Implicit conversion to jclass
  operator jclass() const { return m_cls; }

private:
  Env& m_env;
  jclass m_cls;

  NO_COPY_CLASS(JClass);
};

// ============================================================================
// FieldId<> and MethodId<>
// ============================================================================

// FieldIdAny and MethodIdAny classes are very similar -- so similar, in fact,
// that we use the same template (JavaIdAny<>) to implement both of them and
// simply instantiate for 2 different (dummy) types
//
// Implementation details: as we can't specialize just a single function of a
// class, we need to use a separate template function in JavaIdAny. Moreover,
// due to a bug in VC++ 6 which makes it impossible to explicitly specialize
// function templates (independently of the number of instantiations the first
// one is going to be used!), we have to use a helper class here.

// ----------------------------------------------------------------------------
// JavaIdAny
// ----------------------------------------------------------------------------

/**
    Helper class for JavaIdAny.

    This just allows us to call either GetFieldID() or GetMethodID() in
    JavaIdAny.
 */
template <typename T>
struct JavaId
{
  /// Return the value of the given field or method id
  static T Get(Env& env, jclass cls,
        const char *name, const char *signature);
};

/**
    Helper class for JavaStaticIdAny.

    This just allows us to call either GetFieldID() or GetMethodID() in
    JavaStaticIdAny.
 */
template <typename T>
struct JavaStaticId
{
  /// Return the value of the given field or method id
  static T Get(Env& env, jclass cls,
        const char *name, const char *signature);
};

/// specialize JavaId for jfieldID
template <>
struct JavaId<jfieldID>
{
  /**
      Return the field id for the given field.

      @param env the environment
      @param cls the Java class
      @param name the name of the field
      @param signature the signature describing the type of the field
   */
  static jfieldID Get(Env& env, jclass cls,
            const char *name, const char *signature)
  {
    return env.GetFieldID(cls, name, signature);
  }
};

/// specialize JavaId for jmethodID
template <>
struct JavaId<jmethodID>
{
  /**
      Return the id for the given method.

      @param env the environment
      @param cls the Java class
      @param name the name of the method
      @param signature the methods signature
   */
  static jmethodID Get(Env& env, jclass cls,
            const char *name, const char *signature)
  {
    return env.GetMethodID(cls, name, signature);
  }
};

/// specialize JavaId for jfieldID
template <>
struct JavaStaticId<jfieldID>
{
  /**
      Return the field id for the given field.

      @param env the environment
      @param cls the Java class
      @param name the name of the field
      @param signature the signature describing the type of the field
   */
  static jfieldID Get(Env& env, jclass cls,
            const char *name, const char *signature)
  {
    return env.GetStaticFieldID(cls, name, signature);
  }
};

/// specialize JavaId for jmethodID
template <>
struct JavaStaticId<jmethodID>
{
  /**
      Return the id for the given method.

      @param env the environment
      @param cls the Java class
      @param name the name of the method
      @param signature the methods signature
   */
  static jmethodID Get(Env& env, jclass cls,
            const char *name, const char *signature)
  {
    return env.GetStaticMethodID(cls, name, signature);
  }
};

/**
    This class represents a Java object field or method id for any type.

    For the primitive types and jstring, the corresponding FieldId or MethodId
    template class specialization should be used. For all the other types you
    need to use the base classes FieldIdAny or MethodIdAny directly and
    manually provide the Java field signature.
 */
template <typename T>
class JavaIdAny
{
protected:
  /**
      Get the id for the given class.

      Note that it is safe to call this method multiple times and even
      from multiple threads as we assume that the value returned by
      GetFieldID() will always be the same so the worst which can happen
      to us is that we call it (unnecessarily) twice or more...

      This method will throw JavaException if the id is not found.

      @param env pointer to the JNI environment
      @param cls the Java class whose member we're interested in
   */
  void Init(Env& env, jclass cls)
  {
    m_id = JavaId<T>::Get(env, cls, m_name, m_signature);
  }

  /**
      Get the id for the given object.

      This is just like Init() overload above except that it finds class from
      the object automatically -- this makes it more convenient to use as we
      usually have the object pointer, not the class pointer.

      @param env pointer to the JNI environment
      @param obj the Java object whose member we want to access
   */
  void Init(Env& env, jobject obj)
  {
    Init(env, JClass(env, obj));
  }

public:
  /// the type of ids we're working with
  typedef T IdType;

  /**
      Default ctor, must use Init() later.

      Note that the ctor arguments @b must be literal strings as it
      stores them internally (until Init() is called).

      @param name the method name
      @param signature the method signature as defined by Java rules
   */
  JavaIdAny(const char *name, const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    m_id = 0;
  }

  /**
      Ctor gets the id, throws an exception if failed.

      This constructor retrives the id directly from the class, it is more
      efficient than doing it from the object -- especially if you already
      have the class at hand.

      @param env pointer to the JNI environment
      @param cls the Java class whose method we're interested in
      @param name the member name
      @param signature the member signature as defined by Java rules
   */
  JavaIdAny(Env& env,
       jclass cls,
       const char *name,
       const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    Init(env, cls);
  }

  /**
      Ctor gets the id, throws an exception if failed.

      This constructor first retrieves the class from the object and only
      then gets the id so it is less efficient than the constructor taking
      the class argument directly if you already have the class pointer
      somewhere. But if you don't, this constructor might be more convenient
      to use.

      @param env pointer to the JNI environment
      @param obj the Java object whose method we're interested in
      @param name the member name
      @param signature the member signature as defined by Java rules
   */
  JavaIdAny(Env& env,
       jobject obj,
       const char *name,
       const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    Init(env, obj);
  }

  /**
      Initializes the id if necessary and returns it.

      This method allows to use JavaIdAny objects like this:
      @code
          static JavaIdAny s_fid("name", "signature");
          ...
          env.SetStaticIntField(obj, s_fid(env, cls), ...);
      @endcode
   */
  T operator()(Env& env, jclass cls)
  {
    if ( !m_id )
      Init(env, cls);

    return m_id;
  }

  /**
      Initializes the id if necessary and returns it.

      This method allows to use JavaIdAny objects like this:
      @code
          static JavaIdAny s_fid("name", "signature");
          ...
          env.SetIntField(obj, s_fid(env, obj), ...);
      @endcode
   */
  T operator()(Env& env, jobject obj)
  {
    if ( !m_id )
      Init(env, obj);

    return m_id;
  }

  /**
      Provides direct access to the id.

      This supposes that the id had been already initialized, no checks are
      done to ensure that this is indeed the case in release build but in
      debug build an assert is raised if an uninitialized id is used.
   */
  operator T() const
  {
    ASSERT_MSG( m_id, "Using uninitialized JavaIdAny" );

    return m_id;
  }

private:
  T m_id;
  const char * const m_name;
  const char * const m_signature;

  NO_COPY_CLASS(JavaIdAny);
};

/**
    This class represents a Java static field or method id for any type.

    For the primitive types and jstring, the corresponding FieldId or MethodId
    template class specialization should be used. For all the other types you
    need to use the base classes FieldIdAny or MethodIdAny directly and
    manually provide the Java field signature.
 */
template <typename T>
class JavaStaticIdAny
{
protected:
  /**
      Get the id for the given class.

      Note that it is safe to call this method multiple times and even
      from multiple threads as we assume that the value returned by
      GetStaticFieldID() will always be the same so the worst which can happen
      to us is that we call it (unnecessarily) twice or more...

      This method will throw JavaException if the id is not found.

      @param env pointer to the JNI environment
      @param cls the Java class whose member we're interested in
   */
  void Init(Env& env, jclass cls)
  {
    m_id = JavaStaticId<T>::Get(env, cls, m_name, m_signature);
  }

  /**
      Get the id for the given object.

      This is just like Init() overload above except that it finds class from
      the object automatically -- this makes it more convenient to use as we
      usually have the object pointer, not the class pointer.

      @param env pointer to the JNI environment
      @param obj the Java object whose member we want to access
   */
  void Init(Env& env, jobject obj)
  {
    Init(env, JClass(env, obj));
  }

public:
  /// the type of ids we're working with
  typedef T IdType;

  /**
      Default ctor, must use Init() later.

      Note that the ctor arguments @b must be literal strings as it
      stores them internally (until Init() is called).

      @param name the method name
      @param signature the method signature as defined by Java rules
   */
  JavaStaticIdAny(const char *name, const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    m_id = 0;
  }

  /**
      Ctor gets the id, throws an exception if failed.

      This constructor retrives the id directly from the class, it is more
      efficient than doing it from the object -- especially if you already
      have the class at hand.

      @param env pointer to the JNI environment
      @param cls the Java class whose method we're interested in
      @param name the member name
      @param signature the member signature as defined by Java rules
   */
  JavaStaticIdAny(Env& env,
       jclass cls,
       const char *name,
       const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    Init(env, cls);
  }

  /**
      Ctor gets the id, throws an exception if failed.

      This constructor first retrieves the class from the object and only
      then gets the id so it is less efficient than the constructor taking
      the class argument directly if you already have the class pointer
      somewhere. But if you don't, this constructor might be more convenient
      to use.

      @param env pointer to the JNI environment
      @param obj the Java object whose method we're interested in
      @param name the member name
      @param signature the member signature as defined by Java rules
   */
  JavaStaticIdAny(Env& env,
       jobject obj,
       const char *name,
       const char *signature)
    : m_name(name),
     m_signature(signature)
  {
    Init(env, obj);
  }

  /**
      Initializes the id if necessary and returns it.

      This method allows to use JavaStaticIdAny objects like this:
      @code
          static JavaStaticIdAny s_fid("name", "signature");
          ...
          env.SetStaticIntField(obj, s_fid(env, cls), ...);
      @endcode
   */
  T operator()(Env& env, jclass cls)
  {
    if ( !m_id )
      Init(env, cls);

    return m_id;
  }

  /**
      Initializes the id if necessary and returns it.

      This method allows to use JavaStaticIdAny objects like this:
      @code
          static JavaStaticIdAny s_fid("name", "signature");
          ...
          env.SetIntField(obj, s_fid(env, obj), ...);
      @endcode
   */
  T operator()(Env& env, jobject obj)
  {
    if ( !m_id )
      Init(env, obj);

    return m_id;
  }

  /**
      Provides direct access to the id.

      This supposes that the id had been already initialized, no checks are
      done to ensure that this is indeed the case in release build but in
      debug build an assert is raised if an uninitialized id is used.
   */
  operator T() const
  {
    ASSERT_MSG( m_id, "Using uninitialized JavaStaticIdAny" );

    return m_id;
  }

private:
  T m_id;
  const char * const m_name;
  const char * const m_signature;

  NO_COPY_CLASS(JavaStaticIdAny);
};

// ----------------------------------------------------------------------------
// FieldIdAny and FieldId<>
// ----------------------------------------------------------------------------

/**
    FieldIdAny is JavaIdAny specialized for the fields.
 */
typedef JavaIdAny<jfieldID> FieldIdAny;

/**
    FieldId template determines the Java field signature from its argument.

    This template may only be used for the primitive types, for the other
    ones you have to use FieldIdGeneric directly.

    Note that there is no default implementation of this class, only its
    specializations are defined.

    The template parameter is the Java type of the field.
 */
template <typename T>
class FieldId;

/// This macro is used below to define all FieldIdAny specializations.
#define JNI_DEFINE_FIELD_ID(Type)                                             \
  template <>                                                               \
  class FieldId<JNI_JTYPE(Type)> : public FieldIdAny                        \
  {                                                                         \
  public:                                                                   \
    FieldId(const char *name) : FieldIdAny(name, JNI_SIGNATURE(Type)) { } \
                                       \
    FieldId(Env& env, jclass cls, const char *name)             \
      : FieldIdAny(env, cls, name, JNI_SIGNATURE(Type)) { }             \
    \
  private: \
    NO_COPY_CLASS(FieldId); \
  };

/// And this one is used to define FieldId<> for the array types
#define JNI_DEFINE_FIELD_ID_ARRAY(Type)                                       \
  template <>                                                               \
  class FieldId<JNI_ATYPE(Type)> : public FieldIdAny                        \
  {                                                                         \
  public:                                                                   \
    FieldId(const char *name) : FieldIdAny(name, JNI_ASIGNATURE(Type)) { }\
                                       \
    FieldId(Env& env, jclass cls, const char *name)             \
      : FieldIdAny(env, cls, name, JNI_ASIGNATURE(Type)) { }            \
  };

/// Define FieldId specializations for all primitive types.
JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_FIELD_ID)

/// And for all arrays of primitive types
JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_FIELD_ID_ARRAY)

/// Finally also define FieldId<jstring>
JNI_DEFINE_FIELD_ID(String)

#undef JNI_DEFINE_FIELD_ID
#undef JNI_DEFINE_FIELD_ID_ARRAY


/**
    StaticFieldIdAny is JavaIdAny specialized for the fields.
 */
typedef JavaStaticIdAny<jfieldID> StaticFieldIdAny;

/**
    StaticFieldId template determines the Java field signature from its argument.

    This template may only be used for the primitive types, for the other
    ones you have to use FieldIdGeneric directly.

    Note that there is no default implementation of this class, only its
    specializations are defined.

    The template parameter is the Java type of the field.
 */
template <typename T>
class StaticFieldId;

/// This macro is used below to define all StaticFieldIdAny specializations.
#define JNI_DEFINE_FIELD_ID(Type)                                             \
  template <>                                                               \
  class StaticFieldId<JNI_JTYPE(Type)> : public StaticFieldIdAny                        \
  {                                                                         \
  public:                                                                   \
    StaticFieldId(const char *name) : StaticFieldIdAny(name, JNI_SIGNATURE(Type)) { } \
                                       \
    StaticFieldId(Env& env, jclass cls, const char *name)             \
      : StaticFieldIdAny(env, cls, name, JNI_SIGNATURE(Type)) { }             \
  private: \
    NO_COPY_CLASS(StaticFieldId); \
  };

/// And this one is used to define StaticFieldId<> for the array types
#define JNI_DEFINE_FIELD_ID_ARRAY(Type)                                       \
  template <>                                                               \
  class StaticFieldId<JNI_ATYPE(Type)> : public StaticFieldIdAny                        \
  {                                                                         \
  public:                                                                   \
    StaticFieldId(const char *name) : StaticFieldIdAny(name, JNI_ASIGNATURE(Type)) { }\
                                       \
    StaticFieldId(Env& env, jclass cls, const char *name)             \
      : StaticFieldIdAny(env, cls, name, JNI_ASIGNATURE(Type)) { }            \
  };

/// Define FieldId specializations for all primitive types.
JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_FIELD_ID)

/// And for all arrays of primitive types
JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_FIELD_ID_ARRAY)

/// Finally also define FieldId<jstring>
JNI_DEFINE_FIELD_ID(String)

#undef JNI_DEFINE_FIELD_ID
#undef JNI_DEFINE_FIELD_ID_ARRAY



// ----------------------------------------------------------------------------
// MethodIdAny and MethodId<>
// ----------------------------------------------------------------------------

/**
    This class represents a Java method id.

    It is similar to FieldIdAny but is used for methods instead of members. For
    the primitive types and arrays of them as well as jstring you should use
    the corresponding specialization of MethodId, this class should only be
    used for the other types as it doesn't provide type-safeness unlike
    MethodId.

    @sa FieldIdAny
 */
typedef JavaIdAny<jmethodID> MethodIdAny;

/// @todo: implement MethodId
template <typename T>
class MethodId;

/**
    This class represents a Java satic method id.

    It is similar to StaticFieldIdAny but is used for methods instead of members. For
    the primitive types and arrays of them as well as jstring you should use
    the corresponding specialization of MethodId, this class should only be
    used for the other types as it doesn't provide type-safeness unlike
    MethodId.

    @sa StaticFieldIdAny
 */
typedef JavaStaticIdAny<jmethodID> StaticMethodIdAny;

/// @todo: implement MethodId
template <typename T>
class StaticMethodId;

// ----------------------------------------------------------------------------
// JArray<>
// ----------------------------------------------------------------------------

/**
    Class for working with the arrays.

    This is a template class which may be used with the arrays of primitive
    types as well as with the arbitrary jobject arrays -- but without
    type-safety in the latter case.
 */
template <typename T>
class JArray;

/**
    Array of objects.

    This class can either create a new array or work with an existing one but
    in any case it never destroys the array, the idea being that if we create
    it, it's only to return it to Java side (and if we didn't create it, we
    surely shouldn't destroy it neither).
 */
template <>
class JArray<jobject>
{
public:
  /**
    Creates an array associated to an existing Java array.

    @param env the environment
    @param jarr the Java array
   */
  JArray(Env& env, jobjectArray jarr) : m_env(&env), m_array(jarr) { }

  /**
    Creates a new array of objects of given class.

    @param env the JNI environment
    @param length the array size
    @param cls the class of array elements
    @param init the value to initialize all the objects in the array with
   */
  JArray(Env& env, size_t length, const JClass& cls, jobject init = NULL)
    : m_env(&env)
    {
      m_array = env.NewObjectArray(static_cast<jsize>(length), cls, init);
    }

  /// Returns the number of elements in the array
  jsize GetLength() const
  {
    return m_env->GetArrayLength(m_array);
  }

  /// Gets the object at the given index
  jobject Get(size_t n) const
  {
    return m_env->GetObjectArrayElement(m_array, static_cast<jsize>(n));
  }

  /// Sets the object at the given index
  void Set(size_t n, jobject value)
  {
    m_env->SetObjectArrayElement(m_array, static_cast<jsize>(n), value);
  }

  /// Get the underlying Java array
  jobjectArray GetJArray() const { return m_array; }

private:
  /// The Java environment
  Env *m_env;

  /// The array itself
  jobjectArray m_array;
};

/**
    Base class for all JArray specializations.

    Only this class can access Env array methods because they're made
    private there to increase type-safeness, so the derived classes can access
    them only via the methods of this class.
 */
class JArrayAny
{
public:
  /**
      Default constructor leaves the array uninitialized, use Init() later.
   */
  JArrayAny() : m_env(NULL), m_iMode(JNI_ABORT) { }

  /**
      Constructor sets the mode initially to JNI_ABORT, call Commit() to
      change this.
   */
  JArrayAny(Env& env) : m_env(&env), m_iMode(JNI_ABORT) { }

  /**
      Initialize the object after using the default constructor.
   */
  void Init(Env& env) { m_env = &env; }

  /**
      Tell the array that the changes to its elements should be saved in
      Java.

      By default, any changes to the array elements are discarded and the
      Java array it corresponds to stays unchanged. However if Commit() is
      called, the changes will be copied back to Java but only after this
      object is destroyed and not immeidately.
   */
  void Commit() { m_iMode = 0; }

protected:
  /**
      Declare some methods for the derived JArray classes.
   */
#define JNI_FWD_GETRELEASE_ARRAY(Type)                                        \
  JNI_JTYPE(Type) *Get ## Type ## ArrayElements(JNI_ATYPE(Type) array)      \
  {                                                                         \
      return m_env->Get ## Type ## ArrayElements(array, NULL);          \
  }                                                                         \
                                       \
  void Release ## Type ## ArrayElements(JNI_ATYPE(Type) array,              \
                     JNI_JTYPE(Type) *elems)             \
  {                                                                         \
    m_env->Release ## Type ## ArrayElements(array, elems, m_iMode);       \
  }                                                                         \

  JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_FWD_GETRELEASE_ARRAY)

#undef JNI_FWD_GETRELEASE_ARRAY

  /// The real Java environment
  Env *m_env;

  /// Describes if we should propagate changes to array back to Java
  jint m_iMode;
};

/// Macro used to define JArray for all primitive types
#define JNI_DEFINE_ARRAY(Type)                                                \
  template <>                                                               \
  class JArray<JNI_JTYPE(Type)> : public JArrayAny                          \
  {                                                                         \
  public:                                                                   \
    JArray()                                                              \
    {                                                                     \
      m_data = NULL;                                                    \
    }                                                                     \
                                       \
    JArray(Env& env, JNI_ATYPE(Type) values)                              \
      : JArrayAny(env), m_values(values)                                \
    {                                                                     \
      m_data = Get ## Type ## ArrayElements(values);                    \
    }                                                                     \
                                       \
    void SetData()                                                        \
    {                                                                     \
      m_data = Get ## Type ## ArrayElements(m_values);                  \
    }                                                                     \
                                       \
    JArray(Env& env, size_t n_, JNI_JTYPE(Type) *values = NULL)            \
      : JArrayAny(env)                                                  \
    {                                                                     \
      const jsize n = static_cast<jsize>(n_);                           \
      m_values = env.New ## Type ## Array(n);                           \
      if ( values )                                                     \
      {                                                                 \
       env.Set ## Type ## ArrayRegion(m_values, 0, n, values);         \
      }                                                                 \
                                       \
      SetData();                                                        \
    }                                                                     \
                                       \
    void Init(Env& env, JNI_ATYPE(Type) values)                           \
    {                                                                     \
      JArrayAny::Init(env);                                             \
                                       \
      m_values = values;                                                \
                                       \
      SetData();                                                        \
    }                                                                     \
                                       \
    ~JArray()                                                             \
    {                                                                     \
      if ( m_data )                                                     \
        Release ## Type ## ArrayElements(m_values, m_data);           \
    }                                                                     \
                                       \
    jsize GetLength() const                                               \
    {                                                                     \
      return m_env ? m_env->GetArrayLength(m_values) : 0;               \
    }                                                                     \
                                       \
    operator const JNI_JTYPE(Type) *() const                              \
    {                                                                     \
      return m_data;                                                    \
    }                                                                     \
                                       \
    operator JNI_JTYPE(Type) *()                                          \
    {                                                                     \
      return m_data;                                                    \
    }                                                                     \
                                       \
    JNI_ATYPE(Type) GetJArray() const { return m_values; }                \
                                       \
  private:                                                                  \
    JNI_ATYPE(Type) m_values;                                             \
    JNI_JTYPE(Type) *m_data;                                              \
  };

JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_ARRAY)


/**
    This function creates a Java array from C data.

    A new Java array of the correct type (i.e. jdoubleArray for jdouble
    argument, ...) is created and initialized with the provided data.

    If the array creation fails, it throws.

    @param env pointer to the JNI environment
    @param nItems the number of items in pItems array
    @param pItems pointer to the array of (at least) nItems items
    @return a new Java array of the correct type
 */
#ifdef DOXYGEN

JNI_ATYPE(Type)
CreateFillArray(Env& env, size_t nItems, JNI_JTYPE(Type) *pItems);

#else // !DOXYGEN
#define JNI_DEFINE_CREATEFILLARRAY(Type)                                          \
  inline                                                                    \
  JNI_ATYPE(Type)                                                           \
  CreateFillArray(Env& env,                                                 \
          size_t nItems,                                            \
          JNI_JTYPE(Type) *pItems)                                  \
  {                                                                         \
    JNI_ATYPE(Type) jarray = env.New ## Type ## Array(static_cast<jsize>(nItems));            \
    env.Set ## Type ## ArrayRegion(jarray, 0, static_cast<jsize>(nItems), pItems);            \
    return jarray;                                                        \
  }
#endif // DOXYGEN/!DOXYGEN

JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_CREATEFILLARRAY)

#undef JNI_DEFINE_CREATEFILLARRAY

// ----------------------------------------------------------------------------
// JStringChars
// ----------------------------------------------------------------------------

/**
    This class provides easy and safe access to jstring contents.

    JStringChars implements RIAA policy for jstrings: it acquires the data in
    its constructor and frees it in its destructor, thus making sure that you
    don't forget to do it.
 */
class JStringChars
{
public:
  /**
      Constructor initializes data from the given string.

      May throw if an error occurs while getting data.

      @param env the checked JNIEnv object
      @param jstr the Java string to access
   */
  JStringChars(Env& env, jstring jstr)
    : m_env(env), m_jstr(jstr)
  {
    m_data = m_env.GetStringUTFChars(m_jstr, NULL);
  }

  /// Destructor frees the data.
  ~JStringChars()
  {
    m_env.ReleaseStringUTFChars(m_jstr, m_data);
  }

  /// Implicit conversion to an UTF-8 string.
  operator const char *() const { return m_data; }

private:
  Env& m_env;
  jstring m_jstr;
  const char *m_data;

  NO_ASSIGN_CLASS(JStringChars);
};

// ============================================================================
// inline functions implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Set/Get[Static]<Type>Field
// ----------------------------------------------------------------------------

/**
    This macro is used to define Get/Set[Static]TypeField for the given Type.
 */
#define JNI_DEFINE_SETGET_FIELD(Type)                                         \
  inline JNI_JTYPE(Type)                                                    \
  Env::Get ## Type ## Field(jobject obj,                          \
                    FieldId<JNI_JTYPE(Type)>& fid)        \
  {                                                                         \
    return Get ## Type ## Field(obj, fid(*this, obj));                    \
  }                                                                         \
                                       \
  inline JNI_JTYPE(Type)                                                    \
  Env::GetStatic ## Type ## Field(jclass cls,                     \
                       StaticFieldId<JNI_JTYPE(Type)>& fid)  \
  {                                                                         \
    return GetStatic ## Type ## Field(cls, fid(*this, cls));              \
  }                                                                         \
                                       \
  inline void                                                               \
  Env::Set ## Type ## Field(jobject obj,                          \
                    FieldId<JNI_JTYPE(Type)>& fid,        \
                    JNI_JTYPE(Type) value)                \
  {                                                                         \
    Set ## Type ## Field(obj, fid(*this, obj), value);                    \
  }                                                                         \
                                       \
  inline void                                                               \
  Env::SetStatic ## Type ## Field(jclass cls,                     \
                       StaticFieldId<JNI_JTYPE(Type)>& fid,  \
                       JNI_JTYPE(Type) value)          \
  {                                                                         \
    SetStatic ## Type ## Field(cls, fid(*this, cls), value);              \
  }

JNI_DO_FOR_ALL_PRIMITIVE_TYPES(JNI_DEFINE_SETGET_FIELD)

// ----------------------------------------------------------------------------
// Set/Get[Static]ObjectField
// ----------------------------------------------------------------------------

inline jobject
Env::GetObjectField(jobject obj, FieldIdAny& fid)
{
  return GetObjectField(obj, fid(*this, obj));
}

inline jobject
Env::GetStaticObjectField(jclass cls, StaticFieldIdAny& fid)
{
  return GetStaticObjectField(cls, fid(*this, cls));
}

inline JStringChars
Env::GetStringField(jobject obj, jfieldID fid)
{
  jstring jstr = reinterpret_cast<jstring>(GetObjectField(obj, fid));
  return JStringChars(*this, jstr);
}

// ----------------------------------------------------------------------------
// exception handling
// ----------------------------------------------------------------------------

/* static */ inline
void Env::JavaThrow(JNIEnv *env, const char *clsname, const std::string& msg)
{
  CHECK_VOID( env, "NULL JNIEnv in JavaThrow()" );

  try
  {
    Env envTmp(env);

    envTmp.m_env->ThrowNew(JClass(envTmp, clsname), msg.c_str());
  }
  catch ( ... )
  {
    // ignore
  }
}

} // namespace JNI

} // namespace ito33

#endif // _ITO33_JNIXX_H_

