/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/sharedptr.h
// Purpose:     a ref-counted smart pointer
// Author:      Vadim Zeitlin
// Created:     16.04.03
// RCS-ID:      $Id: sharedptr.h,v 1.51 2006/08/19 18:40:28 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/sharedptr.h
    @brief Provides a ref-counted smart pointer class, and brings 
           boost::shared_ptr into our namespace.

    Ref counted smart pointers provide normal value semantics unlike the
    standard auto_ptr<> (or our AutoPtr<>).

    The name SharedPtr (but not the implementation) were taken from the
    analogous boost smart pointer class.
 */

#ifndef _ITO33_SHAREDPTR_H_
#define _ITO33_SHAREDPTR_H_

#include "ito33/beforestd.h"
#include "ito33/afterstd.h"

#include "ito33/common.h"

#ifndef __CPP2ANY__

#include <boost/shared_ptr.hpp>

#include "ito33/thread.h"

namespace ito33
{

/// using directive to bring boost shared_ptr related class/function into ito33
using boost::shared_ptr;
using boost::static_pointer_cast;
using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;

/**
    This function explicitly converts a raw pointer to shared_ptr.

    The advantage of using it instead of just writing shared_ptr<T>(new T(...))
    is simply that you don't have to repeat (potentially long) name of T.

    The disavantage is that shared_ptr will be instantiated for T even if
    you need only a shared pointer to a base class.

    Normally, the value should be assigned to a variable previously declared
    since using unnamed shared_ptr temporaries might not be exception-safe.
 */
template <typename T>
inline shared_ptr<T> make_ptr(T *ptr)
{
  return shared_ptr<T>(ptr);
}

// forward declarations
template <typename T, class M> class WeakPtr;
namespace Private
{
  template <class M> class SharedPtrCount;
}


/**
    This class is used by SharedPtr to delete the pointer it holds.

    This is a template class but it doesn't have any generic implementation,
    you must specialize it for each type with which you use SharedPtr and
    implement the Free() method for it.
 */
template <class T>
struct SharedPtrDeleter
{
  /// this is called by ~SharedPtr
  static void Free(T *p);
};

/**
    Namespace containing different multithreading policies, i.e. classes
    which are used as template parameters and determine the MT-safe (or not)
    aspect of behaviour of another class.
 */
namespace Policy
{

/**
    MT unsafe policy: no overhead but can't be used in presence of multiple
    threads.
 */
struct MTUnsafe
{
  static ::ito33::Private::SharedPtrCount<MTUnsafe> *CreateSharedCount();
  static void Destroy(::ito33::Private::SharedPtrCount<MTUnsafe> *p);
};

/**
    MT-safe policy using atomic operations: smallest overhead necessary for
    MT-safe operation.
 */
struct MTSafe
{
  static ::ito33::Private::SharedPtrCount<MTSafe> *CreateSharedCount();
  static void Destroy(::ito33::Private::SharedPtrCount<MTSafe> *p);
};

} // namespace Policy


template <class T, class M, class D> class EnableSharedFromThis;

// private stuff, don't document
namespace Private
{

/**
    Object containing the reference count for SharedPtr and WeakPtr.

    All SharedPtr objects using the same pointer store one and the same
    SharedPtrCount object. It is destroyed only when there are no more
    SharedPtr nor WeakPtr using it.

    This class is private, you should never use it directly.
 */
class SharedPtrCountBase
{
public:
  /// Ctor initializes the reference count to 1, making the object alive.
  SharedPtrCountBase() : m_nTotalRef(1), m_nStrongRef(1) { }

protected:
  // total reference count: strong + weak, when it reaches 0 this object itself
  // is deleted
  int m_nTotalRef;

  // strong reference count: when it reaches 0, the pointer is deleted
  int m_nStrongRef;
};

template <class M>
class SharedPtrCount;

// MT-unsafe version of SharedPtrCount
template <>
class SharedPtrCount<Policy::MTUnsafe> : private SharedPtrCountBase
{
public:
  typedef Policy::MTUnsafe MTPolicy;

  SharedPtrCount() : SharedPtrCountBase() { }

  /// increment SharedPtr reference count
  void AddStrongRef() { m_nStrongRef++; m_nTotalRef++; }

  /// increment SharedPtr reference count if the pointer is still alive; do
  /// nothing and return false otherwise
  bool AddStrongRefIfAlive()
  {
    return m_nStrongRef == 0 ? false : (AddStrongRef(), true);
  }

  /// increment WeakPtr reference count
  void AddWeakRef() { m_nTotalRef++; }

  /// decrement strong ref count and returns false if it drops to 0
  ///
  /// also deletes the counter itself if total ref count drops to 0 as well
  bool DecStrongRef()
  {
    const bool keepAlivePtr = --m_nStrongRef != 0;
    if ( !--m_nTotalRef )
      MTPolicy::Destroy(this);

    return keepAlivePtr;
  }

  /// same as DecStrongRef() but only if there is exactly one strong ref left,
  /// otherwise doesn't do anything
  ///
  /// returns false if we were destroyed, true if there is more than one ref
  bool DecStrongRefIfLast()
  {
    return m_nStrongRef == 1 ? DecStrongRef() : true;
  }

  /// decrement SharedPtr reference count and deletes this object if it drops
  /// to 0
  void DecWeakRef()
  {
    if ( !--m_nTotalRef )
      MTPolicy::Destroy(this);
  }
};

// MT-safe version of SharedPtrCount
template <>
class SharedPtrCount<Policy::MTSafe> : private SharedPtrCountBase
{
public:
  typedef Policy::MTSafe MTPolicy;

  SharedPtrCount() : SharedPtrCountBase() { }

  void AddStrongRef()
  {
    Lock<CriticalSection> lock(m_cs);

    m_nStrongRef++;
    m_nTotalRef++;
  }

  bool AddStrongRefIfAlive()
  {
    Lock<CriticalSection> lock(m_cs);

    if ( !m_nStrongRef )
      return false;

    m_nStrongRef++;
    m_nTotalRef++;

    return true;
  }

  void AddWeakRef()
  {
    Lock<CriticalSection> lock(m_cs);

    m_nTotalRef++;
  }

  bool DecStrongRef()
  {
    bool keepAlivePtr,
         destroySelf;

    // we have to unlock m_cs before calling Destroy() which destroys it
    {
      Lock<CriticalSection> lock(m_cs);
      keepAlivePtr = --m_nStrongRef != 0;
      destroySelf = --m_nTotalRef == 0;
    }

    if ( destroySelf )
      MTPolicy::Destroy(this);

    return keepAlivePtr;
  }

  bool DecStrongRefIfLast()
  {
    bool destroySelf;

    {
      Lock<CriticalSection> lock(m_cs);
      if ( m_nStrongRef != 1 )
        return true;

      --m_nStrongRef;
      destroySelf = --m_nTotalRef == 0;
    }

    if ( destroySelf )
      MTPolicy::Destroy(this);

    return false;
  }

  void DecWeakRef()
  {
    // as above, take care to unlock m_cs before calling Destroy()
    bool last;
    {
      Lock<CriticalSection> lock(m_cs);
      last = --m_nTotalRef == 0;
    }

    if ( last )
      MTPolicy::Destroy(this);
  }

private:
  // critical section to protect counters from concurrent modification
  mutable CriticalSection m_cs;

  NO_COPY_TEMPLATE_CLASS(SharedPtrCount, Policy::MTSafe);
};

// parameters for SharedPtr ctors, just to distinguish them
namespace Cast
{
    struct Const;
    struct Static;
    struct Dynamic;
} // namespace Cast

/**
    This class is a friend of both WeakPtr and EnableSharedFromThis and so
    can use their private methods which are only used for EnableSharedFromThis
    initialization.
 */
struct SelfPtrSetter
{
  template <class T, class M, class D>
  static void Assign(EnableSharedFromThis<T, M, D>& p,
                     ::ito33::Private::SharedPtrCount<M>& count)
  {
    p.m_self.InitFromExistingCount(static_cast<T *>(&p), count);
  }
};

/**
    Helper for EnableSharedFromThis: this function does template dispatching
    depending on whether the class T derives from EnableSharedFromThis or not.
 */
template <class T, class M, class D>
inline void
DoEnableSharedFromThis(::ito33::Private::SharedPtrCount<M>& count,
                       EnableSharedFromThis<T, M, D> *p)
{
  if ( p )
    SelfPtrSetter::Assign(*p, count);
}

template <class M>
inline void
DoEnableSharedFromThis(::ito33::Private::SharedPtrCount<M>& /* count */, ...)
{
  // this version does nothing
}

} // namespace Private


/**
    SharedPtr implements smart pointer using reference counting.

    None of the methods of this class throw exceptions, except if T dtor throws
    (it shouldn't!)

    Template parameters are:
        - T is the type of the pointee
        - D is the type used for deleting the objects of type T and is by
            default SharedPtrDeleter which simply calls delete on the pointer
        - M is the multithreading policy
 */
template <typename T,
          class M = ::ito33::Policy::MTUnsafe,
          class D = ::ito33::SharedPtrDeleter<T> >
class SharedPtr
{
private:
  /// hold to the given pointer, if any
  void Hold(T *ptr)
  {
    if ( m_pCount )
      m_pCount->AddStrongRef();

    m_ptr = ptr;
  }

  /// give up ownership of the pointer we currently have, if any
  void Release()
  {
    if ( m_pCount )
    {
      if ( !m_pCount->DecStrongRef() )
      {
        // last reference, destroy the object
        D::Free(m_ptr);
      }
    }
  }

  /// reset both pointers to NULL (without freeing anything)
  void Reinit()
  {
    m_pCount = NULL;
    m_ptr = NULL;
  }

public:
  /// the type of the object we're pointing to
  typedef T element_type;

  /// synonym for element_type
  typedef T value_type;

  /// synonym for pointer to element_type
  typedef T *pointer;

  /// a scalar type which doesn't risk to be converted to anything
  typedef T *(SharedPtr<T, M, D>::*unspecified_bool_type)() const;

  /// our MT policy
  typedef M MTPolicy;


  /// default ctor, we don't have any valid pointer
  SharedPtr() { Reinit(); }

  /// takes ownership of the pointer (which shouldn't normally be @c NULL)
  explicit SharedPtr(T *ptr) : m_pCount(MTPolicy::CreateSharedCount())
  {
    m_ptr = ptr;

    ::ito33::Private::DoEnableSharedFromThis(*m_pCount, m_ptr);
  }

  /// template copy ctor from a compatible shared pointer
  template <typename T2>
  SharedPtr(const SharedPtr<T2, M>& other) : m_pCount(other.GetCount())
  {
    Hold(other.get());
  }

  /// template copy ctor from a compatible weak pointer
  template <typename T2>
  SharedPtr(const WeakPtr<T2, M>& other)
  {
    m_pCount = other.GetCount();
    if ( m_pCount && m_pCount->AddStrongRefIfAlive() )
    {
        m_ptr = other.get();
    }
    else
    {
        Reinit();
    }
  }

  /**
     copy ctor using const_cast<>

     @sa ConstCast()
    */
  template <typename T2>
      SharedPtr(const SharedPtr<T2, M>& other, ::ito33::Private::Cast::Const *)
    : m_pCount(other.GetCount())
  {
    Hold(const_cast<element_type *>(other.get()));
  }

  /**
    copy ctor using static_cast<>

    @sa StaticCast()
    */
  template <typename T2>
  SharedPtr(const SharedPtr<T2, M>& other, ::ito33::Private::Cast::Static *)
    : m_pCount(other.GetCount())
  {
    Hold(static_cast<element_type *>(other.get()));
  }

  /**
    copy ctor using dynamic_cast<>

    @sa DynamicCast()
    */
  template <typename T2>
  SharedPtr(const SharedPtr<T2, M>& other, ::ito33::Private::Cast::Dynamic *)
    : m_pCount(other.GetCount())
  {
    Hold(dynamic_cast<element_type *>(other.get()));
  }

  // VC++ 6 can't compile the template copy ctor/assignment operator
#if !defined(_MSC_VER) || (_MSC_VER >= 1300)
  /// normal copy ctor
  SharedPtr(const SharedPtr<T, M>& other) : m_pCount(other.m_pCount)
  {
    Hold(other.m_ptr);
  }

  /// assignment operator
  SharedPtr<T, M>& operator=(const SharedPtr<T, M>& other)
  {
    if ( this != &other )
    {
      Release();

      m_pCount = other.m_pCount;

      Hold(other.m_ptr);
    }

    return *this;
  }
#endif // VC++ < 7

  /// assignment operator from a raw pointer
  SharedPtr<T, M>& operator=(T *ptr)
  {
    Release();

    m_pCount = MTPolicy::CreateSharedCount();
    m_ptr = ptr;

    ::ito33::Private::DoEnableSharedFromThis(*m_pCount, m_ptr);

    return *this;
  }

  /// reset the pointer to @c NULL value
  void reset()
  {
    *this = NULL;
  }

  /**
      Release the ownership of the pointer if we are the sole owner of it.

      This can be used to transfer ownership of the object to another
      SharedPtr<> (with possibly different MTPolicy or deleter) or to auto_ptr
      or AutoPtr.

      @return the current pointer or @c NULL if we are not the only owner
   */
  T *release()
  {
    if ( m_pCount && !m_pCount->DecStrongRefIfLast() )
    {
      T *ptr = m_ptr;
      Reinit();
      return ptr;
    }

    return NULL;
  }

  /// dtor decrements the ref count
  ~SharedPtr() { Release(); }

  /// returns the stored pointer, can be @c NULL
  T *get() const { return m_ptr; }

  /// allows use of this class as a pointer, must be non @c NULL
  T& operator*() const { return *get(); }

  /// allows use of this class as a pointer, must be non @c NULL
  T *operator->() const { return get(); }

  /**
      Implicit, but safe, conversion to bool.

      Having this conversion ensures that we can work with the objects of
      this type as with the plain pointers but using a unique type (instead
      of bool or int) ensures that we don't risk to implicitly confuse the
      unrelated pointers.
   */
  operator unspecified_bool_type() const // never throws
  {
    return m_ptr ? &SharedPtr<T, M, D>::get : NULL;
  }

  /// comparison
  bool operator==(const SharedPtr<T, M>& other) const
  {
    return m_ptr == other.m_ptr;
  }

  /// comparison
  bool operator!=(const SharedPtr<T, M>& other) const
  {
    return !(*this == other);
  }

  // operator< is needed to be able to use SharedPtrs as std::map keys
  bool operator<(const SharedPtr<T, M>& other) const
  {
    return m_ptr < other.m_ptr;
  }

  /**
      @name Pointer casts.

      As with normal pointers, casting the shared pointers shouldn't be
      normally done but, again, just like with the normal ones, sometimes it
      is necessary so we provide functions to do it safely (DynamicCast()) and
      unconditionally (StaticCast()).

      Unfortunately, VC++ 6 doesn't allow explicitly selecting the template
      function to call (i.e. "foo<T>()" doesn't compile), so we have to add
      dummy extra arguments for all the functions below. To hide this ugliness
      from the user, we also define PTR_CONST/STATIC/DYNAMIC_CAST() macros
      below which you should use instead of calling these functions directly.
   */
  //@{

  /// Cast away cv-qualifier using const_cast<>
  template <typename U>
  SharedPtr<U, M> ConstCast(U * = NULL) const
  {
    return SharedPtr<U, M>(*this, (::ito33::Private::Cast::Const *)NULL);
  }

  /// Convert an existing pointer to another type using static_cast<>
  template <typename U>
  SharedPtr<U, M> StaticCast(U * = NULL) const
  {
    return SharedPtr<U, M>(*this, (::ito33::Private::Cast::Static *)NULL);
  }

  /// Convert an existing pointer to another type using dynamic_cast<>
  template <typename U>
  SharedPtr<U, M> DynamicCast(U * = NULL) const
  {
    return SharedPtr<U, M>(*this, (::ito33::Private::Cast::Dynamic *)NULL);
  }

  //@}


  /// Implementation only, don't use
  ::ito33::Private::SharedPtrCount<MTPolicy> *GetCount() const
  {
    return const_cast<SharedPtr *>(this)->m_pCount;
  }

private:
  /// pointer to the object containing the (shared) reference count
  ::ito33::Private::SharedPtrCount<MTPolicy> *m_pCount;

  /// pointer we keep itself
  element_type *m_ptr;
};

/**
    @name Pointer casts.

    These macros are unfortunately needed to work around VC 6 brokeness.

    @sa ConstCast, StaticCast, DynamicCast
 */
//@{

/// const_cast<> for shared pointers
#define PTR_CONST_CAST(T, p)  (p).ConstCast(static_cast<T *>(NULL))

/// static_cast<> for shared pointers
#define PTR_STATIC_CAST(T, p)  (p).StaticCast(static_cast<T *>(NULL))

/// dynamic_cast<> for shared pointers
#define PTR_DYNAMIC_CAST(T, p)  (p).DynamicCast(static_cast<T *>(NULL))

//@}

/**
    A macro which expands into a normal implementation of SharedPtrDeleter.

    This is similar to ITO33_IMPLEMENT_AUTOPTR.
 */
#define ITO33_IMPLEMENT_SHAREDPTR(T) \
  template <> \
  void SharedPtrDeleter< T >::Free(T *p) { delete p; } \
  /* just to force semicolon after the macro */ struct Dummy

/**
    This function explicitly converts a raw pointer to SharedPtr.

    The advantage of using it instead of just writing SharedPtr<T>(new T(...))
    is simply that you don't have to repeat (potentially long) name of T.

    It works with the standard (MT-unsafe) policy and standard deleter only!
 */
template <typename T>
inline SharedPtr<T> MakePtr(T *ptr)
{
  return SharedPtr<T>(ptr);
}

/**
    Same as MakePtr but creates a MT-safe pointer.
 */
template <typename T>
inline SharedPtr<T, ::ito33::Policy::MTSafe> MakeMTPtr(T *ptr)
{
  return SharedPtr<T, ::ito33::Policy::MTSafe>(ptr);
}


/**
    Deriving from this class allows to obtain a WeakPtr to this object which
    is often useful when the object has to register itself with some other
    class which needs a SharedPtr.

    The name is the same of the equivalent boost class, modulo naming
    conventions.

    Example of usage:
      @code
        class X
        {
        public:
          void Register(const SharedPtr<Foo>& foo) { ... }
        };

        class Foo : public EnableSharedFromThis<Foo>
        {
        public:
          Foo() { }

          void RegisterWithX(X& x)
          {
            x.Register(SelfPtr());
          }
        };
      @endcode
 */
template <class T,
          class M = ::ito33::Policy::MTUnsafe,
          class D = ::ito33::SharedPtrDeleter<T> >
class EnableSharedFromThis
{
public:
  /// the type of the object we're pointing to
  typedef T element_type;

  /// our MT policy
  typedef M MTPolicy;

  /// type of weak pointer to this class itself
  typedef WeakPtr<element_type, MTPolicy> WeakPtrToSelf;

  /// type of strong pointer to this class itself
  typedef SharedPtr<element_type, MTPolicy, D> PtrToSelf;


  /**
      Return a weak pointer to this object itself.
   */
  WeakPtrToSelf Self() const { return m_self; }

  /**
      Return a strong pointer to this object itself.

      Normally this pointer is always valid but it may be not if the object is
      being destroyed.
   */
  PtrToSelf SelfPtr() const { return PtrToSelf(m_self); }

protected:
  // ctor is protected, this class is never created directly
  EnableSharedFromThis() { }

  // default copy ctor, assignment operator and dtor are ok


  WeakPtrToSelf m_self;

  // give it access to our m_self so that it could initialize it
  friend struct ::ito33::Private::SelfPtrSetter;
};

/**
    Same as EnableSharedFromThis except it uses MTSafe policy.
 */
template <class T, class D = ::ito33::SharedPtrDeleter<T> >
class EnableMTSharedFromThis : public EnableSharedFromThis<T, Policy::MTSafe, D>
{
protected:
  EnableMTSharedFromThis() { }
};

} // namespace ito33

#endif // __CPP2ANY__

#endif // _ITO33_SHAREDPTR_H_
