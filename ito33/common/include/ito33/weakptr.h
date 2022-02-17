/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/weakptr.h
// Purpose:     a ref-counted weak pointer
// Author:      Vaclav Slavik
// Created:     2005-01-21
// RCS-ID:      $Id: weakptr.h,v 1.8 2006/02/08 17:28:38 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/weakptr.h
    @brief Provides a weak pointer class for monitoring SharedPtr.

    Ref counted smart pointers provide normal value semantics unlike the
    standard auto_ptr<> (or our AutoPtr<>).

    The name (but not the implementation) were taken from the analogous boost
    weak pointer class.
 */

#ifndef _ITO33_WEAKPTR_H_
#define _ITO33_WEAKPTR_H_

#include "ito33/sharedptr.h"

namespace ito33
{

/**
    WeakPtr stores a weak reference to object managed by SharedPtr. Weak
    reference doesn't increase object's reference count, so it can be destroyed
    even when there are WeakPtrs pointing to it. SharedPtr created from WeakPtr
    dereferences to NULL if the object had been already destroyed.

    Weak pointer cannot be dereferenced directly, because it's unsafe
    operation: unlike SharedPtr, having a WeakPtr pointing to an object doesn't
    ensure that the object won't be deleted. If you need to access the object,
    convert the pointer to SharedPtr first:

      @code
      SharedPtr<Foo> ptr(weakptr);
      if ( ptr ) // ptr is NULL if weakptr expired
      {
        ...
      }
      @endcode

    None of the methods of this class throw exceptions.

    Template parameters are:
        - T is the type of the pointee
        - M is the multithread safety policy, MTUnsafe by default
 */
template <class T, class M = ::ito33::Policy::MTUnsafe>
class WeakPtr
{
private:
  typedef Private::SharedPtrCount<M> Count;


  /// hold to the given pointer, if any
  void Hold(T *ptr)
  {
    if ( m_pCount )
      m_pCount->AddWeakRef();

    m_ptr = ptr;
  }

  /// this function is only used by DoEnableSharedFromThis() via SelfPtrSetter
  void InitFromExistingCount(T *ptr, Count& count)
  {
    m_pCount = &count;
    Hold(ptr);
  }

  friend struct Private::SelfPtrSetter;


  /// give up ownership of the pointer we currently have, if any
  void Release()
  {
    if ( m_pCount )
      m_pCount->DecWeakRef();
  }

public:
  /// the type of the object we're pointing to
  typedef T element_type;

  /// synonym for element_type
  typedef T value_type;

  /// synonym for pointer to element_type
  typedef T *pointer;

  /// our MT policy
  typedef M MTPolicy;


  /// default ctor, we don't have any valid pointer
  WeakPtr() { m_pCount = NULL; m_ptr = NULL; }
  
  /// template copy ctor from a compatible shared pointer
  template <typename U>
  WeakPtr(const SharedPtr<U, M>& other) : m_pCount(other.GetCount())
  {
    Hold(other.get());
  }

  /// template copy ctor from a compatible weak pointer
  template <typename U>
  WeakPtr(const WeakPtr<U, M>& other) : m_pCount(other.GetCount())
  {
    Hold(other.get());
  }

  /// copy ctor using const_cast<>, @sa ConstCast()
  template <typename U>
  WeakPtr(const WeakPtr<U, M>& other, Private::Cast::Const *)
    : m_pCount(other.GetCount())
  {
    Hold(const_cast<element_type *>(other.get()));
  }

  /// copy ctor using static_cast<>, @sa StaticCast()
  template <typename U>
  WeakPtr(const WeakPtr<U, M>& other, Private::Cast::Static *)
    : m_pCount(other.GetCount())
  {
    Hold(static_cast<element_type *>(other.get()));
  }

  /// copy ctor using dynamic_cast<>, @sa DynamicCast()
  template <typename U>
  WeakPtr(const WeakPtr<U, M>& other, Private::Cast::Dynamic *)
    : m_pCount(other.GetCount())
  {
    Hold(dynamic_cast<element_type *>(other.get()));
  }

  /// normal copy ctor
  WeakPtr(const WeakPtr<T, M>& other) : m_pCount(other.m_pCount)
  {
    Hold(other.m_ptr);
  }

  /// assignment operator
  WeakPtr<T, M>& operator=(const WeakPtr<T, M>& other)
  {
    if ( this != &other )
    {
      Release();

      m_pCount = other.m_pCount;

      Hold(other.m_ptr);
    }

    return *this;
  }

  /// assignment operator from a shared pointer
  WeakPtr<T, M>& operator=(const SharedPtr<T, M>& other)
  {
    Release();

    m_pCount = other.GetCount();

    Hold(other.get());

    return *this;
  }

  /// get a (possible @c NULL) strong pointer corresponding to this one
  ///
  /// The name of this function is the same as in boost::weak_ptr (or the new
  /// std::tr1::weak_ptr) which explains why we use it although it's hardly the
  /// most logical one for this operation. In any case, it is rarely used as it
  /// is usually simpler to directly create a SharedPtr from a WeakPtr and is
  /// needed only for some template constructs (see weekalgo.h) where otherwise
  /// the compiler might be unable to perform the type deduction.
  SharedPtr<T, M> lock() const
  {
    return SharedPtr<T, M>(*this);
  }

  /// reset the pointer to @c NULL value
  void reset()
  {
    Release();
    m_pCount = NULL;
    m_ptr = NULL;
  }

  /// dtor decrements the ref count
  ~WeakPtr() { Release(); }

  /// comparison
  bool operator==(const WeakPtr<T, M>& other) const
  {
    // NB: We don't compare the pointers themselves because if there are no
    //     strong references left to the pointer, m_ptr usually points to
    //     unallocated memory, but it could also point to newly allocated
    //     _another_ object and so two different WeakPtrs could theoretically
    //     have same m_ptr. OTOH, m_pCount cannot be same for two different
    //     WeakPtr instances.
    return m_pCount == other.m_pCount;
  }

  /// comparison
  bool operator!=(const WeakPtr<T, M>& other) const
  {
    return !(*this == other);
  }

  // operator< is needed to be able to use WeakPtrs as std::map keys
  bool operator<(const WeakPtr<T, M>& other) const
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
  WeakPtr<U, M> ConstCast(U * = NULL) const
  {
    return WeakPtr<U, M>(*this, (Private::Cast::Const *)NULL);
  }

  /// Convert an existing pointer to another type using static_cast<>
  template <typename U>
  WeakPtr<U, M> StaticCast(U * = NULL) const
  {
    return WeakPtr<U, M>(*this, (Private::Cast::Static *)NULL);
  }

  /// Convert an existing pointer to another type using dynamic_cast<>
  template <typename U>
  WeakPtr<U, M> DynamicCast(U * = NULL) const
  {
    return WeakPtr<U, M>(*this, (Private::Cast::Dynamic *)NULL);
  }

  //@}


private:
  /// Implementation only, don't use
  Count *GetCount() const
  {
    return const_cast<WeakPtr *>(this)->m_pCount;
  }

  /// returns the stored pointer, can be @c NULL (we are keeping it private
  /// on purpose, it's only meant for use by WeakPtr and SharedPtr)
  T *get() const { return m_ptr; }
  
  /// pointer to the object containing the (shared) reference count
  Count *m_pCount;

  /// pointer we keep itself
  element_type *m_ptr;

  template <typename U, class V> friend class WeakPtr;
  template <typename U, class V, class D> friend class SharedPtr;
};

} // namespace ito33

#endif // _ITO33_WEAKPTR_H_
