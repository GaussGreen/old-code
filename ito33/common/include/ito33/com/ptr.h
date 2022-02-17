/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/ptr.h
// Purpose:     Ptr<> is a COM template auto pointer class
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: ptr.h,v 1.16 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/ptr.h
    @brief Ptr<> is std::auto_ptr<> for COM.

    Ptr<> is a smart pointer class using reference counting. Unlike auto_ptr<>,
    whose interface it mimics, it has normal copy semantics however.
 */

#ifndef _ITO33_COM_PTR_H_
#define _ITO33_COM_PTR_H_

#include <stdlib.h>
#include <winerror.h>

#include "ito33/com/utils.h"

namespace ito33
{

namespace COM
{

/**
    Ptr<> is a smart pointer for COM interfaces.
 */
template <class T>
class Ptr
{
public:
  /// a scalar type which doesn't risk to be converted to anything
  typedef T *(Ptr<T>::*unspecified_bool_type)() const;


  /**
      Ctor takes ownership of the pointer passed to it and will Release() it
      if it is not @c NULL.

      @param ptr the initial value of the pointer, NULL by default
   */
  explicit Ptr(T *ptr = NULL);

  /**
      This template ctor does a QueryInterface() on the given pointer for our
      interface.

      Don't confuse this with the non template ctor above which simpyl stores
      the pointer, this one is quite different!

      @param ptr any COM pointer or another Ptr<>
   */
  template <class IfacePtr> explicit Ptr(const IfacePtr& ptr)
  {
    // NB: this function (as all template methods) must be defined inline
    //     for VC++ to be able to compile it

    if ( ptr )
    {
      m_ptr = CastTo<T>(ptr.Get());
    }
    else // NULL in, NULL out
    {
      m_ptr = NULL;
    }
  }


  /**
      Normal copy ctor: unlike auto_ptr<>, it doesn't modify its argument.

      The reference count of the other pointer is incremented here (unless it
      is @c NULL)
   */
  Ptr(const Ptr<T>& other);

  /**
      Normal assignement operator.

      The previously held pointer is Release()d and the reference count of
      the new pointer is incremented here.
   */
  Ptr<T>& operator=(const Ptr<T>& other);

  /**
      Assignment operator from a raw pointer.

      The previously held pointer is Release()d but the reference count of
      the new pointer is not incremented here -- in other words, we take
      ownership of the pointer. Call Assign() if this is not what you want.
   */
  Ptr<T>& operator=(T *pOther);

  /**
      Dtor calls Release() on the pointer if it is not @c NULL.
   */
  ~Ptr();

  /**
      Assigns a raw pointer to this object.

      We will hold onto the pointer by calling an AddRef() on it instead of
      simply taking its ownership as the assignment operator does.
   */
  void Assign(T *ptr);

  /**
      Returns the stored pointer or @c NULL if none.

      Don't call AddRef()/Release() via this pointer as it may lead to memory
      leaks/crashes. Better yet, don't use this directly at all.
   */
  T *Get() const;

  /**
      Release the held pointer, if any.

      Decrements the reference count of our pointer, if non @c NULL, and
      returns the old pointer value. Note that the pointer may be invalid by
      then (if it had reference count of 1 before calling this method) so use
      the returned value at your own risk.
   */
  T *Release();

  /**
      Member access operator.
   */
  T *operator->() const;

  /**
      Implicit, but safe, conversion to bool.

      Having this conversion ensures that we can work with the objects of
      this type as with the plain pointers but using a unique type (instead
      of bool or int) ensures that we don't risk to implicitly confuse the
      unrelated pointers.
   */
  operator unspecified_bool_type() const
  {
    return m_ptr ? &Ptr<T>::Get : NULL;
  }

private:
  // hold to the pointer by calling AddRef() on it
  void Hold(T *ptr);

  // release the pointer if it is not NULL
  void Free();

  T *m_ptr;
};

/**
    Returns the pointer cast to the correct type.

    This is a convenience function which is easier to use than doing
    static_cast<U>(ptr->Get()). Its name comes from the fact that it is most
    often used to cast an IFoo pointer to FooImpl one -- but it's not really
    restricted to just this usage.

    @param p the smart pointer to cast (or anything else with get method)
    @return the resulting pointer of the correct type
 */
template <class U, class T>
inline U *CastToImpl(T p) { return static_cast<U *>(p.Get()); }

/// returns a pointer AddRef()ing it
template <class Iface>
inline HRESULT
ReturnValue(const Ptr<Iface>& ptr, Iface **ppOut)
{
    CHECK( ppOut, E_POINTER, "NULL output pointer" );

    if ( (*ppOut = ptr.Get()) != NULL )
      (*ppOut)->AddRef();

    return S_OK;
}

// ----------------------------------------------------------------------------
// implementation only from now on
// ----------------------------------------------------------------------------

template <class T>
inline Ptr<T>::Ptr(T *ptr)
      : m_ptr(ptr)
{
}

template <class T>
inline void Ptr<T>::Free()
{
  if ( m_ptr )
    m_ptr->Release();
}

template <class T>
inline void Ptr<T>::Hold(T *ptr)
{
  if ( (m_ptr = ptr) != NULL )
    m_ptr->AddRef();
}

template <class T>
inline Ptr<T>::Ptr(const Ptr<T>& other)
{
  Hold(other.m_ptr);
}

template <class T>
inline Ptr<T>& Ptr<T>::operator=(T *pOther)
{
  Free();

  // do not AddRef() it here, we take ownership of the pointer
  m_ptr = pOther;

  return *this;
}

template <class T>
inline Ptr<T>& Ptr<T>::operator=(const Ptr<T>& pOther)
{
  Free();

  Hold(pOther.Get());

  return *this;
}

template <class T>
inline Ptr<T>::~Ptr()
{
  Free();
}

template <class T>
inline void Ptr<T>::Assign(T *ptr)
{
  Free();

  Hold(ptr);
}

template <class T>
inline T *Ptr<T>::Get() const
{
  return m_ptr;
}

template <class T>
inline T *Ptr<T>::Release()
{
  Free();

  T *ptrOld = m_ptr;
  m_ptr = NULL;

  return ptrOld;
}

template <class T>
inline T *Ptr<T>::operator->() const
{
  return Get();
}

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_PTR_H_

