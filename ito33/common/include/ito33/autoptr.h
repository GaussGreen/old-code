/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/autoptr.h
// Purpose:     replacement for the standard auto_ptr<> class
// Author:      Vadim Zeitlin
// Created:     16.04.03
// RCS-ID:      $Id: autoptr.h,v 1.23 2006/06/07 19:23:25 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/autoptr.h
    @brief This header defines a replacement for std::auto_ptr<> template class.

    The standard auto_ptr<> class suffers from many problems but even assuming
    that you're ready to live with its weird semantics there are still two
    things which can make your life miserable:
        - as ~auto_ptr is inline, it will call the delete defined in the user
          code which could be inappropriate for freeing the objects allocated
          with the library new (specifically in the case of VC++ which has many
          different and incompatible versions of the standard library)
        - it is impossible to use auto_ptr with incomplete types

    This class solves both of these problems but at a price of a small
    performance loss (extra function call coming from the fact that its dtor is
    not inlined any more) and some increased complexity of use. Additionally,
    it could be more standard-compliant than the auto_ptr<> class provided by
    the standard library, in particular VC++ 6 auto_ptr<> doesn't have the
    standard reset() method which we do.
 */

#ifndef _ITO33_AUTOPTR_H_
#define _ITO33_AUTOPTR_H_

#include <stddef.h>         // for NULL

#include <memory>

#include "ito33/dlldecl.h"

namespace ito33
{

/**
    This class is used by AutoPtr to delete the pointer it holds.

    This is a template class but it doesn't have any generic implementation,
    you must specialize it for each type with which you use AutoPtr and
    implement the Free() method for it.
 */
template <class T>
struct AutoPtrDeleter
{
  /// this is called by ~AutoPtr
  static void Free(T *p);
};

/**
    A macro which expands into a normal implementation of AutoPtrDeleter.

    Although AutoPtrDeleter may be implemented in many different ways, in 90%
    of cases Free() just deletes the pointer passed to it. This macro allows to
    define the specialization of AutoPtrDeleter which does this with minimal
    typing.

    This macro must always be used inside namespace ito33!
 */
#define ITO33_IMPLEMENT_AUTOPTR(T) \
  template <> \
  void AutoPtrDeleter< T >::Free(T *p) { delete p; } \
  /* just to force semicolon after the macro */ struct Dummy

#ifdef _WIN32
  #define ITO33_IMPLEMENT_DLLDECL_AUTOPTR(T) \
    template <> \
    struct ITO33_DLLDECL AutoPtrDeleter<T> { \
    static void Free(T *p) { delete p; } }; \
    template class ITO33_DLLDECL AutoPtr<T>; \
    /* just to force semicolon after the macro */ struct Dummy
#else //!_WIN32, just the same as ITO33_IMPLEMENT_AUTOPTR 
  #define ITO33_IMPLEMENT_DLLDECL_AUTOPTR(T) \
  ITO33_IMPLEMENT_AUTOPTR(T)
#endif 

/**
    AutoPtr mimicks std::auto_ptr<> except that it factors out the destruction
    of the owned pointer to another class.

    All methods behave identically to their std::auto_ptr<> counterparts except
    when explicitly noticed. The only extra method we add is an implicit
    conversion to bool because it is safe (in the way it is implemented here)
    and very useful.

    None of the methods of this class throw exceptions, except if T dtor throws
    (it shouldn't!)

    Template parameters are:
        - T is the type of the pointee
        - D is the type used for deleting the objects of type T and is by
            default AutoPtrDeleter which simply calls delete on the pointer
 */
template <class T, class D = ::ito33::AutoPtrDeleter<T> >
class AutoPtr
{
public:
  /// same as auto_ptr<>::element_type
  typedef T element_type;

  /// synonym for element_type, doesn't exist in auto_ptr
  typedef T value_type;

  /// synonym for pointer to element_type, doesn't exist in auto_ptr
  typedef T *pointer;

  /// a scalar type which doesn't risk to be converted to anything
  typedef T *(AutoPtr<T, D>::*unspecified_bool_type)() const;


  /// takes ownership of the pointer (which may, and by default is, @c NULL)
  explicit AutoPtr(T *ptr = NULL) : m_ptr(ptr) { }

  /// pseudo copy ctor: it has ownership transfer semantics, as auto_ptr one
  AutoPtr(const AutoPtr<T>& other) : m_ptr(other.release()) { }

  /// another pseudo copy ctor: this time from a real auto_ptr
  AutoPtr(std::auto_ptr<T> other)
    : m_ptr(const_cast< std::auto_ptr<T>& >(other).release()) { }

  /// dtor frees the pointer using the provided delete class
  ~AutoPtr() { D::Free(m_ptr); }

  /// returns the stored pointer, may be @c NULL
  T *get() const { return m_ptr; }

  /// allows use of this class as a pointer, must be non @c NULL
  T& operator*() const { return *get(); }

  /// allows use of this class as a pointer, must be non @c NULL
  T *operator->() const { return get(); }

  /// releases ownership of the held pointer and returns it
  T *release() const // this is a misnomer
  {
    T *p = m_ptr;
    const_cast<AutoPtr *>(this)->m_ptr = NULL;

    return p;
  }

  /// replaces the stored pointer with the given one
  void reset(T *p = NULL)
  {
    if ( p != m_ptr )
    {
      D::Free(m_ptr);

      m_ptr = p;
    }
    //else: don't delete the pointer held if we're reset to the same one!
  }

  /// assignment operator has ownership transfer semantics, as auto_ptr one
  AutoPtr<T>& operator=(const AutoPtr<T>& other)
  {
    reset(other.release());

    return *this;
  }

  /**
      Implicit, but safe, conversion to bool.

      Having this conversion ensures that we can work with the objects of
      this type as with the plain pointers but using a unique type (instead
      of bool or int) ensures that we don't risk to implicitly confuse the
      unrelated pointers.
   */
  operator unspecified_bool_type() const // never throws
  {
    return m_ptr ? &AutoPtr<T, D>::get : NULL;
  }

private:
  T *m_ptr;
};

} // namespace ito33

#endif // _ITO33_AUTOPTR_H_


