/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/selfptr.h
// Purpose:     SelfPointer class declaration
// Author:      Vadim Zeitlin
// Created:     2005-11-04
// RCS-ID:      $Id: selfptr.h,v 1.2 2006/02/06 13:07:41 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/selfptr.h
    @brief Helper class for objects which must contain pointer to themselves.

    It is often needed that an object of a given class keeps a weak pointer
    to itself so that it can call methods of other classes passing them a
    strong shared pointer to itself. As it's impossible to pass the pointer to
    the not yet existing object to the ctor of the object being created, the
    standard solution is to have a static Create() function which first creates
    the object and then initializes it with the pointer pointing to it.

    The SelfPointer template class implements this technique. Example of
    declaring a class using it:

    @code
      class Foo : public SelfPointer<Foo>
      {
        ... no public ctors ...
      };
    @endcode

    To create an object of class Foo you now have to call the static Create()
    method (implemented in the base class). And to pass a SharedPtr<Foo> to
    some other function you may use SharedPtr<Foo>(m_self) in any Foo method
    except its dtor because when the dtor is being executed there are no more
    shared pointers pointing to this object and so trying to create one from
    weak pointer we keep would just result in a NULL pointer.
 */

#ifndef _ITO33_SELFPTR_H_
#define _ITO33_SELFPTR_H_

#include "ito33/sharedptr.h"
#include "ito33/weakptr.h"
#include "ito33/typelist.h"

namespace ito33
{

#ifdef DOXYGEN

/**
    SelfPointer is a base class for all classes which need to keep a weak
    pointer to themselves.

    In reality, there is not one SelfPointer class but a number of classes
    SelfPointer0, SelfPointer1, ..., SelfPointerN where N can be defined when
    compiling this file by predefining MAX_CTOR_ARGS preprocessor constant and
    has sufficiently big default value if not explicitely specified. The class
    SelfPointerM should be used for a class whose ctor takes M parameters.
 */
template <class T,
          class A1 = Type::Null,
          class A2 = Type::Null,
          class A3 = Type::Null,
          class A4 = Type::Null,
          class A5 = Type::Null,
          class A6 = Type::Null,
          class A7 = Type::Null,
          class M = ::ito33::Policy::MTUnsafe>
class SelfPointer
{
public:
  /// Typedef for shared pointer returned by Create()
  typedef SharedPtr<T, M> Ptr;

  /// Typedef for weak pointer stored internally.
  typedef WeakPtr<T, M> WeakPtr;

  /**
      Create the object of type T.

      Throws only if ctor of T throws.
   */
  static Ptr Create(A1 a1, A2 a2, A3 a3, A4 a4, A5 a5, A6 a6, A7 a7)
  {
    Ptr p = new T(a1, a2, a3, a4, a5, a6, a7);
    p->m_self = p;
    return p;
  }

  /**
      Return a strong pointer to this object itself.

      This will only return a @c NULL pointer when called from the destructor,
      otherwise the pointer is always valid.
   */
  Ptr Self() const
  {
    return Ptr(m_self);
  }

protected:
  // we're never supposed to be created directly
  SelfPointer() { }


  WeakPtr m_self;
};

#else // !DOXYGEN

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/arithmetic/sub.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

#ifndef MAX_CTOR_ARGS
  #define MAX_CTOR_ARGS 10
#endif // MAX_CTOR_ARGS

// declare SelfPointerN
#define SELF_PTR_for(z, n, unused)                                            \
  template <class T,                                                          \
            BOOST_PP_ENUM_PARAMS(n, typename A)                               \
            BOOST_PP_COMMA_IF(n)                                              \
            class M = ::ito33::Policy::MTUnsafe>                              \
  class SelfPointer ## n                                                      \
  {                                                                           \
  public:                                                                     \
    typedef SelfPointer ## n<T,                                               \
                             BOOST_PP_ENUM_PARAMS(n, A)                       \
                             BOOST_PP_COMMA_IF(n)                             \
                             M> SelfPointerBase;                              \
    typedef SharedPtr<T, M> Ptr;                                              \
    typedef WeakPtr<T, M> WeakPtr;                                            \
                                                                              \
    static Ptr Create(BOOST_PP_ENUM_BINARY_PARAMS(n, A, a))                   \
    {                                                                         \
      Ptr p(new T(BOOST_PP_ENUM_PARAMS(n, a)));                               \
      p->m_self = p;                                                          \
      return p;                                                               \
    }                                                                         \
                                                                              \
    Ptr Self() const { return Ptr(m_self); }                                  \
                                                                              \
  protected:                                                                  \
    SelfPointer ## n() { }                                                    \
                                                                              \
    WeakPtr m_self;                                                           \
  };

// do it for all values of n from 0 to MAX_CTOR_ARGS-1
BOOST_PP_REPEAT(MAX_CTOR_ARGS, SELF_PTR_for, ~)

#undef SELF_PTR_for

#endif // DOXYGEN/!DOXYGEN

} // namespace ito33

#endif // _ITO33_SELFPTR_H_

