/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/thread/container.h
// Purpose:     thread-safe container class
// Author:      Vadim Zeitlin
// Created:     2005-10-21
// RCS-ID:      $Id: container.h,v 1.6 2006/03/07 12:54:58 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/thread/container.h
    @brief Declaration of thread-safe Container template.
 */

#ifndef _ITO33_THREAD_CONTAINER_H_
#define _ITO33_THREAD_CONTAINER_H_

#include "ito33/vector.h"
#include "ito33/thread.h"
#include "ito33/typenull.h"

#include <boost/function.hpp>

namespace ito33
{

namespace thread
{

namespace Private
{

/**
    The generic callback class.

    This can be used for callbacks which don't have any run-time state at all.
    Otherwise see the specialization for boost::function.

    Also note that we have to define exactly the same class twice because
    otherwise we get a problem with Container<T, U, U> as it would try to
    derive from the same ContainerCallback<U> twice.
 */
#define PRIVATE_DECLARE_CONTAINER_CALLBACK(n)                                 \
                                                                              \
template <typename T>                                                         \
class ContainerCallback##n                                                    \
{                                                                             \
protected:                                                                    \
  typedef T Functor;                                                          \
                                                                              \
  ContainerCallback##n() { }                                                  \
                                                                              \
  static void Call() { Functor(); }                                           \
};                                                                            \
                                                                              \
/**                                                                           \
    When no callback is specified, don't call anything.                       \
 */                                                                           \
template <>                                                                   \
class ContainerCallback##n<Type::Null>                                        \
{                                                                             \
protected:                                                                    \
  ContainerCallback##n() { }                                                  \
                                                                              \
  static void Call() { }                                                      \
};                                                                            \
                                                                              \
/**                                                                           \
    Callback with a state implemented using boost::function.                  \
 */                                                                           \
template <typename T>                                                         \
class ContainerCallback##n< boost::function<T> >                              \
{                                                                             \
protected:                                                                    \
  typedef boost::function<T> Function;                                        \
                                                                              \
  ContainerCallback##n(const Function& func = Function()) : m_func(func) { }  \
                                                                              \
  void Call() { m_func(); }                                                   \
                                                                              \
private:                                                                      \
  Function m_func;                                                            \
}

PRIVATE_DECLARE_CONTAINER_CALLBACK(1);
PRIVATE_DECLARE_CONTAINER_CALLBACK(2);

#undef PRIVATE_DECLARE_CONTAINER_CALLBACK

} // namespace Private


/**
    Thread-safe wrapper around a standard sequence container.

    Operations on standard container are not, and can't be, thread-safe. We
    provide a much reduced interface compared to the full standard container
    API but at least all the operations are automatically thread-safe.

    Template parameters are:
      - T the type of container elements, as with the standard containers it
          must have value semantics
      - U type of functor to be executed when the container becomes empty
          (it shouldn't access the container itself again as this would lead to
          a deadlock because the container is locked when it is called)
      - V type of functor to be executed when the container is no longer empty
          (it shouldn't access the container neither)
      - C the type of container to use internally, must have push_back(),
          empty(), erase() and size() and support standard iteration protocol
          via nested iterator type and begin() and end() methods
      - A the allocator to use for the container C
 */
template <typename T,
          typename U = Type::Null,
          typename V = Type::Null,
          template <typename, typename> class C = std::vector,
          class A = std::allocator<T> >
class Container : private Private::ContainerCallback1<U>,
                  private Private::ContainerCallback2<V>
{
public:
  /// The type of container used internally and returned by MakeCopy()
  typedef C<T, A> Elements;

  /// The type of callback called when the container becomes empty
  typedef Private::ContainerCallback1<U> WhenEmptyCallback;

  /// The type of callback called when the container becomes non empty
  typedef Private::ContainerCallback2<V> WhenNonEmptyCallback;

  
  /**
      Default constructor for containters without stateful callbacks.
   */
  Container() { }

  /**
      Constructor doing extra callbacks initialization.
   */
  template <typename EmptyArg, typename NonEmptyArg>
  Container(const EmptyArg& emptyArg, const NonEmptyArg& nonEmptyArg)
    : WhenEmptyCallback(emptyArg),
      WhenNonEmptyCallback(nonEmptyArg)
  {
  }


  /**
      Returns true if the element exists in the container.

      You may want to call this method before Add() to ensure that the
      container doesn't have this element already.

      @param x the element to look for
   */
  bool Exists(const T& x) const
  {
    Lock<CriticalSection> lock(m_cs);

    return ExistsUnlocked(x);
  }

  /**
      Adds an element to the container.

      The element shouldn't occur more than once in the container as otherwise
      Remove() wouldn't work correctly, hence an assert is raised if the
      element is already present.

      If this is the first container element, WhenNonEmptyCallback is invoked.

      @param x the element to add (a copy is made)

      @see UpdateOrAdd
   */
  void Add(const T& x)
  {
    Lock<CriticalSection> lock(m_cs);

    ASSERT_MSG( !ExistsUnlocked(x), "element can't be added twice to container" );

    AddUnlocked(x);
  }

  /**
      Same as Add() but doesn't check that the element isn't already present in
      the container.
   */
  void AddNonUnique(const T& x)
  {
    Lock<CriticalSection> lock(m_cs);

    AddUnlocked(x);
  }

  /**
      Remove the given element from the container.

      If the container becomes empty as the resut, call the WhenEmptyCallback.

      @param x the element to remove if it's present in the container
      @return true if the element was removed, false if it wasn't found
   */
  bool Remove(const T& x)
  {
    Lock<CriticalSection> lock(m_cs);

    const typename Elements::iterator end = m_elements.end();
    typename Elements::iterator i = std::find(m_elements.begin(), end, x);
    if ( i == end )
      return false;

    m_elements.erase(i);

    if ( m_elements.empty() )
      WhenEmptyCallback::Call();

    return true;
  }

  /**
      Remove all elements from the container.

      This method inconditionally invokes WhenEmptyCallback.
   */
  void Clear()
  {
    Lock<CriticalSection> lock(m_cs);

    m_elements.clear();

    WhenEmptyCallback::Call();
  }

  /**
      Return a snapshot of the container at the current moment.

      This method allows to examine the container without preventing it from
      being modified from the other threads. This is completely safe but not
      very efficient for large containers, of course.

      @see Move
   */
  Elements MakeCopy() const
  {
    Lock<CriticalSection> lock(m_cs);

    return Elements(m_elements);
  }

  /**
      Return a snapshot of the container at the current moment, like MakeCopy
      does, and clears the content of the container.

      This method allows to move data out of the container in MT-safe way.

      This method inconditionally invokes WhenEmptyCallback.

      @see MakeCopy
   */
  Elements Move()
  {
    Lock<CriticalSection> lock(m_cs);

    Elements e(m_elements);
    m_elements.clear();
    
    WhenEmptyCallback::Call();

    return e;
  }

  /**
      Either updates existing element in the array by replacing it with a new
      value or adds a new element to the container, depending on the return
      value of supplied function.

      Function @a func is executed for every element in the container
      sequentially until it returns true or the end of container is reached.
      If @a func returns true for any element, it is replaced with @a x.
      If it returns false for all elements, @a x is added to the end of the
      container.

      As a side effect, when we add the first element to a previously empty
      container, the WhenNonEmptyCallback is invoked.

      @param func function used to decide whether to update an element or not
      @param x    the element to add or update with (a copy is made)

      @see Add
   */
  void UpdateOrAdd(boost::function<bool (const T&)> func, const T& x)
  {
    Lock<CriticalSection> lock(m_cs);

    typename Elements::iterator i = std::find_if
                                    (
                                      m_elements.begin(),
                                      m_elements.end(),
                                      func
                                    );

    if ( i == m_elements.end() )
    {
        m_elements.push_back(x);

        if ( m_elements.size() == 1 )
          WhenNonEmptyCallback::Call();
    }
    else
    {
        *i = x;
    }
  }

private:
  // common part of Exists() and Add() (in debug build): supposes the container
  // is locked by the caller
  bool ExistsUnlocked(const T& x) const
  {
    const typename Elements::const_iterator end = m_elements.end();
    return std::find(m_elements.begin(), end, x) != end;
  }

  // common part of Add() and AddNonUnique(): really adds the element to the
  // container (supposes thay it's locked by the caller)
  void AddUnlocked(const T& x)
  {
    m_elements.push_back(x);

    if ( m_elements.size() == 1 )
      WhenNonEmptyCallback::Call();
  }


  // the real container
  Elements m_elements;

  // and the critical section to protect it
  mutable CriticalSection m_cs;
};

} // namespace thread

} // namespace ito33

#endif // _ITO33_THREAD_CONTAINER_H_

