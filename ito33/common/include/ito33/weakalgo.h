/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/weakalgo.h
// Purpose:     algorithms for working with containers of weak pointers
// Author:      Vadim Zeitlin
// Created:     2006-02-07
// RCS-ID:      $Id: weakalgo.h,v 1.3 2006/02/10 16:58:01 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file include/ito33/weakalgo.h
    @brief Algorithms for working with weak pointers.

    Containers containing weak pointers need special handling because the
    pointers in them can expire and so we have to check for this before
    processing every element. Using standard algorithms is extremely
    inconvenient in this case as it requires several levels of nested lambda
    expressions so we provide our own implementations for a few common cases.
 */

#ifndef _ITO33_WEAKALGO_H_
#define _ITO33_WEAKALGO_H_

#include <utility>

namespace ito33
{

namespace Private
{

// this is a helper of ForEachWeak: as we don't know the type T inside the
// function, we can't create a SharedPtr<T> explicitely (until we have auto
// keyword with C0x semantics...) and so we have to forward to another function
template <typename T, typename M, typename F>
inline void DoForOne(const SharedPtr<T, M>& p, F& func)
{
  if ( p )
    func(*p);
}

// this helper for ForEachWeak calls DoForOne for weak pointer in type the
// iterator dereferences to (e.g. *iterator if it is weak pointer or
// iterator->second if the iterator points to std::pair)
template <typename T, typename F>
inline void ApplyForOne(const T& value, F& func)
{
  DoForOne(value.lock(), func);
}

// overload for iterators into associative containers
template <typename K, typename T, typename F>
inline void ApplyForOne(const std::pair<K, T>& value, F& func)
{
  DoForOne(value.second.lock(), func);
}

} // namespace Private

/**
    Algorithm for applying a given functor to all valid elements of a sequence.

    The expired weak pointers are simply skipped, otherwise this is the same as
    std::for_each().

    The iterators must dereference to either WeakPtr<T> (usual case of
    iterators into sequence containers) or to std::pair<K, WeakPtr<T>> (allows
    to also use the algorithm with iterators into associative containers) for
    some types T, K.

    @param begin the starting iterator of the sequence
    @param end the one past end iterator
    @param func the functor to apply to the elements (it should be callable
                with argument of type T)
    @return the functor after all applications
 */
template <typename I, typename F>
F ForEachWeak(I begin, I end, F func)
{
  for ( ; begin != end; ++begin )
    Private::ApplyForOne(*begin, func);

  return func;
}

} // namespace ito33

#endif // _ITO33_WEAKALGO_H_

