/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/typelist.h
// Purpose:     type lists support
// Author:      Vadim Zeitlin
// Created:     27.03.03
// RCS-ID:      $Id: typelist.h,v 1.8 2005/11/05 00:24:28 zeitlin Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/typelist.h
    @brief  Template-based type lists support.

    The code in this file is based on the ideas from Andrei Alexandrescu book
    "Modern C++ Design: Generic Programming and Design" and reading it (at
    least the beginning) could be helpful to get the big picture of what is
    going on here.
 */

#ifndef _ITO33_TYPE_LIST_H_
#define _ITO33_TYPE_LIST_H_

#include "ito33/typenull.h"

namespace ito33
{

namespace Type
{

/**
    This is a linked list of types which may be of any length.

    The type lists are implemented as a simply linked lists with a head which
    represents the first element of the list and a tail which may be either
    Null (then there are no other elements) or another type list.

    To create the type lists of arbitrary length the convenience macros
    TYPE_LIST_N are provided below.

    The template parameters are:
        - H is the first element of the typelist, a "simple" type usually
        - T is either Null or another type list
 */
template <class H, class T>
struct List
{
  /// the head of the linked list (normally a simple type)
  typedef H Head;

  /// the tail of the linked list (either another type list or Type::Null)
  typedef T Tail;
};

/**
    Append a simple type to the type list.

    Here we provide the forward declaration only, to define this template we
    need the helper types from inside Impl namespace below.

    Template parameters:
        - TL a type list
        - T a simple type
 */
template <class TL, typename T>
struct Append;

// This namespace is used to hide the implementation details only.
namespace Impl
{
  /**
      Append is implemented by recursively shifting the elements from the
      original type list to the new one.

      This class provides the general recursion step.
   */
  template <class TL>
  struct AppendImpl
  {
    /// helper needed to emulate partial specialization with VC++
    template <typename T>
    struct Inner
    {
      /// the type of the resulting list
      typedef List<typename TL::Head,
            typename Append<typename TL::Tail, T> > Result;
    };
  };

  /// This speciailization stops the recursion
  template <>
  struct AppendImpl<Null>
  {
    /// helper needed to emulate partial specialization with VC++
    template <typename T>
    struct Inner
    {
      /// the type of the resulting list
      typedef List<T, Null> Result;
    };
  };
} // namespace Impl

/// now we can define the public template
template <class TL, typename T>
struct Append : Impl::AppendImpl<TL>::template Inner<T>::Result
{
};

} // namespace Type

} // namespace ito33

/**
    Macros for simpler type list declaration.

    We define macros TYPE_LIST_N(T1, ..., TN) for definition of the type lists
    with N elements. N is currently limited to 10 but it is trivial to write
    additional macros if we ever need longer type lists.
 */

/// macro to define type list of 1 element
#define TYPE_LIST_1(T1) ::ito33::Type::List< T1, ::ito33::Type::Null>

/// macro to define type list of 2 elements
#define TYPE_LIST_2(T1, T2) ::ito33::Type::List< T1, TYPE_LIST_1(T2) >

/// macro to define type list of 3 elements
#define TYPE_LIST_3(T1, T2, T3) ::ito33::Type::List< T1, TYPE_LIST_2(T2, T3) >

/// macro to define type list of 4 elements
#define TYPE_LIST_4(T1, T2, T3, T4) \
  ::ito33::Type::List< T1, TYPE_LIST_3(T2, T3, T4) >

/// macro to define type list of 5 elements
#define TYPE_LIST_5(T1, T2, T3, T4, T5) \
  ::ito33::Type::List< T1, TYPE_LIST_4(T2, T3, T4, T5) >

/// macro to define type list of 6 elements
#define TYPE_LIST_6(T1, T2, T3, T4, T5, T6) \
  ::ito33::Type::List< T1, TYPE_LIST_5(T2, T3, T4, T5, T6) >

/// macro to define type list of 7 elements
#define TYPE_LIST_7(T1, T2, T3, T4, T5, T6, T7) \
  ::ito33::Type::List< T1, TYPE_LIST_6(T2, T3, T4, T5, T6, T7) >

/// macro to define type list of 8 elements
#define TYPE_LIST_8(T1, T2, T3, T4, T5, T6, T7, T8) \
  ::ito33::Type::List< T1, TYPE_LIST_7(T2, T3, T4, T5, T6, T7, T8) >

/// macro to define type list of 9 elements
#define TYPE_LIST_9(T1, T2, T3, T4, T5, T6, T7, T8, T9) \
  ::ito33::Type::List< T1, TYPE_LIST_8(T2, T3, T4, T5, T6, T7, T8, T9) >

/// macro to define type list of 10 elements
#define TYPE_LIST_10(T1, T2, T3, T4, T5, T6, T7, T8, T9, T10) \
  ::ito33::Type::List< T1, TYPE_LIST_9(T2, T3, T4, T5, T6, T7, T8, T9, T10) >

#endif // _ITO33_TYPE_LIST_H_

