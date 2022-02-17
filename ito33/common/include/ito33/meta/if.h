/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/meta/if.h
// Purpose:     declares If<> metaprogramming construct
// Author:      Vadim Zeitlin
// Created:     2005-07-25
// RCS-ID:      $Id: if.h,v 1.2 2005/07/26 13:43:37 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/meta/if.h
    @brief Declaration of If template.

    This template allows to select one or another type at compile-time
    depending on the result of (compile-time) boolean constant.
 */

#ifndef _ITO33_META_IF_H_
#define _ITO33_META_IF_H_

namespace ito33
{

namespace meta
{

/**
    If<> template defines a nested Value type which is the same as T if cond is
    true and same as U otherwise.

    Template parameters:
      - cond a boolean compile-time constant
      - T the result type if cond == true
      - U the result type if cond == false
 */
template <bool cond, typename T, typename U>
struct If
{
  // by default we assume cond is true
  typedef T Value;
};

// partial specialization for the "else" part
template <typename T, typename U>
struct If<false, T, U>
{
  typedef U Value;
};

} // namespace meta

} // namespace ito33

#endif // _ITO33_META_IF_H_

