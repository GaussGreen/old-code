/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/typenull.h
// Purpose:     Type::Null declaration
// Author:      Vadim Zeitlin
// Created:     2005-11-04
// RCS-ID:      $Id: typenull.h,v 1.1 2005/11/05 00:24:28 zeitlin Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/typenull.h
    @brief Declaration of a distinct class which doesn't represent any type.
 */

#ifndef _ITO33_TYPENULL_H_
#define _ITO33_TYPENULL_H_

namespace ito33
{

/**
    Contains all metaprogramming helpers for working with types.
 */
namespace Type
{

/**
    This class is for types what @c NULL is for pointers.

    It can be used whenever a distinct nil value is needed.
 */
class Null
{
};

} // namespace Type

} // namespace ito33

#endif // _ITO33_TYPENULL_H_

