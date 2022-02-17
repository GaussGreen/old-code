/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/meta/typetraits.h
// Purpose:     Miscellaneous operations on types.
// Author:      Vadim Zeitlin
// Created:     2005-07-26
// RCS-ID:      $Id: typetraits.h,v 1.3 2005/09/06 07:42:35 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/meta/typetraits.h
    @brief Miscellaneous templates operating on types.
 */

#ifndef _ITO33_META_TYPETRAITS_H_
#define _ITO33_META_TYPETRAITS_H_

namespace ito33
{

/**
    Contains all metaprogramming helper templates.
 */
namespace meta
{

/**
    IsSame defines a nested Value constant which is true if the 2 types are the
    same and false otherwise.
 */
template <typename T, typename U>
struct IsSame
{
  enum { Value = false };
};

// partial specialization for the case of same types
template <typename T>
struct IsSame<T, T>
{
  enum { Value = true };
};

/**
    ArgType defines a nested Value type which is type suitable for passing
    arguments to a function and is either T or const T&.
 */
template <typename T>
struct ArgType
{
  typedef const T& Value;
};

// specializations for primitive types:
#define ITO33_ARGTYPE_PRIMITIVE_TYPE(t)         \
    template <> struct ArgType<t>               \
    {                                           \
      typedef const t Value;                    \
    };

ITO33_ARGTYPE_PRIMITIVE_TYPE(bool)
ITO33_ARGTYPE_PRIMITIVE_TYPE(unsigned char)
ITO33_ARGTYPE_PRIMITIVE_TYPE(signed char)
ITO33_ARGTYPE_PRIMITIVE_TYPE(unsigned int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(signed int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(unsigned short int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(signed short int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(signed long int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(unsigned long int)
ITO33_ARGTYPE_PRIMITIVE_TYPE(float)
ITO33_ARGTYPE_PRIMITIVE_TYPE(double)
ITO33_ARGTYPE_PRIMITIVE_TYPE(long double)
#ifndef _MSC_VER
    // VC++ doesn't have wchar_t as primitive type
    ITO33_ARGTYPE_PRIMITIVE_TYPE(wchar_t)
#endif

#undef ITO33_ARGTYPE_PRIMITIVE_TYPE

} // namespace meta

} // namespace ito33

#endif // _ITO33_META_TYPETRAITS_H_
