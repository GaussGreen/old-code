/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/enum_values_names.h
// Purpose:     mapping between enum elements and the corresponding strings
// Author:      Vadim Zeitlin
// Author:      ZHANG Yunzhi
// Created:     2004-05-10
// RCS-ID:      $Id: enum_values_names.h,v 1.6 2005/10/24 13:58:26 cosmin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/enum_values_names.h
    @brief mapping between enum elements and the corresponding strings
 */

#ifndef _ITO33_ENUM_VALUES_NAMES_H_
#define _ITO33_ENUM_VALUES_NAMES_H_

#include "ito33/beforestd.h"
#include <cstring>
#include <string>
#include "ito33/afterstd.h"

namespace ito33
{

/**
    Structure describing the mapping between enum elements and the
    corresponding strings. It is often used in XML file handling.

    Example usage:
      @code
        static const EnumValuesNames<MyEnum> myEnumValues[] =
        {
          { "foo", MyEnum_Foo },
          { "bar", MyEnum_Bar },
        };

        const char* strName = 
          GetNameFromEnumValue(MyEnum_Bar, SIZEOF(myEnumValues), myEnumValues);
      @endcode
 */
template <typename T, typename N = const char*>
struct EnumValuesNames
{
  /// enum type
  typedef T Type;
  /// string type
  typedef N String;

  String name;
  Type value;
};


/**
    Get the text value corresponding to enum value in given mapping.
    If no text is defined for the enum value, NULL is returned.

    @param value given enum value
    @param nSize number of elements in values array
    @param values the mapping between strings and enum values
    @return text value corresponding to enum value or NULL.
 */
template <typename T>
inline 
const char *GetNameFromEnumValue(T value,
                                 size_t nSize,
                                 const EnumValuesNames<T> values[])
{
  for ( size_t n = 0; n < nSize; ++n )
  {
    if ( values[n].value == value)
      return values[n].name;
  }
  
  return 0;
}
/**
    Get the text value corresponding to enum value in given mapping.
    If no text is defined for the enum value, an empty string is returned.

    @param value given enum value
    @param nSize number of elements in values array
    @param values the mapping between strings and enum values
    @return text value corresponding to enum value or string.
 */
template <typename T>
inline 
std::string GetNameFromEnumValue
            (
              T value,
              size_t nSize,
              const EnumValuesNames<T, std::string> values[]
            )
{
  for ( size_t n = 0; n < nSize; ++n )
  {
    if ( values[n].value == value)
      return values[n].name;
  }
  return std::string();
}

/**
    Restore the enum value corresponding to the text value.
    If no matching key is found, the function returns false.

    @param strName to translate to enum
    @param value (output) enum value corresponding to strName
    @param nSize number of elements in values array
    @param values the mapping between strings and enum values
    @return true if it succeeds false if not
 */
template <typename T>
inline
bool RestoreEnumValueFromName(const char * strName,
                              T& value,
                              size_t nSize,
                              const EnumValuesNames<T> values[])
{
  for ( size_t n = 0; n < nSize; ++n )
  {
    if ( strcmp(values[n].name, strName) == 0 )
    {
      value = values[n].value;
      return true;
    }
  }

  return false;
}

/**
    Restore the enum value corresponding to the text value.
    If no matching key is found, the function returns false.

    @param strName to translate to enum
    @param value (output) enum value corresponding to strName
    @param nSize number of elements in values array
    @param values the mapping between strings and enum values
    @return true if it succeeds false if not
 */
template <typename T>
inline
bool RestoreEnumValueFromName(const std::string& strName,
                              T& value,
                              size_t nSize,
                              const EnumValuesNames<T, std::string> values[])
{
  for ( size_t n = 0; n < nSize; ++n )
  {
    if ( values[n].name == strName )
    {
      value = values[n].value;
      return true;
    }
  }

  return false;
}


} // namespace ito33

#endif // _ITO33_ENUM_VALUES_NAMES_H_

