/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/string.h
// Purpose:     includes standard string header and other useful string stuff
// Author:      Vadim Zeitlin
// Created:     24.12.02
// RCS-ID:      $Id: string.h,v 1.24 2005/12/22 17:10:14 vaclav Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file   ito33/string.h
    @brief  Wrapper for the standard string header and more.

    There are two reasons for having this header:
    - having to include before/afterstd.h is annoying, especially for a
      header which is as common as string
    - we also define a couple of other useful functions here in String
      namespace
 */

#ifndef _ITO33_STRING_H_
#define _ITO33_STRING_H_

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/common.h"

#include <stdarg.h>                 // for va_list

#ifdef HAVE_WCHAR_T
  #include <wchar.h>
#endif

namespace ito33
{

namespace String
{

/**
    Formatted output: this is a safe sprintf() replacement.

    Note that it is safe in the sense that it will never overflow the buffer
    because the buffer is allocated dynamically depending on the string size,
    but it is still as unsafe as any other printf() function if the arguments
    passed to it don't correspond to the format string.

    @param format the printf()-like format
    @param ... the arguments (attention, no compile time checking is done
               except if we're using a new enough gcc!)
 */
extern std::string Printf(const char *format, ...) ATTRIBUTE_PRINTF_1;

/**
    Printf() helper.

    @sa Printf

    @param format the printf()-like format
    @param argptr the other arguments
 */
extern std::string PrintfV(const char *format, va_list argptr);

/// @name Functions supplementing the standard C library string functions
//@{

/**
    Returns true if the given string is not empty.

    This function handles both empty and @c NULL strings correctly.
 */
inline bool IsEmpty(const char *str) { return !str || *str == '\0'; }

/**
    Returns true if the strings are equal regardless of case.

    Unfortunately there is no standard function doing case insensitive
    comparison, so we need our own one.

    @param s1 the first string to compare
    @param s2 the second string to compare
    @return 0 if strings are equal, negative value if s1 precedes s2 in
            lexicographic order and positive value otherwise
 */
inline int Stricmp(const char *s1, const char *s2)
{
#ifdef _MSC_VER
  return _stricmp(s1, s2);
#elif HAVE_STRCASECMP
  return strcasecmp(s1, s2);
#else
  #error "No case-insensitive string comparison function defined"
#endif
}

/**
    Returns true if the beginning of the strings are equal regardless of case.

    This is the same as Stricmp() except that it stops after comparing first n
    characters.

    @param s1 the first string to compare
    @param s2 the second string to compare
    @param n the maximum number of characters to examine
    @return 0 if strings are equal (up to n characters), negative value if s1
            precedes s2 in lexicographic order and positive value otherwise
 */
inline int Strnicmp(const char *s1, const char *s2, size_t n)
{
#ifdef _MSC_VER
  return _strnicmp(s1, s2, n);
#elif HAVE_STRNCASECMP
  return strncasecmp(s1, s2, n);
#else
  #error "Did you forget to define HAVE_STRNCASECMP?"
#endif
}

/// The same as Stricmp() above but for string objects
inline int Stricmp(const std::string& s1, const std::string& s2)
{
  return Stricmp(s1.c_str(), s2.c_str());
}

/// The same as Strnicmp() above but for string objects
inline int Strnicmp(const std::string& s1, const std::string& s2, size_t n)
{
  return Strnicmp(s1.c_str(), s2.c_str(), n);
}

/**
    Remove all ' ' white space characters from the begining and 
    the end of the string.
 */
extern void Trim(std::string& str);

/** 
    Convenience overload, returns the trimmed string as result.
 */
extern std::string Trim(const std::string& str);

//@}


/// @name String to number conversions
//@{

/**
    Convert the string to a signed integer.

    @param s string containing the number
    @param val the pointer to the output value, cannot be @c NULL
    @param base the base for the number, 10 by default (0 is special and means
           to use C conventions)
    @return true if string is a number, false otherwise (including the case
            when the string does contain a number but it is followed by
            something else)
 */
extern bool ToLong(const std::string& s, long *val, int base = 10);

/**
    Convert the string to a unsigned integer.

    @param s string containing the number
    @param val the pointer to the output value, cannot be @c NULL
    @param base the base for the number, 10 by default (0 is special and means
           to use C conventions)
    @return true if string is a number, false otherwise (including the case
            when the string does contain a number but it is followed by
            something else)
 */
extern bool ToULong(const std::string& s, unsigned long *val, int base = 10);

/**
    Convert the string to a real number using current locale's formatting.

    @param s string containing the number
    @param val the pointer to the output value, cannot be @c NULL
    @return true if string is a number, false otherwise (including the case
            when the string does contain a number but it is followed by
            something else)

    @sa ToCDouble
 */
extern bool ToDouble(const std::string& s, double *val);

/**
    Convert the string to a real number using "C" locale's formatting
    (radix character is '.').

    @param s string containing the number
    @param val the pointer to the output value, cannot be @c NULL
    @return true if string is a number, false otherwise (including the case
            when the string does contain a number but it is followed by
            something else)

    @sa ToDouble, FromCDouble
 */
extern bool ToCDouble(const std::string& s, double *val);

/**
    Convert real number to string in decimal notation, with '.' radix character
    (i.e. "C" or English locale).

    @param val        the number
    @param precision  max number of significant digits (both before and
                      after decimal point; -1 means "use smallest needed
                      number of digits"

    @return string containing the value

    @sa ToCDouble
 */
extern std::string FromCDouble(double val, int precision = 15);

//@}


/**
    Strict Weak Ordering functor for case-insensitive string comparision.
 */
struct CmpNoCase
{
  bool operator()(const std::string& s1, const std::string& s2) const
  {
    return Stricmp(s1, s2) < 0;
  }
};


#ifdef HAVE_WCHAR_T

/**
    Translates a narrow string into a wide one.

    This class is used for multibyte to wide char translation using the
    standard C mbstowcs() function -- and hence the current locale.
 */
class MB2WC
{
public:
  /**
      Ctor converts the given string in a wide char one and stores the result.

      The passed in string shouldn't be @c NULL as it doesn't make sense but
      this is handled silently, i.e. nothing bad happens (we simply return @c
      NULL from our accessor to the wide string as well).

      However if there is a problem during conversion, we throw an exception.
      We also throw if we memory allocation fails.

      @param str the multibyte string to convert
   */
  MB2WC(const char *str) { Init(str); }

  /**
      Same as the other ctor but for std::string.
   */
  MB2WC(const std::string& str) { Init(str.c_str()); }

  /**
      Copy ctor has auto_ptr<>-like move semantics.

      To be more precise, the source object here is not const and is modified
      -- made invalid -- after it is copied.
    */
  MB2WC(const MB2WC& other)
  {
    m_wcs = other.m_wcs;
    const_cast<MB2WC &>(other).m_wcs = NULL;
  }

  /// implicit conversion to a wide string
  operator wchar_t *() const { return m_wcs; }

  /**
      Destructor frees the string created in constructor.

      Notice that destructor is not virtual and hence this class must not be
      used polymorphically.
   */
  ~MB2WC() { free(m_wcs); }

private:
  // common part of both ctors
  void Init(const char *str);

  // the wide char string we allocated (and will free)
  wchar_t *m_wcs;

  // can't be reassigned
  MB2WC& operator=(const MB2WC&);
};

/**
    Translates a wide string into a narrow one.

    This class is used for wide to multibyte char translation using the
    standard C wcstombs() function -- and hence the current locale.
 */
class WC2MB
{
public:
  /**
      Ctor converts the given string in a multibyte one and stores the result.

      The passed in string shouldn't be @c NULL as it doesn't make sense but
      this is handled silently, i.e. nothing bad happens (we simply return @c
      NULL from our accessor to the multibyte string as well).

      However if there is a problem during conversion, we throw an exception.
      We also throw if we memory allocation fails.

      @param wcs the wide string to convert
   */
  WC2MB(const wchar_t *wcs) { Init(wcs); }

  /**
      Same as the other ctor but for std::string.
   */
  WC2MB(const std::wstring& wcs) { Init(wcs.c_str()); }

  /**
      Copy ctor has auto_ptr<>-like move semantics.

      To be more precise, the source object here is not const and is modified
      -- made invalid -- after it is copied.
    */
  WC2MB(const WC2MB& other)
  {
    m_str = other.m_str;
    const_cast<WC2MB &>(other).m_str = NULL;
  }

  /// implicit conversion to a C string
  operator const char *() const { return m_str; }

  /**
      Destructor frees the string created in constructor.

      Notice that destructor is not virtual and hence this class must not be
      used polymorphically.
   */
  ~WC2MB() { free(m_str); }

private:
  // common part of both ctors
  void Init(const wchar_t *wcs);

  // the multibyte char string we allocated (and will free)
  char *m_str;

  // can't be copied
  WC2MB& operator=(const WC2MB&);
};

#endif // HAVE_WCHAR_T

} // namespace String

} // namespace ito33

#endif // _ITO33_STRING_H_

