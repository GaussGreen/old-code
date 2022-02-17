/////////////////////////////////////////////////////////////////////////////
// Name:        utils/string.cpp
// Purpose:     implements functions from ito33::String namespace
// Author:      Vadim Zeitlin
// Created:     29.01.03
// RCS-ID:      $Id: string.cpp,v 1.14 2006/05/27 20:05:06 zhang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  // '_vsnprintf', 'mbstowcs' and 'wcstombs' were declared deprecated in VC8.
  // Ignore warning C4996
  #define _CRT_SECURE_NO_DEPRECATE 1
#endif

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"
#include "ito33/debug.h"

#include <stdarg.h>
#include <stdio.h>          // vsnprintf() is sometimes defined here
#include <stdlib.h>

#include <new>

#ifdef _MSC_VER
  #define vsnprintf _vsnprintf
#elif !defined(HAVE_VSNPRINTF)
  #error "vsnprintf() is required: #define HAVE_VSNPRINTF if it is available"
#endif

#ifdef _WIN32
  #include "ito33/win32/winwrap.h"  // for OutputDebugString()
#endif

#include <stdexcept>

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// Printf() &c
// ----------------------------------------------------------------------------

std::string ito33::String::Printf(const char *format, ...)
{
  va_list argptr;
  va_start(argptr, format);

  std::string s = PrintfV(format, argptr);

  va_end(argptr);

  return s;
}

std::string ito33::String::PrintfV(const char *format, va_list argptr)
{
  char *buf = NULL;

  int size = 1024;
  for ( ;; ) {
    char *bufNew = (char *)realloc(buf, (size + 1)*sizeof(char));
    if ( !bufNew ) {
      // out of memory (but still don't leak what we've got)
      free(buf);

      return "";
    }

    buf = bufNew;

    int len = vsnprintf(buf, size, format, argptr);

    // some implementations of vsnprintf() don't NUL terminate the string
    // if there is not enough space for it so always do it manually
    buf[size] = '\0';

    if ( len >= 0 )
    {
      // ok, there was enough space
      break;
    }

    // still not enough, double it again
    size *= 2;
  }

  // copythe buffer contents to the string object (this is grossly
  // inefficient but the standard string class doesn't allow to write
  // directly to its internal buffer)
  std::string s(buf);
  free(buf);

  return s;
}

// ----------------------------------------------------------------------------
// implementation of function supplementing the standard library
// ----------------------------------------------------------------------------
extern void ito33::String::Trim(std::string& str)
{
  std::string::size_type pos = str.find_last_not_of(' ');
  if ( pos != std::string::npos ) 
  {
    str.erase(pos + 1);
    pos = str.find_first_not_of(' ');
    if ( pos != std::string::npos ) 
      str.erase(0, pos);
  }
  else 
    str.clear();
}

extern std::string ito33::String::Trim(const std::string& str)
{
  std::string tmp(str);

  Trim(tmp);

  return tmp;
}

// ----------------------------------------------------------------------------
// debug messages
// ----------------------------------------------------------------------------

#ifndef NDEBUG

void
ito33::DebugPrintf(const char* format, ...)
{
    va_list argptr;

    va_start(argptr, format);
    std::string s = String::PrintfV(format, argptr);
    va_end(argptr);

#ifdef _WIN32
    ::OutputDebugString(s.c_str());
#else
    fprintf(stderr, "%s", s.c_str());
#endif
}

#endif // ifndef NDEBUG

// ----------------------------------------------------------------------------
// conversion to numbers
// ----------------------------------------------------------------------------

bool ito33::String::ToLong(const std::string& s, long *val, int base)
{
  CHECK( val, false, "NULL pointer in String::ToLong" );
  ASSERT_MSG( !base || (base > 1 && base <= 36), "invalid base" );

  const char *start = s.c_str();
  char *end;
  *val = strtol(start, &end, base);

  // return TRUE only if scan was stopped by the terminating NUL and if the
  // string was not empty to start with
  return !*end && (end != start);
}

bool ito33::String::ToULong(const std::string& s, unsigned long *val, int base)
{
  CHECK( val, false, "NULL pointer in String::ToULong" );
  ASSERT_MSG( !base || (base > 1 && base <= 36), "invalid base" );

  const char *start = s.c_str();
  char *end;
  *val = strtoul(start, &end, base);

  // return TRUE only if scan was stopped by the terminating NUL and if the
  // string was not empty to start with
  return !*end && (end != start);
}

bool ito33::String::ToDouble(const std::string& s, double *val)
{
  CHECK( val, false, "NULL pointer in String::ToDouble" );

  const char *start = s.c_str();
  char *end;
  *val = strtod(start, &end);

  // return TRUE only if scan was stopped by the terminating NUL and if the
  // string was not empty to start with
  return !*end && (end != start);
}

// ----------------------------------------------------------------------------
// wide char stuff
// ----------------------------------------------------------------------------

#ifdef HAVE_WCHAR_T

void ito33::String::MB2WC::Init(const char *str)
{
  if ( str )
  {
    // first calculate the size needed
    size_t len = mbstowcs(NULL, str, 0);
    if ( len == (size_t)-1 )
    {
      // error during conversion
      throw std::runtime_error("bad multibyte string");
    }

    // take NUL into account
    ++len;

    m_wcs = (wchar_t *)malloc(sizeof(wchar_t) * len);
    if ( !m_wcs )
    {
      // out of memory
      throw std::bad_alloc();
    }

    mbstowcs(m_wcs, str, len);
  }
  else // NULL pointer given
  {
    m_wcs = NULL;
  }
}

void ito33::String::WC2MB::Init(const wchar_t *wcs)
{
  if ( wcs )
  {
    // first calculate the size needed
    size_t len = wcstombs(NULL, wcs, 0);
    if ( len == (size_t)-1 )
    {
      // error during conversion
      throw std::runtime_error("bad wide character string");;
    }

    // take NUL into account
    ++len;

    m_str = (char *)malloc(sizeof(char) * len);
    if ( !m_str )
    {
      // out of memory
      throw std::bad_alloc();
    }

    wcstombs(m_str, wcs, len);
  }
  else // NULL pointer given
  {
    m_str = NULL;
  }
}

#endif // HAVE_WCHAR_T

