/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/typeinfo.h
// Purpose:     enhanced type info class
// Author:      Vadim Zeitlin
// Created:     Dec 20, 2003
// RCS-ID:      $Id: typeinfo.h,v 1.4 2004/10/05 09:13:35 pedro Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file  ito33/typeinfo.h
    @brief More programmer-friendly replacement for std::type_info.

    The code in this file is very similar to TypeInfo from Andrei Alexandrescu
    book "Modern C++ Design: Generic Programming and Design", please refer
    there for more details.
 */

#ifndef _ITO33_TYPEINFO_H_
#define _ITO33_TYPEINFO_H_

#include "ito33/beforestd.h"
#include <typeinfo>
#include "ito33/afterstd.h"

namespace ito33
{

/**
    TypeInfo wraps a std::type_info but, unlike the latter, can be stored
    inside the containers.

    According to the Standard, type_info objects can't be copied so it is
    impossible to store them in any containers. This is really a pity because
    it is sometimes really useful to do it (e.g. for the double dispatch
    implementation) and so we simply wrap a (pointer to) type_info in a class
    which does provide value semantics and so can be stored in containers.
 */
class TypeInfo
{
public:
  /**
      Default ctor puts object in an invalid state.

      This ctor is required for an object to be used with standard containers.
   */
  TypeInfo()
  {
    // initializing the pointer to some unique but non-NULL value allows us to
    // make the remaining code much simpler as we don't have to ever test for
    // NULL
    class NoType {};
    m_pTI = &typeid(NoType); 
  }

  /**
      Conversion ctor from std::type_info.

      This (non explicit) ctor is useful to be able to mix type_info and
      TypeInfo objects together in the expressions.
   */
  TypeInfo(const std::type_info& ti) : m_pTI(&ti) { }

  // default copy ctor, assignment operator and dtor are ok

  /**
      std::type_info-compatible comparison method.

      This method is not strictly speaking necessary but having it allows to
      use TypeInfoexactly in the same way as type_info (although it is not the
      best/easiest way to use this class).
   */
  bool before(const TypeInfo& other) const
  {
    // type_info::before() returns int, use explicit comparison to suppress a
    // warning
    return (m_pTI->before(*other.m_pTI)) != 0;
  }

  /**
      Get the name of the type we represent the type information for.

      As with before(), the name of this method is the same as in
      std::type_info so, again, TypeInfo objects can be used in the same way as
      type_info ones.
   */
  const char *name() const { return m_pTI->name(); }

  /**
      Comparison with another TypeInfo.
   */
  bool operator==(const TypeInfo& other) const
  {
    // incredibly enough, type_info::operator==() returns int, not bool, so we
    // have to use the comparison below to prevent a warning from VC++ about
    // implicit int to bool conversion
    return (*m_pTI == *other.m_pTI) != 0;
  }

private:
  const std::type_info *m_pTI;
};

// ----------------------------------------------------------------------------
// all the other comparison operators fpr TypeInfo
// ----------------------------------------------------------------------------

inline bool operator<(const TypeInfo& lhs, const TypeInfo& rhs)
  { return lhs.before(rhs); }
inline bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs == rhs); }
inline bool operator>(const TypeInfo& lhs, const TypeInfo& rhs)
  { return rhs < lhs; }
inline bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs > rhs); }
inline bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs)
  { return !(lhs < rhs); }

} // namespace ito33

#endif // _ITO33_TYPEINFO_H_

