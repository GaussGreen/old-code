/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/com/traits.h
// Purpose:     COM::Traits class stores meta information about interfaces
// Author:      Vadim Zeitlin
// Created:     23.12.02
// RCS-ID:      $Id: traits.h,v 1.11 2004/10/05 09:13:36 pedro Exp $
// Copyright:   (c) 2002 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/com/traits.h
    @brief Traits for the COM interfaces.

    The trais is a template class whose specialization for the given interface
    contains useful information about it. It is used by UuidOf() function, for
    example.
 */

#ifndef _ITO33_COM_TRAITS_H_
#define _ITO33_COM_TRAITS_H_

#include "ito33/com/uuid.h"

namespace ito33
{

namespace COM
{

/**
    Traits contains meta information about the given COM interface.

    Traits is used to find the base class(es) of the interface and its IID from
    its name in the template functions. The generic Traits template doesn't
    compile, so Traits @b must be specialized for each interface used. This can
    be achieved with minimal effort using the DEFINE_COM_TRAITS macro below.

    Note that currently only one base class is supported.

    @todo support for multiple inheritance
 */
template <class Iface>
struct Traits
{
  /**
      The type of the interface this traits class corresponds to.
   */
  typedef Iface Interface;

  /**
      The type of the base interface.

      For the generic template we don't know the base class so we don't
      specify it here.
   */
  typedef class NoSuchBase Base;

  /**
      The class which implements this interface.

      The macro DEFINE_COM_IFACE() below supposes that for an interface IFoo
      this class is called FooImpl.
   */
  typedef class NoSuchImpl Impl;

  /**
      The class used by Impl to really implement the interface methods.

      In the usual scheme of things we have FooImpl contains a trivial
      implementation of IFoo which simply forwards all calls to the real class
      Foo.
   */
  typedef class NoSuchClass Class;

  /**
      Returns the IID of the interface Iface.
   */
  static const UUID& Uuid();
};

/**
    Macro to be used to create the traits class for the given interface.

    This macro is called with the interface and its base interface names. It
    supposes that the IID for the given interface has the same name with @c
    IID_ prefix.

    This macro is now superseded with DEFINE_COM_IFACE().
 */
#define DEFINE_COM_TRAITS(Iface, IfaceBase)                               \
  template <>                                                             \
  struct ::ito33::COM::Traits<Iface>                                      \
  {                                                                       \
    typedef Iface Interface;                                              \
    typedef IfaceBase Base;                                               \
    static const UUID& Uuid() { return IID_ ## Iface; }                   \
  }

/**
    Macro to be used to fully define everthing we ever want to know about the
    given interface.

    This macro supersedes DEFINE_COM_TRAITS() and so also supposes that the
    IID has the form IID_IFoo. Moreover, it presumes that the implementation
    class for IFoo is FooImpl and it uses the "real" class Foo to implement the
    methods.

    @param name the name of the interface without the leading @c "I"
    @param cls the fully qualified name of the corresponding real C++ class
    @param base the complete (i.e. usually with "I") name of the base interface
 */
#define DEFINE_COM_IFACE_FULL(name, cls, base)                                \
  class name##Impl;                                                           \
  template <>                                                                 \
  struct ::ito33::COM::Traits<I##name>                                        \
  {                                                                           \
    typedef I##name Interface;                                                \
    typedef base Base;                                                        \
    typedef name##Impl Impl;                                                  \
    typedef cls Class;                                                        \
    static const UUID& Uuid() { return IID_I ## name; }                       \
  }

/**
    This is exactly the same as DEFINE_COM_IFACE_FULL except that it supposes
    that the real class name is the same as name.
 */
#define DEFINE_COM_IFACE(name, base) DEFINE_COM_IFACE_FULL(name, name, base)

} // namespace COM

} // namespace ito33

#endif // _ITO33_COM_TRAITS_H_

