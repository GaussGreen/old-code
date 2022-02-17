/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/com/translateoutput.h
// Purpose:     Specialization for modeloutput translation
// Created:     2006/06/30
// RCS-ID:      $Id: translateoutput.h,v 1.2 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/com/translateoutput.h
    @brief Specialization for modeloutput translation

    Many functions return a ModelOutput which might be a specific modeloutput
    (bondlikeoutput, cboptionoutput etc). Currently, cpp2any can't handle it
    automatically(well, without this file).
 */

#ifndef _ITO33_FINANCE_COM_TRANSLATEOUTPUT_H_
#define _ITO33_FINANCE_COM_TRANSLATEOUTPUT_H_

#include "com/modeloutput_impl.h"        // COM interface implementation wrapper
#include "com/bondlikeoutput_impl.h"     // COM interface implementation wrapper
#include "com/cboptionoutput_impl.h"     // COM interface implementation wrapper

#include "ito33/com/c2a.h"

namespace C2A
{

namespace COM
{

// ----------------------------------------------------------------------------
// Specification for ModelOutput
// ----------------------------------------------------------------------------

template<>
struct Translate<IModelOutput *>
{
  typedef IModelOutput Iface;

  typedef Iface *COMType;

  // there is no CppType typedef because many C++ types map to COM interface
  // pointer type

  typedef ::ito33::COM::Traits<Iface> Traits;
  typedef Traits::Impl Impl;
  typedef Traits::Class Class;

  static ::ito33::shared_ptr<Class> From(Iface *p)
  {
    if ( !p )
    {
      using ito33::Exception;
      throw EXCEPTION_MSG(ITO33_NULL_PARAM, "Object is null.");
    }

    return static_cast<Impl *>(p)->GetImpl();
  }

  static Iface *To(const ::ito33::shared_ptr<Class>& ptr)
  {
    // We may use a visitor in the future if there are many different types
    // of output

    // BondLikeOutput
    if ( dynamic_cast<ito33::finance::BondLikeOutput *>( ptr.get() ) != 0 )
      return Translate<IBondLikeOutput *>::To
             ( ito33::static_pointer_cast<ito33::finance::BondLikeOutput>(ptr) );
    // CBOptionOutput
    else if ( dynamic_cast<ito33::finance::CBOptionOutput *>( ptr.get() ) != 0 )
      return Translate<ICBOptionOutput *>::To
             ( ito33::static_pointer_cast<ito33::finance::CBOptionOutput>(ptr) );
    // basic output - option, cds, par bond, eds etc.
    else
      return ptr ? new Impl(ptr) : NULL;
  }

  static Iface *To(const ::ito33::AutoPtr<Class>& ptr)
  {
    return ptr ? To(ito33::shared_ptr<Class>(ptr.release())) : NULL;
  }

  static Iface *To(const ::std::auto_ptr<Class>& ptr)
  {
    return ptr.get() ? To(ito33::shared_ptr<Class>(
          const_cast<std::auto_ptr<Class>&>(ptr).release())
        )
      : NULL;
  }

  static Iface *To(const Class& obj)
  {
    return To(::ito33::shared_ptr<Class>(new Class(obj)));
  }

  static Iface *To(const Class *p)
  {
    return p ? To(*p) : NULL;
  }
};

} // namespace COM

} // namespace C2A

#endif // _ITO33_FINANCE_COM_TRANSLATEOUTPUT_H_
