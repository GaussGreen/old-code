/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/jni/translateoutput.h
// Purpose:     Specialization for Translate of finance::ModelOutput
// Created:     2006/07/01
// RCS-ID:      $Id: translateoutput.h,v 1.2 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/jni/translateoutput.h
    @brief Specialization for modeloutput translation

    Many functions return a ModelOutput which might be a specific modeloutput
    (bondlikeoutput, cboptionoutput etc). Currently, cpp2any can't handle it
    automatically(well, without this file).
 */

#ifndef _ITO33_FINANCE_JNI_TRANSLATEOUTPUT_H_
#define _ITO33_FINANCE_JNI_TRANSLATEOUTPUT_H_

#include "ito33/jni/c2a.h"

#include "modeloutput.h"
#include "bondlikeoutput.h"
#include "cboptionoutput.h"

namespace C2A
{

namespace Java
{

/**
    Specialization for shared pointer of ModelOutput.
 */
template <>
struct Translate< shared_ptr<ito33::finance::ModelOutput> >
{
  typedef shared_ptr<ito33::finance::ModelOutput> CppType;
  typedef jobject JNIType;

  static CppType From(JNI::Env& env, jobject obj)
  {
    if ( !obj )
      throw EXCEPTION_MSG(ITO33_NULL_PARAM, "Object is null.");

    JNI::FieldId<jlong> idSelf("m_self");
    const jlong ptr = env.GetLongField(obj, idSelf);

    // Java object keeps a shared pointer containing the real object
    return *reinterpret_cast<CppType *>(ptr);
  }

  static jobject To(JNI::Env& env, const CppType& ptr)
  {
    if ( dynamic_cast<ito33::finance::BondLikeOutput *>( ptr.get() ) != 0 )
      return Translate< shared_ptr<ito33::finance::BondLikeOutput> >::To
             ( env, ito33::static_pointer_cast<ito33::finance::BondLikeOutput>(ptr) );
    // CBOptionOutput
    else if ( dynamic_cast<ito33::finance::CBOptionOutput *>( ptr.get() ) != 0 )
      return Translate< shared_ptr<ito33::finance::CBOptionOutput> >::To
             ( env, ito33::static_pointer_cast<ito33::finance::CBOptionOutput>(ptr) );
    // basic output - option, cds, par bond, eds etc.
    else
    {
      JNI::JClass cls(env, TypeInfo<ito33::finance::ModelOutput>::GetFullName());
      JNI::MethodIdAny idCtor(env, cls, "<init>", "(J)V");

      return env.NewObject(cls, idCtor, reinterpret_cast<jlong>(&ptr));
    }
  }
};

} // namespace Java

} // namespace C2A

#endif // _ITO33_FINANCE_JNI_TRANSLATEOUTPUT_H_
