/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/hazardratecallback.cpp
// Purpose:     Callback hazard rate class
// Author:      ZHANG Yunzhi
// Created:     2004/06/05
// RCS-ID:      $Id: hazardratecallback.cpp,v 1.9 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/ihg/hazardratecallback.h"

#include "ito33/xml/write.h"

#include "ihg/xml/hazardrate.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace ihg
{

  void HazardRateCallBack::Dump(ito33::XML::Tag& /* tagParent */) const
{
  // do nothing right now
  // TODO : how to dump user defined function?
}


//-----------------------------------------------------------------------------
// TODO: Visitor pattern for call back volatility
//-----------------------------------------------------------------------------
void HazardRateCallBack::Visit(ito33::ihg::HazardRateVisitor &) const
{
  throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          "We are sorry, but we don't know how to visit a call back class."
        );
}



} // namespace ihg

} // namespace ito33
