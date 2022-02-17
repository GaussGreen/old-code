/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/underlyingprocess.cpp
// Purpose:     Implementation of inhomogeneous underlying process class
// Created:     2006/06/01
// RCS-ID:      $Id: underlyingprocess.cpp,v 1.2 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004-2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/ihg/error.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/volatility.h"
#include "ito33/ihg/underlyingprocess.h"

#include "ito33/xml/write.h"
#include "ito33/xml/write_vector.h"

#include "ito33/xml/finance/underlyingprocess.h"
#include "ihg/xml/underlyingprocess.h"

extern const ito33::ihg::Error ITO33_IHG_NULL_VOLATILITY, 
                               ITO33_IHG_NULL_HAZARDRATE;

namespace ito33
{

namespace ihg
{ 

void UnderlyingProcess::CheckVolatility() const
{
  CHECK_COND(m_pVolatility, ITO33_IHG_NULL_VOLATILITY);
}

void UnderlyingProcess::CheckHazardRate() const
{
  CHECK_COND(m_pHazardRate, ITO33_IHG_NULL_HAZARDRATE);
} 

void UnderlyingProcess::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_UNDERLYING_PROCESS, tagParent);

  // Dump the volatility
  GetVolatility()->Dump(tag);

  // Dump the hazard rate dump
  GetHazardRate()->Dump(tag);

  // Dump the data of finance::UnderlyingProcess
  DumpMe(tag);
}     


} // namespace ihg

} // namespace ito33
