/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/callfixedshare.cpp
// Purpose:     class for standard call of bond
// Author:      Wang
// Created:     2004/08/17 
// RCS-ID:      $Id: callfixedshare.cpp,v 1.14 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file finance/bondlike/callfixedshare.cpp
   @todo replace the error code with specific error
 */

#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/callfixedshare.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/callfixedshare.h"
#include "ito33/xml/finance/bondlike/common.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_CALL_WRONG_START_END,
  ITO33_BONDLIKE_SOFTCALL_BADTRIGGER,
  ITO33_PEPSLIKE_CALLFIXEDSHARE_RATIO_NEGATIVE;

namespace ito33 
{

namespace finance
{


CallFixedShare::CallFixedShare(Date startDate, Date endDate, double dRatio)
                             : m_startDate(startDate), m_endDate(endDate),
                               m_dRatio(dRatio), m_dTriggerRate(0.)
{
  CHECK_COND(startDate <= endDate, ITO33_BONDLIKE_CALL_WRONG_START_END);

  CHECK_COND(dRatio > 0, ITO33_PEPSLIKE_CALLFIXEDSHARE_RATIO_NEGATIVE);
}

void CallFixedShare::SetTrigger(double dTriggerRate) 
{
  CHECK_COND(dTriggerRate > 0, ITO33_BONDLIKE_SOFTCALL_BADTRIGGER);

  m_dTriggerRate = dTriggerRate;
}

bool CallFixedShare::HasMakeWhole() const
{
  return m_makeWholeType != MakeWholeType_Max && IsSoft();
}

XML::Tag CallFixedShare::Dump(XML::Tag& tagParent) const
{

  XML::Tag tagMe(XML_TAG_BONDLIKE_CALLFIXEDSHARE_ROOT, tagParent);

  Call::DumpMe(tagMe);

  tagMe.Element(XML_TAG_BONDLIKE_STARTDATE)(m_startDate);

  tagMe.Element(XML_TAG_BONDLIKE_ENDDATE)(m_endDate);

  tagMe.Element(XML_TAG_BONDLIKE_RATIO)(m_dRatio);

  if (m_dTriggerRate > 0)
    tagMe.Element(XML_TAG_BONDLIKE_TRIGGERRATE)(m_dTriggerRate);

  return tagMe;
}


} // namespace finance

} // namespace ito33
