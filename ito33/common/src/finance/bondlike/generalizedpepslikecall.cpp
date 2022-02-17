/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/generalizedpepslikecall.cpp
// Purpose:     class for standard call of bond
// Created:     2004/08/17 
// RCS-ID:      $Id: generalizedpepslikecall.cpp,v 1.5 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
   @file finance/bondlike/callfixedshare.cpp
   @todo replace the error code with specific error
 */

#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/generalizedpepslikecall.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/xml/finance/bondlike/common.h"

extern const ito33::finance::BondError
  ITO33_BONDLIKE_CALL_WRONG_START_END,
  ITO33_GENERALIZEDPEPSLIKE_CALL_TRIGGER,
  ITO33_GENERALIZEDPEPSLIKE_CALL_TYPE;

namespace ito33 
{

namespace finance
{

GeneralizedPEPSLikeCall::GeneralizedPEPSLikeCall
                          (
                            Date startDate,
                            Date endDate,
                            double dTriggerRate,
                            GeneralizedPEPSLikeCallType type
                          )
                          : m_startDate(startDate),
                            m_endDate(endDate),
                            m_dTriggerRate(dTriggerRate),
                            m_type(type)
{
  CHECK_COND(startDate <= endDate, ITO33_BONDLIKE_CALL_WRONG_START_END);

  CHECK_COND(dTriggerRate >= 1., ITO33_GENERALIZEDPEPSLIKE_CALL_TRIGGER);

  CHECK_COND
  (
    type < GeneralizedPEPSLikeCallType_Max,
    ITO33_GENERALIZEDPEPSLIKE_CALL_TYPE
  );  
}

bool GeneralizedPEPSLikeCall::HasMakeWhole() const
{
  return m_makeWholeType < MakeWholeType_Max;
}


XML::Tag GeneralizedPEPSLikeCall::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_ROOT, tagParent);

  Call::DumpMe(tagMe);

  tagMe.Element(XML_TAG_BONDLIKE_STARTDATE)(m_startDate);

  tagMe.Element(XML_TAG_BONDLIKE_ENDDATE)(m_endDate);

  tagMe.Element(XML_TAG_BONDLIKE_TRIGGERRATE)(m_dTriggerRate);

  tagMe.Element(XML_TAG_BONDLIKE_GENERALIZED_PEPSLIKE_CALL_TYPE)
               (
                   GetNameFromEnumValue(
                     m_type,
                     SIZEOF(g_GeneralizedPEPSLikeCallTypes),
                     g_GeneralizedPEPSLikeCallTypes)
               );

  return tagMe;
}


} // namespace finance

} // namespace ito33
