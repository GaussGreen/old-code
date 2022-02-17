/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/conversionperiod.cpp
// Purpose:     class for standard call of bond
// Created:     2004/08/17 
// RCS-ID:      $Id: conversionperiod.cpp,v 1.12 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/conversionperiod.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/conversionperiod.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

extern const ito33::finance::BondError 
    ITO33_BONDLIKE_CONVERSION_WRONG_START_END,
    ITO33_BONDLIKE_CONVERSION_BADRATIO,
    ITO33_BONDLIKE_COCO_BADTRIGGER,
    ITO33_BONDLIKE_COCO_TYPE,
    ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE,
    ITO33_BONDLIKE_NO_COCO;

namespace ito33
{

namespace finance
{


ConversionPeriod::ConversionPeriod(Date startDate, Date endDate, double dRatio)
    : m_startDate(startDate), m_endDate(endDate), m_dRatio(dRatio),
      m_dCash(0.),
      m_dTriggerRate(-1), 
      m_coCoType(CoCoType_Max),
      m_dChangeInTriggerRate(0.),
      m_dExtremeTriggerRate(-1.)

{
  CHECK_COND(startDate <= endDate, 
              ITO33_BONDLIKE_CONVERSION_WRONG_START_END);

  CHECK_COND(dRatio > 0,
              ITO33_BONDLIKE_CONVERSION_BADRATIO);
}

void ConversionPeriod::SetCoCo
     ( 
       double dTriggerRate, CoCoType coCoType,
       double dChangeRate, double dExtremeTriggerRate
     )
{    
  CHECK_COND(dTriggerRate > 0, ITO33_BONDLIKE_COCO_BADTRIGGER);

  CHECK_COND(coCoType < CoCoType_Max,
              ITO33_BONDLIKE_COCO_TYPE);

  CHECK_COND(dExtremeTriggerRate > 0,
              ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE);

  CHECK_COND( (dChangeRate >= 0 && dExtremeTriggerRate >= dTriggerRate) ||
              (dChangeRate < 0 && dExtremeTriggerRate <= dTriggerRate),
              ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE);

  m_dTriggerRate           = dTriggerRate;
  m_coCoType               = coCoType;
  m_dChangeInTriggerRate   = dChangeRate;
  m_dExtremeTriggerRate    = dExtremeTriggerRate;
}

double ConversionPeriod::GetTrigger() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dTriggerRate; 
}

CoCoType ConversionPeriod::GetCoCoType() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_coCoType; 
}

  double ConversionPeriod::GetChangeRate() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dChangeInTriggerRate; 
}

double ConversionPeriod::GetExtremeTrigger() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dExtremeTriggerRate; 
}

XML::Tag ConversionPeriod::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_CONVERSIONPERIOD_ROOT, tagParent);

  tagMe.Element(XML_TAG_BONDLIKE_STARTDATE)(m_startDate);

  tagMe.Element(XML_TAG_BONDLIKE_ENDDATE)(m_endDate);

  tagMe.Element(XML_TAG_BONDLIKE_RATIO)(m_dRatio);

  tagMe.Element(XML_TAG_BONDLIKE_CONVERSIONPERIOD_CASH)(m_dCash);

  if (m_dTriggerRate > 0)
  {
    tagMe.Element(XML_TAG_BONDLIKE_TRIGGERRATE)(m_dTriggerRate);

    tagMe.Element(XML_TAG_BONDLIKE_CHANGERATE)
                    (m_dChangeInTriggerRate);
    tagMe.Element(XML_TAG_BONDLIKE_EXTREMETRIGGERRATE)
                    (m_dExtremeTriggerRate);

    tagMe.Element(XML_TAG_BONDLIKE_COCOTYPE)
                 (
                   GetNameFromEnumValue(
                     m_coCoType, SIZEOF(g_coCoTypes), g_coCoTypes)
                 );

  }

  return tagMe;
}


} // namespace finance

} // namespace ito33
