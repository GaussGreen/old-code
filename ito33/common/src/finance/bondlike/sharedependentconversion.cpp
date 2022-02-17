/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/sharedependentconversion.cpp
// Purpose:     class for share dependent conversion of a bond
// Author:      ITo 33
// Created:     2004/12/31
// RCS-ID:      $Id: sharedependentconversion.cpp,v 1.6 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/constants.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/sharedependentconversion.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/sharedependentconversion.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

extern const ito33::finance::BondError
    ITO33_BONDLIKE_COCO_BADTRIGGER,
    ITO33_BONDLIKE_COCO_TYPE,
    ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE,
    ITO33_BONDLIKE_COCO_CURRENTLYACTIVE,
    ITO33_BONDLIKE_NO_COCO,
    ITO33_AWCB_INVALID_FUNCTION_CALL_TO_GET_RESET_DATE;

extern const ito33::finance::BondError
    ITO33_BONDLIKE_CONVERSION_WRONG_START_END,
    ITO33_BONDLIKE_CONVERSION_BADRATIO,
    ITO33_BONDLIKE_SHAREDEPENDENT_BADRESET,
    ITO33_BONDLIKE_SHAREDEPENDENT_BADINCREMENTALSHAREFACTOR,
    ITO33_BONDLIKE_SHAREDEPENDENT_BADSTRIKE,
    ITO33_BONDLIKE_SHAREDEPENDENT_BADCAP;

namespace ito33
{

namespace finance
{


ShareDependentConversion::ShareDependentConversion
  (Date startDate, Date endDate, double dBaseRatio,
   double dIncrementalShareFactor)
 : Conversion(),
   m_StartDate(startDate), 
   m_EndDate(endDate),
   m_bHasResetDate(false),
   m_dBaseRatio(dBaseRatio),
   m_dIncrementalShareFactor(dIncrementalShareFactor),   
   m_dCapRatio(-1.),
   m_dFixedStrike(-1.),
   m_dTriggerRate(-1.), 
   m_coCoType(CoCoType_Max),
   m_dChangeInTriggerRate(0.),
   m_dExtremeTriggerRate(-1.),
   m_dCurrentRatio(0.)
{

  CHECK_COND(startDate <= endDate, 
             ITO33_BONDLIKE_CONVERSION_WRONG_START_END);

  CHECK_COND(dBaseRatio > 0,
             ITO33_BONDLIKE_CONVERSION_BADRATIO);

  CHECK_COND(dIncrementalShareFactor > 0,
             ITO33_BONDLIKE_SHAREDEPENDENT_BADINCREMENTALSHAREFACTOR);

}

void ShareDependentConversion::SetCapRatio(double dCapRatio)
{  
 
  CHECK_COND(dCapRatio >= m_dBaseRatio,
             ITO33_BONDLIKE_SHAREDEPENDENT_BADCAP);

  m_dCapRatio  = dCapRatio;

}

void ShareDependentConversion::SetCurrentRatio(double dCurrentRatio)
{
  CHECK_COND( dCurrentRatio > 0., ITO33_BONDLIKE_CONVERSION_BADRATIO);

  m_dCurrentRatio = dCurrentRatio;
}

void ShareDependentConversion::SetResetDate(Date resetDate)
{

  CHECK_COND(resetDate > m_StartDate,
             ITO33_BONDLIKE_SHAREDEPENDENT_BADRESET);

  CHECK_COND(resetDate < m_EndDate,
             ITO33_BONDLIKE_SHAREDEPENDENT_BADRESET);

  m_ResetDate     = resetDate;
  m_bHasResetDate = true;
}

void ShareDependentConversion::SetFixedStrike(double dFixedStrike)
{
  CHECK_COND(dFixedStrike >= 0.,
             ITO33_BONDLIKE_SHAREDEPENDENT_BADSTRIKE);

  m_dFixedStrike = dFixedStrike;
  
}

void ShareDependentConversion::SetCoCo
     ( 
       double dTrigger, CoCoType coCoType,
       double dChangeRate, double dExtremeTrigger,
       bool bIsLastTriggerConditionMet
     )
{    
  CHECK_COND(dTrigger > 0, ITO33_BONDLIKE_COCO_BADTRIGGER);

  CHECK_COND(coCoType < CoCoType_Max,
              ITO33_BONDLIKE_COCO_TYPE);

  CHECK_COND(dExtremeTrigger > 0,
              ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE);

  CHECK_COND( (dChangeRate >= 0 && dExtremeTrigger >= dTrigger) ||
              (dChangeRate < 0 && dExtremeTrigger <= dTrigger),
              ITO33_BONDLIKE_COCO_EXTREMETRIGGERRATE);

  CHECK_COND( bIsLastTriggerConditionMet == false || 
              coCoType != CoCoType_CheckAnyTimeAndConvertOnCheckDate,
              ITO33_BONDLIKE_COCO_CURRENTLYACTIVE);

  m_dTriggerRate           = dTrigger;
  m_coCoType               = coCoType;
  m_dChangeInTriggerRate   = dChangeRate;
  m_dExtremeTriggerRate    = dExtremeTrigger;
  m_bIsLastTriggerConditionMet = bIsLastTriggerConditionMet;
}

Date ShareDependentConversion::GetResetDate() const 
{ 
  CHECK_COND( HasResetDate(), 
    ITO33_AWCB_INVALID_FUNCTION_CALL_TO_GET_RESET_DATE);

  return m_ResetDate; 
}

double ShareDependentConversion::GetTrigger() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dTriggerRate; 
}

CoCoType ShareDependentConversion::GetCoCoType() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_coCoType; 
}

double ShareDependentConversion::GetChangeRate() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dChangeInTriggerRate; 
}

double ShareDependentConversion::GetExtremeTrigger() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_dExtremeTriggerRate; 
}

bool ShareDependentConversion::GetIsLastTriggerConditionMet() const 
{ 
  CHECK_COND(HasCoCo(), ITO33_BONDLIKE_NO_COCO);

  return m_bIsLastTriggerConditionMet; 
}

XML::Tag ShareDependentConversion::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_SHAREDEPENDENT_ROOT, tagParent);

  DumpMe(tagMe);

  tagMe.Element(XML_TAG_BONDLIKE_STARTDATE)(m_StartDate);
  tagMe.Element(XML_TAG_BONDLIKE_ENDDATE)(m_EndDate);
  
  if ( HasResetDate() )
  {
    tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_RESETDATE)(m_ResetDate);

    if ( m_dCurrentRatio > 0)   
      tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_CURRENT_RATIO)
      (m_dCurrentRatio);
  }
  
  tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_BASERATIO)(m_dBaseRatio);

  if ( m_dCapRatio >= 0. )
    tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_CAPRATIO)(m_dCapRatio);
  
  tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_SHAREFACTOR)
        (m_dIncrementalShareFactor);

  if (m_dFixedStrike >= 0.)
    tagMe.Element(XML_TAG_BONDLIKE_SHAREDEPENDENT_FIXEDSTRIKE)(m_dFixedStrike);
  
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

    if ( m_bIsLastTriggerConditionMet )
      tagMe.Element(XML_TAG_BONDLIKE_ISLASTTRIGGERMET)
                   (m_bIsLastTriggerConditionMet);
  } //end if

  return tagMe;
}


} // namespace finance

} // namespace ito33
