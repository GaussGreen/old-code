/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/conversionschedulechedule.cpp
// Purpose:     class for standard conversion provision of a bond
// Author:      ZHANG Yunzhi
// Created:     2004 may 6
// RCS-ID:      $Id: conversionschedule.cpp,v 1.20 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/conversionperiod.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/bondlike/cocotype.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/conversionschedule.h"
#include "ito33/xml/finance/bondlike/common.h"
#include "ito33/xml/finance/bondlike/cocotype.h"

extern const ito33::finance::BondError
    ITO33_BONDLIKE_CONVERSIONS_NULL_PERIOD,
    ITO33_BONDLIKE_CONVERSIONS_WRONGSCHEDULE,
    ITO33_BONDLIKE_COCO_CURRENTLYACTIVE;

namespace ito33
{

namespace finance
{


void ConversionSchedule::AddConversionPeriod
                         ( const shared_ptr<ConversionPeriod>& pConversion )
{
  CHECK_PTR(pConversion, ITO33_BONDLIKE_CONVERSIONS_NULL_PERIOD);

  Elements::iterator 
    iterBefore = m_pConversions.end(),
    iterNext = m_pConversions.begin();

  for( ;
       iterNext != m_pConversions.end()
          && pConversion->GetEndDate() > (*iterNext)->GetStartDate();
       iterBefore = iterNext, iterNext++)
     ;

  // at this moment pConversion->endDate <= (*iterNext)->startDate or iterNext == end()
  // and we compare conversion to iterBefore to check schedule error
  if(
      // wrong schedule with the next one : equal
      iterNext != m_pConversions.end() && pConversion->GetStartDate() == (*iterNext)->GetEndDate() 
    ||
      // wrong schedule with the last one
      (iterBefore != m_pConversions.end() 
      && pConversion->GetStartDate() < (*iterBefore)->GetEndDate())
      // we have always endDate > (*iterBefore)->Start, thus equal is not possible
    )
    throw EXCEPTION_MSG
          ( 
            ITO33_BONDLIKE_CONVERSIONS_WRONGSCHEDULE,
            ITO33_BONDLIKE_CONVERSIONS_WRONGSCHEDULE.GetMessage()
          );
        
  m_pConversions.insert(iterNext, pConversion);
}


void ConversionSchedule::ValidateWith(Date valuationDate) const
{
  // If m_bIsLastTriggerConditionMet is true, verify that the period
  // containing the valuation date is a valid CoCo period.
  if ( m_bIsLastTriggerConditionMet == false || m_pConversions.empty() )
    return;  

  // Find the period containing the valuation date
  bool bIsValidCoCoPeriod = false;

  if ( valuationDate >= m_pConversions.front()->GetStartDate() &&
       valuationDate < m_pConversions.back()->GetEndDate() )
  {
    Elements::const_iterator iter;

    for ( iter = m_pConversions.begin(); iter != m_pConversions.end(); ++iter )
    {
      if ( !(*iter)->HasCoCo() )
        continue;

      if ( valuationDate >= (*iter)->GetStartDate() &&
           valuationDate < (*iter)->GetEndDate() )
      {
        CoCoType cocoType = (*iter)->GetCoCoType();

        if ( cocoType == CoCoType_CheckAnyTimeAndConvertAsOfCheckDate ||
             cocoType == CoCoType_CheckQuarterlyAndConvertAsOfCheckDate || 
             cocoType == CoCoType_CheckQuarterlyAndConvertDuringNextQuarter )
        {
          bIsValidCoCoPeriod = true;    
          break;
        }
      }
    } // loop over periods
  } // is valuation within a period

  CHECK_COND( bIsValidCoCoPeriod, ITO33_BONDLIKE_COCO_CURRENTLYACTIVE );
}


XML::Tag ConversionSchedule::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_ROOT, tagParent);

  DumpMe(tagMe);

  if (m_bIsLastTriggerConditionMet)
    tagMe.Element(XML_TAG_BONDLIKE_ISLASTTRIGGERMET)
                 (m_bIsLastTriggerConditionMet);

  XML::Tag tagPeriods(
             XML_TAG_BONDLIKE_CONVERSIONSCHEDULE_CONVERSIONPERIODS, tagMe);

  Elements::const_iterator iter;
  for(iter = m_pConversions.begin(); iter != m_pConversions.end(); ++iter)
    tagPeriods.Element(*(*iter));

  return tagMe;
}


} // namespace finance

} // namespace ito33
