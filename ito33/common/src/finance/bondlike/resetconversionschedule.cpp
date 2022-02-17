/////////////////////////////////////////////////////////////////////////////
// Name:        finance/bondlike/resetconversionschedule.cpp
// Purpose:     Reset schedule for a resettable convertible bond
// Author:      Yann and David
// Created:     2004/10/19
// RCS-ID:      $Id: resetconversionschedule.cpp,v 1.8 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file finance/bondlike/resetconversionschedule.cpp
*/

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/conversionpricereset.h"
#include "ito33/finance/bondlike/resetconversionschedule.h"
#include "ito33/finance/bondlike/resetflooredby.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/resetconversionschedule.h"
#include "ito33/xml/finance/bondlike/common.h"

extern const ito33::finance::BondError 
    ITO33_BONDLIKE_INVALID_RESET_PERIOD,
    ITO33_BONDLIKE_INVALID_RESET_INITIAL_CONV_PRICE,
    ITO33_BONDLIKE_INVALID_RESET_CURRENT_CONV_PRICE,
    ITO33_BONDLIKE_INVALID_RESET_DATE;

namespace ito33
{

namespace finance
{

ResetConversionSchedule::ResetConversionSchedule
    (Date startDate, Date endDate,
     double dInitialConvPrice, double dCurrentConvPrice, 
     ResetFlooredBy flooredBy)
   : Conversion(),
     m_StartDate(startDate),
     m_EndDate(endDate),             
     m_dInitialConvPrice(dInitialConvPrice),
     m_dCurrentConvPrice(dCurrentConvPrice),
     m_ResetFlooredBy(flooredBy),
     m_dCash(0.)
{
  CHECK_COND
  (
    m_StartDate < m_EndDate,
    ITO33_BONDLIKE_INVALID_RESET_PERIOD
  );

  CHECK_COND
  (
    m_dInitialConvPrice > 0.0,
    ITO33_BONDLIKE_INVALID_RESET_INITIAL_CONV_PRICE
  );

  CHECK_COND
  (
    m_dCurrentConvPrice > 0.0,
    ITO33_BONDLIKE_INVALID_RESET_CURRENT_CONV_PRICE
  );

} 

void ResetConversionSchedule::AddConversionPriceReset
     ( const shared_ptr<ConversionPriceReset>& pConversionPriceReset )
{

  // Basic checks.  Make sure new date is within the conversion period
  CHECK_COND
  (
    pConversionPriceReset->GetDate() >= m_StartDate,
    ITO33_BONDLIKE_INVALID_RESET_DATE
  );

  CHECK_COND
  (
    pConversionPriceReset->GetDate() <= m_EndDate,
    ITO33_BONDLIKE_INVALID_RESET_DATE
  );

  // search through exisiting list to find insertion point
  Elements::iterator iterNext; 

  for( iterNext = m_ResetTermList.begin();
       iterNext != m_ResetTermList.end()
          && pConversionPriceReset->GetDate() > (*iterNext)->GetDate();
       iterNext++)
     ;

  // At this moment pConversionPriceReset->date <= (*iterNext)->date or iterNext== end()
  // Check for schedule error
  if ( iterNext != m_ResetTermList.end() )
  {
    CHECK_COND
    (
      pConversionPriceReset->GetDate() != (*iterNext)->GetDate(),
      ITO33_BONDLIKE_INVALID_RESET_DATE
    );
  }
  
  // valid entry, so insert
  m_ResetTermList.insert(iterNext, pConversionPriceReset);
}


void 
ResetConversionSchedule::SetInitialConversionPrice(double dInitialConvPrice)
{
  m_dInitialConvPrice = dInitialConvPrice;

  CHECK_COND
  (
    m_dInitialConvPrice > 0.0,
    ITO33_BONDLIKE_INVALID_RESET_INITIAL_CONV_PRICE
  );

}

void 
ResetConversionSchedule::SetCurrentConversionPrice(double dCurrentConvPrice)
{
  m_dCurrentConvPrice = dCurrentConvPrice;

  CHECK_COND
  (
    m_dCurrentConvPrice > 0.0,
    ITO33_BONDLIKE_INVALID_RESET_CURRENT_CONV_PRICE
  );

}


XML::Tag ResetConversionSchedule::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMe(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_ROOT, tagParent);

  DumpMe(tagMe);

  tagMe.Element(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_START)(m_StartDate);
  tagMe.Element(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_END)(m_EndDate);
  tagMe.Element(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_INITIAL)(m_dInitialConvPrice);
  tagMe.Element(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CURRENT)(m_dCurrentConvPrice);
  tagMe.Element(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CASH)(m_dCash);

  tagMe.Element (XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_FLOOREDBY)
                    (
                      GetNameFromEnumValue(
                        m_ResetFlooredBy, 2, g_resetFlooredBy)
                    );

  XML::Tag tagTerms(XML_TAG_BONDLIKE_RESETCONVERSIONSCHEDULE_CONVPRICERESET, tagMe);

  Elements::const_iterator iter;
  for(iter = m_ResetTermList.begin(); iter != m_ResetTermList.end(); ++iter)
    tagTerms.Element(*(*iter));

  return tagMe;
}


} // namespace finance

} // namespace ito33
