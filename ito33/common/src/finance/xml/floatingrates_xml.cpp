/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/floatingrates_xml.cpp
// Purpose:     Restore FloatingRates object from XML
// Created:     2005/09/09
// RCS-ID:      $Id: floatingrates_xml.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004-2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"
#include "ito33/date.h"
#include "ito33/useexception.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/frequency.h"
#include "ito33/finance/floatingrates.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/read_frequency.h"
#include "ito33/xml/finance/read_daycountconvention.h"
#include "ito33/xml/finance/lastpaymenttype.h"
#include "ito33/xml/finance/floatingrates.h"

extern const ito33::Error ITO33_UNEXPECTED;

namespace ito33
{

namespace XML
{

shared_ptr<finance::FloatingRates>
GetFloatingRatesFromNode(const xml::node& node)
{
  bool bHasUnknownPayments = false;
  finance::LastPaymentType lastPaymentType = finance::LastPaymentType_Max;

  Date first, last;

  Date startOfAccruedDate = GetDateFromName
                            (node, XML_TAG_FLOATINGRATES_STARTOFACCRUEDDATE);

  xml::node::const_iterator iterSub;

  if ( (iterSub = node.find(XML_TAG_FLOATINGRATES_FIRSTPAYMENTDATE) )
        != node.end() )
  {
    bHasUnknownPayments = true;

    first = GetDateFromNode(*iterSub);
        
    last = GetDateFromName(node, XML_TAG_FLOATINGRATES_LASTDATE);
    
    if ( (iterSub = node.find(XML_TAG_LASTPAYMENTTYPE) ) != node.end() )
    {
      lastPaymentType = GetEnumFromNode
                        (
                          *iterSub,
                          SIZEOF(g_LastPaymentType),
                          g_LastPaymentType
                        );
    }
  }
      
  finance::Frequency
    freq = GetFrequencyFromName(node, XML_TAG_PAYMENTFREQUENCY);

  double dMargin = GetDoubleFromName
                   (node, XML_TAG_FLOATINGRATES_MARGIN);

  shared_ptr<finance::FloatingRates> pFloatingRates;

  if ( !bHasUnknownPayments )
    pFloatingRates = make_ptr( new finance::FloatingRates
                                   ( startOfAccruedDate, freq ) );
  else if ( IsValid(lastPaymentType) )
    pFloatingRates = make_ptr( new finance::FloatingRates 
                                   ( dMargin, startOfAccruedDate, first, last, 
                                     freq, lastPaymentType ) );
  else
    pFloatingRates = make_ptr( new finance::FloatingRates 
                                   ( dMargin, startOfAccruedDate, first, last,
                                     freq ) );

  
  Date::DayCountConvention
    dcm = GetDayCountConventionFromName(node, XML_TAG_DAYCOUNTCONVENTION);

  pFloatingRates->SetDayCountConvention(dcm);
  
  double 
    dMultiplier = GetDoubleFromName(node, XML_TAG_FLOATINGRATES_MULTIPLIER);

  pFloatingRates->SetMultiplier(dMultiplier);
  
  double dCap = GetDoubleFromName(node, XML_TAG_FLOATINGRATES_CAP);

  pFloatingRates->SetCap(dCap);

  double dFloor = GetDoubleFromName
                  (node, XML_TAG_FLOATINGRATES_FLOOR);

  pFloatingRates->SetFloor(dFloor);

  int iFixingDelay = GetLongFromName
                     (node, XML_TAG_FLOATINGRATES_FIXINGDELAY);

  pFloatingRates->SetFixingDelay(iFixingDelay);

  // The known payment stream
  std::vector<Date> pDates;
  std::vector<double> pdValues;

  if (     (iterSub = node.find(XML_TAG_FLOATINGRATES_KNOWNPAYMENTS) )
        != node.end() )
  {
    for ( xml::node::const_iterator iter = iterSub->begin(); 
          iter != iterSub->end(); 
          iter++)
    {
      if( strcmp(iter->get_name(), XML_TAG_FLOATINGRATES_KNOWNPAYMENT) == 0)
      {
        pDates.push_back(
          GetDateFromName(*iter, XML_TAG_FINANCE_DATE));
        pdValues.push_back(
          GetDoubleFromName(*iter, XML_TAG_FINANCE_RATE));
      }
    }

    if ( !pDates.empty() )
      pFloatingRates->SetKnownPaymentStream(pDates, pdValues);
  }

  return pFloatingRates;
}

} // namespace XML

} // namespace ito33
