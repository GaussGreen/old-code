/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/sessiondata.cpp
// Purpose:     Implementation of the SessionData class
// Created:     2006/03/16
// RCS-ID:      $Id: sessiondata.cpp,v 1.4 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/sessiondata.h"
#include "ito33/xml/finance/dividends.h"

extern const ito33::finance::Error
  ITO33_SESSIONDATA_INVALID_RATEDATA,
  ITO33_SESSIONDATA_INVALID_EQUITY;

extern const ito33::Error
  ITO33_BAD_DATE;

namespace ito33
{

namespace finance
{

SessionData::SessionData(const shared_ptr<RateData>& pRateData,
                         const shared_ptr<Equity>& pEquity,
                         Date valuationDate)
{

  CHECK_COND( pRateData, ITO33_SESSIONDATA_INVALID_RATEDATA );

  CHECK_COND( pEquity, ITO33_SESSIONDATA_INVALID_EQUITY );

  CHECK_COND( valuationDate.IsValid(), ITO33_BAD_DATE );

  m_pRateData = pRateData;
  
  m_pEquity = pEquity; 

  m_dValuationDate = valuationDate;
}


void SessionData::Dump(XML::Tag& tagParent) const
{
  // Create tag, dump the valuation date, then let the rate data and
  // equity objects dump themselves
  XML::Tag tagFinance(XML_TAG_SESSIONDATA, tagParent);

  tagFinance.Element(XML_TAG_SESSIONDATA_VALUATION_DATE)(m_dValuationDate);

  m_pRateData->Dump(tagFinance);

  m_pEquity->Dump(tagFinance);

}


} // end namespace finance

} // end namespace ito33 
