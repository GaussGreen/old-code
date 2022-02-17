/////////////////////////////////////////////////////////////////////////////
// Name:        ratedata_xml.cpp
// Purpose:     Restore RateData object from XML document
// Created:     2006/03/23
// RCS-ID:      $Id: ratedata_xml.cpp,v 1.3 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/ratedata.h"
#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/moneymarket.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/ratedata.h"
#include "ito33/xml/finance/spotfxrates.h"
#include "ito33/xml/finance/moneymarket.h"


namespace ito33
{

namespace XML
{
 
  
/**
    Exception thrown when the Equity node is not found in XML document.
 */
class MissingEquityInSessionData : public Exception
{
public:
  /**
      Constructs the exception thrown when an Equity XML node is not found.

      This exception is thrown when equity is not found in session data tag.

      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingEquityInSessionData(const char *filename,
                             size_t line,
                             const char *function)
    : Exception(ITO33_BAD_DATA,
                "Equity data not found in session data tag.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_EQUITY_IN_SESSIONDATA_EXCEPTION                        \
  throw MissingEquityInSessionData(__FILE__, __LINE__, __FUNCTION__)



shared_ptr<finance::RateData> GetRateDataFromNode(const xml::node& node)
{

  xml::node::const_iterator pNodeFound;

  // RateData object is initially empty. Construct and fill with optional data
  shared_ptr<finance::RateData> pRateData( new finance::RateData );

  // SpotFXRates is optional (at least, not in ctr)
  pNodeFound = node.find(XML_TAG_SPOT_FX_RATES_ROOT);

  if ( pNodeFound != node.end() )
  {
    shared_ptr<finance::SpotFXRates> 
      pSpotFXRates = GetSpotFXRatesFromNode(*pNodeFound);

    pRateData->SetSpotFXRates(pSpotFXRates);
  }

  // Read the money markets (yield curves)
  xml::node::const_iterator iter;
  for ( iter = node.begin(); iter != node.end(); ++iter )
  {
    if ( strcmp(iter->get_name(), XML_TAG_MONEYMARKET_ROOT) == 0 )
    {
      shared_ptr<finance::MoneyMarket> 
        pMoneyMarket = GetMoneyMarketFromNode(*iter);

      pRateData->SetYieldCurve( pMoneyMarket->GetNumeraire(),
                                pMoneyMarket->GetYieldCurve() );
    }

  } // loop over nodes

  return pRateData;
}

} // namespace XML

} // namespace ito33

