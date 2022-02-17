/////////////////////////////////////////////////////////////////////////////
// Name:        sessiondata_xml.cpp
// Purpose:     Restore SessionData object from XML document
// Created:     2006/03/23
// RCS-ID:      $Id: sessiondata_xml.cpp,v 1.3 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/date.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/sessiondata.h"
#include "ito33/xml/finance/equity.h"
#include "ito33/xml/finance/ratedata.h"
#include "ito33/xml/finance/common.h"

using namespace ito33;

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


/**
    Exception thrown when the RateData node is not found in XML document.
 */
class MissingRateDataInSessionData : public Exception
{
public:
  /**
      Constructs the exception thrown when a RateData XML node is not found.

      This exception is thrown when RateData is not found in session data tag.

      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingRateDataInSessionData(const char *filename,
                               size_t line,
                               const char *function)
    : Exception(ITO33_BAD_DATA,
                "RateData not found in session data tag.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_RATEDATA_IN_SESSIONDATA_EXCEPTION                        \
  throw MissingRateDataInSessionData(__FILE__, __LINE__, __FUNCTION__)


shared_ptr<finance::SessionData> GetSessionDataFromNode(const xml::node& node)
{

  xml::node::const_iterator pNodeFound;

  // valuation date
  Date valuationDate = 
    GetDateFromName(node, XML_TAG_SESSIONDATA_VALUATION_DATE);

  // Equity 
  pNodeFound = node.find(XML_TAG_EQUITY_ROOT);

  if ( pNodeFound == node.end() )
    THROW_MISSING_EQUITY_IN_SESSIONDATA_EXCEPTION;

  shared_ptr<finance::Equity> pEquity = GetEquityFromNode(*pNodeFound);

  // RateData
  pNodeFound = node.find(XML_TAG_RATEDATA_ROOT);

  if ( pNodeFound == node.end() )
    THROW_MISSING_RATEDATA_IN_SESSIONDATA_EXCEPTION;

  shared_ptr<finance::RateData> pRateData = GetRateDataFromNode(*pNodeFound);

  // Construct and return
  shared_ptr<finance::SessionData> 
    pSessionData( new finance::SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}

} // namespace XML

} // namespace ito33

