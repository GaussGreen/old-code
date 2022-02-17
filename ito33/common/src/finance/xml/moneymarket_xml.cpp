/////////////////////////////////////////////////////////////////////////////
// Name:        session_xml.cpp
// Purpose:     Restore MoneyMarket object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-04
// RCS-ID:      $Id: moneymarket_xml.cpp,v 1.8 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <string>
#include "ito33/afterstd.h"

#include "ito33/finance/moneymarket.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/moneymarket.h"
#include "ito33/xml/finance/yieldcurve.h"

namespace ito33
{

namespace XML
{

/**
    Exception thrown when the YieldCurve node is not found in XML document.
 */
class MissingYieldCurveInMoneyMarket : public Exception
{
public:
  /**
      Constructs the exception thrown when a YieldCurve XML node is not found.

      This exception is thrown when YieldCurve is not found in money market
      data tag.

      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingYieldCurveInMoneyMarket(const char *filename,
                                 size_t line,
                                 const char *function)
    : Exception(ITO33_BAD_DATA,
                "YieldCurve not found in money market data tag.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_YIELDCURVE_IN_MONEYMARKET_EXCEPTION   \
  throw MissingYieldCurveInMoneyMarket(__FILE__, __LINE__, __FUNCTION__)


/**
    Exception thrown when the Numuraire node is not found in XML document.
 */
class MissingNumeraireInMoneyMarket : public Exception
{
public:
  /**
      Constructs the exception thrown when a Numeraire XML node is not found.

      This exception is thrown when Numeraire is not found in money market
      data tag.

      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingNumeraireInMoneyMarket(const char *filename,
                                size_t line,
                                const char *function)
    : Exception(ITO33_BAD_DATA,
                "Numeraire not found in money market data tag.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_NUMERAIRE_IN_MONEYMARKET_EXCEPTION   \
  throw MissingNumeraireInMoneyMarket(__FILE__, __LINE__, __FUNCTION__)


shared_ptr<finance::MoneyMarket> GetMoneyMarketFromNode(const xml::node& node)
{

  xml::node::const_iterator pNodeFound;

  // Yield curve
  pNodeFound = node.find(XML_TAG_MONEYMARKET_YIELDCURVE);
  
  if (pNodeFound == node.end())
    THROW_MISSING_YIELDCURVE_IN_MONEYMARKET_EXCEPTION;

  shared_ptr<finance::YieldCurve> 
    pYieldCurve = GetYieldCurveFromNode(*pNodeFound);

  // Numeraire (currency)
  pNodeFound = node.find(XML_TAG_MONEYMARKET_NUMERAIRE);

  if (pNodeFound == node.end())
    THROW_MISSING_NUMERAIRE_IN_MONEYMARKET_EXCEPTION;

  shared_ptr<finance::Numeraire> 
    pNumeraire(new finance::Numeraire(pNodeFound->get_content() ) );    

  // Construct and return
  shared_ptr<finance::MoneyMarket> pMoneyMarket(
    new finance::MoneyMarket(pNumeraire, pYieldCurve) );

  return pMoneyMarket;
}

} // namespace XML

} // namespace ito33
