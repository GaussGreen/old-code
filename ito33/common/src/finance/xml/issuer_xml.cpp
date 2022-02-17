/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/xml/issuer_xml.cpp
// Purpose:     Restore Issuer object from XML document
// Author:      Nabil
// Created:     2005/03/17
// RCS-ID:      $Id: issuer_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005-2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/finance/issuer.h"
#include "ito33/finance/moneymarket.h"

#include "ito33/xml/read.h"
#include "ito33/xml/read_vector.h"
#include "ito33/xml/finance/issuer.h"
#include "ito33/xml/finance/moneymarket.h"


namespace ito33
{
namespace XML
{
 
  
/**
    Exception thrown when the expected node is not found in XML document.
 */
class MissingMoneyMarketInIssuer : public Exception
{
public:
  /**
      Constructs the exception thrown when an expected XML node is not found.

      This exception is thrown when money_market is not found in session tag.

      @param node the parent node
      @param name the name of child note we didn't find
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingMoneyMarketInIssuer(const char *filename,
                       size_t line,
                       const char *function)
    : Exception(ITO33_BAD_DATA,
                "Money market data is missed.",
                filename,
                line,
                function)
  {
  }
};

#define THROW_MISSING_MONEYMARKET_IN_ISSUER_EXCEPTION                        \
  throw MissingMoneyMarketInIssuer(__FILE__, __LINE__, __FUNCTION__)

shared_ptr<finance::Issuer> GetIssuerFromNode(const xml::node& node)
{ 
  shared_ptr<finance::Issuer> pIssuer(new finance::Issuer);
  
  // Read the fiscal year start date
  xml::node::const_iterator
    pNodeRoot( node.find(XML_TAG_ISSUER_FISCALYEARSTART_DATE) );
  if ( pNodeRoot != node.end() )
  {
    Date fiscalYearStartDate
      = GetDateFromName(node, XML_TAG_ISSUER_FISCALYEARSTART_DATE);
      
    pIssuer->SetFiscalYearStartDate(fiscalYearStartDate);
  }

  // Read the default intensity
  pNodeRoot = node.find(XML_TAG_ISSUER_DEFAULTINTENSITY);
  if ( pNodeRoot != node.end() )
  {
     pIssuer->SetDefaultIntensity
              (
                GetVectorFromNode<Date>  // time array
                          (*pNodeRoot,
                           XML_TAG_FINANCE_DATES,
                           XML_TAG_FINANCE_DATE),
                GetVectorFromNode<double>  // value array
                          (*pNodeRoot,
                           XML_TAG_FINANCE_VALUES,
                           XML_TAG_FINANCE_VALUE)
              );
  }

  return pIssuer;
}

} // namespace XML

} // namespace ito33
