/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/varianceswap_xml.cpp
// Purpose:     Restore variance swap object from XML document
// Created:     2006/03/08
// RCS-ID:      $Id: varianceswap_xml.cpp,v 1.14 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2006 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/varianceswapterms.h"
#include "ito33/finance/varianceswap.h"
#include "ito33/finance/gammavarianceswap.h"
#include "ito33/finance/conditionalvarianceswap.h"
#include "ito33/finance/optionvarianceswap.h"
#include "ito33/finance/varianceswaption.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/varianceswap.h"
#include "ito33/xml/finance/sessiondata.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/returntype.h"
#include "ito33/xml/finance/swaptype.h"
#include "ito33/xml/finance/swappayofftype.h"
#include "ito33/xml/finance/optiontype.h"

using namespace ito33;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(varianceswap_xml);

namespace ito33
{

namespace XML
{

shared_ptr<finance::VarianceSwapTerms>
GetVarianceSwapTermsFromNode(const xml::node& node)
{
  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);
  Date startOfSamplingPeriod = 
    GetDateFromName(node, XML_TAG_VARIANCESWAP_STARTOFPERIOD);
  size_t nNbSamplingReturns = 
    GetLongFromName(node, XML_TAG_VARIANCESWAP_NBSAMPLINGRETURNS);

  finance::SwapType swapType = GetEnumFromName( node,
                                       XML_TAG_VARIANCESWAP_SWAP_TYPE,
                                       SIZEOF(g_swapTypes),
                                       g_swapTypes );

  shared_ptr<finance::VarianceSwapTerms>
    pTerms( new finance::VarianceSwapTerms
                (maturityDate, swapType, 
                 startOfSamplingPeriod, nNbSamplingReturns) );

  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAP_RETURN_TYPE);

  if ( pNodeFound != node.end() )
  {
    finance::ReturnType 
      returnType = GetEnumFromNode(*pNodeFound, 
                                   SIZEOF(g_returnTypes), g_returnTypes);

    pTerms->SetReturnType(returnType);
  }

  // Get the number of sampling returns in a year, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_ANNUALRETURNFREQUENCY);
  if ( pNodeFound != node.end() )
     pTerms->SetAnnualReturnFrequency(GetLongFromNode(*pNodeFound));

  // Get the cap, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_CAP_MULTIPLIER);
  if ( pNodeFound != node.end() )
    pTerms->SetCapMultiplier( GetDoubleFromNode(*pNodeFound) );

  // Get the corridor barriers, if they exist
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_UP_CORRIDOR);
  if ( pNodeFound != node.end() )
    pTerms->SetUpCorridorBarrier( GetDoubleFromNode(*pNodeFound) );

  pNodeFound = node.find(XML_TAG_VARIANCESWAP_DOWN_CORRIDOR);
  if ( pNodeFound != node.end() )
    pTerms->SetDownCorridorBarrier( GetDoubleFromNode(*pNodeFound) );

  return pTerms;
}

} // namespace XML

} // namespace ito33

/**
    Restores a variance swap object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new variance swap object to be deleted by the caller
 */
static
finance::Derivative *ReadVarianceSwap(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  // Read all mandatory parameters
  shared_ptr<finance::VarianceSwapTerms> 
    pTerms( GetVarianceSwapTermsFromNode(node) );

  double dStrike = 
    GetDoubleFromName(node, XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE);

  // Create the variance swap
  finance::VarianceSwap *
    pVarianceSwap = new finance::VarianceSwap(pTerms, dStrike);
   
  GetOptionalDerivativeDataFromNode(node, *pVarianceSwap);

  GetMarketPrice(node, *pVarianceSwap);

  // Get the current volatility, if it exists
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY);

  double dCurrentVolatility = -1.0;
  if ( pNodeFound != node.end() )
    dCurrentVolatility = GetDoubleFromNode(*pNodeFound);

  // Get the number of samples used, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED);

  size_t nNbSamplesUsed = 0;
  if ( pNodeFound != node.end() )
    nNbSamplesUsed = GetLongFromNode(*pNodeFound);  
  
  // Let the finance object throw an exception if the current volatility
  // and number of samples used are not consistent
  if (dCurrentVolatility >= 0.0 || nNbSamplesUsed > 0)
    pVarianceSwap->SetCurrentValues( dCurrentVolatility, nNbSamplesUsed );

  return pVarianceSwap; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_VARIANCESWAP_ROOT, VarianceSwap);

/**
    Restores a gamma variance swap object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new gamma variance swap object to be deleted by the caller
 */
static
finance::Derivative *ReadGammaVarianceSwap(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  // Read all mandatory parameters
  shared_ptr<finance::VarianceSwapTerms> 
    pTerms( GetVarianceSwapTermsFromNode(node) );

  double dStrike = 
    GetDoubleFromName(node, XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE);

  // Create the gamma variance swap
  finance::GammaVarianceSwap*
    pVarianceSwap = new finance::GammaVarianceSwap(pTerms, dStrike);
   
  GetOptionalDerivativeDataFromNode(node, *pVarianceSwap);

  GetMarketPrice(node, *pVarianceSwap);

  // Get the current volatility, if it exists
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY);

  double dCurrentVolatility = -1.0;
  if ( pNodeFound != node.end() )
    dCurrentVolatility = GetDoubleFromNode(*pNodeFound);

  // Get the number of samples used, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED);

  size_t nNbSamplesUsed = 0;
  if ( pNodeFound != node.end() )
    nNbSamplesUsed = GetLongFromNode(*pNodeFound);  
  
  // Get the share price at the start of the sampling period, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_START_SHARE_PRICE);

  double dStartSharePrice = 0.0;
  if ( pNodeFound != node.end() )
    dStartSharePrice = GetDoubleFromNode(*pNodeFound); 

  // Let the finance object throw an exception if the current values
  // are inconsistent
  if (dCurrentVolatility >= 0.0 || nNbSamplesUsed > 0 || dStartSharePrice > 0)
    pVarianceSwap->SetCurrentValues( dCurrentVolatility, 
                                     nNbSamplesUsed, 
                                     dStartSharePrice );

  return pVarianceSwap; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_GAMMAVARIANCESWAP_ROOT, 
                               GammaVarianceSwap);

/**
    Restores a conditional variance swap object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new gamma variance swap object to be deleted by the caller
 */
static
finance::Derivative *ReadConditionalVarianceSwap(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  // Read all mandatory parameters
  shared_ptr<finance::VarianceSwapTerms> 
    pTerms( GetVarianceSwapTermsFromNode(node) );

  double dStrike = 
    GetDoubleFromName(node, XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE);

  // Create the gamma variance swap
  finance::ConditionalVarianceSwap*
    pVarianceSwap = new finance::ConditionalVarianceSwap(pTerms, dStrike);
   
  GetOptionalDerivativeDataFromNode(node, *pVarianceSwap);

  GetMarketPrice(node, *pVarianceSwap);

  // Get the current volatility, if it exists
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY);

  double dCurrentVolatility = -1.0;
  if ( pNodeFound != node.end() )
    dCurrentVolatility = GetDoubleFromNode(*pNodeFound);

  // Get the number of samples used, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED);

  size_t nNbSamplesUsed = 0;
  if ( pNodeFound != node.end() )
    nNbSamplesUsed = GetLongFromNode(*pNodeFound);  

  // Get the current conditional count, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_CURRENT_COUNT);

  int iCurrentCount = -1;
  if ( pNodeFound != node.end() )
    iCurrentCount = GetLongFromNode(*pNodeFound); 

  // Let the finance object throw an exception if the current values
  // are inconsistent
  if (dCurrentVolatility >= 0.0 || nNbSamplesUsed > 0 || iCurrentCount > -1)
    pVarianceSwap->SetCurrentValues( dCurrentVolatility, 
                                     nNbSamplesUsed, 
                                     iCurrentCount );

  return pVarianceSwap; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_CONDITIONALVARIANCESWAP_ROOT, 
                               ConditionalVarianceSwap);

/**
    Restores an option variance swap object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new gamma variance swap object to be deleted by the caller
 */
static
finance::Derivative *ReadOptionVarianceSwap(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  // Read all mandatory parameters
  shared_ptr<finance::VarianceSwapTerms> 
    pTerms( GetVarianceSwapTermsFromNode(node) );

  double dStrike = 
    GetDoubleFromName(node, XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE);

  finance::OptionType optionType = GetEnumFromName
                                   (
                                      node,
                                      XML_TAG_OPTION_TYPE,
                                      SIZEOF(g_optionTypes),
                                      g_optionTypes
                                   );

  // Create the gamma variance swap
  finance::OptionVarianceSwap*
    pVarianceSwap = new finance::OptionVarianceSwap(pTerms, dStrike, optionType);
   
  GetOptionalDerivativeDataFromNode(node, *pVarianceSwap);

  GetMarketPrice(node, *pVarianceSwap);

  // Get the current volatility, if it exists
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY);

  double dCurrentVolatility = -1.0;
  if ( pNodeFound != node.end() )
    dCurrentVolatility = GetDoubleFromNode(*pNodeFound);

  // Get the number of samples used, if it exists
  pNodeFound = node.find(XML_TAG_VARIANCESWAP_NB_SAMPLES_USED);

  size_t nNbSamplesUsed = 0;
  if ( pNodeFound != node.end() )
    nNbSamplesUsed = GetLongFromNode(*pNodeFound);  
  
  // Let the finance object throw an exception if the current values
  // are inconsistent
  if (dCurrentVolatility >= 0.0 || nNbSamplesUsed > 0 )
    pVarianceSwap->SetCurrentValues( dCurrentVolatility, 
                                     nNbSamplesUsed );

  return pVarianceSwap; 
}


ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_OPTIONVARIANCESWAP_ROOT, 
                               OptionVarianceSwap);

/**
    Restores a variance swaption object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new variance swaption object to be deleted by the caller
 */
static
finance::Derivative *ReadVarianceSwaption(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  // Read all mandatory parameters
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_VARIANCESWAPTION_VSTERMS);

  if ( pNodeFound == node.end() )
    throw;

  shared_ptr<finance::VarianceSwapTerms> 
    pTerms( GetVarianceSwapTermsFromNode(*pNodeFound) );

  finance::OptionType optionType = GetEnumFromName
                                   (
                                      node,
                                      XML_TAG_OPTION_TYPE,
                                      SIZEOF(g_optionTypes),
                                      g_optionTypes
                                   );

  double dStrike = GetDoubleFromName(node, XML_TAG_FINANCE_STRIKE);

  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  // Create the variance swaption
  finance::VarianceSwaption *
    pVarianceSwaption = new finance::VarianceSwaption
                            (pTerms, optionType, dStrike, maturityDate);
   
  GetOptionalDerivativeDataFromNode(node, *pVarianceSwaption);

  GetMarketPrice(node, *pVarianceSwaption);
  
  return pVarianceSwaption; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_VARIANCESWAPTION_ROOT, VarianceSwaption);


