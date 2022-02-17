/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/option_xml.cpp
// Purpose:     Restore Option object from XML document
// Author:      Vadim Zeitlin
// Created:     2004-05-08
// RCS-ID:      $Id: option_xml.cpp,v 1.16 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 - 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/option.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/optiontype.h"
#include "ito33/xml/finance/exercisetype.h"

using namespace ito33;
using namespace ito33::XML;
using finance::Option;

ITO33_FORCE_LINK_THIS_MODULE(option_xml);

/**
    Restore an Option object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::Derivative *ReadOption(const xml::node *pNode)
{
  const xml::node& node = *pNode;

  double dStrike = GetDoubleFromName(node, XML_TAG_FINANCE_STRIKE);
  Date expiryDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  finance::OptionType optionType = GetEnumFromName
                                   (
                                      node,
                                      XML_TAG_OPTION_TYPE,
                                      SIZEOF(g_optionTypes),
                                      g_optionTypes
                                   );


  finance::ExerciseType exerciseType = GetEnumFromName
                                       (
                                          node,
                                          XML_TAG_OPTION_EXERCISE_TYPE,
                                          SIZEOF(g_exerciseTypes),
                                          g_exerciseTypes
                                       );

  finance::Option *pOption = 
    new Option(dStrike, expiryDate, optionType, exerciseType);

  GetOptionalDerivativeDataFromNode(node, *pOption);

  GetMarketPrice(node, *pOption);

  //get the implied volatility if it exists
  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_OPTION_IMPLIED_VOL);

  if ( pNodeFound != node.end() )
    pOption->SetImpliedVol( GetDoubleFromNode(*pNodeFound) );

  return pOption; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_OPTION_ROOT, Option);

// This funtion is used when reading an EDS termstructure
bool ito33::XML::Restore(const xml::node& node, shared_ptr<Option>& pOption)
{
  if ( strcmp(node.get_name(), XML_TAG_OPTION_ROOT) != 0 )
    return false;

#ifndef NDEBUG
  pOption = shared_ptr<Option>( dynamic_cast<Option*>(ReadOption(&node)) );
#else
  pOption = shared_ptr<Option>( static_cast<Option*>(ReadOption(&node)) );
#endif

  xml::node::const_iterator 
    pNodeFound = node.find(XML_TAG_FINANCE_MARKETPRICE);

  if ( pNodeFound != node.end() )
    pOption->SetMarketPrice( GetDoubleFromNode(*pNodeFound) );

  return true;
}
