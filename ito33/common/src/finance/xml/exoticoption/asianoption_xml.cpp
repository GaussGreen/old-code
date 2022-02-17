/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/exoticoption/asianoption_xml.cpp
// Purpose:     Restore an Asian Option object from XML document
// Author:      ITO 33 Canada
// Created:     April 13, 2005
// RCS-ID:      $Id: asianoption_xml.cpp,v 1.6 2006/07/28 21:01:10 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/exoticoption/asianoption.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/option.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/exoticoption/asianoption.h"
#include "ito33/xml/finance/optiontype.h"
#include "ito33/xml/finance/exercisetype.h"
#include "ito33/xml/finance/frequency.h"

using namespace ito33;
using namespace ito33::XML;

ITO33_FORCE_LINK_THIS_MODULE(asianoption_xml);

/**
    Restore an Asian Option object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the parent tag in DOM tree
    @return the new option object to be deleted by the caller
 */
static
finance::Derivative *ReadAsianOption(const xml::node *pNode)
{
  const xml::node& node = *pNode;
  
  //read maturity
  Date maturityDate = GetDateFromName(node, XML_TAG_FINANCE_MATURITY);

  finance::ExerciseType exerciseType = GetEnumFromName
                                       (
                                          node,
                                          XML_TAG_OPTION_EXERCISE_TYPE,
                                          SIZEOF(g_exerciseTypes),
                                          g_exerciseTypes
                                       );

  finance::OptionType optionType = GetEnumFromName
                                      (
                                         node,
                                         XML_TAG_OPTION_TYPE,
                                         SIZEOF(g_optionTypes),
                                         g_optionTypes
                                       );

  size_t nNumberOfSamplingAverages 
    = GetLongFromName(node, XML_TAG_ASIAN_OPTION_NB_SAMPLING_AVERAGES);
  
  Date averageStartDate = 
    GetDateFromName(node, XML_TAG_ASIAN_OPTION_AVG_START_DATE);

  //Read the current average optional
  double dCurrentAverage  = -1.;
  size_t nNbSamplesUsed   = 0 ;
  xml::node::const_iterator pNodeCurrentAverage;
  if(  (pNodeCurrentAverage = node.find(XML_TAG_ASIAN_OPTION_CURRENT_AVERAGE)) 
     != node.end() )
  {
    dCurrentAverage = GetDoubleFromNode(*pNodeCurrentAverage);

    nNbSamplesUsed = GetLongFromName(node, XML_TAG_ASIAN_OPTION_NB_SAMPLES_USED);
  }

  //Strike is optional
  double dStrike = -1;
  xml::node::const_iterator pNodeFound;
  
  pNodeFound = node.find(XML_TAG_FINANCE_STRIKE);
  if ( pNodeFound != node.end() )
    dStrike = GetDoubleFromNode(*pNodeFound);
  
  //end average date is optional 
  bool bHasAverageEndDate = false;
  Date averageEndDate;
  xml::node::const_iterator pNodeAverageEndDate;
 
  pNodeAverageEndDate = node.find( XML_TAG_ASIAN_OPTION_AVG_END_DATE);
  if ( pNodeAverageEndDate != node.end() )
  {
    bHasAverageEndDate = true;
    averageEndDate = GetDateFromNode(*pNodeAverageEndDate);
  }

  finance::AsianOption *deriv;
  if ( dStrike > 0 )
    deriv = new finance::AsianOption(dStrike, 
                                     maturityDate, 
                                     optionType, 
                                     exerciseType, 
                                     averageStartDate,
                                     nNumberOfSamplingAverages);
  else
    deriv = new finance::AsianOption(maturityDate, 
                                     optionType, 
                                     exerciseType, 
                                     averageStartDate,
                                     nNumberOfSamplingAverages);

  if ( dCurrentAverage > 0 )
    deriv->SetCurrentAverage(dCurrentAverage, nNbSamplesUsed);

  if ( bHasAverageEndDate )
    deriv->SetAverageEndDate(averageEndDate);

  GetOptionalDerivativeDataFromNode(node, *deriv);

  GetMarketPrice(node, *deriv);

  return deriv; 
}

ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_ASIAN_OPTION_ROOT, AsianOption);
