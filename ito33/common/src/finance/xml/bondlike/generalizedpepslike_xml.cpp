/////////////////////////////////////////////////////////////////////////////
// Name:        generalizedpepslike_xml.cpp
// Purpose:     Restore Generalized PEPS Like object from XML document
// Author:      ZHANG Yunzhi
// Created:     2005-03-29
// RCS-ID:      $Id: generalizedpepslike_xml.cpp,v 1.6 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/link.h"

#include "ito33/finance/bondlike/generalizedpepslike.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/finance/bondlike/pepsaveragingperiod.h"

#include "ito33/xml/read.h"
#include "ito33/xml/finance/derivative.h"
#include "ito33/xml/finance/bondlike/generalizedpepslike.h"
#include "ito33/xml/finance/bondlike/generalizedpepslikecall.h"
#include "ito33/xml/finance/bondlike/bondliketerms_all.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"
#include "ito33/xml/finance/bondlike/callschedule.h"
#include "ito33/xml/finance/bondlike/pepsaveragingperiod.h"

using namespace ito33;
using namespace ito33::XML;
using namespace ito33::finance;

ITO33_FORCE_LINK_THIS_MODULE(generalizedpepslike_xml);

/**
    Restore an Generalized PEPS-like object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the session tag in DOM tree
    @return the new bond object to be deleted by the caller
 */
static
finance::Derivative *ReadGeneralizedPEPSLike(const xml::node *pNode)
{
  shared_ptr<BondLikeTerms> pBondLikeTerms = GetBondLikeTermsFromNode(*pNode);

  double dLowerStrike = GetDoubleFromName
            (*pNode, XML_TAG_GENERALIZED_PEPSLIKE_LOWER_STRIKE);
  double dHigherStrike = GetDoubleFromName
            (*pNode, XML_TAG_GENERALIZED_PEPSLIKE_HIGHER_STRIKE);
  double dDownside = GetDoubleFromName
            (*pNode, XML_TAG_GENERALIZED_PEPSLIKE_DOWNSIDE_CONVERSION_RATIO);
  double dUpside = GetDoubleFromName
           (*pNode, XML_TAG_GENERALIZED_PEPSLIKE_UPSIDE_BASE_CONVERSION_RATIO);

  GeneralizedPEPSLike
    *pPEPSLike = new GeneralizedPEPSLike
                        (
                          pBondLikeTerms,
                          dDownside,
                          dLowerStrike,
                          dUpside,
                          dHigherStrike
                        );

  GetOptionalConvertibleLikeDataFromNode(*pNode, *pPEPSLike);

  // xml::node::const_iterator *pNodeFound;

  // we can forget call_type,
  // pNodeFound = pNode->find(XML_TAG_PEPSLIKE_CALL_TYPE);

  // optional call provision  
  shared_ptr<CallSchedule> pCall;
  if(Restore(*pNode, pCall))
  {
    pPEPSLike->SetCallFixedCash(pCall);
  }
  else
  {
    shared_ptr<GeneralizedPEPSLikeCall> pGeneralizedPEPSLikeCall;

    if(Restore(*pNode, pGeneralizedPEPSLikeCall))
      pPEPSLike->SetGeneralizedPEPSLikeCall(pGeneralizedPEPSLikeCall);
  }

  // conversion provision
  xml::node::const_iterator 
    pNodeFound = pNode->find
                 (XML_TAG_GENERALIZED_PEPSLIKE_HAS_OPTIONAL_CONVERSION);

  if ( pNodeFound != pNode->end() && GetBoolFromNode(*pNodeFound) )
    pPEPSLike->EnableOptionalConversion(); 


  //Look for the averaging part
  xml::node::const_iterator pNodeAveraging;
  pNodeAveraging = pNode->find(XML_TAG_FINANCE_BONDLIKE_PEPS_AVERAGING_PERIOD);

  if (  pNodeAveraging != pNode->end() )
  {
    Date avgStartDate = GetDateFromName(*pNodeAveraging,
      XML_TAG_FINANCE_BONDLIKE_PEPS_START_DATE);
    
    Date avgEndDate = GetDateFromName(*pNodeAveraging,
      XML_TAG_FINANCE_BONDLIKE_PEPS_END_DATE);

    size_t nNbSampling = GetLongFromName(*pNodeAveraging,
      XML_TAG_FINANCE_BONDLIKE_PEPS_NB_SAMPLING_AVERAGES);

    long nIsStockAveraging = GetLongFromName
                             ( *pNodeAveraging,
                               XML_TAG_FINANCE_BONDLIKE_PEPS_STOCK_AVERAGING );

    shared_ptr<PEPSAveragingPeriod> pAvg;

    if ( nIsStockAveraging )
      pAvg = PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nNbSampling);
    else
      pAvg = PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, avgEndDate, nNbSampling);
   

    //look for the current average
    xml::node::const_iterator pNodeCurrentAverage;
    pNodeCurrentAverage = 
      pNodeAveraging->find(XML_TAG_FINANCE_BONDLIKE_PEPS_CURRENT_AVERAGE);

    if ( pNodeCurrentAverage != pNodeAveraging->end() )
    {
      double dCurrentAverage = GetDoubleFromNode(*pNodeCurrentAverage);

      size_t nNbSamplesUsed  = GetLongFromName(*pNodeAveraging,
        XML_TAG_FINANCE_BONDLIKE_PEPS_NB_SAMPLES_USED);
      
      if ( nIsStockAveraging )
        pAvg->SetCurrentStockAverage(dCurrentAverage, nNbSamplesUsed);
      else
        pAvg->SetCurrentConversionRatioAverage(dCurrentAverage, nNbSamplesUsed);

    }

    pPEPSLike->SetAveragingPeriod(pAvg);
  }
 
  return pPEPSLike;
}


ITO33_DEFINE_DERIVATIVE_READER(XML_TAG_GENERALIZED_PEPSLIKE_ROOT, GeneralizedPEPSLike);
