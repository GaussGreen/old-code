/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/varianceswap.h
// Purpose:     Names of elements and attributes in XML for variance swaps
// Created:     2006/02/21
// RCS-ID:      $Id: varianceswap.h,v 1.13 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/varianceswap.h
    @brief Contains the names of the elements used in the XML description of
           variance swaps.
 */

#ifndef _ITO33_XML_FINANCE_VARIANCESWAP_H_
#define _ITO33_XML_FINANCE_VARIANCESWAP_H_

#include "ito33/sharedptr.h" 

//=============================================================================
//=  @name            tag name macros                                              =
//=============================================================================
//@{

#define XML_TAG_VARIANCESWAP_ROOT               "variance_swap"
#define XML_TAG_GAMMAVARIANCESWAP_ROOT          "gamma_variance_swap"
#define XML_TAG_OPTIONVARIANCESWAP_ROOT         "option_variance_swap"
#define XML_TAG_CONDITIONALVARIANCESWAP_ROOT    "conditional_variance_swap"
#define XML_TAG_VARIANCESWAP_VOLATILITYSTRIKE   "volatility_strike"
#define XML_TAG_VARIANCESWAP_STARTOFPERIOD      "start_of_sampling_period"
#define XML_TAG_VARIANCESWAP_NBSAMPLINGRETURNS  "number_of_sampling_returns"
#define XML_TAG_VARIANCESWAP_CURRENT_VOLATILITY "current_volatility"
#define XML_TAG_VARIANCESWAP_NB_SAMPLES_USED    "number_of_samples_used"
#define XML_TAG_VARIANCESWAP_CURRENT_COUNT      "current_conditional_count"
#define XML_TAG_VARIANCESWAP_START_SHARE_PRICE  "start_share_price"
#define XML_TAG_VARIANCESWAP_CAP_MULTIPLIER     "cap_multiplier"
#define XML_TAG_VARIANCESWAP_UP_CORRIDOR        "up_corridor_barrier"
#define XML_TAG_VARIANCESWAP_DOWN_CORRIDOR      "down_corridor_barrier"
#define XML_TAG_VARIANCESWAP_ANNUALRETURNFREQUENCY "annual_return_frequency"

#define XML_TAG_VARIANCESWAPTION_ROOT           "variance_swaption"

#define XML_TAG_VARIANCESWAPTION_VSTERMS        "variance_swap_terms"

 //@}

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL VarianceSwapTerms;
}

namespace XML
{

/// Read variance swap terms from a node
shared_ptr<finance::VarianceSwapTerms>
GetVarianceSwapTermsFromNode(const xml::node& node);

}

} // namespace ito33

#endif // _ITO33_XML_FINANCE_VARIANCESWAP_H_
