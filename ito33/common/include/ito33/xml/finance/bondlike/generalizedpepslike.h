/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/generalizedpepslike.h
// Purpose:     Names of elements and attributes used in XML for 
//              a generalized peps
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: generalizedpepslike.h,v 1.3 2006/03/24 10:18:27 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/generalizedpepslike.h
    @brief Contains the names of the elements used in the XML description of a
           generalized pepslike instrument
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_


/**
   @name Tag name macros
*/
//@{
#define XML_TAG_GENERALIZED_PEPSLIKE_ROOT               "generalized_peps_like"
#define XML_TAG_GENERALIZED_PEPSLIKE_LOWER_STRIKE                "lower_strike"                                              
#define XML_TAG_GENERALIZED_PEPSLIKE_HIGHER_STRIKE              "higher_strike"
#define XML_TAG_GENERALIZED_PEPSLIKE_DOWNSIDE_CONVERSION_RATIO \
                                                    "downside_conversion_ratio"
#define XML_TAG_GENERALIZED_PEPSLIKE_UPSIDE_BASE_CONVERSION_RATIO \
                                                 "upside_base_conversion_ratio"
#define XML_TAG_GENERALIZED_PEPSLIKE_HAS_OPTIONAL_CONVERSION  \
                                                      "has_optional_conversion"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_GENERALIZED_PEPSLIKE_H_

