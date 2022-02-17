/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/pepslike.h
// Purpose:     Names of elements and attributes used in XML for 
//              a convertible bond
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: pepslike.h,v 1.5 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/pepslike.h
    @brief Contains the names of the elements used in the XML description of a
           pepslike instrument

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_PEPSLIKE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_PEPSLIKE_H_


/**
   @name Tag name macros
*/
//@{
#define XML_TAG_PEPSLIKE_ROOT                  "peps_like"

#define XML_TAG_PEPSLIKE_MIN_CONVERSION_RATIO "min_conversion_ratio"
                                              
#define XML_TAG_PEPSLIKE_MAX_CONVERSION_RATIO "max_conversion_ratio"

#define XML_TAG_PEPSLIKE_CALL_TYPE                        "call_type"

#define XML_TAG_PEPSLIKE_HAS_OPTIONAL_CONVERSION      "has_optional_conversion"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_PEPSLIKE_H_

