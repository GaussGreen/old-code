/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/bondliketerms_all.h
// Purpose:     Names of elements and attributes used in XML for 
//              all classes derived from BondLikeTerms
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: bondliketerms_all.h,v 1.28 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/bondliketerms_all.h
    @brief Contains the names of the elements used in the XML description of a
           all possible bond like terms.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_BONDLIKETERMS_ALL_H_
#define _ITO33_XML_FINANCE_BONDLIKE_BONDLIKETERMS_ALL_H_

#include "ito33/sharedptr.h"

namespace xml
{
  class node;
}

//=============================================================================
//=              readers                                                      =
//=============================================================================
namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL BondLikeTerms;
  class ITO33_DLLDECL BondTerms;
}

namespace XML
{

/**
    Restore an BondLikeTerms object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the tag in DOM tree presents the specific terms
    @return the new BondLikeTerms object to be deleted by the caller
 */
shared_ptr<finance::BondLikeTerms>
GetBondLikeTermsFromNode(const xml::node& node);


/**
    Restore an BondTerms object from XML.

    If XML is mal-formed or invalid, an exception is thrown.

    @param node the tag HOLDS bond_terms tag
    @return sharedptr to new BondTerms object 
 */
shared_ptr<finance::BondTerms>
GetBondTermsFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_BONDLIKETERMS_ROOT             "bond_like_terms"

#define XML_TAG_BONDLIKETERMS_CASHDISTRIBUTION "cash_distribution"

//_______Bond Terms____________________________________________________________ 
//
#define XML_TAG_BONDTERMS_ROOT                      "bond_terms"

#define XML_TAG_BONDTERMS_REDEMPTIONPRICE           "redemption_price"

#define XML_TAG_BONDTERMS_OIDYIELD                  "OID_yield"

#define XML_TAG_BONDTERMS_CASHPAYTOZEROACCRETIONRATE "cash_pay_then_OID_yield"

#define XML_TAG_BONDTERMS_COMPOUNDING_FREQUENCY  "compounding_frequency"

#define XML_TAG_BONDTERMS_YIELD_DCC "yield_day_count_convention"

//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_BONDLIKETERMS_ALL_H_
