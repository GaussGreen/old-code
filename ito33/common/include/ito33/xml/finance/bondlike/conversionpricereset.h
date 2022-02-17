/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/conversionpricereset.h
// Purpose:     Names of elements and attributes used in XML for 
//              the data at a given conversion price reset date
// Author:      ITO 33
// Created:     2004/10/13
// RCS-ID:      $Id: conversionpricereset.h,v 1.6 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/conversionpricereset.h
    @brief Contains the names of the elements used in the XML description of 
           the conversion price reset data at a given date

    @sa ito33/xml/finance/bondlike/conversionpricereset.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVPRICERESET_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVPRICERESET_H_

#include "ito33/sharedptr.h"

//=============================================================================
//=              readers                                                      =
//=============================================================================

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ConversionPriceReset;
}

namespace XML
{
  /// Restore conversion price reset from XML node
  bool Restore
        (
          const xml::node& node,
          shared_ptr<finance::ConversionPriceReset>& pConversionPriceReset
        );
}

}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_CONVPRICERESET_ROOT                "conversion_price_reset"

#define XML_TAG_CONVPRICERESET_CAPRATE             "cap_rate"

#define XML_TAG_CONVPRICERESET_FLOORRATE           "floor_rate"

#define XML_TAG_CONVPRICERESET_MULTIPLIER          "multiplier"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVPRICERESET_H_

