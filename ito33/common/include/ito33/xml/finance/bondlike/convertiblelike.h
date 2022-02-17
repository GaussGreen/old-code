/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/bondlike/convertiblelike.h
// Purpose:     Names of elements and attributes used in XML for 
//              a convertible-like bond
// Author:      ITO 33
// Created:     2004/10/13
// RCS-ID:      $Id: convertiblelike.h,v 1.14 2006/04/04 16:29:47 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/bondlike/convertiblelike.h
    @brief Contains the names of the elements used in the XML description of a
           convertible-like bond

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_
#define _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_

#include "ito33/sharedptr.h"

namespace xml
{
  class node;
}

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL ConvertibleLike;
}

namespace XML
{

/**
    Gets optional data of convertilbe Like class from xml node.

    @param node xml node
    @param convertilbeLike object whose optional ConvertibleLike
                           data are required.
  */
void GetOptionalConvertibleLikeDataFromNode
     (const xml::node& node, finance::ConvertibleLike& convertilbeLike);
}

}

/**
   @name Tag name macros
*/
//@{
#define XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_ROOT "exchangeable_feature"

#define XML_TAG_CONVERTIBLELIKE_EXCHANGEABLE_COD  "exchangeable_upon_default"

#define XML_TAG_CONVERTIBLELIKE_NEWSHARE          "new_share"

#define XML_TAG_BONDLIKE_FIXEDFXRATE "fixed_fx_rate"

#define XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_ROOT "fixedquanto_feature"

#define XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_FXVOL "fx_rate_vol"

#define XML_TAG_CONVERTIBLELIKE_FIXEDQUANTO_CORRELATION "correlation"
//@}

#endif // _ITO33_XML_FINANCE_BONDLIKE_CONVERTIBLELIKE_H_
