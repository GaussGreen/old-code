/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/issuer.h
// Purpose:     Names of elements and attributes used in XML for an issuer
// Author:      Nabil
// Created:     2005/03/17
// RCS-ID:      $Id: issuer.h,v 1.6 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004-2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/issuer.h
    @brief Contains the names of the elements used in the XML description of
           the issuer.
 */

#ifndef _ITO33_XML_FINANCE_ISSUER_H_
#define _ITO33_XML_FINANCE_ISSUER_H_

#include "ito33/sharedptr.h"
#include "ito33/xml/finance/common.h"

/**
    @name Names of the elements and values in the XML description of the
          Issuer class.
 */
//@{

#define XML_TAG_ISSUER_ROOT                  "issuer"

#define XML_TAG_ISSUER_FISCALYEARSTART_DATE  "fiscal_year_start_date"

#define XML_TAG_ISSUER_DEFAULTINTENSITY      "default_intensity"

//@}

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  class ITO33_DLLDECL Issuer;
}

namespace XML
{

/**
    Get an Issuer object from XML. 
    
    Exception is thrown if format is wrong.

    @param node the tag in DOM tree
    @return completed Issuer object with information we read from node
 */
shared_ptr<finance::Issuer> GetIssuerFromNode(const xml::node& node);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_ISSUER_H_
