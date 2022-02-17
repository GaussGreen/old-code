/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/modeloutput.h
// Purpose:     Names of elements and attributes used in XML for ModelOutput
// Author:      Vadim Zeitlin
// Created:     2004-05-10
// RCS-ID:      $Id: modeloutput.h,v 1.11 2006/07/11 14:28:11 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/modeloutput.h
    @brief Contains the names of the elements used in the XML description of
           base model output.

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_MODELOUTPUT_H_
#define _ITO33_XML_FINANCE_MODELOUTPUT_H_

#include "ito33/finance/modeloutput.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_OUTPUT_PRICE            "price"
#define XML_TAG_OUTPUT_DELTA            "delta"
#define XML_TAG_OUTPUT_FXDELTA          "fx_delta"
#define XML_TAG_OUTPUT_GAMMA            "gamma"
#define XML_TAG_OUTPUT_THETA            "theta"
#define XML_TAG_OUTPUT_VEGA             "vega"
#define XML_TAG_OUTPUT_RHO              "rho"
#define XML_TAG_OUTPUT_UNDERLYING_RHO   "underlying_rho"
#define XML_TAG_OUTPUT_FUGIT            "fugit"

//@}

namespace xml { class node; }

namespace ito33
{

namespace XML
{

/**
    Restore a ModelOutput object from XML.

    If any unrecognized data is found, the function skips it. If some expected
    data is not found, the function still reads the data in the node, but
    returns false and not true in the end. This behaviour should ensure that we
    can read newer versions of our XML and allows the caller a choice as to
    what to do if we have an older version.

    @param node the session tag in DOM tree
    @param output to fill in with information we read from node
    @return true if we found all the data we expected, false if some expected
            tags or values were not found

    @todo This function returns boolean as we don't care whether the xml file
    has full output data or not. However, we don't have solution yet to handle
    the case where Restore returns false.
 */
bool Restore(const xml::node& node, finance::ModelOutput& output);

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_MODELOUTPUT_H_
