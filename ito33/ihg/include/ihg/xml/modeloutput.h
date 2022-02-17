/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/xml/modeloutput.h
// Purpose:     Restoring ihg::ModelOutput from XML
// Author:      Vadim Zeitlin
// Created:     2004-05-10
// RCS-ID:      $Id: modeloutput.h,v 1.7 2006/01/05 14:47:08 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/xml/modeloutput.h
    @brief Support for restoring ihg::ModelOutput from XML.
 */

#ifndef _ITO33_IHG_XML_MODELOUTPUT_H_
#define _ITO33_IHG_XML_MODELOUTPUT_H_

#define XML_TAG_OUTPUT_VEGA             "vega"
#define XML_TAG_OUTPUT_RHO              "rho"
#define XML_TAG_OUTPUT_UNDERLYING_RHO   "underlying_rho"
#define XML_TAG_OUTPUT_FUGIT            "fugit"


namespace xml { class node; }

namespace ito33
{

namespace ihg
{

class ModelOutput;

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

    @todo should use visitor as ModelOutput may have sub-classes
 */
bool Restore(const xml::node& node, ihg::ModelOutput& output);

} // namespace XML

} // namespace ihg

} // namespace ito33

#endif // _ITO33_IHG_XML_MODELOUTPUT_H_
