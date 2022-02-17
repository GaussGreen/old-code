/////////////////////////////////////////////////////////////////////////////
// Name:        hg/xml/spotcomponent.h
// Purpose:     Names of elements and attributes in XML for parametrizations
// Created:     2005/05/20
// RCS-ID:      $Id: parametrization.h,v 1.3 2006/08/19 23:46:52 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file hg/xml/parametrization.h
    @brief Contains the names of the elements used in the XML description of a
           HG parametrization. 
 */

#ifndef _HG_XML_PARAMETRIZATION_H_
#define _HG_XML_PARAMETRIZATION_H_

#include "ito33/sharedptr.h"

#define XML_TAG_HG_CALIBRATION           "calibration"
#define XML_TAG_PARAMETRIZATION_ROOT     "parametrization"
#define XML_TAG_PARAMETRIZATION_INIT     "initial_guess"

namespace xml { class node; }

namespace ito33
{

namespace hg
{

class Parametrization;

namespace XML
{


/**
   Read a Parametrization from XML.

   @param pNode the parametrization node in DOM tree
   @return pointer to new parametrization object
 */
shared_ptr<Parametrization>
ReadParametrizationHRWithTimeComponent(const xml::node *pNode);


} // namespace xml

} // namespace hg

} // namespace ito33

#endif // #ifndef _HG_XML_PARAMETRIZATION_H_
