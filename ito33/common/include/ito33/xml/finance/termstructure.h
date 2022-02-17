/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/finance/termstructure.h
// Purpose:     Names of elements and attributes used in XML for termstructures
// Author:      ZHANG Yunzhi
// Created:     2004-09-24
// RCS-ID:      $Id: termstructure.h,v 1.14 2006/08/19 22:21:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/finance/termstructure.h

    @sa ito33/xml/finance/sessiondata.h
 */

#ifndef _ITO33_XML_FINANCE_TERMSTRUCTURE_H_
#define _ITO33_XML_FINANCE_TERMSTRUCTURE_H_

#include "ito33/sharedptr.h"

/**
   @name Tag name macros
*/
//@{

#define XML_TAG_TERMSTRUCTURE_CDS         "term_structure_cds"
#define XML_TAG_TERMSTRUCTURE_PARBOND     "term_structure_par_bond"
#define XML_TAG_TERMSTRUCTURE_EDS         "term_structure_eds"
#define XML_TAG_TERMSTRUCTURE_OPTION      "term_structure_option"

//@}

namespace xml { class node; }

namespace ito33
{

namespace finance
{
  /**
     @name Forward declaration
   */
//@{
  class TermStructureEnumerator;
  class ITO33_DLLDECL TermStructureCDS;
  class ITO33_DLLDECL TermStructureParBond;
  class ITO33_DLLDECL TermStructureEDS;
  class ITO33_DLLDECL TermStructureOption;
//@}
}

namespace XML
{
/// Gets all possible term structures from XML
void ReadTermStructures
              (
                const xml::node& node,
                finance::TermStructureEnumerator& termstructures
              );

/**
    Restores a CDS curve from XML node 
    if it contains a TermsStrutureCDS node.
 */
bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureCDS>& pCDSCurve);

/**
    Restores a ParBond curve from XML node 
    if it contains a TermsStrutureParBond node.
 */
bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureParBond>& pParBondCurve);

/**
    Restores a EDS curve from XML node 
    if it contains a TermsStrutureEDS node.
 */
bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureEDS>& pEDSCurve);

/**
    Restores a Option curve from XML node 
    if it contains a TermsStrutureOption node.
 */
bool Restore(const xml::node& node, 
             shared_ptr<finance::TermStructureOption>& pOptionCurve);

/// gets a CDS curve from XML, throws exception if it fails.
shared_ptr<finance::TermStructureCDS> GetTermStructureCDSInNode
        (
        const xml::node& node
        );

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_FINANCE_TERMSTRUCTURE_H_
