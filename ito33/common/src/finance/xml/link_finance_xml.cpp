/////////////////////////////////////////////////////////////////////////////
// Name:        link_finance_xml.cpp
// Purpose:     force linking of xml files related to financial classes
// Author:      ZHANG Yunzhi
// Created:     2004-09-02
// RCS-ID:      $Id: link_finance_xml.cpp,v 1.16 2006/07/19 17:38:55 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file must be included in the main project when using XML parsing code
// in a (static) library because otherwise the definition of the static object
// by the macro below would be discarded by the overzealous MSVC linker

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(asianoption_xml);
ITO33_FORCE_LINK_MODULE(option_xml);

ITO33_FORCE_LINK_MODULE(cds_xml);
ITO33_FORCE_LINK_MODULE(referencecds_xml);
ITO33_FORCE_LINK_MODULE(parbond_xml);
ITO33_FORCE_LINK_MODULE(eds_xml);
ITO33_FORCE_LINK_MODULE(logcontract_xml);
ITO33_FORCE_LINK_MODULE(varianceswap_xml);

ITO33_FORCE_LINK_MODULE(bond_xml);
ITO33_FORCE_LINK_MODULE(attachedwarrantcb_xml);
ITO33_FORCE_LINK_MODULE(convertiblebond_xml);
ITO33_FORCE_LINK_MODULE(reset_xml);
ITO33_FORCE_LINK_MODULE(generalizedpepslike_xml);
ITO33_FORCE_LINK_MODULE(pepslike_xml);
ITO33_FORCE_LINK_MODULE(percslike_xml);
ITO33_FORCE_LINK_MODULE(cboption_xml);

ITO33_FORCE_LINK_MODULE(yieldcurve_xml);
