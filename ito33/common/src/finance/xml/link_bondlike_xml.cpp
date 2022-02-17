/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/link_bondlike_xml.cpp
// Purpose:     force linking of xml files related to bondlike classes
// Created:     2004/11/22
// RCS-ID:      $Id: link_bondlike_xml.cpp,v 1.5 2005/12/27 15:33:47 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// this file must be included in the main project when using XML parsing code
// in a (static) library because otherwise the definition of the static object
// by the macro below would be discarded by the overzealous MSVC linker

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(cds_xml);
ITO33_FORCE_LINK_MODULE(bond_xml);
ITO33_FORCE_LINK_MODULE(attachedwarrantcb_xml);
ITO33_FORCE_LINK_MODULE(convertiblebond_xml);
ITO33_FORCE_LINK_MODULE(reset_xml);
ITO33_FORCE_LINK_MODULE(pepslike_xml);
ITO33_FORCE_LINK_MODULE(generalizedpepslike_xml);
ITO33_FORCE_LINK_MODULE(percslike_xml);
ITO33_FORCE_LINK_MODULE(yieldcurve_xml);
ITO33_FORCE_LINK_MODULE(cboption_xml);
