/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/xml/bondlike/bondlikeoutput_xml.cpp
// Purpose:     Restoring ihg::BondLikeOutput from XML
// Author:      ITO33 Canada 
// Created:     September 27, 2005 
// RCS-ID:      $Id: bondlikeoutput_xml.cpp,v 1.2 2006/02/28 16:13:48 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/xml/read.h"

#include "ito33/finance/bondlike/bondlikeoutput.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/xml/finance/bondlike/bondlikeoutput.h"
#include "ito33/xml/finance/modeloutput.h"

using namespace ito33;

bool
ito33::XML::Restore(const xml::node& node, finance::BondLikeOutput& output)
{
  // read base class fields
  bool rc = XML::Restore(node, static_cast<finance::ModelOutput&>(output));

  const xml::node::const_iterator end = node.end();
  xml::node nodeBond = ito33::XML::GetNodeByName(node, XML_TAG_OUTPUT_BONDLIKE);

  xml::node::const_iterator i = nodeBond.find(XML_TAG_OUTPUT_BONDLIKE_FLOOR);
  if ( i == nodeBond.end() )
  {
    rc = false;
  }
  else // found it
  {
    output.SetBondFloor(ito33::XML::GetDoubleFromNode(*i));
  }

  return rc;
}

