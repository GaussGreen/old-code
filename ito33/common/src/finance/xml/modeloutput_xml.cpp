/////////////////////////////////////////////////////////////////////////////
// Name:        finance/xml/modeloutput_xml.cpp
// Purpose:     Restoring ModelOutput from XML
// Author:      Vadim Zeitlin
// Created:     2004-05-10
// RCS-ID:      $Id: modeloutput_xml.cpp,v 1.3 2006/02/28 16:13:48 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/xml/read.h"

#include "ito33/finance/modeloutput.h"

#include "ito33/xml/finance/modeloutput.h"

using namespace ito33;

// reads an output or set rc to false if it doesn't exist
#define READ_HERE_(tagname, name) \
{                      \
  i = node.find(XML_TAG_OUTPUT_##tagname); \
  if ( i != end ) \
  { \
  output.Set ##name(ito33::XML::GetDoubleFromNode(*i));\
  }\
  else\
  {\
    rc = false;\
  }\
}

// read an output if it exists
#define READ_POSSIBLEVALUE_HERE_(tagname, name) \
{                      \
  i = node.find(XML_TAG_OUTPUT_##tagname); \
  if ( i != end ) \
  { \
  output.Set ##name(ito33::XML::GetDoubleFromNode(*i));\
  }\
}

bool
ito33::XML::Restore(const xml::node& node, finance::ModelOutput& output)
{
  // read base class fields
  xml::node::const_iterator i = node.find(XML_TAG_OUTPUT_PRICE);
  if ( i == node.end() )
    return false;

  output.SetPrice(XML::GetDoubleFromNode(*i));


  const xml::node::const_iterator end = node.end();
  bool rc = true;
  
  READ_HERE_(DELTA, Delta);
  READ_HERE_(GAMMA, Gamma);
  READ_HERE_(THETA, Theta);
  READ_POSSIBLEVALUE_HERE_(VEGA, Vega);
  READ_POSSIBLEVALUE_HERE_(RHO, Rho);
  READ_POSSIBLEVALUE_HERE_(FUGIT, Fugit);

  return rc;
}
