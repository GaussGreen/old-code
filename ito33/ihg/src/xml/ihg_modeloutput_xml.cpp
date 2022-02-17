/////////////////////////////////////////////////////////////////////////////
// Name:        xml/ihg_modeloutput_xml.cpp
// Purpose:     Restoring ihg::ModelOutput from XML
// Author:      Vadim Zeitlin
// Created:     2004-05-10
// RCS-ID:      $Id: ihg_modeloutput_xml.cpp,v 1.7 2005/12/30 10:10:09 nabil Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/xml/read.h"

#include "ito33/ihg/modeloutput.h"

#include "ito33/xml/finance/modeloutput.h"
#include "ihg/xml/modeloutput.h"

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
ito33::ihg::XML::Restore(const xml::node& node, ihg::ModelOutput& output)
{
  // read base class fields
  bool rc = ito33::XML::Restore(node, output);

  const xml::node::const_iterator end = node.end();
  xml::node::const_iterator i;
  
  READ_HERE_(DELTA, Delta);
  READ_HERE_(GAMMA, Gamma);
  READ_HERE_(THETA, Theta);
  READ_POSSIBLEVALUE_HERE_(VEGA, Vega);
  READ_POSSIBLEVALUE_HERE_(RHO, Rho);
  READ_POSSIBLEVALUE_HERE_(FUGIT, Fugit);

  return rc;
}

