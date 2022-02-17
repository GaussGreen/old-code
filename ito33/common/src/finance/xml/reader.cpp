/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/xml/reader.cpp
// Purpose:     code for reading ITO33 XML files
// Created:     2005/04/19
// RCS-ID:      $Id: reader.cpp,v 1.5 2006/08/19 22:39:19 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <fstream>

#include "ito33/xml/read.h"
#include "ito33/xml/reader.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/sessiondata.h"
#include "ito33/xml/finance/derivative.h"

#include "ito33/finance/derivative.h"
#include "ito33/finance/derivative_visitor.h"

#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

using namespace ito33;
using namespace ito33::XML;

// ----------------------------------------------------------------------------
// private functions
// ----------------------------------------------------------------------------

// this function is used to initialize xmlwrapp library
static void InitializeXMLWrapp()
{
  static xml::init s_xmlInit;
}

// ============================================================================
// XML::Reader implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctors/dtor
// ----------------------------------------------------------------------------

Reader::Reader(const char *filename) : m_parser(NULL)
{
  InitializeXMLWrapp();

  std::ifstream inFile(filename, std::ios::in);
  
  if (!inFile)
  {
     throw EXCEPTION_MSG
           ( 
             ITO33_BAD_DATA,
             String::Printf(TRANS("File does not exists: %s"), filename)
           );
  }

  inFile.close();

  m_parser = new xml::tree_parser(filename);
}

Reader::Reader(const char *data, size_t len) : m_parser(NULL)
{
  InitializeXMLWrapp();

  m_parser = new xml::tree_parser(data, len);
}

Reader::~Reader()
{
  delete m_parser;
}

// ----------------------------------------------------------------------------
// reading session data information
// ----------------------------------------------------------------------------

shared_ptr<finance::SessionData> Reader::ReadSessionData() const
{
  return GetSessionDataFromNode(
    GetNodeByName(GetMainNode(), XML_TAG_SESSIONDATA));
}

const xml::node& Reader::GetRootNode() const
{
  return m_parser->get_document().get_root_node();
}

xml::node Reader::GetMainNode() const
{
  return GetRootNode();
}

void 
XML::Reader::ReadDerivatives(finance::DerivativeVisitor& visitor) const
{
  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  xml::node nodeDerivs = GetNodeByName(GetMainNode(), XML_TAG_DERIVATIVES);

  const xml::node::const_iterator end = nodeDerivs.end();
  for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
  {
    shared_ptr<finance::Derivative>
      derivative( DerivativeFactory::Create(i->get_name(), &(*i)) );

    if ( derivative )
    {
      derivative->SetSessionData(pSessionData);

      derivative->Visit(visitor);
    }
    //else: ignore unknown derivatives
  }
}

void
XML::Reader::ReadDerivative(shared_ptr<finance::Derivative>& pDerivative) const
{
  shared_ptr<finance::SessionData> pSessionData( ReadSessionData() );

  xml::node nodeDerivs = GetNodeByName(GetMainNode(), XML_TAG_DERIVATIVES);

  const xml::node::const_iterator end = nodeDerivs.end();
  for ( xml::node::const_iterator i = nodeDerivs.begin(); i != end; ++i )
  {
    shared_ptr<finance::Derivative>
      derivative( DerivativeFactory::Create(i->get_name(), &(*i)) );

    if ( derivative )
    {
      derivative->SetSessionData(pSessionData);

      pDerivative = derivative;

      return;
    }
  }
}
