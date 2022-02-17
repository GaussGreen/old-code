/////////////////////////////////////////////////////////////////////////////
// Name:        xmlreader_fortest.cpp
// RCS-ID:      $Id: xmlreader_fortest.cpp,v 1.3 2004/11/23 15:17:02 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include <fstream>
#include <xmlwrapp/init.h>
#include <xmlwrapp/document.h>
#include <xmlwrapp/tree_parser.h>

#include "ito33/useexception.h"
#include "./xmlreader_fortest.h"

extern const ito33::Error ITO33_BAD_PARAM;

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
// XML::TestingReader implementation
// ============================================================================

// ----------------------------------------------------------------------------
// ctors/dtor
// ----------------------------------------------------------------------------

TestingReader::TestingReader(const char *filename)
      : m_parser(NULL)
{
  InitializeXMLWrapp();

  std::ifstream inFile(filename, std::ios::in);
  if (!inFile)
  {
     throw EXCEPTION_MSG
          ( 
            ITO33_BAD_PARAM,
            TRANS("File does not exists")
          );
  }
  inFile.close();

  m_parser = new xml::tree_parser(filename);

  // TODO: validate?
}

TestingReader::~TestingReader()
{
  delete m_parser;
}
