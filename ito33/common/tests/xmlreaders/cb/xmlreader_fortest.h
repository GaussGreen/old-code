/////////////////////////////////////////////////////////////////////////////
// Name:        tests/xmlreader_fortest.h
// Purpose:     reading XML files
// Author:      ZHANG Yunzhi
// Created:     2004-09-03
// RCS-ID:      $Id: xmlreader_fortest.h,v 1.2 2004/10/05 09:13:53 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file tests/xmlreader_fortest.h
    @brief most simplest class for reading XML-files.
 */

#ifndef _ITO33_XML_READER_FOR_TEST_H_
#define _ITO33_XML_READER_FOR_TEST_H_

#include "ito33/sharedptr.h"

namespace xml
{
  class node;
  class tree_parser;
}

namespace ito33
{

namespace XML
{

/**
    This class reads XML document and provides easy access to the data in it.

    Right now XML data can only be loaded from a file but we could also load it
    from a string if needed. Unfortunately, xmlwrapp doesn't take std::istream
    as input...
 */
class TestingReader
{
public:
  /**
    Load XML from the given file.

    If an error occurs while loading date, i.e. file is not found, couldn't be
    read or doesn't contain a valid XML document with the IHG roto tag, an
    exception is thrown.

    @param filename the name of the file to read XML from
   */
  TestingReader(const char *filename); // throws

  /**
    Destructor cleans up resources used for XML parsing.

    Destructor is not virtual, this class is not supposed to be derived from.
   */
  ~TestingReader();

  xml::tree_parser *m_parser;

private:
  NO_COPY_CLASS(TestingReader);

};

}

} // namespace ito33

#endif // _ITO33_XML_READER_FOR_TEST_H_

