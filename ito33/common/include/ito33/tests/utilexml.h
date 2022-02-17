/////////////////////////////////////////////////////////////////////////////
// Name:        common/include/ito33/tests/utilexml.h
// Purpose:     functions that permits the comparison of two streams
//              this is very handy for testing the dump function 
// Author:      Yann d'Halluin
// Created:     1/10/2004
// RCS-ID:      $Id: utilexml.h,v 1.2 2004/10/05 09:13:39 pedro Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"


#include "ito33/xml/write.h"

// this function simply removes all white space between tags
static std::string SqueezeXML(const std::string& src)
{
  std::string dst;

  const size_t len = src.length();
  dst.reserve(len);

  bool atStartOfLine = true;
  for ( size_t n = 0; n < len; ++n )
  {
    switch ( const char ch = src[n] )
    {
      case '\r':
      case '\n':
        atStartOfLine = true;
        break;

      case ' ':
        if ( atStartOfLine )
          break;
        //else: fall through

      default:
        dst += ch;
    }
  }

  return dst;
}

#define ASSERT_SAME_XML(s1, s2) CPPUNIT_ASSERT( SqueezeXML((s1)) == SqueezeXML((s2)) )

class ExpectedXML
{
public:
  ExpectedXML(const std::ostringstream& oss, const char *expected)
    : m_oss(oss),
      m_expected(expected)
  {
  }

  ~ExpectedXML()
  {
    ASSERT_SAME_XML( m_oss.str(), m_expected );
  }

private:
  const std::ostringstream& m_oss;
  const std::string m_expected;

  NO_COPY_CLASS(ExpectedXML);
};
