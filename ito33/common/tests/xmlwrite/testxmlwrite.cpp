/////////////////////////////////////////////////////////////////////////////
// Name:        test/datadump/main.cpp
// Purpose:     implemtation for XMLWrite test program
// Author:      Vadim Zeitlin
// Created:     22.03.04
// RCS-ID:      $Id: testxmlwrite.cpp,v 1.2 2004/10/05 09:13:53 pedro Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/debug.h"
#include "ito33/string.h"
#include "ito33/cppunit.h"

#include "ito33/xml/write.h"

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/tests/testxmlwrite.h"

using namespace ito33;

// ----------------------------------------------------------------------------
// define a simple "custom" class just to show how it can be dumped too
// ----------------------------------------------------------------------------

struct Person
{
  Person(const char *first, const char *last)
    : m_first(first),
      m_last(last)
  {
  }

  std::string m_first,
              m_last;
};

// defining operator<<() for our class allows not only to print it to cout but
// also to insert it as a text value or as an attribute in the generated XML
std::ostream& operator<<(std::ostream& ostr, const Person& p)
{
  ostr << p.m_last << ", " << p.m_first;

  return ostr;
}

// however if we want to represent it as a complex tag instead of a simple
// text value, we need to specialize this template function instead
namespace ito33
{

namespace XML
{

template <>
Tag MakeNamedTag(const char *name, const Person& p, Tag& parent)
{
  Tag tag(name, parent);
  tag.Element("firstName")(p.m_first);
  tag.Element("lastName")(p.m_last);

  return tag;
}

template <>
Tag MakeTag(const Person& p, Tag& parent)
{
  return MakeNamedTag("person", p, parent);
}

} // namespace XML

} // namespace ito33

// ----------------------------------------------------------------------------
// helper function for checking XML produced by our classes
// ----------------------------------------------------------------------------

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

#define ASSERT_SAME_XML(s1, s2) \
  CPPUNIT_ASSERT( SqueezeXML((s1)) == SqueezeXML((s2)) )

// see comment in XMLWriteTestCase::Empty() to understand why do we need this
// class
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

// ----------------------------------------------------------------------------
// test functions
// ----------------------------------------------------------------------------

void XMLWriteTestCase::Empty()
{
  // this is a bit tricky: we shouldn't check m_oss contents before root is
  // destroyed as the document is only closed in its dtor; so we create a
  // special local object which verifies the condition after destruction of
  // root
  ExpectedXML expected(m_oss, "<?xml version=\"1.0\"?><root/>");

  XML::RootTag root("root", m_oss);
}

void XMLWriteTestCase::Close()
{
  XML::RootTag root("root", m_oss);
  root.Close();

  // here, do check that the document is already ok, before dtor
  ASSERT_SAME_XML( m_oss.str(), "<?xml version=\"1.0\"?><root/>" );
}

void XMLWriteTestCase::Attr()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root foo=\"1\" bar='none'/>"
              );

  XML::RootTag("root", m_oss).Attr("foo", 1).Attr("bar", "none", '\'');
}

void XMLWriteTestCase::AttrQuote()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root nick=\"&quot;foo&quot;\" name='O&apos;Henry'/>"
              );

  XML::RootTag("root", m_oss)
      .Attr("nick", "\"foo\"")
      .Attr("name", "O'Henry", '\'');
}

void XMLWriteTestCase::Value()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "   Hello, XML!\n"
                "   That's all, folks.\n"
                "</root>"
              );

  XML::RootTag("root", m_oss)("Hello, XML!")("That's all, folks.");
}

void XMLWriteTestCase::ValueQuote()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "   &lt;b&gt;Hello&lt;/b&gt; &amp; goodbye, \"XML\"!\n"
                "</root>"
              );

  XML::RootTag("root", m_oss).Value("<b>Hello</b> & goodbye, \"XML\"!");
}

void XMLWriteTestCase::Comment()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "   <!-- generated by XML::Tag -->\n"
                "</root>"
              );

  XML::RootTag("root", m_oss).Comment("generated by XML::Tag");
}

void XMLWriteTestCase::SubElements()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "   <child>\n"
                "       <grandchild>\n"
                "         Foo!\n"
                "       </grandchild>\n"
                "   </child>\n"
                "</root>"
              );

  XML::RootTag("root", m_oss).Element("child").Element("grandchild")("Foo!");
}

void XMLWriteTestCase::SiblingElements()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "   <input>\n"
                "       <value>2</value>\n"
                "       <value>2</value>\n"
                "       <oper>*</oper>\n"
                "   </input>\n"
                "   <output>\n"
                "       <value>5</value>\n"
                "   </output>\n"
                "</root>"
              );

  XML::RootTag root("root", m_oss);
  XML::Tag input("input", root);
  input.Element("value")(2);
  input.Element("value")(2);
  input.Element("oper")('*');
  input.Close();

  root.Element("output").Element("value")(5);
}

void XMLWriteTestCase::RandomExample()
{
  ExpectedXML expected(m_oss,
                "<?xml version=\"1.0\"?>\n"
                "<root format=\"verbose\">\n"
                "    <child>\n"
                "        value\n"
                "    </child>\n"
                "    <child id=\"1\"/>\n"
                "    <another>\n"
                "        <grandchild id=\"2.1\">\n"
                "            first\n"
                "        </grandchild>\n"
                "        <grandchild id=\"2.2\">\n"
                "            second\n"
                "        </grandchild>\n"
                "    </another>\n"
                "    <yet_another>\n"
                "        <!-- another level -->\n"
                "        <grandchild id='3.1'/>\n"
                "    </yet_another>\n"
                "    <author>\n"
                "        Doe, John\n"
                "    </author>\n"
                "    <person>\n"
                "        <firstName>\n"
                "            John\n"
                "        </firstName>\n"
                "        <lastName>\n"
                "            Doe\n"
                "        </lastName>\n"
                "    </person>\n"
                "    <anonymous>\n"
                "        <firstName>\n"
                "            John\n"
                "        </firstName>\n"
                "        <lastName>\n"
                "            Doe\n"
                "        </lastName>\n"
                "    </anonymous>\n"
                "</root>"
              );

  XML::RootTag root("root", m_oss);
  root.Attr("format", "verbose");
  root.Element("child").Value("value");
  root.Element("child").Attr("id", "1");

  XML::Tag child("another", root);
  child.Element("grandchild").Attr("id", "2.1").Value("first");
  child.Element("grandchild").Attr("id", "2.2").Value("second");
  child.Close();

  {
    XML::Tag child2("yet_another", root);
    child2.Comment("another level").Element("grandchild").Attr("id", "3.1", '\'');
  }

  Person person("John", "Doe");
  root.Element("author").Value(person);
  root.Element(person);
  root.Element("anonymous", person);
}
