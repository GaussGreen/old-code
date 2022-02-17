/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/write.h
// Purpose:     Classes for generating XML output
// Author:      Vadim Zeitlin
// Created:     23.03.04
// RCS-ID:      $Id: write.h,v 1.17 2006/08/02 13:29:02 cosmin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/write.h
    @brief Allows to easily output any data in XML.

    Writing XML is simple, at least compared to parsing it, but still, doing
    it manually is not very plesant, if only because you have to carefully
    close all opened tags, quote all attribute values and so on. The classes
    here free you of such mundane choirs.
 */

#ifndef _ITO33_XML_WRITE_H_
#define _ITO33_XML_WRITE_H_

#include "ito33/common.h"
#include "ito33/string.h"
#include "ito33/debug.h"

#include "ito33/beforestd.h"
#include <iostream>
#include <sstream>
#include "ito33/afterstd.h"

namespace ito33
{

/**
    This namespace contains everything related to generating and parsing XML.
 */
namespace XML
{

class Tag;

/**
    This function is used to output an object of any type as a complex XML tag.

    The XML writing classes are extensible in that they allow to output not
    just the values of the fundamental types (int, char *, ...) but also of
    arbitrary user-defined classes. This can be done in two ways:
      - simply defining operator<<(const T&, ostream&) for the type T allows
        to define the text representation of a value of type T in XML and is
        useful for values represented using simple tags, such as those
        defining some primitive values
      - specializing this function for the type T allows to output the value
        as a tag of arbitrary complexity in the XML and so is usually done for
        structures containing several orthogonal elements

    To do the latter you simply have to define template specialization of this
    function for your type T, create a new Tag object there and add as many
    attributes, elements and text to it as needed. Here is an example:
      @code
        struct PersonName
        {
          std::string first,
                      last;
        };

        namespace XML
        {
          template <>
          Tag MakeTag(const PersonName& pn, Tag& parent)
          {
            Tag tag("person", parent);
            tag.Element("firstName")(pn.first);
            tag.Element("lastName")(pn.last);

            return tag;
          }

        } // namespace XML
      @endcode

    Note that there is no default implementation for this template function,
    so if you use Tag::Element(const T&), you must define it.

    Finally, you should never have to call this function yourself,
    Tag::Element() does it when needed.

    @sa Tag::Element()

    @param value the value to output
    @param parent the parent tag
    @return the new tag containing the representation of value in XML
 */
template <typename T>
// intel compiler doesn't need and like this work around specific for msvc
#ifndef __INTEL_COMPILER  
extern inline
#endif

Tag MakeTag(const T& value, Tag& parent);

/**
    Similar to MakeTag() but takes the tags name as extra parameter.

    MakeTag() is used for the objects which only occur once in the output and
    in this case it is convenient that it chooses the tag name itself, but if
    several different objects of the same class may appear at the same level,
    we need to give different names to different objects and this is where
    MakeNamedTag() comes into play.

    If you define MakeNamedTag() for any type, it is trivial to define
    MakeTag() for it as well by calling MakeNamedTag() with some default name.

    @param name the name of the tag to create
    @param value the value to output
    @param parent the parent tag
    @return the new tag containing the representation of value in XML
 */
template <typename T>
// intel compiler doesn't need and like this work around specific for msvc
#ifndef __INTEL_COMPILER  
extern inline
#endif

Tag MakeNamedTag(const char *name, const T& value, Tag& parent);

/**
    This class represents any XML element.

    An XML element, or tag, may have any number of attributes and embedded
    tags (which may contain other tags themselves and so on). This clearly
    calls for a hierarchical structure and we define it using Tag for the
    nodes of this tree. The root of the tree is represented by RootTag which
    is also a Tag, albeit a special one.

    Tag provides methods to add attributes and elements to it, but in the
    simplest case you may not call them at all -- this will create an empty
    element. In general, all Tag methods return the reference to the object
    itself which allows to chain calls to them. For not too complex tags this
    means that you don't even have to name the Tag object but can simply use a
    temporary.

    Here are some examples of using this class. First, here is an empty tag:
      @code
        void AddEmptyTag(const RootTag& root)
        {
          root.Element("empty_tag");
        }
      @endcode
    Although the name of the class doesn't appear at all in the example, it is
    still used because Tag::Element() returns a temporary Tag object.

    Adding attributes is easy:
      @code
        void AddTagWithAttr(const RootTag& root)
        {
          root.Element("attr_tag").Attr("attr_name", "attr_value");
        }
      @endcode

    ... as is adding sub-elements:
      @code
        void AddTagWithSubTags(const RootTag& root)
        {
          root.Element("withsub_tag").Element("sub_tag").Element("sub_sub_tag");
        }
      @endcode

    Of course, both can be combined:
      @code
        void AddComplexTag(const RootTag& root)
        {
          root.Element("complex_tag").Attr("foo", "1.0").Element("bar");
        }
      @endcode
    however you @b must call Attr() before calling Element() in such case.

    As explained in MakeTag() description, it is also possible to use any
    other type and not just strings with both Attr() and Element(). For the
    former, the type object is converted to string using operator<<(ostream&) 
    which should be defined:
      @code
        struct Foo { ... };
        std::ostream& operator<<(std::ostream& o, const Foo& foo) { ... }

        void AddUserAttr(const RootTag& root)
        {
          root.Element("userattr_tag").Attr("foo", Foo(...));
        }
      @endcode

    For the latter, MakeTag() should be defined. Reusing the class from the
    MakeTag example we may write
      @code
        void AddPersonTag(const RootTag& root)
        {
          root.Element(Person("John", "Doe"));
        }
      @endcode

        
    Comment about the generated XML: we make no attempt to ensure that it is
    valid or even well-formed, it is the callers responsability to check for
    potential violations.
 */
class Tag
{
private:
  // common part of all ctors: sets all member variables to initial values
  void Init(unsigned indent)
  {
    m_isEmpty = true;
    m_isClosed = false;
    m_isTerse = false;
    m_indent = indent;
  }

protected:
  /// return the string representing the indent of this tag plus extra indent
  std::string Indent(unsigned extra = 0)
  {
    return std::string((m_indent + extra)*4, ' ');
  }

  /// opens this tag
  void BeginStartTag()
  {
    m_ostr << Indent() << "<" << m_name;
  }

  /// closes this tag if this hadn't been done yet, it is safe to call this
  /// function many times
  void EndStartTag()
  {
    if ( m_isEmpty )
    {
      m_ostr << '>';
      if ( !m_isTerse )
        m_ostr << '\n';
      m_isEmpty = false;
    }
  }

  /// common part of QuoteAttr() and QuoteValue()
  static std::string DoQuote(const std::string& src,
                             const char *orig,
                             const char **replacements)
  {
    // we could do something clever with stringstream and transform() but it
    // would be more complicated than this simple loop...
    std::string dst;
    const size_t len = src.length();
    dst.reserve(len);

    for ( size_t n = 0; n < len; ++n )
    {
      // check the current character
      const char ch = src[n];
      const char *p = strchr(orig, ch);
      if ( p )
      {
        // replace with the corresponding entity
        dst += '&';
        dst += replacements[p - orig];
        dst += ';';
      }
      else // pass through as is
      {
        dst += ch;
      }
    }

    return dst;
  }

  /// returns an attribute value with all quote characters replaced with the
  /// corresponding XML entities
  static std::string QuoteAttr(const std::string& s, char quote)
  {
    static const char *replacements[] = { "quot", "apos" };
    char orig[2];
    orig[0] = quote;
    orig[1] = '\0';

    return DoQuote(s, orig, quote == '"' ? replacements : replacements + 1);
  }

  /// returns an attribute value with '<', '>; and '&' characters replaced with
  /// the corresponding XML entities
  static std::string QuoteValue(const std::string& s)
  {
    static const char *replacements[] = { "lt", "gt", "amp" };

    return DoQuote(s, "<>&", replacements);
  }

public:
  /**
      Create a tag with the given name under parent tag.

      Immediately opens the tag, it will be closed when the object is
      destroyed or when Close() is called. Between them Attr(), Element() and
      Value() may be called to add contents to the tag.

      @param name of the tag, should be a valid XML identifier
      @param parent the parent tag, either another Tag or RootTag
   */
  Tag(const std::string& name, Tag& parent)
    : m_ostr(parent.m_ostr),
      m_nPrecision(parent.m_nPrecision),
      m_name(name)
  {
    Init(parent.m_indent + 1);

    parent.EndStartTag();

    BeginStartTag();
  }

  /**
      Dangerous copy constructor.

      Normally it doesn't make sense to copy the Tag objects because each of
      them correponds to an opening and closing tag in the XML and we don't
      want to close the tag twice if a Tag object is copied to another before
      both the original and new copy are destroyed.

      Unfortunately, due to C++ rules, an object without copy ctor can't be
      returned from a function and we really want to allow this so that
      MakeTag() could work. Thus we're forced to define a dangerous copy
      constructor which @b modifies its "const" argument. To be precise, it
      prevents the source tag from being closed after it had been copied. Of
      course, the source tag must not be used any more after this neither,
      once again this copy constructor only exists to allow returning Tags
      from the functions.

      @param tag the source object which must not be used after this call
   */
  Tag(const Tag& tag)
    : m_ostr(tag.m_ostr),
      m_nPrecision(tag.m_nPrecision),
      m_name(tag.m_name)
  {
    m_isEmpty = tag.m_isEmpty;
    m_isTerse = tag.m_isTerse;
    m_indent = tag.m_indent;

    // we should close the tag exactly once!
    if ( (m_isClosed = tag.m_isClosed) == false )
    {
      const_cast<Tag&>(tag).m_isClosed = true;
    }
  }

  /**
    Set the number of digits to display in a floating-point number

    @param nPrecision number of digits
    */
  void precision(std::streamsize nPrecision)
  {
    m_nPrecision = nPrecision;
    m_ostr.precision(m_nPrecision);
  }

  /**
      Add an attribute to this tag.

      Adds an attribute with the given name and value to the tag. The value
      may be of any type which can be streamed in an std::ostream, i.e. for
      which the <tt>operator<<(std::ostream&, const T&)</tt> is defined.

      This method cannot be called any longer if either Element(), Value() or
      Comment() had been called. If this happens, an assert failure is
      generated in debug build and badly formed XML is generated silently in
      release.

      @param name of the attribute, should be a valid XML identifier
      @param value value of the attribute which is written using operator<<()
      @param quote the character to use for quoting the attribute value,
                   double quote by default but may also be a single quote; no
                   other values are valid
      @return this tag so that calls to Attr() may be chained
   */
  template <typename T>
  Tag& Attr(const char *name, const T& value, char quote = '"')
  {
    ASSERT_MSG( m_isEmpty,
                String::Printf("XML::Tag %s: attr %s added too late",
                               m_name.c_str(), name) );

    std::ostringstream ostrstr;
    if(m_nPrecision > 0)
      ostrstr.precision(m_nPrecision);
    ostrstr << value;
    m_ostr << " " << name << "=" << quote
           << QuoteAttr(ostrstr.str(), quote)
           << quote;

    return *this;
  }

  /**
      Add an empty new subtag to this tag.

      The new tag is created with the given name and this tag as parent and
      returned by this method so that Attr(), Element() or Value() could be
      called on it in turn.

      @param name of the new subtag, should be a valid XML identifier
      @return the new subtag
   */
  Tag Element(const char *name)
  {
    ASSERT_MSG( !m_isClosed,
                String::Printf("XML::Tag %s: element %s added after Close()",
                               m_name.c_str(), name) );

    return Tag(name, *this);
  }

  /**
      Add a new complex subtag to this tag.

      A new tag is created using MakeTag() and added to this tag as subtag.
      The contents of the new tag is defined by the specialization of MakeTag
      for the type T.

      @param value the value passed to MakeTag() to create the subtag
      @return the new subtag so that more elements or text could be added to it,
              however this is usually not done because MakeTag() normally
              fills the tag entirely.
   */
  template <typename T>
  Tag Element(const T& value)
  {
    ASSERT_MSG( !m_isClosed,
                String::Printf("XML::Tag %s: complex element added after Close()",
                               m_name.c_str()) );

    EndStartTag();

    return MakeTag(value, *this);
  }

  /**
      Add a new complex subtag with the given name to this tag.

      A new tag is created using MakeNamedTag() and added to this tag as subtag.
      The contents of the new tag is defined by the specialization of
      MakeNamedTag for the type T.

      @param name of the new tag
      @param value the value passed to MakeNamedTag() to create the subtag
      @return the new subtag
   */
  template <typename T>
  Tag Element(const char *name, const T& value)
  {
    ASSERT_MSG( !m_isClosed,
                String::Printf("XML::Tag %s: complex element added after Close()",
                               m_name.c_str()) );

    EndStartTag();

    return MakeNamedTag(name, value, *this);
  }

  /**
      Adds some raw text to the tag.

      This method may be called multiple times to fill in the tag contents.
      The contents of value is streamed in the generated XML using operator<<(),
      just as for Attr().

      Any special XML characters (e.g. '<') are quoted automatically, there is
      no need to worry about them.

      @param value the value to be written to XML using operator<<() for it
      @return this tag itself so that other calls to Element() or Value()
              could be chained to this one
   */
  template <typename T>
  Tag& Value(const T& value)
  {
    ASSERT_MSG( !m_isClosed,
                String::Printf("XML::Tag %s: value added after Close()",
                                m_name.c_str()) );

    // if there hadn't been anything before the value yet, assume that it is
    // "terse", i.e. should be output on the same line as the element
    //
    // this is certainly a hack but one which works in 90% of cases and makes
    // the generated XML slightly easier on the eyes, so why not
    if ( m_isEmpty )
      m_isTerse = true;
    //else: can't be a terse element if it already has anything

    EndStartTag();

    if ( !m_isTerse )
      m_ostr << Indent(+1);

    m_ostr << QuoteValue(ToXMLValue(value));

    if ( !m_isTerse )
      m_ostr << '\n';

    return *this;
  }

  /**
      Shorter form of Value().

      This is exactly the same as Value() but requires less typing.
   */
  template <typename T>
  Tag& operator()(const T& value)
  {
    return Value(value);
  }

  /**
      Insert a comment inside this tag.

      Comment must @b not contain two consecutive dashes ("--"), this would
      make the generated XML badly formed.

      @param comment the text of the comment
      @return this tag itself
   */
  Tag& Comment(const std::string& comment)
  {
    ASSERT_MSG( !m_isClosed,
                String::Printf("XML::Tag %s: comment added after Close()",
                               m_name.c_str()) );

    EndStartTag();

    m_ostr << Indent(+1) << "<!-- " << comment << " -->\n";

    return *this;
  }

  /**
      Close the current tag explicitly.
   */
  void Close()
  {
    ASSERT_MSG( !m_isClosed, 
                String::Printf("XML::Tag %s closed twice.", m_name.c_str()) );

    m_isClosed = true;

    if ( m_isEmpty )
    {
      // special case of an empty tag
      m_ostr << "/>\n";
    }
    else
    {
      // normal tag already containing something
      if ( !m_isTerse )
        m_ostr << Indent();
      m_ostr << "</" << m_name << ">\n";
    }
  }

  /**
      Destructor closes the tag if Close() hadn't been called.

      Note that the destructor is not virtual and so this class should not be
      used polymorphically.
   */
  ~Tag()
  {
    if ( !m_isClosed )
    {
      Close();
    }
  }

protected:
  /**
     This ctor is only used by RootTag: this ensures that it is impossible to
     construct a Tag without parent accidentally
   */
  Tag(const std::string& name, std::ostream& ostr)
    : m_ostr(ostr),
      m_nPrecision(0),
      m_name(name)
  {
    Init(0);
  }

private:
  template<typename T>
  std::string ToXMLValue(const T& value)
  {
    std::ostringstream ostrstr;
    ostrstr << value;
    return ostrstr.str();
  }

  std::string ToXMLValue(float value)
  {
    return String::FromCDouble(value, m_nPrecision > 0 ? m_nPrecision : -1);
  }

  std::string ToXMLValue(double value)
  {
    return String::FromCDouble(value, m_nPrecision > 0 ? m_nPrecision : -1);
  }

protected:
  /// where the output goes
  std::ostream& m_ostr;

  /// Specifies the number of digits to display in a floating-point number
  std::streamsize m_nPrecision;

private:
  // the name of this tag (needed to close it)
  std::string m_name;

  // the indent of this tag (in indent units, not spaces)
  unsigned m_indent;

  // did we have any contents (subtags, text, ...) already?
  bool m_isEmpty;

  // was the tag already closed?
  bool m_isClosed;

  // should we output this tag entirely on one line, without indentation?
  bool m_isTerse;

  // forbid assigning tags
  Tag& operator=(const Tag&);
};


/**
    RootTag is the root of XML document.

    Each XML document must have a single root and so creating a new RootTag
    creates a new XML document. All other tags in this document must be
    created with this root tag as parent, directly or indirectly (i.e.
    recursively).

    The direction of the XML output stream is also defined when creating the
    root tag: all the tags in this tree will share the stream specified in our
    constructor.
 */
class RootTag : public Tag
{
public:
  /**
      Start building XML document with the given root tag on output stream.

      When the root tag is destroyed, the XML document is finalized. This in
      particular means that the underlying stream object must have a greater
      life time than the root tag, i.e. it can only be destroyed when the root
      tag doesn't exist any more.

      @param name of the root tag of the XML document
      @param ostr the stream where the entire contents of this XML document is
                  sent
   */
  RootTag(const std::string& name, std::ostream& ostr)
    : Tag(name, ostr)
  {
    m_ostr << "<?xml version=\"1.0\"?>\n";

    BeginStartTag();
  }


   /**
      Start building XML document with the given root tag on output stream and
      a specific style sheet
      When the root tag is destroyed, the XML document is finalized. This in
      particular means that the underlying stream object must have a greater
      life time than the root tag, i.e. it can only be destroyed when the root
      tag doesn't exist any more.

      @param name of the root tag of the XML document
      @param ostr the stream where the entire contents of this XML document is
                  sent
      @param stylesheetname to use to transform xml to html

   */
  RootTag(const std::string& name, std::ostream& ostr,const std::string stylesheetname)
    : Tag(name, ostr)
  {
    m_ostr << "<?xml version=\"1.0\"?>\n";
    //<?xml-stylesheet type="text/xsl" href="test_style_sheet.xsl"?>
    m_ostr << "<?xml-stylesheet type=\"text/xsl\" href=\" "<<stylesheetname<<" \"?>\n";
    BeginStartTag();
  }

private:
  // there is no need for a copy ctor for RootTag as there is no need to
  // return it from functions -- it will be almost always created as a local
  // variable
  RootTag(const RootTag&);
  RootTag& operator=(const RootTag&);
};


template <class T>

// intel compiler doesn't need and like this work around specific for msvc
#ifndef __INTEL_COMPILER  
extern inline
#endif

Tag MakeTag(const T& value, Tag& tagParent)
{
  return value.Dump(tagParent);
}

template <class T>

// intel compiler doesn't need and like this work around specific for msvc
#ifndef __INTEL_COMPILER  
extern inline
#endif

Tag MakeNamedTag(const char *name, const T& value, Tag& parent)
{
  return value.Dump(name, parent);
}

} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_WRITE_H_


