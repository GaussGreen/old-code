/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/xml/read.h
// Purpose:     Helpers for parsing XML
// Author:      Vadim Zeitlin
// Created:     2004-05-10
// RCS-ID:      $Id: read.h,v 1.10 2006/06/15 16:45:35 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/xml/read.h
    @brief Helper functions used by XML parsing code.

    The real job of XML parsing is done by xmlwrapp library but it provides
    only low level functions and we sometimes need higher level functions, such
    as defined in this file.
 */

#ifndef _ITO33_XML_READ_H_
#define _ITO33_XML_READ_H_

#include "ito33/error.h"
#include "ito33/date.h"
#include "ito33/gettext.h"
#include "ito33/string.h"

#include "ito33/enum_values_names.h"

#include <xmlwrapp/node.h>

extern const ito33::Error ITO33_BAD_PARAM, ITO33_BAD_PARAM,
                          ITO33_UNDEF_DATA, ITO33_BAD_DATA,
                          ITO33_BAD_DATE;

namespace ito33
{

namespace XML
{

/**
    Base class for all exceptions thrown by our code while parsing XML.

    Note that xmlwrapp itself also throws std::runtime_error which should be
    caught in addition to XML::Exception.
 */
class Exception : public ito33::Exception
{
public:
  /**
      Constructs the exception object.

      The parameter meanings are the same as for the base class constructor.

      @param error the error code
      @param msg the human readable error message (may be empty)
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  Exception(int error,
            const std::string& msg,
            const char *filename,
            size_t line,
            const char *function)
    : ito33::Exception(error, msg, filename, line, function)
  {
  }
};

/**
    Exception thrown when the expected node is not found in XML document.
 */
class MissingNodeException : public Exception
{
public:
  /**
      Constructs the exception thrown when an expected XML node is not found.

      This exception is thrown by GetNodeByName() when it fails to find an
      expected child tag.

      @param node the parent node
      @param name the name of child note we didn't find
      @param filename the name of the file where the exception occured
                      (usually just __FILE__)
      @param line the line number where it occured (usually __LINE__)
      @param function the name of the function where the exception occured
                      (__FUNCTION__ is unfortunately not yet supported by all
                      compilers so this argument is left empty for them)
   */
  MissingNodeException(const xml::node& node,
                       const char *name,
                       const char *filename,
                       size_t line,
                       const char *function)
    : Exception(ITO33_UNDEF_DATA,
                String::Printf(TRANS("Expected tag \"%s\" not found under "
                                     "tag \"%s\""), name, node.get_name()),
                filename,
                line,
                function)
  {
  }
};

/**
    Exception thrown when the node doesn't contain a value of expected type.
 */
class TypeMismatchException : public Exception
{
public:
  /**
      Constructs the exception object.

      The parameter meanings are the same as for the base class constructor.
   */
  TypeMismatchException(int error,
                        const std::string& msg,
                        const char *filename,
                        size_t line,
                        const char *function)
    : Exception(error, msg, filename, line, function)
  {
  }
};



/**
    Returns the child node with the given name or throw if not found.

    Throws if this node doesn't have child element with this name.

    @param node parent node
    @param name the name of the child node
    @return the first child node with the given name
 */
inline
xml::node GetNodeByName(const xml::node& node, const char *name)
{
  xml::node::const_iterator i = node.find(name);
  if ( i == node.end() )
  {
    // name is not exactly message but we can still (ab)use this macro here
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, name);
  }

  return *i;
}

/**
    Returns the child node with the given name or throw if not found.
    Search in the tree is done recursively.

    @param node parent node
    @param name the name of the child node
    @return the first child node with the given name or 
            parent if not found
 */
inline 
xml::node GetNodeRecursive(const xml::node &node, const char *name)
{

  xml::node::const_iterator i = node.find(name);
  if ( i != node.end() )
  {
    return *i;
  }

  for ( i = node.begin(); i!=node.end(); ++i)
  {
    xml::node tmp = GetNodeRecursive(*i,name);

    if ( strcmp(tmp.get_name(),name) == 0)
      return tmp;
  }
  
  return node;

}

/**
    Returns the child node with the given name or throw if not found.
    Search in the tree is done recursively.

    Throws if this node doesn't have any child element with this name.

    @param node parent node
    @param name the name of the child node
    @return the first child node with the given name
 */
inline 
xml::node GetNodeByNameRecursive(const xml::node &node, const char *name)
{

  xml::node nodeFound = GetNodeRecursive(node,name);

  if ( strcmp(nodeFound.get_name(),name) != 0 )
  {
    typedef MissingNodeException Exception;
    throw EXCEPTION_MSG(node, name);
  }
    
  return nodeFound;
}

/**
    Returns the text contents of the given node as date.

    Throws if the node doesn't contain date in ISO 8601 format.
 */
inline
Date GetDateFromNode(const xml::node& node)
{
  Date date;
  if ( !date.Parse(node.get_content(), "%Y-%m-%d") )
  {
    typedef TypeMismatchException Exception;

    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATE,
            String::Printf
            (
              TRANS("Node \"%s\" contains \"%s\" which is not a valid date "
                    "in ISO 8601 format."),
              node.get_name(),
              node.get_content()
            )
          );
  }

  return date;
}

/**
    Returns the date stored in the child with the given name.

    An exception is thrown if the child node is not found or doesn't contain a
    valid date.

    @sa GetDateFromNode(), GetNodeByName()
 */
inline
Date GetDateFromName(const xml::node& node, const char *name)
{
  return GetDateFromNode(GetNodeByName(node, name));
}

/**
    Returns the text contents of the given node as a floating point number.

    Throws if the node doesn't contain a string representation of a double.
 */
inline
double GetDoubleFromNode(const xml::node& node)
{
  double d;
  if ( !String::ToCDouble(node.get_content(), &d) )
  {
    typedef TypeMismatchException Exception;

    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            String::Printf
            (
              TRANS("Node \"%s\" contains \"%s\" which is not a valid "
                    "floating point number."),
              node.get_name(),
              node.get_content()
            )
          );
  }

  return d;
}

/**
    Returns the double stored in the child with the given name.

    An exception is thrown if the child node is not found or doesn't contain a
    valid double value.

    @sa GetDoubleFromNode(), GetNodeByName()
 */
inline
double GetDoubleFromName(const xml::node& node, const char *name)
{
  return GetDoubleFromNode(GetNodeByName(node, name));
}

/**
    Returns the text contents of the given node as number(of type long).

    Throws if the node doesn't contain a string representation of
    a number(of type long).
 */
inline
long GetLongFromNode(const xml::node& node)
{
  long i;
  if ( !String::ToLong(node.get_content(), &i) )
  {
    typedef TypeMismatchException Exception;

    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            String::Printf
            (
              TRANS("Node \"%s\" contains \"%s\" which is not a valid "
                    "number."),
              node.get_name(),
              node.get_content()
            )
          );
  }

  return i;
}

/**
    Returns the number(of type long) stored in the child with the given name.

    An exception is thrown if the child node is not found or doesn't contain a
    valid number(of type long).

    @sa GetLongFromNode(), GetNodeByName()
 */
inline
long GetLongFromName(const xml::node& node, const char *name)
{
  return GetLongFromNode(GetNodeByName(node, name));
}

/**
    Returns the text contents of the given node as boolean.

    Throws if the node doesn't contain a string representation of
    a number(of type long).
 */
inline
bool GetBoolFromNode(const xml::node& node)
{
  return GetLongFromNode(node) != 0;
}

/**
    Returns the boolean stored in the child with the given name.

    An exception is thrown if the child node is not found or doesn't contain a
    valid number(of type long).

    @sa GetBoolFromNode(), GetNodeByName()
 */
inline
bool GetBoolFromName(const xml::node& node, const char *name)
{
  return GetBoolFromNode(GetNodeByName(node, name));
}

/**
    Returns the enum value corresponding to the text value of the given node.

    The mappings between strings in the XML document and enum values is defined
    by the values parameter. If no matching key is found, an exception is
    thrown.

    @param node to translate to enum
    @param count number of elements in values array
    @param values the mapping between strings and enum values
    @return the enum value corresponding to the text content of the node
 */
template <typename T>
inline
T GetEnumFromNode(const xml::node& node,
                  size_t count,
                  const EnumValuesNames<T> values[])
{
  T value;
  if(!RestoreEnumValueFromName(node.get_content(), value, count, values))
  {
    typedef TypeMismatchException Exception;

    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            String::Printf
            (
              TRANS("Node \"%s\" contains \"%s\" which is not a valid "
                    "enum element."),
              node.get_name(),
              node.get_content()
            )
          );
  }

  return value;
}

/**
    Returns the enum value corresponding to the text value of the given child
    node.

    @param node the parent node
    @param name the name of the child node
    @param count number of elements in values array
    @param values the mapping between strings and enum values
    @return the enum value corresponding to the text content of the node
 */
template <typename T>
inline
T GetEnumFromName(const xml::node& node,
                  const char *name,
                  size_t count,
                  const EnumValuesNames<T> values[])
{
  return GetEnumFromNode(GetNodeByName(node, name), count, values);
}


} // namespace XML

} // namespace ito33

#endif // _ITO33_XML_READ_H_
