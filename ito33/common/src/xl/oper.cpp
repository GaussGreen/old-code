/////////////////////////////////////////////////////////////////////////////
// Name:        src/XL/oper.cpp
// Purpose:     implementation of XL::Oper class
// Author:      Vadim Zeitlin
// Created:     2006-05-05
// RCS-ID:      $Id: oper.cpp,v 1.1 2006/05/05 01:58:57 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/XL/oper.h"

namespace ito33
{

namespace XL
{

// ============================================================================
// Oper implementation
// ============================================================================

/* static */
std::string Oper::GetErrorName(WORD err)
{
  std::string s;

  switch ( err )
  {
#define CASE_ERR(e) case e: s = #e; break

    CASE_ERR(xlerrNull);
    CASE_ERR(xlerrDiv0);
    CASE_ERR(xlerrValue);
    CASE_ERR(xlerrRef);
    CASE_ERR(xlerrName);
    CASE_ERR(xlerrNum);
    CASE_ERR(xlerrNA);

#undef CASE_ERR

    default:
      s = String::Printf("unknown (%d)", err);
  }

  return s;
}

std::string Oper::Dump() const
{
  std::string s;

  switch ( GetType() )
  {
    case Type_Num:
      s = String::Printf("num = %g", val.num);
      break;

    case Type_Str:
      s = String::Printf("str = \"%s\"",
                         std::string(val.str + 1, *val.str).c_str());
      break;

    case Type_Bool:
      s = String::Printf("bool = %s", val.boolean ? "TRUE" : "FALSE");
      break;

    case Type_Ref:
      s = String::Printf("ref: sheet %d, %d ranges",
                         val.mref.idSheet, val.mref.lpmref->count);
      break;

    case Type_Err:
      s = String::Printf("err = %s", GetErrorName(val.err).c_str());
      break;

    case Type_Flow:
      s = String::Printf("flow = ...");
      break;

    case Type_Multi:
      s = String::Printf("multi %d*%d", val.array.rows, val.array.columns);
      break;

    case Type_Missing:
      s = "missing";
      break;

    case Type_Nil:
      s = "nil";
      break;

    case Type_SRef:
      {
        const XLREF& ref = val.sref.ref;
        if ( ref.rwLast == ref.rwFirst && ref.colLast == ref.colFirst )
        {
          // just a single cell
          s = String::Printf("sref = (%c, %d)",
                             'A' + ref.colFirst, 1 + ref.rwFirst);
        }
        else // real range
        {
          s = String::Printf("sref = (%c, %d)-(%c, %d)",
                             'A' + ref.colFirst, 1 + ref.rwFirst,
                             'A' + ref.colLast, 1 + ref.rwLast);
        }
      }
      break;

    case Type_Int:
      s = String::Printf("int = %d", val.w);
      break;

    default:
      s = String::Printf("unknown type = %d", GetType());
  }

  return "{" + s + "}";
}

const Oper *Oper::GetAsArray(unsigned& w, unsigned& h) const
{
  if ( GetType() != Type_Multi )
    return NULL;

  w = val.array.columns;
  h = val.array.rows;
  const Oper * const data = static_cast<Oper *>(val.array.lparray);

  // A range is often given with blank cells beneath it because it allows to
  // add values to it easily, without changing the formulae using it. So
  // determine the effective end of range: it stops at the first blank line
  for ( unsigned n = 0; n < h; n++ )
  {
    const Oper *p = data + n*w;
    const Oper * const end = p + w;
    while ( p < end )
    {
      if ( !p->IsMissing() )
        break;

      ++p;
    }

    if ( p == end )
    {
      // this row is empty so stop at it
      h = n;
      break;
    }
  }

  return data;
}
} // namespace XL

} // namespace ito33

