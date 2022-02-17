/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/dividends.cpp
// Purpose:     Implementation of non inline Dividends methods
// Author:      Vadim Zeitlin
// Created:     30.07.03
// RCS-ID:      $Id: dividends.cpp,v 1.23 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <algorithm>
#include "ito33/afterstd.h"

#include "ito33/useexception.h"
#include "ito33/debug.h"

#include "ito33/finance/dividends.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/common.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace finance
{

// ----------------------------------------------------------------------------
// local functions
// ----------------------------------------------------------------------------

// used as predicate for adjacent_find()
static inline bool AreDatesEqual(const Dividend& div1, const Dividend& div2)
{
  return div1.date == div2.date;
}

// used with for_each()
static inline void CheckValidity(const Dividend& div)
{
  if ( (div.type == Dividend::Yield && div.value >= 1.) ||
          (div.type == Dividend::PseudoCash && div.pseudoYield >= 1.) )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("The dividend percentage value cannot be greater than 1.")
          );
  }

  if (div.type == Dividend::Cash && div.value <= 0.)
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("The dividend cash value must be positive.")
          );
}

// ============================================================================
// Dividends implementation
// ============================================================================

void Dividends::DoValidate()
{
  // NB: it is much more efficient to use the member sort() than std::sort()
  //     for std::list<>
  m_elements.sort();

  // check that there are no duplicates now
  const Elements::iterator begin = m_elements.begin(),
                           end = m_elements.end();
  if ( std::adjacent_find(begin, end, AreDatesEqual) != end )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Multiple dividends on same date are not allowed.")
          );
  }

  // finally check that the percent and pseudo cash dividends (if any) have
  // correct values
  std::for_each(begin, end, CheckValidity);

  m_bValidated = true;
}

void Dividends::Add(const Dividend& div)
{
  CheckValidity(div);

  m_elements.push_back(div);

  m_bValidated = false;
}

bool Dividends::HasCashBetween(Date dateBegin, Date dateEnd)
{
  ASSERT_MSG(dateBegin <= dateEnd, "Begin date must be before end date.");

  Dividends::Elements::iterator iter;

  for ( iter = m_elements.begin(); iter != m_elements.end(); ++iter )
  {
    if ( iter->date >= dateBegin && iter->date <= dateEnd )
      if ( iter->type != Dividend::Yield )
        return true;
  } 

  return false;
}

// ----------------------------------------------------------------------------
// XML stuff
// ----------------------------------------------------------------------------

using ito33::XML::Tag;

Tag
Dividend::Dump(Tag& tagParent) const
{
  ASSERT_MSG( type != PseudoCash,
                "pseudo cash dividends not supported here, TODO!" );

  Tag tagDividend(XML_TAG_DIVIDEND, tagParent);
  tagDividend.Element(XML_TAG_FINANCE_DATE)(date);
  tagDividend.Element(XML_TAG_FINANCE_VALUE)(value);
  
  tagDividend.Element(XML_TAG_FINANCE_TYPE)(type == Yield
      ? XML_VALUE_DIVIDEND_TYPE_YIELD
      : XML_VALUE_DIVIDEND_TYPE_CASH
  );

  return tagDividend;
}

Tag Dividends::Dump(Tag& tagParent) const
{ 
  Tag tagDividends(XML_TAG_FINANCE_DIVIDENDS, tagParent);

  const finance::Dividends::Elements& divs = GetAll();  
  finance::Dividends::Elements::const_iterator i;

  for ( i = divs.begin(); i != divs.end(); ++i )
  {
    tagDividends.Element(*i);
  }

  return tagDividends;
}

} // namespace finance

} // namespace ito33
