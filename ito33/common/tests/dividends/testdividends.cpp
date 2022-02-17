/////////////////////////////////////////////////////////////////////////////
// Name:        test/dividends/main.cpp
// Purpose:     Unit test for Dividends class
// Author:      Vadim Zeitlin
// Created:     30.07.03
// RCS-ID:      $Id: testdividends.cpp,v 1.3 2004/11/23 10:45:39 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/exception.h"
#include "ito33/finance/dividends.h"

#include "ito33/debug.h"

#include "ito33/list.h"
#include "ito33/vector.h"

#include "ito33/cppunit.h"

#include "ito33/tests/testdividends.h"

using namespace ito33;

using ito33::finance::Dividend;
using ito33::finance::Dividends;


void
DivsTestCase::DoCheckElement(const Dividend& div,
                             Dividend::Type type,
                             int date,
                             double value,
                             double pseudoYield)
{
  CPPUNIT_ASSERT( div.type == type );
  CPPUNIT_ASSERT( div.date == date );
  CPPUNIT_ASSERT( div.value == value );
  CPPUNIT_ASSERT( div.pseudoYield == pseudoYield );
}

void DivsTestCase::DoCheckElements()
{
  const Dividends::Elements& elements = m_divs.GetAll();

  CPPUNIT_ASSERT( elements.size() == 4 );

  Dividends::Elements::const_iterator pDiv = elements.begin();

  DoCheckElement(*pDiv++, Dividend::Cash, 7, 12);
  DoCheckElement(*pDiv++, Dividend::Yield, 30, 0.05);
  DoCheckElement(*pDiv++, Dividend::Cash, 365, 15);
  DoCheckElement(*pDiv++, Dividend::PseudoCash, 730, 15, 0.2);

  CPPUNIT_ASSERT( pDiv == elements.end() );
}


void DivsTestCase::Add()
{
  m_divs.AddCash(7, 12);
  m_divs.AddYield(30, 0.05);
  m_divs.Add(Dividend::Cash, 365, 15);
  m_divs.AddPseudoCash(730, 15, 0.2);

  DoCheckElements();
}

void DivsTestCase::AddReverse()
{
  m_divs.AddPseudoCash(730, 15, 0.2);
  m_divs.Add(Dividend::Cash, 365, 15);
  m_divs.AddYield(30, 0.05);
  m_divs.AddCash(7, 12);

  DoCheckElements();
}

void DivsTestCase::AddDuplicates()
{
  // having 2 dividends on the same day is not allowed
  m_divs.AddCash(30, 12);
  m_divs.AddYield(30, 0.1);

  m_divs.Validate();
}

void DivsTestCase::SetCash()
{
  static const int dates[] = { 7, 30, 365 };
  static const double values[] = { 1, 12, 15 };

  static const int nDivs = SIZEOF(dates);

  m_divs.SetCash(dates, dates + nDivs, values);

  const Dividends::Elements& elements = m_divs.GetAll();

  CPPUNIT_ASSERT( elements.size() == nDivs );

  Dividends::Elements::const_iterator pDiv = elements.begin();
  for ( size_t n = 0; n < nDivs; ++n )
  {
    DoCheckElement(*pDiv++, Dividend::Cash, dates[n], values[n]);
  }

  CPPUNIT_ASSERT( pDiv == elements.end() );
}

void DivsTestCase::SetYield()
{
  static const int dates[] = { 7, 30, 365 };
  static const double values[] = { 0.01, 0.12, 0.15 };

  static const int nDivs = SIZEOF(dates);

  m_divs.SetYield(dates, dates + nDivs, values);

  const Dividends::Elements& elements = m_divs.GetAll();

  CPPUNIT_ASSERT( elements.size() == nDivs );

  Dividends::Elements::const_iterator pDiv = elements.begin();
  for ( size_t n = 0; n < nDivs; ++n )
  {
    DoCheckElement(*pDiv++, Dividend::Yield, dates[n], values[n]);
  }

  CPPUNIT_ASSERT( pDiv == elements.end() );
}

void DivsTestCase::NegativeCash()
{
  Dividends div;
  div.AddCash(10000, -0.2);
}

void DivsTestCase::LargeYield()
{
  // 1.2 is not a valid percentage value
  static const int dates[] = { 7 };
  static const double values[] = { 1.2 };

  m_divs.SetYield(dates, dates + SIZEOF(dates), values);
}

void DivsTestCase::LargePseudoYield()
{
  // 1.2 is not a valid percentage value
  m_divs.AddPseudoCash(7, 15, 1.2);
  m_divs.Validate();
}

void DivsTestCase::Clear()
{
  m_divs.AddCash(7, 12);

  CPPUNIT_ASSERT( !m_divs.GetAll().empty() );

  m_divs.Clear();
  CPPUNIT_ASSERT( m_divs.GetAll().empty() );
}
