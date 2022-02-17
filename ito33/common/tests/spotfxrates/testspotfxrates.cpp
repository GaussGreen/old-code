/////////////////////////////////////////////////////////////////////////////
// Name:        tests/spotfxrates/testspotfxrates.cpp
// Purpose:     Unit test for SpotFXRates
// Author:      Wang
// Created:     2004/09/02
// RCS-ID:      $Id: testspotfxrates.cpp,v 1.6 2006/08/19 23:21:44 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <cmath>
#include "ito33/beforestd.h"

#include "ito33/finance/spotfxrates.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testspotfxrates.h"

using namespace ito33;
using namespace ito33::finance;


void SpotFXRatesTest::GetWithoutSet()
{
  SpotFXRates rates;

  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );
  shared_ptr<Numeraire> pCurrencyEUR( new Numeraire("EUR") );

  rates.GetFXRate(pCurrencyUSD, pCurrencyEUR ); 
}

void SpotFXRatesTest::OneForSameNumeraires()
{
  SpotFXRates rates;

  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );

  CPPUNIT_ASSERT_EQUAL(rates.GetFXRate(pCurrencyUSD, pCurrencyUSD), 1.);
}

void SpotFXRatesTest::GetWhileSetInversely()
{
  SpotFXRates rates;
  
  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );
  shared_ptr<Numeraire> pCurrencyEUR( new Numeraire("EUR") );

  rates.SetFXRate(pCurrencyUSD, pCurrencyEUR, 0.9);

  CPPUNIT_ASSERT_EQUAL(rates.GetFXRate(pCurrencyEUR, pCurrencyUSD), 1 / 0.9);
}

void SpotFXRatesTest::GetAfterUpdate()
{
  SpotFXRates rates;

  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );
  shared_ptr<Numeraire> pCurrencyEUR( new Numeraire("EUR") );

  rates.SetFXRate(pCurrencyUSD, pCurrencyEUR, 0.9);

  rates.SetFXRate(pCurrencyUSD, pCurrencyEUR, 0.8);

  CPPUNIT_ASSERT_EQUAL(rates.GetFXRate(pCurrencyUSD, pCurrencyEUR), 0.8);
}

void SpotFXRatesTest::MultiValues()
{
  SpotFXRates rates;
  
  shared_ptr<Numeraire> pCurrencyUSD( new Numeraire("USD") );
  shared_ptr<Numeraire> pCurrencyEUR( new Numeraire("EUR") );
  shared_ptr<Numeraire> pCurrencyJPY( new Numeraire("JPY") );

  rates.SetFXRate(pCurrencyUSD, pCurrencyEUR, 0.8);

  rates.SetFXRate(pCurrencyUSD, pCurrencyJPY, 100);

  CPPUNIT_ASSERT_EQUAL(rates.GetFXRate(pCurrencyUSD, pCurrencyEUR), 0.8);

  CPPUNIT_ASSERT_EQUAL(rates.GetFXRate(pCurrencyUSD, pCurrencyJPY), 100.);
}
