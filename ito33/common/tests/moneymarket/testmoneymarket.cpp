/////////////////////////////////////////////////////////////////////////////
// Name:        testmoneymarket.cpp
// Purpose:     Acceptance test for moneymarkets 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testmoneymarket.cpp,v 1.5 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/sharedptr.h"

#include "ito33/finance/moneymarket.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testmoneymarket.h"

#include "ito33/tests/utilexml.h"

using namespace ito33;
using namespace ito33::finance;

void MoneyMarketTest::NoYieldCurve()
{
  shared_ptr<YieldCurve> pYC;

  shared_ptr<Numeraire> pCurrency( new Numeraire("USD") );
  MoneyMarket moneyMarket(pCurrency, pYC);

} //MoneyMarketTest::NoYieldCurve()


void MoneyMarketTest::Dump()
{

  shared_ptr<finance::YieldCurve> pYc( new finance::YieldCurveFlat(.02) );
  
  shared_ptr<Numeraire> pCurrency( new Numeraire("USD") );

  MoneyMarket moneyMarket(pCurrency, pYc);


  std::ostringstream oss;

    ExpectedXML expected(oss,
    "<?xml version=\"1.0\"?>"
    "<root>\n"
    "<money_market>\n"
    "<yield_curve>\n"
    "<yield_curve_flat>\n"
    "<flat>0.02</flat>\n"
    "</yield_curve_flat>\n"
    "</yield_curve>\n"
    "<currency>USD</currency>\n"
    "</money_market>\n"
    "</root>\n");

  ito33::XML::RootTag root("root",oss);

  moneyMarket.Dump(root);

}//MoneyMarketTest::Dump()
