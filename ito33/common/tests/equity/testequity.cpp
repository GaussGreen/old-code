/////////////////////////////////////////////////////////////////////////////
// Name:        testoption.cpp
// Purpose:     Acceptance test for option 
// Author:      Yann d'Halluin
// Created:     17/11/2004
// RCS-ID:      $Id: testequity.cpp,v 1.5 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2004 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/sharedptr.h"

#include "ito33/finance/equity.h"
#include "ito33/finance/numeraire.h"

#include "ito33/tests/testequity.h"

using namespace ito33;
using namespace ito33::finance;

void EquityTest::NegativeSpotSharePrice()
{
 Equity equity(1, shared_ptr<Numeraire>(new Numeraire("USD") ) );

 equity.SetSpotSharePrice(-100);

} //EquityTest::NegativeSpotSharePrice()
