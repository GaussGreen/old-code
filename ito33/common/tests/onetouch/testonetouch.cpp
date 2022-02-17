/////////////////////////////////////////////////////////////////////////////
// Name:        common/tests/onetouch/testonetouch.cpp
// Purpose:     Acceptance test for one touch 
// Created:     2005/12/08
// RCS-ID:      $Id: testonetouch.cpp,v 1.7 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/sharedptr.h"

#include "ito33/finance/exoticoption/onetouch.h"
#include "ito33/finance/exoticoption/fxonetouch.h"
#include "ito33/finance/exoticoption/onetouchutils.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/dividends.h"

#include "ito33/tests/utilexml.h"

#include "ito33/tests/testonetouch.h"

#include "ito33/xml/write.h"

#include "ito33/link.h"

using namespace ito33;
using namespace ito33::finance;

void OneTouchTest::Delta2Barrier()
{
  Date valuationDate(2002, Date::Jul, 1);

  double dSpot = 0.9;

  shared_ptr<Numeraire> pNumeraire(new Numeraire("EUR"));
  shared_ptr<Equity> pEquity(new Equity(dSpot, pNumeraire));

  shared_ptr<Dividends> pDividends( new Dividends() );
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurve> pyf( new YieldCurveFlat(0.00) );
    
  pEquity->SetBorrowCurve(pyf);
    
  double dRate = 0.00;
  shared_ptr<YieldCurve> pyc( new YieldCurveFlat(dRate) );

  shared_ptr<RateData> pRateData(new RateData);
  pRateData->SetYieldCurve(pNumeraire, pyc);
  
  shared_ptr<SessionData> 
    pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  FXOneTouch oneTouch(Date(2003, Date::Jul, 1), 0.2, Barrier_UpAndOut, 0.5);

  oneTouch.SetSessionData(pSessionData);

  double dBarrier = oneTouch.GetBarrier();
  
  CPPUNIT_ASSERT( dBarrier > dSpot );

  FXOneTouch oneTouch2(Date(2003, Date::Jul, 1), 0.2, Barrier_DownAndOut, 0.5);

  oneTouch2.SetSessionData(pSessionData);

  double dBarrier2 = oneTouch2.GetBarrier();
  
  CPPUNIT_ASSERT( dBarrier2 < dSpot );
}
