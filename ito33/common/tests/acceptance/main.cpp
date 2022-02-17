/////////////////////////////////////////////////////////////////////////////
// Name:        tests/acceptance/main.cpp
// Purpose:     Acceptance testing for common code
// Author:      Testing
// Created:     25/06/2004
// RCS-ID:      $Id: main.cpp,v 1.21 2006/02/28 16:13:48 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/cppunit.h"

#include "ito33/tests/testautoptr.h"
#include "ito33/tests/testbisecnewt.h"

#include "ito33/tests/testbondlike.h"
#include "ito33/tests/testyieldtomaturity.h"
#include "ito33/tests/testyieldtoPut.h"
#include "ito33/tests/testnewshare.h"

#include "ito33/tests/testBSut.h"
#include "ito33/tests/testcashflowstream.h"
#include "ito33/tests/const_constraint_case.h"
#include "ito33/tests/testbinarysearch.h"
#include "ito33/tests/testdividends.h"
#include "ito33/tests/testequity.h"
#include "ito33/tests/testissuer.h"
#include "ito33/tests/testmoneymarket.h"
#include "ito33/tests/testerror.h"
#include "ito33/tests/testeventmanager.h"
#include "ito33/tests/testinterput.h"
#include "ito33/tests/testnls.h"
#include "ito33/tests/testnd.h"
#include "ito33/tests/testsharedptr.h"
#include "ito33/tests/testsurfaces.h"
#include "ito33/tests/testtimemesh.h"
#include "ito33/tests/testxmlwrite.h"
#include "ito33/tests/testycurve.h"
#include "ito33/tests/testconversionpricereset.h"
#include "ito33/tests/testresetconversionschedule.h"
#include "ito33/tests/testresets.h"
#include "ito33/tests/testoption.h"
#include "ito33/tests/testcds.h"
#include "ito33/tests/testparbond.h"
#include "ito33/tests/testnewton2d.h"
#include "ito33/tests/testqpminimizer.h"
#include "ito33/tests/testsharedependentconversion.h"
#include "ito33/tests/testattachedwarrantconvertiblebond.h"
#include "ito33/tests/testasianoption.h"
#include "ito33/tests/testmandatory.h"
#include "ito33/tests/testfloatingrates.h"
#include "ito33/tests/testsa.h"
#include "ito33/tests/testspotfxrates.h"

int main()
{
  CppUnit::TextUi::TestRunner runner;

  runner.addTest(AutoPtrTestCase::suite());
  runner.addTest(BinarySearchTest::suite());
  runner.addTest(BSutileTest::suite());
  runner.addTest(ConstConstraintCase::suite());
  runner.addTest(DivsTestCase::suite());
  runner.addTest(IssuerTest::suite());
  runner.addTest(EquityTest::suite());
  runner.addTest(MoneyMarketTest::suite());
  runner.addTest(ErrorCodeTestCase::suite());
  runner.addTest(EventManagerTest::suite());
  runner.addTest(SharedPtrTestCase::suite());
  runner.addTest(DomainFixedSpaceMeshTest::suite());
  runner.addTest(SurfaceGeneralTest::suite());
  runner.addTest(TimeMeshTest::suite());
  runner.addTest(XMLWriteTestCase::suite());
  runner.addTest(YCurveFlatTestCase::suite());
  runner.addTest(YCurveAnnuallyCompoundedTestCase::suite());
  runner.addTest(CashFlowStreamTest::suite());
  runner.addTest(FloatingRatesTest::suite());
  runner.addTest(SpotFXRatesTest::suite());

  //math 
  runner.addTest(NLSTest::suite());
  runner.addTest(NormalDistTest::suite());
  runner.addTest(InterpUtTest::suite());
  runner.addTest(BisecNewtTest::suite());
  runner.addTest(Newton2DTest::suite());
  runner.addTest(QPMinimizerTest::suite());
  runner.addTest(SATest::suite());
  runner.addTest(ASATest::suite());

  //resets
  runner.addTest(ConversionPriceResetTest::suite());
  runner.addTest(ConversionScheduleResetTest::suite());
  runner.addTest(ResetTest::suite());

  //option
  runner.addTest(OptionTest::suite());

  //cds
  runner.addTest(CDSTest::suite());
  
  //parbond
  runner.addTest(ParBondTest::suite());

  //bondlike
  runner.addTest(BondLikeTest::suite());
  runner.addTest(YieldToMaturityTest::suite());
  runner.addTest(YieldToPutTest::suite());
  runner.addTest(NewShareTest::suite());

  //Shared dependent conversion
  runner.addTest(ShareDependentConversionTest::suite());

  //attached warrant convertible bond
   runner.addTest(AttachedWarrantConvertibleBondTest::suite());

  //asian option
  runner.addTest(AsianOptionTest::suite());
  runner.addTest(CurranTest::suite());

  //mandatory
  runner.addTest(PepsAveragingPeriodTest::suite());

  return runner.run("") ? 0 : 1;
}
