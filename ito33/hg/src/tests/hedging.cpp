/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/hedging.cpp
// Purpose:     Implementation of hedging tests
// Created:     2005/10/01
// RCS-ID:      $Id: hedging.cpp,v 1.6 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/hedgeratiodata.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/xml/hedgereader.h"

#include "hg/tests/hedging.h"

namespace ito33
{

namespace hg
{

namespace HedgingTest
{

  void HedgingTest::Setup(const std::string& sFileName)
{

  // Read the hedging data from the specified file
  hg::XML::HedgeReader reader( sFileName.c_str() );
    
  reader.ReadTheoreticalModel(m_pModel);
      
  m_pTarget = reader.ReadTarget();

  CPPUNIT_ASSERT( m_pTarget != NULL );
    
  m_pHedgeInstruments = reader.ReadHedgeInstruments();

  CPPUNIT_ASSERT( m_pHedgeInstruments != NULL );
}


void HedgingTest::BSDeltaTest()
{

  // One regime, no jumps
  std::string sFileName("xmlfiles/hedging01.xml");
  Setup(sFileName);

  shared_ptr<hg::HedgeOutput>
    pHedgeOutput( m_pModel->Hedge(*m_pTarget, *m_pHedgeInstruments) );

  // Get the delta
  shared_ptr<finance::ModelOutput> pTargetOutput = 
    pHedgeOutput->GetTargetModelOutput();

  double dDelta = pTargetOutput->GetDelta();

  // Get the underlying ratio
  double dUnderlyingRatio = pHedgeOutput->GetUnderlyingHedgeRatio();

  // underlying ratio must be the delta if no jumps, no default
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dDelta, dUnderlyingRatio, 1.e-8);

}


void HedgingTest::HedgeWithTargetTest()
{

  // The target contract is also the second hedge contract
  std::string sFileName("xmlfiles/hedging02.xml");
  Setup(sFileName);

  shared_ptr<hg::HedgeOutput>
    pHedgeOutput( m_pModel->Hedge(*m_pTarget, *m_pHedgeInstruments) );

  // Get the hedge ratios
  double dUnderlyingRatio = pHedgeOutput->GetUnderlyingHedgeRatio();
  std::vector< shared_ptr<HedgeRatioData> > ppHedgeRatioData
    = pHedgeOutput->GetHedgeRatioData();  

  // All ratios should be zero, except the second ratio which is 1
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dUnderlyingRatio, 1.e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, ppHedgeRatioData[0]->GetRatio(), 1.e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0, ppHedgeRatioData[1]->GetRatio(), 1.e-8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, ppHedgeRatioData[2]->GetRatio(), 1.e-8);

}


void HedgingTest::GenericTests()
{

  // Just see if the code runs for these examples.

  // Make a list of xml files to test
  std::list<std::string> fileList;
  
  fileList.push_back("xmlfiles/hedging01.xml");
  fileList.push_back("xmlfiles/hedging02.xml");
  fileList.push_back("xmlfiles/hedging03.xml");

  std::list<std::string>::const_iterator iterFileList;

  // Run the hedge code
  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {

    // Read the file and run the example
    Setup( *iterFileList );

    shared_ptr<hg::HedgeOutput>
      pHedgeOutput( m_pModel->Hedge(*m_pTarget, *m_pHedgeInstruments) );
/*        
    // Output hedge ratios
    double dUnderlyingRatio = pHedgeOutput->GetUnderlyingHedgeRatio();
    std::cout << "Underlying ratio = " << dUnderlyingRatio << std::endl;

    std::vector< shared_ptr<HedgeRatioData> > ppHedgeRatioData
      = pHedgeOutput->GetHedgeRatioData();  
    for (size_t nIdx = 0; nIdx < pdHedgeRatios.size(); nIdx++)
    {
      double dRatio = ppHedgeRatioData[nIdx]->GetRatio()
      std::cout << " Ratio " << nIdx << " = " << dRatio 
                << std::endl;
    }
    std::cout << std::endl;
*/
  }

}


} // namespace HedgingTest

} // namespace ito33

} // namespace hg
