/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/tests/hero.cpp
// Purpose:     Implementation of HERO tests
// Created:     2005/10/01
// RCS-ID:      $Id: hero.cpp,v 1.8 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <iostream>
#include <string>
#include <list>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/derivatives.h"

#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/hedgingdata.h"

#include "hg/xml/hedgereader.h"

#include "hg/tests/hero.h"

namespace ito33
{

namespace hg
{

namespace HeroTest
{

void HeroTest::Setup(const std::string& sFileName)
{

  // Read the hedging data from the specified file
  hg::XML::HedgeReader reader( sFileName.c_str() );
    
  reader.ReadTheoreticalModel(m_pModel);
      
  m_pTarget = reader.ReadTarget();

  CPPUNIT_ASSERT( m_pTarget != NULL );
    
  m_pHedgeInstruments = reader.ReadHedgeInstruments();

  CPPUNIT_ASSERT( m_pHedgeInstruments != NULL );
}


void HeroTest::ZeroHeroTest()
{

  // One regime, no jumps
  std::string sFileName("xmlfiles/hero01.xml");
  Setup(sFileName);

  shared_ptr<hg::HedgeOutput> 
    pHedgeOutput( m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments) );

  // Get the hero
  double dHero = pHedgeOutput->GetHERO();

  // if no jumps, should be able to hedge exactly
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dHero, 1.e-6);


  // Hedge with the target contract
  std::string sFileName2("xmlfiles/hero02.xml");
  Setup(sFileName2);

  pHedgeOutput = m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments);

  // Get the hero
  dHero = pHedgeOutput->GetHERO();

  // should be able to hedge exactly with the target contract
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dHero, 1.e-6);

}


void HeroTest::AddHedgingContractsTest()
{

  // Two regimes, zero hedge contracts to start
  std::string sFileName("xmlfiles/hero03.xml");
  Setup(sFileName);

  shared_ptr<hg::HedgeOutput> 
    pHedgeOutput( m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments) );

  // Get the initial hero
  double dOldHero = pHedgeOutput->GetHERO();

  shared_ptr<finance::SessionData> pSessionData = m_pTarget->GetSessionData();

  double dSpot = pSessionData->GetSpotSharePrice();
  Date valuationDate = pSessionData->GetValuationDate();
  Date maturityDate = m_pTarget->GetMaturityDate();
  
  
  const double pdStrikeDiffs[] = 
    {-16.0, 16.0, -8.0, 8.0, -4.0, 4.0, -6.0, 6.0, 12.0};

  //size_t nNbStrikes = sizeof(pdStrikeDiffs) / sizeof(double);
  size_t nNbStrikes = 2;
  if (m_pModel->GetUnderlyingProcess()->GetNbRegimes() == 2)
    nNbStrikes = 6;
  else
    nNbStrikes = 9;

  shared_ptr<finance::Option> opt;
  for (size_t nIdx = 0; nIdx < nNbStrikes; nIdx++)
  {
    opt = make_ptr( new finance::Option
                        (dSpot + pdStrikeDiffs[nIdx], maturityDate, 
                         finance::Option_Put, finance::ExerciseType_European) );

    opt->SetSessionData(pSessionData);

    m_pHedgeInstruments->Add(opt);

    // Compute hero with extra hedge contract
    pHedgeOutput = m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments);

    double dNewHero = pHedgeOutput->GetHERO();

    std::vector< shared_ptr<HedgeRatioData> > ppHedgeRatioData 
      = pHedgeOutput->GetHedgeRatioData();
    //size_t nn = ppHedgeRatioData.size();
    
    //std::cout << nIdx << std::endl;
    //std::cout << "Old hero = " << dOldHero 
    //          << ", New hero = " << dNewHero << std::endl;
    
    CPPUNIT_ASSERT(dNewHero <= dOldHero);
    dOldHero = dNewHero;

    // If we get lucky and get exact match early
    if (dOldHero < 1.e-6)
      break;
  }

  // When we have enough hedge contracts to cover the jumps, the hero 
  // should be zero
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dOldHero, 1.e-6);

}


void HeroTest::CompareSVDandNAGTest()
{

  // Make a list of xml files to test
  std::list<std::string> fileList;
  
  fileList.push_back("xmlfiles/hero01.xml");
  fileList.push_back("xmlfiles/hero02.xml");
  fileList.push_back("xmlfiles/hero03.xml");
  fileList.push_back("xmlfiles/hero04.xml");

  std::list<std::string>::const_iterator iterFileList;

  // Calculate hero by both SVD and NAG
  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {

    Setup( *iterFileList );

    // Compute using SVD (the default)
    shared_ptr<hg::HedgeOutput> 
      pHedgeOutput( m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments) );

    double dSVDHero = pHedgeOutput->GetHERO();

    // Compute using NAG
    HedgingData::bForceNAG = true;

    pHedgeOutput = m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments); 

    double dNAGHero = pHedgeOutput->GetHERO();

    HedgingData::bForceNAG = false;

    // Compare relative error
    double dError = fabs( dSVDHero - dNAGHero );
    if ( fabs(dNAGHero) > 1.e-6 )
      dError /= fabs(dNAGHero);

    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0, dError, 1.e-2);

    //std::cout.precision(12);
    //std::cout << "nag HERO = " << dNAGHero << ", svd hero = " << dSVDHero << std::endl;    

  } // loop over file list

}


void HeroTest::GenericTests()
{

  // Make a list of xml files to test
  std::list<std::string> fileList;
  
  fileList.push_back("xmlfiles/hero01.xml");
  fileList.push_back("xmlfiles/hero02.xml");
  fileList.push_back("xmlfiles/hero03.xml");

  std::list<std::string>::const_iterator iterFileList;

  // Run the hedge code
  for (iterFileList = fileList.begin();
       iterFileList != fileList.end();
       ++iterFileList)
  {
    Setup( *iterFileList );

    shared_ptr<hg::HedgeOutput> 
      pHedgeOutput( m_pModel->ComputeHERO(*m_pTarget, *m_pHedgeInstruments) );

    double dHero = pHedgeOutput->GetHERO();

    CPPUNIT_ASSERT( dHero >= 0.0 );
  }

}

/*
// Fix this test when path-dependent contracts are supported in HG
void HeroTest::MissingDataTest()
{

  // Two regimes, zero hedge contracts to start
  std::string sFileName("xmlfiles/hero03.xml");
  Setup(sFileName);

  // Add path dependent hedge contract with insufficient surface data
  Date maturityDate = m_pTarget->GetMaturityDate();
  maturityDate.AddDays(-1);
  shared_ptr<finance::SessionData> pSessionData = m_pTarget->GetSessionData();
  double dSpot = pSession->GetSpotSharePrice();

  shared_ptr<finance::Option> opt;
  opt = new finance::Option(dSpot, maturityDate, finance::Option_Put, 
                            finance::ExerciseType_European);

  opt->SetSessionData(pSessionData);

  m_pHedgeInstruments->Add(opt);

  // Should throw missing data exception
  shared_ptr<hg::HedgeOutput> pHedgeOutput =  
    m_pModel->Hedge(*m_pTarget, *m_pHedgeInstruments, m_bComputeHERO);

}
*/

} // namespace HeroTest

} // namespace ito33

} // namespace hg
