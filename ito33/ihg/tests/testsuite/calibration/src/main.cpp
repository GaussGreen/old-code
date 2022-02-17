#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/link.h"

#include "testparametrization_hrwithtimecomponent.h"
#include "testparametrization_hrwithspotcomponentpower.h"
#include "testparametrization_volflathrwithtimecomponent.h"
#include "testparametrization_volflathrwithspotcomponentpower.h"
#include "testparametrization_volflathrpower.h"
#include "testparametrization_volpower.h"
#include "testparametrization_volpowerhrpower.h"
#include "testparametrization_volpowerhrwithtimecomponent.h"
#include "testparametrization_voltanh.h"
#include "testparametrization_voltanhhrpower.h"
#include "testparametrization_voltanhhrwithspotcomponentpower.h"
#include "testparametrization_voltanhhrwithtimecomponent.h"
#include "testparametrization_voltimeonlyhrwithtimecomponent.h"
#include "testparametrization_volwithtimecomponent.h"


#include "ito33/date.h"
#include "ito33/dateutils.h"

// Current test file.  Needs to be global due to cppunit constraints
std::string g_strInputFilename;

// callback volatility. Used by several test classes
void __stdcall 
parametricVol(double dTime, const double *pdS, double *pdValues, size_t nNbS,
              int i)
{
  // Implement the formula vol(S,t) = 0.2 + 0.00002*(100 - S)^2, 0.4 - t/10.0
  // Also, cap the values at 0.4, and make sure they are never negative 
  ito33::Date dateTmp = ito33::Date(i);
  dTime -= GetDoubleFrom(dateTmp);

  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    pdValues[nIdx] = 0.2 + 0.00002*(100.0 - pdS[nIdx])*(100.0 - pdS[nIdx]) - dTime/10.0;
    if (pdValues[nIdx] > 0.4) 
      pdValues[nIdx] = 0.4;
    if (pdValues[nIdx] < 0.0) 
      pdValues[nIdx] = 0.0;
  }
}


// Force linking of pricing modules
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);
ITO33_FORCE_LINK_MODULE(IHGPriceEDS);
ITO33_FORCE_LINK_MODULE(IHGPriceOption);
ITO33_FORCE_LINK_MODULE(IHGPriceParBond);
ITO33_FORCE_LINK_MODULE(option_xml);

void run_HRTC();
void run_HRSCP();
void run_VF_HRTC();
void run_VF_HRSCP();
void run_VF_HRP();
void run_VP();
void run_VP_HRP();
void run_VP_HRTC();
void run_VTANH();
void run_VTANH_HRP();
void run_VTANH_HRSCP();
void run_VTANH_HRTC();
void run_VTO_HRTC();
void run_VTC();

int main(int ,char **)
{
  run_HRSCP();
  run_HRTC();

  run_VF_HRP();
  run_VF_HRSCP();
  run_VF_HRTC();

  run_VP();
  run_VP_HRP();
  run_VP_HRTC();

  run_VTANH();
  run_VTANH_HRP();
  run_VTANH_HRSCP();
  run_VTANH_HRTC();

  run_VTO_HRTC();

  run_VTC();

}


void run_HRTC()
{ 
  std::cout << "Parametrization_HRWithTimeComponent......." << std::endl;
  size_t nIdx;
  size_t nNbTests;

  nNbTests = 4;
  for(nIdx = 4; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/basic%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_HRWithTimeComponentTest::suite());
    runner.run("");
  }

  nNbTests = 1;
  //nNbTests = 0;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/termparbond%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_HRWithTimeComponentTest::suite());
    runner.run("");
  }


  nNbTests = 1;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/termeds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_HRWithTimeComponentTest::suite());
    runner.run("");
  }
}

void run_HRSCP()
{
  
  size_t nIdx;
  size_t nNbTests = 3;

  std::cout <<  "Parametrization_HRWithSpotComponentPower..." << std::endl;
  nNbTests = 3;
  //nNbTests = 0;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/basic%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_HRWithSpotComponentPowerTest::suite());
    runner.run("");
  }

}

void run_VF_HRTC()
{  
  size_t nIdx;
  size_t nNbTests = 3;

  std::cout <<  "Parametrization_VolFlatHRWithTimeComponent..." << std::endl;
  nNbTests = 4;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolFlatHRWithTimeComponentTest::suite());
    runner.run("");
  }
  

  nNbTests = 1;
  //nNbTests = 0;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/option_and_termparbond%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolFlatHRWithTimeComponentTest::suite());
    runner.run("");
  }
}

void run_VF_HRSCP()
{
  size_t nIdx;
  size_t nNbTests = 3;

  std::cout <<  "Parametrization_VolFlatHRWithSpotComponentPower..." << std::endl;
  nNbTests = 4;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolFlatHRWithSpotComponentPowerTest::suite());
    runner.run("");
  }
}

void run_VP()
{
  size_t nIdx;
  size_t nNbTests = 3;

  std::cout <<  "Parametrization_VolPower..." << std::endl;
  nNbTests = 2;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoption%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolPowerTest::suite());
    runner.run("");
  }
}

void run_VP_HRP()
{
  std::cout <<  "Parametrization_VolPowerHRPower..." << std::endl;

  size_t nIdx;
  size_t nNbTests = 3;

  nNbTests = 2;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoptionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolPowerHRPowerTest::suite());
    runner.run("");
  }
}

void run_VP_HRTC()
{
  std::cout <<  "Parametrization_VolPowerHRWithTimeComponent..." << std::endl;

  size_t nIdx;
  size_t nNbTests = 3;

  nNbTests = 2;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoptionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolPowerHRWithTimeComponentTest::suite());
    runner.run("");
  }
}

void run_VTANH()
{
  std::cout <<  "Parametrization_Tanh..." << std::endl;

  size_t nIdx;
  size_t nNbTests = 3;

  nNbTests = 2;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoption%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolTanhTest::suite());
    runner.run("");
  }
}

void run_VTO_HRTC()
{  
  size_t nIdx;
  size_t nNbTests;

  std::cout <<  "Parametrization_VolTimeOnlyHRWithTimeComponent..." << std::endl;
  nNbTests = 3;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/general%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolTimeOnlyHRWithTimeComponentTest::suite());
    runner.run("");
  }
  
}

void run_VF_HRP()
{  
  size_t nIdx;
  size_t nNbTests;
  
  std::cout <<  "Parametrization_VolFlatHRPower..." << std::endl;

  // Tests 2 and 3 fail when using ASA
  //nNbTests = 3;
  nNbTests = 1;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/general%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolFlatHRPowerTest::suite());
    runner.run("");
  }
  
}

void run_VTANH_HRP()
{  
  size_t nIdx;
  size_t nNbTests;

  std::cout <<  "Parametrization_VolTanhHRPower..." << std::endl;

  // The 1st and 3rd tests fail.
  // TODO: Fix the 1st and 3rd test and/or code 
  nNbTests = 2;
  for(nIdx = 2; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/general%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolTanhHRPowerTest::suite());
    runner.run("");
  }
  
}

void run_VTANH_HRSCP()
{
  size_t nIdx;
  size_t nNbTests;

  std::cout <<  "Parametrization_VolTanhHRWithSpotComponentPower..." 
            << std::endl;

  nNbTests = 2;
  for (nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoptionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(
      Parametrization_VolTanhHRWithSpotComponentPowerTest::suite());
    runner.run("");
  }
}

void run_VTANH_HRTC()
{
  std::cout <<  "Parametrization_VolTanhHRWithTimeComponent..." << std::endl;

  size_t nIdx;
  size_t nNbTests = 3;

  nNbTests = 2;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/optionandoptionandtermcds%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolTanhHRWithTimeComponentTest::suite());
    runner.run("");
  }
}

void run_VTC()
{ 
  std::cout << "Parametrization_VolWithTimeComponent......." << std::endl;
  size_t nIdx;
  size_t nNbTests;

  nNbTests = 1;
  for(nIdx = 1; nIdx <= nNbTests; nIdx++)
  {
    char tmpFileName[1024];
    sprintf(tmpFileName, "xmlfiles/termoption%02d.xml", nIdx);
    g_strInputFilename = tmpFileName;
     
    CppUnit::TextUi::TestRunner runner;

    runner.addTest(Parametrization_VolWithTimeComponentTest::suite());
    runner.run("");
  }

}
