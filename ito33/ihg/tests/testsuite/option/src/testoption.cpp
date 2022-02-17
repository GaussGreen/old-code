

#include "ito33/exception.h"

using namespace std; 

ito33::ihg::OptionTester optionTester;

void RunAutomatedTests();

class TestClass
{
public:
  TestClass() : m_dValue(0.0) {}

  double GetValue() { return m_dValue; }
  void SetValue(double dValue) { m_dValue = dValue; }

protected:
  double m_dValue;
};

int main()
{
  try {

/*
    ito33::shared_ptr<TestClass> pTest(new TestClass());
    std::cout << pTest->GetValue() << std::endl;

    pTest->SetValue(5.0);
    double dTmp = pTest->GetValue();

    std::cout << pTest->GetValue() << std::endl;

    ito33::shared_ptr<TestClass> pTest2 = pTest;

    std::cout << pTest2->GetValue() << std::endl;


    exit(0);
*/


    //RunAutomatedTests();
    //exit(0);

    string sYieldCurveFile = g_strPath + "common\\constyieldmedium.txt";
    string sForeignCurveFile = g_strPath + "common\\constyieldsmall.txt";
    string sDividendFile = g_strPath + "common\\mixeddividends.txt";
    string sVolatilityFile = g_strPath + "common\\constvolmedium.txt";
    string sHazardRateFile = g_strPath + "common\\consthrzero.txt";
    string sMeshParamFile = g_strPath + "common\\uniformmeshes.txt";
    string sNumParamFile = g_strPath + "common\\implicitnumparams.txt";
    //string sNumParamFile = g_strPath + "common\\cnnumparams.txt";
    //string sNumParamFile = g_strPath + "common\\bdfnumparams.txt";
    string sValuationDateFile = g_strPath + "common\\valuationdate.txt";
    //string sOptionFile = g_strPath + "option\\amercall.txt";
    string sOptionFile = g_strPath + "option\\eurocall.txt";
    
    optionTester.ReadFromMultipleFiles(sYieldCurveFile, sForeignCurveFile,
      sDividendFile, sVolatilityFile, sHazardRateFile, sMeshParamFile,
      sNumParamFile, sValuationDateFile, sOptionFile);

    // To read from a single file, use this code
    //string strFilename = g_strPath + "option\\testAmerPut1.txt";
    //optionTester.ReadInputFile(strFilename);

    //optionTester.SetVolatility(ito33::shared_ptr<ito33::ihg::Volatility>
    //  (new ito33::ihg::VolatilityCallBack(parametricVol,0)));

    //optionTester.SetHazardRate(ito33::shared_ptr<ito33::ihg::HazardRate>
    //  (new ito33::ihg::HazardRateCallBack(parametricHR,0)));

    //int iDate = 0;
    //ito33::Date pricingDate(iDate);
    //optionTester.SetValuationDate(pricingDate);

    optionTester.SetSpotSharePrice(100.0);
  
    size_t nNbTests = 1;
    optionTester.RunConvergenceTest(nNbTests);    

    if (nNbTests >= 3)
      optionTester.ReportConvergenceTest();
    else
    {
      cout << "Price = " << optionTester.GetPrice() << endl;
      cout << "Delta = " << optionTester.GetDelta() << endl;
      cout << "Gamma = " << optionTester.GetGamma() << endl;
      cout << "Vega = " << optionTester.GetVega() << endl;

    }

  }
  catch(const ito33::Exception &e)
  {
    cout << e.GetErrorMessage() << endl;
  }

  return 0;
}




