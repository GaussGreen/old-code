#include "ito33/beforestd.h"
#include <string>
#include <iostream>
#include <fstream>
#include <ctime>
#include "ito33/afterstd.h"

#include "ito33/array.h"
#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/pricing/options.h"

#include "ito33/finance/forwardoption.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hazardratelinear.h"

#include "ihg/forwardoptionpricer.h"

#include "forwardoptiontester.h"

using namespace std;
using namespace ito33::finance;
using namespace ito33::pricing;

using ito33::ihg::ForwardOptionTester;


void ForwardOptionTester::ReadContractParams(std::ifstream& sIn)
{

  size_t nNbOptions;
  int iType;
  std::string sMaturity;
  double dStrike;

  // make sure the exisiting list is empty
  m_listOfOptions.clear();

  // Get the number of option contracts to read
  sIn >> nNbOptions;

  // Get the type of option (put, call)
  sIn >> iType;
  finance::OptionType optionType = finance::Option_Call;
  if (iType == 0)
    optionType = finance::Option_Put;
  else if (iType == 1)
    optionType = finance::Option_Call;
  else 
    ASSERT_MSG(false, "Can only handle calls and puts for forward options");

  // Read in the strike and matrurity of each option, adding to the list
  for (size_t nIdx = 0; nIdx < nNbOptions; nIdx++)
  {
    sIn >> sMaturity >> dStrike;
    Date maturityDate(sMaturity.c_str(), "%Y-%m-%d");

    m_listOfOptions.push_back
                    (
                      shared_ptr<finance::Option>
                      ( new finance::Option
                            (dStrike, maturityDate,
                             optionType, ExerciseType_European)
                      )
                    );
  }
}

void ForwardOptionTester::ReadContractParams(std::string& sContractFile)
{
  ifstream sIn;
  sIn.open(sContractFile.c_str(), std::ios::in);
  ASSERT_MSG(sIn, "Problem opening forward option param file");
  ReadContractParams(sIn);
  sIn.close();
}


void ForwardOptionTester::RunPricer()
{  
  // Specify what needs to be computed
  m_bComputeArrays = true;
  m_bComputeSurfaces = false;

  // Create the options and forwardoptionparams objects, then set everything
  finance::ForwardOption forwardOption(m_listOfOptions);

  Options options(forwardOption);

  ForwardOptionParams forwardOptionParams(options);

  forwardOptionParams.SetValuationTime( GetDoubleFrom(m_ValuationDate) );

  forwardOptionParams.SetSpotSharePrice(m_dSpot);
  
  forwardOptionParams.SetYieldCurve(m_pYieldCurve);
  
  forwardOptionParams.SetYieldCurveForMesh(m_pYieldCurve);
  
  forwardOptionParams.SetForeignCurve(m_pForeignCurve);

  forwardOptionParams.SetDividends(m_pDividends);
 
  forwardOptionParams.SetNumParams(m_pNumParams);

  forwardOptionParams.SetMeshParams(m_pMeshParams);
   
  
  ihg::Model model(m_pVolatility, m_pVolatility, m_pHazardRate);
  

  ComputationalFlags computationalFlags;
 
  computationalFlags.SetComputeVega(false);    
  computationalFlags.SetComputeSurface(m_bComputeSurfaces);
  // The actual pricer
  ForwardOptionPricer pricer(forwardOptionParams, model, computationalFlags);
   
  m_pOutput = pricer.Price()->GetModelOutput();
    
}


void ForwardOptionTester::ReportConvergenceTest()
{
  ASSERT_MSG(m_nNbRuns > 0, 
        "Must call RunConvergenceTest before calling ReportConvergenceTest");

  cout << endl;
  cout << "Data at S=" << m_dSpot << endl << endl;

  string strText = "Price";
  ReportConvergence(m_nNbRuns, m_pdPrices.Get(), m_pdTimes.Get(), strText);
  cout << endl;
}

