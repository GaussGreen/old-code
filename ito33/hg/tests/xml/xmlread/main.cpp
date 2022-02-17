#include "ito33/timer.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/hg/bondlikeoutput.h"
#include "ito33/hg/cboptionoutput.h"
#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/hedgeratiodata.h"

#include "hg/numoutput.h"

#include "ito33/tests/showconvergence.h"
#include "hg/xml/hedgereader.h"

#include <iostream>

using namespace ito33;
using namespace ito33::finance;

class MyOutput
{
public:

  MyOutput(shared_ptr<finance::ModelOutput> pOutput) : m_pOutput(pOutput)
  {
  }

  void Print()
  {
    std::cout.precision(12);
    
    std::cout << "Price " << m_pOutput->GetPrice() << std::endl;
    
    std::cout << "Delta " << m_pOutput->GetDelta() << std::endl;
    std::cout << "Gamma " << m_pOutput->GetGamma() << std::endl;

    std::cout << "RecoveryValue " << m_pOutput->GetValueAfterDefault()
              << std::endl;

    if ( m_pOutput->HasTheta() )
      std::cout << "Theta " << m_pOutput->GetTheta() << std::endl;
    
    if ( m_pOutput->HasVega() )
      std::cout << "Vega " << m_pOutput->GetVega() << std::endl;

    if ( m_pOutput->HasRho() )
      std::cout << "Rho " << m_pOutput->GetRho() << std::endl;

    if ( m_pOutput->HasFugit() )
      std::cout << "Fugit " << m_pOutput->GetFugit() << std::endl;
  }

private:

  shared_ptr<finance::ModelOutput> m_pOutput;
};

int main()
{
  try
  {    
    //hg::XML::HedgeReader reader("c:\\ito33\\output\\hgvarswap.xml");
    hg::XML::HedgeReader reader("c:\\ito33\\output\\hg.xml");
    
    shared_ptr<hg::TheoreticalModel> pModel;
    
    reader.ReadTheoreticalModel(pModel);
    
    //pModel->SetDebugOutputFile("hg_cb.xml"); 
    
    shared_ptr<Derivative> pTarget( reader.ReadTarget() );
    if ( pTarget )
    {
      shared_ptr<Derivatives> 
        pHedgeInstruments( reader.ReadHedgeInstruments() );

      bool bComputeHERO = reader.ReadComputeHERO();

      shared_ptr<hg::HedgeOutput> pHedgeOutput;
      
      if (bComputeHERO)
        pHedgeOutput = pModel->ComputeHERO(*pTarget, *pHedgeInstruments);
      else
        pHedgeOutput = pModel->Hedge(*pTarget, *pHedgeInstruments);

      // Output HERO
      if ( pHedgeOutput->HasHERO() )
        std::cout << "HERO = " 
                  << pHedgeOutput->GetHERO()
                  << std::endl 
                  << std::endl;

      // Output hedge ratios
      double dUnderlyingRatio = pHedgeOutput->GetUnderlyingHedgeRatio();
      std::cout << "Underlying ratio = " << dUnderlyingRatio << std::endl;

      std::vector< shared_ptr< hg::HedgeRatioData> > 
        ppHedgeRatioData( pHedgeOutput->GetHedgeRatioData() );
      for (size_t nIdx = 0; nIdx < ppHedgeRatioData.size(); nIdx++)
      {
        double dRatio = ppHedgeRatioData[nIdx]->GetRatio();
        std::cout << " Ratio " << nIdx << " = " << dRatio 
                  << std::endl;
      }
      std::cout << std::endl;

      return 0;
    }

    // Try to read the contract
    shared_ptr<finance::Derivative> pDerivative;
    reader.ReadDerivative(pDerivative);    

    shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);
    pFlags->SetComputeSurface(true);
    pFlags->SetComputeVega(true);
    pFlags->SetComputeRho(true);
    pFlags->SetComputeFugit(true);
    pFlags->ActivateAllSensitivities(true);

    pDerivative->SetComputationalFlags(pFlags);

    shared_ptr<finance::ModelOutput> pMO(pModel->Compute(*pDerivative));
    
    shared_ptr<hg::NumOutput>
      pNumOutput( dynamic_pointer_cast<hg::NumOutput>( pMO->GetNumOutput() ) );
    
    // Output sensitivities, if computed
    if ( pNumOutput->HasSensitivities() )
    {
      std::cout << std::endl;
      std::cout << "Sensitivities:" << std::endl;
      std::vector<double> pdSensitivities = pNumOutput->GetSensitivities();
      for (size_t nIdx = 0; nIdx < pdSensitivities.size(); nIdx++)
        std::cout << "index = " << nIdx << ", value = " << pdSensitivities[nIdx] << std::endl;
    }
    std::cout << std::endl;

    MyOutput myOutput(pMO);

    myOutput.Print();

    StopWatch sw;
    ShowConvergence(*pModel, *pDerivative, 4);
    std::cout << "Total time = " << sw() << std::endl;
  }
  catch ( ito33::Exception& e )
  {
    printf("ITO33 exception:\n%s\n", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    printf("std exception: %s\n", e.what());
  }
  catch ( ... )
  {
    puts("unknown exception!");
  }

  return 0;
}
