#include "ito33/timer.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/modeloutput.h"
#include "ito33/ihg/bondlikeoutput.h"
#include "ito33/ihg/cboptionoutput.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"
#include "ihg/xml/pricingreader.h"

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
    std::cout << "Price " << m_pOutput->GetPrice() << std::endl;
    
    std::cout << "Delta " << m_pOutput->GetDelta() << std::endl;
    std::cout << "Gamma " << m_pOutput->GetGamma() << std::endl;
    std::cout << "Theta " << m_pOutput->GetTheta() << std::endl;
    std::cout << "RecoveryValue " << m_pOutput->GetValueAfterDefault()
              << std::endl;

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
    std::cout.precision(16);

    ihg::XML::PricingReader reader("c:\\ito33\\output\\ihg_option.xml");
    
    shared_ptr<ito33::ihg::TheoreticalModel> 
      pModel(new ito33::ihg::TheoreticalModel);
    reader.ReadTheoreticalModel(pModel);
    
    // pModel->SetDebugOutputFile("ihg.xml");
    shared_ptr<Derivative> pDerivative;
    reader.ReadDerivative(pDerivative);

    shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);
    pFlags->SetComputeSurface(true);
    pFlags->SetComputeVega(true);
    pFlags->SetComputeRho(true);
    pFlags->SetComputeFugit(true);

    pDerivative->SetComputationalFlags(pFlags);

    shared_ptr<finance::ModelOutput> pMO( pModel->Compute(*pDerivative) );

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

