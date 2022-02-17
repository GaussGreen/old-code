#include "ito33/finance/sessiondata.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"
#include "ito33/finance/bondlike/cboptionoutput.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ihg/xml/pricingreader.h"

#include <iostream>

using namespace ito33;
using namespace ito33::finance;

template<class T>
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
    
    if ( m_pOutput->HasFugit() )
      std::cout << "Fugit " << m_pOutput->GetFugit() << std::endl;

    if ( m_pOutput->HasVega() )      
      std::cout << "Vega " << m_pOutput->GetVega() << std::endl;
    
    if ( m_pOutput->HasRho() )      
      std::cout << "Rho " << m_pOutput->GetRho() << std::endl;
    
    std::cout << "RecoveryValue " << m_pOutput->GetValueAfterDefault() 
              << std::endl;
  }

private:
  shared_ptr<finance::ModelOutput> m_pOutput;
};

template<>
class MyOutput<finance::BondLikeOutput> : public MyOutput<finance::ModelOutput>
{
public:
  MyOutput(shared_ptr<finance::BondLikeOutput> pOutput)
    : MyOutput<finance::ModelOutput>( dynamic_pointer_cast<finance::ModelOutput>(pOutput) ),
      m_pOutput(pOutput)
  {
  }

  void Print()
  {
    MyOutput<finance::ModelOutput>::Print();

    std::cout << "BdFlr " << m_pOutput->GetBondFloor() << std::endl;
  }

private:
  shared_ptr<finance::BondLikeOutput> m_pOutput;
};


template<>
class MyOutput<finance::CBOptionOutput> : public MyOutput<finance::ModelOutput>
{
public:
  MyOutput(shared_ptr<finance::CBOptionOutput> pOutput)
    : MyOutput<finance::ModelOutput>( dynamic_pointer_cast<finance::ModelOutput>(pOutput) ),
      m_pOutput(pOutput)
  {
  }

  void Print()
  {
    std::cout << "CB option results : " << std::endl;
    MyOutput<finance::ModelOutput>::Print();
    
    MyOutput< shared_ptr<finance::BondLikeOutput> > 
      mycboutput( m_pOutput->GetCBOutput());
    
    std::cout << std::endl << "CB results :" << std::endl;

    mycboutput.Print();
  }

private:
  shared_ptr<finance::CBOptionOutput> m_pOutput;
};


template<class T>
void ComputePrice(const ihg::TheoreticalModel& tm, Derivative& contract)
{
  MyOutput<T> tmp( dynamic_pointer_cast<T>( tm.Compute(contract) ) );
  tmp.Print();
}

int main()
{
  try
  {
    std::cout.precision(16);

    ihg::XML::PricingReader reader("ihg_cbcrash.xml");

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
    
    ComputePrice<finance::ModelOutput>(*pModel, *pDerivative );
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

