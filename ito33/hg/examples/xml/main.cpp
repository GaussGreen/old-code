#include "ito33/finance/computationalflags.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/bondlike/bondlikeoutput.h"
#include "ito33/finance/bondlike/cboptionoutput.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/hedgeoutput.h"
#include "ito33/hg/hedgeratiodata.h"

#include "hg/xml/hedgereader.h"

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
    : MyOutput<finance::ModelOutput>(dynamic_pointer_cast<finance::ModelOutput>(pOutput) ),
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


template<class T>
void ComputePrice(const hg::TheoreticalModel& tm, Derivative& contract)
{  
  MyOutput<T> tmp( dynamic_pointer_cast<T>( tm.Compute(contract) ) );
  
  tmp.Print();
}

int main()
{
  try
  {    
    hg::XML::HedgeReader reader("hg_cb.xml");
    
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

      if ( bComputeHERO )
        pHedgeOutput = pModel->Hedge(*pTarget, *pHedgeInstruments);
      else
        pHedgeOutput = pModel->ComputeHERO(*pTarget, *pHedgeInstruments);
        
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

    shared_ptr<Derivative> pDerivative;
        reader.ReadDerivative(pDerivative);
	
    shared_ptr<finance::ComputationalFlags> 
      pFlags(new finance::ComputationalFlags);
    // pFlags->SetComputeSurface(true);
    
    pDerivative->SetComputationalFlags(pFlags);
    
    std::cout.precision(16);
    
    ComputePrice<finance::ModelOutput>(*pModel, *pDerivative);
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

