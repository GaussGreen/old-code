#include "ito33/beforestd.h"
#include <vector>
#include <stdio.h>
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/modeloutput.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/tests/showconvergence.h"
#include "ihg/xml/pricingreader.h"

#include "ihg/tests/analysisdate_surface.h"
#include "ihg/tests/test_greeks.h"


using namespace ito33;
using namespace ito33::finance;

void PrintOutput(finance::ModelOutput& output)
{
  std::cout << "Price " << output.GetPrice() << std::endl;
  std::cout << "Delta " << output.GetDelta() << std::endl;
  std::cout << "Gamma " << output.GetGamma() << std::endl;
  std::cout << "Theta " << output.GetTheta() << std::endl;

  if ( output.HasRho() )
    std::cout << "Rho   " << output.GetRho() << std::endl;

  if ( output.HasVega() )
    std::cout << "Vega  " << output.GetVega() << std::endl;

  if ( output.HasFugit() )
    std::cout << "Fugit " << output.GetFugit() << std::endl;

  std::cout << std::endl;
}



template<class T>
void ComputePrice(const ihg::TheoreticalModel& tm,
             const shared_ptr<SessionData>& pSessionData,
             T& contract
            )
{
  contract.SetSessionData(pSessionData);
  PrintOutput(*tm.Compute(contract));

  ihg::AnalysisDateTest myTests(tm, &contract);

  myTests.Report();
  
  ihg::SurfaceTest surfaceTest(tm, &contract);

  surfaceTest.Report();

  // ShowConvergence(tm, contract, 4);
}

int main()
{
  try
  {
    /*************************************************************************
    * Enter list of xml file you need to test                               
    * file+absolute path                                                    
    * so if file is in testsuite/bondlike/xmlfiles                          
    * filename = "../bondlike/xmlfile/<filename>                            
    /*************************************************************************/
    std::list<std::string> fileList;
    
    fileList.push_back("../bondlike/xmlfiles/basic01.xml");
    fileList.push_back("../bondlike/xmlfiles/basic02.xml");
    fileList.push_back("../bondlike/xmlfiles/basic03.xml");
    fileList.push_back("../bondlike/xmlfiles/basic04.xml");
    fileList.push_back("../bondlike/xmlfiles/basic05.xml");
    fileList.push_back("../bondlike/xmlfiles/basic06.xml");
    fileList.push_back("../bondlike/xmlfiles/basic07.xml");
    fileList.push_back("../bondlike/xmlfiles/basic08.xml");
    fileList.push_back("../bondlike/xmlfiles/basic09.xml");
    fileList.push_back("../bondlike/xmlfiles/basic10.xml");
    fileList.push_back("../bondlike/xmlfiles/basic11.xml");
    fileList.push_back("../bondlike/xmlfiles/basic12.xml");
    fileList.push_back("../bondlike/xmlfiles/basic13.xml");
    fileList.push_back("../bondlike/xmlfiles/bond01.xml");
    fileList.push_back("../bondlike/xmlfiles/bond02.xml");
    fileList.push_back("../bondlike/xmlfiles/cbnewshare01.xml");
    fileList.push_back("../bondlike/xmlfiles/cbnewshare02.xml"); 
    fileList.push_back("../bondlike/xmlfiles/cbnewshare03.xml"); 
    fileList.push_back("../bondlike/xmlfiles/cbnewshare04.xml"); 
    fileList.push_back("../bondlike/xmlfiles/cbnewshare05.xml");
    fileList.push_back("../bondlike/xmlfiles/cbnewshare06.xml");
    fileList.push_back("../bondlike/xmlfiles/cbnewshare07.xml");
    fileList.push_back("../bondlike/xmlfiles/cboption01.xml");  
    fileList.push_back("../bondlike/xmlfiles/cboption02.xml");  
    fileList.push_back("../bondlike/xmlfiles/cboption03.xml");    
    fileList.push_back("../bondlike/xmlfiles/cboption04.xml");    
    fileList.push_back("../bondlike/xmlfiles/cboption05.xml");    
    fileList.push_back("../bondlike/xmlfiles/cboption06.xml");   
    fileList.push_back("../bondlike/xmlfiles/cboption07.xml");              
    fileList.push_back("../bondlike/xmlfiles/conversion01.xml");  
    fileList.push_back("../bondlike/xmlfiles/conversion02.xml");
    fileList.push_back("../bondlike/xmlfiles/conversion03.xml");
    fileList.push_back("../bondlike/xmlfiles/generalizedpeps01.xml");
    fileList.push_back("../bondlike/xmlfiles/generalizedpeps02.xml"); 
    fileList.push_back("../bondlike/xmlfiles/generalizedpeps03.xml");  
    fileList.push_back("../bondlike/xmlfiles/peps01.xml");
    fileList.push_back("../bondlike/xmlfiles/peps02.xml");
    fileList.push_back("../bondlike/xmlfiles/peps03.xml");
    fileList.push_back("../bondlike/xmlfiles/percs01.xml");
    fileList.push_back("../bondlike/xmlfiles/percs02.xml");
    fileList.push_back("../bondlike/xmlfiles/percs03.xml");
    fileList.push_back("../bondlike/xmlfiles/reset01.xml");
    fileList.push_back("../bondlike/xmlfiles/reset02.xml");
    fileList.push_back("../bondlike/xmlfiles/reset03.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D01.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D02.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D03.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D04.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D05.xml");
    fileList.push_back("../bondlike/xmlfiles/reset1D06.xml");
    fileList.push_back("../bondlike/xmlfiles/reset100.xml");
    fileList.push_back("../bondlike/xmlfiles/reset101.xml");
    fileList.push_back("../bondlike/xmlfiles/reset102.xml");
    fileList.push_back("../bondlike/xmlfiles/warrants01.xml");
    fileList.push_back("../bondlike/xmlfiles/warrants02.xml");
    
    

    /**************************************************************************
      First implementation to test the greeks
      Please feel free to add more things or complete if necessary
    ***************************************************************************/
    ihg::TestGreeks testGreeks;
    testGreeks.Run(fileList);

    return 0;

    std::cout.precision(16);
    ihg::XML::PricingReader reader("./ihg.xml");

    // ihg::XML::Reader reader("c:\\ito33\\output\\ihgoption.xml");
    shared_ptr<finance::SessionData> pSessionData(reader.ReadSessionData());

    DerivativeVisitorGoodType visitor;
    reader.ReadDerivatives(visitor);
    
    shared_ptr<ito33::ihg::TheoreticalModel> 
      pModel(new ito33::ihg::TheoreticalModel);
    reader.ReadTheoreticalModel(pModel);
    pModel->SetDebugOutputFile("ihg.xml");
    
    printf("\nSessionData:\n"
           "--------------------------\n"
           "Spot:\t\t%g\n"
           "Val. date:\t%s\n\n",
           pSessionData->GetSpotSharePrice(),
           pSessionData->GetValuationDate().Format().c_str());
    
    finance::ModelOutput output;
    reader.ReadOutput(output);

    printf("Output:\n"
           "--------------------------\n"
           "Price:\t\t%g\n"
           "Delta:\t\t%g\n"
           "\n",
           output.GetPrice(),
           output.GetDelta());
    
    if(visitor.GetOption())
    {
      Option& option = *visitor.GetOption();

      printf("Option:\n"
            "--------------------------\n"
            "Type:\t\t%s\n"
            "Exercise:\t%s\n"
            "Strike:\t\t%g\n"
            "Exp. date:\t%s\n\n",
            option.GetOptionType() == finance::Option_Put ? "put" : "call",
            option.GetExerciseType() == finance::ExerciseType_European
                                          ? "european"
                                          : "american",
            option.GetStrike(),
            option.GetMaturityDate().Format().c_str()
            );

      ComputePrice(*pModel, pSessionData, option);

      /*
      ihg::AnalysisDateTest myTests(*pModel, &option);

      myTests.Report();
      
      ihg::SurfaceTest surfaceTest(*pModel, &option);

      surfaceTest.Report();

      
      if ( ! )
        exit(-1);
      */
    }

    if(visitor.GetCDS())
    {
      CDS& cds = *visitor.GetCDS();

      printf("CDS:\n"
            "--------------------------\n"
            "Issue Date:\t%s\n"
            "MaturityDate Date:\t%s\n"
            "Recovery rate:\t%g\n\n",
            cds.GetIssueDate().Format().c_str(),
            cds.GetMaturityDate().Format().c_str(),
            cds.GetRecoveryRate()
            );
      ComputePrice(*pModel, pSessionData, cds);
    }
    
    
    if(visitor.GetBond())
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetBond());
    }
    
    if(visitor.GetConvertibleBond())
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetConvertibleBond());
    }
    

    if(visitor.GetPEPSLike())
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetPEPSLike());
    }
    
    if(visitor.GetPERCSLike())
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetPERCSLike());
    }


    if ( visitor.GetReset() )
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetReset() );
    }


    if ( visitor.GetAttachedWarrantCB() )
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetAttachedWarrantCB() );
    }


    if ( visitor.GetCBOption() )
    {
      ComputePrice(*pModel, pSessionData, *visitor.GetCBOption() );
    }
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

