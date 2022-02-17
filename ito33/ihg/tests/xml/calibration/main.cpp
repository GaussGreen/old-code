#include <iostream>

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/derivativevisitors/derivative_visitor_goodtype.h"
#include "ito33/finance/termstructure_enumerator.h"
#include "ito33/finance/derivatives.h"

#include "ito33/ihg/modeloutput.h"
#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/tests/showconvergence.h"
#include "ihg/xml/calibrationreader.h"
#include "ihg/xml/parametrization_visitor_goodtype.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(option_xml);
ITO33_FORCE_LINK_MODULE(cds_xml);
ITO33_FORCE_LINK_MODULE(yieldcurve_xml);

ITO33_FORCE_LINK_MODULE(IHGPriceOption);
ITO33_FORCE_LINK_MODULE(IHGPriceCDS);

using namespace ito33;
using namespace ito33::finance;

int main()
{
  try
  {
    std::cout.precision(16);

    ihg::XML::CalibrationReader 
      reader("c:\\ito33\\output\\ihg_calibration.xml");

    // Get the sessiondata
    shared_ptr<SessionData> pSessionData = reader.ReadSessionData();
    
    // Get the parametrization
    ihg::ParametrizationVisitorGoodType param_visitor;
    DerivativeVisitorGoodType deriv_visitor;
    TermStructureEnumerator termstructures;
    shared_ptr<ito33::finance::Derivatives> 
      pDerivatives(new ito33::finance::Derivatives());

    reader.ReadCalibration(param_visitor, termstructures, deriv_visitor, 
                           *pDerivatives);

    if (param_visitor.GetParametrizationHRWithTimeComponent())
    {
      shared_ptr<ihg::ParametrizationHRWithTimeComponent>
        pParametrzation = param_visitor.GetParametrizationHRWithTimeComponent();
    
      // Calibrate a cds term structure
      if ( termstructures.GetTermStructureCDS() )
      {
        shared_ptr<TermStructureCDS> 
          ptsCDS = termstructures.GetTermStructureCDS();
        
        pParametrzation->CalibrateWithCDSs(*ptsCDS);
      }
    }

    if (param_visitor.GetParametrizationHRWithSpotComponentPower())
    {
      shared_ptr<ihg::ParametrizationHRWithSpotComponentPower>
        pParametrzation = param_visitor.GetParametrizationHRWithSpotComponentPower();

      // Expect a cds termstructure 
      shared_ptr<TermStructureCDS> 
        ptsCDS = termstructures.GetTermStructureCDS();
  
      pParametrzation->CalibrateWithCDSs(*ptsCDS);
    }  

    if (param_visitor.GetParametrizationVolFlatHRWithTimeComponent())
    {
      shared_ptr<ihg::ParametrizationVolFlatHRWithTimeComponent>
        pParametrzation = param_visitor.GetParametrizationVolFlatHRWithTimeComponent();
    
      // expect an option
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      // Calibrate a cds term structure
      if ( termstructures.GetTermStructureCDS() )
      {
        shared_ptr<TermStructureCDS> 
          ptsCDS = termstructures.GetTermStructureCDS();

        pParametrzation->CalibrateWithOptionAndCDSs(*pOption, *ptsCDS);
      }

      if ( termstructures.GetTermStructureOption() )
      {
        shared_ptr<TermStructureOption> 
          ptsOption = termstructures.GetTermStructureOption();

        pParametrzation->CalibrateWithOptionAndOptions(*pOption, *ptsOption);
      }
    }

    if (param_visitor.GetParametrizationVolFlatHRWithSpotComponentPower())
    {
      shared_ptr<ihg::ParametrizationVolFlatHRWithSpotComponentPower>
        pParametrzation = param_visitor.GetParametrizationVolFlatHRWithSpotComponentPower();

      // expect an option 
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      // Expect a cds termstructure 
      shared_ptr<TermStructureCDS> 
        ptsCDS = termstructures.GetTermStructureCDS();
  
      pParametrzation->CalibrateWithOptionAndCDSs(*pOption, *ptsCDS);
    } 

    if (param_visitor.GetParametrizationVolPowerHRWithTimeComponent())
    {
      shared_ptr<ihg::ParametrizationVolPowerHRWithTimeComponent>
        pParametrzation = param_visitor.GetParametrizationVolPowerHRWithTimeComponent();
    
      // expect two options
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      shared_ptr<Option> pOption2 = deriv_visitor.GetOption2();
  
      // Calibrate a cds term structure
      if ( termstructures.GetTermStructureCDS() )
      {
        shared_ptr<TermStructureCDS> 
          ptsCDS = termstructures.GetTermStructureCDS();

        pParametrzation->CalibrateWithOptionsAndCDSs(*pOption, *pOption2, *ptsCDS);
      }
    }

    if (param_visitor.GetParametrizationVolPowerHRPower())
    {
      shared_ptr<ihg::ParametrizationVolPowerHRPower>
        pParametrzation = param_visitor.GetParametrizationVolPowerHRPower();

      // expect two options
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      shared_ptr<Option> pOption2 = deriv_visitor.GetOption2();

      // Expect a cds termstructure 
      shared_ptr<TermStructureCDS> 
        ptsCDS = termstructures.GetTermStructureCDS();
  
      pParametrzation->CalibrateWithOptionsAndCDSs(*pOption, *pOption2, *ptsCDS);
    } 

    if (param_visitor.GetParametrizationVolTanh())
    {
      shared_ptr<ihg::ParametrizationVolTanh>
        pParametrzation = param_visitor.GetParametrizationVolTanh();
    
      // expect two options
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      shared_ptr<Option> pOption2 = deriv_visitor.GetOption2();

      pParametrzation->CalibrateWithOptions(*pOption, *pOption2);
    }

    if (param_visitor.GetParametrizationVolTanhHRWithSpotComponentPower())
    {
      shared_ptr<ihg::ParametrizationVolTanhHRWithSpotComponentPower>
        pParametrzation = param_visitor.GetParametrizationVolTanhHRWithSpotComponentPower();

      // expect two options
      shared_ptr<Option> pOption = deriv_visitor.GetOption();

      shared_ptr<Option> pOption2 = deriv_visitor.GetOption2();

      // Expect a cds termstructure 
      shared_ptr<TermStructureCDS> 
        ptsCDS = termstructures.GetTermStructureCDS();

      pParametrzation->CalibrateWithOptionsAndCDSs(*pOption, *pOption2, *ptsCDS);
    } 

    if ( param_visitor.GetParametrizationVolWithTimeComponent() )
    {
      shared_ptr<ihg::ParametrizationVolWithTimeComponent>
        pParametrzation = param_visitor.GetParametrizationVolWithTimeComponent();

      // Expect an option terms structure 
      shared_ptr<TermStructureOption> 
        ptsOption = termstructures.GetTermStructureOption();
  
      pParametrzation->CalibrateWithOptions(*ptsOption);
    } 

    if ( param_visitor.GetParametrizationVolTimeOnlyHRWithTimeComponent() )
    {
      shared_ptr<ihg::ParametrizationVolTimeOnlyHRWithTimeComponent>
        pParametrzation = param_visitor.GetParametrizationVolTimeOnlyHRWithTimeComponent();

      // Expect a general derivative list
        
      pParametrzation->CalibrateWithDerivatives(*pDerivatives);
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

