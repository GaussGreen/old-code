#include <iostream>

#include "ito33/exception.h"

#include "ito33/finance/derivatives.h"

#include "ito33/hg/parametrization.h"
#include "ito33/hg/theoreticalmodel.h"

#include "hg/xml/calibrationreader.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(option_xml);
ITO33_FORCE_LINK_MODULE(cds_xml);
ITO33_FORCE_LINK_MODULE(eds_xml);
ITO33_FORCE_LINK_MODULE(varianceswap_xml);
ITO33_FORCE_LINK_MODULE(yieldcurve_xml);

ITO33_FORCE_LINK_MODULE(HGPriceOption);
ITO33_FORCE_LINK_MODULE(HGPriceCDS);
ITO33_FORCE_LINK_MODULE(HGPriceEDS);
ITO33_FORCE_LINK_MODULE(HGPriceVarianceSwap);

using namespace ito33;
using namespace ito33::finance;

int main()
{
  try
  {
    std::cout.precision(16);

    hg::XML::CalibrationReader 
      reader("c:\\ito33\\output\\hg_calibration.xml");
    
    // Get the parametrization
    shared_ptr<hg::Parametrization>
      pParametrization( reader.ReadParametrization() );

    shared_ptr<Derivatives> pDerivatives( reader.ReadDerivatives() );

    shared_ptr<hg::UnderlyingProcess> pProcess;
    
    try
    {
      pProcess = pParametrization->Calibrate(*pDerivatives);
    }
    catch( const ito33::Exception& e)
    {
      pProcess = pParametrization->GetCalibratedUnderlyingProcess();

      std::cout << "calibration failed:" << e.GetErrorMessage() << std::endl;
    }

    hg::TheoreticalModel model(pProcess);

    Derivatives::Elements elements = pDerivatives->GetAll();
    Derivatives::Elements::const_iterator iter;
    for (iter = elements.begin(); iter != elements.end(); ++iter)
    {
      shared_ptr<finance::ModelOutput> pOutput = model.Compute(*iter->first);

      double dComputedPrice = pOutput->GetPrice();

      double dMarketPrice = iter->first->GetMarketPrice();

      std::cout << "Market price = " << dMarketPrice
                << ", computed price = " << dComputedPrice
                << std::endl;

    }
  }
  catch ( const ito33::Exception& e )
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
