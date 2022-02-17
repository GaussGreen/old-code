#include "ito33/beforestd.h"
#include <iostream>
#include <time.h>
#include <stdlib.h>
#include "ito33/afterstd.h"

#include "ito33/date.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/modeloutput.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/computationalflags.h"
#include "ito33/finance/forwardoption.h"
#include "ito33/finance/dividends.h"

#include "ito33/hg/theoreticalmodel.h"
#include "ito33/hg/multioutput.h"

#include "hg/priceforwardoption.h"
#include "hg/translator.h"
#include "hg/numoutput.h"

#include "ito33/tests/showconvergence.h"

#include "ito33/link.h"

ITO33_FORCE_LINK_MODULE(HGPriceOption);

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;


std::list< shared_ptr<Option> > 
ConstructOptionList( const shared_ptr<SessionData> &pSessionData )
{
  //seed the random number generator with time
  srand( (unsigned)time( NULL ) );

  // the return derivatives
  std::list< shared_ptr<Option> > optionList;

  // Get session data
  double dSpot = pSessionData->GetSpotSharePrice();
  Date valuationDate = pSessionData->GetValuationDate();

  // Construct some options
  double dStrikeMin = dSpot / 2.0;
  double dStrikeDelta = 2.0 * (dSpot / 10.0);
  double dStrikeMax = dSpot * 1.5;
  //OptionType optionType = Option_Put;
  OptionType optionType = Option_Call;

  size_t nNbYears = 3;
  Date maturityDate = valuationDate;
  for (size_t nIdx = 0; nIdx < nNbYears; nIdx++)
  {
    maturityDate.AddYears(1);
   
    double dStrike = dStrikeMin;
    while (dStrike <= dStrikeMax + 1.e-8)
    {

      optionType = (rand()/double(RAND_MAX) >= .5)? Option_Call:Option_Put;

      shared_ptr<Option> 
        pOption( new Option(dStrike, maturityDate, 
                            optionType, ExerciseType_European) );
      
      pOption->SetSessionData(pSessionData);

      pOption->SetMarketPrice(1);

      optionList.push_back(pOption);

      dStrike += dStrikeDelta;

    } // loop over strike

  } // loop over maturities

  return optionList;
}


int main()
{
  try
  {
    // Make the session
    Date valuationDate(2003, Date::Feb, 1);
    shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );

    shared_ptr<Equity> pEquity(new Equity(45.0, pCurrency));

    shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.01));

    pEquity->SetBorrowCurve(pyf);

    shared_ptr<Dividends> pDividends( new Dividends );
    pDividends->AddCash(Date(2003, Date::Jul, 1), 2.0);
    pEquity->SetDividends(pDividends);

    shared_ptr<YieldCurve> pyc( new YieldCurveFlat(0.05) );

    shared_ptr<RateData> pRateData(new RateData);
    pRateData->SetYieldCurve(pCurrency, pyc);
    
    shared_ptr<SessionData>
      pSessionData( new SessionData(pRateData, pEquity, valuationDate) );

    // Make the model
    size_t nNbRegimes = 2;

    std::vector<double> pdVols;
    pdVols.resize(nNbRegimes, 0.5);

    std::vector<double> pdDefaultIntensities;
    pdDefaultIntensities.resize(nNbRegimes, 0.1);

    shared_ptr<hg::UnderlyingProcess>
      pUnderlyingProcess( new hg::UnderlyingProcess
                              (nNbRegimes, pdVols, pdDefaultIntensities) );

    hg::Jumps jumps;

    jumps.push_back(hg::Jump(0.2, -0.3));
    pUnderlyingProcess->SetJumps(0, 0, jumps);

    if (nNbRegimes > 1)
    {
      jumps.clear();
      jumps.push_back(hg::Jump(0.5, -0.4));
      pUnderlyingProcess->SetJumps(0, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.15, -0.1));
      pUnderlyingProcess->SetJumps(1, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.3, -0.2));
      pUnderlyingProcess->SetJumps(1, 0, jumps); 
    }
    /**/

    if (nNbRegimes > 2)
    {
      jumps.clear();
      jumps.push_back(hg::Jump(0.1, -0.2));
      pUnderlyingProcess->SetJumps(0, 2, jumps); 
 
      jumps.clear();
      jumps.push_back(hg::Jump(0.1, -0.3));
      pUnderlyingProcess->SetJumps(1, 2, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.3, -0.2));
      pUnderlyingProcess->SetJumps(2, 0, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.2, -0.4));
      pUnderlyingProcess->SetJumps(2, 1, jumps); 

      jumps.clear();
      jumps.push_back(hg::Jump(0.05, -0.25));
      pUnderlyingProcess->SetJumps(2, 2, jumps); 
    }

    shared_ptr<hg::TheoreticalModel> 
      pModel(new hg::TheoreticalModel(pUnderlyingProcess));
    
    shared_ptr<ComputationalFlags> pFlags(new ComputationalFlags);

    pFlags->ActivateAllSensitivities(true);
    pFlags->SetSensitivityMethod(2);
    //pFlags->SetSensitivityMethod(1);

    pModel->SetExternalFlags(pFlags);

    // make the option list
    std::list< shared_ptr<Option> > 
      optionList = ConstructOptionList(pSessionData);

    // Make the forward option and then price
    shared_ptr<ForwardOption> forwardOpt(new ForwardOption(optionList));
 
    shared_ptr<hg::MultiOutput> 
      forwardOutput( hg::PriceForwardOption(*pModel, *forwardOpt) );
    
    double dObjectif = forwardOutput->GetObjectif();

    std::vector<double> pdPrices;
    pdPrices = forwardOutput->GetPrices();

    std::cout << "Forward prices: " << std::endl;
    size_t nNbOptions = optionList.size();

    for (size_t nIdx = 0; nIdx < nNbOptions; nIdx++)
    {
      std::cout << pdPrices[nIdx] << std::endl;
    }
    std::cout << std::endl;

    shared_ptr<hg::NumOutput>
      pNumOutput( dynamic_pointer_cast<hg::NumOutput>
                  ( forwardOutput->GetNumOutput() ) );

    // Output sensitivities, if computed
    if ( pNumOutput->HasSensitivities() )
    {
      std::cout << std::endl;
      std::cout << "Sensitivities:" << std::endl;

      if ( forwardOutput->IsSensitivityOnObjectif() )
      {
        std::vector<double> 
          pdSensitivities = pNumOutput->GetSensitivities();

        for (size_t nIdxD = 0; nIdxD < pdSensitivities.size(); nIdxD++)
          std::cout << pdSensitivities[nIdxD] << std::endl;
      }
      else
      {
        std::vector< std::vector<double> > 
          ppdSensitivities = forwardOutput->GetMultiSensitivities();
      
        size_t nNbParams = ppdSensitivities[0].size();
        for (size_t nIdxOption = 0; nIdxOption < nNbOptions; nIdxOption++)
        {
          for (size_t nIdxD = 0; nIdxD < nNbParams; nIdxD++)
            std::cout << "  index = " << nIdxD 
                      << ", value = " << ppdSensitivities[nIdxOption][nIdxD] 
                      << std::endl;
        }
      }
      
      std::cout << std::endl;
    }

    std::cout << std::endl;

    // Compute by finite difference in case of sensitivity by adjoint
    if (    pNumOutput->HasSensitivities() 
         && forwardOutput->IsSensitivityOnObjectif() )
    {
      // use a clone of the model
      hg::TheoreticalModel myModel(*pModel);

      shared_ptr<ComputationalFlags> pFlagsTmp(new ComputationalFlags);

      myModel.SetExternalFlags(pFlagsTmp);

      std::vector<double> pdMarketPrices;
      std::list< shared_ptr<Option> >::const_iterator i;
      for ( i = optionList.begin(); i != optionList.end(); ++i)
        pdMarketPrices.push_back( (*i)->GetMarketPrice() );

      const double dShift = 1.e-8;

      hg::Translator translator( *myModel.GetUnderlyingProcess() );
      std::vector<double> parameters( translator.GetParameters() );
      
      std::vector<double> sensitivities( parameters.size() );
      for (size_t nIdx = 0; nIdx < parameters.size(); nIdx++)
      {
        // Setup the perturbed parameters
        std::vector<double> newParameters(parameters);
        newParameters[nIdx] += dShift;

        // Get the perturbed underlying process
        shared_ptr<hg::UnderlyingProcess> pUP = translator(newParameters);

        myModel.SetUnderlyingProcess(pUP);

        shared_ptr<hg::MultiOutput> 
          forwardOutputTmp( hg::PriceForwardOption(myModel, *forwardOpt) );
      
        sensitivities[nIdx] = (forwardOutputTmp->GetObjectif() - dObjectif)
                            / dShift;
      }

      for (size_t nIdxD = 0; nIdxD < sensitivities.size(); nIdxD++)
        std::cout << sensitivities[nIdxD] << std::endl;
    }

    std::list< shared_ptr<Option> >::const_iterator iterOptions;

    std::cout << "Option prices:" << std::endl;
    size_t nCounter = 0;
    for (iterOptions = optionList.begin();
         iterOptions != optionList.end();
         ++iterOptions)
    {

      shared_ptr<finance::ModelOutput> output = pModel->Compute(**iterOptions);

      double dPrice = output->GetPrice();

      
      std::cout << "T = " << (**iterOptions).GetMaturityDate() 
                << ", K = " << (**iterOptions).GetStrike();

      if ( (**iterOptions).GetOptionType() == Option_Call )
        std::cout << ", (c)";
      else
        std::cout << ", (p)";

      std::cout << ", price = " << dPrice 
                << ", fwd price = " << pdPrices[nCounter++]
                << std::endl;      

      // Output sensitivities, if computed
      shared_ptr<hg::NumOutput>
        pNumOutput( dynamic_pointer_cast<hg::NumOutput>
                    ( output->GetNumOutput() ) );
    
      if ( pNumOutput->HasSensitivities() )
      {
        std::cout << "Sensitivities:" << std::endl;
        std::vector<double> pdSensitivities = pNumOutput->GetSensitivities();
        for (size_t nIdx = 0; nIdx < pdSensitivities.size(); nIdx++)
          std::cout << "  index = " << nIdx << ", value = " << pdSensitivities[nIdx] << std::endl;

        std::cout << std::endl;
      }

    }

    return 0;
  }
  catch ( const ::ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
              << e.GetFullMessage() << std::endl;

    return 1;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught.\n";

    return 2;
  }
}
