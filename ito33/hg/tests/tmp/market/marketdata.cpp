#include <iostream>
#include <fstream>
#include <cmath>

#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"
#include "ito33/autoptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/ratedata.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/finance/option.h"
#include "ito33/finance/cds.h"
#include "ito33/finance/eds.h"

#include "ito33/hg/underlyingprocess.h"
#include "hg/translator.h"

#include "marketdata.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::hg;

extern const ito33::Error ITO33_UNEXPECTED;

std::string GetToken(std::string& sBuffer)
{
  std::string sTemp;
  size_t nPos = sBuffer.find(",", 0);
  if ( nPos == std::string::npos )
  {
    sTemp = sBuffer;
    sBuffer.erase();
  }
  else
  {
    sTemp = sBuffer.substr(0, nPos);
    sBuffer.erase(0, nPos + 1);
  }

  return sTemp;
}

std::list< shared_ptr<Option> > ReadOptionData(std::string sFilename,
                                              shared_ptr<SessionData> pSessionData)
{
  // return list
  std::list< shared_ptr<Option> > optionList;

  // Open the file
  std::ifstream dataFile;
  dataFile.open(sFilename.c_str(), std::ios::in);

  if ( dataFile.fail() )
  {
    std::cout << "Cannot open file " << sFilename << std::endl;
    exit(1);
  }

  // temp string used to tokenize each line
  std::string temp;

  // Read each line separately
  std::string sBuffer;
  while ( getline(dataFile, sBuffer) )
  {

    // Filter out blank lines
    if (sBuffer[0] == NULL)
      continue;

    // First token is the (Bloomberg) description. Not needed.
    temp = GetToken(sBuffer);

    // Next token is valuation date
    temp = GetToken(sBuffer);    
    Date valuationDate(temp.c_str(), "%d/%m/%Y");

    // Next is the market price
    temp = GetToken(sBuffer);
    double dMarketPrice = atof(temp.c_str());

    // European or American
    ExerciseType exerciseType;
    temp = GetToken(sBuffer);
    if ( strcmp(temp.c_str(), "European") == 0 )
      exerciseType = ExerciseType_European;
    else
      exerciseType = ExerciseType_American;

    // Maturity date
    temp = GetToken(sBuffer);    
    Date maturityDate(temp.c_str(), "%d/%m/%Y");

    // Call or put
    OptionType optionType;
    temp = GetToken(sBuffer);
    if ( strcmp(temp.c_str(), "Call") == 0 )
      optionType = Option_Call;
    else
      optionType = Option_Put;

    // Finally, the strike
    temp = GetToken(sBuffer); 
    double dStrike = atof(temp.c_str());

    // Create the option and add to list
    shared_ptr<Option> pOption(new 
    Option(dStrike, maturityDate, optionType, exerciseType) );
    pOption->SetMarketPrice(dMarketPrice);
    pOption->SetSessionData(pSessionData);

    optionList.push_back(pOption);
    
  } // loop over lines in data file

  return optionList;
}

std::list< shared_ptr<Option> > ReadOptionData2(std::string sFilename,
                                               shared_ptr<SessionData> pSessionData)
{
  // return list
  std::list< shared_ptr<Option> > optionList;

  // Open the file
  std::ifstream dataFile;
  dataFile.open(sFilename.c_str(), std::ios::in);

  if ( dataFile.fail() )
  {
    std::cout << "Cannot open file " << sFilename << std::endl;
    exit(1);
  }

  // temp string used to tokenize each line
  std::string temp;

  // Read each line separately
  std::string sBuffer;
  while ( getline(dataFile, sBuffer) )
  {

    // Filter out blank lines
    if (sBuffer[0] == NULL)
      continue;

    // First token is the description. Not needed.
    temp = GetToken(sBuffer);

    // Next token is valuation date
    temp = GetToken(sBuffer);    
    Date valuationDate(temp.c_str(), "%m/%d/%Y");

    // Next is the spot price (not needed, but must be parsed)
    temp = GetToken(sBuffer);

    // Maturity date
    temp = GetToken(sBuffer);    
    Date maturityDate(temp.c_str(), "%m/%d/%Y");

    // Next is the forward price???? Not used
    temp = GetToken(sBuffer);

    // Next is the implied vol
    temp = GetToken(sBuffer); 
    double dStrike = atof(temp.c_str());

    // Finaly, the implied vol
    temp = GetToken(sBuffer); 
    double dImpliedVol = atof(temp.c_str());

    // Unspecified data
    ExerciseType exerciseType = ExerciseType_European;
    OptionType optionType = Option_Call;

    // Create the option and add to list
    shared_ptr<Option> pOption(new 
    Option(dStrike, maturityDate, optionType, exerciseType) );
    pOption->SetSessionData(pSessionData);
    //pOption->SetMarketPrice(dImpliedVol);
    pOption->SetImpliedVol(dImpliedVol);
    
    optionList.push_back(pOption);
    
  } // loop over lines in data file

  return optionList;
}



shared_ptr<SessionData> InitSessionDataMarket(Company company)
{

  double dSpot;  
  double dDividendAmount;
  Date dividendDate;
  double dDividendRate;

  // Spot and dividend info comes from finance.yahoo.com (found on
  // April 28, 2005).  Spot info is consistent with the data provided by 
  // polygon.  Note also that stock symbols are given with the polygon data.
  switch (company)
  {
  case Accor:
    dSpot = 38.31;
    dividendDate = Date(2005, Date::May, 17);
    dDividendAmount = 1.05;
    // might become 1.30 (proposed exceptional adds 0.25)
    dDividendRate = 0.025;
    break;

  case Carrefour:
    dSpot = 40.18;
    dividendDate = Date(2005, Date::Apr, 2);
    dDividendAmount = 0.94;
    dDividendRate = 0.025;
    break;

  case FranceTelecom:
    dSpot = 22.70;
    dividendDate = Date(2005, Date::Jun, 3);
    dDividendAmount = 0.48;
    dDividendRate = 0.02;
    break;

  case Lafarge:
    dSpot = 76.95;
    dividendDate = Date(2005, Date::Jul, 1);
    dDividendAmount = 2.40;
    dDividendRate = 0.03;
    break;

  case Valeo:
    dSpot = 35.58;
    dividendDate = Date(2005, Date::May, 16);
    dDividendAmount = 1.10;
    dDividendRate = 0.03;
    break;

  case BNBParibas:
    dSpot = 55.10;
    dividendDate = Date(2005, Date::May, 30);
    dDividendAmount = 2.00;
    dDividendRate = 0.035;
    break;

  case Endesea:
    dSpot = 17.06;
    dividendDate = Date(2005, Date::Jul, 1);
    dDividendAmount = 0.74;
    dDividendRate = 0.04;
    break;

  //case Iberdrola:
  //  dSpot = 19.57;
  //  dividendDate = Date(2005, Date::Jul, 1);
  //  dDividendAmount = 0.77;
  //  dDividendRate = 0.04;
  //  break;

  case Renault:
    dSpot = 68.75;
    dividendDate = Date(2005, Date::May, 13);
    dDividendAmount = 1.8;
    dDividendRate = 0.025;
    break;

  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Company not known");
    break;
  }


  // market data reported on march 10, 2005
  Date valuationDate(2005, Date::Mar, 10);

  shared_ptr<Numeraire> pCurrency( new Numeraire("EUR") );
  shared_ptr<Equity> pEquity(new Equity(dSpot, pCurrency));

  // What to do for dividends?  Ans: A dividend yield was computed
  // by dividing the dividend amount by the stock price, and rounding
  // to the nearest 0.005 (assuming that the companies aim for a constant
  // dividend yield, while the stock price fluctuates). This dividend rate
  // will be used after the discrete dividend.
  shared_ptr<Dividends> pDividends(new Dividends());
  pDividends->AddCash(dividendDate, dDividendAmount);
  pEquity->SetDividends(pDividends);

  shared_ptr<YieldCurveAnnuallyCompounded> pyf(new
    YieldCurveAnnuallyCompounded( valuationDate, 3));

  size_t nNbDays = dividendDate.GetExcel() - valuationDate.GetExcel();
  pyf->AddLeg( nNbDays-1, 0.0);
  pyf->AddLeg( nNbDays, dDividendRate);
  pyf->AddLeg( 20*365, dDividendRate);

  //size_t nNbDays = dividendDate.GetExcel() - valuationDate.GetExcel();
  //pyf->AddLeg( 1, dDividendRate);
  //pyf->AddLeg( nNbDays, dDividendRate);
  //pyf->AddLeg( 20*365, dDividendRate);
  
  //shared_ptr<YieldCurve> pyf(new YieldCurveFlat(0.0));
     
  pEquity->SetBorrowCurve(pyf);
    
  // US yield curve, april 1, 2005
  // http://futures.fxstreet.com/Futures/content/100750/content.asp?data=3/27/2005
  //  1m    3m    6m    1y    2y    3y    5y    7y   10y    20y
  // 2.66  2.80  3.13  3.34  3.75  3.90  4.13  4.29  4.46  4.85

  // European yield curve
  // http://epp.eurostat.cec.eu.int/cache/ITY_PUBLIC/EYC/EN/page4.htm#yields
  //   1y    2y    3y    4y    5y    6y    7y    8y    9y    10y
  //  2.11  2.30  2.50  2.69  2.87  3.03  3.18  3.31  3.43  3.53
  
  //double dRate = 0.03;
  //shared_ptr<YieldCurve> pyc = new YieldCurveFlat(dRate);

  shared_ptr<YieldCurveAnnuallyCompounded> pyc(new
    YieldCurveAnnuallyCompounded( Date(2005, Date::Mar, 1), 10));

  // ignore leap years
  pyc->AddLeg(1*365, 0.0211);
  pyc->AddLeg(2*365, 0.0230);
  pyc->AddLeg(3*365, 0.0250);
  pyc->AddLeg(4*365, 0.0269);
  pyc->AddLeg(5*365, 0.0287);
  pyc->AddLeg(6*365, 0.0303);
  pyc->AddLeg(7*365, 0.0318);
  pyc->AddLeg(8*365, 0.0331);
  pyc->AddLeg(9*365, 0.0343);
  pyc->AddLeg(10*365, 0.0353);

  shared_ptr<RateData> pRateData( new RateData() );
  pRateData->SetYieldCurve(pCurrency, pyc );

  shared_ptr<SessionData> pSessionData(new SessionData(pRateData, pEquity, valuationDate));

  return pSessionData;
}


std::string GetInputFile(Company company, bool bIsAdi)
{
  std::string sInputFile;

  switch (company)
  {
  case Accor:
    if (bIsAdi)
      sInputFile = "accor_adi.txt";
    else
      sInputFile = "accor.euro";

    break;

  case Carrefour:
    if (bIsAdi)
      sInputFile = "carrefour_adi.txt";
    else
      sInputFile = "carrefour.txt";
    
    break;

  case FranceTelecom:
    if (bIsAdi)
      sInputFile = "francetelecom_adi.txt";
    else
      sInputFile = "francetelecom.euro";
    
    break;

  case Lafarge:
    if (bIsAdi)
      sInputFile = "lafarge_adi.txt";
    else
      sInputFile = "lafarge.euro";

    break;

  case Valeo:
    if (bIsAdi)
      sInputFile = "valeo_adi.txt";
    else
      sInputFile = "valeo.txt";
    
    break;

  case BNBParibas:
    if (bIsAdi)
      sInputFile = "bnbparibas_adi.txt";
    else
      sInputFile = "bnbparibas.txt";
    
    break;

  case Endesea:
    if (bIsAdi)
      sInputFile = "endesea_adi.txt";
    else
      sInputFile = "endesea.txt";
    
    break;

  //case Iberdrola:
  //  sInputFile = "iberdrola.txt";
  //  break;

  case Renault:
    if (bIsAdi)
      sInputFile = "renault_adi.txt";
    else
      sInputFile = "renault.euro";
    
    break;

  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Company not known");
    break;
  }

  return sInputFile;
}


shared_ptr<hg::UnderlyingProcess> GetUnderlyingProcess(Company company)
{

  // Setup model structure
  size_t nNbRegimes = 2;
  std::vector<double> pdVols(nNbRegimes);
  std::vector<double> pdIntensities(nNbRegimes);

  size_t nIdx;
  for (nIdx = 0; nIdx < nNbRegimes; nIdx++)
  {
    pdVols[nIdx] = 0.1;
    pdIntensities[nIdx] = 0.001;

    //pdVols[nIdx] = 0.0;
    //pdIntensities[nIdx] = 0.0;
  }

  shared_ptr<hg::UnderlyingProcess> 
    pUnderlyingProcessInitial( new hg::UnderlyingProcess
                                   (nNbRegimes, pdVols, pdIntensities) );

  size_t nIdx1, nIdx2;
  for (nIdx1 = 0; nIdx1 < nNbRegimes; nIdx1++)
  {
    for (nIdx2 = 0; nIdx2 < nNbRegimes; nIdx2++)
    {
      if (nIdx1 == nIdx2)
        continue;

      //if (nIdx1 > 0)
      //  continue;

      // Make new jumps with default values
      Jumps newJumps;
      newJumps.push_back( Jump(0.01, -0.01) );
      //newJumps.push_back( Jump(0.2, -0.4) );      

      pUnderlyingProcessInitial->SetJumps(nIdx1, nIdx2, newJumps);
    }
  }

  // Use a translator to translate a known good guess into a new
  // underlying process
  AutoPtr<Translator> pTranslator(new Translator(*pUnderlyingProcessInitial));

  // By default, return what was constructed above
  shared_ptr<hg::UnderlyingProcess> pReturnProcess = pUnderlyingProcessInitial;

  // vector of parameters
  size_t nNbParams = pTranslator->GetNbParameters();
  std::vector<double> pdX( nNbParams );

  switch (company)
  {
  case Accor:
    /*
    if ( nNbParams == 8 )
    {
      // assume 2 regimes, no internal jumps
      pdX[0] = 0.141374;
      pdX[1] = 0.493682;
      pdX[2] = 0.0;
      pdX[3] = 0.0;
      pdX[4] = 0.243057;
      pdX[5] = -0.153369;
      pdX[6] = 2.0;
      pdX[7] = -0.163968;
    }
    else if ( nNbParams == 12 )
    {
      pdX[0] = 1.7939e-001;
      pdX[1] =  1.0000e-006;     
      pdX[2] =  1.0000e-006; 
      pdX[3] =  1.6340e-002; 
      pdX[4] =  1.0229e+000; 
      pdX[5] = -1.1025e-001; 
      pdX[6] =  3.4061e-001; 
      pdX[7] = -1.8322e-002; 
      pdX[8] =  0.0000e+000;  
      pdX[9] = -9.9000e-001;  
      pdX[10] =  6.3081e-002;   
      pdX[11] = -9.9000e-001;    

      //pdX[0] = 0.0000e+000;
      //pdX[1] = 0.0000e+000;
      //pdX[2] = 0.0000e+000;
      //pdX[3] = 4.0440e-002;
      //pdX[4] = 5.4683e-001;
      //pdX[5] = -3.1735e-001;
      //pdX[6] = 0.0000e+000;
      //pdX[7] = -3.4521e-001;
      //pdX[8] = 0.0000e+000;
      //pdX[9] = -3.3122e-001;
      //pdX[10] = 0.0000e+000;
      //pdX[11] = -4.7332e-001;
    }
    else if ( nNbParams == 18 )
    {
      // assume 3 regimes, no internal jumps
      pdX[0] = 0.1311191;
      pdX[1] = 0.0;
      pdX[2] = 0.01546597;
      pdX[3] = 0.001088863;
      pdX[4] = 0.6813949;
      pdX[5] = 0.0;
      pdX[6] = 0.0;          // 0 -> 1
      pdX[7] = -0.2360307;   // 0 -> 1
      pdX[8] = 0.2606002;    // 0 -> 2
      pdX[9] = -0.1192591;   // 0 -> 2
      pdX[10] = 0.6221086;   // 1 -> 0
      pdX[11] = -0.99;       // 1 -> 0
      pdX[12] = 0.0;         // 1 -> 2
      pdX[13] = 0.3357349;   // 1 -> 2
      pdX[14] = 0.6966612;   // 2 -> 0
      pdX[15] = 0.5173301;   // 2 -> 0
      pdX[16] = 0.03479919;  // 2 -> 1
      pdX[17] = -0.4081513;  // 2 -> 1
    }
    
    pReturnProcess = (*pTranslator)(&pdX[0]);
    */
    break;

  case Carrefour:
    break;

  case FranceTelecom:
    break;

  case Lafarge:
    break;

  case Valeo:
    break;

  case BNBParibas:
    break;

  case Endesea:
    break;

  //case Iberdrola:
  //  break;

  case Renault:
    break;

  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Company not known");
    break;
  }

  return pReturnProcess;

}

std::list< shared_ptr<CDS> > GetCDSData(Company company, 
                                       shared_ptr<SessionData> pSessionData)
{

  // The return list
  std::list< shared_ptr<CDS> > listCDS;

  // The market cds apread (ask)
  double dSpread3year;
  double dSpread5year;

  switch (company)
  {
  case Accor:
    dSpread3year = 24.0;
    //dSpread3year = 29.0;
    dSpread5year = 43.0;
    break;

  case Carrefour:
    dSpread3year = 12.0;
    dSpread5year = 18.0;
    break;

  case FranceTelecom:
    dSpread3year = 13.0;
    dSpread5year = 27.0;
    break;

  case Lafarge:
    dSpread3year = 20.0;
    dSpread5year = 34.0;
    break;

  case Valeo:
    dSpread3year = 37.0;
    dSpread5year = 54.0;
    break;

  case BNBParibas:
    dSpread3year = 6.0;
    dSpread5year = 9.0;
    break;

  case Endesea:
    dSpread3year = 18.0;
    dSpread5year = 26.0;
    break;

  //case Iberdrola:
  //  dSpread3year = 13.0;
  //  dSpread5year = 18.0;
  //  break;

  case Renault:
    dSpread3year = 17.0;
    dSpread5year = 30.0;
    break;

  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Company not known");
    break;
  }

  // The spreads are quoted in basis points
  dSpread3year /= 10000.0;
  dSpread5year /= 10000.0;

  Date valuationDate = pSessionData->GetValuationDate();
  Date issueDate = valuationDate;
  Date firstDate = valuationDate;
  firstDate.AddMonths(3);

  // Create the 3 year cds
  Date lastDate = valuationDate;
  lastDate.AddYears(3);
 
  double dRecoveryRate = 0.5;

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread3year,
                          Date::DayCountConvention_Act365,
                          //Frequency_SemiAnnual
                          Frequency_Quarterly
                        )
                  );

  shared_ptr<finance::CDS>
    pCDS( new ito33::finance::CDS(dRecoveryRate, pSpreadStream) );

  //shared_ptr<finance::CDS> pCDS(new ito33::finance::CDS(issueDate, 
  //                               3*12, 
  //                               dSpread3year, 
  //                               Date::DayCountConvention_Act360, 
  //                               Frequency_Quarterly, 
  //                               dRecoveryRate));

  pCDS->SetSessionData(pSessionData);
  pCDS->SetMarketPrice(0.0);

  listCDS.push_back(pCDS);

  // create the 5 year cds
  lastDate.AddYears(2);

  pSpreadStream = make_ptr( new CashFlowStreamUniform
                                (
                                  issueDate,
                                  firstDate,
                                  lastDate,
                                  dSpread5year,
                                  Date::DayCountConvention_Act365,
                                  //Frequency_SemiAnnual
                                  Frequency_Quarterly
                                ) );

  pCDS = make_ptr( new CDS(dRecoveryRate, pSpreadStream) );
 
  pCDS->SetSessionData(pSessionData);
  pCDS->SetMarketPrice(0.0);

  listCDS.push_back(pCDS);

  return listCDS;

}


std::list< shared_ptr<EDS> > GetEDSData(Company company, 
                                       shared_ptr<SessionData> pSessionData)
{

  // The return list
  std::list< shared_ptr<EDS> > listEDS;

  // The market eds apread (ask)
  double dSpread3year;
  double dSpread5year;

  switch (company)
  {
  case Accor:
    dSpread3year = (49.0 + 61.0)/2.0;
    //dSpread3year = 64.0;
    dSpread5year = (78.0 + 99.0)/2.0;
    break;

  case Carrefour:
    dSpread3year = (53.0 + 66.0)/2.0;
    dSpread5year = (74.0 + 92.0)/2.0;
    break;

  case FranceTelecom:
    dSpread3year = (69.0 + 84.0)/2.0;
    dSpread5year = (98.0 + 119.0)/2.0;
    break;

  case Lafarge:
    dSpread3year = (60.0 + 74.0)/2.0;
    dSpread5year = (90.0 + 110.0)/2.0;
    break;

  case Valeo:
    dSpread3year = 55.0;  // bid
    dSpread5year = 78.0;  // bid
    break;

  case BNBParibas:
    dSpread3year = (74.0 + 88.0)/2.0;
    dSpread5year = (117.0 + 138.0)/2.0;
    break;

  case Endesea:
    dSpread3year = (35.0 + 46.0)/2.0;
    dSpread5year = (53.0 + 68.0)/2.0;
    break;

  //case Iberdrola:
  //  dSpread3year = 0.0; // no data
  //  dSpread5year = 0.0; // no data
  //  break;

  case Renault:
    dSpread3year = (58.0 + 71.0)/2.0;
    dSpread5year = (81.0 + 100.0)/2.0;
    break;

  default:
    throw EXCEPTION_MSG(ITO33_UNEXPECTED, "Company not known");
    break;
  }
  
  // The spreads are quoted in basis points
  dSpread3year /= 10000.0;
  dSpread5year /= 10000.0;

  Date valuationDate = pSessionData->GetValuationDate();
  Date issueDate = valuationDate;
  Date firstDate = valuationDate;
  firstDate.AddMonths(3);

  // Create the 3 year eds
  Date lastDate = valuationDate;
  lastDate.AddYears(3);
 
  double dRecoveryRate = 0.5;
  double dBarrier = 0.3 * pSessionData->GetSpotSharePrice();

  shared_ptr<CashFlowStreamUniform>
    pSpreadStream( new CashFlowStreamUniform
                        (
                          issueDate,
                          firstDate,
                          lastDate,
                          dSpread3year,
                          Date::DayCountConvention_Act365,
                          //Frequency_SemiAnnual
                          Frequency_Quarterly
                        )
                  );

  shared_ptr<finance::EDS>
    pEDS( new ito33::finance::EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  pEDS->SetMarketPrice(0.0);

  listEDS.push_back(pEDS);

  // create the 5 year cds
  lastDate.AddYears(2);

  pSpreadStream = make_ptr( new CashFlowStreamUniform
                                (
                                  issueDate,
                                  firstDate,
                                  lastDate,
                                  dSpread5year,
                                  Date::DayCountConvention_Act365,
                                  //Frequency_SemiAnnual
                                  Frequency_Quarterly
                                ) );

  pEDS = make_ptr( new EDS(dRecoveryRate, pSpreadStream, dBarrier) );

  pEDS->SetSessionData(pSessionData);
  pEDS->SetMarketPrice(0.0);

  listEDS.push_back(pEDS);

  return listEDS;

}

