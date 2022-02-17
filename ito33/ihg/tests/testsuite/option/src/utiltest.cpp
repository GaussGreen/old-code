/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/option/testsuite/utiltest.cpp
// Purpose:     util function for tests
// Author:      Yann d'Halluin David Pooley
// Created:     2004/05/02
// RCS-ID:      $Id: utiltest.cpp,v 1.5 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/tests/option/testsuite/utiltest.cpp
    @brief set of functions to help when writting tests results
**/

#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"

#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/sharedptr.h"
#include "ito33/autoptr.h"
#include "ito33/array.h"

#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/volatility.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/hazardrateflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratedecay.h"
#include "ito33/ihg/hazardratelinear.h"
#include "ito33/ihg/hazardratepower.h"

#include "ito33/numeric/numparams.h"
#include "ito33/numeric/meshparams.h"

#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/yieldcurve_annuallycompounded.h"
#include "ito33/finance/dividends.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/option.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/numeraire_list.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/ratedata.h"

#include "ito33/pricing/option.h"
#include "ito33/pricing/optionparams.h"

#include "ihg/model.h"
#include "ihg/optionpricer.h"


#include "utiltest.h"
#include "testparam.h"


namespace ito33
{
  
namespace XML 
{
      
  template <> Tag MakeNamedTag(const char *name, 
    const ParameterOneValue &parameterOneValue, Tag& parent)  
  {  
    Tag tag(name, parent);
     
    std::vector<double>::const_iterator iterParam;
    std::vector<double>::const_iterator iterValue;

    for ( iterParam  = parameterOneValue.m_pdParameter.begin(), 
          iterValue  = parameterOneValue.m_pdValue.begin();
          iterParam != parameterOneValue.m_pdParameter.end(); ++iterParam, ++iterValue )
    {
       tag.Element("value").Attr("paramval",(*iterParam)).Attr("val",(*iterValue));
    }

   return tag;
  }


       
  template <> Tag MakeNamedTag(const char *name, 
    const ParameterTwoValues &parameterTwoValues,  Tag& parent)  
  {
    Tag tag(name, parent);
     
    size_t nIdx = 0;
    size_t nSize = parameterTwoValues.m_pdParameter.size();

    std::vector<double>::const_iterator iterParam;
    std::vector<double>::const_iterator iterValue1;
    std::vector<double>::const_iterator iterValue2;

    for ( nIdx = 0; nIdx < nSize; nIdx++)
    {
      double dParamVal = parameterTwoValues.m_pdParameter[nIdx];
      double dVal1     = parameterTwoValues.m_pdValue1[nIdx];
      double dVal2     = parameterTwoValues.m_pdValue2[nIdx];

      double dError = fabs(dVal1/dVal2 - 1.);
      if ( dVal1 < 1.0 )
        dError = fabs( dVal1 - dVal2 );

      tag.Element("value").Attr("paramval",dParamVal).Attr("val1",dVal1).Attr("val2",dVal2).Attr("error",dError);
    }

    return tag;    
  }

  } //end XML

namespace ihg
{

namespace test
{

//-----------------------------------------------------------------------------
// Generate several files to test that the
// xml real works fine
// all the files created are installed ino
// ihg/tests/testsuite/testing
//-----------------------------------------------------------------------------
void GenerateXMLTestingFiles(size_t &nFileCounter)
{
  std::vector<int> year;   //number of years to add to the pricing date
  std::vector<int> day;    // number of days to add to the pricing date
  std::vector<int> month;  // number of  months to add to the pricing date
  char fileName[1024];     // allocate space for the file name
  nFileCounter = 0;        // file ihg{FileCounter}.xml
  
  
  Date pricingDate(2003, Date::Feb, 1); //pricing date

  std::vector<double> spot; //spot price
  std::vector<double> strike; // strike price
  
  std::vector<finance::ExerciseType> exerciseType; //American,European
  std::vector<finance::OptionType>   optionType; //call,put ....
  std::vector<HazardRateType>               hazardRateType; //DECAY,FLAT,LINEAR,POLYNOM,TIMEONLY;
  std::vector<YieldType>                    yieldType; //{FLAT, ZEROCOUPON}
  std::vector<YieldType>                    foreignType;//{FLAT,ZEROCOUPON}
  std::vector<DividendType>                 dividendType;// {NODIVIDEND,CASH,YIELD};
  std::vector<VolatilityType>               volatilityType; //VolatilityType {FLATVOL};

  //Create the different vectors
  year.push_back(1);
  year.push_back(2);
  day.push_back(0);
  month.push_back(6);


  spot.push_back(100);
  spot.push_back(50);
  spot.push_back(150);
  
  strike.push_back(100);
  strike.push_back(150);
  strike.push_back(50);
    
  exerciseType.push_back(finance::ExerciseType_American);
  exerciseType.push_back(finance::ExerciseType_European);

  optionType.push_back(finance::Option_Call);
  optionType.push_back(finance::Option_Put);

  hazardRateType.push_back(FLATHR);
  //hazardRateType.push_back(POWERHR);
  // hazardRateType.push_back(TIMEONLYHR);

  yieldType.push_back(FLAT);
  foreignType.push_back(FLAT);
  //foreignType.push_back(ZEROCOUPON);

  dividendType.push_back(NODIVIDEND);
  //dividendType.push_back(CASH);

  volatilityType.push_back(FLATVOL);





 
  size_t nFiles = year.size()*day.size()*month.size()*spot.size()*strike.size()
    *exerciseType.size()*optionType.size()*hazardRateType.size()
    *yieldType.size()*foreignType.size()*dividendType.size()
    *volatilityType.size();

  std::cout <<"Generating " << nFiles << " file(s) " << std::endl;
  
  try
  {
  size_t nIdxYear;

  for ( nIdxYear = 0 ; nIdxYear < year.size() ; nIdxYear++ )
    { //loop over the years
     size_t nIdxMonth;
     
     for ( nIdxMonth = 0; nIdxMonth < month.size(); nIdxMonth++)
     { //loop over the month
       size_t nIdxDay;
       
       for ( nIdxDay = 0 ; nIdxDay < day.size(); nIdxDay++)
       {//loop over the day
          size_t nIdxSpot;
          
          for ( nIdxSpot = 0; nIdxSpot < spot.size(); nIdxSpot++)  
          { //loop over the spot    
            size_t nIdxYield;

            for ( nIdxYield = 0; nIdxYield < yieldType.size() ; nIdxYield++)
            { //loop over yield curve
              size_t nIdxForeign; 

            for ( nIdxForeign = 0; nIdxForeign < foreignType.size(); 
              nIdxForeign++)
              { //loop over foreign curve
                size_t nIdxStrike;
                 
                for ( nIdxStrike = 0; nIdxStrike < strike.size();  
                  nIdxStrike++)
                  { //loop over strike
                  size_t nIdxOptionType;

                  for ( nIdxOptionType = 0; nIdxOptionType< optionType.size();
                    nIdxOptionType++)
                  { //loop over option type
                    size_t nIdxExerciseType;

                    for ( nIdxExerciseType = 0; nIdxExerciseType< exerciseType.size();
                      nIdxExerciseType++)        
                    { //loop over exercise type
                     size_t nIdxDividend;

                     for ( nIdxDividend = 0 ; nIdxDividend < dividendType.size();
                       nIdxDividend++)
                     {//loop over dividends

                      size_t nIdxVolatility;

                     
                      for ( nIdxVolatility = 0; nIdxVolatility< volatilityType.size();
                        nIdxVolatility++)
                      { //loop over volatility
                        size_t nIdxHazardRate;

                        for ( nIdxHazardRate = 0; nIdxHazardRate< hazardRateType.size(); 
                          nIdxHazardRate++)                 
                        { //loop over hazadrate
                          //generate a unique file name; 
                          fileName[0]=0;          
                          sprintf(fileName,"xmlfiles/ihg_option_temp%d.xml",nFileCounter);      
                          nFileCounter++;
              
                          //std::cout << "Creating File " << fileName << std::endl;

                          shared_ptr<finance::Numeraire>
                            pNumeraire(new finance::Numeraire("EUR"));

                          shared_ptr<finance::Equity> 
                            pEquity(new finance::Equity(spot[nIdxSpot], 
                                                        pNumeraire));

                          //avoid to have a borrow curve and a dividend
                          //at the same time
                          if ( foreignType[nIdxForeign] == NOYIELD )
                          {
                          pEquity->SetDividends( 
                            GetDividend(dividendType[nIdxDividend],pricingDate,spot[nIdxSpot]));
                          }
                          else
                          {
                          pEquity->SetDividends( 
                            GetDividend(NODIVIDEND,pricingDate,spot[nIdxSpot]));
                          }
                          pEquity->SetBorrowCurve(
                            GetForeignCurve(foreignType[nIdxForeign],pricingDate) );

                          shared_ptr<finance::RateData> pRateData(new finance::RateData);

                          pRateData->SetYieldCurve(pNumeraire, 
                            GetYieldCurve(yieldType[nIdxYield],pricingDate));
                                                   
                          shared_ptr<finance::SessionData>
                            pSessionData( new finance::SessionData(pRateData, 
                                                                   pEquity, 
                                                                   pricingDate) );

                          Date maturityDate(pricingDate);
                          maturityDate.AddDays(day[nIdxDay]);
                          maturityDate.AddYears(year[nIdxYear]);
                          maturityDate.AddMonths(day[nIdxMonth]);

                          shared_ptr<finance::Option> 
                            opt(new finance::Option(
                                     strike[nIdxStrike], 
                                     maturityDate,
                                     optionType[nIdxOptionType],
                                     exerciseType[nIdxExerciseType]
                                     ));
                          opt->SetMarketPrice(4.35);
                          opt->SetSessionData(pSessionData);

                          shared_ptr<ihg::TheoreticalModel>
                            pModel(new ihg::TheoreticalModel);



                          pModel->SetVolatility( 
                            GetVolatility(volatilityType[nIdxVolatility]) );

                          pModel->SetHazardRate(
                            GetHazardRate(hazardRateType[nIdxHazardRate],spot[nIdxSpot],pricingDate));

                          pModel->SetDebugOutputFile(fileName);
                          
                          shared_ptr<finance::ModelOutput> 
                            output = pModel->Compute(*opt);
                          //std::cout.precision(10);
                          //std::cout << "The option price is: " << output->GetPrice() << std::endl;

                        } //end hazard rate
                      }//end volatility
                     } //end dividend
                    } //end exercise type
                  } // end option type
                } //end strike
              } //end foreign curve
            } //end yield curve
          } //end spot
       } //end day
     } //end month
   } //end year

  } //end try
  catch ( const ito33::Exception& e )
  {
    std::cerr << "Exception caught:\n"
      << e.GetFullMessage() << std::endl;
  }
  catch ( ... )
  {
    std::cerr << "Unexpected exception caught.\n";
  }


} //end GenerateXMLTestingFiles()


//-----------------------------------------------------------------------------
// Return a shared pointer to a particular hazard rate type.
//-----------------------------------------------------------------------------
shared_ptr<Volatility> GetVolatility(VolatilityType volatilityType)
{  
  if ( volatilityType == FLATVOL )
  {
    return shared_ptr<Volatility>(new VolatilityFlat(.2));
  }
  else
  {
    std::cerr << "No such type of volatility" << std::endl;
    exit(0);
  }//end if

  return shared_ptr<Volatility>();

} // shared_ptr<ihg::Volatility>

//-----------------------------------------------------------------------------
// return a shared pointer to a particular type of dividends
//-----------------------------------------------------------------------------
shared_ptr<finance::Dividends> 
GetDividend(DividendType dividendType, Date date,double Spot)
{
 // enum DividendType {NODIVIDEND,DIVIDENDS};

  if ( dividendType == NODIVIDEND) 
  {
    return shared_ptr<finance::Dividends>(new finance::Dividends());
  } 
  else if ( dividendType == CASH )
  {
    shared_ptr<finance::Dividends> div(new finance::Dividends());

    size_t nLeg = 10;
    size_t nIdx;
    int nMonthCash = 3;
    double dCash = 1./100. * Spot / double(nMonthCash) / 12; //3% of spot

    for (nIdx = 0 ; nIdx < nLeg; nIdx++)
    {
      Date loopDate = date;
      div->AddCash(loopDate.AddMonths(nIdx+nMonthCash), dCash);
    }

    return div;
  }
  else if ( dividendType == YIELD )
  {
     shared_ptr<finance::Dividends> 
      div(new finance::Dividends());

    size_t nLeg     = 10;
    size_t nIdx;
    int nMonthYield = 2;
    double dPercent = 1.5/100./12.;

    for (nIdx = 0 ; nIdx < nLeg; nIdx++)
    {
      Date loopDate = date;
      div->AddYield(loopDate.AddMonths(nIdx+nMonthYield),dPercent);
    }

    return div;

  }
  else 
  {
    std::cerr << "No such type of dividends" <<std::endl;
    exit(0);
  } //end if
} // shared_ptr<finance::Dividends>

//-----------------------------------------------------------------------------
// return a shared pointer to a yield curve
//-----------------------------------------------------------------------------
shared_ptr<finance::YieldCurve> 
GetYieldCurve(YieldType yieldType,Date date)
{
    //enum YieldType {FLAT, RATIO, FIXED};

  if ( yieldType == FLAT )
  {
    double flatRate = .05;
    return shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(flatRate));
  } 
  else if ( yieldType == ZEROCOUPON )
  {  
    size_t nLeg = 10;
    size_t nIdx;
    int nYearShift     = 0;
    double dRate       = .05;

    shared_ptr<finance::YieldCurveAnnuallyCompounded> 
      yfc(new finance::YieldCurveAnnuallyCompounded(date));

    for (nIdx = 0 ; nIdx < nLeg; nIdx++)
    {
      Date loopDate = date;
      yfc->AddLeg(int(nIdx+nYearShift)*INVERSEONEDAY,dRate/(double(nIdx)+1.));
    }
    
    return yfc;

  }
  else if ( yieldType == NOYIELD)
  {
    double flatRate = .0;
    return shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(flatRate));
  }
  else
  {
    std::cerr << "No such type of yield curve"<< std::endl;
    exit(0);
  } //end if 
} // shared_ptr<finance::YieldCurve>

//-----------------------------------------------------------------------------
// return a shared pointer to a foreign curve
//-----------------------------------------------------------------------------
shared_ptr<finance::YieldCurve> 
GetForeignCurve(YieldType foreignType,Date date)
{
  if ( foreignType == FLAT )
  {
    double flatRate = .01;
    return shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(flatRate));
  } 
  else if ( foreignType == ZEROCOUPON )
  { 
    size_t nLeg = 10;
    size_t nIdx;
    int nYearShift     = 0;
    double dRate       = .02;

    shared_ptr<finance::YieldCurveAnnuallyCompounded> 
      yfc(new finance::YieldCurveAnnuallyCompounded(date));

    for (nIdx = 0 ; nIdx < nLeg; nIdx++)
    {
      Date loopDate = date;
      yfc->AddLeg(int(nIdx+nYearShift)*INVERSEONEDAY,dRate/(double(nIdx)+1.));
    }
    
    return yfc;
  }
  else if ( foreignType == NOYIELD)
  {
    double flatRate = .0;
    return shared_ptr<finance::YieldCurve>(new finance::YieldCurveFlat(flatRate));
  }
  else
  {
    std::cerr << "No such type of foreign curve"<< std::endl;
    exit(0);
  } //end if
} //shared_ptr<finance::YieldCurve>

//-----------------------------------------------------------------------------
// Return the appropriate hazard rate shared ptr
//-----------------------------------------------------------------------------
shared_ptr<HazardRate> 
GetHazardRate(HazardRateType hazardRateType, double dSpot, Date ValuationDate)
{
  if ( hazardRateType == DECAYHR )
  {
    double dAlpha = .02;
    double dS0   = (1.5)*dSpot;
     
    return shared_ptr<HazardRate>( new HazardRateDecay(dAlpha,dS0) );
  }
  else if ( hazardRateType == FLATHR)
  {    
    double dValue = .02;
    return shared_ptr<HazardRate>( new HazardRateFlat(dValue) ); 
  }
  else if ( hazardRateType == POWERHR )
  {
    double dAlpha = .02;
    double dBeta  = 0.9;
    double dS0    = (1.5)*dSpot;
    return shared_ptr<HazardRate>( new HazardRatePower(dAlpha,dBeta,dS0) );
  }
  else if ( hazardRateType == TIMEONLYHR)
  {
    Date times[4] = { ValuationDate, ValuationDate.AddYears(1),
                      ValuationDate.AddYears(2), ValuationDate.AddYears(3) };
    double values[4] = { 0.01, .01, .01, .01 };
    return shared_ptr<HazardRate>(new HazardRateTimeOnly(times, values, 4));
  }
  else
  {
   std::cerr<< "No such Hazard rate type: " <<std::endl;
   exit(0);
  }
} //end shared_ptr<HazardRate>
} //end namespace test
} //end namespace ihg
} //endnamespace ito33
