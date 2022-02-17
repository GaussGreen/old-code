/////////////////////////////////////////////////////////////////////////////
// Name:        tests/bondlike/testmandatory.cpp
// Purpose:     test file mandatory
// Author:      Ito33 Canada
// Created:     May 11, 2005
// RCS-ID:      $Id: testmandatory.cpp,v 1.5 2006/08/19 23:22:40 wang Exp $
// Copyright:   (c) 2005- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------
#include "ito33/beforestd.h"
#include <iostream>
#include "ito33/afterstd.h"
#include "ito33/date.h"

#include "ito33/common.h"
#include "ito33/debug.h"
#include "ito33/exception.h"
#include "ito33/cppunit.h"
#include "ito33/sharedptr.h"

#include "ito33/finance/bondlike/pepsaveragingperiod.h"

#include "ito33/finance/cashflowstream_general.h"

#include "ito33/tests/testmandatory.h"

#include "ito33/tests/utilexml.h"

using namespace ito33;

using namespace ito33::finance;

  
void PepsAveragingPeriodTest::AveragingPeriodTooLarge()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2006, Date::Mar, 1) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );

}
  
void PepsAveragingPeriodTest::NegativeNumberOfSampling()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 0;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );
}

void PepsAveragingPeriodTest::TooManySamplesUsed()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod   
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );

 pAveragingPeriod->SetCurrentStockAverage(100., 35);
}


void PepsAveragingPeriodTest::SetNegativeStockAverage()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );

  pAveragingPeriod->SetCurrentStockAverage(-100., 10);
}

void PepsAveragingPeriodTest::SetNegativeConversionRatioAverage()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
    PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, avgEndDate, nSampling)
    );

  pAveragingPeriod->SetCurrentConversionRatioAverage(-1., 10);

}

void PepsAveragingPeriodTest::SetConversionRatioAveragingWhenStockAveragingPeriod()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );

  pAveragingPeriod->SetCurrentConversionRatioAverage(1., 5);

}
 
void PepsAveragingPeriodTest::SetStockAveragingWhenConversionAveragingPeriod()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
    PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, avgEndDate, nSampling)
    );

  pAveragingPeriod->SetCurrentStockAverage(100., 10);
}
 
void PepsAveragingPeriodTest::GetStockAveragingWhenConversionAveragingPeriod()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
    PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, avgEndDate, nSampling)
    );

  double dCurrentAvg = 0.;
  dCurrentAvg = pAveragingPeriod->GetCurrentStockAverage();
}
 
void PepsAveragingPeriodTest::GetConversionRatioAveragingWhenStockAveragingPeriod()
{
  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );

 double dCurrentAvg = 0.;
 dCurrentAvg = pAveragingPeriod->GetCurrentConversionRatioAverage();

}

void PepsAveragingPeriodTest::DumpStockAverage()
{
  std::ostringstream oss;

  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
     PEPSAveragingPeriod::CreateWithStock(avgStartDate, avgEndDate, nSampling)
    );
  
  pAveragingPeriod->SetCurrentStockAverage(100.,5);

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<peps_averaging_period>\n"
                "<start_date>2005-01-01</start_date>\n"
                "<end_date>2005-01-20</end_date>\n"
                "<stock_averaging>1</stock_averaging>\n"
                "<number_of_sampling_averages>20</number_of_sampling_averages>\n"
                "<current_average>100</current_average>\n"
                "<number_of_samples_used>5</number_of_samples_used>\n"
                "</peps_averaging_period>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

  pAveragingPeriod->Dump(root);

}

void PepsAveragingPeriodTest::DumpConversionRatioAverage()
{
  std::ostringstream oss;

  Date avgStartDate(2005, Date::Jan, 1);
  Date avgEndDate(2005, Date::Jan, 20) ; 
  size_t nSampling = 20;

  shared_ptr<PEPSAveragingPeriod> pAveragingPeriod 
    (
    PEPSAveragingPeriod::CreateWithConversionRatio(avgStartDate, avgEndDate, nSampling)
    );
  
  pAveragingPeriod->SetCurrentConversionRatioAverage(1.23, 10);

  ExpectedXML expected(oss,
                "<?xml version=\"1.0\"?>"
                "<root>\n"
                "<peps_averaging_period>\n"
                "<start_date>2005-01-01</start_date>\n"
                "<end_date>2005-01-20</end_date>\n"
                "<stock_averaging>0</stock_averaging>\n"
                "<number_of_sampling_averages>20</number_of_sampling_averages>\n"
                "<current_average>1.23</current_average>\n"
                "<number_of_samples_used>10</number_of_samples_used>\n"
                "</peps_averaging_period>\n"
                "</root>"
              );


  ito33::XML::RootTag root("root",oss);

  pAveragingPeriod->Dump(root);

}
