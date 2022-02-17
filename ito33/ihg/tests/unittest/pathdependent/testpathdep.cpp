/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/tests/common/pathdependent/testpathdep.cpp
// Purpose:     Acceptance test for path dependent 
// Author:      ITO 33
// Created:     04/02/2005
// RCS-ID:      $Id: testpathdep.cpp,v 1.11 2006/08/20 09:49:27 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/beforestd.h"
#include <vector>
#include "ito33/afterstd.h"

#include "ito33/cppunit.h"
#include "ito33/exception.h"
#include "ito33/date.h"
#include "ito33/dateutils.h"

#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/conversionschedule.h"
#include "ito33/finance/qualitycontrol.h"
#include "ito33/finance/modeloutput.h"
#include "ito33/finance/computationalflags.h"

#include "ito33/ihg/theoreticalmodel.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardrateflat.h"


#include "ihg/cbpathdepstructure.h"
#include "ihg/cbnumoutput.h"
#include "ihg/cbinstdata.h"
#include "ihg/model.h"

#include "ito33/pricing/model.h"
#include "ito33/pricing/pathdeppricer.h"
#include "ito33/pricing/pathdepevent.h"
#include "ito33/pricing/cbmeshmanager.h"
#include "ito33/pricing/cbparams.h"

#include "ito33/numeric/meshparams.h"
#include "ito33/numeric/numparams.h"
#include "ito33/numeric/extrapolationmode.h"
#include "ito33/numeric/interpolation.h"

#include "ihg/tests/testpathdep.h"
#include "utils.h"
#include "continuousevent.h"
#include "discreteevent.h"

using namespace ito33;
using namespace ito33::finance;
using namespace ito33::ihg;
using namespace ito33::numeric;
using namespace ito33::pricing;

void PathDepTest::TestPathDepNoEvents()
{
  size_t nPath     = 5;
  Date issueDate   = Date(2003, Date::Feb, 1);
  Date maturityDate= Date(2004, Date::May, 1);

  shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);
  shared_ptr<ConvertibleBond> pCB = InitCB(pSessionData,issueDate,maturityDate);
  
  double dVol    = 0.5;  
  double dLambda = .02;

  double dPriceRegularCB = Price(pCB,dVol,dLambda,issueDate);
   
//________________________ Path dependent part ______________________________//
  shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
  shared_ptr<HazardRate> pHR( new HazardRateFlat(dLambda) ); 
  ihg::Model model( pVol, pVol, pHR );
  shared_ptr<MeshParams> pMeshParams(new MeshParams);
  shared_ptr<QualityControl> pQualityControl(new QualityControl);

  shared_ptr<NumParams>
  pNumParams( new NumParams
                  (
                    *pQualityControl,
                    GetDoubleFrom( pCB->GetMaturityDate() ) 
                  - GetDoubleFrom( pSessionData->GetValuationDate() )
                  )
            );

  //Solve the different path dependent problem   
  std::vector<double> pdGridY;
  std::vector< pricing::CB > pPricingCB(nPath);
  std::vector< AutoPtr<pricing::CBLikeParams> > pCBPathDep(nPath);

  size_t nIdx;
  for ( nIdx = 0 ; nIdx < nPath ; nIdx++)
  { 
    pdGridY.push_back(nIdx);

    shared_ptr<ConvertibleBond>     
      pCb = InitCB(pSessionData,issueDate,maturityDate);

    pPricingCB[nIdx] = pricing::CB( *pCb );
    
    AutoPtr<pricing::CBLikeParams> pParamsCB(  
        new pricing::CBParams (pPricingCB[nIdx]) ) ;
    
    pParamsCB->SetNumParams(     pNumParams);
    pParamsCB->SetMeshParams(    pMeshParams);
    pParamsCB->SetYieldCurve(    pSessionData->GetYieldCurve() );
    pParamsCB->SetForeignCurve(  pSessionData->GetForeignCurve() );
    pParamsCB->SetDividends(     pSessionData->GetDividends() );
    pParamsCB->SetValuationTime( GetDoubleFrom(pSessionData->GetValuationDate()) );
    pParamsCB->SetSpotSharePrice(pSessionData->GetSpotSharePrice() );

    pCBPathDep[nIdx] = pParamsCB;
    
  } //end for loop
 
   
  //empty list of events
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
  ComputationalFlags flags;
  ihg::CBPathDepStructure cbPath(pdGridY,pCBPathDep,model,flags,pathDepEvents);
     
  size_t nNPathToSave = 0;
  cbPath.PrepareForTimestepping( nNPathToSave );

  PathDepPricer pathDepPricer;
  pathDepPricer.Price(cbPath, pathDepEvents);

  // Get and compare the path dep price
  AutoPtr<CBNumOutput> pNumOutput = cbPath.GetOutput();

  shared_ptr<ModelOutput> pOutput = pNumOutput->GetOutput();

  double dValue = pOutput->GetPrice();
 
  CPPUNIT_ASSERT_DOUBLES_EQUAL(dValue, dPriceRegularCB, 1.e-8);

}//PathDepTest::TestPathDepNoEvents()


void PathDepTest::TestPathDepConstraints()
{
 
  //--------------------------------------------------------------------------
  // Create the different windows of constraints
  //--------------------------------------------------------------------------
  Date issueDate     = Date(2003,Date::Feb, 1);
  Date maturityDate  = Date(2004,Date::May, 1);
  Date startCallDate = Date(2003,Date::Mar, 1);
  Date endCallDate   = Date(2003,Date::Jul, 1);
  Date startConvDate = Date(2003,Date::Apr, 1);
  Date endConvDate   = Date(2003,Date::Nov, 1);
  Date putDate       = Date(2003,Date::Jul, 1);

  //risk parameters
  double dVol    = 0.5;
  double dLambda = .02;

  std::vector<double> dPrice;

  //--------------------------------------------------------------------------
  //  0. no calls, no conversion, no puts     PDE0
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB0 
    = CreateCB(0, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB0, dVol, dLambda, issueDate )  );

  //--------------------------------------------------------------------------
  //  1. calls, no conversion, no puts        PDE1
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB1 
    = CreateCB(1, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB1, dVol, dLambda, issueDate ) );

  //--------------------------------------------------------------------------
  //  2. no calls, conversions, no puts       PDE2
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB2 
    = CreateCB(2, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB2, dVol, dLambda, issueDate ) );

  //--------------------------------------------------------------------------
  //  3. no calls, no conversion, put         PDE3
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB3 
    = CreateCB(3, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB3, dVol, dLambda, issueDate ) );

  //--------------------------------------------------------------------------
  //  4. calls, conversions, puts             PDE4
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB4 
    = CreateCB(4, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB4, dVol, dLambda, issueDate ) );

  //--------------------------------------------------------------------------
  //  5. calls notice, no conversion, no puts PDE5
  //--------------------------------------------------------------------------
  shared_ptr<finance::ConvertibleBond> pCB5 
    = CreateCB(5, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate,putDate);
  dPrice.push_back( Price( pCB5, dVol, dLambda, issueDate ) );

  //____________________Path Dependent Code____________________________________
  shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
  shared_ptr<HazardRate> pHR( new HazardRateFlat(dLambda) ); 
  ihg::Model model( pVol, pVol, pHR );

  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);
  shared_ptr<finance::QualityControl> m_pQualityControl(new QualityControl);

  shared_ptr<numeric::NumParams>
    pNumParams( new numeric::NumParams
                  (
                    *m_pQualityControl,
                    GetDoubleFrom( maturityDate ) 
                  - GetDoubleFrom( issueDate )
                  )
            );

  //create the y grid direction
  std::vector<double> pdGridY;
  size_t nIdx = 0;
  size_t nPath  = dPrice.size();
  std::vector< AutoPtr<pricing::CBLikeParams> > pCBPathDep(nPath);
  std::vector< pricing::CB > pPricingCB(nPath);

  shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);

  for ( nIdx = 0; nIdx < nPath ; nIdx++ )
  {
    pdGridY.push_back(nIdx);

    shared_ptr<finance::ConvertibleBond> 
      pCB = CreateCB(nIdx, issueDate, maturityDate, startCallDate, endCallDate,
      startConvDate, endConvDate, putDate);   

    pPricingCB[nIdx] = pricing::CB( *pCB );
    
    AutoPtr<pricing::CBLikeParams> pParamsCB( 
      new pricing::CBParams (pPricingCB[nIdx]) ) ;
      
    pParamsCB->SetNumParams(     pNumParams);
    pParamsCB->SetMeshParams(    pMeshParams);
    pParamsCB->SetYieldCurve(    pSessionData->GetYieldCurve() );
    pParamsCB->SetForeignCurve(  pSessionData->GetForeignCurve() );
    pParamsCB->SetDividends(     pSessionData->GetDividends() );
    pParamsCB->SetValuationTime(GetDoubleFrom(pSessionData->GetValuationDate()));
    pParamsCB->SetSpotSharePrice(pSessionData->GetSpotSharePrice() );

    pCBPathDep[nIdx] = pParamsCB;
  }

  //create an empty list of events
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
 
  // Create the path dep structure
  ComputationalFlags flag;
  ihg::CBPathDepStructure cbPath(pdGridY,pCBPathDep,model,flag,pathDepEvents);   

  //***************************************************************************
  // Compare results
  //***************************************************************************

  for (size_t nIdx = 0; nIdx < pdGridY.size(); nIdx++)
  {
    size_t nNPathToSave = nIdx;
    cbPath.PrepareForTimestepping( nNPathToSave );

    // Solve
    PathDepPricer pathDepPricer;
    pathDepPricer.Price(cbPath, pathDepEvents);

    AutoPtr<CBNumOutput> pNumOutput = cbPath.GetOutput();

    shared_ptr<ModelOutput> pOutput = pNumOutput->GetOutput();

    double dValue = pOutput->GetPrice();

    CPPUNIT_ASSERT_DOUBLES_EQUAL( dValue,dPrice[nIdx],1.e-8);

  } //end for loop


}//PathDepTest::TestPathDepConstraints()


void PathDepTest::TestContinuousEvent()
{

  //--------------------------------------------------------------------------
  // Create the different windows of constraints
  //--------------------------------------------------------------------------
  Date issueDate     = Date(2003,Date::Feb, 1);
  Date maturityDate  = Date(2004,Date::May, 1);
  Date startCallDate = Date(2003,Date::Mar, 1);
  Date endCallDate   = Date(2003,Date::Sep, 1);
  Date startConvDate = Date(2003,Date::Oct, 1);
  Date endConvDate   = Date(2003,Date::Dec, 1);

  //--------------------------------------------------------------------------
  // Risk parameters
  //--------------------------------------------------------------------------
  double dVol = 0.5;
  double dLambda = .02;
 
  shared_ptr<finance::ConvertibleBond> pCB 
    = CreateCB(6, issueDate, maturityDate, startCallDate, endCallDate,
              startConvDate, endConvDate, issueDate);

  double dPrice = Price( pCB, dVol, dLambda, issueDate );

  //______________path Dependent structure_____________________________________ 
  shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
  shared_ptr<HazardRate> pHR( new HazardRateFlat(dLambda) ); 
  ihg::Model model( pVol, pVol, pHR );


  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);
  shared_ptr<finance::QualityControl> m_pQualityControl(new QualityControl);

   shared_ptr<numeric::NumParams>
     pNumParams( new numeric::NumParams
                    (
                      *m_pQualityControl,
                      GetDoubleFrom( maturityDate ) 
                    - GetDoubleFrom( issueDate )
                    )
              );


  shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);

  size_t nPath = 2;
  size_t nIdx;

  std::vector<double> pdGridY;
  std::vector< AutoPtr<pricing::CBLikeParams> > pCBPathDep(nPath);
  std::vector< pricing::CB > pPricingCB(nPath);

  for ( nIdx = 0; nIdx < nPath ; nIdx++)
  {

    pdGridY.push_back(nIdx);
    
    //pde0 no call no conversion
    //pde1 call and conversion
    size_t nPDENum = 0;
    if ( nIdx != 0 )
     nPDENum = 6;

    shared_ptr<finance::ConvertibleBond> 
      pCB = CreateCB(nPDENum, issueDate, maturityDate, 
                     startCallDate, endCallDate,
                     startConvDate, endConvDate, issueDate);  

     pPricingCB[nIdx] = pricing::CB( *pCB );

     AutoPtr<pricing::CBLikeParams> 
       pParamsCB(  new pricing::CBParams (pPricingCB[nIdx]) ) ;
     
    pParamsCB->SetNumParams(     pNumParams);
    pParamsCB->SetMeshParams(    pMeshParams);
    pParamsCB->SetYieldCurve(    pSessionData->GetYieldCurve() );
    pParamsCB->SetForeignCurve(  pSessionData->GetForeignCurve() );
    pParamsCB->SetDividends(     pSessionData->GetDividends() );
    pParamsCB->SetValuationTime( GetDoubleFrom(pSessionData->GetValuationDate()) );
    pParamsCB->SetSpotSharePrice(pSessionData->GetSpotSharePrice() );

    pCBPathDep[nIdx] = pParamsCB;
    
  } //end for

   //create list of events
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;

  double startCallTime = GetDoubleFrom(startCallDate);
  double endCallTime   = GetDoubleFrom(endCallDate);

  shared_ptr<pricing::PathDepEvent> pStartCallEvent( 
    new ContinuousEvent(startCallTime, endCallTime, startCallTime) );
    
  shared_ptr<pricing::PathDepEvent> pEndCallEvent( 
    new ContinuousEvent(startCallTime, endCallTime, endCallTime) );

  pathDepEvents.push_back(pStartCallEvent);
  pathDepEvents.push_back(pEndCallEvent);

  double startConvTime = GetDoubleFrom(startConvDate);
  double endConvTime   = GetDoubleFrom(endConvDate);

  shared_ptr<pricing::PathDepEvent> pStartConvEvent( 
    new ContinuousEvent(startConvTime, endConvTime, startConvTime) );
   
  shared_ptr<pricing::PathDepEvent> pEndConvEvent( 
    new ContinuousEvent(startConvTime, endConvTime, endConvTime) );

  pathDepEvents.push_back(pStartConvEvent);
  pathDepEvents.push_back(pEndConvEvent); 
  
   // Create the path dep structure
  ComputationalFlags flag;
  ihg::CBPathDepStructure cbPath(pdGridY, pCBPathDep, model, flag, pathDepEvents);
     
  size_t nNPathToSave = 0;
  cbPath.PrepareForTimestepping( nNPathToSave );
  
  //the first path is set to be inactive
  cbPath.m_pbIsActive[1] = false;

  // Solve
  PathDepPricer pathDepPricer;
  pathDepPricer.Price(cbPath, pathDepEvents);

   //***************************************************************************
   // Compare results
   //***************************************************************************

   for (size_t nIdx = 0; nIdx < pdGridY.size(); nIdx++)
  {
    double dSpot = pSessionData->GetSpotSharePrice();
    double dValue = 0.0;

    if ( cbPath.m_pbIsActive[nIdx] == true ) 
    {
     Interpolate(
      cbPath.m_path[nIdx].meshes->GetS(), 
      cbPath.m_path[nIdx].instdata->m_pdPrices.Get(), 
      cbPath.m_path[nIdx].meshes->GetNbS(),
      &dSpot, 
      &dValue, 
      1, 
      numeric::ExtrapolationMode_Linear, 
      numeric::ExtrapolationMode_Linear);

     CPPUNIT_ASSERT_DOUBLES_EQUAL( dValue,dPrice,1.e-8);
    }

  } //end for loop

    
}//PathDepTest::TestContinuousEvent()
 

void PathDepTest::TestDiscreteEvent()
{
  Date issueDate     = Date(2003,Date::Feb, 1);
  Date maturityDate  = Date(2004,Date::May, 1);

  std::vector<Date> eventDate;

  eventDate.push_back( Date(2003,Date::Mar, 1) );
  eventDate.push_back( Date(2003,Date::Apr, 1) );
  eventDate.push_back( Date(2003,Date::Jun, 1) );
  eventDate.push_back( Date(2003,Date::Jul, 1) );
  eventDate.push_back( Date(2003,Date::Sep, 1) );

  //--------------------------------------------------------------------------
  // Risk parameters
  //--------------------------------------------------------------------------
  double dVol = 0.5;
  double dLambda = .02;
 
  shared_ptr<finance::ConvertibleBond> 
    pCB ( CreateCB(0, issueDate, maturityDate, issueDate, issueDate,
              issueDate, issueDate, issueDate) );

    
  shared_ptr<CallSchedule> pCall( new CallSchedule() );
  shared_ptr<ConversionSchedule> pConv( new ConversionSchedule() );
  shared_ptr<PutSchedule> pPut( new PutSchedule() );

  size_t nIdEvent = 0;
  size_t nEvent = eventDate.size();

  for ( nIdEvent ; nIdEvent < nEvent ; nIdEvent++ )
  {
    //call events
    shared_ptr<CallPeriod> pCallPeriod( CallPeriod::CreateWithStrike
                      (eventDate[nIdEvent], eventDate[nIdEvent], .9) );
    pCall->AddCallPeriod( pCallPeriod );

    // conv events
    shared_ptr<ConversionPeriod> 
      pConvPeriod(new ConversionPeriod(eventDate[nIdEvent], eventDate[nIdEvent], 1.0) );
    pConv->AddConversionPeriod( pConvPeriod );

    //put events
    pPut->AddPutWithStrike(eventDate[nIdEvent], .9);
  } //end loop

  pCB->SetCallSchedule(pCall);
  pCB->SetConversionSchedule(pConv);
  pCB->SetPutSchedule(pPut);

  double dPrice = Price( pCB, dVol, dLambda, issueDate );
  
  //_______________path Dependent structure____________________________________
  shared_ptr<Volatility> pVol( new VolatilityFlat(dVol) );
  shared_ptr<HazardRate> pHR( new HazardRateFlat(dLambda) );
  ihg::Model model( pVol, pVol, pHR );

  shared_ptr<numeric::MeshParams> pMeshParams(new numeric::MeshParams);
  shared_ptr<finance::QualityControl> m_pQualityControl(new QualityControl);

  shared_ptr<numeric::NumParams>
     pNumParams( new numeric::NumParams
                    (
                      *m_pQualityControl,
                      GetDoubleFrom( maturityDate ) 
                    - GetDoubleFrom( issueDate )
                    )
              );

  shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);

  size_t nPath = 2;
  size_t nIdx;

  std::vector< AutoPtr<pricing::CBLikeParams> > pCBPathDep(nPath);
  std::vector< pricing::CB > pPricingCB(nPath);
  std::vector<double> pdGridY;

  for ( nIdx = 0 ; nIdx < nPath ; nIdx++)
  {
    pdGridY.push_back(nIdx);

    shared_ptr<finance::ConvertibleBond> 
      pCB =  CreateCB(0, issueDate, maturityDate, issueDate, issueDate,
                     issueDate, issueDate, issueDate);
  
    
    if ( nIdx == 1 )
    {
     pCB->SetCallSchedule(pCall);
     pCB->SetConversionSchedule(pConv);
     pCB->SetPutSchedule(pPut);
    }
   

    pPricingCB[nIdx] = pricing::CB( *pCB );
    
    AutoPtr<pricing::CBLikeParams> 
      pParamsCB(  new pricing::CBParams (pPricingCB[nIdx]) ) ;
     
    pParamsCB->SetNumParams(     pNumParams);
    pParamsCB->SetMeshParams(    pMeshParams);
    pParamsCB->SetYieldCurve(    pSessionData->GetYieldCurve() );
    pParamsCB->SetForeignCurve(  pSessionData->GetForeignCurve() );
    pParamsCB->SetDividends(     pSessionData->GetDividends() );
    pParamsCB->SetValuationTime( GetDoubleFrom(pSessionData->GetValuationDate()) );
    pParamsCB->SetSpotSharePrice(pSessionData->GetSpotSharePrice() );

    pCBPathDep[nIdx] = pParamsCB;

  }

  //Create the list of events
  std::list< shared_ptr<pricing::PathDepEvent> > pathDepEvents;
 
  size_t nEvents = eventDate.size();

  for ( nIdx = 0 ; nIdx < nEvents ; nIdx++ )
  {
    double eventTime = GetDoubleFrom(eventDate[nIdx]);

    shared_ptr<pricing::PathDepEvent> 
      pCallEvent( new DiscreteEvent(eventTime,eventTime,eventTime) );
  
    shared_ptr<pricing::PathDepEvent> 
      pConvEvent( new DiscreteEvent(eventTime,eventTime,eventTime) );

    shared_ptr<pricing::PathDepEvent> 
      pPutEvent( new DiscreteEvent(eventTime,eventTime,eventTime) );

    pathDepEvents.push_back(pCallEvent);
    pathDepEvents.push_back(pConvEvent);
    pathDepEvents.push_back(pPutEvent);
  }

   // Create the path dep structure
  ComputationalFlags flag;
  ihg::CBPathDepStructure cbPath(pdGridY,pCBPathDep,model,flag,pathDepEvents);
 
  size_t nNPathToSave = 0;
  cbPath.PrepareForTimestepping( nNPathToSave );

   // Solve
   PathDepPricer pathDepPricer;
   pathDepPricer.Price(cbPath, pathDepEvents);

   //***************************************************************************
   // Compare results
   //***************************************************************************

   for (size_t nIdx = 0; nIdx < pdGridY.size(); nIdx++)
  {
    double dSpot = pSessionData->GetSpotSharePrice();
    double dValue = 0.0;

    if ( cbPath.m_pbIsActive[nIdx] == true ) 
    {
     Interpolate(
      cbPath.m_path[nIdx].meshes->GetS(), 
      cbPath.m_path[nIdx].instdata->m_pdPrices.Get(), 
      cbPath.m_path[nIdx].meshes->GetNbS(),
      &dSpot, 
      &dValue, 
      1, 
      numeric::ExtrapolationMode_Linear, 
      numeric::ExtrapolationMode_Linear);

     CPPUNIT_ASSERT_DOUBLES_EQUAL( dValue,dPrice,1.e-8);

    }

  } //end for loop

}//PathDepTest::TestDiscreteEvent()


//--------------------------------------------------------------------------
// price CB using the classic method
//--------------------------------------------------------------------------
double PathDepTest::Price(shared_ptr<finance::ConvertibleBond> &pCB,
                          double dVol,double dLambda,Date issueDate)
{
 
   shared_ptr<ihg::TheoreticalModel> pModel(new ihg::TheoreticalModel);
  
   shared_ptr<ihg::Volatility> pVol( new ihg::VolatilityFlat(dVol) );
   pModel->SetVolatility( pVol );
   shared_ptr<ihg::HazardRate> pHR( new ihg::HazardRateFlat(dLambda) ); 
   pModel->SetHazardRate( pHR );

   shared_ptr<ComputationalFlags> flags(new ComputationalFlags);
   flags->SetAnalysisDate( issueDate );

   pCB->SetComputationalFlags(flags);

   shared_ptr<finance::ModelOutput> output = pModel->Compute(*pCB);

   return ( output->GetPrice() );
}// PathDepTest::Price


//--------------------------------------------------------------------------
// Create CB with some constraints
//--------------------------------------------------------------------------
shared_ptr<ito33::finance::ConvertibleBond> 
 PathDepTest::CreateCB(size_t nPDENum,
             ito33::Date issueDate, ito33::Date maturityDate,
             ito33::Date startCallDate, ito33::Date endCallDate,
             ito33::Date startConvDate, ito33::Date endConvDate,
             ito33::Date putDate)
{

   shared_ptr<SessionData> pSessionData = InitSessionData(issueDate);
   shared_ptr<finance::ConvertibleBond> 
     pCB = InitCB(pSessionData,issueDate,maturityDate);

    shared_ptr<CallSchedule>       pCall( new CallSchedule() );
    shared_ptr<PutSchedule>        pPut( new PutSchedule() );
    shared_ptr<ConversionSchedule> pConv( new ConversionSchedule() );
   
    SetConstraints(nPDENum,pCall,pPut,pConv,startCallDate,
      endCallDate,startConvDate,endConvDate,putDate);

    pCB->SetCallSchedule(pCall);
    pCB->SetPutSchedule(pPut);
    pCB->SetConversionSchedule(pConv);

   return pCB;
}


//--------------------------------------------------------------------------
// Set the constraints
//--------------------------------------------------------------------------
void PathDepTest::SetConstraints(size_t nPDENum,
    ito33::shared_ptr<ito33::finance::CallSchedule> &pCall,
    ito33::shared_ptr<ito33::finance::PutSchedule>  &pPut,
    ito33::shared_ptr<ito33::finance::ConversionSchedule> &pConv,
    ito33::Date startCallDate, ito33::Date endCallDate,
    ito33::Date startConvDate, ito33::Date endConvDate,
    ito33::Date putDate)
{
  Date earlyDateStart(1900, Date::Jan, 1);
  Date earlyDateEnd(1900, Date::Jan, 2);

   if ( nPDENum  == 0 )
   {
     //  0. no calls, no conversion, no puts     PDE0
     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(earlyDateStart,earlyDateEnd,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );
   }
   else if ( nPDENum == 1)
   {
     shared_ptr<CallPeriod> pCallPeriod( CallPeriod::CreateWithStrike
                           (startCallDate, endCallDate,1.0) );
     pCall->AddCallPeriod( pCallPeriod );

     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(earlyDateStart,earlyDateEnd,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );
   }
   else if ( nPDENum == 2 )
   {
     //  2. no calls, conversions, no puts       PDE2
     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(startConvDate,endConvDate,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );
       
   }
   else if ( nPDENum == 3 )
   {
     //  3. no calls, no conversion, put         PDE3
     pPut->AddPutWithStrike(putDate,1.0);

     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(earlyDateStart,earlyDateEnd,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );
   }
   else if ( nPDENum == 4 )
   {
    //  4. calls, conversions, puts             PDE4
     shared_ptr<CallPeriod> pCallPeriod( CallPeriod::CreateWithStrike
                                       (startCallDate, endCallDate, 1.0) ); 
     pCall->AddCallPeriod( pCallPeriod );

     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(startConvDate,endConvDate,1.0) ); 
     pConv->AddConversionPeriod( pConvPeriod );
       

     pPut->AddPutWithStrike(putDate, 1.0);
   }
   else if ( nPDENum == 5 )
   {
    //  5. calls notice, no conversion, no puts PDE5
     pCall->SetNoticePeriod(20);
     shared_ptr<CallPeriod> 
       pCallPeriod( CallPeriod::CreateWithStrike
                                       (startCallDate, endCallDate, 1.0)  );
     pCall->AddCallPeriod( pCallPeriod );

     shared_ptr<ConversionPeriod> 
       pConvPeriod( new ConversionPeriod(earlyDateStart,earlyDateEnd,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );
   }
   else if ( nPDENum == 6 )
   {
     //call and conversion period
     shared_ptr<CallPeriod>
       pCallPeriod( CallPeriod::CreateWithStrike
                                       (startCallDate, endCallDate, 1.0) );
     pCall->AddCallPeriod( pCallPeriod );

     shared_ptr<ConversionPeriod> 
       pConvPeriod( 
         new ConversionPeriod(startConvDate,endConvDate,1.0) );
     pConv->AddConversionPeriod( pConvPeriod );

   }
   else
   {
    FAIL("Should not reached this statement.");
   }

}//PathDepTest::SetConstraints
