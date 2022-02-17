/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_hrwithtimecomponent.cpp
// Purpose:     implement for ParametrizationHRWithTimeComponent class
// Created:     03/11/04
// RCS-ID:      $Id: parametrization_hrwithtimecomponent.cpp,v 1.26 2006/08/22 12:29:04 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/termstructurederivative.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureeds.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/sessiondata.h"

#include "ihg/hrtimecomponentcalibrator.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/version.h"

#include "ito33/ihg/parametrization_hrwithtimecomponent.h"

#include "ito33/xml/write.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/spotcomponent.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/common.h"

extern const ito33::Error
  ITO33_BAD_PARAM,
  ITO33_BAD_DATA;


extern const ito33::finance::Error
  ITO33_EMPTY_PARBONDCURVE,
  ITO33_EMPTY_CDSCURVE,  
  ITO33_EMPTY_EDSCURVE,  
  ITO33_EMPTY_OPTIONCURVE,
  ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{

void ParametrizationHRWithTimeComponent::CheckVolatility()
{
  if (!m_pVolatility)
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Calibration on hazard rate having a spot component needs"
                  " that volatility is set.")
          );
}

void ParametrizationHRWithTimeComponent::SetSpotComponent
     ( const shared_ptr<SpotComponent>& pSpotComponent )
{
  m_pSpotComponent = pSpotComponent;

  if ( !m_pSpotComponent )
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Setting invalid spot component for hazard rate")
          );
}

shared_ptr<HazardRateWithTimeComponent>
ParametrizationHRWithTimeComponent::CalibrateWithTermStructurePossibleTimeOnly
(const finance::TermStructureDerivative& tsDeriv)
{
  shared_ptr<HazardRateWithTimeComponent> pHazardRate;
  shared_ptr<Volatility> pVolatility;

  if ( !m_pSpotComponent ) // hazard rate is time only
  {
    // generate a temporary volatility. The value does not matter.
    pVolatility = make_ptr( new VolatilityFlat(0.2) );
  }
  else
  {
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 

    CheckVolatility();

    pVolatility = m_pVolatility;
  } 
  
  shared_ptr<HazardRateWithTimeComponent> pHazardRateTmp;

  try
  {
    HazardRateTimeComponentCalibrator calibrator;

    pHazardRateTmp = calibrator.Calibrate(tsDeriv, pVolatility, pHazardRate);
  }
  catch(ito33::numeric::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return pHazardRateTmp;
}

shared_ptr<HazardRateWithTimeComponent>
ParametrizationHRWithTimeComponent::CalibrateWithParBonds
(const finance::TermStructureParBond& tsParBond)
{
  CHECK_COND(!tsParBond.GetAll().empty(), ITO33_EMPTY_PARBONDCURVE);
  
  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0
      ( tsParBond.GetAll().front()->GetSessionData()->GetSpotSharePrice() );

  if ( IsDebugOutputEnabled() )
  {    
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {      
      tsParBond.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      tsParBond.Dump(tagParam);
    }
  }

  return CalibrateWithTermStructurePossibleTimeOnly(tsParBond);
}

shared_ptr<HazardRateWithTimeComponent>
ParametrizationHRWithTimeComponent::CalibrateWithCDSs
(const finance::TermStructureCDS& tsCDS)
{
  CHECK_COND(!tsCDS.GetAll().empty(), ITO33_EMPTY_CDSCURVE);
  
  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0
      ( tsCDS.GetAll().front()->GetSessionData()->GetSpotSharePrice() );

  if ( IsDebugOutputEnabled() )
  {    
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      tsCDS.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      tsCDS.Dump(tagParam);
    }
  }

  return CalibrateWithTermStructurePossibleTimeOnly(tsCDS);
}


shared_ptr<HazardRateWithTimeComponent>
ParametrizationHRWithTimeComponent::CalibrateWithEDSs
(const finance::TermStructureEDS& tsEDS)
{
  CHECK_COND(!tsEDS.GetAll().empty(), ITO33_EMPTY_EDSCURVE);

  // An eds always needs a volatility
  if (!m_pVolatility)
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Calibration to an EDS term structure requires that"
                  " the volatility be set.")
          );

  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0
      ( tsEDS.GetAll().front()->GetSessionData()->GetSpotSharePrice() );

  if ( IsDebugOutputEnabled() )
  {    
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      tsEDS.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      tsEDS.Dump(tagParam);
    }
  }


  shared_ptr<HazardRateWithTimeComponent> pHazardRate;

  if ( m_pSpotComponent ) // hazard rate is not time only
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 
  
  shared_ptr<HazardRateWithTimeComponent> pHazardRateTmp;

  try
  {
    HazardRateTimeComponentCalibrator calibrator;

    pHazardRateTmp = calibrator.Calibrate(tsEDS, m_pVolatility, pHazardRate);
  }
  catch(ito33::numeric::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return pHazardRateTmp;
}

shared_ptr<HazardRateWithTimeComponent>
ParametrizationHRWithTimeComponent::CalibrateWithOptions
(const finance::TermStructureOption& tsOption)
{ 
  CHECK_COND(!tsOption.GetAll().empty(), ITO33_EMPTY_OPTIONCURVE);

  // Calibration with option always needs a volatility
  if (!m_pVolatility)
    throw EXCEPTION_MSG
          (
            ITO33_BAD_DATA,
            TRANS("Calibration to an option term structure requires that"
                  " the volatility be set.")
          );

  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0
      ( tsOption.GetAll().front()->GetSessionData()->GetSpotSharePrice() );

  if ( IsDebugOutputEnabled() )
  {    
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      tsOption.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      tsOption.Dump(tagParam);
    }
  }

  shared_ptr<HazardRateWithTimeComponent> pHazardRate;

  if ( m_pSpotComponent ) // hazard rate is not time only
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 

  shared_ptr<HazardRateWithTimeComponent> pHazardRateTmp;

  try
  {
    HazardRateTimeComponentCalibrator calibrator;

    pHazardRateTmp = calibrator.Calibrate
                                ( tsOption, m_pVolatility, pHazardRate );
  }
  catch(ito33::numeric::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return pHazardRateTmp;
}

void ParametrizationHRWithTimeComponent::Dump(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_HRWITHTIMECOMPONENT_ROOT, 
                        tagParent);
  
  // only output data if it is set
  if (m_pVolatility)
    m_pVolatility->Dump(tagParametrization);

  if (m_pSpotComponent)
    tagParametrization.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);
}

void 
ParametrizationHRWithTimeComponent::Visit(ParametrizationVisitor& visit) const
{
  visit.OnParametrizationHRWithTimeComponent(*this);
}

} // namespace ihg

} // namespace ito33
