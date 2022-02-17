/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volpower_hrwithtimecomponent.cpp
// Purpose:     implement ParametrizationVolPowerHRWithTimeComponent class
// Author:      ITO 33
// Created:     2005/01/05
// RCS-ID:      $Id: parametrization_volpower_hrwithtimecomponent.cpp,v 1.25 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/basket_goodtype.h"

#include "ihg/volpowerhrtimecomponentcalibrator.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_volpower_hrwithtimecomponent.h"

#include "ito33/xml/write.h"
#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/spotcomponent.h"

extern const ito33::Error
  ITO33_BAD_PARAM;

extern const ito33::finance::Error
  ITO33_EMPTY_CDSCURVE, 
  ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{

void ParametrizationVolPowerHRWithTimeComponent::SetSpotComponent
     ( const shared_ptr<SpotComponent>& pSpotComponent )
{
  if ( !pSpotComponent )
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Setting invalid spot component for hazard rate.")
          );

  m_pSpotComponent = pSpotComponent;
}

void ParametrizationVolPowerHRWithTimeComponent::CalibrateWithOptionsAndCDSs
     (const finance::Option& option1, 
      const finance::Option& option2,
      const finance::TermStructureCDS& tsCDS)
{
  CHECK_COND(!tsCDS.GetAll().empty(), ITO33_EMPTY_CDSCURVE);

  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0( option1.GetSessionData()->GetSpotSharePrice() );

  // Write the calibration parameters to XML file when debug requested
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
      option1.Dump(tagParam);
      option2.Dump(tagParam);
      tsCDS.Dump(tagParam);      
    }
  }

  shared_ptr<HazardRateWithTimeComponent> pHazardRate;

  if ( m_pSpotComponent ) // hazard rate is not time only
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 

  VolPowerHRTimeComponentCalibrator calibrator;
  
  try
  {
    // The first option passed to the calibrator is used to obtain a
    // guess for volatility. The guess will be better if the option is
    // close to at-the-money
    double dStrike1 = option1.GetStrike();
    double dStrike2 = option2.GetStrike();
    double dS0 = option1.GetSessionData()->GetSpotSharePrice();
    if ( fabs(dStrike1-dS0) < fabs(dStrike2-dS0) )
      calibrator.Calibrate(option1, option2, tsCDS, pHazardRate);
    else
      calibrator.Calibrate(option2, option1, tsCDS, pHazardRate);      
    
  }
  catch(ito33::numeric::Exception)
  {
    // Even if the calibration fails, keep the last calibrated values
    m_pVolatility = calibrator.GetVolatility();

    m_pHazardRate = calibrator.GetHazardRate();

    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  // Save the calibrated values
  m_pVolatility = calibrator.GetVolatility();

  m_pHazardRate = calibrator.GetHazardRate();
}

void ParametrizationVolPowerHRWithTimeComponent::Calibrate
     (const finance::BasketGoodType& basket)
{
  CalibrateWithOptionsAndCDSs(*basket.GetOption(),
                              *basket.GetOption2(),
                              *basket.GetTermStructureCDS());
}

shared_ptr<finance::TheoreticalModel> 
  ParametrizationVolPowerHRWithTimeComponent::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pVolatility, m_pHazardRate));

  return pTM;
}

void ParametrizationVolPowerHRWithTimeComponent::Dump
    (ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLPOWERHRWITHTIMECOMPONENT_ROOT,
                       tagParent);
  
  // only output data if it is set
  if (m_pSpotComponent)
    tagParametrization.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);
}

void ParametrizationVolPowerHRWithTimeComponent::Visit
    (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolPowerHRWithTimeComponent(*this);
}

} // namespace ihg

} // namespace ito33
