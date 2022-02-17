/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volwithtimecomponent.cpp
// Purpose:     implement for ParametrizationVolWithTimeComponent class
// Author:      Wang
// Created:     03/11/04
// RCS-ID:      $Id: parametrization_volwithtimecomponent.cpp,v 1.8 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/sessiondata.h"

#include "ihg/voltimecomponentcalibrator.h"

#include "ito33/ihg/version.h"
#include "ito33/ihg/hazardrate.h"
#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/volatilitycombo.h"

#include "ito33/ihg/parametrization_volwithtimecomponent.h"

#include "ito33/xml/write.h"

#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/spotcomponent.h"

#include "ihg/xml/parametrization_visitor.h"

extern const ito33::Error
    ITO33_BAD_PARAM,
    ITO33_NULL_PARAM;

extern const ito33::finance::Error
    ITO33_EMPTY_OPTIONCURVE,
    ITO33_CALIBRATION_FAIL;

namespace ito33
{

namespace ihg
{


ParametrizationVolWithTimeComponent::ParametrizationVolWithTimeComponent
    (const shared_ptr<HazardRate>& pHR)
{
  if ( !pHR )
    throw EXCEPTION_MSG
          (
            ITO33_NULL_PARAM,
            TRANS("Setting null hazard rate")
          );

  m_pHazardRate = pHR;
}

void ParametrizationVolWithTimeComponent::SetSpotComponent
     ( const shared_ptr<SpotComponent>& pSpotComponent )
{
  m_pSpotComponent = pSpotComponent;

  if ( !m_pSpotComponent )
    throw EXCEPTION_MSG
          (
            ITO33_NULL_PARAM,
            TRANS("Setting invalid spot component for volatility")
          );
}

shared_ptr<VolatilityWithTimeComponent>
ParametrizationVolWithTimeComponent::CalibrateWithOptions
(const finance::TermStructureOption& tsOption)
{ 
  CHECK_COND(!tsOption.GetAll().empty(), ITO33_EMPTY_OPTIONCURVE);

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

  shared_ptr<VolatilityWithTimeComponent> pVol;

  if ( m_pSpotComponent ) // volatility is not time only
  {
    pVol = make_ptr( new VolatilityCombo(m_pSpotComponent) ); 
  } 
  
  shared_ptr<VolatilityWithTimeComponent> pVolTmp;

  try
  {
    VolTimeComponentCalibrator calibrator(m_pHazardRate);

    pVolTmp = calibrator.Calibrate( tsOption, pVol );
  }
  catch(ito33::numeric::Exception)
  {
    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  return pVolTmp;
}

void ParametrizationVolWithTimeComponent::Dump(ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLWITHTIMECOMPONENT_ROOT, 
                        tagParent);

  m_pHazardRate->Dump(tagParametrization);

  if (m_pSpotComponent)
    tagParametrization.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);
}

void 
ParametrizationVolWithTimeComponent::Visit(ParametrizationVisitor& visit) const
{
  visit.OnParametrizationVolWithTimeComponent(*this);
}


} // namespace ihg

} // namespace ito33
