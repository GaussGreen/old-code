/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volflat_hrwithtimecomponent.cpp
// Purpose:     implement for ParametrizationVolFlatHRWithTimeComponent class
// Created:     03/11/04
// RCS-ID:      $Id: parametrization_volflat_hrwithtimecomponent.cpp,v 1.39 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/option.h"
#include "ito33/finance/cdslike.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/termstructureparbond.h"
#include "ito33/finance/termstructureoption.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/termstructurederivative.h"

#include "ito33/finance/basket_goodtype.h"

#include "ito33/xml/finance/termstructure.h"
#include "ito33/xml/write.h"

#include "ihg/volflathrtimecomponentcalibrator.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_volflat_hrwithtimecomponent.h"

#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/spotcomponent.h"

extern const ito33::Error
  ITO33_BAD_PARAM,
  ITO33_NULL_PARAM;

extern const ito33::finance::Error
  ITO33_CALIBRATION_FAIL,
  ITO33_EMPTY_DERIVATIVECURVE;

namespace ito33
{

namespace ihg
{


void ParametrizationVolFlatHRWithTimeComponent::SetSpotComponent
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

shared_ptr<finance::TheoreticalModel> 
  ParametrizationVolFlatHRWithTimeComponent::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pVolatility, m_pHazardRate));

  return pTM;
}

void ParametrizationVolFlatHRWithTimeComponent::CalibrateWithOptionAndParBonds
      (const finance::Option& option,
       const finance::TermStructureParBond& tsParBond)
{
  Calibrate(option, tsParBond);
}


void ParametrizationVolFlatHRWithTimeComponent::CalibrateWithOptionAndCDSs
     (const finance::Option& option, const finance::TermStructureCDS& tsCDS)
{
  Calibrate(option, tsCDS); 
}

void 
ParametrizationVolFlatHRWithTimeComponent::CalibrateWithOptionAndOptions
     (const finance::Option& option, 
      const finance::TermStructureOption& tsOption)
{
  Calibrate(option, tsOption);
}

void ParametrizationVolFlatHRWithTimeComponent::Calibrate(
    const finance::Derivative& deriv, 
    const finance::TermStructureDerivative& tsDerivs)
{
    
  CHECK_COND( tsDerivs.Validate(), ITO33_BAD_PARAM);
  
  CHECK_COND(!tsDerivs.GetAll().empty(), ITO33_EMPTY_DERIVATIVECURVE);

   std::string sXmlTagName;
  if ( dynamic_pointer_cast<finance::ParBond>( tsDerivs.GetAll().front() ) )
    sXmlTagName = XML_TAG_TERMSTRUCTURE_PARBOND;
  else if ( dynamic_pointer_cast<finance::Option>( tsDerivs.GetAll().front() ) )
    sXmlTagName = XML_TAG_TERMSTRUCTURE_OPTION;
  else if ( dynamic_pointer_cast<finance::CDSLike>( tsDerivs.GetAll().front() ) )
    sXmlTagName = XML_TAG_TERMSTRUCTURE_CDS;
  else
     throw EXCEPTION_MSG
      (
        ITO33_BAD_PARAM,
        TRANS("VolFlatHRTimeOnly: The derivatives passed in are not supported.")
       );


  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0( deriv.GetSessionData()->GetSpotSharePrice() );

  // Write the calibration parameters to XML file when debug requested
  if ( IsDebugOutputEnabled() )
  {
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      tsDerivs.GetAll().front()->GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      DumpTermStructureDerivative(tsDerivs, tagParam, sXmlTagName);
      deriv.Dump(tagParam);
    }
  }

  shared_ptr<HazardRateWithTimeComponent> pHazardRate;


  if ( !m_pSpotComponent ) // hazard rate is time only
  {    
    // unlike with cds and zero coupons, even if the hazard rate is time only
    // we still can't do decoupled calibration, so we just create an empty
    // hazard rate, the calibrator checks on the passed pointer to decide
    // on the algorithm
    if ( dynamic_pointer_cast<finance::Option>( tsDerivs.GetAll().front() ) )
      pHazardRate = make_ptr( new HazardRateTimeOnly() );
  }
  else
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 
  
  VolFlatHRTimeComponentCalibrator calibrator;

  try
  {
    calibrator.Calibrate(deriv, tsDerivs, pHazardRate);
  }
  catch (const ito33::numeric::Exception& )
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

void ParametrizationVolFlatHRWithTimeComponent::Calibrate
     (const finance::BasketGoodType& basket)
{
  shared_ptr<finance::Derivative> pDerivative = basket.GetOption();
  if ( !pDerivative )
    pDerivative = basket.GetBondLike();

  ASSERT(pDerivative);

  ASSERT( !basket.GetTermStructure().GetAll().empty() );

  Calibrate(*pDerivative, basket.GetTermStructure());
}

void ParametrizationVolFlatHRWithTimeComponent::Dump
    (ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_VOLFLATHRWITHTIMECOMPONENT_ROOT,
                        tagParent);
  
  // only output data if it is set
  if (m_pSpotComponent)
    tagParametrization.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);
}

void ParametrizationVolFlatHRWithTimeComponent::Visit
    (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolFlatHRWithTimeComponent(*this);
}


} // namespace ihg

} // namespace ito33
