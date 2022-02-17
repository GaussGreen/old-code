/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_hrwithspotcomponentpower.cpp
// Purpose:     implement for HRSpotComponentPower class
// Author:      Ito33
// Created:     2004/11/23
// RCS-ID:      $Id: parametrization_hrwithspotcomponentpower.cpp,v 1.12 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2003 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include <cmath>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/termstructurecds.h"
#include "ito33/finance/sessiondata.h"

#include "ihg/hrspotcomponentpowercalibrator.h"

#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/hazardratewithtimecomponent.h"
#include "ito33/ihg/parametrization_hrwithspotcomponentpower.h"
#include "ito33/ihg/parametrization_hrwithtimecomponent.h"
#include "ito33/ihg/volatility.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"

#include "ito33/xml/write.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/common.h"

extern const ito33::Error
  ITO33_BAD_PARAM;
extern const ito33::finance::Error
  ITO33_CALIBRATION_FAIL,
  ITO33_EMPTY_CDSCURVE;

extern const ito33::ihg::Error
  ITO33_IHG_CALIBRATION_DEGENERATE;

namespace ito33
{

namespace ihg
{


ParametrizationHRWithSpotComponentPower::ParametrizationHRWithSpotComponentPower
(shared_ptr<Volatility> pVolatility)
{  
  if ( !pVolatility )
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Calibration on hazard rate having a spot component needs"
                  " that volatility is set.")
          );

  m_pVolatility = pVolatility;
}

shared_ptr<HazardRateCombo> 
ParametrizationHRWithSpotComponentPower::CalibrateWithCDSs
(const finance::TermStructureCDS& tsCDS)
{
  CHECK_COND(!tsCDS.GetAll().empty(), ITO33_EMPTY_CDSCURVE);

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

  // Init output objects to default values in case calibration fails
  double dS0 = tsCDS.GetAll().front()->GetSessionData()->GetSpotSharePrice();
  m_pHRSpotComponentPower = make_ptr( new HRSpotComponentPower(0.0, dS0) );

  HazardRateSpotComponentPowerCalibrator calibrator;

  try
  {
    // check that the cds term structure contains at least two elements
    CheckSize(tsCDS,2);
    
    shared_ptr<HazardRateCombo>
      pHRCombo = calibrator.Calibrate(tsCDS, m_pVolatility);

    m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower> 
                              ( pHRCombo->GetSpotComponent() );

    return pHRCombo;
  }  
  catch(ito33::numeric::Exception)
  {
    // Save the results, even if it did not converge
    shared_ptr<HazardRateCombo>
      pHRCombo = calibrator.GetHazardRate();
    m_pHRSpotComponentPower = static_pointer_cast<HRSpotComponentPower>
                              ( pHRCombo->GetSpotComponent() );

    // If beta is close to zero, the code likely failed because the data 
    // is degenerate, and a time-only parametrization should be used
    if ( fabs(m_pHRSpotComponentPower->GetBeta()) < 1.e-3)
      throw EXCEPTION(ITO33_IHG_CALIBRATION_DEGENERATE);
    else
      throw EXCEPTION(ITO33_CALIBRATION_FAIL);

  }
}

void 
ParametrizationHRWithSpotComponentPower::Dump(ito33::XML::Tag& tagParent) const
{
  {
  // create the root tag
  ito33::XML::Tag 
    tagParametrization(XML_TAG_PARAMETRIZATION_HRWITHSPOTCOMPONENTPOWER_ROOT,
                       tagParent);
  
  // dump the vol (set in the constructor, so must be valid)
  m_pVolatility->Dump(tagParametrization);

  }
}

void ParametrizationHRWithSpotComponentPower::Visit
     (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationHRWithSpotComponentPower(*this);
}


} // namespace ihg

} // namespace ito33
