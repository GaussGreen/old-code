/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_volflat_hrpower.cpp
// Purpose:     implement ParametrizationVolFlatHRower class
// Author:      ITO 33
// Created:     2005/07/29
// RCS-ID:      $Id: parametrization_volflat_hrpower.cpp,v 1.15 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/beforestd.h"
#include <fstream>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"
#include "ito33/numeric/exception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/derivative.h"
#include "ito33/finance/derivatives.h"
#include "ito33/finance/sessiondata.h"
#include "ito33/finance/basket_goodtype.h"

#include "ihg/calibratorgeneral.h"
#include "ihg/translator.h"

#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"
#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_volflat_hrpower.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"

extern const ito33::finance::Error  
  ITO33_CALIBRATION_FAIL,
  ITO33_EMPTY_DERIVATIVELIST;


namespace ito33
{

namespace ihg
{


void ParametrizationVolFlatHRPower::CalibrateWithDerivatives
     (const finance::Derivatives& derivatives) 
{

  finance::Derivatives::Elements elements = derivatives.GetAll();
  size_t nNbDerivatives = elements.size();

  CHECK_COND(nNbDerivatives > 0, ITO33_EMPTY_DERIVATIVELIST);

  // Write the calibration parameters to XML file when debug requested
  if ( IsDebugOutputEnabled() )
  {
    std::ofstream ofs(GetDebugOutputFile().c_str());
    ito33::XML::RootTag tagRoot(XML_TAG_IHG_ROOT, ofs);
    tagRoot.precision(10);
    tagRoot.Attr(XML_ATTR_IHG_ROOT_VERSION, ITO33_IHG_VERSION_DOT_STRING);

    {
      derivatives.GetSessionData()->Dump(tagRoot);
    }    

    {
      ito33::XML::Tag tagParam(XML_TAG_IHG_CALIBRATION, tagRoot);
      Dump(tagParam);
      // dump the derivatives to be calibrated
      {
        ito33::XML::Tag tagDerivatives(XML_TAG_DERIVATIVES, tagParam);
      
        finance::Derivatives::Elements::const_iterator iter;
      
        for (iter = elements.begin(); iter != elements.end(); ++iter)
          iter->first->Dump(tagDerivatives);
      }
    }
  } // if debug output

  CalibratorGeneral calibrator(*m_pProcess);

  shared_ptr<finance::SessionData> 
    pSessionData = derivatives.GetSessionData();

  shared_ptr<ihg::Translator> pTranslator(new ihg::Translator(VolType_flat, 
    HRType_power, pSessionData->GetSpotSharePrice()));

  shared_ptr<ihg::TheoreticalModel> pModel;
  
  try
  {
    pModel = calibrator.Calibrate(pTranslator.get(), derivatives);   
  }
  catch(ito33::numeric::Exception)
  {
    // Even if the calibration fails, keep the last calibrated values
    pModel = calibrator.GetLastCalibratedProcess();

    m_pVolatility = static_pointer_cast<VolatilityFlat>
                    ( pModel->GetVolatility() ); 
    m_pHazardRate = static_pointer_cast<HazardRatePower>
                    ( pModel->GetHazardRate() ); 

    throw EXCEPTION(ITO33_CALIBRATION_FAIL);
  }

  // Save the calibrated values
  m_pVolatility = static_pointer_cast<VolatilityFlat>
                  ( pModel->GetVolatility() ); 
  m_pHazardRate = static_pointer_cast<HazardRatePower>
                  ( pModel->GetHazardRate() ); 

}

void ParametrizationVolFlatHRPower::Calibrate
     (const finance::BasketGoodType& basket)
{
  CalibrateWithDerivatives( *basket.GetDerivatives() );
}

void ParametrizationVolFlatHRPower::Dump
    (ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag tagParametrization
    (XML_TAG_PARAMETRIZATION_VOLFLATHRPOWER_ROOT, tagParent);
}

void ParametrizationVolFlatHRPower::Visit
    (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolFlatHRPower(*this);
}

shared_ptr<finance::TheoreticalModel> 
  ParametrizationVolFlatHRPower::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pVolatility, m_pHazardRate) );

  return pTM;
}

} // namespace ihg

} // namespace ito33
