/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/calibration/parametrization_voltimeonly_hrwithtimecomponent.cpp
// Purpose:     implement ParametrizationVolTimeOnlyHRWithTimeComponent class
// Author:      ITO 33
// Created:     2005/07/14
// RCS-ID:      $Id: parametrization_voltimeonly_hrwithtimecomponent.cpp,v 1.19 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
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

#include "ihg/voltimeonlyhrtimecomponentcalibrator.h"

#include "ito33/ihg/hazardratetimeonly.h"
#include "ito33/ihg/hazardratecombo.h"
#include "ito33/ihg/spotcomponent.h"
#include "ito33/ihg/hrspotcomponentpower.h"
#include "ito33/ihg/volatilitytimeonly.h"
#include "ito33/ihg/version.h"
#include "ito33/ihg/error.h"
#include "ito33/ihg/theoreticalmodel.h"

#include "ito33/ihg/parametrization_voltimeonly_hrwithtimecomponent.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/common.h"
#include "ihg/xml/common.h"
#include "ihg/xml/parametrization.h"
#include "ihg/xml/parametrization_visitor.h"
#include "ihg/xml/spotcomponent.h"

extern const ito33::Error
  ITO33_BAD_PARAM;

extern const ito33::finance::Error
  ITO33_CALIBRATION_FAIL,
  ITO33_EMPTY_DERIVATIVELIST;

extern const ito33::ihg::Error
  ITO33_IHG_INCORRECT_DERIVATIVELIST;

namespace std
{

template<>
bool greater< ito33::shared_ptr<ito33::finance::Derivative> >::operator()
     (const ito33::shared_ptr<ito33::finance::Derivative> & pDeriv1, 
      const ito33::shared_ptr<ito33::finance::Derivative> & pDeriv2) const
{
  return pDeriv1->GetMaturityDate() < pDeriv2->GetMaturityDate();
}

} // namespace std

namespace ito33
{

namespace ihg
{

void ParametrizationVolTimeOnlyHRWithTimeComponent::SetSpotComponent
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

void ParametrizationVolTimeOnlyHRWithTimeComponent::CalibrateWithDerivatives
     (const finance::Derivatives& derivatives) 
{

  const finance::Derivatives::Elements& elements( derivatives.GetAll() );
  size_t nNbDerivatives = elements.size();

  CHECK_COND(nNbDerivatives > 0, ITO33_EMPTY_DERIVATIVELIST);

  // Make sure the derivative list is sorted, and that maturities define
  // period pairs:
  //   maturity[0] <= maturity[1] < maturity[3] <= maturity[4] < etc  
  // The period pairs are then 0/1, 2/3, 4/5, etc and the second
  // maturities of each pair are increasing. Alternatively, we
  // need an even number of derivatives, and no 3 (or more) can have the
  // same maturity date.
  CHECK_COND(nNbDerivatives % 2 == 0,
             ITO33_IHG_INCORRECT_DERIVATIVELIST);

  std::list< shared_ptr<finance::Derivative> > derivListTmp;

  // since derivatives is const, must copy before sorting
  finance::Derivatives::Elements::const_iterator iter;
  for (iter = elements.begin(); iter != elements.end(); ++iter)
    derivListTmp.push_back(iter->first);

  derivListTmp.sort( std::greater< shared_ptr<finance::Derivative> >() );

  // Check ordering
  std::list< shared_ptr<finance::Derivative> >::const_iterator iterTmp;
  Date lastDate = derivListTmp.front()->GetMaturityDate();
  lastDate.AddYears(-1);
  size_t nCounter;
  for (iterTmp = derivListTmp.begin(), nCounter = 0; 
       iterTmp != derivListTmp.end(); 
       ++iterTmp, nCounter++)
  {
    // pairs are 0/1, 2/3, 4/5, 6/7, etc
    // 'even' numbers must be strictly greater
    // 'odd' numbers must be greater or equal (redundant check due to sort)
    if (nCounter % 2 == 0)
    {
      CHECK_COND( (*iterTmp)->GetMaturityDate() > lastDate, 
                  ITO33_IHG_INCORRECT_DERIVATIVELIST );
    }
    else
    {
      CHECK_COND( (*iterTmp)->GetMaturityDate() >= lastDate, 
                  ITO33_IHG_INCORRECT_DERIVATIVELIST );
    }

    lastDate = (*iterTmp)->GetMaturityDate();

  } // checking that derivs form valid pairs


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

  if ( m_pSpotComponent )
    m_pSpotComponent->TrySetS0
      ( derivatives.GetSessionData()->GetSpotSharePrice() );

  // Create the combo hazard rate if spot component is defined
  shared_ptr<HazardRateWithTimeComponent> pHazardRate;

  if ( m_pSpotComponent ) // hazard rate is not time only
    pHazardRate = make_ptr( new HazardRateCombo(m_pSpotComponent) ); 

  VolTimeOnlyHRTimeComponentCalibrator calibrator;
  
  try
  {
    calibrator.Calibrate(derivListTmp, pHazardRate);
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

void ParametrizationVolTimeOnlyHRWithTimeComponent::Calibrate
     (const finance::BasketGoodType& basket)
{
  ASSERT( basket.GetDerivatives() );

  CalibrateWithDerivatives( *basket.GetDerivatives() );
}

void ParametrizationVolTimeOnlyHRWithTimeComponent::Dump
    (ito33::XML::Tag& tagParent) const
{
  // create the root tag
  ito33::XML::Tag tagParametrization
    (XML_TAG_PARAMETRIZATION_VOLTIMEONLYHRWITHTIMECOMPONENT_ROOT, tagParent);
  
  // only output data if it is set
  if (m_pSpotComponent)
    tagParametrization.Element(XML_TAG_SPOTCOMPONENT_ROOT, *m_pSpotComponent);
}

void ParametrizationVolTimeOnlyHRWithTimeComponent::Visit
    (ParametrizationVisitor& visitor) const
{
  visitor.OnParametrizationVolTimeOnlyHRWithTimeComponent(*this);
}

shared_ptr<finance::TheoreticalModel> 
ParametrizationVolTimeOnlyHRWithTimeComponent::GetTheoreticalModel()
{
  shared_ptr<TheoreticalModel> 
    pTM( new TheoreticalModel(m_pVolatility, m_pHazardRate) );

  return pTM;
}


} // namespace ihg

} // namespace ito33
