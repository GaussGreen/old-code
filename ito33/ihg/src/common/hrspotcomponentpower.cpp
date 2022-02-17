/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hrspotcomponentpower.cpp
// Purpose:     hazard rate as a power function of spot
// Author:      Wang 
// Created:     2004/06/11
// RCS-ID:      $Id: hrspotcomponentpower.cpp,v 1.19 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/useexception.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/hrspotcomponentpower.h"

#include "ito33/xml/write.h"

#include "ihg/xml/spotcomponent.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace ihg
{

void HRSpotComponentPower::DoValidate()
{
  if ( m_dBeta < 0.0 )
  {
    throw EXCEPTION_MSG
      (
       ITO33_BAD_PARAM,
       TRANS("Parameter beta must be non-negative in HazardSpotComponentPower.")
      );
  }

  if ( m_dBeta > 2.0 )
  {
    throw EXCEPTION_MSG
      (
       ITO33_BAD_PARAM,
       TRANS("Parameter beta must less than 2 in HazardSpotComponentPower.")
      );
  }  
}

HRSpotComponentPower::HRSpotComponentPower(double dBeta, double dS0)
    : m_dBeta(dBeta), SpotComponent(dS0)
{
  DoValidate();

  if ( m_dS0 < 0.0 )
  {
    throw EXCEPTION_MSG
      (
       ITO33_BAD_PARAM,
       TRANS("Parameter S0 must be non-negative in HazardSpotComponentPower.")
      );
  }
}

HRSpotComponentPower::HRSpotComponentPower(double dBeta)
    : m_dBeta(dBeta), SpotComponent(-100)
{
  DoValidate(); 
}

void HRSpotComponentPower::GetValues
     (const double *pdS, double *pdValues, size_t nNbS) const
{
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    double dTmp;

    if (pdS[nIdx] < 1.e-16)
      dTmp = 1.e16;
    else
      dTmp = 1. / pdS[nIdx];

    pdValues[nIdx] = pow(m_dS0 * dTmp, m_dBeta);
  }
}

ito33::XML::Tag
HRSpotComponentPower::Dump(const char *name, ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tagRoot(name, tagParent);
  ito33::XML::Tag tag(XML_TAG_SPOTCOMPONENT_POWER_ROOT, tagRoot);
  tag.Element(XML_TAG_SPOTCOMPONENT_BETA)(m_dBeta);
  tag.Element(XML_TAG_SPOTCOMPONENT_S0)(m_dS0);

  return tagRoot;
}
 
void HRSpotComponentPower::GetModelParameters
  (finance::ModelParametersConsumer& visitor, const char* componentName) const
{
    std::vector<finance::ScalarModelParameter> pScalarModelParam(2);
  
  pScalarModelParam[0].name = SCALAR_MODEL_PARAM_NAME_BETA;
  pScalarModelParam[0].value = m_dBeta;
  pScalarModelParam[0].expressedInPercents = false;

  pScalarModelParam[1].name = SCALAR_MODEL_PARAM_NAME_S0;
  pScalarModelParam[1].value = m_dS0;
  pScalarModelParam[1].expressedInPercents = false;

  visitor.OnScalarValues(componentName, pScalarModelParam);

}


} // namespace ihg

} // namespace ito33
