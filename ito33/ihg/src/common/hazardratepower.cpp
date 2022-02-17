/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/hazardratepower.cpp
// Purpose:     hazard rate as a power function of spot
// Author:      Wang 
// Created:     2004/06/11
// RCS-ID:      $Id: hazardratepower.cpp,v 1.19 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/useexception.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/hazardratepower.h"
#include "ito33/ihg/hazardrate_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/hazardrate.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace ihg
{

HazardRatePower::HazardRatePower(double dAlpha, double dBeta, double dS0)
  : m_dAlpha(dAlpha), m_dBeta(dBeta), m_dS0(dS0)
{

  if ( dAlpha < 0.0 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Parameter alpha must be non-negative in HazardRatePower.")
          );
  }

  if ( dS0 < 0.0 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Parameter S0 must be non-negative in HazardRatePower.")
          );
  }

  if ( dBeta < 0.0 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Parameter beta must be non-negative in HazardRatePower.")
          );
  }

  if ( dBeta > 2.0 )
  {
        throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Parameter beta must less than 2 in HazardRatePower.")
          );
  }
}


void HazardRatePower::GetHazardRates
     (const double *pdS, double *pdValues, size_t nNbS) const
{
  size_t nIdx;

  for (nIdx = 0; nIdx < nNbS; nIdx++)
  {
    double dTmp;

    if (pdS[nIdx] < 1.e-16)
      dTmp = 1.e16;
    else
      dTmp = 1. / pdS[nIdx];

    pdValues[nIdx] = m_dAlpha * pow(m_dS0 * dTmp, m_dBeta);
  }
}

void HazardRatePower::GetHazardRates(double /* dTime */, 
                                       const double* pdS,
                                       double *pdValues, 
                                       size_t nNbS) const
{
  GetHazardRates(pdS, pdValues, nNbS);
}

void HazardRatePower::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_HAZARDRATEPOWER_ROOT, tagParent);
  /* ctor: HazardRatePolynom(double dAlpha, double dBeta, double dS0) */
  /* \lambda (S,t) = \alpha (S_0 / S)^\beta    */
  tag.Element(XML_TAG_HAZARDRATE_ALPHA)(m_dAlpha);
  tag.Element(XML_TAG_HAZARDRATE_BETA)(m_dBeta);
  tag.Element(XML_TAG_HAZARDRATE_S0)(m_dS0);
}


void HazardRatePower::Visit(ito33::ihg::HazardRateVisitor& visitor) const
{
  visitor.OnHazardRatePower(*this);
}

void HazardRatePower::GetModelParameters
                      (finance::ModelParametersConsumer& visitor) const
{
  std::vector<finance::ScalarModelParameter> pScalarModelParam(3);
 
  pScalarModelParam[0].name = SCALAR_MODEL_PARAM_NAME_ALPHA;
  pScalarModelParam[0].value = m_dAlpha;

  pScalarModelParam[1].name = SCALAR_MODEL_PARAM_NAME_BETA;
  pScalarModelParam[1].value = m_dBeta;
  pScalarModelParam[1].expressedInPercents = false;

  pScalarModelParam[2].name = SCALAR_MODEL_PARAM_NAME_S0;
  pScalarModelParam[2].value = m_dS0;
  pScalarModelParam[2].expressedInPercents = false;

  visitor.OnScalarValues(MODEL_PARAM_NAME_HR_POWER, pScalarModelParam);
}


} // namespace ihg

} // namespace ito33
