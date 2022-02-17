/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatiliytanh.cpp
// Purpose:     implementations of volatility as a variant of the tanh function of spot
// Created:     2005/02/04
// RCS-ID:      $Id: volatilitytanh.cpp,v 1.12 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/useexception.h"

#include "ito33/ihg/volatilitytanh.h"
#include "ito33/ihg/volatility_visitor.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/xml/write.h"
#include "ihg/xml/volatility.h"

extern const ito33::Error ITO33_BAD_DATA;

namespace ito33
{

namespace ihg
{


VolatilityTanh::VolatilityTanh(double dLeft, double dRight, 
                               double dScale, double dS0)
                             : m_dLeft(dLeft), m_dRight(dRight),
                               m_dScale(dScale), m_dS0(dS0)
{
  CHECK_COND_MSG(dLeft >= 0.,
                 ITO33_BAD_DATA,
                 "The left limit of the tanh function must not be negative.");

  CHECK_COND_MSG(dRight >= 0.,
                 ITO33_BAD_DATA,
                 "The right limit of the tanh function must not be negative.");

  CHECK_COND_MSG(dScale >= 0.,
                 ITO33_BAD_DATA,
                 "The scale must be positive.");

  CHECK_COND_MSG(m_dS0 > 0,
                 ITO33_BAD_DATA,
                 "The shift to the tanh function must be positive.");
}

void VolatilityTanh::GetVols(double /*  dTime  */, 
                             const double* pdS, 
                             double *pdVols,
                             size_t nNbS) const
{
  for (size_t nIdx = 0; nIdx < nNbS; nIdx++)
  {
    pdVols[nIdx] = 0.5 * (m_dRight - m_dLeft) 
                       * ( 1. + tanh( m_dScale * (pdS[nIdx] - m_dS0) ) )
                 + m_dLeft; 
  }
}

void VolatilityTanh::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_VOLATILITYTANH_ROOT,tagParent);
 
  tag.Element(XML_TAG_VOLATILITYTANH_LEFT)(m_dLeft);
  tag.Element(XML_TAG_VOLATILITYTANH_RIGHT)(m_dRight);
  tag.Element(XML_TAG_VOLATILITYTANH_SCALE)(m_dScale);
  tag.Element(XML_TAG_VOLATILITYTANH_S0)(m_dS0);
}

void VolatilityTanh::Visit(VolatilityVisitor& visitor) const
{
  visitor.OnVolatilityTanh(*this);
}

void VolatilityTanh::GetModelParameters
                     (finance::ModelParametersConsumer& visitor) const
{
  std::vector<finance::ScalarModelParameter> pScalarModelParam(4);
  
  pScalarModelParam[0].name = SCALAR_MODEL_PARAM_NAME_LEFT_LIMIT;
  pScalarModelParam[0].value = m_dLeft;

  pScalarModelParam[1].name = SCALAR_MODEL_PARAM_NAME_RIGHT_LIMIT;
  pScalarModelParam[1].value = m_dRight;

  pScalarModelParam[2].name = SCALAR_MODEL_PARAM_NAME_SCALE;
  pScalarModelParam[2].value = m_dScale;
  pScalarModelParam[2].expressedInPercents = false;

  pScalarModelParam[3].name = SCALAR_MODEL_PARAM_NAME_S0;
  pScalarModelParam[3].value = m_dS0;
  pScalarModelParam[3].expressedInPercents = false;

  visitor.OnScalarValues(MODEL_PARAM_NAME_VOL_TANH, pScalarModelParam);

}


} // namespace ihg

} // namespace ito33
