/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatiliypower.cpp
// Purpose:     implementations of spot power  volatility
// Author:      Ito33
// Created:     2004/12/13
// RCS-ID:      $Id: volatilitypower.cpp,v 1.10 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/common/volatilitypower.cpp
 */

#include <cmath>

#include "ito33/useexception.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/volatilitypower.h"
#include "ito33/ihg/volatility_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/volatility.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{

namespace ihg
{

VolatilityPower::VolatilityPower(double dAlpha,double dBeta,double dS0)
                               : Volatility()
{

  if ( dAlpha < 0.0 )		 
    throw EXCEPTION_MSG
        (
          ITO33_BAD_PARAM,
          TRANS("Parameter alpha must be non-negative.")
        );

  if ( dS0 < 0.0 )
    throw EXCEPTION_MSG
        (
          ITO33_BAD_PARAM,
          TRANS("Parameter S0 must be non-negative.")
        );


  m_dAlpha = dAlpha;
  m_dBeta  = dBeta;	
  m_dS0    = dS0;
}

void VolatilityPower::GetVols(double /*  dTime  */, 
                              const double* pdS, 
                              double *pdVols,
                              size_t nNbS) const
{
  size_t nIdx;


  for (nIdx = 0; nIdx < nNbS; nIdx++)
  {  
    double dTmp;
   
    if (pdS[nIdx] < 1.e-16)
      dTmp = 1.e16;
    else
      dTmp = 1. / pdS[nIdx];

    pdVols[nIdx] = m_dAlpha * pow(m_dS0 * dTmp, m_dBeta);
  }//end for loop

} // VolatilitySpotPower::GetVols


void VolatilityPower::GetVolsSquared(double /*  dTime  */, 
                                    const double*   pdS  , 
                                    double *pdVolsSquared,
                                    size_t nNbS) const
{	

  double dAlphaSquared = m_dAlpha * m_dAlpha;
  size_t nIdx;

  for (nIdx = 0; nIdx < nNbS; nIdx++)
	{
    double dTmp;

    if (pdS[nIdx] < 1.e-16)
			dTmp = 1.e16;
		else
			dTmp = 1. / pdS[nIdx];

    pdVolsSquared[nIdx] = dAlphaSquared * pow(m_dS0 * dTmp, 2.*m_dBeta);
  } //end for
    
}//VolatilitySpotPower::GetVolsSquared

void VolatilityPower::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_VOLATILITYPOWER_ROOT,tagParent);
  tag.Element(XML_TAG_VOLATILITYPOWER_ALPHA)(m_dAlpha);
  tag.Element(XML_TAG_VOLATILITYPOWER_BETA)(m_dBeta);
  tag.Element(XML_TAG_VOLATILITYPOWER_S0)(m_dS0);
}//VolatilitySpotPower::Dump


void VolatilityPower::Visit(VolatilityVisitor& visitor) const
{
  visitor.OnVolatilityPower(*this);
}//VolatilitySpotPower::Visit

void VolatilityPower::GetModelParameters
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

  visitor.OnScalarValues(MODEL_PARAM_NAME_VOL_POWER, pScalarModelParam);
}


} // namespace ihg
}// namespace ito33
