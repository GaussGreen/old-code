/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatiliyflat.cpp
// Purpose:     implementations of flat volatility class
// Author:      Wang
// Created:     2004/03/17
// RCS-ID:      $Id: volatilityflat.cpp,v 1.18 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/common/volatilityflat.cpp
 */

#include "ito33/useexception.h"

#include "ito33/finance/modelparametersconsumer.h"

#include "ito33/ihg/volatilityflat.h"
#include "ito33/ihg/volatility_visitor.h"

#include "ito33/xml/write.h"

#include "ihg/xml/volatility.h"

extern const ito33::Error ITO33_BAD_PARAM;

namespace ito33
{
  
namespace ihg
{


VolatilityFlat::VolatilityFlat(double dValue) : Volatility()
{
  if ( dValue < 0 || dValue > 5 )
  {
    throw EXCEPTION_MSG
          (
            ITO33_BAD_PARAM,
            TRANS("Flat volatility can't be negative or greater than 500%!")
          );
  }

  m_dValue = dValue;
}

void VolatilityFlat::GetVols(double /*  dTime  */, 
                             const double*  /*  pdS  */, 
                             double *pdVols, 
                             size_t nNbS) const
{
  size_t nIdxS;
  for ( nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdVols[nIdxS] = m_dValue;
}

void VolatilityFlat::GetVolsSquared(double /*  dTime  */, 
                                    const double* /*  pdS  */, 
                                    double *pdVolsSquared,
                                    size_t nNbS) const
{
  double dVolSquared = m_dValue * m_dValue;
  size_t nIdxS;

  for ( nIdxS = 0; nIdxS < nNbS; nIdxS++)
    pdVolsSquared[nIdxS] = dVolSquared;
}

void VolatilityFlat::Dump(ito33::XML::Tag& tagParent) const
{
  ito33::XML::Tag tag(XML_TAG_VOLATILITYFLAT_ROOT,tagParent);
  tag.Element(XML_TAG_VOLATILITY_FLAT)(m_dValue);
}


void VolatilityFlat::Visit(ito33::ihg::VolatilityVisitor& visitor) const
{
  visitor.OnVolatilityFlat(*this);
}

void VolatilityFlat::GetModelParameters
                     (finance::ModelParametersConsumer& visitor) const
{
  std::vector<finance::ScalarModelParameter> pScalarModelParam(1);
  pScalarModelParam[0].name = SCALAR_MODEL_PARAM_NAME_VALUE;
  pScalarModelParam[0].value = m_dValue;

  visitor.OnScalarValues(MODEL_PARAM_NAME_VOL_FLAT, pScalarModelParam);
}


  
} // namespace ihg

} // namespace ito33
