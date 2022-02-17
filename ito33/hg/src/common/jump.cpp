/////////////////////////////////////////////////////////////////////////////
// Name:        hg/src/common/jump.cpp
// Purpose:     implementation of HG Jump class
// RCS-ID:      $Id: jump.cpp,v 1.7 2006/08/19 23:44:27 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/hg/error.h"
#include "ito33/hg/jump.h"

#include "ito33/xml/write.h"

extern const ito33::hg::Error
  ITO33_HG_INTENSITY,
  ITO33_HG_AMPLITUDE;

namespace ito33
{

namespace hg
{


Jump::Jump(double dIntensity, double dAmplitude)
         : m_dIntensity(dIntensity), m_dAmplitude(dAmplitude)
{
  /* CHECK_COND(m_dIntensity >= 0, ITO33_HG_INTENSITY); */

  CHECK_COND(m_dAmplitude > -1., ITO33_HG_AMPLITUDE);
}


} // namespace hg

} // namespace ito33
