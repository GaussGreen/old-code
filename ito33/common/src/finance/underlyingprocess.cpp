/////////////////////////////////////////////////////////////////////////////
// Name:        src/finance/underlyingprocess.cpp
// Purpose:     Base class for underlying process
// Created:     2006/06/01
// RCS-ID:      $Id: underlyingprocess.cpp,v 1.2 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004-2006 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/underlyingprocess.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error ITO33_POST_DEFAULT_VOLATILITY_NEGATIVE,
                                   ITO33_POST_DEFAULT_VOLATILITY_TOO_LARGE;

namespace ito33
{

namespace finance
{

UnderlyingProcess::UnderlyingProcess()
                  :m_dPostDefaultVolatility(0.0)
{
}

void UnderlyingProcess::SetPostDefaultVolatility(double dPostDefaultVolatility)
{
  CHECK_COND( dPostDefaultVolatility >= 0., 
    ITO33_POST_DEFAULT_VOLATILITY_NEGATIVE);
  
  CHECK_COND( dPostDefaultVolatility < 5, 
    ITO33_POST_DEFAULT_VOLATILITY_TOO_LARGE);

  m_dPostDefaultVolatility = dPostDefaultVolatility;
}

void UnderlyingProcess::DumpMe(ito33::XML::Tag& tagParent) const
{
  if ( m_dPostDefaultVolatility > 0.0 )
    tagParent.Element(XML_TAG_POST_DEFAULT_VOLATILITY)
        (m_dPostDefaultVolatility);
}

} // namespace finance

} // namespace ito33
