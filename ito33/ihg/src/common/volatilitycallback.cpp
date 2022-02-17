/////////////////////////////////////////////////////////////////////////////
// Name:        ihg/src/common/volatiliycallback.cpp
// Purpose:     implementations of call back volatility class
// Author:      Wang
// Created:     2004/03/18
// RCS-ID:      $Id: volatilitycallback.cpp,v 1.13 2006/08/20 09:31:04 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ihg/src/common/volatilitycallback.cpp
 */

#include "ito33/useexception.h"

#include "ito33/ihg/volatilitycallback.h"

extern const ito33::Error ITO33_UNEXPECTED; 

using ito33::ihg::VolatilityCallBack;
using ito33::ihg::VolatilityCallBackBase;

//-----------------------------------------------------------------------------
// CallBackVolatility
//-----------------------------------------------------------------------------

void
VolatilityCallBackBase::GetVols(double dTime,
                                const double *pdS, 
                                double *pdVols, 
                                size_t nNbS) const
{
  DoGetValues(dTime, pdS, pdVols, nNbS);
  size_t nIdS;

  for ( nIdS = 0; nIdS < nNbS; nIdS++)
    pdVols[nIdS] += m_dShift;
}

//-----------------------------------------------------------------------------
// TODO: Dump for volatility call back 
//-----------------------------------------------------------------------------
void VolatilityCallBackBase::Dump(ito33::XML::Tag& ) const
{
  
}

//-----------------------------------------------------------------------------
// TODO: Visitor pattern for call back volatility
//-----------------------------------------------------------------------------
void VolatilityCallBackBase::Visit(ito33::ihg::VolatilityVisitor&) const
{
  throw EXCEPTION_MSG
        (
          ITO33_UNEXPECTED,
          "We are sorry, but we don't know how to visit a call back class."
        );
}

