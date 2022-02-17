/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cboption.cpp
// Author:      Nabil
// Created:     2004/10/13
// RCS-ID:      $Id: cboption.cpp,v 1.4 2006/07/26 21:24:12 wang Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"

#include "ito33/pricing/cboption.h"

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::CBOption);
}

namespace ito33
{

namespace pricing
{

CBOption::CBOption(const finance::CBOption& cboption)
{
  m_cb = CB( *cboption.GetConvertibleBond() );

  // Initialize the contract data
  // Use the maturity time of the cb but not the one from cboption!
  m_dMaturityTime = m_cb.GetMaturityTime();
  m_pDerivativeCurve = m_cb.GetDerivativeCurve();

  // Gets the cb option data
  m_pCBOptionData = AutoPtr<CBOptionData>( new CBOptionData(cboption) );
}

} //namespace pricing

} //namespace ito33
