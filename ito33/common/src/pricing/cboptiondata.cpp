/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/pricing/cboptiondata.cpp
// Author:      Nabil
// Created:     2005/06/14
// RCS-ID:      $Id: cboptiondata.cpp,v 1.12 2006/08/19 23:18:27 wang Exp $
// Copyright:   (c) 2005 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/autoptr.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/predicatetime.h"

#include "ito33/finance/numeraire.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/cboption.h"

#include "ito33/pricing/cboptiondata.h"

using namespace ito33::finance;

namespace ito33
{
  ITO33_IMPLEMENT_AUTOPTR(pricing::CBOptionData);
}

namespace ito33
{

namespace pricing
{

CBOptionData::CBOptionData(const finance::CBOption& cboption)
{
  m_dMaturityTime = GetDoubleFrom( cboption.GetMaturityDate() );
  
  m_dRedemptionRate = cboption.GetASWRedemptionRate();
  
  m_dCbFaceValue = cboption.GetCbNominal();
  
  m_dCbOptionFaceValue = cboption.GetASWNotional();

  m_dBalloonCoupon = cboption.GetBalloonCoupon();

  m_pFloatingRates = cboption.GetASWFloatingRates();

  // initialize the fixed cash flows
  shared_ptr<CashFlowStream> 
    pFixedCashFlowStream = cboption.ComputeASWFixedPayments();
  
  if ( pFixedCashFlowStream )
    m_pFixedCashFlows = make_ptr( new CashFlows( pFixedCashFlowStream,
                                                m_dCbFaceValue) );
  
  // initialize the floating cash flows
  m_pFloatingCashFlows = make_ptr( new CashFlows
                                      ( cboption.ComputeASWFloatingPayments(),
                                        m_dCbOptionFaceValue ) );
}

} //namespace pricing

} //namespace ito33
