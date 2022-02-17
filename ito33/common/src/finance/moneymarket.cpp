/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/moneymarket.cpp
// Purpose:     do the necessary for moneymarket class
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: moneymarket.cpp,v 1.16 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/moneymarket.h"
#include "ito33/finance/yieldcurve.h"
#include "ito33/finance/numeraire.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/moneymarket.h"

extern const ito33::finance::Error
  ITO33_MONEYMARKET_INVALID_YIELD_CURVE,
  ITO33_MONEYMARKET_INVALID_NUMERAIRE;

namespace ito33
{

namespace finance
{


MoneyMarket::MoneyMarket(const shared_ptr<Numeraire>& pCurrency,
                         const shared_ptr<YieldCurve>& pYieldCurve)
{
  // standard checks before setting
  CHECK_COND( pYieldCurve, ITO33_MONEYMARKET_INVALID_YIELD_CURVE);
  m_pYieldCurve = pYieldCurve;

  CHECK_COND( pYieldCurve, ITO33_MONEYMARKET_INVALID_NUMERAIRE);
  m_pNumeraire = pCurrency;
}


void MoneyMarket::SetYieldCurve(const shared_ptr<YieldCurve>& pYieldCurve)
{
  // standard check before setting
  CHECK_COND( pYieldCurve, ITO33_MONEYMARKET_INVALID_YIELD_CURVE);
  m_pYieldCurve = pYieldCurve;
}


void MoneyMarket::SetYieldCurveForMesh(const shared_ptr<YieldCurve>& yc) 
{ 
  // internal function, so assert instead of exception
  ASSERT_MSG( yc, "Invalid yield curve." );

  m_pYieldCurveForMesh = yc;
}


const shared_ptr<YieldCurve>& MoneyMarket::GetYieldCurveForMesh() const 
{ 
  // internal function, so assert instead of exception
  ASSERT_MSG( m_pYieldCurveForMesh, "Invalid yield curve for mesh.");

  return m_pYieldCurveForMesh; 
}


void MoneyMarket::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagMarket(XML_TAG_MONEYMARKET_ROOT, tagParent);

  tagMarket.Element(XML_TAG_MONEYMARKET_YIELDCURVE, *GetYieldCurve());

  tagMarket.Element(XML_TAG_MONEYMARKET_NUMERAIRE)(m_pNumeraire->GetCode());
}

} // namespace finance

} // namespace ito33
