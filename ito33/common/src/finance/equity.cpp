/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/equity.cpp
// Purpose:     Implemenation of the Equity class
// Author:      ZHANG Yunzhi
// Created:     Feb 09, 2004
// RCS-ID:      $Id: equity.cpp,v 1.11 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2003 - 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/dividends.h"
#include "ito33/finance/yieldcurve_flat.h"
#include "ito33/finance/equity.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/error.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/equity.h"
#include "ito33/xml/finance/dividends.h"
#include "ito33/xml/finance/common.h"

extern const ito33::finance::Error
  ITO33_EQUITY_INVALID_ISSUER,
  ITO33_EQUITY_INVALID_NUMERAIRE,
  ITO33_EQUITY_INVALID_DIVIDENDS,
  ITO33_EQUITY_INVALID_BORROWCURVE,
  ITO33_NEG_SPOT,
  ITO33_UNDEF_SPOT,
  ITO33_NEG_PREVIOUS_SHARE_PRICE;

namespace ito33 
{

namespace finance
{

// by default, equity has a zero borrow curve
Equity::Equity(double dSpotSharedPrice,
               const shared_ptr<Numeraire>& pCurrency)
  : m_dSpot(-1.),
    m_pBorrowCurve(new YieldCurveFlat(0.)),
    m_dPreviousSharePrice(-1.0)
{
  CHECK_COND( pCurrency, ITO33_EQUITY_INVALID_NUMERAIRE );

  m_pNumeraire = pCurrency;

  SetSpotSharePrice(dSpotSharedPrice);

  // Set a default issuer
  m_pIssuer = shared_ptr<Issuer>(new Issuer);
}

// by default, equity has a zero borrow curve
Equity::Equity(const shared_ptr<Numeraire>& pCurrency)
  : m_dSpot(-1.),
    m_pBorrowCurve(new YieldCurveFlat(0.)),
    m_dPreviousSharePrice(-1.0)
{
  CHECK_COND( pCurrency, ITO33_EQUITY_INVALID_NUMERAIRE );

  m_pNumeraire = pCurrency;

  // Set a default issuer
  m_pIssuer = shared_ptr<Issuer>(new Issuer);
}



void Equity::SetIssuer(const shared_ptr<Issuer>& pIssuer)
{
  CHECK_COND( pIssuer, ITO33_EQUITY_INVALID_ISSUER );

  m_pIssuer = pIssuer;  
}


void Equity::SetSpotSharePrice(double dSpot) 
{ 
  CHECK_COND( dSpot > 0, ITO33_NEG_SPOT );

  m_dSpot = dSpot; 
}

double Equity::GetSpotSharePrice() const
{ 
  CHECK_COND( m_dSpot > 0, ITO33_UNDEF_SPOT );

  return m_dSpot;
}

void Equity::SetPreviousSharePrice(double dPreviousSharePrice)
{
  CHECK_COND( dPreviousSharePrice > 0., ITO33_NEG_PREVIOUS_SHARE_PRICE);

  m_dPreviousSharePrice = dPreviousSharePrice;
}

void Equity::SetDividends(const shared_ptr<Dividends>& pDividends) 
{
  CHECK_COND( pDividends, ITO33_EQUITY_INVALID_DIVIDENDS );
 
  pDividends->Validate();

  m_pDividends = pDividends; 
}


void Equity::SetBorrowCurve(const shared_ptr<YieldCurve>& pBorrowCurve) 
{ 
  CHECK_COND( pBorrowCurve, ITO33_EQUITY_INVALID_BORROWCURVE );

  m_pBorrowCurve = pBorrowCurve; 
}


void Equity::Dump(XML::Tag& tagParent) const
{
  XML::Tag tagEquity(XML_TAG_EQUITY_ROOT, tagParent);

  tagEquity.Element(XML_TAG_EQUITY_SPOTSHAREPRICE)(m_dSpot);  

  tagEquity.Element(XML_TAG_CURRENCY)(m_pNumeraire->GetCode());

  tagEquity.Element(XML_TAG_EQUITY_BORROWCURVE, *m_pBorrowCurve);

  if ( m_dPreviousSharePrice > 0)
    tagEquity.Element(XML_TAG_EQUITY_PREVIOUS_SHARE_PRICE)
                     (m_dPreviousSharePrice);

  if ( m_pIssuer )
    m_pIssuer->Dump(tagEquity);

  if ( m_pDividends )
    tagEquity.Element(*m_pDividends);
  
}

} // namespace finance

} // namespace ito33
