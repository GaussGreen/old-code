/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/bondlike/convertiblelike.cpp
// Purpose:     financial convertible-like bond class
// Author:      ITO33
// Created:     2004/10/13
// RCS-ID:      $Id: cb_base.cpp,v 1.34 2006/08/19 22:43:51 wang Exp $
// Copyright:   (c) 2004 - 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/sharedptr.h"
#include "ito33/useexception.h"

#include "ito33/finance/error.h"
#include "ito33/finance/issuer.h"
#include "ito33/finance/numeraire.h"
#include "ito33/finance/cashflowstream.h"
#include "ito33/finance/sessiondata.h"

#include "ito33/finance/bondlike/bonderror.h"
#include "ito33/finance/bondlike/callschedule.h"
#include "ito33/finance/bondlike/putschedule.h"
#include "ito33/finance/bondlike/convertiblelike.h"
#include "ito33/finance/bondlike/cb_base.h"
#include "ito33/finance/bondlike/bond.h"
#include "ito33/finance/bondlike/utils.h"

#include "ito33/xml/write.h"
#include "ito33/xml/finance/bondlike/convertiblelike.h"

extern const ito33::finance::Error  ITO33_INVALID_FREQUENCY;

extern const ito33::finance::BondError
  ITO33_BONDLIKE_NULL_BONDTERMS,
  ITO33_BONDLIKE_NULL_CALLSCHEDULE,
  ITO33_BONDLIKE_NULL_PUTSCHEDULE,
  ITO33_BONDLIKE_NO_CALL,
  ITO33_BONDLIKE_NO_PUT,
  ITO33_BONDLIKE_CANNOT_READ_BONDLIKETERMS,
  ITO33_BONDLIKE_STRAIGHTBOND_FOR_EXCHANGEABLE;

namespace ito33
{

namespace finance
{


CBBase::CBBase(const shared_ptr<BondTerms>& pBondTerms)
  : ConvertibleLike(CHECK_PTR(pBondTerms, ITO33_BONDLIKE_NULL_BONDTERMS)),
    m_pBondTerms(pBondTerms)
{
}

void CBBase::SetCallSchedule(const shared_ptr<CallSchedule>& pCallSchedule)
{
  m_pCallSchedule = CHECK_PTR(pCallSchedule, ITO33_BONDLIKE_NULL_CALLSCHEDULE);
}

void CBBase::SetPutSchedule(const shared_ptr<PutSchedule>& pPutSchedule)
{
  m_pPutSchedule = CHECK_PTR(pPutSchedule, ITO33_BONDLIKE_NULL_PUTSCHEDULE);
}

double 
CBBase::ComputeYieldToMaturity
(double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const
{
  CheckSessionData();

  CHECK_COND( IsValid(cmpFrequency), ITO33_INVALID_FREQUENCY );

  Bond bond( m_pBondTerms );
  bond.SetSessionData( GetSessionData() );

  return bond.ComputeYieldToMaturity(dPrice, cmpFrequency, dcc);
}

double 
CBBase::ComputeYieldToPut
(double dPrice, Frequency cmpFrequency, Date::DayCountConvention dcc) const
{
  CheckSessionData();

  CHECK_COND( IsValid(cmpFrequency), ITO33_INVALID_FREQUENCY );

  CHECK_COND_MSG
  (
    m_pPutSchedule,
    ITO33_BONDLIKE_NULL_PUTSCHEDULE,
    "Can't compute a yield-to-put without a put schedule"
  );

  Bond bond( m_pBondTerms );
  bond.SetSessionData( GetSessionData() );  
  
  bond.SetPutSchedule( m_pPutSchedule );

  return bond.ComputeYieldToPut(dPrice, cmpFrequency, dcc);
}

double CBBase::ComputePutPrice(Date putDate) const
{
  CheckSessionData();

  CHECK_COND(m_pPutSchedule, ITO33_BONDLIKE_NO_PUT);

  Bond bond( m_pBondTerms );
  bond.SetSessionData( GetSessionData() );
  
  bond.SetPutSchedule( m_pPutSchedule );

  return bond.ComputePutPrice(putDate);
}

double CBBase::ComputePutPrice(Date putDate, double& dCouponAmount) const
{
  CheckSessionData();

  CHECK_COND(m_pPutSchedule, ITO33_BONDLIKE_NO_PUT);

  Bond bond( m_pBondTerms );
  bond.SetSessionData( GetSessionData() );
  
  bond.SetPutSchedule( m_pPutSchedule );

  return bond.ComputePutPrice(putDate, dCouponAmount);
}

double CBBase::ComputeCallPrice(Date callDate) const
{ 
  CheckSessionData();

  CHECK_COND(m_pCallSchedule, ITO33_BONDLIKE_NO_CALL);

  return ComputeBondLikeCallPrice( callDate, *m_pBondTerms, *GetSessionData(), 
                                   *m_pCallSchedule, m_pNumeraire );
}

double CBBase::GetAccruedInterestValue() const
{
  CheckSessionData();

  shared_ptr<CashFlowStream>
    pCoupons = ComputeCouponRates(*m_pBondTerms, *GetSessionData(), m_pNumeraire);

  if (pCoupons)
    return pCoupons->GetAccrued(GetSessionData()->GetValuationDate())
          * GetBondTerms()->GetNominal();
  else
    return 0;
}

shared_ptr<Bond> CBBase::GetStraightBond() const
{
  CHECK_COND(!IsExchangeable(), ITO33_BONDLIKE_STRAIGHTBOND_FOR_EXCHANGEABLE);

  shared_ptr<Bond> pBond( new Bond(m_pBondTerms) );

  if ( m_pSessionData )
    pBond->SetSessionData(m_pSessionData);

  return pBond;
}

void CBBase::Validate() const
{
  Derivative::Validate();

  m_pBondTerms->Validate();
  
  if ( m_pCallSchedule )
    m_pCallSchedule->Validate();

  if ( m_pPutSchedule )
    m_pPutSchedule->Validate();

  CheckYieldCompoundingFrequency
      (m_pPutSchedule, m_pCallSchedule, *m_pBondTerms);

  CheckYieldDayCountConvention
      (m_pPutSchedule, m_pCallSchedule, *m_pBondTerms);
}

void CBBase::DumpMe(XML::Tag& tagParent) const
{
  ConvertibleLike::DumpMe(tagParent);

  if ( m_pCallSchedule && !( m_pCallSchedule->GetAll().empty() ) )
    tagParent.Element(*m_pCallSchedule);

  if ( m_pPutSchedule && !( m_pPutSchedule->GetAll().empty() ) )
    tagParent.Element(*m_pPutSchedule);
}


} // namespace finance

} // namespace ito33

