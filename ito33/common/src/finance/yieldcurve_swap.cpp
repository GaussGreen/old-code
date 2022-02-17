/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/yieldcurve_swap.cpp
// Purpose:     YieldCurveSwap class implementation
// Author:      ZHANG Yunzhi
// Created:     2005/04/11
// RCS-ID:      $Id: yieldcurve_swap.cpp,v 1.14 2006/08/23 10:23:05 wang Exp $
// Copyright:   (c) 2004 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include <cmath>

#include "ito33/beforestd.h"
#include <algorithm>
#include <list>
#include "ito33/afterstd.h"

#include "ito33/enum_values_names.h"
#include "ito33/useexception.h"
#include "ito33/constants.h"
#include "ito33/dateutils.h"

#include "ito33/numeric/interpolation.h"
#include "ito33/numeric/extrapolationmode.h"

#include "ito33/finance/error.h"
#include "ito33/finance/yieldcurve_swap.h"
#include "ito33/finance/yieldcurve_visitor.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/yieldcurve.h"
#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/timeunit.h"
#include "ito33/xml/finance/frequency.h"
#include "ito33/xml/finance/daycountconvention.h"

extern const ito33::Error
  ITO33_UNEXPECTED, ITO33_BAD_PARAM;

extern const ito33::finance::Error
  ITO33_SWAP_RATE_MATURITY_UNIT,
  ITO33_SWAP_RATE_WRONG_FREQUENCY,
  ITO33_YIELDCURVELEG_MATURITY,
  ITO33_YIELDCURVE_SWAP_NODATA,
  ITO33_YIELDCURVE_SWAP_CASHRATES_MATURITES,
  ITO33_YIELDCURVE_SWAP_SWAPRATES_MATURITES,
  ITO33_YIELDCURVE_SWAP_INCONSISTENCY_SWAPCASH,
  ITO33_BOOTSTRAPPING_FAILED;

namespace ito33
{

namespace finance
{

extern void DumpPart(const YieldCurveLeg& part, XML::Tag& tagObject);

extern XML::Tag 
Dump(const YieldCurveLeg& leg, const char* legName, XML::Tag& tagParent);

SwapRate::SwapRate( double dRate,
            size_t nMaturityDuration,
            TimeUnit maturityUnit,
            Frequency paymentFrequency
            )
            : YieldCurveLeg(dRate, nMaturityDuration, maturityUnit),
              m_paymentFrequency(paymentFrequency)
{
  switch (maturityUnit)
  {
  case TimeUnit_Month:
    m_nNumberMonths = nMaturityDuration;
    break;
  case TimeUnit_Year:
    m_nNumberMonths = nMaturityDuration * 12;
    break;
  default:
    throw EXCEPTION(ITO33_SWAP_RATE_MATURITY_UNIT);
  }

  CHECK_COND( (m_nNumberMonths * paymentFrequency) % 12 == 0,
              ITO33_SWAP_RATE_WRONG_FREQUENCY);
}


// this function is used for equality of two dates expressed as doubles (we
// can't compare them using ==, of course...)
static inline bool
IsEqual(const SwapRate& leg1, const SwapRate& leg2)
{
  return leg1.GetNumberMonths() == leg2.GetNumberMonths();
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

inline double GetITO33YearsDiff(const Date& date1, const Date& date2)
{
  return GetDoubleFrom(date2) - GetDoubleFrom(date1);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

/**
   Helps to calculate do bootstrapping.

   In fact, this class stocks and calculates, for a certain term, the sum of
   discount factors at the terms before the former and with a fixed frequency.

           m_dValue = sum_{term < nNbMonths} { DF(term) }
 */
class SumDiscountFactor
{
public:

  SumDiscountFactor(Frequency paymentFrequency)
    : m_paymentFrequency(paymentFrequency),
      m_nNbMonthsBase(12 / m_paymentFrequency),
      m_nNbMonths(0),
      m_dValue(0)
  {
  }

  Frequency GetPaymentFrequency() const
  {
    return m_paymentFrequency;
  }

  double GetValue() const
  {
    return m_dValue;
  }

  // update values (sum, DF, zero-coupon rate) for new maturity term
  void UpdateValues
            (
              std::vector<size_t>& pnTerms,
              std::vector<double>& pdRates,
              Date dateReference,
              size_t nNbMonths,
              double dPaymentRate
            )
  {
    ASSERT_MSG( nNbMonths > m_nNbMonths,
                "SumDiscountFactor::UpdateValues(): "
                "Given maturity term must be longer than known one");

    ASSERT_MSG( nNbMonths % m_nNbMonthsBase == 0,
                "SumDiscountFactor::UpdateValues(): "
                "Given maturity term is not consistent with frequency");

    if ( (nNbMonths - m_nNbMonthsBase) > m_nNbMonths)
    {
      // we have to interpolate to have the discount factor at missing points

      // allocation
      size_t
        nIdx,
        nNbNewTerms = (nNbMonths - m_nNbMonths) / m_nNbMonthsBase - 1;

      YieldCurveAnnuallyCompounded yieldCurve(dateReference, pnTerms.size());
      yieldCurve.SetLegs(pnTerms, pdRates);

      std::vector<double>
        pdTermDays(nNbNewTerms),
        pdTermDFs(nNbNewTerms);

      for(m_nNbMonths += m_nNbMonthsBase, nIdx = 0;
          m_nNbMonths <= (nNbMonths - m_nNbMonthsBase);
          m_nNbMonths += m_nNbMonthsBase, nIdx++)
      {
        pdTermDays[nIdx]
            = GetDoubleFrom
                ( AddMonthsAdjustedForEndOfMonth  // term date
                    ( dateReference, static_cast<int>(m_nNbMonths) ) );
      }
      // at this point, m_nNbMonths is nNbMonths+m_nNbMonthsBase, so the real
      // nNbMonths input

      yieldCurve.GetDiscountFactor(&pdTermDays[0], &pdTermDFs[0],
                                   nNbNewTerms);

      for (nIdx = 0; nIdx < nNbNewTerms; nIdx++)
        m_dValue += pdTermDFs[nIdx];
    }

    // update members
    m_nNbMonths = nNbMonths;
    double dDF = (1 - m_dValue * dPaymentRate) / (1 + dPaymentRate);

    // we should throw this exception if bootstrapping fails because of a
    // problem with the values used and if it still makes sense to retry
    // bootstrapping with different values -- if it is (ever) not the case,
    // another error code should be used
    CHECK_COND( dDF > 0, ITO33_BOOTSTRAPPING_FAILED );

    m_dValue += dDF;

    // update pdTerms and pdRates
    Date dateNew = AddMonthsAdjustedForEndOfMonth(dateReference, m_nNbMonths);
    size_t nDiff = size_t(Date::DaysDiff(dateReference, dateNew));
    double dFractionITO33 = nDiff * ONEDAY;

    pnTerms.push_back(nDiff);
    pdRates.push_back( pow(dDF, -1/dFractionITO33) - 1);
  }

private:

  Frequency m_paymentFrequency;

  size_t m_nNbMonthsBase;

  size_t m_nNbMonths;

  double m_dValue;

};



class MyCashRate : public finance::CashRate
{
public:
  MyCashRate( const CashRate& cashRate, Date referenceDate)
      : CashRate(cashRate)
  {
    switch (m_maturityUnit)
    {
    case TimeUnit_Day:
      m_maturityDate = referenceDate.AddDays(m_nMaturityDuration);
      break;
    case TimeUnit_Month:
      m_maturityDate = AddMonthsAdjustedForEndOfMonth
                          (referenceDate, m_nMaturityDuration);
      break;
    case TimeUnit_Year:
      m_maturityDate = AddMonthsAdjustedForEndOfMonth
                          (referenceDate, m_nMaturityDuration * 12);
      break;
    }
  }

  // the maturity date
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }

  /// compare two legs
  bool operator<(const MyCashRate& other) const
    { return m_maturityDate < other.m_maturityDate; }
  /// compare two legs
  bool operator>(const MyCashRate& other) const
    { return m_maturityDate > other.m_maturityDate; }
  /// compare two legs
  bool operator<=(const MyCashRate& other) const
    { return m_maturityDate <= other.m_maturityDate; }
  /// compare two legs
  bool operator>=(const MyCashRate& other) const
    { return m_maturityDate >= other.m_maturityDate; }

private:

  Date m_maturityDate;

};

class MySwapRate : public finance::SwapRate
{
public:
  MySwapRate( const SwapRate& swapRate, Date referenceDate)
      : SwapRate(swapRate)
  {
    switch (m_maturityUnit)
    {
    case TimeUnit_Day:
      m_maturityDate = referenceDate.AddDays(m_nMaturityDuration);
      break;
    case TimeUnit_Month:
      m_maturityDate = AddMonthsAdjustedForEndOfMonth
                          (referenceDate, m_nMaturityDuration);
      break;
    case TimeUnit_Year:
      m_maturityDate = AddMonthsAdjustedForEndOfMonth
                          (referenceDate, m_nMaturityDuration * 12);
      break;
    }
  }

  MySwapRate( double dRate,
            size_t nMaturityDuration,
            TimeUnit maturityUnit,
            Frequency paymentFrequency,
            Date maturityDate
            )
        : SwapRate(dRate, nMaturityDuration, maturityUnit, paymentFrequency),
          m_maturityDate(maturityDate)
  {
  }

  // the maturity date
  Date GetMaturityDate() const
  {
    return m_maturityDate;
  }

  /// compare two legs
  bool operator<(const MySwapRate& other) const
    { return m_maturityDate < other.m_maturityDate; }
  /// compare two legs
  bool operator>(const MySwapRate& other) const
    { return m_maturityDate > other.m_maturityDate; }
  /// compare two legs
  bool operator<=(const MySwapRate& other) const
    { return m_maturityDate <= other.m_maturityDate; }
  /// compare two legs
  bool operator>=(const MySwapRate& other) const
    { return m_maturityDate >= other.m_maturityDate; }

private:

  Date m_maturityDate;

};

bool YieldCurveSwap::CheckData()
{
  CheckReferenceDate();

  // should have at least one cash rates or one swap rates.
  CHECK_COND( m_cashRates.GetAll().size() || m_swapRates.GetAll().size(),
              ITO33_YIELDCURVE_SWAP_NODATA);

  //_____________ fill in local CashRate array
  std::vector<MyCashRate> pCashRates;
  std::vector<MySwapRate> pSwapRates;
  size_t nRate = 0;

  for (nRate = 0; nRate < m_cashRates.GetAll().size(); nRate++)
  {
    pCashRates.push_back(
                  MyCashRate(m_cashRates.GetAll()[nRate], GetReferenceDate()));
  }

  for (nRate = 0; nRate < m_swapRates.GetAll().size(); nRate++)
  {
    pSwapRates.push_back(
                  MySwapRate(m_swapRates.GetAll()[nRate], GetReferenceDate()));
  }



  // we need the two following vector as the yieldcurve will not stop
  // being updated and validated, while YieldCurveAnnalCompounded::AddLeg()
  // doesn't work once yeild curve has been validated.
  std::vector<size_t>
    pnFinalTerms;
  std::vector<double>
    pdFinalRates;

  Date dateLast = GetReferenceDate();
  dateLast.AddDays(-1);

  //___________________________________________________________________________
  // treat cash rates
  // udpate pnFinalTerms and pdFinalRates
  if ( !pCashRates.empty() )
  {
    const std::vector<MyCashRate>::iterator
      pCashBegin = pCashRates.begin(),
      pCashEnd = pCashRates.end();

    std::sort(pCashBegin, pCashEnd);

    std::vector<MyCashRate>::iterator pCash;
    for (pCash = pCashRates.begin();
         pCash != pCashRates.end();
         pCash++)
    {
      if (dateLast >= pCash->GetMaturityDate())
      {
        throw EXCEPTION_MSG
          (
            ITO33_YIELDCURVE_SWAP_CASHRATES_MATURITES,
            ito33::String::Printf
              ( TRANS("%s %d %s"),
                  ITO33_YIELDCURVE_SWAP_CASHRATES_MATURITES.GetMessage(),
                  pCash->GetMaturityDuration(),
                  GetNameFromEnumValue( pCash->GetMaturityUnit(),
                                        SIZEOF(g_timeUnits),
                                        g_timeUnits )
              )
          );
      }

      dateLast = pCash->GetMaturityDate();

      size_t nDiffs = size_t(Date::DaysDiff(GetReferenceDate(), dateLast));
      double
        dFracLeg = Date::YearsDiff(GetReferenceDate(), dateLast,
                                   m_cashRates.GetBasis() ),
        dFracITO = nDiffs * ONEDAY;

      // (1+r)^ti = 1 + r*tl
      pnFinalTerms.push_back( nDiffs );
      pdFinalRates.push_back(pow( 1 + pCash->GetRate() * dFracLeg,
                                  1 / dFracITO ) - 1);

    }
  }
  /* end of treatment of cash rates ******************************************/
  /***************************************************************************/


  //___________________________________________________________________________
  // treat swap rates
  // udpate pnFinalTerms and pdFinalRates
  if ( !pSwapRates.empty() )
  {
    const std::vector<MySwapRate>::iterator
      iterBegin = pSwapRates.begin(),
      iterEnd = pSwapRates.end();

    std::sort(iterBegin, iterEnd);

    // now we can check for duplicates
    if ( std::adjacent_find(iterBegin, iterEnd, IsEqual) != iterEnd )
      throw EXCEPTION(ITO33_YIELDCURVE_SWAP_SWAPRATES_MATURITES);

    std::list<MySwapRate> pNewSwapRates;
    std::vector<MySwapRate>::iterator
      iterCurrent = pSwapRates.begin(),
      iterBefore = pSwapRates.end();

    CHECK_COND(iterCurrent->GetMaturityDate() >= dateLast,
               ITO33_YIELDCURVE_SWAP_INCONSISTENCY_SWAPCASH);

    if (iterCurrent->GetMaturityDate() > dateLast)
    {
      Date dateNewLeg = iterCurrent->GetMaturityDate();
      size_t nNumberOfMonth = iterCurrent->GetNumberMonths();
      size_t nFraction = 12 / iterCurrent->GetPaymentFrequency();

      pNewSwapRates.push_back(*iterCurrent);

      nNumberOfMonth -= nFraction;
      dateNewLeg = AddMonthsAdjustedForEndOfMonth
                        (GetReferenceDate(), nNumberOfMonth);

      while(dateNewLeg > dateLast)
      {
        // suppose that the swap rate is constant for the missing
        // terms at left hand
        MySwapRate newRate( iterCurrent->GetRate(),
                          nNumberOfMonth,
                          TimeUnit_Month,
                          iterCurrent->GetPaymentFrequency(),
                          dateNewLeg);

        pNewSwapRates.push_back(newRate);

        nNumberOfMonth -= nFraction;
        dateNewLeg = AddMonthsAdjustedForEndOfMonth
                        (GetReferenceDate(), nNumberOfMonth);
      }

      dateLast = iterCurrent->GetMaturityDate();
    }

    iterBefore = iterCurrent;
    iterCurrent++;

    // do interpolation for missing terms
    for(; iterCurrent != pSwapRates.end(); ++iterCurrent)
    {
      Date dateNewLeg;
      size_t
        nNumberOfMonth = iterCurrent->GetNumberMonths(),
        nFraction = 12 / iterCurrent->GetPaymentFrequency();

      double
        dRate,
        dRateBefore = iterBefore->GetRate(),
        dRateFraction = (iterCurrent->GetRate() - dRateBefore)
                      / Date::YearsDiff
                            (iterBefore->GetMaturityDate(),
                             iterCurrent->GetMaturityDate(),
                             m_swapRates.GetBasis());

      pNewSwapRates.push_back(*iterCurrent);

        // calculate the rate by interpolation
      nNumberOfMonth -= nFraction;
      dateNewLeg = AddMonthsAdjustedForEndOfMonth
                      (GetReferenceDate(), nNumberOfMonth);
      dRate = dRateBefore
            + dRateFraction * Date::YearsDiff(iterBefore->GetMaturityDate(),
                                              dateNewLeg,
                                              m_swapRates.GetBasis());

      while (dateNewLeg > dateLast)
      {
        MySwapRate newRate( dRate,
                            nNumberOfMonth,
                            TimeUnit_Month,
                            iterCurrent->GetPaymentFrequency(),
                            dateNewLeg);

        pNewSwapRates.push_back(newRate);
        nNumberOfMonth -= nFraction;
        dateNewLeg = AddMonthsAdjustedForEndOfMonth
                        (GetReferenceDate(), nNumberOfMonth);
        dRate = dRateBefore
              + dRateFraction * Date::YearsDiff(iterBefore->GetMaturityDate(),
                                                dateNewLeg,
                                                m_swapRates.GetBasis());
      }

      dateLast = iterCurrent->GetMaturityDate();
    }

    pNewSwapRates.sort();

    std::list<MySwapRate>::iterator
      pNewCurrent;

    SumDiscountFactor
      helperAnnual(Frequency_Annual),
      helperSemiAnnual(Frequency_SemiAnnual),
      helperQuarterly(Frequency_Quarterly),
      helperBiMonthly(Frequency_BiMonthly),
      helperMonthly(Frequency_Monthly);

    // go through all swap rates
    for (pNewCurrent = pNewSwapRates.begin();
         pNewCurrent != pNewSwapRates.end();
         ++pNewCurrent)
    {
      double
        dRate = pNewCurrent->GetRate()
              / double(pNewCurrent->GetPaymentFrequency());

      switch(pNewCurrent->GetPaymentFrequency())
      {
      case Frequency_Annual:
        helperAnnual.UpdateValues(pnFinalTerms, pdFinalRates,
                                     GetReferenceDate(),
                                     pNewCurrent->GetNumberMonths(),
                                     dRate);
        break;
      case Frequency_SemiAnnual:
        helperSemiAnnual.UpdateValues(pnFinalTerms, pdFinalRates,
                                         GetReferenceDate(),
                                         pNewCurrent->GetNumberMonths(),
                                         dRate);
        break;
      case Frequency_Quarterly:
        helperQuarterly.UpdateValues(pnFinalTerms, pdFinalRates,
                                        GetReferenceDate(),
                                        pNewCurrent->GetNumberMonths(),
                                        dRate);
        break;
      case Frequency_BiMonthly:
        helperBiMonthly.UpdateValues(pnFinalTerms, pdFinalRates,
                                        GetReferenceDate(),
                                        pNewCurrent->GetNumberMonths(),
                                        dRate);
        break;
      case Frequency_Monthly:
        helperMonthly.UpdateValues(pnFinalTerms, pdFinalRates,
                                      GetReferenceDate(),
                                      pNewCurrent->GetNumberMonths(),
                                      dRate);
        break;
      }
    }

  }

  // yieldCurveAnnuallyCompounded need at least 2 legs
  if (pnFinalTerms.size() == 1)
  {
    pnFinalTerms.push_back(pnFinalTerms[0]);
    pdFinalRates.push_back(pdFinalRates[0]);
  }

  // create yieldcurve
  m_pYieldCurve = make_ptr( new YieldCurveAnnuallyCompounded
                                (GetReferenceDate(), pnFinalTerms.size()) );
  m_pYieldCurve->SetLegs(pnFinalTerms, pdFinalRates);

  return true;
}


void YieldCurveSwap::DoValidate()
{
  CheckData();
}


void
YieldCurveSwap::GetZeroRates( const double *pdDays,
                              double *pdValues,
                              size_t nPoints) const
{
  Validate();

  m_pYieldCurve->GetZeroRates(pdDays, pdValues, nPoints);
}

void YieldCurveSwap::Visit(YieldCurveVisitor& visitor) const
{
  visitor.OnYieldCurveSwap(*this);
}

/// writes CashRate object in tagParent and return tag of the object
XML::Tag Dump(const CashRate& object, XML::Tag& tagParent)
{
  return Dump(object, XML_TAG_CASHRATE_ROOT, tagParent);
}

/// writes SwapRate object in tagParent and return tag of the object
XML::Tag Dump(const SwapRate& object, XML::Tag& tagParent)
{
  XML::Tag tagSwapRate(XML_TAG_SWAPRATE_ROOT, tagParent);

  DumpPart(object, tagSwapRate);

  tagSwapRate.Element(XML_TAG_PAYMENTFREQUENCY)
                     (XML::GetNameOfFrequency(object.GetPaymentFrequency()));
                  
  return tagSwapRate;
}

XML::Tag
YieldCurveSwap::Dump(const char *name, XML::Tag& tagParent) const
{
  XML::Tag tagYCVirtual(name, tagParent);
  XML::Tag tagYC(XML_TAG_YIELDCURVESWAP_ROOT, tagYCVirtual);

  tagYC.Element(XML_TAG_YIELDCURVE_REFERENCE_DATE)(GetReferenceDate());
  tagYC.Element(XML_TAG_YIELDCURVE_SWAP_CASHBASIS)
               (XML::GetNameOfDayCountConvention(m_cashRates.GetBasis()));
  tagYC.Element(XML_TAG_YIELDCURVE_SWAP_SWAPBASIS)
               (XML::GetNameOfDayCountConvention(m_swapRates.GetBasis()));

  size_t n;
  for (n = 0; n < m_cashRates.GetAll().size(); n++)
    finance::Dump(m_cashRates.GetAll()[n], tagYC);

  for (n = 0; n < m_swapRates.GetAll().size(); n++)
    finance::Dump(m_swapRates.GetAll()[n], tagYC);

  return tagYCVirtual;
}


shared_ptr<YieldCurveAnnuallyCompounded> YieldCurveSwap::BootStrap() const
{
  Validate();

  const std::vector<ZeroRate>& pLegs(m_pYieldCurve->GetAll());

  shared_ptr<YieldCurveAnnuallyCompounded>
    pyc(new YieldCurveAnnuallyCompounded(m_pYieldCurve->GetReferenceDate(),
                                         pLegs.size()));

  for(std::vector<ZeroRate>::const_iterator iter = pLegs.begin();
      iter != pLegs.end();
      ++iter)
    pyc->AddLeg(iter->GetDuration(), iter->GetRate());

  return pyc;
}


shared_ptr<YieldCurve>
YieldCurveSwap::Perturb(double dShift) const
{
  Validate();

  size_t n;

  CashRates cashRates;
  cashRates.SetBasis(m_cashRates.GetBasis());
  const std::vector<CashRate>& pCashRates(m_cashRates.GetAll());
  for ( n = 0; n < pCashRates.size(); ++n )
    cashRates.AddLeg( pCashRates[n].GetRate() + dShift,
                      pCashRates[n].GetMaturityDuration(),
                      pCashRates[n].GetMaturityUnit()
                    );


  SwapRates swapRates;
  swapRates.SetBasis(m_swapRates.GetBasis());
  const std::vector<SwapRate>& pSwapRates(m_swapRates.GetAll());
  for ( n = 0; n < pSwapRates.size(); ++n )
    swapRates.AddLeg( pSwapRates[n].GetRate() + dShift,
                      pSwapRates[n].GetMaturityDuration(),
                      pSwapRates[n].GetMaturityUnit(),
                      pSwapRates[n].GetPaymentFrequency()
                    );

  shared_ptr<YieldCurveSwap>
    pyc( new YieldCurveSwap(GetReferenceDate()) );
  pyc->SetSwapRates(swapRates);
  pyc->SetCashRates(cashRates);

  return pyc;
}

} // namespace finance

} // namespace ito33
