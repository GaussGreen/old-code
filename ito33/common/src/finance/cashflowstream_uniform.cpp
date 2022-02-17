/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/finance/cashflowstream_uniform.cpp
// Purpose:     implementation of CashFlowStreamUniform
// Author:      Pedro Ferreira
// Created:     Mar 25, 2004
// RCS-ID:      $Id: cashflowstream_uniform.cpp,v 1.44 2006/08/19 23:06:55 wang Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/dateutils.h"
#include "ito33/arraycheckers.h"
#include "ito33/error.h"

#include "ito33/finance/error.h"
#include "ito33/finance/cashflowstream_uniform.h"

#include "ito33/xml/write.h"

#include "ito33/xml/finance/common.h"
#include "ito33/xml/finance/cashflowstream_all.h"
#include "ito33/xml/finance/daycountconvention.h"
#include "ito33/xml/finance/frequency.h"

using namespace ito33;
using namespace ito33::finance;

extern const ito33::Error
  ITO33_INVALID_DAYCOUNTCONVENTION;

extern const ito33::finance::Error
  ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE,
  ITO33_CASHFLOWSTREAM_FIRST_PAYMENT_DATE_AFTER_LAST,
  ITO33_CASHFLOWSTREAM_INCONSISTENT_DATES,
  ITO33_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT;

// Used to generate a period from a date
enum Direction
{
  Direction_Backward,
  Direction_Forward
};

/**
    Generates \a nNbPeriods regular periods from the date \a date in the
    direction \a direction.

    @param date the start date of the generation.
    @param freq the frequency of the payments.
    @param dcc the day count convention.
    @param nNbPeriods the number of periods needed.
    @param direction the direction of the generation.

    @return The found date.
 */
Date 
GenerateRegularPeriodsFrom(Date date, Frequency freq, 
                           Date::DayCountConvention dcc, size_t nNbPeriods, 
                           Direction direction)
{ 
  Date newDate;

  // regular period length in months is (12 / freq)
  int iInterval = 12 / freq * nNbPeriods;

  if (direction == Direction_Backward)
    iInterval = - iInterval;

  // End of month case
  if ( IsEOM(dcc) )
    newDate = AddMonthsAdjustedForEndOfMonth(date, iInterval);
  else // Non end of month case
    newDate = date.AddMonths(iInterval);

  return newDate;
}

CashFlowStreamUniform::CashFlowStreamUniform
    ( Date contractingDate, size_t nMonths, Frequency freq, double annualAmt,
      Date::DayCountConvention dcc)
    : CashFlowStream(contractingDate, dcc, freq),
      m_dAnnualPr(annualAmt),
      m_bAdjusted(false),
      m_lastPaymentType(LastPaymentType_Max),
      m_nMonths(nMonths)
{
  CHECK_COND( m_dAnnualPr > 0, ITO33_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT);

  // End of month case
  if ( IsEOM(m_dcc) )
  {
    m_firstPaymentDate = 
      AddMonthsAdjustedForEndOfMonth(contractingDate, 12 / freq);

    m_lastPaymentDate = 
      AddMonthsAdjustedForEndOfMonth(contractingDate, nMonths);
  }
  else // Non end of month case
  {
    m_firstPaymentDate = contractingDate.AddMonths(12 / freq);

    m_lastPaymentDate = contractingDate.AddMonths(nMonths);
  }

  Validate();

  Init();
} 

CashFlowStreamUniform::CashFlowStreamUniform
                       ( Date contractingDate,
                         const std::vector<Date>& paymentDates,
                         const std::vector<double>& paymentRates,
                         double annualAmt,
                         Date::DayCountConvention dcc,
                         Frequency freq )
                       : CashFlowStream(contractingDate, dcc, freq),
                         m_dAnnualPr(annualAmt),
                         m_bAdjusted(true),
                         m_lastPaymentType(LastPaymentType_Max),
                         m_nMonths(0)
{
  CheckIncreasingOrder(paymentDates, "Cash flow payment dates");

  CheckNonNegativity(paymentRates, "Cash flow payment rates");
 
  CHECK_COND
  (
    contractingDate < paymentDates[0],
    ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE
  );

  CHECK_COND( m_dAnnualPr > 0, ITO33_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT); 

  m_lastPaymentDate = paymentDates.back();
  m_firstPaymentDate = paymentDates.front();

  for (size_t n = 0; n < paymentDates.size(); n++)
    AddCashFlow(paymentDates[n], paymentRates[n]);
}

CashFlowStreamUniform::CashFlowStreamUniform( Date contractingDate,
                                              Date firstPaymentDate,
                                              Date lastPaymentDate,
                                              double annualAmt,
                                              Date::DayCountConvention dcc,
                                              Frequency freq,
                                              LastPaymentType lastPaymentType )
    : CashFlowStream(contractingDate, dcc, freq),
      m_lastPaymentDate(lastPaymentDate),
      m_firstPaymentDate(firstPaymentDate),
      m_dAnnualPr(annualAmt),
      m_bAdjusted(false),
      m_lastPaymentType(lastPaymentType),
      m_nMonths(0)
{
  finance::Validate(lastPaymentType);

  Validate();

  Init();
}
  
void CashFlowStreamUniform::Validate() const
{ 
  CHECK_COND( m_firstPaymentDate > m_contractingDate,
              ITO33_CASHFLOWSTREAM_BAD_FIRST_PAYMENT_DATE );

  CHECK_COND( m_firstPaymentDate <= m_lastPaymentDate,
              ITO33_CASHFLOWSTREAM_FIRST_PAYMENT_DATE_AFTER_LAST );

  CHECK_COND( m_dAnnualPr > 0, ITO33_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT); 
}

void CashFlowStreamUniform::Init()
{
  // value of each normal coupon (first and last one require extra
  // adjustments)
  const double valueCouponNormal = m_dAnnualPr / m_freq;

  Date date;  
  double dCoupon;

  // first coupon: we need to adjust the coupon value 
  dCoupon = GenerateCouponAmount(m_contractingDate, m_firstPaymentDate);

  AddCashFlow(m_firstPaymentDate, dCoupon);

  // Nothing more to generate if this is the only coupon
  if ( m_firstPaymentDate == m_lastPaymentDate )
    return;

  // regular coupon dates
  std::list<Date>
    regularDates( GenerateRegularDates(m_firstPaymentDate, m_lastPaymentDate, 
                                       12 / m_freq, IsEOM(m_dcc)) );

  // Check if last coupon is irregular
  bool bLastIrregular = false;
  if ( regularDates.empty() || regularDates.back() != m_lastPaymentDate )
  {
    CHECK_COND( IsValid(m_lastPaymentType),
                ITO33_CASHFLOWSTREAM_INCONSISTENT_DATES );
    
    bLastIrregular = true;
    if (  m_lastPaymentType == LastPaymentType_Long )
      if ( !regularDates.empty() )
        regularDates.erase(--regularDates.end());
  }

  std::list<Date>::const_iterator i;
  for ( i = regularDates.begin(); i != regularDates.end(); ++i )
    AddCashFlow(*i, valueCouponNormal);

  // Last irregular coupon
  if ( bLastIrregular )
  {
    Date lastButOneDate;
    if ( regularDates.empty() )
      lastButOneDate = m_firstPaymentDate;
    else
      lastButOneDate = regularDates.back();

    dCoupon = GenerateCouponAmount(lastButOneDate, m_lastPaymentDate);

    AddCashFlow(m_lastPaymentDate, dCoupon);
  }
}

bool CashFlowStreamUniform::IsRegularPeriod(Date startDate, Date endDate) const
{
  ASSERT_MSG( startDate < endDate, "endDate should be greater than startDate");

  if ( endDate == GenerateRegularPeriodsFrom(startDate, m_freq, m_dcc, 1, 
                                             Direction_Forward) )
    return true;
  else
    return false;
}

double CashFlowStreamUniform::GenerateCouponAmount(Date previousDate, 
                                                   Date currentDate) const
{
  ASSERT_MSG( previousDate < currentDate, 
              "currentDate should be greater than previousDate");
  
  ASSERT_MSG( currentDate == m_firstPaymentDate || 
              currentDate == m_lastPaymentDate, 
              "currentDate should be either the first or the last "
              "payment date");

  double dNormalCoupon = m_dAnnualPr / m_freq;

  if ( IsRegularPeriod(previousDate, currentDate) )
    return dNormalCoupon;

  // Irregular period
  double dCouponAmount;
  
  long 
    iNbDaysDiff = Date::DaysDiffWithDayCount(previousDate, currentDate, m_dcc);
  
  switch ( m_dcc ) 
  {
    case Date::DayCountConvention_30360:
    case Date::DayCountConvention_30E360:
    case Date::DayCountConvention_30U360:
    case Date::DayCountConvention_Act360:
    // Day count conventions with NO EOM (Non End Of Month)
    case Date::DayCountConvention_30360_NO_EOM:
    case Date::DayCountConvention_30E360_NO_EOM:
    case Date::DayCountConvention_30U360_NO_EOM:
    case Date::DayCountConvention_Act360_NO_EOM:
      {
        dCouponAmount = m_dAnnualPr * ( iNbDaysDiff / 360.);
      }
      break;
    
    case Date::DayCountConvention_Act365:
    case Date::DayCountConvention_Act365_NO_EOM:
      {
        dCouponAmount = m_dAnnualPr * ( iNbDaysDiff / 365.);
      }
      break;
    
    case Date::DayCountConvention_Act365L:
    case Date::DayCountConvention_Act365L_NO_EOM:
      {
        long iNbDays;

        double Y = 365.;

        if ( m_freq == Frequency_Annual )
        {
          // see footnote 7 page 9 of the documentation
          // accrued_interest_swx.pdf in mangue\ito33\docs\CashFlows
          iNbDays = Date::DaysDiffWithDayCount(previousDate, currentDate, 
                                               m_dcc);
          if (iNbDays == 365 || iNbDays == 366)
            Y = (double) iNbDays;
          else
          {
            Date dateTmp;

            long
              iIdx,
              iYear1 = previousDate.GetYear(),
              iYear2 = currentDate.GetYear();

            for (iIdx = iYear1; iIdx <= iYear2; ++iIdx)
            {
              if ( Date::IsLeapYear(iIdx) )
              {
                dateTmp.Set(iIdx, Date::Feb, 29);

                if (previousDate < dateTmp && dateTmp <= currentDate)
                  Y = 366.;
              }
            }
          }          
        }
        else
        {
          if ( currentDate.IsLeap() )
            Y = 366.;
        }

        dCouponAmount = m_dAnnualPr * (iNbDaysDiff / Y);
      }
      break;

    case Date::DayCountConvention_ActAct:
    case Date::DayCountConvention_ActAct_NO_EOM:
      {
        Direction direction;
        
        Date
          beginP,
          endP;

        dCouponAmount = 0.;
        if (currentDate == m_firstPaymentDate)
        {
          direction = Direction_Backward;
          endP = m_firstPaymentDate;
          beginP = GenerateRegularPeriodsFrom(endP, m_freq, m_dcc, 1, 
                                              direction);
          // create notional periods.
          while (beginP > m_contractingDate)
          {
            dCouponAmount += dNormalCoupon;
            endP = beginP;
            beginP = GenerateRegularPeriodsFrom(endP, m_freq, m_dcc, 1, 
                                                direction);
          }
          dCouponAmount += dNormalCoupon 
            * Date::DaysDiffWithDayCount(m_contractingDate, endP, m_dcc)
            / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);
        }
        else // current date is the last payment date
        {
          direction = Direction_Forward;
          beginP = previousDate; // Which is here the last but one payment date
          endP = GenerateRegularPeriodsFrom(beginP, m_freq, m_dcc, 1, 
                                            direction);
          // create notional periods.
          while (endP < m_lastPaymentDate)
          {
            dCouponAmount += dNormalCoupon;
            beginP = endP;
            endP = GenerateRegularPeriodsFrom(beginP, m_freq, m_dcc, 1, 
                                              direction);
          }
          dCouponAmount += dNormalCoupon 
            * Date::DaysDiffWithDayCount(beginP, m_lastPaymentDate, m_dcc)
            / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);
        }
      }
      break;

    default:
      throw EXCEPTION(ITO33_INVALID_DAYCOUNTCONVENTION);
  }

  return dCouponAmount;
}

double CashFlowStreamUniform::GetAccrued(Date date) const
{
  if ( date <= GetContractingDate() || date >= GetLastPaymentDate() )
    return 0;
  
  if (m_bAdjusted) // This is like the CashFlowStreamGeneral case
    return CashFlowStream::GetAccrued(date);

  double dAccrued;
  
  switch ( m_dcc ) 
  {
    case Date::DayCountConvention_30360:
    case Date::DayCountConvention_30E360:
    case Date::DayCountConvention_30U360:
    case Date::DayCountConvention_Act360:
    case Date::DayCountConvention_Act365:
    case Date::DayCountConvention_Act365L:
    // Day count conventions with NO EOM (Non End Of Month)
    case Date::DayCountConvention_30360_NO_EOM:
    case Date::DayCountConvention_30E360_NO_EOM:
    case Date::DayCountConvention_30U360_NO_EOM:
    case Date::DayCountConvention_Act360_NO_EOM:
    case Date::DayCountConvention_Act365_NO_EOM:
    case Date::DayCountConvention_Act365L_NO_EOM:
      {
        dAccrued = CashFlowStream::GetAccrued(date);
      }
      break;
    
    case Date::DayCountConvention_ActAct:
    case Date::DayCountConvention_ActAct_NO_EOM:
      {         
        double dNormalCoupon = m_dAnnualPr / m_freq;

        Direction direction;
        
        Date
          beginP,
          endP;        

        Elements::const_iterator iter;

        beginP = m_contractingDate;
        for (iter = m_elements.begin(); 
             iter != m_elements.end() && iter->first < date; 
             ++iter)
             beginP = iter->first;

        dAccrued = 0.;
        
        if (iter == m_elements.begin())
        {
          // i.e: contracting date < date <= first payment date

          direction = Direction_Backward;
          endP = m_firstPaymentDate;
          beginP = GenerateRegularPeriodsFrom(endP, m_freq, m_dcc, 1, 
                                              direction);
          // create notional periods.
          while (beginP > date)
          {
            endP = beginP;
            beginP = GenerateRegularPeriodsFrom(endP, m_freq, m_dcc, 1, 
                                                direction);
          }
          
          if (beginP <= m_contractingDate)
          {
            dAccrued += dNormalCoupon 
              * Date::DaysDiffWithDayCount(m_contractingDate, date, m_dcc)
              / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);
          }
          else // beginP > m_contractingDate
          {
            dAccrued += dNormalCoupon 
              * Date::DaysDiffWithDayCount(beginP, date, m_dcc)
              / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);

            while (beginP > m_contractingDate)
            {
              endP = beginP;
              beginP = GenerateRegularPeriodsFrom(endP, m_freq, m_dcc, 1, 
                                                  direction);
              if (beginP > m_contractingDate)
                dAccrued += dNormalCoupon;
            }
            dAccrued += dNormalCoupon 
              * Date::DaysDiffWithDayCount(m_contractingDate, endP, m_dcc)
              / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);     
          }
        }
        else if (iter == (--m_elements.end()))
        {
          // i.e: before laste payment date < date <= last payment date.
          // Here, beginP = last but one payment date
          direction = Direction_Forward;
          endP = GenerateRegularPeriodsFrom(beginP, m_freq, m_dcc, 1, 
                                            direction);
          // create notional periods.
          while (endP < date)
          {
            dAccrued += dNormalCoupon;
            beginP = endP;
            endP = GenerateRegularPeriodsFrom(beginP, m_freq, m_dcc, 1, 
                                              direction);
          }
          dAccrued += dNormalCoupon 
            * Date::DaysDiffWithDayCount(beginP, date, m_dcc)
            / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);
        }
        else
        {
          // i.e first payment date < date <= last but one payment date
          // Here, beginP < date <= iter->first
          endP = iter->first;
          dAccrued = dNormalCoupon 
            * Date::DaysDiffWithDayCount(beginP, date, m_dcc)
            / Date::DaysDiffWithDayCount(beginP, endP, m_dcc);
        }
      }
      break;

    default:
      throw EXCEPTION(ITO33_INVALID_DAYCOUNTCONVENTION);
  }

  return dAccrued;
}

shared_ptr<CashFlowStreamUniform>
CashFlowStreamUniform::ChangeSpread(double annualAmt) const
{
  if ( IsAdjusted() )
  {
    double dRatio = annualAmt / m_dAnnualPr;
    std::vector<Date> dates;
    std::vector<double> values;

    for (const_iterator i = begin(); i != end(); i++)
    {
      dates.push_back(i->first);
      values.push_back(i->second * dRatio);
    }

    return make_ptr( new CashFlowStreamUniform
                         (
                           GetContractingDate(),
                           dates,
                           values,
                           annualAmt,
                           GetDayCountConvention(),
                           GetPaymentFrequency()
                         ) );
  }

  return make_ptr( m_lastPaymentType != LastPaymentType_Max
                   ? new CashFlowStreamUniform
                         (
                           GetContractingDate(),
                           GetFirstPaymentDate(),
                           GetLastPaymentDate(),
                           annualAmt, // the only thing which changes
                           GetDayCountConvention(),
                           GetPaymentFrequency(),
                           m_lastPaymentType
                         )
                  : new CashFlowStreamUniform
                        (
                          GetContractingDate(),
                          GetFirstPaymentDate(),
                          GetLastPaymentDate(),
                          annualAmt, // the only thing which changes
                          GetDayCountConvention(),
                          GetPaymentFrequency()
                        ));
}

ito33::XML::Tag 
CashFlowStreamUniform::Dump(const char *name, ito33::XML::Tag& tagParent) const
{
  XML::Tag tagName(name, tagParent); 
  XML::Tag tagHere(XML_TAG_CASHFLOWSTREAMUNIFORM_ROOT, tagName);

  CashFlowStream::DumpMe(tagHere);

  if (m_bAdjusted)
  {
    Elements::const_iterator iter;

    for (iter = m_elements.begin(); iter != m_elements.end(); ++iter)
    {
      XML::Tag tagElement(XML_TAG_CASHFLOWSTREAMGENERAL_CASHFLOW, tagHere);

      tagElement.Element(XML_TAG_FINANCE_DATE)(iter->first);
      tagElement.Element(XML_TAG_FINANCE_RATE)(iter->second);
    }
  }
  else if (m_nMonths)
  {
    tagHere.Element(XML_TAG_CASHFLOWSTREAMUNIFORM_DURATION)(m_nMonths);
  }
  else
  {
    tagHere.Element(XML_TAG_CASHFLOWSTREAMUNIFORM_FIRSTDATE)
                   ( GetFirstPaymentDate() );

    if (m_lastPaymentType != LastPaymentType_Max)
    {
      tagHere.Element(XML_TAG_LASTPAYMENTTYPE)
                     (
                      GetNameFromEnumValue(
                        m_lastPaymentType,
                        SIZEOF(g_LastPaymentType),
                        g_LastPaymentType)
                     );
    }

    tagHere.Element(XML_TAG_CASHFLOWSTREAMUNIFORM_LASTDATE)
                  ( GetLastPaymentDate() );
  }

  tagHere.Element(XML_TAG_CASHFLOWSTREAMUNIFORM_ANNUALAMOUNT)
                 ( GetAnnualPaymentAmount() );
  
  return tagName;
}
