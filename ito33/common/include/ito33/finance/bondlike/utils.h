/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/finance/bondlike/utils.h
// Purpose:     useful functions for using bondlike classes
// Author:      ZHANG Yunzhi
// Created:     2005/02/24
// RCS-ID:      $Id: utils.h,v 1.13 2006/08/19 19:11:50 wang Exp $
// Copyright:   (c) 2004- Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/finance/bondlike/utils.h
    @brief useful functions for using bondlike classes
 */

#ifndef _ITO33_FINANCE_BONDLIKE_UTILS_H_
#define _ITO33_FINANCE_BONDLIKE_UTILS_H_

#include "ito33/sharedptr.h"
#include "ito33/date.h"
#include "ito33/finance/frequency.h"

namespace ito33
{

namespace finance
{

class ITO33_DLLDECL BondTerms;
class ITO33_DLLDECL CallSchedule;
class ITO33_DLLDECL CashFlowStream;
class ITO33_DLLDECL SessionData;
class ITO33_DLLDECL PutSchedule;
class ITO33_DLLDECL CallSchedule;
class ITO33_DLLDECL Numeraire;


/**
    Computes Coupons rates of given BondTerms espcially when it
    has floating coupon.

    @param terms bond terms
    @param sessionData session data containing yield curve information
    @param pNumeraire currency of the bond

    @return fixed cash distribution in normal situation. if bond terms has
            floating coupons, returned cash flow stream is a guess of 
            expected coupons payment in the future.
 */
shared_ptr<CashFlowStream> 
ComputeCouponRates(const BondTerms& terms, const SessionData& sessionData,
                   const shared_ptr<Numeraire>& pNumeraire);

/**
    @brief Computes the call price at the given date.
    
    Throw an exception if there are no call at the given date.

    @param callDate the date at which the call price needs to be computed.
    @param bondTerms bond terms
    @param sessionData session data containing yield curve information
    @param calls the call schedule
    @param pNumeraire currency of the bond

    @return the call price at the given date if there are one at this date.
 */
double ComputeBondLikeCallPrice(Date callDate, 
                                const BondTerms& bondTerms,
                                const SessionData& sessionData,
                                const CallSchedule& calls,
                                const shared_ptr<Numeraire>& pNumeraire);

/**
    Gets the claim value at the given date.

    @param date the date at which the claim needs to be computed
    @param bondTerms bond terms
    @param sessionData session data containing yield curve information
    @param pNumeraire currency of the bond
  
    @return the claim value at the given date.
 */
double GetClaim(Date date, 
                const BondTerms& bondTerms,
                const SessionData& sessionData,
                const shared_ptr<Numeraire>& pNumeraire);

/**
    Gets the claim for a call or a put.
  
    @param date the date at which the claim needs to be computed.
    @param dYield the yield (must be positive)
    @param bondTerms bond terms
    @param sessionData session data containing yield curve information
    @param pNumeraire currency of the bond

    @return the claim value at the given date.
 */
double GetClaimFromYield(Date date, 
                         double dYield,
                         const BondTerms& bondTerms,
                         const SessionData& sessionData,
                         const shared_ptr<Numeraire>& pNumeraire);

/**
    @brief Gets yield compounding frequency from BondTerms object even when
           it is not explicitly specified.

    In fact, if the bond pays coupons, the function returns the payment 
    frequency.

    @param terms given bond terms object

    @return reasonable yield compounding frequency value.
            Otherwise Frequency_Undefined will be returned.
 */
Frequency GetYieldCompoundingFrequencyEvenIfUndefined(const BondTerms& terms);

/**
    Gets yield day count convention from BondTerms object even when
    it is not explicitly specified.

    In fact, if the bond pays coupons, the function returns the day count
    convention in the cash distribution.

    @param terms given bond terms object

    @return reasonable yield day count convention value.
            Otherwise DayCountConvention_Max will be returned.
 */

Date::DayCountConvention
GetYieldDayCountConventionEvenIfUndefined(const BondTerms& terms);

/**
    Checks whether yield compounding frequency has been well defined and throws
    exception if it is not the case.

    @param pPutSchedule gives yield to put information
    @param pCallSchedule gives yield to call information
    @param bondTerms gives compounding frequency information
 */
void CheckYieldCompoundingFrequency
     (
       const shared_ptr<PutSchedule>& pPutSchedule, 
       const shared_ptr<CallSchedule>& pCallSchedule,
       const BondTerms& bondTerms
     );

/**
    Checks whether yield day count convention has been well defined and throws
    exception if it is not the case.

    @param pPutSchedule gives yield to put information
    @param pCallSchedule gives yield to call information
    @param bondTerms gives yield day count convention information
 */
void CheckYieldDayCountConvention
     (
       const shared_ptr<PutSchedule>& pPutSchedule, 
       const shared_ptr<CallSchedule>& pCallSchedule,
       const BondTerms& bondTerms
     );

} // namespace finance

} // namespace ito33

#endif // #ifndef _ITO33_FINANCE_BONDLIKE_UTILS_H_
