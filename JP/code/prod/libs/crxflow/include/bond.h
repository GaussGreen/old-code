/*********************************************************************************
 *    risky bond pricer
 *
 ********************************************************************************/

#ifndef __BOND_H__
#define __BOND_H__

#ifdef __cplusplus
extern "C"{
#endif

#include <alib/yield.h>
#include <alib/duration.h>

#include "crxutil.h"
#include "crcrv.h"
    
int RiskyBondPV_O(
    double              *result,          /* (O) bond price or duration         */
    TDate               valDate,          /* (I) value date                     */
    TDate               issueDate,        /* (I) issue date                     */
    TDate               maturityDate,     /* (I) maturity date                  */
    TDateInterval       cpnInterval,      /* (I) coupon interval                */
    TBoolean            stubAtEnd,        /* (I) T=Stub at end; F=Stub at beg.  */
    long                stubConv,         /* (I) stub conv                      */
    TDateInterval       delay,            /* (I) default payment delay          */
    long                DCC,              /* (I) payment daycount convention    */
    double              notional,         /* (I) notional                       */
    double              coupon,           /* (I) coupon rate                    */
    double              recovery,         /* (I) recovery                       */
    KAccrualConv        accrualConv,      /* (I) accrual conv                   */
    KStubType           priceConv,        /* (I) price conv                     */
    TCurve              *discCurve,       /* (I) ir curve                       */
    TCurve              *spdCurve,        /* (I) cr curve                       */
    char                choice);          /* (I) 'P': Price, 'D': Duration      */

#ifdef __cplusplus
}
#endif

#endif
