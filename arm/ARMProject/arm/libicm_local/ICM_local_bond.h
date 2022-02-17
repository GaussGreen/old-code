#ifndef ICM_LOCAL_BOND_H
#define ICM_LOCAL_BOND_H

class ARM_result;
#include "ICMKernel/glob/icm_enums.h"


extern long ICMLOCAL_BOND   (double	CouponRateIn,
					  double	Int_Accrual_DateIn,
					  double	MaturityIn,
					  int		frequencyIn,
					  double	First_Period_Reference_DateIn ,
					  double	NotionalAmountIn ,
					  qPAYMENT_PREMIUM_LEG       AccOnDef,
 					  CCString	Currency,
					  CCString  PaymentCalendar,
					  int		DayCount,
					  int		AccruedDayCount,
					  int		SettlementGap,
					  double	RedemptionValue,	
					  ARM_result& result,
					  long objId = -1);


extern long ICMLOCAL_YTOPRICE (long bondId, double settlement, double yield, ARM_result& result);

extern long ICMLOCAL_PTOYIELD (long bondId, double settlement, double price, ARM_result& result);


#endif	// ARM_LOCAL_SWAP_H

/*----End Of File ----*/
// EOF %M%