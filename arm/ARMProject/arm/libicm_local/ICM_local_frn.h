#ifndef ICM_LOCAL_FRN_H
#define ICM_LOCAL_FRN_H

#include "ICMKernel/glob/icm_enums.h"
class ARM_result;



extern long ICMLOCAL_FRN   (double	Spread,
					 double	Int_Accrual_DateIn,
					 double	MaturityIn,
					 long IrIndex,
					 double InitialRate,
					 double LastIndexFixing,
					 qPAYMENT_PREMIUM_LEG AccOnDef,
					 double	NotionalAmountIn ,
					 int		DayCount,
					 int		AccruedDayCount,
					 int		SettlementGap,	
					 CCString	DiscCurrency,
					 CCString   ResetCalendar,
					 CCString   PayCalendar,
					 double	First_Period_Reference_DateIn,
					 ARM_result& result,
					 long objId = -1);

extern long ICMLOCAL_CLN(const double& startdate,
				  const double& enddate,
				  const double& refdate,
				  const double& fstcpnrefdate,
				  const long&	IrIndex,
				  const double& spread,
				  const double&	Notional,
				  qPAYMENT_PREMIUM_LEG AccOnDef,
				  const int&	DayCount,
				  const int&	decompfreq,
				  const int&	stubrule,
				  const int&	resetgap,
				  CCString	DiscCurrency,
				  CCString   ResetCalendar,
				  CCString   PayCalendar,
				  const int&	nxchange,
				  const bool&   includematurity,
				  const int&	adjstartdate,
				  const int&	LegType,
				  const double& binary,
				  ARM_result& result,
				  long objId = -1);

#endif	// ARM_LOCAL_SWAP_H

/*----End Of File ----*/
// EOF %M%