

#ifndef ARM_XL_UTILITIES_H
#define ARM_XL_UTILITIES_H

#include <libCCxll\CCxll.h>

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_CQSO_Create(
	LPXLOPER XL_startDate,
	LPXLOPER XL_endDate,
	LPXLOPER XL_FundingCurrency,
	LPXLOPER XL_UnderlyingCurrency,
	LPXLOPER XL_1stIndex,
	LPXLOPER XL_2ndIndex,
	LPXLOPER XL_ResetDates,
	LPXLOPER XL_Notionals,
	LPXLOPER XL_Margins,
	LPXLOPER XL_Strikes,
	LPXLOPER XL_Leverages1,
	LPXLOPER XL_Leverages2,
	LPXLOPER XL_CpnMin,
	LPXLOPER XL_CpnMax,
	LPXLOPER XL_Fees,
	LPXLOPER XL_ScheduleArguments
	);

#endif