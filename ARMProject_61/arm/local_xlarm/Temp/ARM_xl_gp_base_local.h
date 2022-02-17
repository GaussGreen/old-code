/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calculators_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARM_XL_GP_BASE_LOCAL_H
#define ARM_XL_GP_BASE_LOCAL_H


#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"
//////////////////////////////////
/// 
/// In order to factorise code
/// many function are calling
/// common function
/// these functions are not
/// declared here 
/// 
/// only exported functions are 
/// included here to facilitate 
/// the reading of the header
/// 
//////////////////////////////////



__declspec(dllexport) LPXLOPER WINAPI Local_DateStripCreate(LPXLOPER XL_startDate,
															LPXLOPER XL_endDate,
															LPXLOPER XL_resetFreq,
															LPXLOPER XL_dayCount,
															LPXLOPER XL_resetCalendar,
															LPXLOPER XL_fwdRule,
															LPXLOPER XL_intRule,
															LPXLOPER XL_stubRule,
															LPXLOPER XL_resetGap,
															LPXLOPER XL_payFreq,
															LPXLOPER XL_payGap,
															LPXLOPER XL_payCalendar,
															LPXLOPER XL_resetTiming,
															LPXLOPER XL_payTiming,
															LPXLOPER XL_adjFirstDate,
															LPXLOPER XL_refDate,
                                                            LPXLOPER XL_accruedFull,
															LPXLOPER XL_firstDateFwdRule);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripCreate(LPXLOPER XL_startDate,
																LPXLOPER XL_endDate,
																LPXLOPER XL_resetFreq,
																LPXLOPER XL_dayCount,
																LPXLOPER XL_resetCalendar,
																LPXLOPER XL_fwdRule,
																LPXLOPER XL_intRule,
																LPXLOPER XL_stubRule,
																LPXLOPER XL_resetGap,
																LPXLOPER XL_payFreq,
																LPXLOPER XL_payGap,
																LPXLOPER XL_payCalendar,
																LPXLOPER XL_resetTiming,
																LPXLOPER XL_payTiming,
																LPXLOPER XL_adjFirstDate,
																LPXLOPER XL_refDate,
                                                                LPXLOPER XL_accruedFull,
																LPXLOPER XL_firstDateFwdRule);

__declspec(dllexport) LPXLOPER WINAPI Local_DataFromDateStrip(LPXLOPER XL_DateStripId,
															  LPXLOPER XL_DataType );


__declspec(dllexport) LPXLOPER WINAPI Local_DateStripVecCreate(LPXLOPER XL_FlowStartDatesV,
															   LPXLOPER XL_FlowEndDatesV,
															   LPXLOPER XL_FwdStartDatesV,
															   LPXLOPER XL_FwdEndDatesV,
															   LPXLOPER XL_ResetDatesV,
															   LPXLOPER XL_PaymentDatesV,
															   LPXLOPER XL_InterestDaysV,
															   LPXLOPER XL_InterestTermsV );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DateStripVecCreate(LPXLOPER XL_FlowStartDatesV,
																   LPXLOPER XL_FlowEndDatesV,
																   LPXLOPER XL_FwdStartDatesV,
																   LPXLOPER XL_FwdEndDatesV,
																   LPXLOPER XL_ResetDatesV,
																   LPXLOPER XL_PaymentDatesV,
																   LPXLOPER XL_InterestDaysV,
																   LPXLOPER XL_InterestTermsV );


__declspec(dllexport) LPXLOPER WINAPI Local_Warning_SetPopUpFlag( LPXLOPER XL_Flag );
__declspec(dllexport) LPXLOPER WINAPI Local_Warning_Activate( LPXLOPER XL_Flag );


#endif