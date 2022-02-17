/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calib.h,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


/*! \file ARM_local_gp_calib.h
 *
 *  \brief file for addins of the generic calibration in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef ARMLOCAL_GP_CALIB_H
#define ARMLOCAL_GP_CALIB_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"


extern long ARMLOCAL_DateStrip_Create(double startDate,
									  double endDate,
									  long resetFreq,			
									  long dayCount,
									  const CCString& resetCalendar,
									  long fwdRule,
									  long longRule,
									  long stubRule,
									  long resetGap,
									  long payFreq,
									  long payGap,
									  const CCString& payCalendar,
									  long resetTiming,
									  long payTiming,
									  long adjFirstdate,
									  double refDate,
                                      long accruedfull,
									  long firstDateFwdRule,
									  ARM_result& result,                                      
									  long objId = -1 );

extern long ARMLOCAL_DateStripGetData(long dateStripId,
									  long dataType,
									  VECTOR<double>& Data,
									  ARM_result& result );


extern long ARMLOCAL_DateStripFromVec_Create(const VECTOR<double>&  FlowStartDatesV,
											 const VECTOR<double>&  FlowEndDatesV,
											 const VECTOR<double>&  FwdStartDatesV,
											 const VECTOR<double>&  FwdEndDatesV,
											 const VECTOR<double>&  ResetDatesV,
											 const VECTOR<double>&  PaymentDatesV,
											 const VECTOR<double>&  InterestDaysV,
											 const VECTOR<double>&  InterestTermsV,
											 ARM_result&	result,
											 long objId = -1 );


extern long ARMLOCAL_SetWarningFlag( 
	const bool val,
	ARM_result&	result );

extern long ARMLOCAL_Warning_Activate( 
	const bool val,
	ARM_result&	result );


#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
