/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  ppbase.h
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Davis (Shuenn-Tyan) Lee
 *			   Robert Lenk
 *			   (Derivatives Research)
 *	Code version    :  1.10
 *	Extracted	:  7/12/96 at 18:41:12
 *	Last Updated	:  7/12/96 at 08:08:45
 ***************************************************************************
 *      Common utilities for prepay code
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#ifndef __ppbase_h
#define __ppbase_h

#include "ppconst.h"


/*****************************************************************************
 *   Public functions
 ****************************************************************************/

/************************************************************************
 *  AgencyOf
 *  Returns agency code (e.g., MBS_AGENCY_GNMAI) given
 *  the name of agency:
       'FNMA'  --> MBS_AGENCY_FNMA
       'FHLMC' --> MBS_AGENCY_FHLMC
       'GOLD'  --> MBS_AGENCY_GOLD
       'GNMAI' --> MBS_AGENCY_GNMAI
       'GNMAII'--> MBS_AGENCY_GNMAII
       'WHOLE' --> MBS_AGENCY_WHOLE
 ************************************************************************/
int AgencyOf(char *agencyString); /* (I) name of agency  */

/************************************************************************
 *  is_valid_prepay_form
 *  Boolean function to test if prepayForm is a valid type
 *  (one of MBS_PP_SPD_PSA, MBS_PP_SPD_CPR, MBS_PP_SPD_SMM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_prepay_form(long prepayForm);

/************************************************************************
 *  is_valid_prepay_model
 *  Boolean function to test if prepayModel is a valid type
 *  (one of  MBS_PP_MODEL_NONE, MBS_PP_MODEL_CONST, MBS_PP_MODEL_VECTOR,
 *  MBS_PP_MODEL_FIX_MRGP, MBS_PP_MODEL_ARM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_prepay_model(long prepayModel);

/************************************************************************
 *  is_valid_mbs_type
 *  Boolean function to test if mbsType is a valid type
 *  (one of  MBS_MBSTYPE_FIXED, MBS_MBSTYPE_ARM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_mbs_type(long mbsType);

/************************************************************************
 *  is_valid_agency_type
 *  Boolean function to test if agency is a valid type
 *  (one of  MBS_AGENCY_FNMA, MBS_AGENCY_FHLMC, MBS_AGENCY_GOLD,
 *  MBS_AGENCY_GNMAI, MBS_AGENCY_GNMAII, MBS_AGENCY_WHOLE)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_agency_type(long agency);

/************************************************************************
 *  is_valid_agency_type
 *  Boolean function to test if mbsTerm is a valid type
 *  (one of  MBS_PP_MBSTERM_30, MBS_PP_MBSTERM_15)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_mbs_term(long mbsTerm);

/************************************************************************
 *  is_valid_refi_index_rate_type
 *  Boolean function to test if refi_index_rate_type is a valid type
 *  (one of  MBS_PP_REFI_TYPE_FH30, MBS_PP_REFI_TYPE_FH15,
 *  MBS_PP_REFI_TYPE_CMT10)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_refi_index_rate_type(long refiIndexRateType);

/************************************************************************
 *  CheckCapFloorNumber
 *  Boolean function to test if floor/cap spread has wrong sign or
 *  life cap/floor has wrong number.
 *  RETURNS true/false
 ************************************************************************/
int
CheckCapFloorNumber
   (double floorSpread,
    double capSpread,
    double lifeFloor,
    double lifeCap);

/************************************************************************
 *  get_DOM_for_ARM
 *  Utility to get Day-of-month for cashflow for this TYPE of ARM
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int get_DOM_for_ARM
   (long  agency,       /* (I) agency type */
    long *cashflowDOM); /* (O) Day-of-month */

/************************************************************************
 *  check_const_months_interval
 *  Utility to check time series has constant month(s) apart in adjacent
 *  time points
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int check_const_months_interval
   (long    monthsApart, /* (I) agency type */
    long    numTimePts,  /* (I) Number of time points */
    TDate *timeSeries); /* (I) Array of time series */

/************************************************************************
 *  ConvertPrepays
 *  Utility to convert prepay rate from curr_form (cpr, smm, PSA) 
 *  to new_form (same choices)
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int
ConvertPrepays
   (long     currForm, /* (I) curr form: MBS_PP_SPD_PSA, MBS_PP_SPD_CPR,
                       * or MBS_PP_SPD_SMM */
    long     newForm,  /* (I) new form for prepay rate */
    double  currRate, /* (I) current prepay rate (in decimal, so that
                       * 6% CPR would be 0.06) */
    double  age,      /* (I) (if PSA involved) age, in months */
    double *newRate); /* (O) new prepay rate, in decimal */


/************************************************************************
 *  PointAdjustCommitRate()
 *  Utility to adjust an array commitment rates (w/varying points)
 *  to an array of adjusted commitment rates corresponding to
 *  a constant points value
 ************************************************************************/
int PointAdjustCommitRate
   (long      numFhCommit,      /* (I) # rates supplied */
    long      fhCommitRateType, /* (I) type of FHLMC commitment rate used;
                                * should be either MBS_PP_REFI_TYPE_FH30 
                                * or MBS_PP_REFI_TYPE_FH15 */
    double  *fhCommit,         /* (I) array of monthly commit rates */
    double  *fhCommitPts,      /* (I) array of corresponding (monthly)
                                *  points (in decimal) */
    double   desiredPoints,    /* (I) desired percentage points 
                                *  for adjusted rate */
    double   ioMultiplier,     /* (I) IO multiplier */ 
    double  *fhAdjCommit);     /* (O) adjusted array of commitment rates
                                *  corresponding to desired_points */

#endif
