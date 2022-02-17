/***************************************************************************
 *	SCCS Keyword Information
 *	------------------------
 *	Module name	:  ppbase.c
 *	Company Name	:  JP Morgan Securities Inc.
 *	Authors  	:  Davis (Shuenn-Tyan) Lee
 *			   (Derivatives Research)
 *	Code version    :  1.6
 *	Extracted	:  5/19/97 at 13:54:16
 *	Last Updated	:  5/19/97 at 13:54:05
 ***************************************************************************
 *      Common utilities for prepay code
 *
 *      Copyright 1996 J.P. Morgan & Co. Incorporated. All rights reserved.
 ***************************************************************************/
#include <math.h>
extern "C" {
#include <cdate.h>
#include <ldate.h>
#include <convert.h>
}
#include "mbsbase.h"

#include "ppconst.h"
#include "mbsconst.h"
#include "ppbase.h"

/* Convenience macro to eval CPR of 100% PSA (guards against bad age) */
#define CPR_OF_PSA100(a) ((0.06)*(MIN(MAX((a),1.e-6),30.)/30.))

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
int AgencyOf(char *agencyString) /* (I) name of agency  */
{
    F_INIT("AgencyOf");
    char  aStr[SHORTSTRLEN];


    XONTRUE( agencyString IS NULL, "No string supplied" );
    strncpy(aStr,agencyString, SHORTSTRLEN-1);
    XONFAIL( upcase(aStr) );

    if( (!strcmp(aStr,"FNMA")) )
    {
        return MBS_AGENCY_FNMA;
    }
    else if( (!strcmp(aStr,"FHLMC")) )
    {
        return MBS_AGENCY_FHLMC;
    }
    else if( (!strcmp(aStr,"GOLD")) || (!strcmp(aStr,"GOLDPC")) )
    {
        return MBS_AGENCY_GOLD;
    }
    else if( (!strcmp(aStr,"GNMA")) || (!strcmp(aStr,"GNMAI")) )
    {
        return MBS_AGENCY_GNMAI;
    }
    else if( (!strcmp(aStr,"GNMAII")) || (!strcmp(aStr,"GNMA2")) )
    {
        return MBS_AGENCY_GNMAII;
    }
    else if( (!strcmp(aStr,"WHOLE")) )
    {
        return MBS_AGENCY_WHOLE;
    }
    else
    {
        XONTRUE( TRUE, "Invalid agency name" );
    }
    
    F_END;
}


/************************************************************************
 *  WARNING: The input age (in months) must be positive!!!
 *  The input CPR is .06 for 6% CPR; The output PSA is 1 for 100% PSA.
 *  is_valid_prepay_form
 *  Boolean function to test if prepayForm is a valid type
 *  (one of MBS_PP_SPD_PSA, MBS_PP_SPD_CPR, MBS_PP_SPD_SMM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_prepay_form(long prepayForm)
{
    int status = FAILURE;
    static char routine[] = "is_valid_prepay_form";

    switch (prepayForm)
    {
    case MBS_PP_SPD_PSA:
    case MBS_PP_SPD_CPR:
    case MBS_PP_SPD_SMM:
         break;
    default: GtoErrMsg("%s: Invalid prepay form %d\n", routine,prepayForm);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  is_valid_prepay_model
 *  Boolean function to test if prepayModel is a valid type
 *  (one of  MBS_PP_MODEL_NONE, MBS_PP_MODEL_CONST, MBS_PP_MODEL_VECTOR,
 *  MBS_PP_MODEL_FIX_MGRP, MBS_PP_MODEL_ARM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_prepay_model(long prepayModel)
{
    int status = FAILURE;
    static char routine[] = "is_valid_prepay_model";

    switch (prepayModel)
    {
    case MBS_PP_MODEL_NONE:
    case MBS_PP_MODEL_CONST:
    case MBS_PP_MODEL_VECTOR:
    case MBS_PP_MODEL_FIX_MGRP:
    case MBS_PP_MODEL_ARM:
    case MBS_PP_MODEL_ARMSTATIC:
         break;
    default: GtoErrMsg("%s: Invalid prepay model %d\n", routine,prepayModel);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  is_valid_mbs_type
 *  Boolean function to test if mbsType is a valid type
 *  (one of  MBS_MBSTYPE_FIXED, MBS_MBSTYPE_ARM)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_mbs_type(long mbsType)
{
    int status = FAILURE;
    static char routine[] = "is_valid_mbs_type";

    switch (mbsType)
    {
    case MBS_MBSTYPE_FIXED:
    case MBS_MBSTYPE_ARM:
         break;
    default: GtoErrMsg("%s: Invalid mbs type %d\n", routine,mbsType);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  is_valid_agency_type
 *  Boolean function to test if agency is a valid type
 *  (one of  MBS_AGENCY_FNMA, MBS_AGENCY_FHLMC, MBS_AGENCY_GOLD,
 *  MBS_AGENCY_GNMAI, MBS_AGENCY_GNMAII, MBS_AGENCY_WHOLE)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_agency_type(long agency)
{
    int status = FAILURE;
    static char routine[] = "is_valid_agency_type";

    switch (agency)
    {
    case MBS_AGENCY_FNMA:
    case MBS_AGENCY_FHLMC:
    case MBS_AGENCY_GOLD:
    case MBS_AGENCY_GNMAI:
    case MBS_AGENCY_GNMAII:
    case MBS_AGENCY_WHOLE:
         break;
    default: GtoErrMsg("%s: Invalid agency type %d\n", routine,agency);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  is_valid_agency_type
 *  Boolean function to test if mbs term is a valid type
 *  (one of  MBS_PP_MBSTERM_30, MBS_PP_MBSTERM_15)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_mbs_term(long mbsTerm)
{
    int status = FAILURE;
    static char routine[] = "is_valid_mbs_term";

    switch (mbsTerm)
    {
    case MBS_PP_MBSTERM_15:
    case MBS_PP_MBSTERM_30:
         break;
    default: GtoErrMsg("%s: Invalid mbs term %d\n", routine,mbsTerm);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  is_valid_refi_index_rate_type
 *  Boolean function to test if refi_index_rate_type is a valid type
 *  (one of  MBS_PP_REFI_TYPE_FH30, MBS_PP_REFI_TYPE_FH15,
 *  MBS_PP_REFI_TYPE_CMT10)
 *  RETURNS true/false
 ************************************************************************/
int is_valid_refi_index_rate_type(long refiIndexRateType)
{
    int status = FAILURE;
    static char routine[] = "is_valid_refi_index_rate_type";

    switch (refiIndexRateType)
    {
    case MBS_PP_REFI_TYPE_FH30:
    case MBS_PP_REFI_TYPE_FH15:
    case MBS_PP_REFI_TYPE_CMT10:
         break;
    default: GtoErrMsg("%s: Invalid refi index rate type %d\n",
                 routine,refiIndexRateType);  
         goto done;  
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

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
    double lifeCap)
{
    int status = FAILURE;
    static char routine[] = "CheckCapFloorNumber";

    if (floorSpread > 0.)
    {
        GtoErrMsg("%s: Periodic floor spread(%lf) should be negative\n",
                  routine, floorSpread);
        goto done;
    }
    if (capSpread < 0.)
    {
        GtoErrMsg("%s: Periodic floor spread(%lf) should be positive\n",
                  routine, capSpread);
        goto done;
    }
    if (lifeFloor < 0.)
    {
        GtoErrMsg("%s: Minimum life coupon (%f) < 0.\n",
                  routine, lifeFloor);
        goto done;
    }
    if (lifeCap < 0.)
    {
        GtoErrMsg("%s: Maximum life coupon (%f) < 0.\n",
                  routine, lifeCap);
        goto done;
 
    }
    if (lifeFloor >= lifeCap)
    {
        GtoErrMsg("%s: Min life coupon (%f) >= Max life coupon (%f).\n",
                  routine, lifeFloor, lifeCap);
        goto done;
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
        return FALSE;
    }
    else
    {
        return TRUE;
    }
}

/************************************************************************
 *  get_DOM_for_ARM
 *  Utility to get Day-of-month for cashflow for this TYPE of ARM
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int get_DOM_for_ARM
   (long  agency,       /* (I) agency type */
    long *cashflowDOM)  /* (O) Day-of-month */
{
    int status = FAILURE;
    static char routine[] = "get_DOM_for_ARM";

    switch (agency)
    {
    case MBS_AGENCY_GNMAII:
        (*cashflowDOM) = 20;
        break;
        /* For now, refuse to recognize other agency types */
    case MBS_AGENCY_FNMA:
    case MBS_AGENCY_FHLMC:
    case MBS_AGENCY_GOLD:
    case MBS_AGENCY_GNMAI:
    case MBS_AGENCY_WHOLE:
    default:
        GtoErrMsg("%s: For now, cannot price ARMs for this agency\n",
                  routine);
        goto done;
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}

/************************************************************************
 *  check_const_months_interval
 *  Utility to check time series has constant month(s) apart in adjacent
 *  time points
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int check_const_months_interval
   (long    monthsApart, /* (I) agency type */
    long    numTimePts,  /* (I) Number of time points */
    TDate *timeSeries)  /* (I) Array of time series */
{
    int status = FAILURE;
    static char routine[] = "check_const_months_interval";
    long i,
        numOfMonths;
    TMonthDayYear
        mdy1,
        mdy2;

    if (numTimePts < 2)
    {
        GtoErrMsg("%s: Time series is less than two points\n", routine);
        goto done;
    }
    for (i = 1; i < numTimePts; i++)
    {
        GtoDateToMDY(timeSeries[i-1],&mdy1);
        GtoDateToMDY(timeSeries[i],&mdy2);

        numOfMonths = 12*(mdy2.year-mdy1.year)+(mdy2.month-mdy1.month);
        if (numOfMonths != monthsApart)
        {
            GtoErrMsg("%s: Time series is not constantly increased by [%d] months in %s and %s\n",
                routine,monthsApart,GtoFormatDate(timeSeries[i-1]),
                GtoFormatDate(timeSeries[i]));
            goto done;
        }
    }

    status = SUCCESS;
done:
    if (status IS FAILURE)
    {
        GtoErrMsg("%s : Failed\n",routine);
    }
    return(status);
}


/************************************************************************
 *  ConvertPrepays
 *  Utility to convert prepay rate from currForm (cpr, smm, PSA) 
 *  to newForm (same choices)
 *  Returns SUCCESS/FAILURE
 ************************************************************************/
int ConvertPrepays
   (long     currForm,  /* (I) curr form: MBS_PP_SPD_PSA, MBS_PP_SPD_CPR,
                        * or MBS_PP_SPD_SMM */
    long     newForm,   /* (I) new form for prepay rate */
    double  currRate,  /* (I) current prepay rate (in decimal, so that 6%
                        * CPR would be 0.06) */
    double  age,       /* (I) (if PSA involved) age, in months */
    double *newRate)   /* (O) new prepay rate, in decimal */
{
    F_INIT("ConvertPrepays");

    /* if no change needed, just quit */
    if( currForm == newForm )
	goto done;

    XONTRUE( !is_valid_prepay_form(currForm), "Bad currForm" );
    XONTRUE( !is_valid_prepay_form(newForm), "Bad newForm" );
    XONTRUE( currRate < 0., "Cannot have negative prepays" );
    if( currForm == MBS_PP_SPD_PSA ||
        newForm == MBS_PP_SPD_PSA )
    {
	XONTRUE( age <= 0., "Cannot convert to/from PSA with age<=0.0" );
    }

    switch (currForm) {
    case MBS_PP_SPD_PSA:  /* From PSA to ... */
        switch (newForm)
        {
        case MBS_PP_SPD_CPR:
             *newRate = currRate*CPR_OF_PSA100(age);
             break;
        case MBS_PP_SPD_SMM:
             *newRate = 1.-mbs_pow(1.-currRate*CPR_OF_PSA100(age),1./12.);
             break;
        default:
             XONTRUE( TRUE, "Fatal error 1" );
             break;
	}
	break;
    case MBS_PP_SPD_CPR:   /* From CPR to ... */
	switch (newForm)
        {
        case MBS_PP_SPD_PSA:
              *newRate = currRate/CPR_OF_PSA100(age);
              break;
        case MBS_PP_SPD_SMM:
             *newRate = 1.-mbs_pow(1.-currRate,1./12.);
             break;
        default:
             XONTRUE( TRUE, "Fatal error 2" );
             break;
        }
	break;
    case MBS_PP_SPD_SMM:   /* From SMM to ... */
	switch (newForm)
        {
        case MBS_PP_SPD_PSA:
             *newRate = (1.-mbs_pow(1.-currRate,12.))/CPR_OF_PSA100(age);
             break;
        case MBS_PP_SPD_CPR:
             *newRate = 1.-mbs_pow(1.-currRate,12.);
             break;
        default:
             XONTRUE( TRUE, "Fatal error 2" );
             break;
	}
	break;
    default:
	XONTRUE( TRUE, "Fatal error 0" );
	break;
    }


    F_END;
}

/************************************************************************
 *  PointAdjustCommitRate()
 *  Utility to adjust an array commitment rates (w/varying points)
 *  to an array of adjusted commitment rates corresponding to
 *  a constant points value
 ************************************************************************/
int PointAdjustCommitRate
   (long     numFhCommit,      /* (I) # rates supplied */
    long     fhCommitRateType, /* (I) type of FHLMC commitment rate used;
                               * should be either MBS_PP_REFI_TYPE_FH30 
                               * or MBS_PP_REFI_TYPE_FH15 */
    double *fhCommit,         /* (I) array of monthly commit rates */
    double *fhCommitPts,      /* (I) array of corresponding (monthly)
                               * points (in decimal) */
    double  desiredPoints,    /* (I) desired percentage points 
                               * for adjusted rate */
    double  ioMultiplier,     /* (I) IO multiplier */
    double *fhAdjCommit)      /* (O) adjusted array of commitment rates
                               * corresponding to desiredPoints */
{
    long irate;

    /* Make the adjustment */
    for (irate = 0; irate < numFhCommit; irate++)
    {
        fhAdjCommit[irate] = fhCommit[irate] -
            ( (desiredPoints - fhCommitPts[irate]) / ioMultiplier );
    }
    return(SUCCESS);
}
