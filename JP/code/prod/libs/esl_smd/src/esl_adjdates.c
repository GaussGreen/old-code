/****************************************************************************/
/*      Code for adjustable and flexible dates  .                           */
/****************************************************************************/
/*      ESL_ADJDATES_C                                                      */
/****************************************************************************/

/* 
$Id$
*/

#include "esl_adjdates.h"
#include "esl_paryield.h"


/*****************************************
 *
 *  COUPON STRUCTURE ADJUSTMENT FUNCTIONS
 *
 *****************************************/

/* if coupon is exposed in any way to past fixing logic, it is done in this function */
static int ADJ_AdjustCouponPastResets(
    ADJ_COUPON_DATA  *coupon_data,            /* (I/O) coupon structure            */
    T_CURVE          *zeroCurve,              /* (I) zero curve for reset          */
    long             *newResetEffDates,       /* (I) array of adjusted reset dates */
    long             valueDate,               /* (I) value date                    */
    int              offset,                  /* (I) adjustment required in days   */
    int              cpn,                     /* (I) current coupon number         */
    int              *pastAdjustmentRequired)  /* (O) TRUE/FALSE                   */
{
    int status = FAILURE;
    int dealPastReset = FALSE;
    int modelPastReset = FALSE;
    char *routine = "ADJ_AdjustPastResets";

    /* check if this current coupon has any exposure to past fixings */
    if (coupon_data->ResetEffDates[cpn] <= valueDate)
        dealPastReset = TRUE;
    if (newResetEffDates[cpn] <= valueDate)
        modelPastReset = TRUE;

    if (dealPastReset == TRUE && modelPastReset == TRUE)
    {
        /* if both deal and model reset dates are past fixings, nothing to do with past 
           fixings except change the past fixing date from the deal date to the model date
           if it is a multiple past fixing coupon type. */
        int j;

        *pastAdjustmentRequired = TRUE;

        if (coupon_data->PastFixingType == Multiple)
        {
            /* find the past fixing entry associated with this effective reset date and 
               move it for multiple past fixings */
            for (j = 0; j < *(coupon_data->NbPastFixings); j++)
            {
                if (coupon_data->PastResets.ResetEffDate[j] == coupon_data->ResetEffDates[cpn])
                {
                    coupon_data->PastResets.ResetEffDate[j] = newResetEffDates[cpn];
                    break;
                }
            }
        }
        coupon_data->YieldRatio[cpn] = 0.0;
    }
    else if (dealPastReset == TRUE && offset > 0) /* deal reset is before model reset */
    {
        int j;
        *pastAdjustmentRequired = TRUE;

        /* this handles deal reset <= valueDate < model reset */
        /* we must forward the entered past deal reset to the new model reset date, which
           is still in the future and hence a tree slice.  To account for this we set
           the ConstantReset field to the value of the entered deal past reset, which thus
           sets the reset slice as constant in the calc function */
        if (coupon_data->PastFixingType == Single)
        {
            coupon_data->ConstantReset[cpn] = *(coupon_data->PastResets.ResetRate);
            *(coupon_data->PastResets.ResetRate) = 0.0;  /* no longer past reset in model */
        }
        else  /* multiple */
        {
            int resetFound = FALSE;
            /* find the past fixing entry associated with this effective reset date */
            for (j = 0; j < *(coupon_data->NbPastFixings); j++)
            {
                if (coupon_data->PastResets.ResetEffDate[j] == coupon_data->ResetEffDates[cpn])
                {
                    coupon_data->ConstantReset[cpn] = coupon_data->PastResets.ResetRate[j];
                    resetFound = TRUE;
                    break;
                }
            }
            if (resetFound == FALSE)
            {
                DR_Error("%s unable to find entry for past reset effective date %ld",
                         routine,
                         coupon_data->ResetEffDates[cpn]);
                goto RETURN;
            }
            /* This is no longer a past reset from the model's point of view, but for simplicity
               leave it in the past reset array for the meantime as generally products don't 
               mind if there are unrequired past fixings entered, as long as they are on correct 
               reset dates */
            coupon_data->PastResets.ResetEffDate[j] = newResetEffDates[cpn];
        }
        coupon_data->YieldRatio[cpn] = 0.0;
    }
    else if (modelPastReset == TRUE && offset < 0)  /* model reset is before deal reset */
    {
        double YieldModel, Annuity;

        *pastAdjustmentRequired = TRUE;

        /* create a sythetic past reset by estimating the forward fixing at the deal
           reset date */
        if (Par_Yield(&YieldModel,
                      &Annuity,
                      zeroCurve->NbZero,
                      zeroCurve->Zero,
                      zeroCurve->ZeroDate,
                      zeroCurve->ValueDate,
                      coupon_data->ResetEffDates[cpn],
                      coupon_data->IndexMaturity,
                      coupon_data->IndexDayCount,
                      coupon_data->IndexFrequency) == FAILURE)
        {
            goto RETURN;
        }
            
        /* add this as a new past fixing to the collection of past fixings */
        if (coupon_data->PastFixingType == Single)
        {
            *(coupon_data->PastResets.ResetRate) = YieldModel;
        }
        else
        {
            /* if analytics support multiple past fixings we have to add an additional one
               to the list. */
            if (*(coupon_data->NbPastFixings) == 0)
            {
                *(coupon_data->NbPastFixings) = 1;
                coupon_data->PastResets.ResetEffDate[0] = newResetEffDates[cpn];
                /* leave the PastReset.ResetDate as is */
                coupon_data->PastResets.ResetRate[0] = YieldModel;
            }
            else
            {
                /* if >0 fixings already present we have to add this to the list in the
                   correct place as past fixings must be in ascending order */
                int j,k;
                /* if synthetic past reset is before all existing past resets */
                if (newResetEffDates[cpn] < coupon_data->PastResets.ResetEffDate[0])
                {
                    /* move all the existing fixings up one place */
                    for (k = *(coupon_data->NbPastFixings);  k > 0 ; k--)
                    {
                        coupon_data->PastResets.ResetEffDate[k] = coupon_data->PastResets.ResetEffDate[k-1];
                    }
                    /* now insert the new fixing in the correct place */
                    coupon_data->PastResets.ResetDate[0] = 0;  /* to signify synthetic fixing */
                    coupon_data->PastResets.ResetEffDate[0] = newResetEffDates[cpn];
                    coupon_data->PastResets.ResetRate[0] = YieldModel;
                    *(coupon_data->NbPastFixings)++;
                }
                else  
                {
                    /* synthetic reset is not first synthetic reset, so have to find the
                       appropriate place in the list to add it to */
                    for (j = *(coupon_data->NbPastFixings)-1; j >= 0; j--)
                    {
                        if (newResetEffDates[cpn] > coupon_data->PastResets.ResetEffDate[j])
                        {
                            /* move all the existing fixings up one place */
                            for (k = *(coupon_data->NbPastFixings); k > j+1; k--)
                            {
                                coupon_data->PastResets.ResetEffDate[k] = coupon_data->PastResets.ResetEffDate[k-1];
                                coupon_data->PastResets.ResetRate[k] = coupon_data->PastResets.ResetRate[k-1];
                            }
                            /* now insert the new fixing in the correct place */
                            coupon_data->PastResets.ResetDate[k] = 0;  /* to signify synthetic fixing */
                            coupon_data->PastResets.ResetEffDate[k] = newResetEffDates[cpn];
                            coupon_data->PastResets.ResetRate[k] = YieldModel;
                            break;
                        }
                    }
                    *(coupon_data->NbPastFixings) += 1;
                }
                /* double check past fixings are still in ascending order after inserting
                   the synthetic reset */
                for (j = 1; j < *(coupon_data->NbPastFixings); j++)
                {
                    if (coupon_data->PastResets.ResetEffDate[j] <= coupon_data->PastResets.ResetEffDate[j-1])
                    {
                        DR_Error("%s: Past Reset Dates (including new synthetic reset) not in ascending order",
                                 routine);
                        goto RETURN;               
                    }
                }
            }
        }
        coupon_data->YieldRatio[cpn] = 0.0;
    }

    status = SUCCESS;

    RETURN:

    return status;
}

int ADJ_AdjustCouponDates(
    ADJ_COUPON_DATA  *coupon_data,  /* (I/O) coupon structure     */
    long             valueDate,     /* (I) value date             */
    T_CURVE          *zeroCurve,    /* (I) zero curve for reset   */
    int              includeValueDateEvents,    /* (I) true/false */
    const char       *errorMsg)     /* (I) err msg help string    */
{
    int  status = FAILURE;
    int  i;

    int              adjustmentRequired = FALSE;
    char             *routine = "ADJ_AdjustCouponDates";
    ADJ_RESET_TYPE   resetType = UnDefined;
    int              pastAdjustmentRequired = FALSE;
    int              effectiveValueDate;

    long newResetEffDates[MAXNBDATE];

    for (i = 0; i < *(coupon_data->NbCoupons); i++)
    {
        coupon_data->YieldRatio   [i] = 1.0;
        coupon_data->ConstantReset[i] = 0.0;
    }

    /* Adjust accrual start/end dates to match the associated payment dates if required
       - Coupon 1 accrual start date must be assumed to be a good business day (ie. payment date)
         as principle exchanges/payments often occur on this day */

    /* Accrual start dates (apart from first coupon) must equal payment date of previous coupon */
    if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons)-1,
                         coupon_data->AccStartDates+1,
                         *(coupon_data->NbCoupons),
                         coupon_data->PaymentDates) == FAILURE)
    {
        for (i = 1; i < *(coupon_data->NbCoupons); i++)
        {
            coupon_data->AccStartDates[i] = coupon_data->PaymentDates[i-1];
        }
    }

    /* Accrual end dates must equal payment dates */
    if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons),
                         coupon_data->AccEndDates,
                         *(coupon_data->NbCoupons),
                         coupon_data->PaymentDates) == FAILURE)
    {
        for (i = 0; i < *(coupon_data->NbCoupons); i++)
        {
            coupon_data->AccEndDates[i] = coupon_data->PaymentDates[i];
        }
    }
    
    /* if all reset effective dates coincide with accrual start or end dates, nothing to do. */

    /* check if associated with accrual start dates - ie. fixing in advance */
    if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons),
                     coupon_data->ResetEffDates,
                     *(coupon_data->NbCoupons),
                     coupon_data->AccStartDates) == FAILURE)
    {
        /* check if associated with accrual end dates - ie. fixing in arrears */
        if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons),
                         coupon_data->ResetEffDates,
                         *(coupon_data->NbCoupons),
                         coupon_data->AccEndDates) == FAILURE)
        {
            adjustmentRequired = TRUE;
        }
        else
        {
            resetType = InArrears;
        }
    }
    else
    {
        resetType = InAdvance;
    }

    if (adjustmentRequired == FALSE)
    {
        /* simply set the adjustment factors to default, and return */
        for (i = 0; i < *(coupon_data->NbCoupons); i++)
        {
            coupon_data->YieldRatio[i] = 1.0;
            coupon_data->ConstantReset[i] = 0.0;
        }

        /* if reset offset type was undefined but the effective reset dates are
           aligned with either accrual start or end, then redefine this flag to
           reflect the status */
        if (coupon_data->ResetType == UnDefined)
        {
            coupon_data->ResetType = resetType;
        }

        status = SUCCESS;
        goto RETURN;
    }

    /********************************
      Reset Date Adjustment Required 
     ********************************/

    /* if we include the events on the value date reset on today will not be past fixings,
       so an effective valueDate for the past fixings logic used below has to be defined as
       it assumes a reset today is a past fixing */
    if (includeValueDateEvents == TRUE)
        effectiveValueDate = Nxtday(valueDate, -1);
    else
        effectiveValueDate = valueDate;

    if (coupon_data->ResetType == UnDefined)
    {
        DR_Error("%s, %s: Unable to determine whether resets are to be adjusted to inAdvance or inArrears",
                 routine, errorMsg);
        goto RETURN;
    }

    /* Adjust all the coupons to either in advance or arrears and calculate the YieldRatio */
    for (i = 0; i < *(coupon_data->NbCoupons); i++)
    {
        long   offset;

        pastAdjustmentRequired = FALSE;

        /* must adjust all coupons, even if they pay before the value date */
        if (coupon_data->ResetType == InArrears)
        {
            double Annuity;
            double YieldModel;
            double YieldDeal;

            /* if in arrears, we can only move forwards - any negative direction would mean the unadjusted 
               reset date > payment date, which is invalid */
            if (coupon_data->ResetEffDates[i] > coupon_data->PaymentDates[i])
            {
                DR_Error("%s. %s: When fixing in arrears, effective reset date number %d (%ld) must be <= "
                         "corresponding payment date %ld",
                         routine, errorMsg, i+1, coupon_data->ResetEffDates[i], coupon_data->PaymentDates[i]);
                goto RETURN;
            }
            /* calculate number of calendar days we have to move the reset by */
            offset = Daysact(coupon_data->ResetEffDates[i], coupon_data->PaymentDates[i]);
            
            newResetEffDates[i] = Nxtday(coupon_data->ResetEffDates[i], offset);

            /* take account of past resets */
            if (ADJ_AdjustCouponPastResets(coupon_data,
                                           zeroCurve,
                                           newResetEffDates,
                                           effectiveValueDate,
                                           offset,
                                           i,
                                           &pastAdjustmentRequired) == FAILURE)
            {
                goto RETURN;
            }

            /* if current coupon exposed to past fixing adjustment, then above
               function has done the job, so we move to the next coupon */
            if (pastAdjustmentRequired == TRUE)
                continue;

    
            /* useful information during debugging - small compute overhead */
            Par_Yield(&YieldModel,
                      &Annuity,
                      zeroCurve->NbZero,
                      zeroCurve->Zero,
                      zeroCurve->ZeroDate,
                      zeroCurve->ValueDate,
                      coupon_data->PaymentDates[i],
                      coupon_data->IndexMaturity,
                      coupon_data->IndexDayCount,
                      coupon_data->IndexFrequency);

            Par_Yield(&YieldDeal,
                      &Annuity,
                      zeroCurve->NbZero,
                      zeroCurve->Zero,
                      zeroCurve->ZeroDate,
                      zeroCurve->ValueDate,
                      coupon_data->ResetEffDates[i],
                      coupon_data->IndexMaturity,
                      coupon_data->IndexDayCount,
                      coupon_data->IndexFrequency);

            /* calculate the yield curve ratio from the deterministic t_curve */
            if (ParYieldRatio(&(coupon_data->YieldRatio[i]),
                                coupon_data->ResetEffDates[i],
                                coupon_data->PaymentDates[i],
                                0.0,  /* no spread */
                                zeroCurve->NbZero,
                                zeroCurve->ZeroDate,
                                zeroCurve->Zero,
                                zeroCurve->ValueDate,
                                coupon_data->IndexMaturity,
                                coupon_data->IndexDayCount,
                                coupon_data->IndexFrequency) == FAILURE)
            {
                goto RETURN;
            }
        }
        else  /* in advance */
        {
            double Annuity;
            double YieldModel;
            double YieldDeal;

            /* if in advance, move offset must be relative the payment date of the previous coupon */
            if (i == 0)                    
            {
                offset = Daysact(coupon_data->ResetEffDates[0], coupon_data->AccStartDates[0]);
            }
            else
            {
                offset = Daysact(coupon_data->ResetEffDates[i], coupon_data->PaymentDates[i-1]);
            }
            
            newResetEffDates[i] = Nxtday(coupon_data->ResetEffDates[i], offset);

            /* take account of past resets */
            if (ADJ_AdjustCouponPastResets(coupon_data,
                                           zeroCurve,
                                           newResetEffDates,
                                           effectiveValueDate,
                                           offset,
                                           i,
                                           &pastAdjustmentRequired) == FAILURE)
            {
                goto RETURN;
            }

            /* if current coupon exposed to past fixing adjustment, then above
               function has done the job, so we move to the next coupon */
            if (pastAdjustmentRequired == TRUE)
                continue;

            /* if reached here, the coupon has no exposure to the past reset logic, so we can
               simply calculate the yield adjustment factor */

            /* useful information during debugging - small compute overhead */
            Par_Yield(&YieldModel,
                      &Annuity,
                      zeroCurve->NbZero,
                      zeroCurve->Zero,
                      zeroCurve->ZeroDate,
                      zeroCurve->ValueDate,
                      coupon_data->AccStartDates[i],
                      coupon_data->IndexMaturity,
                      coupon_data->IndexDayCount,
                      coupon_data->IndexFrequency);

            Par_Yield(&YieldDeal,
                      &Annuity,
                      zeroCurve->NbZero,
                      zeroCurve->Zero,
                      zeroCurve->ZeroDate,
                      zeroCurve->ValueDate,
                      coupon_data->ResetEffDates[i],
                      coupon_data->IndexMaturity,
                      coupon_data->IndexDayCount,
                      coupon_data->IndexFrequency);

            /* calculate the yield curve ratio from the deterministic t_curve */
            if (ParYieldRatio(&(coupon_data->YieldRatio[i]),
                                coupon_data->ResetEffDates[i],
                                coupon_data->AccStartDates[i],
                                0.0,  /* no spread */
                                zeroCurve->NbZero,
                                zeroCurve->ZeroDate,
                                zeroCurve->Zero,
                                zeroCurve->ValueDate,
                                coupon_data->IndexMaturity,
                                coupon_data->IndexDayCount,
                                coupon_data->IndexFrequency) == FAILURE)
            {
                goto RETURN;
            }
        }
    }

        
    /* finally populate the new adjusted reset effective dates into the trade structure */
    for (i = 0; i < *(coupon_data->NbCoupons); i++)
    {
        coupon_data->ResetEffDates[i] = newResetEffDates[i];
        /* leave ResetDates as the original ones as they're not used in the analytics and it
           is essentially impossible for us to backward engineer what the calendar offset
           should be to reconstruct the "adjusted reset dates" */
    }

    /* double check all dates were adjusted correctly to coincide exactly with the accrual start or end dates */
    if (coupon_data->ResetType == InArrears)
    {
        if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons),
                         coupon_data->ResetEffDates,
                         *(coupon_data->NbCoupons),
                         coupon_data->AccEndDates) == FAILURE)
        {
            DR_Error("%s, %s: Adjusted reset in arrears dates do not coincide with the accrual end dates",
                     routine, errorMsg);
            goto RETURN;
        }
    }
    else
    {
        if (ADJ_DrDatesInSet(*(coupon_data->NbCoupons),
                         coupon_data->ResetEffDates,
                         *(coupon_data->NbCoupons),
                         coupon_data->AccStartDates) == FAILURE)
        {
            DR_Error("%s, %s: Adjusted reset in advance dates do not coincide with the accrual end dates",
                     routine, errorMsg);
            goto RETURN;
        }
    }


    status = SUCCESS;

    RETURN:
    
    return status;

}



/***********************
 *
 *  OPTION ADJUSTMENTS
 *
 ***********************/


int ADJ_AdjustOptionDates(
    ADJ_OPTION_DATA  *option_data,    /* (I/O) adapter structure           */
    long             valueDate)       /* (I) value date                    */
{
    int i;
    int status = FAILURE;
    char *routine = "ADJ_AdjustOptionDates";
    int removeExDateRow=-1;

    /* check if there are "past fixing" option notification dates (ie notif <= valueDate < Exercise)
       The middle office procedure should avoid this situation as analytics can not support this.
       We simply remove the offending option notif/exercise date from the list */
    for (i = 0; i < *(option_data->NbExercise); i++)
    {
        if (option_data->ExerDates[i] > valueDate &&
            option_data->NotifEffDates[i] <= valueDate)
        {
            removeExDateRow = i;
            /* print warning message,but do not return as error  */
            DR_Error("%s: Option Notification date (%ld) can not be on or before value date (%ld) "
                     "when exercise date (%ld) is after value date",
                     routine, option_data->NotifEffDates[i], valueDate, option_data->ExerDates[i]);
            break;
        }
    }

    /* if past fixing option notification, remove date from list and shuffle up remaining dates */
    if (removeExDateRow != -1)
    {
        for (i = 0; i < *(option_data->NbExercise); i++)
        {
            if (i >= removeExDateRow)
            {
                option_data->ExerDates[i] = option_data->ExerDates[i+1];
                option_data->NotifDates[i] = option_data->NotifDates[i+1];
                option_data->NotifEffDates[i] = option_data->NotifEffDates[i+1];
            }
        }
        *(option_data->NbExercise) -= 1;
    }

    /* for the meantime, if any call effective notification dates do not 
       coincide with the exercise date, simply move the notification date
       to the exercise date.  We are mispicing the call and applying no 
       adjustment factor at this stage */
    for (i = 0; i < *(option_data->NbExercise); i++)
    {
        /* this should have been checked in manager of product, but just in case */
        if (option_data->NotifEffDates[i] > option_data->ExerDates[i])
        {
            DR_Error("%s: Option effective notification date (%ld) must be <= exercise date (%ld) "
                     "Removing exercise date from the option schedule",
                     routine, option_data->NotifEffDates[i], option_data->ExerDates[i]);
            goto RETURN;
        }
        option_data->NotifEffDates[i] = option_data->ExerDates[i];
        /* leave the notifDates unchanged - they are not used in analysis and no point trying to
           determine the equivalent adjusted notifDate */
    }


    status = SUCCESS;

    RETURN:

    return status;
}




/****************************************************************************/
/*                                                                          */
/*   FUNCTION    ADJ_DrDatesInSet                                           */
/*                                                                          */
/*   "Silent" version of DRDatesInSet that doesn't write any error          */
/*   messages to DR_Error, simply returns FAILURE if the dates do not       */
/*   match   .                                                              */
/*                                                                          */
/*   DatesIn must be an arrays of ordered, strictly increasing dates.       */
/*                                                                          */
/****************************************************************************/

int ADJ_DrDatesInSet(
    int        NbDates1,      /* (I) Number of dates to check   */
    long      *Dates1,        /* (I) Dates to check             */
    int        NbDates2,      /* (I) Number of org dates        */ 
    long      *Dates2)        /* (I) Org dates                  */
{

    int     i, j;
    int     status = FAILURE;


    /* Pre-checks */
    if (NbDates1 == 0) return(SUCCESS); /* Nothing to check */

    /* Perform the comparison */
    i = 0;
    for (j = 0; j < NbDates2; j++)
    {
        if (Dates2[j] == Dates1[i]) i++;

        if (i>NbDates1-1) break;
    }

    if (i < NbDates1)
    {
        goto RETURN;
    }


    status = SUCCESS;

    RETURN:

    return(status);

} 
