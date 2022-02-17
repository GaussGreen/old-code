/****************************************************************************/
/*      Construct the callable ladder time line.                            */
/****************************************************************************/
/*      TIME39.c                                                            */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fix123head.h"
#include "template39.h"

#ifndef LADDER
#define LADDER(f,bl,bh,d,m,u)  ( ((f)>(bh)) ? (u) : (((f)>(bl)) ? (m) : (d)) )
#endif

/*****  Ladder_Schedule ************************************************/
/**
 *       Sets up the time line. 
 *
 *       In this product, there are eleven 'sources' of critical dates:
 *
 *           0 - The exercise dates
 *                   (supp value: (strike - 1))
 *           1 - The ladder known pmt dates
 *                   (supp value: amt, rate, dcf, outs, factor )
 *     2,6,7,8 - The sticky reset dates
 *                   (supp value: lfloor, lcap, outs, dcf,
 *                                down, mid, up rate, lo, barrier hi,
 *                                wt1, wt2, sprd, stic wt,
 *                                lev, ifloor, icap
 *                    supp date : pmt date, acs date) 
 *           3 - The floating known pmt dates
 *                   (supp value: amt, rate, dcf, outs)
 *           4 - The floating reset dates
 *                   (supp value: floor, cap, stepup, outs, dcf,
 *                    supp date : pmt date, acs date) 
 *           5 - The state variable dates
 *                   (supp value: min, max)
 *           9 - The Rib Observation dates
 *           10- The corresponding complex payment dates
 *           11- Event Stat
 *
 *       This routine deals with the detailed date generation from input
 *       parameters and the subsequent addition of dates to the timeline
 *       critical date list.
 */
int     Ladder_Schedule 
            (long              ValueDate,         /**< (I) Value date         */
             T_CURVE          *t_curve,           /**< (I) Zero curve data    */
             MKTVOL_DATA      *mktvol_data,       /**< (I) volatility data    */
             LADDER_DATA      *ladder_data,       /**< (I) Structure of deal  */
             FIX3_TREE_DATA   *tree_data)         /**< (O) Tree data          */
{

    CRIT_DATE   *CritDate = NULL;        /* Critical date list             */
    int         NbCritDate = 0;          /* Number of critical dates       */

    /* Date lists to define zero banks */
    long       *ZeroMatDL[3] = {NULL,NULL,NULL};
    long       *ZeroUseDL[3] = {NULL,NULL,NULL};
    int         NbZeroMatDates[3] = {0,0,0};
    int         NbZeroUseDates[3] = {0,0,0};
    int         DCurve;

    /* Zero bank datelists for rate indices before optimisation */
    long       *IdxZMatDL[3] = {NULL,NULL,NULL};
    long       *IdxZUseDL[3] = {NULL,NULL,NULL};
    int         NbIdxZMatDates[3] = {0,0,0};
    int         NbIdxZUseDates[3] = {0,0,0};
    int         NbResetSt = 0;
    int         NbResetFl = 0;
    int         NbResetRib= 0;
    long       *ResetDLSt  = NULL;
    long       *ResetDLFl  = NULL;
    long       *ResetDLRib = NULL;
    int         IdxCurveSt [2] = {0, 0};
    int         IdxCurveRib[2] = {0, 0};
    int         IdxCurveFl  = 0;

    /* Exercise variables */
    EVENT_LIST  *ExerEventList = NULL;   /* Event list for exercises       */
    
    /* Ladder reset & pmt variables */
    EVENT_LIST  *StEventList = NULL;     /* Event list for fixed pmts      */
    int         FirstResetI = 0;         /* First reset with state var     */
    int         EndResetI = 0;           /* Last + 1 reset with state var  */
    int         StResetI = 0;            /* Reset date offset in EventList */ 
    int         stI = 0;                 /* Index for step up array        */
    int         payI = 0;                /* Index for step up array        */
    double      X;                       /* quantity to calc the strike    */
    double      CouponRate;              /* Rate of prossesed coupon       */
    double      CouponAmt;               /* Amount of prossesed coupon     */
    double      CouponFact;              /* Zero coupon factor for fixing  */
    double      MinState,  
                MaxState; 
    double      MinYield[2] = {0., 0.};  
    double      MaxYield[2] = {0., 0.};  /* bounds for yield at each reset */
    double      MinYieldTot, MaxYieldTot;/* bounds for the total yield     */
    double      RateLoObs[2], RateHiObs[2];/* variables for min,max  yields*/
    double      RateLoPmt[2],RateHiPmt[2];/* variables for min,max  yields*/
    double      StepN, StepX;            /* variables for min and max state*/
    double      InitOuts;                /* tmp var for cumulated notional */
    double      lb,hb,dr,mr,ur;          /* tmp var for barriers and rates */
    double      spd;                     /* tmp var for weight and spread  */
    double      wt[2] = {0., 0.};
    double      obswt[2] = {0., 0.};     /* tmp var for observation weight */
    double      stC;                     
    /* tmp var for sticky coef        */
    double      lev, idxcap, idxfloor;   /* tmp var for lev, icap & ifloor */
    long        MatDate[2]={0,0};        /* tmp var for calc of max var    */

    /* Floating reset & pmt variables */
    EVENT_LIST  *FlEventList = NULL;     /* Event list for float pmts      */
    int         suI = 0;                 /* Index for step up array        */

    /* Common variables */
    long        CurrResetDate;           /* Date of processed reset        */
    long        CurrCouponDate;          /* Date of prossesed coupon       */
    long        CurrExerciseDate;        /* Date of prossesed exercise     */
    int         NbFloatingPmts   = 0;    
    int         NbFloatingResets = 0;
    int         amI = 0;                 /* Index for amortization array   */
    int         i,j,event;               /* Index for convenience          */
    int         status = FAILURE;        /* Status = FAILURE initially     */
    double      Outstanding;             /* Outstanding notional           */
    double      DCF;                     /* Day count of prossesed coupon  */
    double      Fixing[2] ={0, 0};       /* Past fixing                    */
    double      FixingFl;

    int         s;                       /* Index for state variables      */
    int         RibIdx = 0;              /* Index for Rib observation index*/
    

    /* State variable MC simulation */
    IR_SIM      ir_sim[2];                  /* State variable path simulation */
    int         p;                       /* path */
    double      StepUp;
    int         NbPaths,
                NbEqSt;    /* Nb of equal spaced state variables between 
                            * the Max and Min */
    double      deltaS;
    double      YieldSlice[2][MAXNBSTATES];   /* Yield slice        */
    double      Yield[MAXNBSTATES];



    for(i = 0; i < 2; i++)
    {
        IdxCurveSt[i] = ladder_data->IdxIoDSt[i];
    }
    for(i = 0; i < 2; i++)
    {
        IdxCurveRib[i] = ladder_data->RibIdxIoD[i];
    }
    IdxCurveFl = ladder_data->IdxIoDFl;
    DCurve     = tree_data->CvDisc;


    /* Initialise critical date info in tree data */
    for (i=0; i<ZbkEVENT; i++)
    {
        tree_data->CritType[i] = 'D';
        tree_data->NbZeros[i]  =  0;
    }
 
    /* Allocate empty critical date list and then add dates successively. */
    CritDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (CritDate == NULL)
    {
        DR_Error("Ladder_Schedule: unable to allocate memory for critical "
                 "dates array !");
        goto FREE_MEM_AND_RETURN;
    }

    /* Always add value date to critical date list. */
    if (Add_To_DateList (&NbCritDate,
                         &CritDate,              
                         ValueDate,
                         NBCRITDATE,             /* No specific type */
                         0, 0, 0, 0, 0,              
                         0, 0, 0) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /********************************/
    /*  TYPE 0: EXERCISE DATES      */ 
    /********************************/
        
    tree_data->CritType[0] = 'D';    /* Discrete */
    tree_data->NbZeros[0]  = 0;

    /* Process the exercise input dates to generate the appropriate    */
    /* event list which will be a temporary EVENT_LIST structure.      */
    ExerEventList = DrNewEventListFromFreq 
                        (ladder_data->NbExer,
                         ladder_data->ExerDate,
                         ladder_data->Style,
                         'N',            /* Stub not allowed            */
                         'Y',            /* Input dates must be in list */
                         ladder_data->Strike,
                         NULL, NULL, NULL, NULL);
    
    if (ExerEventList == NULL)
    {
        goto FREE_MEM_AND_RETURN;
    }
           
    /* And finally add exercise dates to the critical date */
    /* list. Include strike as the supplementary values.   */
    amI  = 0;
    Outstanding = ladder_data->OrgNotPerc;
    for (i = 0; i < ExerEventList->NbEntries; i++)          
    {
        CurrExerciseDate = ExerEventList->Dates[i];

        /* Search for the amortization date falling before the current */
        /* date in the array of the amortization dates. Note that the  */
        /* outstanding notional on the amortization date should be     */
        /* decreased by that amortization amount.                      */
        while((amI < ladder_data->NbAmort) &&
              (CurrExerciseDate >= ladder_data->AmortDate[amI]))
        {
            Outstanding -= ladder_data->Amort[amI];
            amI ++;
        }

        if (ExerEventList->Dates[i] >= ValueDate)
        {
            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 ExerEventList->Dates[i],
                                 0,
                                 (ExerEventList->Curve[0][i] - 1.0)
                                 * ladder_data->NotionalSign * Outstanding, 
                                 0, 0, 0, 0, 
                                 0, 0, 0) == FAILURE)        
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
    }  /* for i */
    

    /*************************************/
    /*  TYPE 1: STICKY PMT DATES         */ 
    /*  TYPE 2: STICKY RESET DATES       */ 
    /*************************************/
        
    tree_data->CritType[1] = 'D';      /* Discrete */      
    tree_data->CritType[2] = 'D';      /* Discrete */      
    tree_data->NbZeros[1]  =  0;
    tree_data->NbZeros[2]  =  0;

    /* Process the swap inputs to generate the appropriate ladder coupon */
    /* payment dates which will be placed in a temporary EVENT_LIST.     */
    {
        long TempDates[2];             /* Only accrual start and final mat */

        TempDates[0] = ladder_data->AccStDate;
        TempDates[1] = ladder_data->MatDate;

        StEventList = DrNewEventListFromFreq 
                              (2,
                               TempDates,
                               ladder_data->PayFreqSt,
                               ladder_data->StubConv,
                               'N',    /* 'Dates in' check not required */
                               NULL, NULL, NULL, NULL, NULL);
        
        if (StEventList == NULL)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    /* Count floating pmt and set zerobank maturities when needed */
    NbFloatingPmts = 0;
    for (i = 1; i < StEventList->NbEntries; i++)            
    {
        if (StEventList->Dates[i] >= ValueDate)
        {
            /* For RIBs, discount factors from reset to payment
               are not needed. discount is done through tree DEV */
            if (ladder_data->ArrearsSt == 'N' &&
                !ladder_data->CplxIsRib)
            {
                /* Add pmt and reset pair to zerobank date list */
                if (AddDateToList(&(NbZeroMatDates[DCurve]),
                                  &(ZeroMatDL[DCurve]),
                                  StEventList->Dates[i]
                                 ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (AddDateToList(&(NbZeroUseDates[DCurve]),
                                  &(ZeroUseDL[DCurve]),
                                  StEventList->Dates[i-1]
                                 ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }
            NbFloatingPmts++; 
        } 
    }  

    /* Resets will either be on accrual starts or on payment    */
    /* dates of ladder so that the same event list can be used  */
     
    /* Add to the critical date list reset date and ladder info.     */
    /* In-arrears: no known pmt and no need for zeros to pmt.        */
    /* In-advnce:  possibly one known pmt and need for zeros to pmt. */
    /* Always calculate X for first reset: given or first lookback   */
    amI  = 0;
    stI = 0;
    Outstanding       = ladder_data->OrgNotPerc;
    X                 = ladder_data->FirstLevel;
    NbFloatingResets  = 0;
    ladder_data->FixingGivenSt = FALSE;
    InitOuts = ladder_data->Notional;

    if (ladder_data->ArrearsSt == 'Y')
    {
        EndResetI = StEventList->NbEntries;
        StResetI  = 1;
        FirstResetI = StResetI;

        ladder_data->FirstResetDateSt = StEventList->Dates[1];

        for (i = 1; i < StEventList->NbEntries; i++)            
        {
            CurrCouponDate = StEventList->Dates[i]; /* payment=reset */
            
            /* Search for ladder date falling before the current date */
            while((stI < ladder_data->NbStepUpSt - 1) &&
                  (CurrCouponDate > ladder_data->StDate[stI+1]))
            {
                stI++;
            }

            /* Search for amortization date falling before the current date */
            while((amI < ladder_data->NbAmort) &&
                  (CurrCouponDate > ladder_data->AmortDate[amI]))
            {
                Outstanding -= ladder_data->Amort[amI];
                amI++;
            }

            /* Day count fraction */
            if (DrDayCountFraction(StEventList->Dates[i-1],
                                   StEventList->Dates[i], 
                                   ladder_data->DayCountSt,
                                   &DCF) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
            if (StEventList->Dates[i] < ValueDate)
            {
                if (ladder_data->FixingDate[i-1] != StEventList->Dates[i])
                {
                    DR_Error("Ladder_Schedule: date of fixing != coupon!");
                    goto FREE_MEM_AND_RETURN;
                }

                /* calculate coupon rate for payment on Dates[i] */
                for(j = 0; j < 2; j++)
                {
                   Fixing[j] = ladder_data->FixingSt[j][i-1];
                }

                if( ladder_data->AoM == 'A')
                {
                    /* Binary */
                    X = ladder_data->StickyCoef[stI] * X 
                        + LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                 ladder_data->IdxObsWeightSt[1] * Fixing[1],
                                 ladder_data->BarrierLo[stI],
                                 ladder_data->BarrierHi[stI],
                                 ladder_data->DownRateSt[stI],
                                 ladder_data->MidRateSt[stI],
                                 ladder_data->UpRateSt[stI]);

                     /* Linear combination of index rate and spread */   
                    X += ladder_data->Leverage[stI] * 
                         COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] +
                                ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                                
                                ladder_data->SprdSt[stI],                                
                                ladder_data->IdxCap[stI],                                
                                ladder_data->IdxFloor[stI]);

                    X =  COLLAR(X,
                                ladder_data->CapSt[stI],                             
                                ladder_data->FloorSt[stI]);
                }
                else if( ladder_data->AoM == 'B')
                {
                    X = ladder_data->StickyCoef[stI] * X +
                        ( LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                 ladder_data->IdxObsWeightSt[1] * Fixing[1],
                                 ladder_data->BarrierLo[stI],
                                 ladder_data->BarrierHi[stI],
                                 ladder_data->DownRateSt[stI],
                                 ladder_data->MidRateSt[stI],
                                 ladder_data->UpRateSt[stI]) *
                          ladder_data->Leverage[stI] * 
                          COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] +
                                 ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                                
                                 ladder_data->SprdSt[stI],                                
                                 ladder_data->IdxCap[stI],                                
                                 ladder_data->IdxFloor[stI]) );

                    X =  COLLAR(X,
                                ladder_data->CapSt[stI],                             
                                ladder_data->FloorSt[stI]);
                }
                else
                {
                    X = ladder_data->StickyCoef[stI] * X +
                        ladder_data->Leverage[stI] * 
                        COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] + 
                               ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                               
                               ladder_data->SprdSt[stI],                               
                               ladder_data->IdxCap[stI],                               
                               ladder_data->IdxFloor[stI]);

                    X = COLLAR(X,
                               ladder_data->CapSt[stI],                            
                               ladder_data->FloorSt[stI]);

                    /* Linear combination of index rate and spread */                   
                    X *= LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                ladder_data->IdxObsWeightSt[1] * Fixing[1],                                
                                ladder_data->BarrierLo[stI],
                                ladder_data->BarrierHi[stI],
                                ladder_data->DownRateSt[stI],
                                ladder_data->MidRateSt[stI],
                                ladder_data->UpRateSt[stI]);
                }

                /* in Rib case, multiply by percentage of Rib obs in range */
                if (ladder_data->CplxIsRib)
                {
                    /* find last Rib observation in complex period */
                    while ((RibIdx < ladder_data->NbRibObsDates - 1) &&
                           (ladder_data->RibObsEffDate[RibIdx+1] < CurrCouponDate))
                    {
                        RibIdx++;
                    }

                    X *= ladder_data->RibInRangeWeight[RibIdx]  *
                         ladder_data->FixingRibPerc[i-1] +
                         ladder_data->RibOutRangeWeight[RibIdx] *
                         (1. - ladder_data->FixingRibPerc[i-1]);
                }
                
                /* accrue notional for zero coupon swap */
                if (ladder_data->SoZ=='Z') 
                    InitOuts *= (1 + ACC_FN(X,DCF,ladder_data->CompSt=='S'));
                
                /* first reset index */
                FirstResetI = i + 1;
            }
            else
            {
                /* Add to critical list */
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     2,
                                     ladder_data->FloorSt[stI],
                                     ladder_data->CapSt[stI],
                                     ladder_data->NotionalSign * Outstanding,
                                     DCF,
                                     0,
                                     StEventList->Dates[i],
                                     StEventList->Dates[i-1], 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     6,
                                     ladder_data->DownRateSt[stI],
                                     ladder_data->MidRateSt[stI],
                                     ladder_data->UpRateSt[stI],
                                     ladder_data->BarrierLo[stI],
                                     ladder_data->BarrierHi[stI],
                                     StEventList->Dates[i],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     7,
                                     ladder_data->IdxWeightSt[0][stI],
                                     ladder_data->IdxWeightSt[1][stI],
                                     ladder_data->SprdSt[stI],
                                     ladder_data->StickyCoef[stI],
                                     0,
                                     StEventList->Dates[i],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                 if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     8,
                                     ladder_data->Leverage[stI],
                                     ladder_data->IdxFloor[stI],
                                     ladder_data->IdxCap[stI],
                                     0,
                                     0,
                                     StEventList->Dates[i],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Prepare reset datelist for idx zero mats calculation */
                if (AddDateToList (&NbResetSt,
                                   &ResetDLSt,
                                   StEventList->Dates[i]
                                  ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                NbFloatingResets++;
            }  
        } 
    }
    else   /*Reset is NOT in arrears */
    {
        EndResetI = StEventList->NbEntries - 1;
        StResetI  = 0;
        FirstResetI = StResetI;

        ladder_data->FirstResetDateSt = StEventList->Dates[0];

        for (i = 0; i < StEventList->NbEntries-1; i++)          
        {
            CurrCouponDate = StEventList->Dates[i+1]; /* payment=reset+1*/
            
            /* Search for ladder date falling before the current date */
            while((stI < ladder_data->NbStepUpSt - 1) &&
                  (CurrCouponDate > ladder_data->StDate[stI+1]))
            {
                stI ++;
            }

            /* Search for amortization date falling before the current date */
            while((amI < ladder_data->NbAmort) &&
                  (CurrCouponDate > ladder_data->AmortDate[amI]))
            {
                Outstanding -= ladder_data->Amort[amI];
                amI ++;
            }

            /* Day count fraction */
            if (DrDayCountFraction(StEventList->Dates[i],
                                   StEventList->Dates[i+1], 
                                   ladder_data->DayCountSt,
                                   &DCF) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
            if (StEventList->Dates[i] < ValueDate)
            {
                if (ladder_data->FixingDate[i] != StEventList->Dates[i])
                {
                    DR_Error("Ladder_Schedule: date for fixing != coupon!");
                    goto FREE_MEM_AND_RETURN;
                }

                /* calculate coupon rate for payment on Dates[i+1] */
                for (j = 0; j < 2; j++)
                    Fixing[j] = ladder_data->FixingSt[j][i];

                if( ladder_data->AoM == 'A')
                {
                    /* Binary */
                    X = ladder_data->StickyCoef[stI] * X 
                        + LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                 ladder_data->IdxObsWeightSt[1] * Fixing[1],
                                 ladder_data->BarrierLo[stI],
                                 ladder_data->BarrierHi[stI],
                                 ladder_data->DownRateSt[stI],
                                 ladder_data->MidRateSt[stI],
                                 ladder_data->UpRateSt[stI]);

                     /* Linear combination of index rate and spread */
                     X += ladder_data->Leverage[stI] * 
                          COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] + 
                                 ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                                 
                                 ladder_data->SprdSt[stI],                                 
                                 ladder_data->IdxCap[stI],                                 
                                 ladder_data->IdxFloor[stI]);
                     X =  COLLAR(X,
                                 ladder_data->CapSt[stI],                            
                                 ladder_data->FloorSt[stI]);
                }
                else if( ladder_data->AoM == 'B')
                {
                    X = ladder_data->StickyCoef[stI] * X +
                        ( LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                 ladder_data->IdxObsWeightSt[1] * Fixing[1],
                                 ladder_data->BarrierLo[stI],
                                 ladder_data->BarrierHi[stI],
                                 ladder_data->DownRateSt[stI],
                                 ladder_data->MidRateSt[stI],
                                 ladder_data->UpRateSt[stI]) *
                          ladder_data->Leverage[stI] * 
                          COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] + 
                                 ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                                 
                                 ladder_data->SprdSt[stI],                                 
                                 ladder_data->IdxCap[stI],                                 
                                 ladder_data->IdxFloor[stI]) );
                     X =  COLLAR(X,
                                 ladder_data->CapSt[stI],                            
                                 ladder_data->FloorSt[stI]);
                }
                else
                {
                    X = ladder_data->StickyCoef[stI] * X +
                        ladder_data->Leverage[stI] * 
                        COLLAR(ladder_data->IdxWeightSt[0][stI] * Fixing[0] +
                               ladder_data->IdxWeightSt[1][stI] * Fixing[1] +                               
                               ladder_data->SprdSt[stI],                               
                               ladder_data->IdxCap[stI],                               
                               ladder_data->IdxFloor[stI]);
                    X = COLLAR(X,
                               ladder_data->CapSt[stI],                            
                               ladder_data->FloorSt[stI]);

                    X *= LADDER(ladder_data->IdxObsWeightSt[0] * Fixing[0] +
                                ladder_data->IdxObsWeightSt[1] * Fixing[1],
                                ladder_data->BarrierLo[stI],
                                ladder_data->BarrierHi[stI],
                                ladder_data->DownRateSt[stI],
                                ladder_data->MidRateSt[stI],
                                ladder_data->UpRateSt[stI]);
                }

                /* in Rib case, multiply by percentage of Rib obs in range */
                if (ladder_data->CplxIsRib && CurrCouponDate < ValueDate)
                {
                    /* find last Rib observation in complex period */
                    while ((RibIdx < ladder_data->NbRibObsDates - 1) &&
                           (ladder_data->RibObsEffDate[RibIdx+1] < CurrCouponDate))
                    {
                        RibIdx++;
                    }

                    X *= ladder_data->RibInRangeWeight[RibIdx]  *
                         ladder_data->FixingRibPerc[i] +
                         ladder_data->RibOutRangeWeight[RibIdx] *
                         (1. - ladder_data->FixingRibPerc[i]);
                }
                
                /* accrue notional for zero coupon swap */
                if (ladder_data->SoZ=='Z' &&  CurrCouponDate < ValueDate) 
                    InitOuts *= (1 + ACC_FN(X,DCF,ladder_data->CompSt=='S'));
                
                /* first reset index */
                FirstResetI = i + 1;

                /* add known pmt date and amount if date > value date */        
                if (StEventList->Dates[i+1] >= ValueDate) 
                {
                    ladder_data->FixingGivenSt = TRUE;

                    CouponAmt  = ACC_FN(X,DCF,ladder_data->CompSt=='S');
                    CouponFact = 1. + CouponAmt;
                    CouponAmt *= ladder_data->NotionalSign * Outstanding;

                    if (Add_To_DateList (&NbCritDate,
                                         &CritDate,
                                         StEventList->Dates[i+1],
                                         1,
                                         CouponAmt,
                                         X,
                                         DCF,
                                         ladder_data->NotionalSign * Outstanding,
                                         CouponFact, 
                                         0, 0, 0) == FAILURE)     
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                } /* if pay date >= value date */
            }
            else
            {
                /* Add to critical list */
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     2,
                                     ladder_data->FloorSt[stI],
                                     ladder_data->CapSt[stI],
                                     ladder_data->NotionalSign * Outstanding,
                                     DCF,
                                     0,
                                     StEventList->Dates[i+1],
                                     StEventList->Dates[i], 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     6,
                                     ladder_data->DownRateSt[stI],
                                     ladder_data->MidRateSt[stI],
                                     ladder_data->UpRateSt[stI],
                                     ladder_data->BarrierLo[stI],
                                     ladder_data->BarrierHi[stI],
                                     StEventList->Dates[i+1],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     7,
                                     ladder_data->IdxWeightSt[0][stI],
                                     ladder_data->IdxWeightSt[1][stI],
                                     ladder_data->SprdSt[stI],
                                     ladder_data->StickyCoef[stI],
                                     0,
                                     StEventList->Dates[i+1],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     StEventList->Dates[i],
                                     8,
                                     ladder_data->Leverage[stI],
                                     ladder_data->IdxFloor[stI],
                                     ladder_data->IdxCap[stI],
                                     0,
                                     0,
                                     StEventList->Dates[i+1],
                                     0, 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Prepare reset datelist for idx zero mats calculation */
                if (AddDateToList (&NbResetSt,
                                   &ResetDLSt,
                                   StEventList->Dates[i]
                                  ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                NbFloatingResets++;

            }
        } 
    }

    /* record bounds for X : MinState, MaxState  */        
    MinState = X;
    MaxState = X;
    ladder_data->InitState = X;
    ladder_data->InitOuts  = InitOuts;
    
    
    /********************************/
    /*  CHECK OF THE FIRST FIXING   */ 
    /********************************/

    if (ladder_data->FixingGivenSt)
    {
        if ( NbFloatingPmts != (NbFloatingResets+1) )
        {
            DR_Error("If 1st fixing is given, nb of ladder payments \n"
                     "must be exactly nb of resets plus one!");
            goto FREE_MEM_AND_RETURN;
        }
    }
    else
    {
        if (NbFloatingPmts != NbFloatingResets)
        {
            DR_Error("Input dates have resulted in nb of ladder payments \n"
                     "different from nb of floating resets!");
            goto FREE_MEM_AND_RETURN;
        }
    }  

    /**********************************/
    /*  TYPE 3: FLOATING PMT DATES    */ 
    /*  TYPE 4: FLOATING RESET DATES  */ 
    /**********************************/
        
    tree_data->CritType[3] = 'D';    /* Discrete */      
    tree_data->CritType[4] = 'D'; 
    tree_data->NbZeros[3]  =  0;
    tree_data->NbZeros[4]  =  0;
        
    /* Analogous to the payment dates for the ladder side. */
    {
        long TempDates[2];     /* Only accrual start and final mat */
        
        TempDates[0] = ladder_data->AccStDate;
        TempDates[1] = ladder_data->MatDate;
        
        FlEventList = DrNewEventListFromFreq
                              (2,     /* Nb of user dates */
                               TempDates,
                               ladder_data->PayFreqFl,
                               ladder_data->StubConv,
                               'N', /* 'Dates in' check not req */
                               NULL, NULL, NULL, NULL, NULL);

        if (FlEventList == NULL)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }
                                            
    /* Calculate number of flt pmt and set zerobank maturities when needed */
    NbFloatingPmts = 0;
    for (i = 1; i < FlEventList->NbEntries; i++)            
    {
        if (FlEventList->Dates[i] >= ValueDate)
        {
            if (ladder_data->ArrearsFl == 'N')
            {
                /* Add pmt and reset pair to zerobank date list */
                if (AddDateToList(&(NbZeroMatDates[DCurve]),
                                  &(ZeroMatDL[DCurve]),
                                  FlEventList->Dates[i]
                                 ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
                if (AddDateToList(&(NbZeroUseDates[DCurve]),
                                  &(ZeroUseDL[DCurve]),
                                  FlEventList->Dates[i-1]
                                 ) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }
            NbFloatingPmts++; 
        }
    }  

    /* Add to the critical date list reset date and info.            */
    /* In-arrears: no known pmt and no need for zeros to pmt.        */
    /* In-advnce:  possibly one known pmt and need for zeros to pmt. */
    amI  = 0;
    suI = 0;
    Outstanding = ladder_data->OrgNotPerc;
    NbFloatingResets = 0;

    if (ladder_data->ArrearsFl == 'Y')
    {
        ladder_data->FirstResetDateFl = FlEventList->Dates[1];

        for (i = 1; i < FlEventList->NbEntries; i++)            
        {
            CurrCouponDate = FlEventList->Dates[i]; /* payment=reset */
            
            /* Search for the amortization date falling before current date */
            while((amI < ladder_data->NbAmort) &&
                  (CurrCouponDate > ladder_data->AmortDate[amI]))
            {
                Outstanding -= ladder_data->Amort[amI];
                amI ++;
            }

            /* Search for stepup date falling before the current date */
            while((suI < ladder_data->NbStepUpFl - 1) &&
                  (CurrCouponDate > ladder_data->FlDate[suI+1]))
            {
                suI ++;
            }

            /* Day count fraction */
            if (DrDayCountFraction(FlEventList->Dates[i-1],
                                   FlEventList->Dates[i], 
                                   ladder_data->DayCountFl,
                                   &DCF) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }

            if (FlEventList->Dates[i] >= ValueDate)
            {
                /* Add to critical list */
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     FlEventList->Dates[i],
                                     4,
                                     ladder_data->FloorFl[suI],
                                     ladder_data->CapFl[suI],
                                     ladder_data->UpRateFl[suI],
                                     ladder_data->NotionalSign * Outstanding,
                                     DCF,
                                     FlEventList->Dates[i],
                                     FlEventList->Dates[i-1], 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Prepare reset datelist for idx zero mats calculation */
                if (ladder_data->IdxOnFl)
                {
                    if (AddDateToList (&NbResetFl,
                                       &ResetDLFl,
                                       FlEventList->Dates[i]
                                      ) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }

                NbFloatingResets++;
            }  
        } 
    }
    else   /*Reset is NOT in arrears */
    {
        ladder_data->FirstResetDateFl = FlEventList->Dates[0];

        for (i = 0; i < FlEventList->NbEntries-1; i++)          
        {
            CurrCouponDate = FlEventList->Dates[i+1]; /* payment=reset+1*/
            
            /* Search for the amortization date falling before current date */
            while((amI < ladder_data->NbAmort) &&
                  (CurrCouponDate > ladder_data->AmortDate[amI]))
            {
                Outstanding -= ladder_data->Amort[amI];
                amI ++;
            }

            /* Search for stepup date falling before the current date */
            while((suI < ladder_data->NbStepUpFl - 1) &&
                  (CurrCouponDate > ladder_data->FlDate[suI+1]))
            {
                suI ++;
            }

            /* Day count fraction including the abs(notional) */
            if (DrDayCountFraction(FlEventList->Dates[i],
                                   FlEventList->Dates[i+1], 
                                   ladder_data->DayCountFl,
                                   &DCF) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
                
            if (FlEventList->Dates[i] < ValueDate)
            {
                /* add known pmt date and amount if date > value date */    
                FixingFl= ladder_data->FixingFl;

                if (FlEventList->Dates[i+1] >= ValueDate)
                {
                    CouponRate = COLLAR(FixingFl * ladder_data->IdxWeightFl + 
                                        ladder_data->UpRateFl[suI],
                                        ladder_data->CapFl[suI],
                                        ladder_data->FloorFl[suI]);

                    CouponAmt  = ACC_FN(CouponRate,DCF,ladder_data->CompFl=='S');
                    CouponAmt *= ladder_data->NotionalSign * Outstanding;

                    if (Add_To_DateList (&NbCritDate,
                                         &CritDate,
                                         FlEventList->Dates[i+1],
                                         3,
                                         CouponAmt,
                                         CouponRate,
                                         DCF,
                                         ladder_data->NotionalSign * Outstanding,
                                         0, 
                                         0, 0, 0) == FAILURE)     
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                } /* if pay date >= value date */
            }
            else
            {
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     FlEventList->Dates[i],
                                     4,
                                     ladder_data->FloorFl[suI],
                                     ladder_data->CapFl[suI],
                                     ladder_data->UpRateFl[suI],
                                     ladder_data->NotionalSign * Outstanding,
                                     DCF,
                                     FlEventList->Dates[i+1],
                                     FlEventList->Dates[i], 0) == FAILURE)     
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Prepare reset datelist for idx zero mats calculation */
                if (ladder_data->IdxOnFl)
                {
                    if (AddDateToList (&NbResetFl,
                                       &ResetDLFl,
                                       FlEventList->Dates[i]
                                      ) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }
                NbFloatingResets++;
            }
        } 
    }
    
    /********************************/
    /*  CHECK OF THE FIRST FIXING   */ 
    /********************************/

    if (ladder_data->FixingGivenFl)
    {
        if ( NbFloatingPmts != (NbFloatingResets+1) )
        {
            DR_Error("If 1st fixing is specified, nb of floating payments \n"
                     "must be exactly nb of resets plus one!");
            goto FREE_MEM_AND_RETURN;
        }
    }
    else
    {
        if (NbFloatingPmts != NbFloatingResets)
        {
            DR_Error("Input dates have resulted in nb of floating payments \n"
                     "different from nb of floating resets!");
            goto FREE_MEM_AND_RETURN;
        }
    }  
    
    /*********************************/
    /*  TYPE 5: STATE VARIABLE DATES */ 
    /*********************************/
        
    tree_data->CritType[5] = 'D';       /* Discrete type */      
    tree_data->NbZeros[5]  =  0;

    MinYieldTot = MaxYieldTot = 0.;

    /* Bootstrap swaption volatility to get spot volatility curve */
    /* This will be used in Fix3_IndexVol().                           */
    if (Fix3_SpotVol (mktvol_data,
                      &t_curve[tree_data->CvDiff]) != SUCCESS)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Add state variable dates and bounds                              */
    /* state var date(i) = reset eff date(i)  (only those >= value date) */
    /* bounds(i) are calculated using bounds(i-1)                       */
    stI = 0;


    /* Dividing state variables into two groups:
     * 1) Equal spaced state varibales between the max and min of the cone
     * 2) MC simulated states - equal spaced in probability space
     * Require that NbStates >= 4
     *
     * GENERAL RULE:
     * 1) if the relationship between the state variable and payoff, e.g.
     * coupon vs. stickLeg, is LINEAR, then NbEqSt = NbStates, and NbPaths = 0
     * 2) if the relationship is NON-LINEAR, then allocate more MC generated
     * state variables will ensure proper sampling in the probability space.
     *
     * We choose half-half here to cover the non-linear 'C'ompounding 
     * vs. linear 'S'imple rate conventions (albeit small non-linear effect) 
     */
    if (ladder_data->CompSt == 'C' && !(ladder_data->CplxIsRib))
    {
        NbPaths = ladder_data->NbStates / 2;
    }
    else
    {
        NbPaths = 0;
    }

    NbEqSt = ladder_data->NbStates - NbPaths;


    /* Initialize the random variables and base date */ 
    for(j = 0; j < 2; j++)
    {
        if (Fix3_Initialize_IR_SIM (&ir_sim[j],
                           ValueDate,
                           NbPaths) == FAILURE)
        goto FREE_MEM_AND_RETURN;
        
    }
    ladder_data->State       = NULL;
    ladder_data->FirstResetI = FirstResetI;                                            
    ladder_data->EndResetI   = EndResetI;

    ladder_data->State = (double **) DR_Matrix (DOUBLE,
                                                FirstResetI-1, (EndResetI - 1L),
                                                0, ladder_data->NbStates - 1L);
    
    if (ladder_data->State == NULL)
    {
        DR_Error("Ladder_Schedule: unable to allocate memory for slices.");
        goto FREE_MEM_AND_RETURN;
    }
    else
    {
        /* Initialize state variables */
        for (p=0; p<ladder_data->NbStates; p++)
            ladder_data->State[FirstResetI-1][p] = ladder_data->InitState;
    }



    /* Add Min/Max states to the boundary.
     * They are estimated based on the NbStdev input,
     * but the simple algorithm here should be in general much more 
     * conservative than NbStdev specified for the state variable.
     */
    for (i = FirstResetI; i < EndResetI; i++)
    {
        /* add to critical date list values found in ladder reset */

        if (Add_To_DateList(&NbCritDate,
                            &CritDate,
                            StEventList->Dates[i],
                            5,
                            i,        /* reset index */
                            0,
                            0,
                            0, 0, 
                            0, 0, 0) == FAILURE)     
        {
            goto FREE_MEM_AND_RETURN;
        }


        /* calculate the avg vol of index from val date to reset */
        for (j = 0; j < 2; j++)
        {
            MatDate[j] = Nxtmth(StEventList->Dates[i],
                         ladder_data->IdxMatSt[j],
                         1L);

            /* MC simulation of yield distribution */
            if (Fix3_SwapYield_MC(
                     YieldSlice[j],
                     &MaxYield[j],
                     &MinYield[j],
                     ladder_data->NbStDev,
                     StEventList->Dates[i],
                     MatDate[j],
                     ladder_data->IdxFreqSt[j],
                     ladder_data->IdxBaseSt[j],
                     &ir_sim[j],
                     mktvol_data,
                     &t_curve[ladder_data->IdxIoDSt[j]]) != SUCCESS)
            {
               goto FREE_MEM_AND_RETURN;
            }
        }



        /* work out current ladder levels */
        /* find index for pay and reset dates and tenor */
        if (ladder_data->ArrearsSt == 'Y')
        {
            payI = i;
        }
        else
        {
            payI = i+1;
        }

        /* work out current loan levels and day count */
        CurrCouponDate = StEventList->Dates[payI]; 
        while((stI < ladder_data->NbStepUpSt - 1) &&
              (CurrCouponDate > ladder_data->StDate[stI+1]))
        {
            stI ++;
        }

        /* Current ladder levels */
        lb = ladder_data->BarrierLo[stI];
        hb = ladder_data->BarrierHi[stI];
        dr = ladder_data->DownRateSt[stI];
        mr = ladder_data->MidRateSt[stI];
        ur = ladder_data->UpRateSt[stI];
        for( j = 0; j < 2; j++)
        {
           wt[j] = ladder_data->IdxWeightSt[j][stI];
           obswt[j] = ladder_data->IdxObsWeightSt[j];
        }
        spd      = ladder_data->SprdSt[stI];
        stC      = ladder_data->StickyCoef[stI];
        lev      = ladder_data->Leverage[stI];
        idxcap   = ladder_data->IdxCap[stI];
        idxfloor = ladder_data->IdxFloor[stI];

        for (j = 0; j < 2; j++)
        {
            if ( obswt[j] > 0)
            {
                RateLoObs[j] = MinYield[j];
                RateHiObs[j] = MaxYield[j];
            }
            else
            {
                RateLoObs[j] = MaxYield[j];
                RateHiObs[j] = MinYield[j];
            }
        }

        MaxYieldTot = obswt[0] * RateHiObs[0] + obswt[1] * RateHiObs[1];                         
        MinYieldTot = obswt[0] * RateLoObs[0] + obswt[1] * RateLoObs[1];

        for (j = 0; j < 2; j++)
        {
            if ( wt[j] > 0)
            {
                RateLoPmt[j] = MinYield[j];
                RateHiPmt[j] = MaxYield[j];
            }
            else
            {
                RateLoPmt[j] = MaxYield[j];
                RateHiPmt[j] = MinYield[j];
            }
        }


        /* ------    
         * 1. Work out the boundary of the state variables in 
         * the most conservative scenario.
           ------*/

        if (MaxYieldTot <= lb)
        {
            StepN = StepX = dr; 
        }
        else if (MaxYieldTot <= hb)
        {
            if (MinYieldTot <= lb)
            {
                StepN = MIN(dr,mr); StepX = MAX(dr,mr);
            }
            else
            {
                StepN = StepX = mr; 
            }
        }
        else
        {
            if (MinYieldTot <= lb)
            {
                StepN = MIN(dr,MIN(mr,ur)); StepX = MAX(dr,MAX(mr,ur));
            }
            else if (MinYieldTot <= hb)
            {
                StepN = MIN(mr,ur); StepX = MAX(mr,ur);
            }
            else
            {
                StepN = StepX = ur; 
            }
        }



        if(ladder_data->AoM == 'A')
        {
            /* Binary */
            MinState = MIN(stC*MinState+StepN,stC*MinState+StepX);
            MaxState = MAX(stC*MaxState+StepN,stC*MaxState+StepX);

           /* Collared linear combination of index and spread */
            MinState += lev * COLLAR(wt[0] * RateLoPmt[0] + wt[1] * RateLoPmt[1] + spd,
                                     idxcap,
                                     idxfloor);            
            MaxState += lev * COLLAR(wt[0] * RateHiPmt[0] + wt[1] * RateHiPmt[1] + spd,
                                     idxcap,
                                     idxfloor);

            MinState = COLLAR(MinState,
                          ladder_data->CapSt[stI], 
                          ladder_data->FloorSt[stI]);
            MaxState = COLLAR(MaxState,
                          ladder_data->CapSt[stI], 
                          ladder_data->FloorSt[stI]);        
        }
        else if (ladder_data->AoM == 'B')
        {
            double MinProd, MaxProd, MinI, MaxI;

            MinI    = lev * COLLAR(wt[0] * RateLoPmt[0] + wt[1] * RateLoPmt[1] + spd,
                                   idxcap,
                                   idxfloor);
            MaxI    = lev * COLLAR(wt[0] * RateHiPmt[0] + wt[1] * RateHiPmt[1] + spd,
                                   idxcap,
                                   idxfloor);
            MinProd = MIN( MIN( MinI * StepN, MinI * StepX), 
                           MIN( MaxI * StepN, MaxI * StepX) );
            MaxProd = MAX( MAX( MinI * StepN, MinI * StepX), 
                           MAX( MaxI * StepN, MaxI * StepX) );

            MinState = stC*MinState + MinProd;
            MaxState = stC*MaxState + MaxProd;

            MinState = COLLAR(MinState,
                              ladder_data->CapSt[stI], 
                              ladder_data->FloorSt[stI]);
            MaxState = COLLAR(MaxState,
                              ladder_data->CapSt[stI], 
                              ladder_data->FloorSt[stI]);
        }
        else
        {
            double MinI, MaxI;

            MinI  = stC * MinState;
            MaxI  = stC * MaxState; 
             
            MinI += lev * COLLAR(wt[0] * RateLoPmt[0] + wt[1] * RateLoPmt[1] + spd,
                                 idxcap,
                                 idxfloor);            
            MaxI += lev * COLLAR(wt[0] * RateHiPmt[0] + wt[1] * RateHiPmt[1] + spd,
                                 idxcap,
                                 idxfloor);

            MinI = COLLAR(MinI,
                          ladder_data->CapSt[stI], 
                          ladder_data->FloorSt[stI]);
            MaxI = COLLAR(MaxI,
                          ladder_data->CapSt[stI], 
                          ladder_data->FloorSt[stI]);

            MinState = MIN( MIN( MinI * StepN, MinI * StepX), 
                            MIN( MaxI * StepN, MaxI * StepX) );
            MaxState = MAX( MAX( MinI * StepN, MinI * StepX), 
                            MAX( MaxI * StepN, MaxI * StepX) );
        }


        /* RIB case: add {0,1} bounds for observation variable 
         * Assume MinState <= Max State                       */
        if (ladder_data->CplxIsRib)
        {
            double RibWeight = 0;
            int    i;

            for (i = 0; i < ladder_data->NbRibObsDates; i++)
            {
                RibWeight = MAX (RibWeight,
                            MAX (fabs (ladder_data->RibInRangeWeight[i]),
                                 fabs (ladder_data->RibOutRangeWeight[i])));
            }

            if (MinState - MaxState > TINY)
            {
                DR_Error ("MinState > MaxState\n");
                goto FREE_MEM_AND_RETURN;
            }
            if (MinState > 0) 
            {
                MinState = 0.;
                MaxState *= RibWeight;
            }
            if (MaxState < 0)
            {
                MinState *= RibWeight;
                MaxState = 0.;
            }
        }


        /* Generate the equal spaced state variables between cone
         * boundary.
         */
        /* Special case on Value Date - needs review */
        if (StEventList->Dates[i] == ValueDate)
        {
            for (p=0; p<ladder_data->NbStates; p++)
                ladder_data->State[FirstResetI][p] = MinState;
        } 
        else
        {
            deltaS = (MaxState - MinState) / (NbEqSt - 1);
            for (s=0; s<NbEqSt; s++)
            {
                ladder_data->State[i][s] = MinState + deltaS * s;
            }
        }


        
        /* ----- 2. MC simulation -----*/
        for (p=0; p<NbPaths; p++)
        {
            /* compute YieldSlice */
            Yield[p] = YieldSlice[0][p] * obswt[0] + YieldSlice[1][p] * obswt[1];
        
            /* Binary */
            if (Yield[p] <= lb)
            {
                StepUp = dr; 
            }
            else if (Yield[p] <= hb)
            {
                StepUp = mr; 
            }
            else
            {
                StepUp = ur; 
            }

            if (ladder_data->AoM == 'A')
            {
                /* Linear rate and spread */
                StepUp += lev * COLLAR(wt[0] * YieldSlice[0][p] + wt[1] * YieldSlice[1][p] +
                                       spd,
                                       idxcap,
                                       idxfloor);

                /* Life cap/floor */
                ladder_data->State[i][p+NbEqSt]
                    = COLLAR (ladder_data->State[i-1][p+NbEqSt] * stC + StepUp, 
                                ladder_data->CapSt[stI], 
                                ladder_data->FloorSt[stI]); 
            }
            else if (ladder_data->AoM == 'B')
            {
                /* Linear rate and spread */
                StepUp *= lev * COLLAR(wt[0] * YieldSlice[0][p] + wt[1] * YieldSlice[1][p] +
                                       spd,
                                       idxcap,
                                       idxfloor);

                /* Life cap/floor */
                ladder_data->State[i][p+NbEqSt]
                    = COLLAR (ladder_data->State[i-1][p+NbEqSt] * stC + StepUp, 
                                ladder_data->CapSt[stI], 
                                ladder_data->FloorSt[stI]);
            }
            else
            {
                ladder_data->State[i][p+NbEqSt]
                      = COLLAR (ladder_data->State[i-1][p+NbEqSt] * stC + 
                                lev * COLLAR(wt[0] * YieldSlice[0][p] + wt[1] * YieldSlice[1][p] + spd,
                                             idxcap,
                                             idxfloor),                     
                                ladder_data->CapSt[stI],
                                ladder_data->FloorSt[stI]) * StepUp;
            }
        }


    } /* for each reset date >= value date */


    
    /* Sort the state variable in acsending order. 
     * Performed only after the MC simulation 
     */
    for (i = FirstResetI; i < EndResetI; i++)
    {
        if (Fix3_DoubleVectSort(ladder_data->State[i], ladder_data->NbStates) 
                == FAILURE)
            goto FREE_MEM_AND_RETURN;
    }


    /******************************************/
    /*  EVT_STATS : Event stats date          */
    /******************************************/

    for (i = 0; i < ladder_data->OptNbStats; i++)
    {
        int idx = GetDLOffset(ExerEventList->NbEntries,
                ExerEventList->Dates, ladder_data->OptStatDates[i], CbkLOWER);

        if (idx < 0)
            continue;

        if (ExerEventList->Dates[idx] > ValueDate)
        {
            if (Add_To_DateList(&NbCritDate, &CritDate,
                                ExerEventList->Dates[idx],
                                11,
                                0,0,0,0,0,0,0,0) == FAILURE) 
                goto FREE_MEM_AND_RETURN;
        }
    }

    /****************************************/
    /*  TYPE 9: RIB RESET DATES             */ 
    /****************************************/

    for (i = 0; i < StEventList->NbEntries - 1; i++)
    {
        ladder_data->NbRibObsInPer[i] = 0;
    }

    if (ladder_data->CplxIsRib)
    {
        int CplxPer = 0, RibIdxInPer = 0;    /* index of current Rib observation */

        for (i = 0; i < ladder_data->NbRibObsDates; i++)
        {
            CurrResetDate = ladder_data->RibObsEffDate[i];

            if (CurrResetDate >= ladder_data->MatDate) break;

            /* identify corresponding complex coupon */
            while ((CplxPer < StEventList->NbEntries - 1) &&
                   (CurrResetDate >= StEventList->Dates[CplxPer+1]))
            {
                CplxPer++;
            }

            (ladder_data->NbRibObsInPer[CplxPer])++;

            if (CurrResetDate >= ValueDate)
            {
                if (Add_To_DateList (&NbCritDate,
                                     &CritDate,
                                     CurrResetDate,
                                     9,
                                     i,
                                     RibIdxInPer,
                                     0, 0, 0,
                                     0, 0, 0) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Prepare reset datelist for idx zerobank if needed */
                if (AddDateToList (&NbResetRib,
                                   &ResetDLRib,
                                   CurrResetDate) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }

            else if (StEventList->Dates[CplxPer+1] >= ValueDate) 
            {
                (ladder_data->NbPastRibObsDates)++;
            }
            /* Also check at this point that Rib in / out weights
             * do not step up in the middle of a complex period */
            if (RibIdxInPer > 0)   /* Not the first observation */
            {
                if ((ABS (ladder_data->RibInRangeWeight[i] - 
                          ladder_data->RibInRangeWeight[i-1]) > TINY) ||
                    (ABS (ladder_data->RibOutRangeWeight[i] -
                          ladder_data->RibOutRangeWeight[i-1]) > TINY))
                {
                    DR_Error ("Rib in/out range weights cannot step up in the "
                              "middle of a complex period");
                    goto FREE_MEM_AND_RETURN;
                }
            }

            /* Advance or reset Rib index in period */
            if ((i < ladder_data->NbRibObsDates - 1) &&
                (ladder_data->RibObsEffDate[i+1] < StEventList->Dates[CplxPer+1]))
            {
                RibIdxInPer++;
            }
            else
            {
                RibIdxInPer = 0;
            }
        }

        /* Check at least one Rib observation in each period and
         * add critical dates */
        ladder_data->MaxNbRib = 0;
        for (i = 0; i < StEventList->NbEntries - 1; i++)
        {
            int RibIdx = 0;   /* first Rib obs in this period */

            CurrResetDate = StEventList->Dates[i+1];
            if (CurrResetDate < ValueDate) continue;
          
            if (ladder_data->NbRibObsInPer[i] == 0)
            {
                DR_Error ("Must have at least one Rib observation date "
                          "in ladder period #%d\n", i);
                goto FREE_MEM_AND_RETURN;
            }

            while ((RibIdx < ladder_data->NbRibObsDates - 1) &&
                   (ladder_data->RibObsEffDate[RibIdx+1] < CurrResetDate))
            {
                RibIdx++;
            }

            ladder_data->MaxNbRib = 
                MAX (ladder_data->MaxNbRib, ladder_data->NbRibObsInPer[i]);


            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 CurrResetDate,
                                 10,
                                 ladder_data->NbRibObsInPer[i],
                                 ladder_data->RibInRangeWeight [RibIdx],
                                 ladder_data->RibOutRangeWeight[RibIdx],
                                 0, 0,
                                 0, 0, 0) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
        }
    } /* If RIB */

    /********************************************************************/
    /*  TYPE ZbkEVENT to ZbkEVENT + 2: CURVE 0,1,2 ZERO MATURITY DATES  */
    /********************************************************************/

    /* First generate zero mats for ladder index and flt index*/
    for(j = 0; j < 2; j++)
    {
        if (ZbkDLFromIdx(NbResetSt, 
                     ResetDLSt,
                     ResetDLSt,
                     ladder_data->IdxMatSt[j],
                     ladder_data->IdxFreqSt[j],
                     &(NbIdxZMatDates[IdxCurveSt[j]]),
                     &(IdxZMatDL[IdxCurveSt[j]]),
                     &(NbIdxZUseDates[IdxCurveSt[j]]),
                     &(IdxZUseDL[IdxCurveSt[j]])) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    if (ZbkDLFromIdx(NbResetFl, 
                     ResetDLFl,
                     ResetDLFl,
                     ladder_data->IdxMatFl,
                     ladder_data->IdxFreqFl,
                     &(NbIdxZMatDates[IdxCurveFl]),
                     &(IdxZMatDL[IdxCurveFl]),
                     &(NbIdxZUseDates[IdxCurveFl]),
                     &(IdxZUseDL[IdxCurveFl])) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }
    
    if (ladder_data->CplxIsRib)
    {
        /* ... Rib index, */
        for (j = 0; j < 2; j++)
        {
            if (ladder_data->RibIdxOn[j])
            {
                if (ZbkDLFromIdx (NbResetRib,
                                  ResetDLRib,
                                  ResetDLRib,
                                  ladder_data->RibIdxMat[j],
                                  ladder_data->RibIdxFreq[j],
                                  &(NbIdxZMatDates[ladder_data->RibIdxIoD[j]]),
                                  &(IdxZMatDL[ladder_data->RibIdxIoD[j]]),
                                  &(NbIdxZUseDates[ladder_data->RibIdxIoD[j]]),
                                  &(IdxZUseDL[ladder_data->RibIdxIoD[j]])) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }
            }
        }
    }

    /* optimise the zerobank dates */
    for(j = 0; j < 2; j++)
    {
        if (ZbkOptDates(NbIdxZMatDates[j],
                    IdxZMatDL[j],
                    NbIdxZUseDates[j],
                    IdxZUseDL[j],
                    NbCritDate,
                    CritDate,
                    ValueDate,
                    &(NbZeroMatDates[j]),
                    &(ZeroMatDL[j]),
                    &(NbZeroUseDates[j]),
                    &(ZeroUseDL[j])) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
    }

    /* Do curve 0,1,2 in turn */

    for (j=0, event=ZbkEVENT; j<3; j++, event++)
    {
        tree_data->CritType[event] = 'D';

        /* sort/merge/check the date list */
        if (CbkProcessDL(&(NbZeroMatDates[j]),
                         &(ZeroMatDL[j]),
                         &(NbZeroUseDates[j]),
                         &(ZeroUseDL[j])) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

        /* Finally we can add these to the critical date list */
        /* and calc the zerobank size                         */
        for (i=0; i<NbZeroMatDates[j]; i++)
        {
            if (Add_To_DateList(&NbCritDate,
                                &CritDate,
                                ZeroMatDL[j][i],
                                event,
                                0,0,0,0,0,
                                ZeroUseDL[j][i],
                                0,0) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        } 

        tree_data->NbZeros[event] = Fix3_CbkCalcSize
                                    (NbCritDate,
                                     CritDate, /* critical datelist         */
                                     event,    /* eval event type           */
                                     0);       /* earliest use suppdate idx */
    } /* for j = 0,1,2 */

    /* Finally construct the time line using 'I'ncreasing time steps. */
    if (Fix3_Time_Line (ValueDate,                      
                   NbCritDate,
                   CritDate,
                   'I',                            
                   tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    
    status = SUCCESS;

    FREE_MEM_AND_RETURN:                            

    Free_DR_Array (CritDate,  CRITDATE,   0, NbCritDate-1);
    Free_DR_Array (ResetDLRib,    LONG,   0, NbResetRib-1);
    Free_DR_Array (ResetDLSt,     LONG,   0, NbResetSt-1);
    Free_DR_Array (ResetDLFl,     LONG,   0, NbResetFl-1);

    for (j=0;j<3;j++)
    {
        Free_DR_Array(ZeroMatDL[j],LONG,0,NbZeroMatDates[j]-1);
        Free_DR_Array(ZeroUseDL[j],LONG,0,NbZeroUseDates[j]-1);

        Free_DR_Array(IdxZMatDL[j],LONG,0,NbIdxZMatDates[j]-1);
        Free_DR_Array(IdxZUseDL[j],LONG,0,NbIdxZUseDates[j]-1);
    }
    
    DrFreeEventList(ExerEventList);
    DrFreeEventList(StEventList);
    DrFreeEventList(FlEventList);
        
    if (status ==FAILURE)
    {
        DR_Error("Ladder_Schedule: Failed.\n");
    }
    
    return (status);

}  /* Ladder_Schedule */

