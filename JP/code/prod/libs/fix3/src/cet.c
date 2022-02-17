/****************************************************************************/
/*      Calibration Enhancement Tool                                        */
/****************************************************************************/
/*      CET.C                                                               */
/****************************************************************************/



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix123head.h"

#ifndef  MIN_CET_VOL
#define  MIN_CET_VOL  0.0001
#endif

#ifndef  CET_CONVERGENCE_CRITERIA
#define  CET_CONVERGENCE_CRITERIA  1E-4  /* in vegas */
#endif

/* CET error that gives failure. Change from 2.0 t1o 0.5 requested by DA Jan03
*  JBL 31/1/03.*/
#ifndef  CET_ERROR_TOLERANCE
#define  CET_ERROR_TOLERANCE  0.5        /* in vegas */
#endif

#ifndef  CET_WARNING_TOLERANCE
#define  CET_WARNING_TOLERANCE  0.05      /* in vegas */
#endif

#ifndef  MASKLENGTH
#define  MASKLENGTH   24L   /* in months */
#endif




/*****  Fix3_Cet_Classic  ********************************************************/
/*
*       Main subroutine.
*/
int     Fix3_Cet_Classic 
            (int                    CetOutputFlag,  /* (I) On TRUE write to file */          
             T_CURVE const*         t_curve,        /* (I) Zero curve data       */
             MKTVOL_DATA*           mktvol_dataProd,/* (I/0) Vol  data           */
             FIX3_TREE_DATA const*  tree_dataProd)  /* (I) Tree data             */
{


    MKTVOL_DATA          mktvol_data;   /* Structure of vol data             */
    FIX3_TREE_DATA       tree_data;     /* Structure of tree data            */
    CET_OUT_DATA         cet_out_data;  /* Output: the "enhanced" mkt vols   */




    long           ValueDate;          /* From zero curve                   */
    int            NbInstr;            /* Nb of benchmark swaptions         */
    int            IterIdx;            /* Iteration counter                 */
    double         TimeToExp;          /* Time to exp of each b'mark option */
    double         slopeFO;            /* 1st order 1st derv in NRaphson    */
    int            i;

    int            isFirstWarning = TRUE;  /* Prints warning header info    */

    int            status = FAILURE;   /* Status */ 
    char           ErrorMsg[MAXBUFF];
                 

    /* Return original vols if not enhanced */
    if ((mktvol_dataProd->CetNbIter == 0) || 
        (mktvol_dataProd->CalibFlag == FALSE))
    {
        status = SUCCESS;
        return (status);
    }
    if (mktvol_dataProd->CetNbIter < 0)
    {
        DR_Error("Number of cetIterations is invalid (%d)",
                 mktvol_dataProd->CetNbIter);
        DR_Error("%s failed.", "Fix3_Cet_Classic");
        status = FAILURE;
        return (status);
    }

    /* Each vol point becomes a critical event type on tree_data */
    /* and therefore total nb vols cannot exceed limit           */
    if (mktvol_dataProd->NbVol > NBCRITDATE)
    {
        sprintf (ErrorMsg, "Nb of vols (%d) exceeds maximum limit "
                 "of %d!", mktvol_dataProd->NbVol, NBCRITDATE);
        DR_Error (ErrorMsg);
        goto FREE_MEM_AND_RETURN;

    }

    /*  Initialize tree data structure to NULL. */
    Fix3_Tree_Init (&tree_data);

    /*  I/O manager: Initialize according to product data. */
    if (Fix3_Cet_Manager(mktvol_dataProd,
                    tree_dataProd,
                    &mktvol_data, 
                    &tree_data) == FAILURE)
    {                                      
        goto FREE_MEM_AND_RETURN;
    }
    ValueDate = (t_curve[tree_data.CvDiff]).ValueDate;
    NbInstr = mktvol_data.NbVol;


    /* Black & Scholes prices */
    for (i=0; i<NbInstr; i++)
    {
        if (ParYieldFromDates(&(cet_out_data.ParYield[i]),
                              &(cet_out_data.Annuity[i]), 
                              mktvol_data.SwapSt[i], 
                              mktvol_data.SwapMat[i],       
                              mktvol_data.DCC,           
                              mktvol_data.Freq,          
                              'F',  /* If ever a stub, it will be front */    
                              &(t_curve[tree_data.CvDiff])) != SUCCESS)
        {
           goto FREE_MEM_AND_RETURN;
        }


        /* Price with B&S */
        TimeToExp = Daysact(ValueDate,mktvol_data.VolDate[i])/365.;
        if (TimeToExp < TINY)
        {
            sprintf (ErrorMsg, "B'mark swap start %ld before value date %ld!",
                     YMDDateFromIRDate(mktvol_data.SwapSt[i]), ValueDate);
            DR_Error (ErrorMsg);
            goto FREE_MEM_AND_RETURN;
        }
        
        if (mktvol_data.VolUnit == 1)
        {
            cet_out_data.BSPrice[i] = Put_BS(cet_out_data.ParYield[i],
                                             cet_out_data.ParYield[i],
                                             TimeToExp,
                                             0.0,  /* risk free rate set to 0 */
                                             mktvol_data.Vol[i]);

            cet_out_data.BSVega[i] = Vega_BS(cet_out_data.ParYield[i],
                                             cet_out_data.ParYield[i],
                                             TimeToExp,
                                             0.0,
                                             mktvol_data.Vol[i]);
        }
        else
        {
            cet_out_data.BSPrice[i] = Option_Normal(cet_out_data.ParYield[i],
                                                    cet_out_data.ParYield[i],
                                                    'P',
                                                    TimeToExp,
                                                    0.0,  /* risk free rate set to 0 */
                                                    mktvol_data.Vol[i]); 

            cet_out_data.BSVega[i] = Vega_Normal(cet_out_data.ParYield[i],
                                                 cet_out_data.ParYield[i],
                                                 TimeToExp,
                                                 0.0,
                                                 mktvol_data.Vol[i]); 
        }
        
        
    
        /* Annuity includes the discounting */
        cet_out_data.BSPrice[i] *= cet_out_data.Annuity[i];
        cet_out_data.BSVega[i]  *= cet_out_data.Annuity[i];

    }



    for (IterIdx = 1; IterIdx <= mktvol_dataProd->CetNbIter; IterIdx++)
    {
        
        /* 1 - Build time line with dates of b'mark instruments.      */
        /*  Out_Data is passed so that paryields are stored alongside */
        /*  the timeline; tree_dataProd is passed in order to match   */
        /*  the product & CET timelines (for short-dated swaptions)   */
        if (Fix3_Cet_Schedule(ValueDate,                                 
                         &cet_out_data,                             
                         &mktvol_data,                              
                         tree_dataProd,
                         &tree_data) == FAILURE)                    
        {                                                           
            goto FREE_MEM_AND_RETURN;                               
        }                              
                                                                    

        /*  2 - Build tree from current market vols */
        if (Fix3_Build_Tree (t_curve,
                        &mktvol_data,
                        &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

        
        /*  3 - Price benchmarks on tree */
        if (Fix3_Cet_Calc(&mktvol_data,                                  
                     &tree_data,
                     &cet_out_data) == FAILURE)
        {
            DR_Error ("Could not run enhancement algorithm!");
            goto FREE_MEM_AND_RETURN;
        }


        /*  4 -  Update target diffs and improve vols on used vols */
        for (i=0; i<NbInstr; i++)
        {
            if (mktvol_data.VolUsed[i]) 
            {
                /* CET tree is effectively built with filtered vol */
                /* so use filtered vols in NR routine              */
                if (Fix3_IndexVol(&(cet_out_data.FilterVol[i]),
                             mktvol_data.VolDate[i],
                             mktvol_data.SwapSt[i],
                             mktvol_data.SwapMat[i],
                             mktvol_data.Freq,
                             mktvol_data.DCC,
                             &mktvol_data,
                             &(t_curve[tree_data.CvDiff])) != SUCCESS)
                {
                    DR_Error ("Could not calculate filtered CET vol");
                    goto FREE_MEM_AND_RETURN;
                }

                cet_out_data.VolPrev[i]             = mktvol_data.Vol[i]; /* just for display in CET.prn */  
                cet_out_data.PriceDiff[i]     = cet_out_data.BSPrice[i] 
                                              - cet_out_data.TreePrice[i];
                cet_out_data.TreePricePrev[i] = cet_out_data.TreePrice[i];

                /* aprroximate NR derivative */
                slopeFO  = cet_out_data.TreePrice[i];
                slopeFO /= cet_out_data.FilterVol[i];
                cet_out_data.FOVega[i] = slopeFO;

                /* perform NR iteration */
                mktvol_data.Vol[i] = cet_out_data.FilterVol[i] + 
                                   cet_out_data.PriceDiff[i]/slopeFO;

                /* ensure CET vol is positive */
                if (IterIdx < mktvol_dataProd->CetNbIter)
                {
                    mktvol_data.Vol[i] = MAX(mktvol_data.Vol[i], MIN_CET_VOL);

                }
                else
                {
                    if (mktvol_data.Vol[i] < MIN_CET_VOL)
                    {
                        sprintf (ErrorMsg, "Final CET vol %ld is too low",
                                 mktvol_data.SwapSt[i]);
                        DR_Error (ErrorMsg);
                        goto FREE_MEM_AND_RETURN;
                    }
                }

                cet_out_data.PriceDiffInVega[i] = 100. * cet_out_data.PriceDiff[i]/
                                                  cet_out_data.BSVega[i];
            }
        }

        /*  5 - Print to VOL.prn file */
        if (CetOutputFlag)
        {
            if (Fix3_Cet_Print_Vol_File(&cet_out_data,
                                   &mktvol_data,
                                   NbInstr,
                                   IterIdx) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            } 
        }
        {
            int NbFactor = tree_data.NbFactor;
            int NbTP = tree_data.NbTP;
            int NbEDevDates = tree_data.NbEDevDates;
            int JumpPpy = tree_data.JumpPpy;
            int NbDailyPts = tree_data.NbDailyPts;
            int NbNmr = tree_data.NbNmr;
            Fix3_Tree_Free(&tree_data);
            tree_data.NbFactor = NbFactor;
            tree_data.NbTP = NbTP;
            tree_data.NbEDevDates = NbEDevDates;
            tree_data.JumpPpy = JumpPpy;
            tree_data.NbDailyPts = NbDailyPts;
            tree_data.NbNmr = NbNmr;
        }
        /* re-initialize pointers back to null */
        Fix3_Cet_Tree_Reset(&tree_data);
    }

    /* Check final CET vols are tolerable; report error/warning as appropriate; record new vols */
    /* Do this AFTER call to Fix3_Cet_Print_Vol_File for accesss to full debug info via CET.prn      */ 
    for (i=0; i<NbInstr; i++)
    {        
        if (mktvol_data.VolUsed[i])
        {
            if (fabs(cet_out_data.PriceDiffInVega[i]) > CET_ERROR_TOLERANCE)
            {
                /* print error to stdout and exit */
                DR_Error("Final CET vol has not converged to within %4.2lf vega!\n"
                         "Expiry = %ld SwapSt = %ld SwapMat = %ld VegaDiff = %4.2lf vega\n",
                         CET_ERROR_TOLERANCE, 
                         YMDDateFromIRDate(mktvol_data.VolDate[i]),
                         YMDDateFromIRDate(mktvol_data.SwapSt[i]), 
                         YMDDateFromIRDate(mktvol_data.SwapMat[i]),
                         cet_out_data.PriceDiffInVega[i]);

                goto FREE_MEM_AND_RETURN;
            }
            else
            {
                if (fabs(cet_out_data.PriceDiffInVega[i]) > CET_WARNING_TOLERANCE)
                {
                    if (isFirstWarning)
                    {
                        /* print header info */
                        DR_Warning(TRUE, /* open clean file */
                                   "Fix3_Cet_Main Warning: Final CET vol has not converged to within %6.4lf vega!\n"
                                   "Expiry  \t SwapSt  \t SwapMat \t VegaDiff", 
                                   CET_WARNING_TOLERANCE);

                        isFirstWarning = FALSE;
                    }
                    
                    /* print warning to stdout and warning.log but do not exit */
                    DR_Warning(FALSE, /* append to file */
                               "%ld\t %ld\t %ld\t %6.4lf", 
                               mktvol_data.VolDate[i],
                               mktvol_data.SwapSt[i], 
                               mktvol_data.SwapMat[i],
                               cet_out_data.PriceDiffInVega[i]);
                }
            }
        }

        /* record new vols since we do not exit */
        mktvol_dataProd->Vol[i]     = mktvol_data.Vol[i];
        mktvol_dataProd->VolUsed[i] = mktvol_data.VolUsed[i];
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Fix3_Tree_Free (&tree_data);

    if (status != SUCCESS)
        DR_Error("%s failed", "Fix3_Cet_Classic");
                                               
    return (status);

}  /* Fix3_Cet_Classic */


/*****  Fix3_Cet_Calc ************************************************************/
/*
*       Main calculation routine: discounted expected value of cash-flows
*       going backward in the tree.
*/
int     Fix3_Cet_Calc
             (MKTVOL_DATA         *mktvol_data,    /* (I) Vol data          */
              FIX3_TREE_DATA      *tree_data,      /* (I) Tree data         */
              CET_OUT_DATA        *cet_out_data)   /* (I) Output data       */
{


    FIX3_DEV_DATA    dev_data;          /* Fix3_Dev data structure        */

    /* Variables */                   
    double      *BMarkPrice[MAX_INST];  /* Prices of b'mark swaps         */
    double      *Swaption = NULL;       /* Prices of b'mark swaptions     */
    double      *AuxSlice = NULL;       /* Used for smoothing             */

    /* Flags & amounts */                       
    long        ExerFlag[MAX_INST];     /* Indicates an exercise date     */
    long        FixCpnFlag[MAX_INST];   /* Indicates a fixed payment date */
    double      FixCpnAmt[MAX_INST];
    int         Smoothing;              /* Use smoothing if TRUE */

    /* Numeric variables */
    long        CurrentDate;       /* Current date                        */
    long        ValueDate;         /* Value date                          */
    int         ICurve;
    int         DCurve;            /* Discount curve number               */
    int         NbInstr;           /* Nb of instruments to price          */
   
    double      Aux;
    int         T;                 /* Last time point                     */
    int         t;                 /* Current time point                  */
    int         OffsetAt0;         /* Node offset at t=0                  */
    int         k;                 /* Convenience index                   */
    int         status = FAILURE;  /* Error status = FAILURE initially    */

    
    /* initialise pointers to NULL */

    Fix3_Dev_Init(&dev_data);

    for (k=0; k<MAX_INST; k++) BMarkPrice[k] = NULL;


    /* Total number of time points and value date */
    T   = tree_data->NbTP;        
    ValueDate = tree_data->TPDate[0];
    
    /* Number of market instruments (swaptions) to price */
    NbInstr = mktvol_data->NbVol;

    /* Assigment of discount and index curves in the engine */
    ICurve = tree_data->CvDiff;
    DCurve = tree_data->CvDisc;

    /* Smoothing flag */
    Smoothing = (mktvol_data->SmoothingFlag == 'Y');

    OffsetAt0 = Fix3_Node_Offset(tree_data->NbFactor, 0, 0, 0, tree_data);


    /*   Allocation of variables for tree pricing and initialisation */
    for (k=0; k<NbInstr; k++)
    {
        BMarkPrice[k] = Fix3_Alloc_Slice(tree_data);
        if (BMarkPrice[k] == NULL)
        {
            DR_Error("could not allocate memory for b'mark swaps!");
            goto FREE_MEM_AND_RETURN;        
        }
    }

    Swaption = Fix3_Alloc_Slice(tree_data);
    AuxSlice = Fix3_Alloc_Slice(tree_data);
    if ( (Swaption == NULL) || (AuxSlice == NULL) ) 
    {
        DR_Error("could not allocate memory for b'mark swaptions!");
        goto FREE_MEM_AND_RETURN;        
    }

    if (Fix3_Dev_Alloc (&dev_data, tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    
    /* One single walk back on the tree for all instruments */
    for (t = T; t >= 0; t--)
    {

        CurrentDate = tree_data->TPDate[t];

        /* Set flags for coupon pmts (final principal included) and exercise */
        for (k=0; k<NbInstr; k++)
        {
            FixCpnFlag[k] = tree_data->TPtype[k][t]; 
            if (FixCpnFlag[k])
            {
                FixCpnAmt[k] = (tree_data->CritDate[k][t]).Value[0]; 
            }
        
            ExerFlag[k] = FALSE;
            if (CurrentDate == mktvol_data->VolDate[k])
            {
                ExerFlag[k] = TRUE;
            }
        } 


        /*  'Update' tree. */
        if (Fix3_Lattice (&dev_data,
                     t,                                      
                     T,
                     mktvol_data,                                   
                     tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }


        for (k=0; k<NbInstr; k++)
        {
            /* No need to price the instrument */
            if (!(mktvol_data->VolUsed[k]))
                continue;


            if ((mktvol_data->SwapMat[k] >= CurrentDate) &&
                (mktvol_data->SwapSt[k]  <= CurrentDate))
            {
                /* DEV benchmark price */
                if (Fix3_Dev(BMarkPrice[k],    
                        t,        
                        T,        
                        ICurve,   
                        &dev_data,
                        tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Add coupon (final principal included) */
                if (FixCpnFlag[k])
                {
                    if (Fix3_AddScalar(BMarkPrice[k],
                                  FixCpnAmt[k],
                                  t,
                                  tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }

                /* If Exercise date, evaluate pay-off (par strike) */
                if (ExerFlag[k])
                {
                    /* clean option slice */
                    if (Fix3_Set_Slice(Swaption,
                                  0.,
                                  t,
                                  tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                    /* evaluate option claim */
                    if (Fix3_OptionPlus_t(Swaption,
                                     NULL,          /* Statistics not needed */
                                     NULL,          /* Statistics not needed */
                                     NULL,          /* Statistics not needed */
                                     BMarkPrice[k],
                                     1.,            /* Notional   */
                                     0.,            /* Strike     */
                                     ExerFlag[k],
                                     1,             /* always price call = 1 */
                                     Smoothing,
                                     AuxSlice,
                                     t,
                                     t,             /* no discounting */
                                     DCurve,        /* => not used    */
                                     &dev_data,
                                     tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                    /* Perform express DEV to finalise */
                    if (CurrentDate != tree_data->EDevDate[k])
                    {
                        DR_Error("Program bug: express DEV date incorrectly "
                                 "processed.\nState prices not available.\n");
                        goto FREE_MEM_AND_RETURN;

                    }

                    if (Fix3_MultTwoSlicesAddAll(&Aux,
                                            Swaption,
                                            tree_data->EDevStPrice[k],
                                            t,
                                            tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                    (BMarkPrice[k]+OffsetAt0)[0] = Aux;


                } /* If Exercise date */

            } /* If before final maturity of k-th benchmark */ 

        } /* For k */
        
    }  /* for t */



    /* Record prices. */
    
    for (k=0; k<NbInstr; k++)
    {
        cet_out_data->TreePrice[k] = (BMarkPrice[k] + OffsetAt0)[0];
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:
    
    for (k=0; k<NbInstr; k++)
    {
        Fix3_Free_Slice(BMarkPrice[k],tree_data);
    }    

    Fix3_Free_Slice(Swaption,tree_data);
    Fix3_Free_Slice(AuxSlice,tree_data);


    Fix3_Dev_Free (&dev_data, tree_data);

    if (status != SUCCESS)
        DR_Error("%s failed", "Fix3_Cet_Calc");

    return (status);

}  /* Calc_Cet */




/*****  Fix3_Cet_Schedule ***********************************************/
/*
*       Sets up the time line. 
*
*       For the calibration tool, the dates are all the payment
*       dates of the  instruments (swap or money market) chosen
*       as benchmarks for vol.
*          
*       tree_dataProd is passed in order to match the CET and
*       product timelines (for short-dated swaptions) out to 1y.
*
*       This routine deals with the detailed date generation from input
*       parameters and the subsequent addition of dates to the timeline
*       critical date list.
*
*/

int     Fix3_Cet_Schedule 
            (long                   ValueDate,      /* (I) Value date             */
             CET_OUT_DATA*          cet_out_data,
             MKTVOL_DATA*           mktvol_data,
             FIX3_TREE_DATA const*  tree_dataProd, /* (I) Tree data for product  */
             FIX3_TREE_DATA*        tree_data)     /* (O) Tree data for CET      */

{


    int         NbCritDate = 0;
    CRIT_DATE   *CritDate = NULL;        /* Critical date list             */
    
    EVENT_LIST  *PmtEventListFix[NBCRITDATE]; /* Event list for fixed pmts */
    
    int         NbInstr = 0;             /* Nb of instruments to price     */

    int         i, k;                    /* Index for convenience          */
    int         status = FAILURE;        /* Status = FAILURE initially     */
    
    long        CurrCouponDate;          /* Date of prossesed coupon       */
    double      CurrCouponDayCount;      /* Day count of prossesed coupon  */
    double      CurrCouponAmt;
 
    /* Only CET those vols that are needed to build the product tree */ 
    NbInstr = mktvol_data->NbCetVol;

    /* Check to ensure that there are enough critical dates */
    /*       for all chosen benchmark instruments           */
    if (NbInstr > NBCRITDATE)
    {
        DR_Error("The nb of swaptions to calibrate (%d) exceeds the pre-set\n"
                 "memory allocated for it (%d).\n", NbInstr, NBCRITDATE);
        goto FREE_MEM_AND_RETURN;
    }

    /* Initialise all event lists to NULL */
    for (k=0; k<NbInstr; k++) 
    {
        PmtEventListFix[k] = NULL;
    }

    /* EXPRESS DEV DATE LIST */
    /* The memory for the express DEV  dates gets allocated here, */
    /* but is freed via Fix3_Tree_Free if failure occurs at some point */
    tree_data->NbEDevDates = NbInstr;
    tree_data->EDevDate = (long *) DR_Array(LONG, 0, NbInstr-1);
    if (tree_data->EDevDate == NULL)
    {
        goto FREE_MEM_AND_RETURN;
    }


    /* CRITICAL DATE LIST */
    /*  Allocate empty critical date list and then add dates successively.*/
    CritDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (CritDate == NULL)
    {
        DR_Error("Unable to allocate memory for critical dates array !");
        goto FREE_MEM_AND_RETURN;
    }



    /*   Always add value date to critical date list.  */
    if (Add_To_DateList (&NbCritDate,
                         &CritDate,              
                         ValueDate,
                         NBCRITDATE,             /* No specific type */
                         0, 0, 0, 0, 0,              
                         0, 0, 0) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Match CET and product timeline out to 1 year from value date */
    {
        int    t, MaskIdx = 0;
        long   MaskEndDate;

        MaskEndDate = Nxtmth(ValueDate, MASKLENGTH, 1L);

        while ( (MaskIdx < tree_dataProd->NbTP) &&
                (tree_dataProd->TPDate[MaskIdx] <= MaskEndDate) )
        {
            MaskIdx++;
        }

        for (t=0; t<MaskIdx; t++)
        {
            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,              
                                 tree_dataProd->TPDate[t],
                                 NBCRITDATE,             /* No specific type */
                                 0, 0, 0, 0, 0,              
                                 0, 0, 0) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
    }


    /********************************/
    /*                              */
    /*      FIXED COUPON DATES      */
    /*    (Final princ included     */
    /********************************/

    for (k=0; k<NbInstr; k++)
    {

        tree_data->CritType[k] = 'D';      /* Always of discrete type */      
        tree_data->NbZeros[k]  =  0;


        /* Process the swap inputs to generate the appropriate fixed cp  */
        /* payment dates which will be placed in a temporary EVENT_LIST. */
        {
            long TempDates[2];  /* Only accrual start and final mat */
        
            TempDates[0] = mktvol_data->SwapSt[k];
            TempDates[1] = mktvol_data->SwapMat[k];
        
            PmtEventListFix[k] = DrNewEventListFromFreq 
                                  (2,
                                   TempDates,
                                   mktvol_data->Freq,
                                   'F',    /* Forward stub                  */
                                   'N',    /* 'Dates in' check not required */
                                   NULL, NULL, NULL, NULL, NULL);
            
            if (PmtEventListFix[k] == NULL)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
                                                
        
        for (i = 1; i < PmtEventListFix[k]->NbEntries; i++)            
        {
           
            CurrCouponDate = PmtEventListFix[k]->Dates[i];
        
            if (DrDayCountFraction (PmtEventListFix[k]->Dates[i-1],
                                    PmtEventListFix[k]->Dates[i], 
                                    mktvol_data->DCC,
                                    &CurrCouponDayCount) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        
            CurrCouponAmt = CurrCouponDayCount 
                          * cet_out_data->ParYield[k];
            if (i == PmtEventListFix[k]->NbEntries-1)
            {
                CurrCouponAmt += 1.0; /* Final principal */
            }
        
            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 PmtEventListFix[k]->Dates[i],
                                 k,
                                 CurrCouponAmt,
                                 0, 0, 0, 0,
                                 0, 0, 0) == FAILURE)     
            {
                goto FREE_MEM_AND_RETURN;
            }
        
        }  /* for i */

        /* Initial notional */
        if (Add_To_DateList (&NbCritDate,                  
                             &CritDate,                    
                             PmtEventListFix[k]->Dates[0], 
                             k,                            
                             -1.0,                
                             0, 0, 0, 0,                   
                             0, 0, 0) == FAILURE)          
        {                                                  
            goto FREE_MEM_AND_RETURN;                      
        }  


        /* Express DEV date */
        tree_data->EDevDate[k] = PmtEventListFix[k]->Dates[0];

    } /* For k */




    /************************************/
    /*                                  */
    /*  TYPES NbVol TO NBCRITDATE:  N/A */
    /*                                  */ 
    /************************************/       
                    
    for (i = NbInstr; i < NBCRITDATE; i++)
    {
        tree_data->CritType[i] = 'D';
        tree_data->NbZeros[i]  =  0;
    }

    /* 
     *   Finally construct the time line using 'J' time steps
     *   i.e. increasing time steps with modified jump size
     */
    if (Fix3_Time_Line (ValueDate,                      
                   NbCritDate,
                   CritDate,
                   'J',                            
                   tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }



    status = SUCCESS;

    FREE_MEM_AND_RETURN:                            

    Free_DR_Array (CritDate, CRITDATE, 0, NbCritDate-1);

    for (k=0; k<NbInstr; k++)
    {
        DrFreeEventList(PmtEventListFix[k]);
    }

    if (status != SUCCESS)
        DR_Error("%s Failed.\n", "Fix3_Cet_Schedule");

    return (status);

}  /* Fix3_Cet_Schedule */


/*****  Fix3_Cet_Manager  ********************************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Fix3_Cet_Manager 
            (MKTVOL_DATA const*     mktvol_dataProd,  /* (I) Vol  data             */
             FIX3_TREE_DATA const*  tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA*           mktvol_data,      /* (O) Vol  data             */
             FIX3_TREE_DATA*        tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   

    /*  Standard environment elements: curve and vols */
    tree_data->CvDiff      = tree_dataProd->CvDiff;
    tree_data->CvIdx1      = tree_dataProd->CvIdx1;
    tree_data->CvIdx2      = tree_dataProd->CvIdx2;
    tree_data->CvDisc      = tree_dataProd->CvDisc;
    mktvol_data->CalibFlag = mktvol_dataProd->CalibFlag;  
    mktvol_data->SkipFlag  = mktvol_dataProd->SkipFlag;
    mktvol_data->BaseDate  = mktvol_dataProd->BaseDate;
    mktvol_data->Freq      = mktvol_dataProd->Freq;
    mktvol_data->DCC       = mktvol_dataProd->DCC;
    mktvol_data->NbVol     = mktvol_dataProd->NbVol;
    mktvol_data->NbCetVol  = mktvol_dataProd->NbVol;

    for (i = 0; i < mktvol_data->NbVol; i++)
    {                                                                       
        mktvol_data->VolDate[i] = mktvol_dataProd->VolDate[i];
        mktvol_data->SwapSt[i]  = mktvol_dataProd->SwapSt[i];
        mktvol_data->SwapMat[i] = mktvol_dataProd->SwapMat[i];
        mktvol_data->Vol[i]     = mktvol_dataProd->Vol[i];
        mktvol_data->VolUsed[i] = mktvol_dataProd->VolUsed[i];
    }

    mktvol_data->FilterSpotVolFlag = TRUE;
    mktvol_data->SmoothingFlag = mktvol_dataProd->SmoothingFlag;
    mktvol_data->TraceFlag     = mktvol_dataProd->TraceFlag;


    /*  Model specific data: numerics  */
    tree_data->NbSigmaMax = tree_dataProd->NbSigmaMax;
    mktvol_data->NbFactor = mktvol_dataProd->NbFactor;
    tree_data->NbFactor   = tree_dataProd->NbFactor;
    tree_data->JumpPpy    = tree_dataProd->Ppy;
    tree_data->NbDailyPts = tree_dataProd->NbDailyPts;

    /* Efficiency features */

    /* Reduce PPY. Fix3_Cet_Schedule sets CET = product timeline over first 2y */
    if (tree_data->NbFactor == 1)
    {
        tree_data->Ppy = MIN (tree_dataProd->Ppy, tree_dataProd->PpyCet[0]);
    }
    else  if (tree_data->NbFactor == 2)
    {
        tree_data->Ppy = MIN (tree_dataProd->Ppy, tree_dataProd->PpyCet[1]);
    }
    else if (tree_data->NbFactor == 3) 
    {
        tree_data->Ppy = tree_dataProd->Ppy;
    }


    /* Only CET those vols that are used to build the product tree */ 
    if ((tree_dataProd->NbTP != 0) && (tree_dataProd->TPDate != NULL))
    {
        long    T = tree_dataProd->TPDate[tree_dataProd->NbTP];
        int     LastVolUsedIdx = 0;

        if (mktvol_data->VolDate[mktvol_data->NbVol-1] <= T)
        {
            mktvol_data->NbCetVol = mktvol_data->NbVol;
        }
        else
        {
            while ( (LastVolUsedIdx < mktvol_data->NbVol) && 
                    (mktvol_data->VolDate[LastVolUsedIdx] < T) )
            {
                LastVolUsedIdx++;
            }

            mktvol_data->NbCetVol = LastVolUsedIdx+1;

            for (i = mktvol_data->NbCetVol; i < mktvol_data->NbVol; i++)
            {
                mktvol_data->VolUsed[i] = FALSE;
            }
        }
    }

    /* Model parameters */
    mktvol_data->ModelChoice = mktvol_dataProd->ModelChoice; 
    mktvol_data->IsNmrModel  = mktvol_dataProd->IsNmrModel;

    return (Fix3_Cet_Model_Manager (mktvol_dataProd, 
                                    tree_dataProd, 
                                    mktvol_data, 
                                    tree_data));
}


/*****  Fix3_Cet_Model_Manager_Classic  ********************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Fix3_Cet_Model_Manager_Classic
            (MKTVOL_DATA const*     mktvol_dataProd,  /* (I) Vol  data             */
             FIX3_TREE_DATA const*  tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA*           mktvol_data,      /* (O) Vol  data             */
             FIX3_TREE_DATA*        tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   


    /*  Model specific data  */
    mktvol_data->QRight   = mktvol_dataProd->QRight;
    mktvol_data->QLeft    = mktvol_dataProd->QLeft;
    mktvol_data->FwdShift = mktvol_dataProd->FwdShift;
    mktvol_data->Bbq      = mktvol_dataProd->Bbq;
    mktvol_data->VolUnit  = mktvol_dataProd->VolUnit;
    mktvol_data->VolNorm  = mktvol_dataProd->VolNorm;
    mktvol_data->VolLogn  = mktvol_dataProd->VolLogn;

    for (i = 0; i < 3; i++)
    {
        mktvol_data->Alpha[i] = mktvol_dataProd->Alpha[i]; 
        mktvol_data->Beta[i]  = mktvol_dataProd->Beta[i]; 
        mktvol_data->Rho[i]   = mktvol_dataProd->Rho[i];
    }
    
    status = SUCCESS;
        
    return (status);

}  /* Fix3_Cet_Model_Manager_Classic */


/*****  Fix3_Cet_Model_Manager_TimeDep  ********************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Fix3_Cet_Model_Manager_TimeDep
            (MKTVOL_DATA const*     mktvol_dataProd,  /* (I) Vol  data             */
             FIX3_TREE_DATA const*  tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA*           mktvol_data,      /* (O) Vol  data             */
             FIX3_TREE_DATA*        tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i,k;                              /* Convenience index            */   


    /*  Model specific data  */
    mktvol_data->QRight   = mktvol_dataProd->QRight;
    mktvol_data->QLeft    = mktvol_dataProd->QLeft;
    mktvol_data->FwdShift = mktvol_dataProd->FwdShift;
    mktvol_data->Bbq      = mktvol_dataProd->Bbq;
    mktvol_data->VolUnit        = mktvol_dataProd->VolUnit;
    mktvol_data->VolNorm  = mktvol_dataProd->VolNorm;
    mktvol_data->VolLogn  = mktvol_dataProd->VolLogn;

    mktvol_data->NbFactor     = mktvol_dataProd->NbFactor;
    mktvol_data->NbTDInp      = mktvol_dataProd->NbTDInp;
    mktvol_data->NbSmileDates = mktvol_dataProd->NbSmileDates;

    for (i = 0; i < 3; i++)
    {
        for ( k = 0; k < mktvol_data->NbTDInp; k++)
        {
            mktvol_data->AlphaTD[i][k] = mktvol_dataProd->AlphaTD[i][k]; 
            mktvol_data->BetaTD[i][k]  = mktvol_dataProd->BetaTD[i][k]; 
            mktvol_data->RhoTD[i][k]   = mktvol_dataProd->RhoTD[i][k];
        }
    }


    for ( k = 0; k < mktvol_data->NbTDInp; k++)
    {
        mktvol_data->TDInpDate[k] = mktvol_dataProd->TDInpDate[k];
    }


    for (k = 0; k < mktvol_data->NbSmileDates; k++)
    {
        mktvol_data->SmileDate[k]  = mktvol_dataProd->SmileDate[k];
        mktvol_data->QLeftTD[k]    = mktvol_dataProd->QLeftTD[k];
        mktvol_data->QRightTD[k]   = mktvol_dataProd->QRightTD[k];
        mktvol_data->FwdShiftTD[k] = mktvol_dataProd->FwdShiftTD[k];
    }
    
    status = SUCCESS;
        
    return (status);

}  /* Fix3_Cet_Model_Manager_TimeDep */


/*****  Fix3_Cet_Model_Manager_Smd  ************************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Fix3_Cet_Model_Manager_Smd
            (MKTVOL_DATA const*     mktvol_dataProd,  /* (I) Vol  data             */
             FIX3_TREE_DATA const*  tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA*           mktvol_data,      /* (O) Vol  data             */
             FIX3_TREE_DATA*        tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   


    /*  Model specific data  */
    mktvol_data->QRight   = mktvol_dataProd->QRight;
    mktvol_data->QLeft    = mktvol_dataProd->QLeft;
    mktvol_data->FwdShift = mktvol_dataProd->FwdShift;
    mktvol_data->Bbq      = mktvol_dataProd->Bbq;
    mktvol_data->VolUnit        = mktvol_dataProd->VolUnit;
    mktvol_data->VolNorm  = mktvol_dataProd->VolNorm;
    mktvol_data->VolLogn  = mktvol_dataProd->VolLogn;

    for (i = 0; i < 3; i++)
    {
        mktvol_data->Alpha[i] = mktvol_dataProd->Alpha[i]; 
        mktvol_data->Beta[i]  = mktvol_dataProd->Beta[i]; 
        mktvol_data->Rho[i]   = mktvol_dataProd->Rho[i];
    }

    mktvol_data->Afac   = mktvol_dataProd->Afac;
    mktvol_data->Bfac   = mktvol_dataProd->Bfac;
    mktvol_data->Cfac   = mktvol_dataProd->Cfac;
    mktvol_data->Dfac   = mktvol_dataProd->Dfac;
    
    status = SUCCESS;
        
    return (status);

}  /* Fix3_Cet_Model_Manager_Smd */



/*****  Fix3_Cet_Model_Manager_E2Q  ********************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Fix3_Cet_Model_Manager_E2Q
            (MKTVOL_DATA const*     mktvol_dataProd,  /* (I) Vol  data             */
             FIX3_TREE_DATA const*  tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA*           mktvol_data,      /* (O) Vol  data             */
             FIX3_TREE_DATA*        tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   


    /*  Model specific data  */
    mktvol_data->QRight         = mktvol_dataProd->QRight;
    mktvol_data->QLeft          = mktvol_dataProd->QLeft;
    mktvol_data->Amap           = mktvol_dataProd->Amap;
    mktvol_data->Bmap           = mktvol_dataProd->Bmap;
    mktvol_data->VolUnit        = mktvol_dataProd->VolUnit;
    mktvol_data->FwdShift       = mktvol_dataProd->FwdShift;
    mktvol_data->Bbq            = mktvol_dataProd->Bbq;
    mktvol_data->VolNorm        = mktvol_dataProd->VolNorm;
    mktvol_data->VolLogn        = mktvol_dataProd->VolLogn;

    for (i = 0; i < 3; i++)
    {
        mktvol_data->Alpha[i] = mktvol_dataProd->Alpha[i]; 
        mktvol_data->Beta[i]  = mktvol_dataProd->Beta[i]; 
        mktvol_data->Rho[i]   = mktvol_dataProd->Rho[i];
    }
    
    status = SUCCESS;
        
    return (status);

}  /* Fix3_Cet_Model_Manager_E2Q */

/*****  Fix3_Cet_Tree_Reset  ******************************************************/
/*
*       Re-initialize tree pointers to NULL.
*/
void    Fix3_Cet_Tree_Reset(FIX3_TREE_DATA *tree_data) /* Tree building data structure */
{
    int
            i;


    tree_data->NbTP = 0;

    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->CritDate[i] = NULL;
        tree_data->TPtype[i]   = NULL;
    }

    tree_data->Top1    = NULL;
    tree_data->Bottom1 = NULL;
    tree_data->Top2    = NULL;
    tree_data->Bottom2 = NULL;
    tree_data->Top3    = NULL;
    tree_data->Bottom3 = NULL;

    tree_data->OutTop1    = NULL;
    tree_data->OutBottom1 = NULL;
    tree_data->OutTop2    = NULL;
    tree_data->OutBottom2 = NULL;
    tree_data->OutTop3    = NULL;
    tree_data->OutBottom3 = NULL;

    for (i = 0; i < 3; i++)
    {
        tree_data->ZeroCoupon[i] = NULL;
        tree_data->ZeroRate[i]   = NULL;
        tree_data->FwdRate[i]    = NULL;
        tree_data->TermZero[i]   = NULL;
    }

    for (i = 0; i < 6; i++)
    {
        tree_data->Aweight[i] = NULL;
    }

    for (i = 0; i < 3; i++)
    {
        tree_data->BetaTD[i] = NULL;
    }

    tree_data->QLeft    = NULL;
    tree_data->QRight   = NULL;
    tree_data->FwdShift = NULL;



    
    tree_data->TPDate  = NULL;
    tree_data->ZCenter = NULL;
    tree_data->Length  = NULL;
    tree_data->LengthJ = NULL;

    tree_data->NbEDevDates = 0;
    tree_data->EDevDate    = NULL;
    tree_data->EDevStPrice = NULL;

    tree_data->NbNmr    = 0;
    tree_data->NmrInv   = NULL;
    tree_data->LastZero = NULL;
    tree_data->Libor    = NULL;

    return;

}  /* Fix3_Cet_Tree_Reset */



/*****  Fix3_Cet_Print_Vol_File  ********************************************/
/*
*       Print debug information to a file of choice
*/
int     Fix3_Cet_Print_Vol_File(CET_OUT_DATA   *cet_out_data,
                           MKTVOL_DATA    *mktvol_data,
                           int             NbInstr,
                           int             IterIdx)


{

    int     i; 
    
    
    int     status = FAILURE; /* Error status = FAILURE initially      */
    char    ErrorMsg[MAXBUFF];
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;

    if (mktvol_data->TraceFlag == 'N') return (SUCCESS);

    strcpy (FileName, "CET.prn");

    if (IterIdx == 1)
    {
        stream = fopen (FileName, "w");
    }
    else
    {
        stream = fopen (FileName, "a");
    }
    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s!",FileName);
        DR_Error (ErrorMsg);
        goto FREE_MEM_AND_RETURN;
    } 

    
    fprintf (stream, "###################################################################################\n");
    fprintf (stream, "ITERATION:  %d\n\n",IterIdx);


    fprintf (stream, "       Instrument           Tree       BS     PrDiff   Vega     Vol     FiltVol    NxtVol \n");
    fprintf (stream, "           #                (bp)      (bp)    (Vega)   (%%)     (%%)       (%%)       (%%) \n");


    for (i = 0; i < NbInstr; i++)
    {
        fprintf (stream,
                 "%2d (%8ld/%8ld)   %8.2f  %8.2f  %6.2f  %6.2f  %8.4f  %8.4f (%8.4f)\n",
                 i+1,
                 YMDDateFromIRDate(mktvol_data->SwapSt[i]),
                 YMDDateFromIRDate(mktvol_data->SwapMat[i]),
                 (mktvol_data->VolUsed[i]) ? 10000.0 * cet_out_data->TreePrice[i]: 0.,
                 (mktvol_data->VolUsed[i]) ? 10000.0 * cet_out_data->BSPrice[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? cet_out_data->PriceDiffInVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->FOVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->VolPrev[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->FilterVol[i] : 0.,
                 100.0 * mktvol_data->Vol[i]
                 );
    } 

    fprintf (stream, "\n\n");


    status = SUCCESS;

    FREE_MEM_AND_RETURN:
 
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
        DR_Error("%s failed", "Fix3_Cet_Print_Vol_File");
    
    return (status);

}  /* Fix3_Cet_Print_Vol_File */
