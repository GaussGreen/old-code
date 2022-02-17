/****************************************************************************/
/*      Calibration Enhancement Tool                                        */
/****************************************************************************/
/*      CET.C                                                               */
/****************************************************************************/

/*
$Header$
*/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cupslib.h"

        
/*****  Hyb3_Cet_Main  ***********************************************************/
/*
*       Main subroutine.
*/
int     Hyb3_Cet_Main 
            (int           CetOutputFlag,    /* (I) On TRUE write to file   */ 
             int           CcyIdx,           /* (I) 0,1 for ccy being CET'd */
             T_CURVE       t_curveProd[2][3],/* (I) 2 ccy's,three ZC's each */
             MKTVOL_DATA  *mktvol_dataProd,  /* (I/0) Vol  data             */  
             HYB3_TREE_DATA    *tree_dataProd)    /* (I) Tree data               */
{


    T_CURVE             t_curve[2][3]; /* Structure of zero curve data      */
    MKTVOL_DATA         mktvol_data;   /* Structure of vol data             */
    HYB3_TREE_DATA           tree_data;     /* Structure of tree data            */

    CET_OUT_DATA        cet_out_data;  /* Output: the "enhanced" mkt vols   */




    long           ValueDate;          /* From zero curve                   */
    int            NbInstr;            /* Nb of benchmark swaptions         */
    int            IterIdx;            /* Iteration counter                 */
    double         TimeToExp;          /* Time to exp of each b'mark option */
    double         slopeFO;            /* 1st derivative to use in NRaphson */
    double         MaxDiff=0.;         /* Max price diff between tree & BS  */
    int            i;

    int            status = FAILURE;   /*             Status                */ 

    /* Return original vols if not enhanced */
    if ((mktvol_dataProd[CcyIdx].CetNbIter == 0) || 
        (mktvol_dataProd[CcyIdx].CalibFlag == FALSE))
    {
        status = SUCCESS;
        return (status);
    }

    /* check that nb of vols is within limit */
    if (mktvol_dataProd[CcyIdx].NbVol > NBCRITDATE)
    {
        DR_Error("Nb of vols exceeds maximum limit of %d (Hyb3_Cet_Main)!", NBCRITDATE);                 
        goto FREE_MEM_AND_RETURN;

    }
    /*  Initialize local tree data structure to NULL. */
    switch (mktvol_dataProd[CcyIdx].NbFactor)
    {
    case 1:
        tree_data.TreeType = TTYPE_1IR;
        break;
    case 2:
        tree_data.TreeType = TTYPE_1IR2F;
        break;
    default:
        DR_Error("CET not implemented for "
                    "more than 2 factors! (Hyb3_Cet_Main)\n");
        goto FREE_MEM_AND_RETURN;
    }


    Hyb3_Tree_Init (&tree_data);


    /*  I/O manager: Initialize according to product data. */
    ValueDate = (t_curveProd[CcyIdx][tree_dataProd->CvDiff[CcyIdx]]).ValueDate;
    NbInstr = mktvol_dataProd[CcyIdx].NbVol;
    t_curve[0][0] = t_curveProd[CcyIdx][0]; /* For local use */
    t_curve[0][1] = t_curveProd[CcyIdx][1];
    t_curve[0][2] = t_curveProd[CcyIdx][2];
    
    if (Hyb3_Cet_Manager(CcyIdx,
                    mktvol_dataProd,   
                    tree_dataProd,
                    &mktvol_data,
                    &tree_data) == FAILURE)
    {                                      
        goto FREE_MEM_AND_RETURN;
    }


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
                              &t_curve[0][tree_dataProd->CvDiff[CcyIdx]]) != SUCCESS)
        {
           goto FREE_MEM_AND_RETURN;
        }

        /* Price with B&S */
        TimeToExp = Daysact(ValueDate,mktvol_data.SwapSt[i])/365.;
        if (TimeToExp < TINY)
        {
            DR_Error("B'mark swap start %ld before value date !",
                     YMDDateFromIRDate(mktvol_data.SwapSt[i]));                 
            goto FREE_MEM_AND_RETURN;
        } 

        cet_out_data.BSPrice[i] = Put_BSAS(cet_out_data.ParYield[i],
                                           cet_out_data.ParYield[i],
                                           TimeToExp,
                                           0.0,  /* risk free rate set to 0 */
                                           mktvol_data.Vol[i]);

        cet_out_data.BSVega[i] = Vega_BS(cet_out_data.ParYield[i],
                                         cet_out_data.ParYield[i],
                                         TimeToExp,
                                         0.0,
                                         mktvol_data.Vol[i]);

        /* Annuity includes the discounting */
        cet_out_data.BSPrice[i] *= cet_out_data.Annuity[i];
        cet_out_data.BSVega[i]  *= cet_out_data.Annuity[i];

    }


    /* Iterate on mkt vols for convergence on prices */
    for (IterIdx = 1; IterIdx <= mktvol_dataProd[CcyIdx].CetNbIter; IterIdx++)
    {
        /*******************************************************************/
        /* need to re-initialise SigmaDBL because CET calls tree_init which*/
        /* sets SigmaDBL to -999                                           */
        /*******************************************************************/

        tree_data.NbSigmaDBL = tree_dataProd->NbSigmaDBL;

        /* 1 - Build time line with dates of b'mark instruments */  
        /* (Out_Data is  passed so that  paryields are stored   */  
        /*                alongside the timeline)               */  
        if (Hyb3_Cet_Schedule(ValueDate,   
                         &cet_out_data,                             
                         &mktvol_data,                              
                         &tree_data) == FAILURE)                    
        {                                                           
            goto FREE_MEM_AND_RETURN;                               
        }                              


        /* Recall that local nb CET iterations has been set to 0 */
        /* so that CET is not called recursively                 */         

        /*  2 - Build tree from current market vols */
        if (Hyb3_Build_Tree (FALSE,     /* CetOutputFlag */
                        t_curve,   /* Struct with relevant 3 ZC's as first */
                        &mktvol_data,
                        NULL,  /* No FX data */
                        NULL,  /* No Eq data */
                        &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

        
        /*  3 - Price benchmarks on tree */
        if (Hyb3_Cet_Calc(&mktvol_data,
                     &tree_data,
                     &cet_out_data) == FAILURE)
        {
            DR_Error ("Could not run enhancement algorithm! (main)");
            goto FREE_MEM_AND_RETURN;
        }


        /*  4 -  Update target diffs and improve vols on used vols */
        MaxDiff = 0.;
        for (i=0; i<NbInstr; i++)
        {
            if (mktvol_data.VolUsed[i]) 
            {
                cet_out_data.PriceDiff[i] = cet_out_data.BSPrice[i] 
                                          - cet_out_data.TreePrice[i];


                slopeFO  = cet_out_data.TreePrice[i];
                slopeFO /= mktvol_data.Vol[i];
            
                cet_out_data.FOVega[i] = slopeFO;
        
                cet_out_data.VolPrev[i]       = mktvol_data.Vol[i];
                mktvol_data.Vol[i]           += cet_out_data.PriceDiff[i]/slopeFO;
                cet_out_data.TreePricePrev[i] = cet_out_data.TreePrice[i];
                cet_out_data.PriceDiffInVega[i]  = 100. * cet_out_data.PriceDiff[i]/
                                                 cet_out_data.BSVega[i];
                MaxDiff = MAX(fabs(cet_out_data.PriceDiffInVega[i]),MaxDiff);
            }
        }


        /*  5 - Print to VOL.prn file */
        if (CetOutputFlag)
        {
            if (Hyb3_Cet_Print_Vol_File(CcyIdx,
                                   &cet_out_data,
                                   &mktvol_data,
                                   NbInstr,
                                   IterIdx) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            } 
        }
        

        Hyb3_Tree_Free(&tree_data);
        Hyb3_Tree_Init(&tree_data);
    }

    /* Record new vols and its difference from market    */
    mktvol_dataProd[CcyIdx].CetVegaError = MaxDiff;
    for (i=0; i<NbInstr; i++)
    {
        mktvol_dataProd[CcyIdx].Vol[i]     = mktvol_data.Vol[i];
        mktvol_dataProd[CcyIdx].VolUsed[i] = mktvol_data.VolUsed[i];
    }

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Hyb3_Tree_Free (&tree_data);

    return (status);

}  /* Hyb3_Cet_Main */





/*****  Hyb3_Cet_Manager  ********************************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Hyb3_Cet_Manager 
            (int             CcyIdx,
             MKTVOL_DATA    *mktvol_dataProd,  /* (I) Vol  data             */
             HYB3_TREE_DATA *tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA    *mktvol_data,      /* (O) Vol  data             */
             HYB3_TREE_DATA *tree_data)        /* (O) Tree data             */
{




    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   



    /*  Model specific data  */
    tree_data->NbSigmaMax = tree_dataProd->NbSigmaMax;

    tree_data->Ppy        = tree_dataProd->Ppy;
    mktvol_data->QRight   = mktvol_dataProd[CcyIdx].QRight;
    mktvol_data->QLeft    = mktvol_dataProd[CcyIdx].QLeft;
    mktvol_data->FwdShift = mktvol_dataProd[CcyIdx].FwdShift;
    mktvol_data->NbFactor = mktvol_dataProd[CcyIdx].NbFactor;
    for (i = 0; i < 3; i++)
    {
        mktvol_data->Alpha[i] = mktvol_dataProd[CcyIdx].Alpha[i]; 
        mktvol_data->Beta[i]  = mktvol_dataProd[CcyIdx].Beta[i]; 
        mktvol_data->Rho[i]   = mktvol_dataProd[CcyIdx].Rho[i];
    }

    /*  Standard environment elements */

    tree_data->CvDiff[0]      = tree_dataProd->CvDiff[CcyIdx];
    tree_data->CvIdx1[0]      = tree_dataProd->CvIdx1[CcyIdx];
    tree_data->CvIdx2[0]      = tree_dataProd->CvIdx2[CcyIdx];
    tree_data->CvDisc[0]      = tree_dataProd->CvDisc[CcyIdx];
    mktvol_data->CalibFlag = mktvol_dataProd[CcyIdx].CalibFlag;  
    mktvol_data->SkipFlag  = mktvol_dataProd[CcyIdx].SkipFlag;
    mktvol_data->BaseDate  = mktvol_dataProd[CcyIdx].BaseDate;
    mktvol_data->Freq      = mktvol_dataProd[CcyIdx].Freq;
    mktvol_data->DCC       = mktvol_dataProd[CcyIdx].DCC;
    mktvol_data->NbVol     = mktvol_dataProd[CcyIdx].NbVol;
    mktvol_data->NbCetVol  = mktvol_dataProd[CcyIdx].NbVol;

    for (i = 0; i < mktvol_data->NbVol ; i++)
    {                                                                       
        mktvol_data->VolDate[i] = mktvol_dataProd[CcyIdx].VolDate[i];
        mktvol_data->SwapSt[i]  = mktvol_dataProd[CcyIdx].SwapSt[i];
        mktvol_data->SwapMat[i] = mktvol_dataProd[CcyIdx].SwapMat[i];
        mktvol_data->Vol[i]     = mktvol_dataProd[CcyIdx].Vol[i];
        mktvol_data->VolUsed[i] = mktvol_dataProd[CcyIdx].VolUsed[i];
    }


    /* Recall that local nb CET iterations has been set to 0 */
    /* so that CET is not called recursively                 */
    mktvol_data->CetNbIter = 0;



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



    status = SUCCESS;
        
    return (status);

}  /* Hyb3_Cet_Manager */




/*****  Hyb3_Cet_Calc ************************************************************/
/*
*       Main calculation routine: discounted expected value of cash-flows
*       going backward in the tree.
*/
int     Hyb3_Cet_Calc
             (MKTVOL_DATA         *mktvol_data,    /* (I) Vol data          */
              HYB3_TREE_DATA      *tree_data,      /* (I) Tree data         */
              CET_OUT_DATA        *cet_out_data)   /* (I) Output data       */
{


    HYB3_DEV_DATA    dev_data;               /* Hyb3_Dev data structure             */



    /* Variables */                   
    double      *BMarkPrice[MAX_INST];  /* Prices of b'mark swaptions     */


    /* Flags & amounts */                       
    long        FixCpnFlag[MAX_INST];   /* Indicates a fixed payment date */
    double      FixCpnAmt[MAX_INST];
    

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
    
    
    


    /* initliase Hyb3_Dev structure */
    Hyb3_Dev_Init(&dev_data);


    /* initialise memory for BMarkPrices */
    for (k = 0; k < MAX_INST; k++)
    {
        BMarkPrice[k] = NULL;
    }



    /* Total number of time points and value date */
    T   = tree_data->NbTP;        
    ValueDate = tree_data->TPDate[0];
    
    /* Number of market instruments (swaptions) to price */
    NbInstr = mktvol_data->NbVol;

    /* Assigment of discount and index curves in the engine */
    ICurve = tree_data->CvDiff[0]; /* Running in single IR mode so */
    DCurve = tree_data->CvDisc[0]; /* only data in TREE is [0]     */

    

    OffsetAt0 = Hyb3_Node_Offset(mktvol_data->NbFactor, 0, 0, 0, tree_data); 


    /*   Allocation of variables for tree pricing and initialisation */
    for (k = 0; k < NbInstr; k++)
    {

        BMarkPrice[k] = Hyb3_Alloc_Slice(tree_data,mktvol_data->NbFactor); 
        if (BMarkPrice[k] == NULL)
        {
            DR_Error("Calc_Cet: could not allocate memory for b'mark prices!");
            goto FREE_MEM_AND_RETURN;
        }
    }


    if (Hyb3_Dev_Alloc (&dev_data, tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    
    /* One single walk back on the tree for all instruments */
    for (t = T; t >= 0; t--)
    {


        CurrentDate = tree_data->TPDate[t];

        /* Set flags for coupon pmts (final principal included) */
        for (k=0; k<NbInstr; k++)
        {
            FixCpnFlag[k] = tree_data->TPType[k][t]; 
            if (FixCpnFlag[k])
            {
                FixCpnAmt[k] = (tree_data->CritDate[k][t]).Value[0]; 
            }

        } 
      
      

        /*  'Update' tree. */
        if (Hyb3_Lattice (&dev_data,
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
                if (Hyb3_Dev(BMarkPrice[k],    
                        t,        
                        T,        
                        ICurve,   
                        (mktvol_data->NbFactor ==1)? DISC_1D_NOCUPS:
                                                    DISC_2D_1IR2F_NOCUPS,
                        &dev_data,
                        tree_data) == FAILURE)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* Add coupon (final principal included) */
                if (FixCpnFlag[k])
                {
                    if (Hyb3_AddScalar(BMarkPrice[k],
                                  mktvol_data->NbFactor,
                                  FixCpnAmt[k],
                                  t,
                                  tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }
                }

                /* If Exercise date, evaluate pay-off (par strike) */
                if (CurrentDate == mktvol_data->SwapSt[k])
                {

                    if (Hyb3_MinOnSlice(BMarkPrice[k],
                                   mktvol_data->NbFactor,
                                   0.0,
                                   t,
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

                    if (Hyb3_MultTwoSlicesAddAll(&Aux,
                                            mktvol_data->NbFactor,
                                            BMarkPrice[k],
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
    
    for (k = 0; k < NbInstr; k++)
    {
        Hyb3_Free_Slice(BMarkPrice[k],tree_data, mktvol_data->NbFactor);
    }


    Hyb3_Dev_Free (&dev_data, tree_data);

    return (status);

}  /* Calc_Cet */




/*****  Hyb3_Cet_Schedule ***********************************************/
/*
*       Sets up the time line. 
*
*       For the calibration tool, the dates are all the payment
*       dates of the  instruments (swap or money market) chosen
*       as benchmarks for vol.
*          
*
*       This routine deals with the detailed date generation from input
*       parameters and the subsequent addition of dates to the timeline
*       critical date list.
*
*/

int     Hyb3_Cet_Schedule 
            (long            ValueDate,      /* (I) Value date             */
             CET_OUT_DATA    *cet_out_data,
             MKTVOL_DATA     *mktvol_data,
             HYB3_TREE_DATA  *tree_data)     /* (O) Structure of tree data */
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





    /* Initialise all event lists to NULL */
    for (k = 0; k < NBCRITDATE; k++)
    {
        PmtEventListFix[k] = NULL;
    }


    /* Check to ensure that there are enough critical dates */
    /*       for all chosen benchmark instruments           */
    if (mktvol_data->NbVol > NBCRITDATE)
    {
        DR_Error("The nb of swaptions to calibrate exceeds the pre-set\n"
                 "memory allocated for it. (Hyb3_Cet_Schedule)\n");
        goto FREE_MEM_AND_RETURN;
    }

    NbInstr = mktvol_data->NbCetVol;


    /* CRITICAL DATE LIST */
    /*  Allocate empty critical date list and then add dates successively.*/
    CritDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    if (CritDate == NULL)
    {
        DR_Error("Hyb3_Cet_Schedule: unable to allocate memory for critical "
                    "dates array !");
        goto FREE_MEM_AND_RETURN;
    }


    /* EXPRESS DEV DATE LIST */
    /* The memory for the express DEV  dates gets allocated here, */
    /* but is freed via Hyb3_Tree_Free if failure occurs at some point */
    tree_data->NbEDevDates = NbInstr;
    tree_data->EDevDate = (long *)DR_Array(LONG, 0, NbInstr-1);
    if (tree_data->EDevDate == NULL)
    {
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
     *   Finally construct the time line using 'I'ncreasing time steps.
     */
    if (Hyb3_Time_Line (ValueDate,                      
                   NbCritDate,
                   CritDate,
                   'I',                            
                   tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    


    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Free_DR_Array (CritDate, CRITDATE, 0, NbCritDate-1);
    
    for (k=0; k < NbInstr; k++)
    {
        DrFreeEventList(PmtEventListFix[k]);
    }
   
        
    if (status == FAILURE)
    {
        DR_Error("Hyb3_Cet_Schedule: Failed.\n");
    }
    
    return (status);

}  /* Hyb3_Cet_Schedule */




/*****  Hyb3_Cet_Print_Vol_File  ********************************************/
/*
*       Print debug information to a file of choice
*/
int     Hyb3_Cet_Print_Vol_File(int        CcyIdx,
                           CET_OUT_DATA   *cet_out_data,
                           MKTVOL_DATA    *mktvol_data,
                           int             NbInstr,
                           int             IterIdx)


{

    int     i; 
    
    
    int     status = FAILURE; /* Error status = FAILURE initially      */
    char    FileName[MAXBUFF];
    FILE    *stream = NULL;



    sprintf (FileName, "CET%d.prn", CcyIdx);
    

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
        DR_Error ("Could not open file %s! (Print_Vol)",FileName);                 
        goto FREE_MEM_AND_RETURN;
    } 

    
    fprintf (stream, "###################################################################################\n");
    fprintf (stream, "ITERATION:  %d\n\n",IterIdx);


    fprintf (stream, "       Instrument           Tree Price      BS Price     Price Diff   Vega    Vol      NxtVol\n");
    fprintf (stream, "           #                   (bp)           (bp)         (Vega)      (%%)    (%%)        (%%)\n");

    for (i = 0; i < NbInstr; i++)
    {
       

        fprintf (stream,
                 "%2d (%8ld/%8ld)      %10.6f     %10.6f   %10.6f  %6.2f  %8.4f   (%5.2f)\n",
                 i+1,
                 YMDDateFromIRDate(mktvol_data->SwapSt[i]),
                 YMDDateFromIRDate(mktvol_data->SwapMat[i]),
                 10000.0 *cet_out_data->TreePrice[i],
                 (mktvol_data->VolUsed[i]) ? 10000.0 *cet_out_data->BSPrice[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? cet_out_data->PriceDiffInVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100 * cet_out_data->FOVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->VolPrev[i] : 0.,
                 100.0 * mktvol_data->Vol[i]);
    } 

    fprintf (stream, "\n\n");

    
                   

    status = SUCCESS;

    FREE_MEM_AND_RETURN:
 
    if (stream != NULL)
    {
        fclose (stream);
    }
    
    return (status);

}  /* Hyb3_Cet_Print_Vol_File */



