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

/* CET error that gives failure. Change from 2.0 to 0.5 requested by DA Jan03
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




/*****  Cet_Main  ***********************************************************/
/*
*       Main subroutine.
*/
int     Cet_Main 
            (int            CetOutputFlag,     /* (I) On TRUE write to file */          
             T_CURVE        *t_curveProd,      /* (I) Zero curve data       */
             MKTVOL_DATA    *mktvol_dataProd,  /* (I/0) Vol  data           */
             TREE_DATA      *tree_dataProd)    /* (I) Tree data             */
{


    T_CURVE             *t_curve;      /* Structure of zero curve data      */
    MKTVOL_DATA         mktvol_data;   /* Structure of vol data             */
    TREE_DATA           tree_data;     /* Structure of tree data            */

    CET_OUT_DATA        cet_out_data;  /* Output: the "enhanced" mkt vols   */




    long           ValueDate;          /* From zero curve                   */
    int            NbInstr;            /* Nb of benchmark swaptions         */
    int            IterIdx;            /* Iteration counter                 */
    double         TimeToExp;          /* Time to exp of each b'mark option */
    double         slopeFO;            /* 1st order 1st derv in NRaphson    */
    int            i;

    int            isFirstWarning = TRUE;  /* Prints warning header info    */

    int            status = FAILURE;   /* Status */ 
    char           ErrorMsg[MAXBUFF];

    double         TtoStart;
    double         vol1, vol2;
    int            rIdx;
                 

    /* Return original vols if not enhanced */
    if ((mktvol_dataProd->CetNbIter == 0) || 
        (mktvol_dataProd->CalibFlag == FALSE))
    {
        status = SUCCESS;
        return (status);
    }

    /* Each vol point becomes a critical event type on tree_data */
    /* and therefore total nb vols cannot exceed limit           */
    if (mktvol_dataProd->NbVol > NBCRITDATE)
    {
        sprintf (ErrorMsg, "Nb of vols exceeds maximum limit "
                 "of %ld (Cet_Main)!", NBCRITDATE);                 
        DR_Error (ErrorMsg);
        goto FREE_MEM_AND_RETURN;

    }

    /*  Initialize tree data structure to NULL. */
    Tree_Init (&tree_data);

    /*  I/O manager: Initialize according to product data. */
    t_curve = t_curveProd;
    if (Cet_Manager(mktvol_dataProd,
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
        if (Par_Yield_From_Dates(&(cet_out_data.ParYield[i]),
                                 &(cet_out_data.Annuity[i]), 
                                 mktvol_data.SwapSt[i], 
                                 mktvol_data.SwapMat[i],       
                                 mktvol_data.DCC,           
                                 mktvol_data.Freq,          
                                 'F',  /* If ever a stub, it will be front */    
                                 t_curve[tree_data.CvDiff].NbZero,        
                                 t_curve[tree_data.CvDiff].Zero,         
                                 t_curve[tree_data.CvDiff].ZeroDate,     
                                 ValueDate) == FAILURE)
        {
           goto FREE_MEM_AND_RETURN;
        }

        if (Par_Yield_From_Dates(&(cet_out_data.ParYield2[i]),
                                 &(cet_out_data.Annuity2[i]), 
                                 mktvol_data.SwapSt[i], 
                                 mktvol_data.SwapMat2[i],       
                                 mktvol_data.DCC,           
                                 mktvol_data.Freq,          
                                 'F',  /* If ever a stub, it will be front */    
                                 t_curve[tree_data.CvDiff].NbZero,        
                                 t_curve[tree_data.CvDiff].Zero,         
                                 t_curve[tree_data.CvDiff].ZeroDate,     
                                 ValueDate) == FAILURE)
        {
           goto FREE_MEM_AND_RETURN;
        }

        /* Price with B&S */
        TimeToExp = Daysact(ValueDate,mktvol_data.SwapSt[i])/365.;
        if (TimeToExp < TINY)
        {
            sprintf (ErrorMsg, "B'mark swap start %ld before value date !",
                     mktvol_data.SwapSt[i]);                 
            DR_Error (ErrorMsg);
            goto FREE_MEM_AND_RETURN;
        } 

        cet_out_data.BSPrice[i] = Put_BS(cet_out_data.ParYield[i],
                                         cet_out_data.ParYield[i],
                                         TimeToExp,
                                         0.0,  /* risk free rate set to 0 */
                                         mktvol_data.TargetVol[i]);
        
        cet_out_data.BSVega[i] = Vega_BS(cet_out_data.ParYield[i],
                                         cet_out_data.ParYield[i],
                                         TimeToExp,
                                         0.0,
                                         mktvol_data.TargetVol[i]);

        /* Annuity includes the discounting */
        cet_out_data.BSPrice[i] *= cet_out_data.Annuity[i];
        cet_out_data.BSVega[i]  *= cet_out_data.Annuity[i];

    }

    /* Adjust initial vol input using CetVolShift */  
    for (i = 0; i < mktvol_dataProd->NbVol; i++)
        mktvol_data.Vol[i] *= mktvol_dataProd->CetVolShiftArray[i];

    for (IterIdx = 1; IterIdx <= mktvol_dataProd->CetNbIter; IterIdx++)
    {
        
        /* 1 - Build time line with dates of b'mark instruments.      */
        /*  Out_Data is passed so that paryields are stored alongside */
        /*  the timeline; tree_dataProd is passed in order to match   */
        /*  the product & CET timelines (for short-dated swaptions)   */
        if (Cet_Schedule(ValueDate,                                 
                         &cet_out_data,                             
                         &mktvol_data,                              
                         tree_dataProd,
                         &tree_data) == FAILURE)                    
        {                                                           
            goto FREE_MEM_AND_RETURN;                               
        }                              
                                                                    

        /*  2 - Build tree from current market vols */
        if (Build_Tree (t_curve,
                        &mktvol_data,
                        &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

        
        /*  3 - Price benchmarks on tree */
        if (Cet_Calc(&mktvol_data,                                  
                     &tree_data,
                     &cet_out_data) == FAILURE)
        {
            DR_Error ("Could not run enhancement algorithm! (main)");
            goto FREE_MEM_AND_RETURN;
        }


        /*  4 -  Update target diffs and improve vols on used vols */
        for (i=0; i<NbInstr; i++)
        {
            if (mktvol_data.VolUsed[i]) 
            {
                /* CET tree is effectively built with filtered vol */
                /* so use filtered vols in NR routine              */
                if (IndexVol(&(cet_out_data.FilterVol[i]),
                             mktvol_data.SwapSt[i],
                             mktvol_data.SwapMat[i],
                             mktvol_data.Freq,
                             mktvol_data.DCC,
                             mktvol_data.CalibFlag,
                             mktvol_data.NbVol,
                             mktvol_data.BaseDate,
                             mktvol_data.VolDate,
                             mktvol_data.Aweight, /* has been filtered  */
                             mktvol_data.QLeft,
                             mktvol_data.QRight,
                             mktvol_data.FwdShift,
                             mktvol_data.Bbq,
                             mktvol_data.VolNorm,
                             mktvol_data.VolLogn,
                             tree_data.NbFactor,
                             mktvol_data.NbMr,
                             mktvol_data.MrDate,
                             mktvol_data.BetaTD,
                             mktvol_data.BetaBmk,
                             t_curve[tree_data.CvDiff].NbZero,
                             t_curve[tree_data.CvDiff].Zero,
                             t_curve[tree_data.CvDiff].ZeroDate,
                             ValueDate) == FAILURE)
                {
                    DR_Error ("Could not calculate filtered CET vol");
                    goto FREE_MEM_AND_RETURN;
                } 

                cet_out_data.Vol[i]             = mktvol_data.Vol[i]; /* just for display in CET.prn */  
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
            if (Cet_Print_Vol_File(&cet_out_data,
                                   &mktvol_data,
                                   NbInstr,
                                   IterIdx) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            } 
        }
        
        Tree_Free(&tree_data);

        /* re-initialize pointers back to null */
        Cet_Tree_Reset(&tree_data);
    }

    /* Check final CET vols are tolerable; report error/warning as appropriate; record new vols */
    /* Do this AFTER call to Cet_Print_Vol_File for accesss to full debug info via CET.prn      */ 
    for (i=0; i<NbInstr; i++)
    {        
        if (mktvol_data.VolUsed[i])
        {
            if (fabs(cet_out_data.PriceDiffInVega[i]) > CET_ERROR_TOLERANCE)
            {
                /* print error to stdout and exit */
                DR_Error("Cet_Main Error: Final CET vol has not converged to within %4.2lf vega!\n"
                         "Expiry = %ld SwapSt = %ld SwapMat = %ld VegaDiff = %4.2lf vega\n",
                         CET_ERROR_TOLERANCE, 
                         mktvol_data.VolDate[i],
                         mktvol_data.SwapSt[i], 
                         mktvol_data.SwapMat[i],
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
                                   "Cet_Main Warning: Final CET vol has not converged to within %6.4lf vega!\n"
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

    /* If mean-reversion is calibrated,
	   calculate bmark vols and ratios with final vols */
    if ((mktvol_dataProd->MrCalibFlag == 1)||(mktvol_dataProd->MrCalibFlag == 3))
	{
        /* 1 - Build time line with dates of b'mark instruments.      */
        /*  Out_Data is passed so that paryields are stored alongside */
        /*  the timeline; tree_dataProd is passed in order to match   */
        /*  the product & CET timelines (for short-dated swaptions)   */
        if (Cet_Schedule(ValueDate,                                 
                         &cet_out_data,                             
                         &mktvol_data,                              
                         tree_dataProd,
                         &tree_data) == FAILURE)                    
        {                                                           
            goto FREE_MEM_AND_RETURN;                               
        }                              
                                                                    

        /*  2 - Build tree from current market vols */
        if (Build_Tree (t_curve,
                        &mktvol_data,
                        &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

        
        /*  3 - Price benchmarks on tree */
        if (Cet_Calc(&mktvol_data,                                  
                     &tree_data,
                     &cet_out_data) == FAILURE)
        {
            DR_Error ("Could not run enhancement algorithm! (main)");
            goto FREE_MEM_AND_RETURN;
        }


        /*  4 -  Update target diffs */
        for (i=0; i<NbInstr; i++)
        {
            if (mktvol_data.VolUsed[i]) 
            {
                cet_out_data.Vol[i]             = mktvol_data.Vol[i]; /* just for display in CET.prn */  
                cet_out_data.PriceDiff[i]     = cet_out_data.BSPrice[i] 
                                              - cet_out_data.TreePrice[i];
                cet_out_data.TreePricePrev[i] = cet_out_data.TreePrice[i];

                cet_out_data.PriceDiffInVega[i] = 100. * cet_out_data.PriceDiff[i]/
                                                  cet_out_data.BSVega[i];
            }
        }

        /*  5 - Print to VOL.prn file */
        if (CetOutputFlag)
        {
            if (Cet_Print_Vol_File(&cet_out_data,
                                   &mktvol_data,
                                   NbInstr,
                                   IterIdx) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            } 
        }
        
        Tree_Free(&tree_data);

        /* re-initialize pointers back to null */
        Cet_Tree_Reset(&tree_data);
	

    /* record vol ratios if needed*/
        rIdx = 0;
        for (i = 0; i < mktvol_dataProd->NbVol; i++)
        {
            if (mktvol_dataProd->RatioFlag[i] == TRUE)
            {
                TtoStart = Daysact(ValueDate, mktvol_dataProd->SwapSt[i])/365.;
                vol2 = ImpVol_BS2Q (cet_out_data.ParYield2[i],
                                cet_out_data.ParYield2[i],
                                TtoStart,                
                                cet_out_data.TreePrice2[i]/cet_out_data.Annuity2[i],
                                1,             
                                1.,           
                                1.,           
                                0.,           
                                mktvol_dataProd->TargetVol[i]/*(mktvol_dataProd->VNFMRatio[rIdx])*/);
                if (vol2 < 0.)
                {
                    goto FREE_MEM_AND_RETURN;
                } 

                vol1 = ImpVol_BS2Q (cet_out_data.ParYield[i],
                                cet_out_data.ParYield[i],
                                TtoStart,                
                                cet_out_data.TreePrice[i]/cet_out_data.Annuity[i],
                                1,             
                                1.,           
                                1.,           
                                0.,           
                                mktvol_dataProd->TargetVol[i]);

                if (vol1 < TINY)
                {
                    goto FREE_MEM_AND_RETURN;
                }

                /* printf("\n%ld, %ld, %20.15lf", mktvol_dataProd->SwapSt[i], mktvol_dataProd->SwapMat[i], vol1);
                printf("\n%ld, %ld, %20.15lf", mktvol_dataProd->SwapSt[i], mktvol_dataProd->SwapMat2[i], vol2);
                printf("\nratio = %20.15lf", vol2/vol1); */
                //printf("\n%20.15lf, %20.15lf, %20.15lf", vol1, vol2, vol2/vol1);
            
                mktvol_dataProd->VolRatio[rIdx] = vol2/vol1;
                rIdx++;
            }
        } /* for i */

        if (rIdx != mktvol_dataProd->NbMr)
        {
            goto FREE_MEM_AND_RETURN;
        }

    } /* if */

    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Tree_Free (&tree_data);
                                               
    return (status);

}  /* Cet_Main */


/*****  Cet_Calc ************************************************************/
/*
*       Main calculation routine: discounted expected value of cash-flows
*       going backward in the tree.
*/
int     Cet_Calc
             (MKTVOL_DATA         *mktvol_data,    /* (I) Vol data          */
              TREE_DATA           *tree_data,      /* (I) Tree data         */
              CET_OUT_DATA        *cet_out_data)   /* (I) Output data       */
{


    DEV_DATA    dev_data;               /* Dev data structure             */



    /* Variables */                   
    double      *BMarkPrice[MAX_INST];  /* Prices of b'mark swaps         */
    double      *BMarkPrice2[MAX_INST];  /* Prices of b'mark swaps         */
    double      *Swaption = NULL;       /* Prices of b'mark swaptions     */
    double      *AuxSlice = NULL;       /* Used for smoothing             */

    /* Flags & amounts */                       
    long        ExerFlag[MAX_INST];     /* Indicates an exercise date     */
    long        FixCpnFlag[MAX_INST];   /* Indicates a fixed payment date */
    double      FixCpnAmt[MAX_INST];
    long        FixCpnFlag2[MAX_INST];   /* Indicates a fixed payment date */
    double      FixCpnAmt2[MAX_INST];
    int         Smoothing;              /* Use smoothing if TRUE */

    /* Numeric variables */
    long        CurrentDate;       /* Current date                        */
    long        ValueDate;         /* Value date                          */
    int         ICurve;
    int         DCurve;            /* Discount curve number               */
    int         NbInstr;           /* Nb of instruments to price          */
    int         NbVol;           /* Nb of instruments to price          */
   
    double      Aux;
    int         T;                 /* Last time point                     */
    int         t;                 /* Current time point                  */
    int         OffsetAt0;         /* Node offset at t=0                  */
    int         k;                 /* Convenience index                   */
    int         status = FAILURE;  /* Error status = FAILURE initially    */

    
    /* initialise pointers to NULL */

    Dev_Init(&dev_data);

    for (k=0; k<MAX_INST; k++) 
    {
        BMarkPrice[k] = NULL;
        BMarkPrice2[k] = NULL;
    }

    /* Total number of time points and value date */
    T   = tree_data->NbTP;        
    ValueDate = tree_data->TPDate[0];
    
    /* Number of market instruments (swaptions) to price */
    NbInstr = mktvol_data->NbVol;
    NbVol = mktvol_data->NbVol;

    /* Assigment of discount and index curves in the engine */
    ICurve = tree_data->CvDiff;
    DCurve = tree_data->CvDisc;

    /* Smoothing flag */
    Smoothing = (mktvol_data->SmoothingFlag == 'Y');

    OffsetAt0 = Node_Offset(tree_data->NbFactor, 0, 0, 0, tree_data);


    /*   Allocation of variables for tree pricing and initialisation */
    for (k=0; k<NbInstr; k++)
    {
        BMarkPrice[k] = Alloc_Slice(tree_data);
        BMarkPrice2[k] = Alloc_Slice(tree_data);  
        if ((BMarkPrice[k] == NULL)||(BMarkPrice2[k] == NULL))
        {
            DR_Error("Calc_Cet: could not allocate memory for b'mark swaps!");
            goto FREE_MEM_AND_RETURN;        
        }
    }

    Swaption = Alloc_Slice(tree_data);
    AuxSlice = Alloc_Slice(tree_data);
    if ( (Swaption == NULL) || (AuxSlice == NULL) ) 
    {
        DR_Error("Calc_Cet: could not allocate memory for b'mark swaptions!");
        goto FREE_MEM_AND_RETURN;        
    }

    if (Dev_Alloc (&dev_data, tree_data) == FAILURE)
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

            FixCpnFlag2[k]   = tree_data->TPtype[k+NbVol][t]; 
            if (FixCpnFlag2[k])
            {
                FixCpnAmt2[k]    = (tree_data->CritDate[k+NbVol][t]).Value[0];
            }
        
            ExerFlag[k] = FALSE;
            if (CurrentDate == mktvol_data->SwapSt[k])
            {
                ExerFlag[k] = TRUE;
            }
        } 


        /*  'Update' tree. */
        if (Lattice (&dev_data,
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


            if (mktvol_data->SwapSt[k]  <= CurrentDate)
            {
                if (mktvol_data->SwapMat[k] >= CurrentDate)
                {
                    /* DEV benchmark price */
                    if (Dev(BMarkPrice[k],    
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
                        if (AddScalar(BMarkPrice[k],
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
                        if (Set_Slice(Swaption,
                                  0.,
                                  t,
                                  tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }

                        /* evaluate option claim */
                        if (OptionPlus_t(Swaption,
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

                        if (MultTwoSlicesAddAll(&Aux,
                                            Swaption,
                                            tree_data->EDevStPrice[k],
                                            t,
                                            tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }

                        (BMarkPrice[k]+OffsetAt0)[0] = Aux;
                    } /* if ExerFlag */
                } /* if before final maturity of k-th benchmark (index 1) */ 

                if(mktvol_data->RatioFlag[k] == TRUE)
                {
                if (mktvol_data->SwapMat2[k] >= CurrentDate)
                {
                    /* DEV benchmark price */
                    if (Dev(BMarkPrice2[k],    
                        t,        
                        T,        
                        ICurve,   
                        &dev_data,
                        tree_data) == FAILURE)
                    {
                        goto FREE_MEM_AND_RETURN;
                    }

                    /* Add coupon (final principal included) */
                    if (FixCpnFlag2[k])
                    {
                        if (AddScalar(BMarkPrice2[k],
                                  FixCpnAmt2[k],
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
                        if (Set_Slice(Swaption,
                                  0.,
                                  t,
                                  tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }

                        /* evaluate option claim */
                        if (OptionPlus_t(Swaption,
                                     NULL,          /* Statistics not needed */
                                     NULL,          /* Statistics not needed */
                                     NULL,          /* Statistics not needed */
                                     BMarkPrice2[k],
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

                        if (MultTwoSlicesAddAll(&Aux,
                                            Swaption,
                                            tree_data->EDevStPrice[k],
                                            t,
                                            tree_data) == FAILURE)
                        {
                            goto FREE_MEM_AND_RETURN;
                        }

                        (BMarkPrice2[k]+OffsetAt0)[0] = Aux;
                    } /* If Exercise date */

                } /* If before final maturity of k-th benchmark (index 2) */

                } /* if RatioFlag */
            } /* if after SwapSt[k] */

        } /* For k */
        
    }  /* for t */



    /* Record prices. */
    
    for (k=0; k<NbInstr; k++)
    {
        cet_out_data->TreePrice[k] = (BMarkPrice[k] + OffsetAt0)[0];
        cet_out_data->TreePrice2[k] = (BMarkPrice2[k] + OffsetAt0)[0];
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:
    
    for (k=0; k<NbInstr; k++)
    {
        Free_Slice(BMarkPrice[k],tree_data);
        Free_Slice(BMarkPrice2[k],tree_data);
    }    

    Free_Slice(Swaption,tree_data);
    Free_Slice(AuxSlice,tree_data);


    Dev_Free (&dev_data, tree_data);



    return (status);

}  /* Calc_Cet */




/*****  Cet_Schedule ***********************************************/
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

int     Cet_Schedule 
            (long            ValueDate,      /* (I) Value date             */
             CET_OUT_DATA    *cet_out_data,
             MKTVOL_DATA     *mktvol_data,
             TREE_DATA       *tree_dataProd, /* (I) Tree data for product  */
             TREE_DATA       *tree_data)     /* (O) Tree data for CET      */

{


    int         NbCritDate = 0;
    CRIT_DATE   *CritDate = NULL;        /* Critical date list             */
    
    EVENT_LIST  *PmtEventListFix[NBCRITDATE]; /* Event list for fixed pmts */
    EVENT_LIST  *PmtEventListFix2[NBCRITDATE];/* Event list for fixed pmts (vol idx 2) */
    
    int         NbInstr = 0;             /* Nb of instruments to price     */

    int         i, k;                    /* Index for convenience          */
    int         status = FAILURE;        /* Status = FAILURE initially     */
    
    long        CurrCouponDate;          /* Date of prossesed coupon       */
    double      CurrCouponDayCount;      /* Day count of prossesed coupon  */
    double      CurrCouponAmt;
    int         NbVol = mktvol_data->NbVol;
    

    /* Only CET those vols that are needed to build the product tree */ 
    NbInstr = mktvol_data->NbCetVol;

    /* Check to ensure that there are enough critical dates */
    /*       for all chosen benchmark instruments           */
    if (NbInstr > NBCRITDATE)
    {
        DR_Error("The nb of swaptions to calibrate exceeds the pre-set\n"
                 "memory allocated for it. (Cet_Schedule)\n");
        goto FREE_MEM_AND_RETURN;
    }

    /* Initialise all event lists to NULL */
    for (k=0; k<NbInstr; k++) 
    {
        PmtEventListFix[k] = NULL;
        PmtEventListFix2[k] = NULL;
    }

    /* EXPRESS DEV DATE LIST */
    /* The memory for the express DEV  dates gets allocated here, */
    /* but is freed via Tree_Free if failure occurs at some point */
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
        DR_Error("Cet_Schedule: unable to allocate memory for critical "
                    "dates array !");
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
        tree_data->CritType[k + NbVol] = 'D';
        tree_data->NbZeros[k]  =  0;
        tree_data->NbZeros[k + NbVol]  =  0;


        /* Process the swap inputs to generate the appropriate fixed cp  */
        /* payment dates which will be placed in a temporary EVENT_LIST. */
        {
            long TempDates[2];  /* Only accrual start and final mat */
            long TempDates2[2];  /* Only accrual start and final mat */
            
            TempDates[0] = mktvol_data->SwapSt[k];
            TempDates[1] = mktvol_data->SwapMat[k];

            TempDates2[0] = mktvol_data->SwapSt[k];
            TempDates2[1] = mktvol_data->SwapMat2[k];
        
            PmtEventListFix[k] = DrNewEventListFromFreq 
                                  (2,
                                   TempDates,
                                   mktvol_data->Freq,
                                   'F',    /* Forward stub                  */
                                   'N',    /* 'Dates in' check not required */
                                   NULL, NULL, NULL, NULL, NULL);
            
            PmtEventListFix2[k] = DrNewEventListFromFreq 
                                  (2,
                                   TempDates2,
                                   mktvol_data->Freq,
                                   'F',    /* Forward stub                  */
                                   'N',    /* 'Dates in' check not required */
                                   NULL, NULL, NULL, NULL, NULL);
            
            if (PmtEventListFix[k] == NULL || PmtEventListFix2[k] == NULL )
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
                                                
        /* index 1 */
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

        /* index 2 */
        if (mktvol_data->RatioFlag[k] == TRUE)
        {
        for (i = 1; i < PmtEventListFix2[k]->NbEntries; i++)            
        {
           
            CurrCouponDate = PmtEventListFix2[k]->Dates[i];
        
            if (DrDayCountFraction (PmtEventListFix2[k]->Dates[i-1],
                                    PmtEventListFix2[k]->Dates[i], 
                                    mktvol_data->DCC,
                                    &CurrCouponDayCount) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        
            CurrCouponAmt = CurrCouponDayCount 
                          * cet_out_data->ParYield2[k];
            if (i == PmtEventListFix2[k]->NbEntries-1)
            {
                CurrCouponAmt += 1.0; /* Final principal */
            }
        
            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 PmtEventListFix2[k]->Dates[i],
                                 k + NbVol,
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
                             k + NbVol,                            
                             -1.0,                
                             0, 0, 0, 0,                   
                             0, 0, 0) == FAILURE)          
        {                                                  
            goto FREE_MEM_AND_RETURN;                      
        }
        }
        
        /* Express DEV date */
        tree_data->EDevDate[k] = PmtEventListFix[k]->Dates[0];

    } /* For k */




    /************************************/
    /*                                  */
    /*  TYPES NbVol TO NBCRITDATE:  N/A */
    /*                                  */ 
    /************************************/       
                    
    for (i = 2*NbVol; i < NBCRITDATE; i++)
    {
        tree_data->CritType[i] = 'D';
        tree_data->NbZeros[i]  =  0;
    }

    /* 
     *   Finally construct the time line using 'J' time steps
     *   i.e. increasing time steps with modified jump size
     */
    if (Time_Line (ValueDate,                      
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
        DrFreeEventList(PmtEventListFix2[k]);
    }


    if (status ==FAILURE)
    {
        DR_Error("Cet_Schedule: Failed.\n");
    }

    return (status);

}  /* Cet_Schedule */




/*****  Cet_Manager  ********************************************************/
/*
*       Manage input & output, reading from ascii files into data structures.
*/
int     Cet_Manager 
            (MKTVOL_DATA    *mktvol_dataProd,  /* (I) Vol  data             */
             TREE_DATA      *tree_dataProd,    /* (I) Tree data             */
             MKTVOL_DATA    *mktvol_data,      /* (O) Vol  data             */
             TREE_DATA      *tree_data)        /* (O) Tree data             */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */
    int     i;                              /* Convenience index            */   


    /*  Model specific data  */

    tree_data->NbSigmaMax = tree_dataProd->NbSigmaMax;
    mktvol_data->QRight   = mktvol_dataProd->QRight;
    mktvol_data->QLeft    = mktvol_dataProd->QLeft;
    mktvol_data->FwdShift = mktvol_dataProd->FwdShift;
    mktvol_data->Bbq      = mktvol_dataProd->Bbq;
    mktvol_data->VolNorm  = mktvol_dataProd->VolNorm;
    mktvol_data->VolLogn  = mktvol_dataProd->VolLogn;
    tree_data->NbFactor   = tree_dataProd->NbFactor;
    for (i = 0; i < 3; i++)
    {
        mktvol_data->Alpha[i] = mktvol_dataProd->Alpha[i]; 
        mktvol_data->Beta[i]  = mktvol_dataProd->Beta[i]; 
        mktvol_data->Rho[i]   = mktvol_dataProd->Rho[i];
    }

    /*  Standard environment elements */

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
        mktvol_data->SwapMat2[i] = mktvol_dataProd->SwapMat2[i];
        mktvol_data->Vol[i]     = mktvol_dataProd->Vol[i];
        mktvol_data->TargetVol[i] = mktvol_dataProd->TargetVol[i];
        mktvol_data->VolUsed[i] = mktvol_dataProd->VolUsed[i];
    }

    mktvol_data->FilterSpotVolFlag = TRUE;
    mktvol_data->SmoothingFlag = mktvol_dataProd->SmoothingFlag;
    mktvol_data->TraceFlag     = mktvol_dataProd->TraceFlag;

    tree_data->JumpPpy = tree_dataProd->Ppy;

    /* Efficiency features */

    /* Reduce PPY. Cet_Schedule sets CET = product timeline over first 2y */
    if (tree_data->NbFactor == 1)
    {
        tree_data->Ppy = MIN (tree_dataProd->Ppy, 24);
        //tree_data->Ppy = tree_dataProd->Ppy;
    }
    else if (tree_data->NbFactor == 2)
    {
        tree_data->Ppy = MIN (tree_dataProd->Ppy, 4);
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

            /* if MR calibrated, make sure that the last ratio date is covered */
			if ((mktvol_dataProd->MrCalibFlag == 1)||(mktvol_dataProd->MrCalibFlag == 3)||
				(mktvol_dataProd->MrCalibFlag == 5)||(mktvol_dataProd->MrCalibFlag == 7))
			{
			    if (mktvol_dataProd->LastRatioIdx > LastVolUsedIdx)
                    LastVolUsedIdx = mktvol_dataProd->LastRatioIdx;
			}

            mktvol_data->NbCetVol = LastVolUsedIdx+1;

            for (i = mktvol_data->NbCetVol; i < mktvol_data->NbVol; i++)
            {
                mktvol_data->VolUsed[i] = FALSE;
            }
        }
    }

    /* mean-reversion info */
    mktvol_data->MrVNFM = mktvol_dataProd->MrVNFM; 
    mktvol_data->NbMr = mktvol_dataProd->NbMr;
    mktvol_data->MrCalibFlag = mktvol_dataProd->MrCalibFlag;

    for (i = 0; i < mktvol_dataProd->NbMr; i++)
    {
        mktvol_data->MrDate[i] = mktvol_dataProd->MrDate[i];
        mktvol_data->MrInput[i] = mktvol_dataProd->MrInput[i];
    }

    for (i = 0; i < mktvol_dataProd->NbVol; i++)
    {
        mktvol_data->RatioFlag[i] = mktvol_dataProd->RatioFlag[i];
    }
    
    status = SUCCESS;
        
    return (status);

}  /* Cet_Manager */



/*****  Cet_Tree_Reset  ******************************************************/
/*
*       Re-initialize tree pointers to NULL.
*/
void    Cet_Tree_Reset(TREE_DATA *tree_data) /* Tree building data structure */
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
        tree_data->BetaTD[i] = NULL;
        tree_data->BetaInt[i] = NULL;
    }

    for (i = 0; i < 6; i++)
    {
        tree_data->Aweight[i] = NULL;
    }

    tree_data->TPDate  = NULL;
    tree_data->ZCenter = NULL;
    tree_data->Length  = NULL;
    tree_data->LengthJ = NULL;

    tree_data->NbEDevDates = 0;
    tree_data->EDevDate    = NULL;
    tree_data->EDevStPrice = NULL;

    return;

}  /* Cet_Tree_Reset */



/*****  Cet_Print_Vol_File  ********************************************/
/*
*       Print debug information to a file of choice
*/
int     Cet_Print_Vol_File(CET_OUT_DATA   *cet_out_data,
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
        sprintf (ErrorMsg, "Could not open file %s! (Print_Vol)",FileName);                 
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
                 mktvol_data->SwapSt[i],
                 mktvol_data->SwapMat[i],
                 (mktvol_data->VolUsed[i]) ? 10000.0 * cet_out_data->TreePrice[i]: 0.,
                 (mktvol_data->VolUsed[i]) ? 10000.0 * cet_out_data->BSPrice[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? cet_out_data->PriceDiffInVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->FOVega[i] : 0.,
                 (mktvol_data->VolUsed[i]) ? 100.0 * cet_out_data->Vol[i] : 0.,
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
    
    return (status);

}  /* Cet_Print_Vol_File */
