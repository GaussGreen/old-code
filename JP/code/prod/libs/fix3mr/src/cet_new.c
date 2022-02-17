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
#define  CET_ERROR_TOLERANCE  0.05        /* in vegas */
#endif

#ifndef  CET_WARNING_TOLERANCE
#define  CET_WARNING_TOLERANCE  0.01      /* in vegas */
#endif

#ifndef  MASKLENGTH
#define  MASKLENGTH   24L   /* in months */
#endif

#define  MAX_NR_INIT  1
#define  MIN_NR_JUMP -0.95
#define  MAX_NR_JUMP  0.95



/*****  Cet_Main_New  ***********************************************************/
/*
*       Main subroutine.
*/
int     Cet_Main_New 
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

    double         prevNRstep[MAXNBDATE];
    double         NRstep;

    int            isFirstWarning = TRUE;  /* Prints warning header info    */

    int            status = FAILURE;   /* Status */ 
    char           ErrorMsg[MAXBUFF];

    double         TtoStart;
    double         vol1, vol2;
    int            rIdx;
	double         prevVol[MAXNBDATE];
	double         prevPrice[MAXNBDATE];

	double memVol;
	int    j, volIdx;
                 

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

        /* initialize old NR step */
        prevNRstep [i] = MAX_NR_INIT;
    }

    /* Adjust initial vol input using CetVolShift */  
    for (i = 0; i < mktvol_dataProd->NbVol; i++)
        mktvol_data.Vol[i] *= mktvol_dataProd->CetVolShiftArray[i];

	/* Initialize prevVol and prevPrice to zeros */
	for (i = 0; i < mktvol_dataProd->NbVol; i++)
	{
		prevVol[i]   = 0.;
		prevPrice[i] = 0.;
	}

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
				TtoStart = Daysact(ValueDate, mktvol_dataProd->SwapSt[i])/365.;
				if ((IterIdx > 1) && 
					((cet_out_data.BSPrice[i] - prevPrice[i])*(cet_out_data.BSPrice[i] - cet_out_data.TreePrice[i]) < -TINY) &&
					 (3.*fabs(cet_out_data.BSPrice[i] - cet_out_data.TreePrice[i]) > fabs(cet_out_data.BSPrice[i] - prevPrice[i])))
				{
					NRstep = (prevVol[i] - cet_out_data.FilterVol[i])*0.5;
					slopeFO = 0.;
					cet_out_data.FOVega[i] = slopeFO;
				}
				else
				{
                    slopeFO  = cet_out_data.TreePrice[i];
                    slopeFO /= cet_out_data.FilterVol[i];
                    cet_out_data.FOVega[i] = slopeFO;

					NRstep = cet_out_data.PriceDiff[i]/slopeFO;
					NRstep = MAX(MIN_NR_JUMP * cet_out_data.FilterVol[i],NRstep);
				}
					

                /* perform NR iteration */
                mktvol_data.Vol[i] = cet_out_data.FilterVol[i] + NRstep;

                prevNRstep[i] = fabs (NRstep);

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

				/* record prevPrice and prevVol */
				prevVol[i] = cet_out_data.FilterVol[i];
				prevPrice[i] = cet_out_data.TreePrice[i];
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
    if ((mktvol_dataProd->MrCalibFlag == 5)||(mktvol_dataProd->MrCalibFlag == 7))
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


