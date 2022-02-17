/****************************************************************************/
/*      Calibration Enhancement Tool - MR                                   */
/****************************************************************************/
/*      CET.C                                                               */
/****************************************************************************/



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix123head.h"

#ifndef  MASKLENGTH
#define  MASKLENGTH   24L   /* in months */
#endif

#ifndef  MIN_CET_MR
#define  MIN_CET_MR   0.0010   
#endif

#ifndef  MAX_CET_MR
#define  MAX_CET_MR   0.5   
#endif



/*****  Cet_MR_Main  ***********************************************************/
/*
*       Main subroutine.
*/
int     Cet_MR_Main 
            (int            CetOutputFlag,     /* (I) On TRUE write to file */          
             T_CURVE        *t_curveProd,      /* (I) Zero curve data       */
             MKTVOL_DATA    *mktvol_dataProd,  /* (I/0) Vol  data           */
             TREE_DATA      *tree_dataProd)    /* (I) Tree data             */
{
    int i, j;
    long *RatioDate = mktvol_dataProd->RatioDate;
    int *RatioFlag = mktvol_dataProd->RatioFlag;
    long *SwapSt = mktvol_dataProd->SwapSt;
    long *SwapMat = mktvol_dataProd->SwapMat;
    long *SwapMat2 = mktvol_dataProd->SwapMat2;
    int NbVol = mktvol_dataProd->NbVol;
    int NbMr = mktvol_dataProd->NbMr;
    double memVol[MAX_INST];
    int rIdx;
    int currIter;
    int nrDim; 
    double B1[3];
    double B2[3];
    double *targetRatio = NULL;
    double *ratioL = NULL;
    double *errorL = NULL;
    double *mrv = NULL;
    double *mrMem = NULL;
    double **JM = NULL;
    double **JMInv = NULL;
    int exitFlag;
    double myErrTol = 0.0001;
    double delMr = 0.0001;
    double tmpSum;
    double BetaL[3];
    char   FileName[MAXBUFF];
    FILE    *stream = NULL;
    int status = FAILURE;

    strcpy (FileName, "MNR.prn");

    stream = fopen (FileName, "w");

    /* store target vols */
    for (i = 0; i < NbVol; i++)
        mktvol_dataProd->TargetVol[i] = mktvol_dataProd->Vol[i]; 

	/* calibration flag: MrCalibFlag
	   Original CET algorithm:
	   0 - FIX3 w/constant MR;
	   1 - input VNFM MR, calibrate MR to VNFM ratios;
	   2 - FIX3 w/ input MR term structure;
	   3 - (not for production) calibrate MR to given vol ratios
	   Enhanced CET algorithm:
	   4 - FIX3 w/constant MR;
	   5 - input VNFM MR, calibrate MR to VNFM ratios;
	   6 - FIX3 w/ input MR term structure;
	   7 - (not for production) calibrate MR to given vol ratios*/

    
    if ((mktvol_dataProd->MrCalibFlag == 0)||(mktvol_dataProd->MrCalibFlag == 2)||
		((mktvol_dataProd->MrCalibFlag == 1)||(mktvol_dataProd->MrCalibFlag == 3))&&(mktvol_dataProd->NbMrCet == 0))
    {
        /* initialize all flags to false */
        for (i = 0; i < NbVol; i++)
        {
            RatioFlag[i] = FALSE;
        }

        /*
        *   CET 
        */
        if (Cet_Main (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
        {
            goto RETURN;
        }

        status = SUCCESS;

        goto RETURN;
    }

	if ((mktvol_dataProd->MrCalibFlag == 4)||(mktvol_dataProd->MrCalibFlag == 6)||
		((mktvol_dataProd->MrCalibFlag == 5)||(mktvol_dataProd->MrCalibFlag == 7))&&(mktvol_dataProd->NbMrCet == 0))
    {
        /* initialize all flags to false */
        for (i = 0; i < NbVol; i++)
        {
            RatioFlag[i] = FALSE;
        }

        /*
        *   CET 
        */
        if (Cet_Main_New (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
        {
            goto RETURN;
        }

        status = SUCCESS;

        goto RETURN;
    }
   
    /* allocate memory */
    nrDim = NbMr;
    targetRatio = DR_Array(DOUBLE, 0, nrDim - 1);
    ratioL = DR_Array(DOUBLE, 0, nrDim - 1);
    errorL = DR_Array(DOUBLE, 0, nrDim - 1);
    mrv = DR_Array(DOUBLE, 0, nrDim - 1);
    mrMem = DR_Array(DOUBLE, 0, nrDim - 1);
    JM = DR_Matrix(DOUBLE, 0, nrDim - 1, 0, nrDim - 1);
    JMInv = DR_Matrix(DOUBLE, 0, nrDim - 1, 0, nrDim - 1);
    if ((targetRatio == NULL) ||
        (ratioL == NULL) ||
        (errorL == NULL) ||
        (mrv == NULL) ||
        (mrMem == NULL) ||
        (JM == NULL) || 
        (JMInv == NULL))
    {
        goto RETURN;
    }


    /* set beta for VNFM ratios */
    BetaL[0] = mktvol_dataProd->MrVNFM;

    /* set ratio flags and (if needed) generate VNFM vol ratios */
    for (i = 0; i < NbVol; i++)
        RatioFlag[i] = FALSE;


    //printf("\nB-coeff:");
    rIdx = 0;
    i = 0;
    while (rIdx < NbMr)
    {
        while(RatioDate[rIdx] != SwapSt[i])
        {
            i++;
            if (i == NbVol)
            {
                goto RETURN;
            }
        }

        RatioFlag[i] = TRUE;
        mktvol_dataProd->LastRatioIdx = i;

        if ((mktvol_dataProd->MrCalibFlag == 1)||(mktvol_dataProd->MrCalibFlag == 5)) /* compute VNFM ratios */
        {
            if (BFactor (   &(B1[0]),
                        mktvol_dataProd->SwapSt[i],
                        mktvol_dataProd->SwapMat[i],
                        mktvol_dataProd->DCC,
                        mktvol_dataProd->Freq,
                        tree_dataProd->NbFactor,
                        BetaL,
                        mktvol_dataProd->Bbq,
                        mktvol_dataProd->VolNorm,
                        mktvol_dataProd->VolLogn,
                        (t_curveProd[tree_dataProd->CvDiff]).NbZero,
                        (t_curveProd[tree_dataProd->CvDiff]).Zero,
                        (t_curveProd[tree_dataProd->CvDiff]).ZeroDate,
                        (t_curveProd[tree_dataProd->CvDiff]).ValueDate) == FAILURE)
            {
                goto RETURN;
            }

            if (fabs(B1[0]) < TINY)
            {
                goto RETURN;
            }

            if (BFactor (   &(B2[0]),
                            mktvol_dataProd->SwapSt[i],
                            mktvol_dataProd->SwapMat2[i],
                            mktvol_dataProd->DCC,
                            mktvol_dataProd->Freq,
                            tree_dataProd->NbFactor,
                            BetaL,
                            mktvol_dataProd->Bbq,
                            mktvol_dataProd->VolNorm,
                            mktvol_dataProd->VolLogn,
                            (t_curveProd[tree_dataProd->CvDiff]).NbZero,
                            (t_curveProd[tree_dataProd->CvDiff]).Zero,
                            (t_curveProd[tree_dataProd->CvDiff]).ZeroDate,
                            (t_curveProd[tree_dataProd->CvDiff]).ValueDate) == FAILURE)
            {
                goto RETURN;
            }

            mktvol_dataProd->VNFMRatio[rIdx] = B2[0]/B1[0];
            targetRatio[rIdx] = B2[0]/B1[0];
            //printf("\ni = %d, B1 = %20.15lf, B2 = %20.15lf, ratio = %20.15lf", i, B1[0], B2[0], B2[0]/B1[0]);
            //printf("\n%20.15lf, %20.15lf, %20.15lf", B1[0], B2[0], B2[0]/B1[0]);
            
            fprintf(stream, "\ntargetRatio[%d] = %20.15lf", rIdx, B2[0]/B1[0]); 
			//printf("\ntargetRatio[%d] = %20.15lf", rIdx, B2[0]/B1[0]); 
        }
        else /* use input ratios */
        {
            mktvol_dataProd->VNFMRatio[rIdx] = mktvol_dataProd->VolRatioInput[rIdx];
            targetRatio[rIdx] = mktvol_dataProd->VolRatioInput[rIdx];
        } /* if MrCalibFlag */

        rIdx++;
    } /* while rIdx */

    /* store input vols */
    for (i = 0; i < NbVol; i++)
    {
        memVol[i] = mktvol_dataProd->Vol[i];
    }

    /* initialize mean-reversions to initial guess */
    for (i = 0; i < nrDim; i++)
    {
        mrv[i] = mktvol_dataProd->MrGuess;
    }

    /* Newton-Raphson */
    currIter = 0;
    while (currIter < mktvol_dataProd->NbMrCet)
    {
        fprintf(stream, "\nIter %d:", currIter);
        fprintf(stream, "\nMean-reversions: \n");
        for (i = 0; i < nrDim; i++)
        {
            fprintf(stream, "%20.15lf,", mrv[i]);
        }

		/* printf("\nIter %d:", currIter);
        printf("\nMean-reversions: \n");
        for (i = 0; i < nrDim; i++)
        {
            printf("%20.15lf,", mrv[i]);
        } */

        /* Initial point */
        for (i = 0; i < nrDim; i++)
        {
            mrMem[i] = mrv[i];
            mktvol_dataProd->BetaTD[0][i] = mrv[i];
            mktvol_dataProd->MrInput[i] = mrv[i];
        }

        for (i = 0; i < NbVol; i++)
        {
            mktvol_dataProd->Vol[i] = memVol[i];
        }

        /*
        *   CET 
        */
		if (mktvol_dataProd->MrCalibFlag < 4)
		{
            if (Cet_Main (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
			{
                goto RETURN;
			}
		}
		else
		{
            if (Cet_Main_New (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
			{
                goto RETURN;
			}
		}

        
        /* extract vol ratios from cet_out_data */
        fprintf(stream, "\nErrors:");
        exitFlag = TRUE;
        for (i = 0; i < nrDim; i++)
        {
            ratioL[i] = mktvol_dataProd->VolRatio[i];
            errorL[i] = fabs(ratioL[i] - targetRatio[i]);
            exitFlag = exitFlag && (errorL[i] < myErrTol);
            fprintf(stream, "\n%d : %15.10lf", i, ratioL[i] - targetRatio[i]); 
        }

		/* printf("\nErrors:");
		for (i = 0; i < nrDim; i++)
        {
            printf("\n%d : %15.10lf", i, ratioL[i] - targetRatio[i]); 
        } */

        /* shift in mrv[i] direction */
        for (i = 0; i < nrDim; i++)
        {
            for (j = 0; j < nrDim; j++)
            {
                mrv[j] = mrMem[j];
            }

            mrv[i] += delMr;

            for (j = 0; j < nrDim; j++)
            {
                mktvol_dataProd->BetaTD[0][j] = mrv[j];
                mktvol_dataProd->MrInput[j] = mrv[j];
            }

            for (j = 0; j < NbVol; j++)
            {
                mktvol_dataProd->Vol[j] = memVol[j];
            } 

            /*
            *   CET 
            */
            if (mktvol_dataProd->MrCalibFlag < 4)
			{
                if (Cet_Main (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
				{
                    goto RETURN;
				}
			}
		    else
			{
                if (Cet_Main_New (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
				{
                    goto RETURN;
				}
			}

            /* derivatives */
            for (j = 0; j < nrDim; j++)
            {
                JM[j][i] = (mktvol_dataProd->VolRatio[j] - ratioL[j])/delMr;
            }
        } /* for i */

        fprintf(stream, "\nJacobi matrix:");

        for (i = 0; i < nrDim; i++)
        {
            fprintf(stream, "\n");
            for (j = 0; j < nrDim; j++)
            {
                fprintf(stream, "%20.15lf, ", JM[i][j]);
            }
        }

		/* printf("\nJacobi matrix:");

        for (i = 0; i < nrDim; i++)
        {
            printf("\n");
            for (j = 0; j < nrDim; j++)
            {
                printf("%20.15lf, ", JM[i][j]);
            }
        } */

        /* invert Jacobi matrix */
        if (MatrixInverse(nrDim,
                          JM,   
                          JMInv) == FAILURE)
        {
            DR_Error ("Calibration failure: could not invert Jacobi matrix");
            goto RETURN;
        }

        /* NR iteration */
        for (i = 0; i < nrDim; i++)
        {
            tmpSum = 0.;
            for (j = 0; j < nrDim; j++)
            {
                tmpSum += JMInv[i][j] * (ratioL[j] - targetRatio[j]);
            }
            mrv[i] = mrMem[i] - tmpSum;

            if (mrv[i] < MIN_CET_MR)
            {
                if (currIter < mktvol_dataProd->NbMrCet - 1)
                {
                    mrv[i] = MIN_CET_MR;
                }
                else
                {
                    DR_Error ("Calibration failure: negative mean-reversion");
                    goto RETURN;
                }
            }

            if (mrv[i] > MAX_CET_MR)
            {
                if (currIter < mktvol_dataProd->NbMrCet - 1)
                {
                    mrv[i] = MAX_CET_MR;
                }
                else
                {
                    DR_Error ("Calibration failure: mean-reversion too high");
                    goto RETURN;
                }
            }
               
        }

        currIter++;
    }
    
    /* record final mean-reversions */
    for (i = 0; i < NbMr; i++)
    {
        mktvol_dataProd->MrInput[i] = mrv[i];
        mktvol_dataProd->BetaTD[0][i] = mrv[i];
    }

    /* still need to produce vols with final mean-reversions! */
    {
         fprintf(stream, "\nIter %d:", currIter);
         fprintf(stream, "\nFinal mean-reversions: \n");
         for (i = 0; i < nrDim; i++)
         {
             fprintf(stream, "%20.15lf,", mrv[i]);
         }

		 /* printf("\nIter %d:", currIter);
         printf("\nFinal mean-reversions: \n");
         for (i = 0; i < nrDim; i++)
         {
             printf("%20.15lf,", mrv[i]);
         } */

        for (i = 0; i < NbVol; i++)
        {
            mktvol_dataProd->Vol[i] = memVol[i];
        }

        /*
        *   CET 
        */
        if (mktvol_dataProd->MrCalibFlag < 4)
		{
            if (Cet_Main (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
			{
                goto RETURN;
			}
		}
		else
		{
            if (Cet_Main_New (CetOutputFlag,          /* output to CET.prn */
                  t_curveProd,
                  mktvol_dataProd,
                  tree_dataProd) == FAILURE)
			{
                goto RETURN;
			}
		}
        
        /* extract vol ratios from cet_out_data */
        fprintf(stream, "\nErrors:");
        exitFlag = TRUE;
        for (i = 0; i < nrDim; i++)
        {
            ratioL[i] = mktvol_dataProd->VolRatio[i];
            errorL[i] = fabs(ratioL[i] - targetRatio[i]);
            exitFlag = exitFlag && (errorL[i] < myErrTol);
            fprintf(stream, "\n%d : %15.10lf", i, ratioL[i] - targetRatio[i]); 
        }

		/* printf("\nErrors:");
		for (i = 0; i < nrDim; i++)
        {
            printf("\n%d : %15.10lf", i, ratioL[i] - targetRatio[i]); 
        } */
    }

    /* print NR results */
    fprintf(stream, "\nNR results:   %d iterations \n", currIter);
    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "Ratio%d,", i);
    }
    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "Error%d,", i);
    }
    fprintf(stream, "\n");

    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "%20.15lf,", ratioL[i]);
    }

    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "%20.15lf,", ratioL[i] - targetRatio[i]);
    }

    fprintf(stream, "\n");

    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "mr%d,", i);
    }

    fprintf(stream, "\n");

    for (i = 0; i < nrDim; i++)
    {
        fprintf(stream, "%20.15lf,", mrv[i]);
    }           

    fflush(stream); 

	/* print NR results into stdout */
    /* printf("\nNR results:   %d iterations \n", currIter);
    for (i = 0; i < nrDim; i++)
    {
        printf("Ratio%d,", i);
    }
    for (i = 0; i < nrDim; i++)
    {
        printf("Error%d,", i);
    }
    printf("\n");

    for (i = 0; i < nrDim; i++)
    {
        printf("%20.15lf,", ratioL[i]);
    }

    for (i = 0; i < nrDim; i++)
    {
        printf("%20.15lf,", ratioL[i] - targetRatio[i]);
    }

    printf("\n");

    for (i = 0; i < nrDim; i++)
    {
        printf("mr%d,", i);
    }

    printf("\n");

    for (i = 0; i < nrDim; i++)
    {
        printf("%20.15lf,", mrv[i]);
    }           

    fflush(stdout); */ 


    status = SUCCESS;

RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }


    Free_DR_Array(targetRatio, DOUBLE, 0, nrDim - 1);
    Free_DR_Array(ratioL, DOUBLE, 0, nrDim - 1);
    Free_DR_Array(errorL, DOUBLE, 0, nrDim - 1);
    Free_DR_Array(mrv, DOUBLE, 0, nrDim - 1);
    Free_DR_Array(mrMem, DOUBLE, 0, nrDim - 1);
    Free_DR_Matrix(JM, DOUBLE, 0, nrDim - 1, 0, nrDim - 1);
    Free_DR_Matrix(JMInv, DOUBLE, 0, nrDim - 1, 0, nrDim - 1);

    return (status);
}