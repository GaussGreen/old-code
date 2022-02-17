/****************************************************************************/
/*      General I/O routine for callable ladder swap (LADDER)               */
/****************************************************************************/
/*      MANAGER39.c                                                         */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "fix123head.h" 
#include "template39.h"


/*****  Ladder_Manager  ************************************************/
/**
*       Manage input & output, reading from ascii files into data structures.
*/
int     Ladder_Manager 
            (T_CURVE         *t_curve,         /**< (O) Zero curve data      */
             MKTVOL_DATA     *mktvol_data,     /**< (O) Vol  data            */
             LADDER_DATA     *ladder_data,     /**< (O) Deal data            */
             FIX3_TREE_DATA  *tree_data)       /**< (O) Tree data            */
{

    int        i, j;                        /* Convenience indices         */ 
    int        readerror;                   /* Reading error status        */
    int        status = FAILURE;            /* Status = FAILURE initially  */
    char       Index[2][8];                 /* Calibration indices         */
    char       *getserror;                  /* Reading error for fgets     */
    char       OverWriteString[6][MAXBUFF]; /* Overwrite strings           */
    char       IoD;                         /* Local curve number          */
    char       ErrorMsg[MAXBUFF];
    char       FileName[MAXBUFF];
    FILE       *stream = NULL;
    FILE       *streamRisk = NULL;          /* File var for riskzero.dat   */

    strcpy (FileName, "ladder_t.dat");

    /* 
     *   Open the deal data file (see ladder.h) 
     */
    stream = fopen (FileName, "r"); 
                
    if (stream == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open file %s! (Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /*
     * Deal specific data
     */
    
    /* Title */
    if (FindAndSkipComLine (stream,
                            "Title line",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Long or short the option */
    if (FindAndSkipComLine (stream,
                            "Long or short",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->LoS));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read Long/Short in file %s! (Ladder_Manager)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->LoS = (char) toupper (ladder_data->LoS);

    /* Option style */
    if (FindAndSkipComLine (stream, 
                            "Option style",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->Style));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read option style in file %s!(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->Style = (char) toupper (ladder_data->Style);

    /* Number of exercise dates plus dates and strikes */
    if (FindAndSkipComLine (stream, 
                            "number of exercise dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbExer));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of exerc in file %s! Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Exercise dates and strikes",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for (i = 0; i < ladder_data->NbExer; i++)
    {
        readerror = fscanf (stream, "%ld %lf\n",
                          &(ladder_data->ExerDate[i]),
                          &(ladder_data->Strike[i]));
        ladder_data->ExerDate[i] = IRDateFromYMDDate(ladder_data->ExerDate[i]);

        if (readerror != 2)
        {
            sprintf (ErrorMsg,
                     "Could not read exer date/strike #%d in file %s! "
                     "(Ladder_Manager)",
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        ladder_data->Strike[i] /= 100.0;   /* Strike entered in % */
    }  /* for i */

    if (FindAndSkipComLine (stream, "Exercise stats style", "Ladder_Manager", FileName) == FAILURE)
        goto RETURN;

    readerror = fscanf (stream, "%c \n", &(ladder_data->OptStatStyle));

    if (readerror != 1)
    {        
        DR_Error ("Could not read event statistics schedule style %s", FileName);
        goto RETURN;
    }

    ladder_data->OptStatStyle = (char)toupper(ladder_data->OptStatStyle);

    if (FindAndSkipComLine (stream, "Number of event stats dates", "Ladder_Manager", FileName) == FAILURE)
        goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(ladder_data->OptNbStats));

    if (readerror != 1)
    {        
        DR_Error ("Could not read nb of event statistics dates in file %s", FileName);
        goto RETURN;
    }

    if (FindAndSkipComLine(stream, "Event statistics dates", "Ladder_Manager", FileName) == FAILURE)
        goto RETURN;
    
    for (i = 0; i < ladder_data->OptNbStats; i++)
    {
        readerror = fscanf (stream, "%ld \n", &ladder_data->OptStatDates[i]);
        ladder_data->OptStatDates[i] = IRDateFromYMDDate(ladder_data->OptStatDates[i]);
        if (readerror != 1)
        {
            DR_Error ("Could not read event statistics date #%d in file %s",
                     i+1, FileName);
            goto RETURN;
        }
     }

    /* Notional amount */
    if (FindAndSkipComLine (stream,
                            "Notional",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%lf \n", &(ladder_data->Notional));
    if (readerror != 1)
    {
        sprintf (ErrorMsg, 
                 "Could not read notional in file %s! (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->NotionalSign = (ladder_data->Notional > 0) ? 1. : -1.;
    ladder_data->Notional = ABS(ladder_data->Notional);

    /* Original notional percentage */
    if (FindAndSkipComLine (stream,
                            "Original notional percentage",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%lf \n", &(ladder_data->OrgNotPerc));
    if (readerror != 1)
    {
        sprintf (ErrorMsg, 
                 "Could not read original notional percentage in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->OrgNotPerc /= 100.0;

    /* Number of amortization dates plus dates and amortization percentages */
    if (FindAndSkipComLine (stream, 
                            "number of amortization dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbAmort));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of amortizations in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Amortization dates and percentages",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for (i = 0; i < ladder_data->NbAmort; i++)
    {
        readerror = fscanf (stream, "%ld %lf\n",
                          &(ladder_data->AmortDate[i]),
                          &(ladder_data->Amort[i]));
        ladder_data->AmortDate[i] = IRDateFromYMDDate(ladder_data->AmortDate[i]);
        if (readerror != 2)
        {
            sprintf (ErrorMsg,
                     "Could not read amort date/percentage  #%d in file %s! "
                     "(Ladder_Manager)",
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        ladder_data->Amort[i] /= 100.0;   /* Strike entered in % */
    }  /* for i */

    /* Stub convention */
    if (FindAndSkipComLine (stream, 
                            "stub convention", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(ladder_data->StubConv));
    if (readerror != 1)
    {
        sprintf (ErrorMsg, 
                 "Could not read stub convention in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->StubConv = (char) toupper (ladder_data->StubConv);

    /* Swap or zero coupon */
    if (FindAndSkipComLine (stream, 
                            "Swap or zero coupon",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->SoZ));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read swap or zero in file %s! "
                "(Ladder_Manager)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->SoZ = (char) toupper (ladder_data->SoZ);

    /* Smoothing */
    if (FindAndSkipComLine (stream, 
                            "Smoothing",
                            "Ladder_Manager",
                             FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->Smoothing));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read smoothing flag in file %s! "
                "(Ladder_Manager)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->Smoothing = (char) toupper (ladder_data->Smoothing);

    /* -------------------------------------------------------*
     *                  Rib Event                               *
     * -------------------------------------------------------*/

    /* Rib observation event */

    /* Number of rib range dates plus dates, low and high barriers */
    if (FindAndSkipComLine (stream, 
                            "number of rib range dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbRibObsDates));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of observation dates in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Range dates, low and high barriers",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

  

    ladder_data->RibObsDate = 
        (long *) DR_Array (LONG, 0, ladder_data->NbRibObsDates-1);
    ladder_data->RibObsEffDate = 
        (long *) DR_Array (LONG, 0, ladder_data->NbRibObsDates-1);
    ladder_data->RibLoBarrier = 
        (double *) DR_Array (DOUBLE, 0, ladder_data->NbRibObsDates-1);
    ladder_data->RibHiBarrier = 
        (double *) DR_Array (DOUBLE, 0, ladder_data->NbRibObsDates-1);
    ladder_data->RibInRangeWeight = 
        (double *) DR_Array (DOUBLE, 0, ladder_data->NbRibObsDates-1);
    ladder_data->RibOutRangeWeight = 
        (double *) DR_Array (DOUBLE, 0, ladder_data->NbRibObsDates-1);

    if (ladder_data->RibObsDate        == NULL ||
        ladder_data->RibObsEffDate     == NULL ||
        ladder_data->RibLoBarrier      == NULL ||
        ladder_data->RibHiBarrier      == NULL ||
        ladder_data->RibInRangeWeight  == NULL ||
        ladder_data->RibOutRangeWeight == NULL)
    {        
        sprintf (ErrorMsg, 
                 "Could not allocate memory for rib observation events! "
                 "(Ladder_Manager)"); 
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    for (i = 0; i < ladder_data->NbRibObsDates; i++)
    {
        readerror = fscanf (stream, "%ld %ld %lf %lf %lf %lf\n",
                            &(ladder_data->RibObsDate[i]),
                            &(ladder_data->RibObsEffDate[i]),
                            &(ladder_data->RibLoBarrier[i]),
                            &(ladder_data->RibHiBarrier[i]),
                            &(ladder_data->RibInRangeWeight[i]),
                            &(ladder_data->RibOutRangeWeight[i]));
        ladder_data->RibObsDate[i] = IRDateFromYMDDate(ladder_data->RibObsDate[i]);
        ladder_data->RibObsEffDate[i] = IRDateFromYMDDate(ladder_data->RibObsEffDate[i]);
        if (readerror != 6)
        {
            sprintf (ErrorMsg,
                     "Could not read range date/lo/high  #%d in file %s! "
                     "(Ladder_Manager)",
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        ladder_data->RibLoBarrier[i] /= 100.0;   /* Entered in % */
        ladder_data->RibHiBarrier[i] /= 100.0;   /* Entered in % */
    }  /* for i */



    /* Smoothing of range profile */
    if (FindAndSkipComLine (stream,
                            "Smoothing",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->RibSmoothing));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read smoothing flag in file %s! "
                "(Ladder_Manager)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->RibSmoothing = 
                     (char) toupper (ladder_data->RibSmoothing);



    /* Two lines for rib observation indices
     * (mat, freq, dcc, curve, and weights) 
     */
    if (FindAndSkipComLine(stream, 
                           "rib observation index spec", 
                           "Ladder_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for(j = 0; j < 2; j++)
    {
        readerror = fscanf(stream,"%d %c %c %c %lf\n", 
                           &(ladder_data->RibIdxMat[j]),
                           &(ladder_data->RibIdxFreq[j]),
                           &(ladder_data->RibIdxDCC[j]),
                           &IoD,
                           &(ladder_data->RibIdxWeight[j]));
        if (readerror != 5)
        {        
            sprintf (ErrorMsg,
                     "Could not read rib obs index spec in file %s! "
                     "(Ladder_Manager)", 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        IoD = (char)toupper(IoD);
        if ((IoD == 'I') || (IoD == '0'))
        {
            ladder_data->RibIdxIoD[j] = 0;
        }
        else if ((IoD == 'D') || (IoD == '1'))
        {
            ladder_data->RibIdxIoD[j] = 1;
        } 
        else if (IoD == '2')
        {
            ladder_data->RibIdxIoD[j] = 2;
        } 
        else  /* not valid entry, checked later */
        {
            ladder_data->RibIdxIoD[j] = IoD;
        } 

    }

    /* Past Rib observations in range */
    if (FindAndSkipComLine (stream,
                            "Past Rib obs in range",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%lf \n", &(ladder_data->RibPastObsPerc));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read past realized Rib percentage in file %s! "
                "(Ladder_Manager)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->RibPastObsPerc /= 100.0;   /* Entered in % */
                                                     
    /* St leg */

    /* St index, curves and weight. */
    if (FindAndSkipComLine(stream, 
                           "Ladder index spec", 
                           "Ladder_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for(i = 0; i < 2; i++)
    {
        readerror = fscanf(stream,"%d %c %c %c %lf\n",
                       &(ladder_data->IdxMatSt[i]),
                       &(ladder_data->IdxFreqSt[i]),
                       &(ladder_data->IdxBaseSt[i]),
                       &IoD,
                       &(ladder_data->IdxObsWeightSt[i]));
 
        if (readerror != 5)
        {        
           sprintf (ErrorMsg,
                 "Could not read ladder index info in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
           DR_Error (ErrorMsg);
           goto RETURN;
        }

        IoD = (char)toupper(IoD);
        if ((IoD == 'I') || (IoD == '0'))
        {
            ladder_data->IdxIoDSt[i] = 0;
        }
        else if ((IoD == 'D') || (IoD == '1'))
        {
            ladder_data->IdxIoDSt[i] = 1;
        } 
        else if (IoD == '2')
        {
            ladder_data->IdxIoDSt[i] = 2;
        } 
        else  /* not valid entry, checked later */
        {
            ladder_data->IdxIoDSt[i] = IoD;
        } 
    }

    /* Arrears reset (yes or no) */
    if (FindAndSkipComLine (stream, 
                            "Ladder arrears reset", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->ArrearsSt));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read arrears reset of ladder leg in file %s!"
                 " (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->ArrearsSt = (char)toupper(ladder_data->ArrearsSt);
    
    /* Comp convention */
    if (FindAndSkipComLine (stream, 
                            "Ladder comp convention", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->CompSt));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read comp convention of ladder leg in file %s!"
                 " (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->CompSt = (char)toupper(ladder_data->CompSt);
    
    /* Number of fixing dates plus dates and ladder rates IN PERCENT */
    if (FindAndSkipComLine (stream, 
                            "number of fixing dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbFixing));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of fixing in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Fixing dates and percentages",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    
    for (i = 0; i < ladder_data->NbFixing; i++)
    {
        readerror = fscanf (stream, "%ld %lf %lf %lf\n",
                          &(ladder_data->FixingDate[i]),
                          &(ladder_data->FixingSt[0][i]),
                          &(ladder_data->FixingSt[1][i]),
                          &(ladder_data->FixingRibPerc[i]));
        ladder_data->FixingDate[i] = IRDateFromYMDDate(ladder_data->FixingDate[i]);
         if (readerror != 4)
         {
            sprintf (ErrorMsg,
                     "Could not read exer date/percentage  #%d in file %s! "
                     "(Ladder_Manager)",
                     i+1, 
                     FileName);
                DR_Error (ErrorMsg);
                goto RETURN;
        }
            ladder_data->FixingSt[0][i] /= 100.0;   /* Rate entered in % */
            ladder_data->FixingSt[1][i] /= 100.0;   /* Rate entered in % */
            ladder_data->FixingRibPerc[i] /= 100.0; /* Rate entered in % */
    }  /* for i */
    

    /* First lookback value */
    if (FindAndSkipComLine (stream,
                            "First lookback value",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%lf \n", &(ladder_data->FirstLevel));
    if (readerror != 1)
    {
        sprintf (ErrorMsg, 
                 "Could not read first level in file %s! (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->FirstLevel /= 100.0;



    /* Number of ladder dates plus dates infor */
    if (FindAndSkipComLine (stream, 
                            "number of ladder dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbStepUpSt));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of ladders in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Ladders dates and rates",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for (i = 0; i < ladder_data->NbStepUpSt; i++)
    {
        readerror = fscanf (stream, "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                            &(ladder_data->StDate[i]),
                            &(ladder_data->FloorSt[i]),
                            &(ladder_data->CapSt[i]),
                            &(ladder_data->DownRateSt[i]),
                            &(ladder_data->MidRateSt[i]),
                            &(ladder_data->UpRateSt[i]),
                            &(ladder_data->BarrierLo[i]),
                            &(ladder_data->BarrierHi[i]),
                            &(ladder_data->IdxWeightSt[0][i]),
                            &(ladder_data->IdxWeightSt[1][i]),
                            &(ladder_data->SprdSt[i]),
                            &(ladder_data->StickyCoef[i]),
                            &(ladder_data->Leverage[i]),
                            &(ladder_data->IdxFloor[i]),
                            &(ladder_data->IdxCap[i]));
        ladder_data->StDate[i] = IRDateFromYMDDate(ladder_data->StDate[i]);
        if (readerror != 15)
        {
            sprintf (ErrorMsg,
                     "Could not read ladder date/cap/floor/barrier/spread "
                     "#%d in file %s! (Ladder_Manager)",
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        ladder_data->FloorSt[i]    /= 100.0;  /* in percentage */
        ladder_data->CapSt[i]      /= 100.0;  /* in percentage */
        ladder_data->UpRateSt[i]   /= 100.0;  /* in percentage */
        ladder_data->MidRateSt[i]  /= 100.0;  /* in percentage */
        ladder_data->DownRateSt[i] /= 100.0;  /* in percentage */
        ladder_data->BarrierLo[i]  /= 100.0;  /* in percentage */
        ladder_data->BarrierHi[i]  /= 100.0;  /* in percentage */
        ladder_data->SprdSt[i]     /= 100.0;  /* in percentage */
        ladder_data->IdxFloor[i]   /= 100.0;  /* in percentage */
        ladder_data->IdxCap[i]     /= 100.0;  /* in percentage */

    }  /* for i */

   /* St payoff flag */
    if (FindAndSkipComLine (stream,
                            "additive or multiplicative", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->AoM));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read ladder pay frequency in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->AoM = (char)toupper(ladder_data->AoM);

    /* St payment frequency */
    if (FindAndSkipComLine (stream,
                            "payment frequency", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->PayFreqSt));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read ladder pay frequency in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->PayFreqSt = (char)toupper(ladder_data->PayFreqSt);

    /* St payment day count convention */
    if (FindAndSkipComLine (stream,
                            "ladder day count convention", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->DayCountSt));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read ladder day count convention "
                 "in file %s! (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->DayCountSt = (char)toupper(ladder_data->DayCountSt);

    /* Funding leg */

    /* Fling index name and curve and weight */
    if (FindAndSkipComLine (stream, 
                            "Fl index name", 
                            "Ladder_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf(stream,"%d %c %c %c %lf\n",
                   &(ladder_data->IdxMatFl),
                   &(ladder_data->IdxFreqFl),
                   &(ladder_data->IdxBaseFl),
                   &IoD,
                   &(ladder_data->IdxWeightFl));
    if (readerror != 5)
    {        
        sprintf (ErrorMsg,
                 "Could not read floating index name/curve/weight "
                 "in file %s! (Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    IoD = (char)toupper(IoD);
    if ((IoD == 'I') || (IoD == '0'))
    {
        ladder_data->IdxIoDFl = 0;
    }
    else if ((IoD == 'D') || (IoD == '1'))
    {
        ladder_data->IdxIoDFl = 1;
    } 
    else if (IoD == '2')
    {
        ladder_data->IdxIoDFl = 2;
    } 
    else  /* not valid entry, checked later */
    {
        ladder_data->IdxIoDFl = IoD;
    } 

    /* Arrears reset (yes or no) */
    if (FindAndSkipComLine (stream, 
                            "float arrears reset", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->ArrearsFl));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read arrears reset in file %s!"
                 " (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->ArrearsFl = (char)toupper(ladder_data->ArrearsFl);
    
    /* Fling pmt simple or compound convention */
    if (FindAndSkipComLine (stream, 
                            "flt pmt comp conv", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->CompFl));
    if (readerror != 1)
    {
        sprintf (ErrorMsg, 
                 "Could not read float compound conv in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->CompFl = (char) toupper (ladder_data->CompFl);

    /* First fixing of floating leg IN PERCENT */
    if (FindAndSkipComLine (stream, 
                            "Fl first fixing", 
                            "Ladder_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    
    readerror = fscanf (stream, "%lf \n", &(ladder_data->FixingFl));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read float first fixing "
                 "in file %s! (Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (fabs(ladder_data->FixingFl) < TINY)
    {
        ladder_data->FixingGivenFl = FALSE;
    }
    else
    {
        ladder_data->FixingGivenFl = TRUE;
    }
    ladder_data->FixingFl /= 100.0;

    /* Number of stepup dates plus dates infor */
    if (FindAndSkipComLine (stream, 
                            "number of stepup dates", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbStepUpFl));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read nb of step up dates in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream,
                            "Stepup dates and rates",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for (i = 0; i < ladder_data->NbStepUpFl; i++)
    {
        readerror = fscanf (stream, "%ld %lf %lf %lf\n",
                            &(ladder_data->FlDate[i]),
                            &(ladder_data->FloorFl[i]),
                            &(ladder_data->CapFl[i]),
                            &(ladder_data->UpRateFl[i]));
        ladder_data->FlDate[i] = IRDateFromYMDDate(ladder_data->FlDate[i]);
        if (readerror != 4)
        {
            sprintf (ErrorMsg,
                     "Could not read stepup date/rate  #%d in file %s! "
                     "(Ladder_Manager)",
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
        ladder_data->FloorFl[i]  /= 100.0;   /* Rate entered in % */
        ladder_data->CapFl[i]    /= 100.0;   /* Rate entered in % */
        ladder_data->UpRateFl[i] /= 100.0;   /* Rate entered in % */

    }  /* for i */
    
    /* Fling payment frequency */
    if (FindAndSkipComLine (stream,
                            "payment frequency", 
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->PayFreqFl));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read pay frequency in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->PayFreqFl = (char)toupper(ladder_data->PayFreqFl);

    /* Fling payment day count convention */
    if (FindAndSkipComLine (stream,
                            "payment day count convention", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(ladder_data->DayCountFl));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,
                 "Could not read pay day count convention "
                 "in file %s! (Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    ladder_data->DayCountFl = (char)toupper(ladder_data->DayCountFl);




    /*
     * Model specific data
     */

    /* First index to be used for volatility calibration */
    /* The index name is put directly into tree_data */
    if (FindAndSkipComLine (stream,
                            "first volatility index",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%s \n", Index[0]);
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read first volatility index in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /* Second index to be used for volatility calibration */
    if (FindAndSkipComLine (stream,
                            "second volatility index",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%s \n", Index[1]);
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read second volatility index in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Zero curve to use as discount curve */
    if (FindAndSkipComLine (stream,
                            "curve for discounting",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(tree_data->CvDisc));
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read index of curve for discounting in file %s! "
                 "(Cmsopt_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Number of standard deviations at which to cut the tree */
    /* NbSigmaMax is also put directly into tree_data         */
    if (FindAndSkipComLine (stream,
                            "number of standard deviations",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(tree_data->NbSigmaMax));
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read number of standard deviations in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    
    /* Number of state variables */
    if (FindAndSkipComLine (stream,
                            "number of state variables",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbStates));
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read number of state variables in file %s!"
                 " (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
 
    /* Number of standard deviations for state variable */
    if (FindAndSkipComLine (stream,
                            "number of standard deviations for state",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(ladder_data->NbStDev));
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read number of standard deviations for state "
                 "in file %s! (Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
  
    /* Number of periods per year in the tree */
    if (FindAndSkipComLine (stream,
                            "periods per year",
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[0], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not read periods per year in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    
    /* Model */
    /* This is read as an overwrite string    */
    if (FindAndSkipComLine (stream,
                            "model type",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[1], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, 
                 "Could not read model type in file %s! "
                 "(Amortfloat_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Number of factors */
    if (FindAndSkipComLine (stream,
                            "number of factors",
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(tree_data->NbFactor));
    if (readerror != 1)
    {
        sprintf (ErrorMsg,
                 "Could not read number of factors in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Factor weights */
    if (FindAndSkipComLine (stream,
                            "factor weights",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[2], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, 
                 "Could not read factor weight in file %s! "
                 "(Ladder_Manager)",
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Mean reversion coefficients */
    if (FindAndSkipComLine (stream,
                            "mean reversion", 
                            "Ladder_Manager", 
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[3], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, 
                 "Could not read mean reversion in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Correlations */
    if (FindAndSkipComLine (stream, 
                            "correlation",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[4], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, 
                 "Could not read correlation in file %s! "
                 "(Ladder_Manager)", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Backbone */
    if (FindAndSkipComLine (stream,
                            "bone",
                            "Ladder_Manager",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }
    getserror = fgets (OverWriteString[5], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, "Could not read bone in file %s! (Ladder_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /* ------------------------------------------------------------------*/
    /*   STANDARD ENVIRONMENT ELEMENTS (i.e. zero.dat, basevol.dat etc.) */
    /*-------------------------------------------------------------------*/

    /* Hardcoded assigment of zero curves for the engine */
    tree_data->CvDiff = 0;  /* zero.dat is index curve */
    tree_data->CvIdx1 = 1;
    tree_data->CvIdx2 = 2;
    /* discount curve is defined in the deal file */

    /* Set smoothing flag for CET */
    mktvol_data->SmoothingFlag = ladder_data->Smoothing;

    if (Fix3_MktAndModel_Input (mktvol_data,
                                tree_data,
                                t_curve,
                                Index,
                                OverWriteString) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;
        
    RETURN:

    if (status == FAILURE)
        Ladder_Data_Free (ladder_data);

    if (stream     != NULL) fclose (stream);
    if (streamRisk != NULL) fclose (streamRisk);
    
    return (status);

}  /* Ladder_Manager */



/*****  Ladder_PreProcess  ************************************************/
/**
*       Manage input & output, reading from ascii files into data structures.
*/
int     Ladder_PreProcess
            (T_CURVE         *t_curve,         /**< (I) Zero curve data      */
             MKTVOL_DATA     *mktvol_data,     /**< (I/O) Vol  data          */
             long            ValueDate,        /**< (I) Value date           */
             LADDER_DATA     *ladder_data,     /**< (I/O) Deal data          */
             FIX3_TREE_DATA  *tree_data)       /**< (I) Tree data            */
{

    int     i, j;
    int     status = FAILURE;               /* Status = FAILURE initially  */
    double  norm;

    /* Set accrual start date */
    ladder_data->AccStDate = ladder_data->StDate[0];

   /* Set final maturity date  */
    ladder_data->MatDate = ladder_data->AmortDate[ladder_data->NbAmort-1];

    /* Statistics calculation */
    ladder_data->CalcStats = (ladder_data->OptNbStats > 0 ? 'Y': 'N');

    /* Floating index on or not */
    ladder_data->IdxOnFl = (ABS(ladder_data->IdxWeightFl) > TINY ? TRUE : FALSE);

    /* ... number of Rib observation dates, */
    ladder_data->CplxIsRib = (ladder_data->NbRibObsDates > 0 ? TRUE : FALSE);
    ladder_data->MaxNbRib  = 0; /* will be determined later */
    ladder_data->NbPastRibObsDates = 0;

    /* ... Ladder index weights, */
    for (j = 0; j < 2; j++)
    {
        ladder_data->IdxOnSt[j] = 
            ((ABS (ladder_data->IdxObsWeightSt[j]) < TINY) ? FALSE : TRUE);
 
        /* set to true nevertheless if the PAYMENT weights are non-zero */
        for (i = 0; i < ladder_data->NbStepUpSt; i++)
        {
            if (ABS (ladder_data->IdxWeightSt[j][i] ) > TINY)
            {
                ladder_data->IdxOnSt[j] = TRUE;
            }
        }
    }

    /* ... Rib index weights, */
    for (j = 0; j < 2; j++)
    {
        ladder_data->RibIdxOn[j] = (ABS(ladder_data->RibIdxWeight[j])<TINY ? FALSE:TRUE);
    }

    if (ladder_data->AccStDate > ValueDate || !(ladder_data->CplxIsRib))
    {
        ladder_data->RibPastObsWeight         = 0.;
        ladder_data->RibPastObsInRangeWeight  = 1.;
        ladder_data->RibPastObsOutRangeWeight = 0.;
    }
    else
    {
        EVENT_LIST  *StEventList = NULL;     /* Event list for fixed pmts      */

        int CplxIdx = 1, RibIdx = 0, NbPastObs = 0, NbCurrObs = 0;
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
    
        if (StEventList == NULL) goto RETURN;

        /* determine ratio of past Rib obs / total Rib obs in 
         * complex period including ValueDate */

        while ((CplxIdx < StEventList->NbEntries) &&
               (ValueDate > StEventList->Dates[CplxIdx]))
        {
            CplxIdx ++;
        }

        while ((RibIdx < ladder_data->NbRibObsDates) &&
               (ladder_data->RibObsEffDate[RibIdx] < StEventList->Dates[CplxIdx-1]))
        {
            RibIdx ++;
        }

        ladder_data->RibPastObsInRangeWeight  = ladder_data->RibInRangeWeight[RibIdx];
        ladder_data->RibPastObsOutRangeWeight = ladder_data->RibOutRangeWeight[RibIdx];

        while ((RibIdx < ladder_data->NbRibObsDates) &&
               (ladder_data->RibObsEffDate[RibIdx] < StEventList->Dates[CplxIdx]))
        {
            if (ladder_data->RibObsEffDate[RibIdx] < ValueDate) 
            {
                NbPastObs ++;
            }

            NbCurrObs ++;
            RibIdx ++;
        }

        ladder_data->RibPastObsWeight = ((double) NbPastObs) / NbCurrObs;

        DrFreeEventList(StEventList);
    }
 
    /* Set indices for index curves according to diff curve */
    
    if (tree_data->CvDiff == 0)
    {
        tree_data->CvIdx1 = 1;
        tree_data->CvIdx2 = 2;
    }
    else if (tree_data->CvDiff == 1)
    {
        tree_data->CvIdx1 = 0;
        tree_data->CvIdx2 = 2;
    }
    else /* (tree_data->CvDiff == 2) */
    {
        tree_data->CvIdx1 = 0;
        tree_data->CvIdx2 = 1;
    }

    /* Set total vol constants in mktvol_data */
    norm = 0.;

    for (i = 0; i < tree_data->NbFactor; i++) 
        norm += mktvol_data->Alpha[i] * mktvol_data->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha (%lf) is too small !", norm);
        goto RETURN;
    }       
    if (IS_EQUAL(mktvol_data->Bbq,1))
    {
        mktvol_data->VolNorm = 0.;
        mktvol_data->VolLogn = norm;
    }
    else
    if (IS_EQUAL(mktvol_data->Bbq,0))
    {
        mktvol_data->VolNorm = norm;
        mktvol_data->VolLogn = 0.;
    }
    else
    {
        DR_Error("Bbq parameter must be either 0 or 1 !");
        goto RETURN;
    }       

    /* Check input deal data */
    if (Ladder_Check (ladder_data,
                        t_curve[0].ValueDate,
                        tree_data) == FAILURE)
    {
        goto RETURN;
    }

    /* No modifications of tree_data */
 
    /* Check input market data */
    for (i = 0; i < 3; i++)
        if (Term_Check_W (&(t_curve[i])) == FAILURE)
        {        
            goto RETURN;
        }

    if (MktVol_Check_W (mktvol_data)  == FAILURE)
    {
        goto RETURN;
    }
    
    if (Fix3_Param_Check (tree_data->NbFactor,
                     mktvol_data,
                     tree_data) == FAILURE)
    {        
        goto RETURN;
    }


    status = SUCCESS;
        
RETURN:

    return (status);

} /* Ladder_PreProcess */



/*****  Ladder_Check  **************************************************/
/**
*       Check the inputs.
*/
int     Ladder_Check 
           (LADDER_DATA      *ladder_data,    /**< (I) Ladder data            */
            long             ValueDate,       /**< (I) Value date             */
            FIX3_TREE_DATA   *tree_data)      /**< (I) Structure of tree data */
{

    int         i, j;
    int         status = FAILURE;       /* Error status = FAILURE initially */
    double      TotalAmort;
    int         NbDates = 0;            /* Local variables for exercise list*/  
    int         DateListStart;
    long       *DateList = NULL;


    /* OPTION */

    /* Long or short */
    if ((ladder_data->LoS != 'L') && (ladder_data->LoS != 'S') ) 
    {
        DR_Error ("Option must be long or short (L or S)!");
        goto RETURN;
    } 

    /* European, American or frequency of exercise */
    if ((ladder_data->Style != 'E') && 
        (ladder_data->Style != 'M') && 
        (ladder_data->Style != 'Q') &&  
        (ladder_data->Style != 'S') &&
        (ladder_data->Style != 'A') )   
    {
        DR_Error ("Invalid entry for exercise type!");
        goto RETURN;
    } 

    /* Nb of exercise dates can only be 1 for European options */
    if (ladder_data->NbExer < 1)
    {
        DR_Error ("At least one exercise date is required!");
        goto RETURN;
    }                                               
    
    if ((ladder_data->NbExer == 1) && (ladder_data->Style != 'E'))
    {
        DR_Error ("A minimum of two dates is required, "
                  "except for European options!");
        goto RETURN;
    }   

    /* Exercise date format and ascending order */
    for (i = 0; i < ladder_data->NbExer; i++)
        if (Dateok(ladder_data->ExerDate[i]))
        {
            DR_Error ("Incorrect format for exercise date!");
            goto RETURN;
        }

    for (i = 1; i < ladder_data->NbExer; i++)
        if (ladder_data->ExerDate[i] <= ladder_data->ExerDate[i-1])
        {
            DR_Error ("Exercise dates must be entered in ascending order!");
            goto RETURN;
        }
    
    /* Unless exercise on given dates only, dates must match frequency */
    if (ladder_data->Style != 'E')
    {
        if (DrDatesInSchedule(ladder_data->NbExer,
                              ladder_data->ExerDate,
                              ladder_data->ExerDate[0],
                              ladder_data->ExerDate[ladder_data->NbExer-1],
                              ladder_data->Style,
                              'N') == FAILURE)
        {
            DR_Error("Specified exercise dates are not subset of dates\n"
                     "obtained from exercise frequency.\n");
            goto RETURN;
        }
    }
    
    /* Check that exercise dates are on or before final maturity of swap */
    if (ladder_data->ExerDate[ladder_data->NbExer-1] 
        > ladder_data->MatDate)
    {
        DR_Error("Exercise window must not extend further than swap maturity!");
        goto RETURN;
    }
       
    /* Exercise dates after accrual start only possible on pmt dates */
    /* To check this, create exercise date list starting after acc st*/
    if (ladder_data->Style != 'E') 
    {
        if (DateListFromFreq(ladder_data->ExerDate[0],
                             ladder_data->ExerDate[ladder_data->NbExer-1],
                             ladder_data->Style,
                             'N',             /* no stub in exercise list */
                             &NbDates,    
                             &DateList) == FAILURE)
        {
            DR_Error("Error creating local exercise list. \n");
            goto RETURN;
        }
    }
    else
    {
        NbDates = ladder_data->NbExer;
        DateList = (long *) DR_Array (LONG, 0, NbDates-1);
        if (DateList == NULL)
        {
            DR_Error("Error creating local date list. \n");
            goto RETURN;
        }
        for (i = 0; i < NbDates; i++)
        {
            DateList[i] = ladder_data->ExerDate[i];
        }
    }
    DateListStart = 0;
    for (i = 0; i < NbDates; i++)
    {
        if (DateList[i] < ladder_data->AccStDate) DateListStart = i+1;
    }
    
    /* Exercise dates after accrual start only possible on pmt dates */
    if (DateListStart < NbDates )
    {
        if (DrDatesInSchedule(NbDates - DateListStart,
                              &(DateList[DateListStart]),
                              ladder_data->AccStDate,
                              ladder_data->MatDate,
                              ladder_data->PayFreqSt,
                              ladder_data->StubConv) == FAILURE)
        {
            DR_Error("Exercise dates are not a subset of ladder pmt dates !");
            goto RETURN;
        }
        if (DrDatesInSchedule(NbDates - DateListStart,
                              &(DateList[DateListStart]),
                              ladder_data->AccStDate,
                              ladder_data->MatDate,
                              ladder_data->PayFreqFl,
                              ladder_data->StubConv) == FAILURE)
        {
            DR_Error("Exercise dates are not a subset of floating pmt dates !");
            goto RETURN;
        }
    }

    /* Exercise stats style */
    if ((ladder_data->OptStatStyle != 'I') && 
        (ladder_data->OptStatStyle != 'D') && 
        (ladder_data->OptStatStyle != 'W') && 
        (ladder_data->OptStatStyle != 'M') && 
        (ladder_data->OptStatStyle != 'Q') &&  
        (ladder_data->OptStatStyle != 'S') &&
        (ladder_data->OptStatStyle != 'A') )   
    {
        DR_Error ("Invalid entry for event schedule type");
        goto RETURN;
    } 

    if ((ladder_data->OptNbStats == 1) && (ladder_data->OptStatStyle != 'I'))            
    {
        DR_Error ("A minimum of two dates is required, "
                 "except for 'I' schedule type");
        goto RETURN;
    }   

    for (i = 0; i < ladder_data->OptNbStats; i++)
    {
        if (Dateok(ladder_data->OptStatDates[i]))
        {
            DR_Error ("Incorrect format for event statistics date");
            goto RETURN;
        }  
        if (i > 0 && ladder_data->OptStatDates[i] < ladder_data->OptStatDates[i-1])
        {
            DR_Error("Event dates must be given in ascending order\n");
            goto RETURN;
        }
    }

    /* SWAP */

    /* Nb of amortization dates can must be at least 1. */
    if (ladder_data->NbAmort < 1)
    {
        DR_Error ("At least one amortization date is required!");
        goto RETURN;
    }                                               
    
    /* Amortization date format and ascending order */
    for (i = 0; i < ladder_data->NbAmort; i++)
        if (Dateok(ladder_data->AmortDate[i]))
        {
            DR_Error ("Incorrect format for amortization date!");
            goto RETURN;
        }

    for (i = 1; i < ladder_data->NbAmort; i++)
        if (ladder_data->AmortDate[i] <= ladder_data->AmortDate[i-1])
        {
            DR_Error ("Amortization dates must be entered in ascending order!");
            goto RETURN;
        }
    
    /* Amortization dates only possible on fix & flt coupon dates */
    if (DrDatesInSchedule(ladder_data->NbAmort,
                          ladder_data->AmortDate,
                          ladder_data->AccStDate,
                          ladder_data->MatDate,
                          ladder_data->PayFreqSt,
                          ladder_data->StubConv) == FAILURE)
    {
        DR_Error("Amortisation dates are not subset of ladder pmt dates.\n");
        goto RETURN;
    }
    if (DrDatesInSchedule(ladder_data->NbAmort,
                          ladder_data->AmortDate,
                          ladder_data->AccStDate,
                          ladder_data->MatDate,
                          ladder_data->PayFreqFl,
                          ladder_data->StubConv) == FAILURE)
    {
        DR_Error("Amortisation dates are not subset of flt pmt dates.\n");
        goto RETURN;
    }
    if (ladder_data->AmortDate[0] == ladder_data->AccStDate)
    {
        DR_Error("Amortization on accrual start is not supported. \n");
        goto RETURN;
    }
 
    /* Sum of original percentage and all amortizations must be zero */
    TotalAmort = ladder_data->OrgNotPerc;
    for (i = 0; i < ladder_data->NbAmort; i++)
    {
        TotalAmort -= ladder_data->Amort[i];
        if (TotalAmort < -ERROR)
        {
            DR_Error ("Agregate amortization is not positive !");
                    goto RETURN;
        }
    }
    if (fabs (TotalAmort) > ERROR)
    {
        DR_Error ("Original percentage and amortization rates don't add up "
                  "to 0%%! total amortization = %f (percent)\n", 
                  100 * TotalAmort);
        goto RETURN;
    }

    /* No amortization schedule for zero coupon swap */
    if ((ladder_data->NbAmort > 1) && (ladder_data->SoZ == 'Z'))
    {
        DR_Error("Amort schedule for zero coupon swap is not supported");
        goto RETURN;
    }

    /* For Stub location */
    if ((ladder_data->StubConv != 'F') && 
        (ladder_data->StubConv != 'B') )   
    {
        DR_Error("Invalid entry for stub location ('F' or 'B')!");
        goto RETURN;
    } 

    /* Accrual start from both legs */
    if (ladder_data->StDate[0] != ladder_data->FlDate[0])
    {
        DR_Error("Accrual start on ladder and floating leg is not the same !");
        goto RETURN;
    } 

    /* Swap or zero coupon */
    if ((ladder_data->SoZ != 'S') && (ladder_data->SoZ != 'Z') ) 
    {
        DR_Error ("Underlying must be swap or zero copoun (S or Z)!");
        goto RETURN;
    } 

    /* Smoothing */
    if ((ladder_data->Smoothing != 'Y') && (ladder_data->Smoothing != 'N') ) 
    {
        DR_Error ("Smoothing flag must be either Y or N!");
        goto RETURN;
    } 

    /* For consistency, the ladder smoothing flag should be the same
     * as rib smoothing 
     */
    if (ladder_data->RibSmoothing != ladder_data->Smoothing)
    {
        DR_Error("Option and rib smoothing flags must be the same.\n");
        goto RETURN;
    }

    /* Rib Observation events */

    if (ladder_data->CplxIsRib)
    {
        /* Range date format and ascending order, lo < high barrier */
        for (i = 0; i < ladder_data->NbRibObsDates; i++)
        {
            if (Dateok(ladder_data->RibObsDate[i]))
            {
                DR_Error ("Incorrect format for rib observation date!");
                goto RETURN;
            } 

            if (Dateok(ladder_data->RibObsEffDate[i]))
            {
                DR_Error ("Incorrect format for rib observation eff date!");
                goto RETURN;
            } 
        }

        for (i = 1; i < ladder_data->NbRibObsDates; i++)
        {
            if (ladder_data->RibObsDate[i] <= 
                ladder_data->RibObsDate[i-1])
            {
                DR_Error ("Rib Obs dates must be entered in ascending order!");
                goto RETURN;
            } 

            if (ladder_data->RibObsEffDate[i] <= 
                ladder_data->RibObsEffDate[i-1])
            {
                DR_Error ("Rib Obs eff dates must be entered in ascending order!");
                goto RETURN;
            } 
        }
    
        for (i = 0; i < ladder_data->NbRibObsDates; i++)
        {
            if ((ladder_data->RibHiBarrier[i] + TINY) <= 
                 ladder_data->RibLoBarrier[i])
            {
                DR_Error ("Low barrier must not be bigger than high barrier !");
                goto RETURN;
            } 

            /* Both rib event weights for indices cannot be equal to zero */
            if ((ABS(ladder_data->RibInRangeWeight[i])< TINY) &&
                (ABS(ladder_data->RibOutRangeWeight[i])< TINY))
            {
                DR_Error("Both event weights cannot be equal to zero !");
                goto RETURN;
            }
        }    




        /*  Rib observation indices */
        for (j=0; j<2; j++)
        {
            /* maturity */
            if(    (ladder_data->RibIdxMat[j] != 1)
                && (ladder_data->RibIdxMat[j] != 3)
                && (ladder_data->RibIdxMat[j] != 6)
                && (ladder_data->RibIdxMat[j] != 12)
                && (ladder_data->RibIdxMat[j] != 24)
                && (ladder_data->RibIdxMat[j] != 36)
                && (ladder_data->RibIdxMat[j] != 48)
                && (ladder_data->RibIdxMat[j] != 60)
                && (ladder_data->RibIdxMat[j] != 72)
                && (ladder_data->RibIdxMat[j] != 84)
                && (ladder_data->RibIdxMat[j] != 96)
                && (ladder_data->RibIdxMat[j] != 108)
                && (ladder_data->RibIdxMat[j] != 120)
                && (ladder_data->RibIdxMat[j] != 144)
                && (ladder_data->RibIdxMat[j] != 180)
                && (ladder_data->RibIdxMat[j] != 240)
                && (ladder_data->RibIdxMat[j] != 360))
            {
                DR_Error("Invalid rib obs index (%d) maturity (%d)!\n",
                                 j+1,
                                 ladder_data->RibIdxMat[j]);
                goto RETURN;
            }

            /* frequency */
            if ((ladder_data->RibIdxFreq[j] != 'A') && 
                (ladder_data->RibIdxFreq[j] != 'S') &&
                (ladder_data->RibIdxFreq[j] != 'Q') &&
                (ladder_data->RibIdxFreq[j] != 'M'))
            {
                DR_Error("Rib obs index frequency must be A, S, Q or M!\n");
                goto RETURN;
            }

            /* frequency vs maturity */
            if (   ((ladder_data->RibIdxMat[j] == 1) && 
                    (ladder_data->RibIdxFreq[j] != 'M'))
                || ((ladder_data->RibIdxMat[j] == 3) && 
                    (ladder_data->RibIdxFreq[j] != 'Q'))
                || ((ladder_data->RibIdxMat[j] == 6) && 
                    (ladder_data->RibIdxFreq[j] != 'S')))
            {
                DR_Error("Rbi obs index frequency does not agree with "
                         "its maturity!\n");
                goto RETURN;
            }

            /* day count */
            if ((ladder_data->RibIdxDCC[j] != '0') && 
                (ladder_data->RibIdxDCC[j] != '5') && 
                (ladder_data->RibIdxDCC[j] != '3') )
            {
                DR_Error ("Rib obs index day count conv must be 0, 3 or 5!\n");
                goto RETURN;
            } 

            /* curve */
            if ((ladder_data->RibIdxIoD[j] != 0) && 
                (ladder_data->RibIdxIoD[j] != 1) &&
                (ladder_data->RibIdxIoD[j] != 2))
            {
                DR_Error("Curve for index in rib obs must be 0, 1, 2.\n");
                goto RETURN;
            }

    } /* j = 0,1 */


    /* Not allowed to have both indices set to zero weights */
    if ((ABS(ladder_data->RibIdxWeight[0])< TINY) &&
        (ABS(ladder_data->RibIdxWeight[1])< TINY)  )
    {
        DR_Error("Weights for rib obs indices cannot be both equal zero !");
        goto RETURN;
    }
    } /* If RIB */

    /* STICKY LEG */


    /* Ladder Indices */
    for (j=0; j<2; j++)
    {
        /* maturity */
        if(    (ladder_data->IdxMatSt[j] != 1)
            && (ladder_data->IdxMatSt[j] != 3)
            && (ladder_data->IdxMatSt[j] != 6)
            && (ladder_data->IdxMatSt[j] != 12)
            && (ladder_data->IdxMatSt[j] != 24)
            && (ladder_data->IdxMatSt[j] != 36)
            && (ladder_data->IdxMatSt[j] != 48)
            && (ladder_data->IdxMatSt[j] != 60)
            && (ladder_data->IdxMatSt[j] != 72)
            && (ladder_data->IdxMatSt[j] != 84)
            && (ladder_data->IdxMatSt[j] != 96)
            && (ladder_data->IdxMatSt[j] != 108)
            && (ladder_data->IdxMatSt[j] != 120)
            && (ladder_data->IdxMatSt[j] != 144)
            && (ladder_data->IdxMatSt[j] != 180)
            && (ladder_data->IdxMatSt[j] != 240)
            && (ladder_data->IdxMatSt[j] != 360))
        {
            DR_Error("Ladder rib obs index (%d) maturity (%d)!\n",
                             j+1,
                             ladder_data->IdxMatSt[j]);
            goto RETURN;
        }

        /* frequency */
        if ((ladder_data->IdxFreqSt[j] != 'A') && 
            (ladder_data->IdxFreqSt[j] != 'S') &&
            (ladder_data->IdxFreqSt[j] != 'Q') &&
            (ladder_data->IdxFreqSt[j] != 'M'))
        {
            DR_Error("Ladder index frequency must be A, S, Q or M!\n");
            goto RETURN;
        }

        /* frequency vs maturity */
        if (   ((ladder_data->IdxMatSt[j] == 1) && 
                (ladder_data->IdxFreqSt[j] != 'M'))
            || ((ladder_data->IdxMatSt[j] == 3) && 
                (ladder_data->IdxFreqSt[j] != 'Q'))
            || ((ladder_data->IdxMatSt[j] == 6) && 
                (ladder_data->IdxFreqSt[j] != 'S')))
        {
            DR_Error("Ladder index frequency does not agree with "
                     "its maturity!\n");
            goto RETURN;
        }

        /* day count */
        if ((ladder_data->IdxBaseSt[j] != '0') && 
            (ladder_data->IdxBaseSt[j] != '5') && 
            (ladder_data->IdxBaseSt[j] != '3') )
        {
            DR_Error ("Ladder index day count conv must be 0, 3 or 5!\n");
            goto RETURN;
        } 
    } 

    /* Curves for indices must be 0,1,2 */
    for(i = 0; i <2; i++)
    {
       if ( (ladder_data->IdxIoDSt[i] != 0) && 
            (ladder_data->IdxIoDSt[i] != 1) &&
            (ladder_data->IdxIoDSt[i] != 2) )
       {
           DR_Error("Curve for ladder index estimation must be 0,1,2.\n");
           goto RETURN;
       }
    }

    /* For Arrears reset */
    if ((ladder_data->ArrearsSt != 'Y') && 
        (ladder_data->ArrearsSt != 'N') )   
    {
        DR_Error("Invalid entry for ladder arrears reset ('Y' or 'N')!");
        goto RETURN;
    } 

    /* For comp ladder conv */
    if ((ladder_data->CompSt != 'S') && 
        (ladder_data->CompSt != 'C') )   
    {
        DR_Error("Invalid entry for comp ladder ('S' or 'C')!");
        goto RETURN;
    } 

    /* Check ladder past fixing dates */
    for (i = 0; i < ladder_data->NbFixing; i++)
        if (Dateok(ladder_data->FixingDate[i]))
        {
            DR_Error ("Incorrect format for past fixing date!");
            goto RETURN;
        }

    for (i = 1; i < ladder_data->NbFixing; i++)
        if (ladder_data->FixingDate[i] <= ladder_data->FixingDate[i-1])
        {
            DR_Error ("Past fixing dates must be entered in ascending order!");
            goto RETURN;
        }

    /* Check agreement between given and created fixing date list */
    if ( DrSameDateSchedules
              (&(ladder_data->NbFixing),
               ladder_data->FixingDate,
               ladder_data->AccStDate,
               ladder_data->MatDate,
               ValueDate, 
               ladder_data->PayFreqSt,
               ladder_data->StubConv,
               ladder_data->ArrearsSt) == FAILURE )
    {
        DR_Error("Past fixings do not agree with ladder schedule !");
        goto RETURN;
    }

    /* Nb of stepup dates can must be at least 1. */
    if (ladder_data->NbStepUpSt < 1)
    {
        DR_Error ("At least one ladder date is required!");
        goto RETURN;
    }                                               
    
    /* Stepup date format, ascending order, rate order */
    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (Dateok(ladder_data->StDate[i]))
        {
            DR_Error ("Incorrect format for ladder date!");
            goto RETURN;
        }

    for (i = 1; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->StDate[i] <= ladder_data->StDate[i-1])
        {
            DR_Error ("Ladder dates must be entered in ascending order!");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->CapSt[i] < ladder_data->FloorSt[i])
        {
            DR_Error ("Cap rate must be bigger than floor rate !");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->IdxCap[i] < ladder_data->IdxFloor[i])
        {
            DR_Error ("Index Cap rate must be bigger than Floor rate !");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->BarrierHi[i] + TINY < ladder_data->BarrierLo[i])
        {
            DR_Error ("High barrier can not be below low barrier !");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->StickyCoef[i] < - TINY)
        {
            DR_Error ("Sticky coefficients need to be positive !");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpSt; i++)
        if (ladder_data->Leverage[i] < - TINY)
        {
            DR_Error ("Leverage needs to be non-negative!\n");
            goto RETURN;
        }

    
    /* Stepup dates only possible on ladder leg dates */
    if (DrDatesInSchedule(ladder_data->NbStepUpSt,
                          ladder_data->StDate,
                          ladder_data->AccStDate,
                          ladder_data->MatDate,
                          ladder_data->PayFreqSt,
                          ladder_data->StubConv) == FAILURE)
    {
        DR_Error("Ladder dates are not a subset of ladder pmt dates.\n");
        goto RETURN;
    }
    if (ladder_data->StDate[ladder_data->NbStepUpSt-1] 
        == ladder_data->MatDate)
    {
        DR_Error("Deal cannot ladder on maturity date.\n");
        goto RETURN;
    }


    /* Additive or multiplicative step up */
    if( (ladder_data->AoM != 'A') &&
        (ladder_data->AoM != 'B') &&
        (ladder_data->AoM != 'M'))
    {
        DR_Error ("Invalid entry for swap ladder step up type!");
        goto RETURN;
    }


    /* For ladder coupon payments, the usual frequencies */
    if ((ladder_data->PayFreqSt != 'M') && 
        (ladder_data->PayFreqSt != 'Q') && 
        (ladder_data->PayFreqSt != 'S') &&
        (ladder_data->PayFreqSt != 'A') )   
    {
        DR_Error ("Invalid entry for swap ladder payment frequency!");
        goto RETURN;
    } 
                                                  
    /* For fix coupon day count convention */
    if ((ladder_data->DayCountSt != '0') && 
        (ladder_data->DayCountSt != '5') && 
        (ladder_data->DayCountSt != '3') )   
    {
        DR_Error ("Invalid entry for swap ladder day count convention!");
        goto RETURN;
    } 



    /* FLOATING LEG */
    
    /* Floating Index */

    /* maturity */
    if(    (ladder_data->IdxMatFl != 1)
        && (ladder_data->IdxMatFl != 3)
        && (ladder_data->IdxMatFl != 6)
        && (ladder_data->IdxMatFl != 12)
        && (ladder_data->IdxMatFl != 24)
        && (ladder_data->IdxMatFl != 36)
        && (ladder_data->IdxMatFl != 48)
        && (ladder_data->IdxMatFl != 60)
        && (ladder_data->IdxMatFl != 72)
        && (ladder_data->IdxMatFl != 84)
        && (ladder_data->IdxMatFl != 96)
        && (ladder_data->IdxMatFl != 108)
        && (ladder_data->IdxMatFl != 120)
        && (ladder_data->IdxMatFl != 144)
        && (ladder_data->IdxMatFl != 180)
        && (ladder_data->IdxMatFl != 240)
        && (ladder_data->IdxMatFl != 360))
    {
        DR_Error("Invalid floating index (%d) maturity (%d)!\n",
                         j+1,
                         ladder_data->IdxMatFl);
        goto RETURN;
    }

    /* frequency */
    if ((ladder_data->IdxFreqFl != 'A') && 
        (ladder_data->IdxFreqFl != 'S') &&
        (ladder_data->IdxFreqFl != 'Q') &&
        (ladder_data->IdxFreqFl != 'M'))
    {
        DR_Error("Floating index frequency must be A, S, Q or M!\n");
        goto RETURN;
    }

    /* frequency vs maturity */
    if (   ((ladder_data->IdxMatFl == 1) && 
            (ladder_data->IdxFreqFl != 'M'))
        || ((ladder_data->IdxMatFl == 3) && 
            (ladder_data->IdxFreqFl != 'Q'))
        || ((ladder_data->IdxMatFl == 6) && 
            (ladder_data->IdxFreqFl != 'S')))
    {
        DR_Error("Floating index frequency does not agree with "
                 "its maturity!\n");
        goto RETURN;
    }

    /* day count */
    if ((ladder_data->IdxBaseFl != '0') && 
        (ladder_data->IdxBaseFl != '5') && 
        (ladder_data->IdxBaseFl != '3') )
    {
        DR_Error ("Floating index day count conv must be 0, 3 or 5!\n");
        goto RETURN;
    } 

    if ((ladder_data->IdxIoDFl != 0) && 
        (ladder_data->IdxIoDFl != 1) &&
        (ladder_data->IdxIoDFl != 2))
    {
        DR_Error("Curve for float index estimation must be 0,1,2.\n");
        goto RETURN;
    }

    /* For Arrears reset */
    if ((ladder_data->ArrearsFl != 'Y') && 
        (ladder_data->ArrearsFl != 'N') )   
    {
        DR_Error("Invalid entry for arrears reset ('Y' or 'N')!");
        goto RETURN;
    } 

    /* For comp floating conv */
    if ((ladder_data->CompFl != 'S') && 
        (ladder_data->CompFl != 'C') )   
    {
        DR_Error("Invalid entry for comp float ('S' or 'C')!");
        goto RETURN;
    } 

    /* First fixing of the flt leg can only be given if not arrears */
    if ( (ladder_data->FixingGivenFl) &&
         (ladder_data->ArrearsFl == 'Y') )
    {
        DR_Error("Fl first fixing can only be specified if float reset "
                "is in advance!");
        goto RETURN;
    }

    /* Not allowed to have floating zero weight and advance reset */
    if (!(ladder_data->IdxOnFl) &&
         (ladder_data->ArrearsFl == 'N') )
    {
        DR_Error("Fl reset in advance invalid for zero weight float index !");
        goto RETURN;
    }

    /* Nb of stepup dates can must be at least 1. */
    if (ladder_data->NbStepUpFl < 1)
    {
        DR_Error ("At least one stepup date is required!");
        goto RETURN;
    }                                               
    
    /* Stepup date format, ascending order, rate order */
    for (i = 0; i < ladder_data->NbStepUpFl; i++)
        if (Dateok(ladder_data->FlDate[i]))
        {
            DR_Error ("Incorrect format for stepup date!");
            goto RETURN;
        }

    for (i = 1; i < ladder_data->NbStepUpFl; i++)
        if (ladder_data->FlDate[i] <= ladder_data->FlDate[i-1])
        {
            DR_Error ("Stepup dates must be entered in ascending order!");
            goto RETURN;
        }

    for (i = 0; i < ladder_data->NbStepUpFl; i++)
        if (ladder_data->CapFl[i] < ladder_data->FloorFl[i])
        {
            DR_Error ("Cap rate must be bigger than floor rate !");
            goto RETURN;
        }
    
    /* Stepup dates only possible on ladder leg dates */
    if (DrDatesInSchedule(ladder_data->NbStepUpFl,
                          ladder_data->FlDate,
                          ladder_data->AccStDate,
                          ladder_data->MatDate,
                          ladder_data->PayFreqFl,
                          ladder_data->StubConv) == FAILURE)
    {
        DR_Error("Stepup dates are not a subset of float pmt dates.\n");
        goto RETURN;
    }
    if (ladder_data->FlDate[ladder_data->NbStepUpFl-1] 
        == ladder_data->MatDate)
    {
        DR_Error("Deal cannot step up on maturity date.\n");
        goto RETURN;
    }

    /* For flt coupon payments, the usual frequencies */
    if ((ladder_data->PayFreqFl != 'M') && 
        (ladder_data->PayFreqFl != 'Q') && 
        (ladder_data->PayFreqFl != 'S') &&
        (ladder_data->PayFreqFl != 'A') )   
    {
        DR_Error ("Invalid entry for swap flt payment frequency!");
        goto RETURN;
    } 
                                                  
    /* For flt coupon day count convention */
    if ((ladder_data->DayCountFl != '0') && 
        (ladder_data->DayCountFl != '5') && 
        (ladder_data->DayCountFl != '3') )   
    {
        DR_Error ("Invalid entry for flt swap day count convention!");
        goto RETURN;
    } 

    /* Flt frequency must be >= ladder frequency for zero version */
    if (ladder_data->SoZ == 'Z')
    {
        if (  Conv_Freq(ladder_data->PayFreqSt) 
            > Conv_Freq(ladder_data->PayFreqFl))
        {
            DR_Error ("Ladder leg cannot be more frequent than floating leg!\n"); 
            goto RETURN;
        } 
    }

    /* MODEL CHECKS */

    /* Diffused curve number */
    if ((tree_data->CvDiff != 0) &&
        (tree_data->CvDiff != 1) &&
        (tree_data->CvDiff != 2) )
    {
        DR_Error("Diffused curve number must be 0, 1, 2 !");
        goto RETURN;
    } 

    /* Index curve number */
    if ((tree_data->CvIdx1 != 0) &&
        (tree_data->CvIdx1 != 1) &&
        (tree_data->CvIdx1 != 2) )
    {
        DR_Error("Index1 curve number must be 0, 1, 2 !");
        goto RETURN;
    } 
    if ((tree_data->CvIdx2 != 0) &&
        (tree_data->CvIdx2 != 1) &&
        (tree_data->CvIdx2 != 2) )
    {
        DR_Error("Index2 curve number must be 0, 1, 2 !");
        goto RETURN;
    } 

    /* All three curves have different indexes in tree */
    if ((tree_data->CvDiff == tree_data->CvIdx1) || 
        (tree_data->CvDiff == tree_data->CvIdx2) || 
        (tree_data->CvIdx1 == tree_data->CvIdx2))
    {
        DR_Error("Curve indexing is wrong in the tree !");
        goto RETURN;
    } 

    /* Discount curve number */
    if ((tree_data->CvDisc != 0) &&
        (tree_data->CvDisc != 1) &&
        (tree_data->CvDisc != 2) )
    {
        DR_Error("Discount curve number must be 0, 1, 2 !");
        goto RETURN;
    } 
    

    if (ladder_data->NbStates < 3) 
    {
        DR_Error ("Number of state variables must be >= 3");
        goto RETURN;
    }
   

    if (tree_data->NbSigmaMax != ladder_data->NbStDev)
    {
        DR_Error("The state variable stdev (%d) != tree stdev (%d)!",
                 ladder_data->NbStDev,
                 tree_data->NbSigmaMax);
        goto RETURN;
    }

    if (ladder_data->NbStDev < 1)
    {
        DR_Error("Can't state variable at less than one std devs!");
        goto RETURN;
    }
    
            
    status = SUCCESS;

    RETURN:

    Free_DR_Array (DateList, LONG, 0, NbDates-1);

    return (status);

}  /* Ladder_Check */


/*****  Print_Ladder   *************************************************/
/**
*       Print debug information in an ascii file.
*/
int     Print_Ladder 
            (T_CURVE          *t_curve,          /**< (O) Zero curve data    */
             MKTVOL_DATA      *mktvol_data,      /**< (O) Vol  data          */
             LADDER_DATA      *ladder_data,      /**< (O) Deal data          */
             FIX3_TREE_DATA        *tree_data)        /**< (O) Tree data          */
{

    int     i, j, k;  
    int     status = FAILURE;               /* Status = FAILURE initially  */
    char    ErrorMsg[MAXBUFF];
    char    LegName[MAXBUFF];
    FILE    *stream = NULL;
    double  days;
    int     IdxReset, p;

    char    VolIndex[MAXBUFF];
    int     idxMat;


    stream = fopen ("TERM.prn", "w");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, 
                 "Could not open file TERM.prn! "
                 "(Print_Ladder)");                 
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* DEAL DATA FIRST */

    fprintf (stream, "### LADDER DEAL INFO ###\n");

    /* OPTION */

    strcpy (LegName, "# Option :");

    /* Long or short */
    fprintf (stream, "%s Long or short \n%c \n", 
                     LegName,
                     (ladder_data->LoS));

    /* Option style */
    fprintf (stream, "%s Exercise style \n%c \n", 
                     LegName,
                     (ladder_data->Style));

    /* Number of exercise date, dates and strikes */
    fprintf (stream, "%s Number of exercises \n%d \n", 
                     LegName, 
                     (ladder_data->NbExer));
    
    fprintf (stream, "%s Date and strike \n", 
                     LegName);
    for (i = 0; i < ladder_data->NbExer; i++)
    {
        fprintf(stream, "%8ld %12.6f\n",
                        YMDDateFromIRDate(ladder_data->ExerDate[i]),
                        100 *(ladder_data->Strike[i]));
    }   

    /* exercise out statistics */
    fprintf (stream, "%s Event statistics schedule style I,D,W,M,Q,S,A\n%1c \n", 
                     LegName,
                     ladder_data->OptStatStyle);

    fprintf (stream, "%s Event statistics number of dates \n%ld \n", 
                     LegName,
                     ladder_data->OptNbStats);

    fprintf (stream, "%s Event statistics dates\n", LegName);

    for (i = 0; i < ladder_data->OptNbStats; i++)
        fprintf(stream, "%8ld\n", YMDDateFromIRDate(ladder_data->OptStatDates[i]));


    /* SWAP */

    strcpy (LegName, "# Swap :");

    /* Notional */
    fprintf (stream, "%s Notional & sign \n%-14.4f\n", 
                     LegName, 
                     (ladder_data->Notional) * (ladder_data->NotionalSign));

    /* Original percentage */
    fprintf (stream, "%s Original percentage \n%-12.6f\n", 
                     LegName, 
                     100 * ladder_data->OrgNotPerc);

    /* Number of amort */
    fprintf (stream, "%s Number of amort dates \n%d \n", 
                     LegName, 
                     (ladder_data->NbAmort));
    
    fprintf (stream, "%s Date and rate \n", 
                     LegName);
    for (i = 0; i < ladder_data->NbAmort; i++)
    {
        fprintf(stream, "%8ld %12.6f\n",
                        YMDDateFromIRDate(ladder_data->AmortDate[i]),
                        100 *(ladder_data->Amort[i]));
    }   

    /* Stub conv */
    fprintf (stream, "%s Stub conv \n%c \n", 
                     LegName, 
                     (ladder_data->StubConv));

    /* Swap or zero coupon */
    fprintf (stream, "%s Swap or zero coupon \n%c \n", 
                     LegName, 
                     (ladder_data->SoZ));

    /* Smoothing */
    fprintf (stream, "%s Smoothing \n%c \n", 
                     LegName, 
                     (ladder_data->Smoothing));

    /* Number of rib range dates plus dates, low and high barriers */
    fprintf(stream, "# Rib: Nb of rib observation dates\n");
    fprintf(stream, "%d\n", ladder_data->NbRibObsDates);

    fprintf(stream, "# Rib: observation, obs effective date, low and high "
            "barriers, inside and outside weights\n");
    for (i = 0; i < ladder_data->NbRibObsDates; i++)
    {
        fprintf(stream, "%ld\t %ld\t %lf\t %lf\t %lf\t %lf\n",
                YMDDateFromIRDate(ladder_data->RibObsDate[i]),
                YMDDateFromIRDate(ladder_data->RibObsEffDate[i]),
                ladder_data->RibLoBarrier[i] * 100.,
                ladder_data->RibHiBarrier[i] * 100.,
                ladder_data->RibInRangeWeight[i],
                ladder_data->RibOutRangeWeight[i]);    
    }


    fprintf(stream, "# Rib: Smoothing (Y or N)\n");
    fprintf(stream, "%c\n", ladder_data->RibSmoothing);

    /* Two lines for rib observation indices, curves and weights. */
    fprintf(stream, "# Rib: Observation index maturity, freq, DCC, curve id, "
            "and weight\n");

    for(j = 0; j < 2; j++)
    {
        fprintf(stream, "%-3d\t %c\t %c\t %1d\t %lf\n",
                ladder_data->RibIdxMat[j],
                ladder_data->RibIdxFreq[j],
                ladder_data->RibIdxDCC[j],
                ladder_data->RibIdxIoD[j],
                ladder_data->RibIdxWeight[j]);
    }

    /* STICKY LEG */
    
    /* St index */
    
     fprintf(stream, "%s St indices (mat, freq, dcc, curve)\n", LegName);
     for (i = 0 ; i < 2; i++)
     {
         fprintf(stream,"%-3d %1c %1c %1d %lf\n",           
                    ladder_data->IdxMatSt[i],
                    ladder_data->IdxFreqSt[i],
                    ladder_data->IdxBaseSt[i],
                    ladder_data->IdxIoDSt[i],
                    ladder_data->IdxObsWeightSt[i]);
    }
    

    /* Arrears reset ladder */
    fprintf (stream, "%s Arrears ladder \n%c \n", 
                     LegName, 
                     (ladder_data->ArrearsSt));

    /* Comp ladder */
    fprintf (stream, "%s Comp ladder \n%c \n", 
                     LegName, 
                     (ladder_data->CompSt));

    /* Previous refixes of ladder leg. */
    fprintf (stream, "%s Number of previous ladder fixings \n%d \n", 
                     LegName, 
                     (ladder_data->NbFixing));

    fprintf (stream, "%s Previous ladder fixings \n", LegName);
    for (i=0; i<ladder_data->NbFixing; i++)
    {
        fprintf(stream, "%8ld %14.4f %14.4f %14.4f\n",
                YMDDateFromIRDate(ladder_data->FixingDate[i]),
                100 * ladder_data->FixingSt[0][i],
                100 * ladder_data->FixingSt[1][i],
                100 * ladder_data->FixingRibPerc[i]);
    }

    /* First lookback */
    fprintf (stream, "%s First lookback value\n%.6f\n", 
                     LegName, 
                     100 * ladder_data->FirstLevel);

    /* St date and info.*/
    fprintf (stream, "%s Number of ladder dates \n%d \n", 
                     LegName, 
                     (ladder_data->NbStepUpSt));

    fprintf (stream, "%s St dates and info \n", LegName);
    for (i=0; i<ladder_data->NbStepUpSt; i++)
    {
        fprintf(stream, "%8ld %12.6f %12.6f %12.6f %12.6f %12.6f "
                        "%12.6f %12.6f %12.6f %12.6f %12.6f %12.6f "
                        "%12.6f %12.6f %12.6f\n",
                YMDDateFromIRDate(ladder_data->StDate[i]),
                100 * ladder_data->FloorSt[i],
                100 * ladder_data->CapSt[i],
                100 * ladder_data->DownRateSt[i],
                100 * ladder_data->MidRateSt[i],
                100 * ladder_data->UpRateSt[i],
                100 * ladder_data->BarrierLo[i],
                100 * ladder_data->BarrierHi[i],
                      ladder_data->IdxWeightSt[0][i],
                      ladder_data->IdxWeightSt[1][i],
                100 * ladder_data->SprdSt[i],
                      ladder_data->StickyCoef[i],
                      ladder_data->Leverage[i],
                100 * ladder_data->IdxFloor[i],
                100 * ladder_data->IdxCap[i]);
    }

    /* Additive or multiplicative flag */
    fprintf (stream, "%s St add / mult \n%c \n", 
                     LegName, 
                     (ladder_data->AoM));

    
    /* St freq */
    fprintf (stream, "%s St freq \n%c \n", 
                     LegName, 
                     (ladder_data->PayFreqSt));

    /* St dcc */
    fprintf (stream, "%s St dcc \n%c \n", 
                     LegName, 
                     (ladder_data->DayCountSt));

    /* FLOATING LEG */

    /* Fl index */
    fprintf(stream, "%s Fl idx (mat, freq, dcc, curve) \n%-3d %1c %1c %1d %lf\n",
                    LegName, 
                    ladder_data->IdxMatFl,
                    ladder_data->IdxFreqFl,
                    ladder_data->IdxBaseFl,
                    ladder_data->IdxIoDFl,
                    ladder_data->IdxWeightFl);

    /* Arrears reset float */
    fprintf (stream, "%s Arrears float \n%c \n", 
                     LegName, 
                     (ladder_data->ArrearsFl));

    /* Comp float */
    fprintf (stream, "%s Comp float \n%c \n", 
                     LegName, 
                     (ladder_data->CompFl));

    /* Previous refixes of float leg. */
    fprintf (stream, "%s Previous float fixing \n%-12.6f \n", 
                     LegName, 
                     100 * (ladder_data->FixingFl));

    /* Fl date and info.*/
    fprintf (stream, "%s Number of float dates \n%d \n", 
                     LegName, 
                     (ladder_data->NbStepUpFl));

    fprintf (stream, "%s Fl dates and info \n", LegName);
    for (i=0; i<ladder_data->NbStepUpFl; i++)
    {
        fprintf(stream, "%8ld %12.6f %12.6f %12.6f\n",
                YMDDateFromIRDate(ladder_data->FlDate[i]),
                100 * ladder_data->FloorFl[i],
                100 * ladder_data->CapFl[i],
                100 * ladder_data->UpRateFl[i]);
    }
    
    /* Fl freq */
    fprintf (stream, "%s Fl freq \n%c \n", 
                     LegName, 
                     (ladder_data->PayFreqFl));

    /* Fl dcc */
    fprintf (stream, "%s Fl dcc \n%c \n", 
                     LegName, 
                     (ladder_data->DayCountFl));

   /*
     * Model specific data
     */

    /* First index to be used for volatility calibration */
    if (!mktvol_data->CalibFlag)
    {
        strcpy(VolIndex, "nil");
    }
    else
    {
        idxMat = Months360(mktvol_data->SwapSt[0], mktvol_data->SwapMat[0]);

        /* Fixed maturity date */
        if (ABS(Daysact(mktvol_data->SwapMat[0], 
                        mktvol_data->SwapMat[mktvol_data->NbVol-1])) < 87) 
        {
            sprintf(VolIndex, "%dyFix", idxMat / 12);
        }
        else if (idxMat > 12)  /* CMS */
        {
            sprintf(VolIndex, "%dyCms", idxMat / 12);
        }
        else   /* base vol */
        {
            sprintf(VolIndex, "%dm", idxMat);
        }
        
    }

    fprintf(stream, "# Model: Index for first volatility calibration\n");
    fprintf(stream, "%s\n", VolIndex);

    fprintf(stream, "# Model: Index for second volatility calibration\n");
    fprintf(stream, "%s\n", VolIndex);


    /* Diffuse and discount curve */
    fprintf(stream, "# Model: Curve to discount\n");
    fprintf(stream, "%d\n", tree_data->CvDisc);

    /* Number of standard deviations at which to cut the tree */
    fprintf(stream, "# Model: Nb of standard deviations used to cut the tree\n");
    fprintf(stream, "%d\n", tree_data->NbSigmaMax);

    /* Number of state variables */
    fprintf (stream, "Model: Number of state variables \n%-3d \n", 
                     (ladder_data->NbStates));

    /* Number of st dev for state variable */
    fprintf (stream, "Model: Number of st dev for state var \n%-3d \n", 
                     (ladder_data->NbStDev));


    
    /* Number of periods per year in the tree */
    fprintf(stream, "# Model: Nb of periods per year used to calculate the time step\n");
    fprintf(stream, "%d\n", tree_data->Ppy);

    /* Smile    */
    fprintf(stream, "# Model: Smile parameters (qL, qR, FwdShift, CetIter)\n");
    fprintf(stream, "%lf\t %lf\t %lf %d\n", 
            1.0 - mktvol_data->QLeft,
            1.0 - mktvol_data->QRight,
            mktvol_data->FwdShift,
            mktvol_data->CetNbIter);

    /* Number of factors */
    fprintf(stream, "# Model: Number of factors\n");
    fprintf(stream, "%d\n", tree_data->NbFactor);

    /* Factor weight */
    fprintf(stream, "# Model: Factor weights (one per factor)\n");
    for (i=0; i<tree_data->NbFactor; i++)
        fprintf(stream, "%lf\t", mktvol_data->Alpha[i]);

    fprintf(stream, "\n");

    /* Mean reversion */
    fprintf(stream, "# Model: Mean reversion coefficient (one per factor)\n");
    for (i=0; i<tree_data->NbFactor; i++)
        fprintf(stream, "%lf\t", mktvol_data->Beta[i]);

    fprintf(stream, "\n");

    /* Correlations */
    fprintf(stream, "# Model: correlation (0, 1 or 3 numbers)\n");
    for (i=0; i<(tree_data->NbFactor)*(tree_data->NbFactor-1)/2; i++)
        fprintf(stream, "%lf\t", mktvol_data->Rho[i]);

    fprintf(stream, "\n");

    /* Backbone */
    fprintf(stream, "# Model: Backbone\n");
    fprintf(stream, "%lf\n", 1. - mktvol_data->Bbq);




    /* ------------------------------------------------------------------*/
    /*   STANDARD ENVIRONMENT ELEMENTS (i.e. zero.dat, basevol.dat etc.) */
    /*-------------------------------------------------------------------*/

    /* Zero curves */
    for (j = 0; j < 3; j++)
    {
        if (j==0)      strcpy (LegName, "zero.dat");
        else if (j==1) strcpy (LegName, "disczero.dat");
        else           strcpy (LegName, "riskzero.dat");
        fprintf (stream, "\n### ZERO CURVE NB %1d (%s) ### \n",j,LegName); 
        EslPrintZeroCurve(&t_curve[j], stream);
    }

    /* Market vol data */

    fprintf (stream, "\n### MARKET VOL DATA ### \n"); 

    fprintf (stream, "# Calib flag \n%1d \n\n", (mktvol_data->CalibFlag));

    if (mktvol_data->CalibFlag)
    {
        fprintf (stream, "# Base date \n%8ld \n", YMDDateFromIRDate(mktvol_data->BaseDate));
        fprintf (stream, "# Nb of vol points \n%-3d \n", (mktvol_data->NbVol));
        fprintf (stream, "# Vol dates and rates \n");
        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            fprintf (stream, "%8ld %12.6f \n",
                             YMDDateFromIRDate(mktvol_data->VolDate[i]),
                             100*(mktvol_data->Vol[i]));
        }
        fprintf (stream, "# Freq \n%1c \n", (mktvol_data->Freq));
        fprintf (stream, "# Dcc \n%1c \n", (mktvol_data->DCC));
        fprintf (stream, "# Swap starts and mats \n");
        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            fprintf (stream, "%8ld %8ld \n",
                             YMDDateFromIRDate(mktvol_data->SwapSt[i]),
                             YMDDateFromIRDate(mktvol_data->SwapMat[i]));
        }
        fprintf (stream, "# Skip flag \n%1d \n", (mktvol_data->SkipFlag));
    }

    fprintf (stream, "# Ql, Qr, Fs, CetNb \n%12.6f %12.6f %12.6f %2d\n",
                      (mktvol_data->QLeft), 
                      (mktvol_data->QRight), 
                      (mktvol_data->FwdShift),
                      (mktvol_data->CetNbIter));
    fprintf (stream, "# Alpha \n%12.6f %12.6f %12.6f \n",
                      (mktvol_data->Alpha[0]), 
                      (mktvol_data->Alpha[1]), 
                      (mktvol_data->Alpha[2]));
    fprintf (stream, "# Beta \n%12.6f %12.6f %12.6f \n",
                      (mktvol_data->Beta[0]), 
                      (mktvol_data->Beta[1]), 
                      (mktvol_data->Beta[2]));
    fprintf (stream, "# Rho \n%12.6f %12.6f %12.6f \n",
                      (mktvol_data->Rho[0]), 
                      (mktvol_data->Rho[1]), 
                      (mktvol_data->Rho[2]));

    /* Tree data */
    
    fprintf (stream, "\n### TREE DATA ### \n"); 
    fprintf (stream, "# Ppy \n%-3d \n", (tree_data->Ppy));
    fprintf (stream, "# Nb of factor \n%1d \n", (tree_data->NbFactor));
    fprintf (stream, "# Nb of sigma max \n%1d \n", (tree_data->NbSigmaMax));
    fprintf (stream, "# Diffused curve \n%1d \n", (tree_data->CvDiff));
    fprintf (stream, "# Idx1 curve \n%1d \n", (tree_data->CvIdx1));
    fprintf (stream, "# Idx2 curve \n%1d \n", (tree_data->CvIdx2));
    fprintf (stream, "# Discount curve \n%1d \n", (tree_data->CvDisc));

    /*************************************************************************/
    
    /* Schedule of the OPTION */
    fprintf (stream, "\nLADDER SCHEDULE:\n");
    fprintf (stream,
             "\nTimePt   Date     | Exer    Strike \n");
    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[0][i]) 
        {
            fprintf (stream,
                    "[%4d]   %8ld |  (%d)  %8.6f \n",
                    i,
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    tree_data->TPtype[0][i],
                    tree_data->CritDate[0][i].Value[0]);
        }
    
    }  /* for i */

    fprintf (stream, "\n");

    /* Schedule of the underlying SWAP (STICKY LEG) */
    fprintf (stream, "\nSTICKY LEG RESET SCHEDULE:\n");
    fprintf (stream,   
             "\nTimePt   Date     | Reset    Floor       Cap "
             " DownRate   MidRate    UpRate    BarrLo    BarrHi " 
             " Idx Wgt1   Idx Wgt2   Spread    St Wgt    Lever  " 
             " Idx Cap    Idx Floor  " 
             "    Outs       DCF   PmtDate \n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[2][i])
        {
            fprintf (stream,
                    "[%4d]   %8ld |  (%d) %9.6f %9.6f "
                    "%9.6f %9.6f %9.6f %9.6f %9.6f "
                    "%9.6f %9.6f %9.6f %9.6f %9.6f "                     
                    "%9.6f %9.6f "
                    "%14.2f %9.6f  %8ld\n",
                    i,
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    tree_data->TPtype[2][i],
                    tree_data->CritDate[2][i].Value[0],
                    tree_data->CritDate[2][i].Value[1],
                    tree_data->CritDate[6][i].Value[0],
                    tree_data->CritDate[6][i].Value[1],
                    tree_data->CritDate[6][i].Value[2],
                    tree_data->CritDate[6][i].Value[3],
                    tree_data->CritDate[6][i].Value[4],
                    tree_data->CritDate[7][i].Value[0],
                    tree_data->CritDate[7][i].Value[1],
                    tree_data->CritDate[7][i].Value[2],
                    tree_data->CritDate[7][i].Value[3],
                    tree_data->CritDate[8][i].Value[0],
                    tree_data->CritDate[8][i].Value[1],
                    tree_data->CritDate[8][i].Value[2],
                    tree_data->CritDate[2][i].Value[2],
                    tree_data->CritDate[2][i].Value[3],
                    YMDDateFromIRDate(tree_data->CritDate[2][i].SuppDate[0]));
        }
    }

    fprintf (stream, "\n");

    /* Schedule of the underlying SWAP (FLT LEG) */
    fprintf (stream, "\nFLOATING LEG SCHEDULE:\n");
    fprintf (stream,   
             "\nTimePt   Date     | Reset    Floor       Cap    StepUp"
             "            Outs       DCF   PmtDate \n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[4][i])
        {
            fprintf (stream,
                    "[%4d]   %8ld |  (%d) %9.6f %9.6f %9.6f  "
                    "%14.2f %9.6f  %8ld\n",
                    i,
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    tree_data->TPtype[4][i],
                    tree_data->CritDate[4][i].Value[0],
                    tree_data->CritDate[4][i].Value[1],
                    tree_data->CritDate[4][i].Value[2],
                    tree_data->CritDate[4][i].Value[3],
                    tree_data->CritDate[4][i].Value[4],
                    YMDDateFromIRDate(tree_data->CritDate[4][i].SuppDate[0]));
        }
    }

    fprintf (stream, "\n");

    /* State variable dates */
    fprintf (stream, "\nSTATE VARIABLE DISTRIBUTION    :\n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[5][i])
        {
            IdxReset = (int) (tree_data->CritDate[5][i]).Value[0];

            fprintf (stream,
                    "[%3d]   %8ld  |  ",
                    IdxReset,
                    YMDDateFromIRDate(tree_data->TPDate[i]));

            for (p=0; p<ladder_data->NbStates; p++)
                fprintf(stream, "%9.6f  ", ladder_data->State[IdxReset][p]);

            fprintf(stream, "\n");
        }
    }


    /* Known ladder coupon payments */
    fprintf (stream, "\nKNOWN STICKY COUPON PAYMENTS:\n");
    fprintf (stream, "\nTimePt   Date     | Known       Amount       Rate"
                     "        DCF          Outs    Factor\n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[1][i])
        {
            fprintf (stream,
                    "[%4d]   %8ld |  (%d)  %12.4f  %9.6f  %9.6f  %12.4f %9.6f \n",
                    i,
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    tree_data->TPtype[1][i],
                    tree_data->CritDate[1][i].Value[0],
                    tree_data->CritDate[1][i].Value[1],
                    tree_data->CritDate[1][i].Value[2],
                    tree_data->CritDate[1][i].Value[3],
                    tree_data->CritDate[1][i].Value[4]);
        }
    }

    fprintf (stream, "\n");

    /* Known floating coupon payments */
    fprintf (stream, "\nKNOWN FLOATING COUPON PAYMENTS:\n");
    fprintf (stream, "\nTimePt   Date     | Known       Amount       Rate"
                     "        DCF          Outs\n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[3][i])
        {
            fprintf (stream,
                    "[%4d]   %8ld |  (%d)  %12.4f  %9.6f  %9.6f  %12.4f\n",
                    i,
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    tree_data->TPtype[3][i],
                    tree_data->CritDate[3][i].Value[0],
                    tree_data->CritDate[3][i].Value[1],
                    tree_data->CritDate[3][i].Value[2],
                    tree_data->CritDate[3][i].Value[3]);
        }
    }

    fprintf (stream, "\n");

    /* Zero bank schedules for the product */

    for (j=0; j<3; j++)
    {
        fprintf (stream, "\nZERO BANK (%d) SCHEDULE:\n", j);
        fprintf (stream,"\nTimePt   Date       ZeroMat   EarliestUse\n");
        for (i = 0; i <= tree_data->NbTP; i++)
        {
            if (tree_data->TPtype[ZbkEVENT+j][i])
            {
                fprintf (stream,
                        "[%4d]   %8ld     (%d)        %8ld\n",
                        i,
                        YMDDateFromIRDate(tree_data->TPDate[i]),
                        tree_data->TPtype[ZbkEVENT+j][i],
                        YMDDateFromIRDate(tree_data->CritDate[ZbkEVENT+j][i].SuppDate[0]));
            }
        }
        fprintf(stream, "\n");
        fprintf (stream, "Max Bank Size: %d\n\n", 
                         tree_data->NbZeros[ZbkEVENT+j]);

    } /* for each zero bank j */
    
    fprintf (stream, "\nWidth1 Width2 Width3 \n");

    fprintf(stream, "%5d %5d %5d \n\n",
            tree_data->Width[0],
            tree_data->Width[1],
            tree_data->Width[2]);

    fprintf (stream, "Node       Date     days  forward   "
                     "zero0  discount0  "
                     "zero1  discount1  "
                     "zero2  discount2  "
                     " Drift    Aweigths\n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        days = Daysact (tree_data->TPDate[0], tree_data->TPDate[i]);

        fprintf(stream, "[%4d] \t%8ld  %5.0f  %7.4f  %7.4f   %8.6f%7.4f   "
                        "%6.4f  %7.4f   %6.4f  %9.6f  %7.4f %7.4f %7.4f %7.4f "
                        "%7.4f %7.4f \n",
                i,
                YMDDateFromIRDate(tree_data->TPDate[i]),
                days,
                tree_data->FwdRate[0][i] * 100.,
                tree_data->ZeroRate[0][i] * 100.,
                tree_data->ZeroCoupon[0][i],
                tree_data->ZeroRate[1][i] * 100.,
                tree_data->ZeroCoupon[1][i],
                tree_data->ZeroRate[2][i] * 100.,
                tree_data->ZeroCoupon[2][i],
                tree_data->ZCenter[i],
                tree_data->Aweight[0][i],
                tree_data->Aweight[1][i],
                tree_data->Aweight[2][i],
                tree_data->Aweight[3][i],
                tree_data->Aweight[4][i],
                tree_data->Aweight[5][i]);
    }                                               

    for (k = 0; k < 3; k++)
    {
        fprintf (stream, "\nCURVE NB %1d\n", k);
        EslPrintZeroCurve(&t_curve[k], stream);
    }                                               
      
    status = SUCCESS;

RETURN:
 
    if (stream != NULL)
    {
        fclose (stream);
    }
    
    return (status);

}  /* Print_Ladder */



/*****  Print_Flows_Deal   *************************************************/
/*
*       Print debug information in an ascii file.
*/
int     Print_Flows_Deal 
            (MKTVOL_DATA      *mktvol_data,
             LADDER_DATA      *ladder_data,
             FIX3_TREE_DATA   *tree_data)
{
    char const* routine = "Print_Flows_Deal";

    char        buff[128];

    int         i, j, size;  
    char const* LegName;

    FILE* stream = fopen ("deal.dat", "w");
    if (stream == NULL)
    {
        DR_Error("%s: could not open run.log", routine);
        return FAILURE;
    }

    /* DEAL DATA FIRST */
    fprintf (stream, "### LADDERFLOWS DEAL INFO ###\n");

    /*--- OPTION ---*/
    LegName = "# Option:";

    /* Long or short */
    fprintf (stream, "%s (L)ong or (S)hort \n%c \n", 
                     LegName,
                     ladder_data->LoS);

    /* Option style */
    fprintf (stream, "%s Exercise style \n%c \n", 
                     LegName,
                     ladder_data->Style);

    /* Number of exercise date, dates and strikes */
    fprintf (stream, "%s Number of exercises \n%d \n", 
                     LegName, 
                     ladder_data->NbExer);
    
    fprintf (stream, "%s Not, Eff, Exe, Pct Strike \n", 
                     LegName);
    for (i = 0; i < ladder_data->NbExer; i++)
    {
        fprintf(stream, "%8ld %8ld %8ld %12.2f\n",
                        YMDDateFromIRDate(ladder_data->ExerDate[i]),
                        YMDDateFromIRDate(ladder_data->ExerDate[i]),
                        YMDDateFromIRDate(ladder_data->ExerDate[i]),
                        100.0 * ladder_data->Strike[i]);
    }   

    /*--- Amortization ---*/
    LegName = "# Amortization:";

    /* Number of amortizations */
    fprintf (stream, "%s Number of amortizations \n%d \n", 
                     LegName, 
                     ladder_data->NbAmort);
    
    fprintf (stream, "%s Date, Amount \n", 
                     LegName);
    for (i = 0; i < ladder_data->NbAmort; i++)
    {
        fprintf(stream, "%8ld %12.2f\n",
                        YMDDateFromIRDate(ladder_data->AmortDate[i]),
                        ladder_data->Amort[i] * ladder_data->Notional);
    }   

    /* Coupon style, current notional */
    fprintf (stream, "%s Coupon style, current notional \n%c %f \n", 
                     LegName,
                     ladder_data->SoZ,
                     ladder_data->NotionalSign * ladder_data->Notional);

    /* Smoothing */
    fprintf (stream, "%s Binary smoothing \n%c \n", 
                     LegName,
                     ladder_data->Smoothing);

    /* Number of rib range dates plus dates, low and high barriers */
    fprintf(stream, "# Rib: Nb of rib observation dates\n");
    fprintf(stream, "%d\n", ladder_data->NbRibObsDates);

    fprintf(stream, "# Rib: observation, obs effective date, low and high "
            "barriers, inside and outside weights\n");
    for (i = 0; i < ladder_data->NbRibObsDates; i++)
    {
        fprintf(stream, "%ld\t %ld\t %lf\t %lf\t %lf\t %lf\n",
                YMDDateFromIRDate(ladder_data->RibObsDate[i]),
                YMDDateFromIRDate(ladder_data->RibObsEffDate[i]),
                ladder_data->RibLoBarrier[i] * 100.,
                ladder_data->RibHiBarrier[i] * 100.,
                ladder_data->RibInRangeWeight[i],
                ladder_data->RibOutRangeWeight[i]);    
    }


    fprintf(stream, "# Rib: Smoothing (Y or N)\n");
    fprintf(stream, "%c\n", ladder_data->RibSmoothing);

    /* Two lines for rib observation indices, curves and weights. */
    fprintf(stream, "# Rib: Observation index maturity, freq, DCC, curve id, "
            "and weight\n");

    for(j = 0; j < 2; j++)
    {
        fprintf(stream, "%-3d\t %c\t %c\t %1d\t %lf\n",
                ladder_data->RibIdxMat[j],
                ladder_data->RibIdxFreq[j],
                ladder_data->RibIdxDCC[j],
                ladder_data->RibIdxIoD[j],
                ladder_data->RibIdxWeight[j]);
    }

    /* Past observations in range */
    fprintf (stream, "%s past obs in range percentage\n%.6f\n", LegName,
             100 * ladder_data->RibPastObsPerc);
   
    /*--- Ladder Leg ---*/
    LegName = "# Ladder Leg:";

    /* Number of past refixes */
    fprintf (stream, "%s Number of ladder index rate refixes \n%d \n", 
                     LegName, 
                     ladder_data->NbFixing);
    
    fprintf (stream, "%s Indices rate refix date, rate 0, rate 1\n", 
                     LegName);
    for (i = 0; i < ladder_data->NbFixing; i++)
    {
        fprintf(stream, "%8ld %12.2f %12.2f\n",
                        YMDDateFromIRDate(ladder_data->FixingDate[i]),
                        ladder_data->FixingSt[0][i] * 100,
                        ladder_data->FixingSt[1][i] * 100);
    }   

    /* First ladder rate lavel */
    fprintf (stream, "%s Ladder rate initial refix \n%.6f \n", 
                     LegName, 
                     100 * ladder_data->FirstLevel);
    
    /* Date adjustment offset */
    fprintf (stream, "%s Reset adjusted relatice to 'S'tart, 'E'nd of period or 'U'ndefined\n%c\n",
                     LegName, 
                     'U');
    
    /* Number of dates and levels */
    size = 0;
    for (i = 0; i <= tree_data->NbTP; i++)
        if (tree_data->TPtype[2][i])
            ++size;

    fprintf (stream, "%s Number of dates and levels \n%d \n", 
                     LegName, size);
    
    /* Dates and levels */
    fprintf (stream, "%s res, eff, acs, ace, pay  "
             "dcf, notional, floor, cap, down, mid, up, low, high, weight 0, weight 1, spread, ladder\n", 
             LegName);
    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (!tree_data->TPtype[2][i])
            continue;

        fprintf(stream, "%8ld %8ld %8ld %8ld %8ld "
                        "%12.8f %12.6f %12.6f %12.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n",
                        YMDDateFromIRDate(tree_data->TPDate[i]), // res
                        YMDDateFromIRDate(tree_data->TPDate[i]), // eff
                        YMDDateFromIRDate(tree_data->CritDate[2][i].SuppDate[1]), // acs
                        YMDDateFromIRDate(tree_data->CritDate[2][i].SuppDate[0]), // ace
                        YMDDateFromIRDate(tree_data->CritDate[2][i].SuppDate[0]), // pay
                        tree_data->CritDate[2][i].Value[3], // dcf
                        tree_data->CritDate[2][i].Value[2] * ladder_data->NotionalSign * ladder_data->Notional, // not
                        tree_data->CritDate[2][i].Value[0] * 100, // floor
                        tree_data->CritDate[2][i].Value[1] * 100, // cap
                        tree_data->CritDate[6][i].Value[0] * 100, // down
                        tree_data->CritDate[6][i].Value[1] * 100, // mid
                        tree_data->CritDate[6][i].Value[2] * 100, // up
                        tree_data->CritDate[6][i].Value[3] * 100, // low
                        tree_data->CritDate[6][i].Value[4] * 100, // high
                        tree_data->CritDate[7][i].Value[0], // weight 0
                        tree_data->CritDate[7][i].Value[1], // weight 1
                        tree_data->CritDate[7][i].Value[2] * 100, // spread
                        tree_data->CritDate[7][i].Value[3]); // ladder
    }   

    /* Add/mul flag */
    fprintf (stream, "%s Additive/multiplicative flag\n%c\n", LegName, ladder_data->AoM);

    /* Ladder index */
    fprintf (stream, "%s Indices maturity, frequency, DCC, curve, weight\n"
                     "%d %c %c %d %lf      %d %c %c %d %lf\n", 
                     LegName, 
                     ladder_data->IdxMatSt[0],  // maturity
                     ladder_data->IdxFreqSt[0], // frequency
                     ladder_data->IdxBaseSt[0], // dcc
                     ladder_data->IdxIoDSt[0],  // curve
                     ladder_data->IdxObsWeightSt[0],  // weight

                     ladder_data->IdxMatSt[1],  // maturity
                     ladder_data->IdxFreqSt[1], // frequency
                     ladder_data->IdxBaseSt[1], // dcc
                     ladder_data->IdxIoDSt[1],  // curve
                     ladder_data->IdxObsWeightSt[1]); // weight

    /*--- Funding leg ---*/
    LegName = "# Funding Leg:";

    /* Number of past refixes */
    size = 0;
    for (i = 0; i <= tree_data->NbTP; i++)
        if (tree_data->TPtype[4][i])
            ++size;

    fprintf (stream, "%s Number of past index refixes \n%d \n", 
                     LegName, size); 
    
    fprintf (stream, "%s Past index refix Date, Rate \n", 
                     LegName);
    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[4][i])
        {
            fprintf(stream, "%8ld %12.2f\n",
                    YMDDateFromIRDate(tree_data->TPDate[i]),
                    ladder_data->FixingFl * 100);
        }
    }

    /* Date adjustment offset */
    fprintf (stream, "%s Date adjustment offset \n%c\n",
                     LegName, 
                     'U');
    
    /* Number of dates and levels */
    fprintf (stream, "%s Nb of dates and levels \n%d \n", 
                     LegName, 
                     size);
    
    fprintf (stream, "%s res, eff, acs, ace, pay - "
                     "dcf, notional, floor, cap, spread\n", 
             LegName);
    for (i = 0; i <= tree_data->NbTP; i++)
    {
        if (tree_data->TPtype[4][i])
        {
            fprintf(stream, "%8ld %8ld %8ld %8ld %8ld "
                        "%12.8f %12.6f %12.6f %12.6f %8.6f\n",
                        YMDDateFromIRDate(tree_data->TPDate[i]), // res
                        YMDDateFromIRDate(tree_data->TPDate[i]), // eff
                        YMDDateFromIRDate(tree_data->CritDate[4][i].SuppDate[1]), // acs
                        YMDDateFromIRDate(tree_data->CritDate[4][i].SuppDate[0]), // ace
                        YMDDateFromIRDate(tree_data->CritDate[4][i].SuppDate[0]), // pay
                        tree_data->CritDate[4][i].Value[4], // dcf
                        tree_data->CritDate[4][i].Value[3] * ladder_data->NotionalSign * ladder_data->Notional, // not
                        tree_data->CritDate[4][i].Value[0] * 100, // floor
                        tree_data->CritDate[4][i].Value[1] * 100, // cap
                        tree_data->CritDate[4][i].Value[2] * 100);// spread
        }
    }

    /* Funding index */
    fprintf (stream, "%s Index maturity, frequency, DCC, curve, weight\n%d %c %c %d %lf\n", 
                     LegName, 
                     ladder_data->IdxMatFl,
                     ladder_data->IdxFreqFl,
                     ladder_data->IdxBaseFl,
                     ladder_data->IdxIoDFl,
                     ladder_data->IdxWeightFl);


    /*--- Model ---*/
    LegName = "# Model:";

    /* First index to be used for volatility calibration */
    if (!mktvol_data->CalibFlag)
    {
        strcpy(buff, "nil");
    }
    else
    {
        int idxMat = Months360(mktvol_data->SwapSt[0], mktvol_data->SwapMat[0]);

        /* Fixed maturity date */
        if (ABS(Daysact(mktvol_data->SwapMat[0], 
                        mktvol_data->SwapMat[mktvol_data->NbVol-1])) < 87) 
        {
            sprintf(buff, "%dyFix", idxMat / 12);
        }
        else if (idxMat > 12)  /* CMS */
        {
            sprintf(buff, "%dyCms", idxMat / 12);
        }
        else   /* base vol */
        {
            sprintf(buff, "%dm", idxMat);
        }
    }

    /* Calibration index 1 */
    fprintf (stream, "%s Calibration index 1\n%s \n", LegName, buff);

    /* Calibration index 2 */
    fprintf (stream, "%s Calibration index 2\n%s \n", LegName, buff);

    /* Discount curve */
    fprintf (stream, "%s Discount curve\n%d \n", 
                     LegName,
                     tree_data->CvDisc);

    /* Number of st dev for tree */
    fprintf (stream, "%s Number of St Dev for tree \n%-3d \n", 
                     LegName,
                     tree_data->NbSigmaMax);

    /* Number of state variables */
    fprintf (stream, "%s Number of state variables \n%-3d \n", 
                     LegName,
                     ladder_data->NbStates);

    /* Number of st dev for state variable */
    fprintf (stream, "%s Number of St Dev for state variables \n%-3d \n", 
                     LegName,
                     ladder_data->NbStDev);

    /* Number of st dev for state variable */
    fprintf (stream, "%s PPY override\n%-3d \n", 
                     LegName,
                     tree_data->Ppy);

    /* Smile    */
    fprintf(stream, "%s Smile parameters (qL, qR, FwdShift, CetIter)\n%lf %lf %lf %d\n", 
            LegName,
            1.0 - mktvol_data->QLeft,
            1.0 - mktvol_data->QRight,
            mktvol_data->FwdShift,
            mktvol_data->CetNbIter);

    /* Number of factors */
    fprintf(stream, "%s Number of factors\n%d\n", 
            LegName,
            tree_data->NbFactor);

    /* Factor weight */
    fprintf(stream, "%s Factor weights (one per factor)\n", LegName);
    for (i=0; i<tree_data->NbFactor; i++)
        fprintf(stream, "%lf ", mktvol_data->Alpha[i]);
    fprintf(stream, "\n");

    /* Mean reversion */
    fprintf(stream, "%s Mean reversion coefficient (one per factor)\n", LegName);
    for (i=0; i<tree_data->NbFactor; i++)
        fprintf(stream, "%lf ", mktvol_data->Beta[i]);
    fprintf(stream, "\n");

    /* Correlations */
    fprintf(stream, "%s correlation (0, 1 or 3 numbers)\n", LegName);
    for (i=0; i<(tree_data->NbFactor)*(tree_data->NbFactor-1)/2; i++)
        fprintf(stream, "%lf ", mktvol_data->Rho[i]);
    fprintf(stream, "\n");

    /* Backbone */
    fprintf(stream, "%s Backbone\n%lf\n", LegName, 1. - mktvol_data->Bbq);

    /* Stats and trace flags */
    fprintf(stream, "%s Stats and trace flags\n%c %c\n", LegName, 'N', 'N');

    fclose(stream);
    return SUCCESS;
}



/*****  Ladder_Data_Init  ********************************************/
/*
*       Initialize the data 
*/
void     Ladder_Data_Init(    
    LADDER_DATA    *ladder_data)    /* (I) Deal data       */
{

    ladder_data->RibObsDate        = NULL;
    ladder_data->RibObsEffDate     = NULL;
    ladder_data->RibLoBarrier      = NULL;
    ladder_data->RibHiBarrier      = NULL;
    ladder_data->RibInRangeWeight  = NULL;
    ladder_data->RibOutRangeWeight = NULL;
    
}


/*****  Ladder_Data_Free  ********************************************/
/*
*       Free the data 
*/
void     Ladder_Data_Free(    
    LADDER_DATA    *ladder_data)    /* (I) Deal data       */
{

    long    NbDates = ladder_data->NbRibObsDates;

    Free_DR_Array (ladder_data->RibObsDate,        LONG,   0, NbDates-1);
    Free_DR_Array (ladder_data->RibObsEffDate,     LONG,   0, NbDates-1);
    Free_DR_Array (ladder_data->RibLoBarrier,      DOUBLE, 0, NbDates-1);
    Free_DR_Array (ladder_data->RibHiBarrier,      DOUBLE, 0, NbDates-1);
    Free_DR_Array (ladder_data->RibInRangeWeight,  DOUBLE, 0, NbDates-1);
    Free_DR_Array (ladder_data->RibOutRangeWeight, DOUBLE, 0, NbDates-1);

    ladder_data->RibObsDate        = NULL;
    ladder_data->RibObsEffDate     = NULL;
    ladder_data->RibLoBarrier      = NULL;
    ladder_data->RibHiBarrier      = NULL;
    ladder_data->RibInRangeWeight  = NULL;
    ladder_data->RibOutRangeWeight = NULL;

}


