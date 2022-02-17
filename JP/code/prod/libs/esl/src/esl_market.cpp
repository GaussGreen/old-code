#include "esl_market.h"
#include "esl_error.h"
#include "esl.h"

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>

/********************************************
 *
 *  Populate MKTVOL_DATA structure
 *
 *  Takes model parameters and overwrites
 *  from the deal file and populates
 *  MKTVOL_DATA, FX_DATA, TREE_DATA
 *  see also: MktVol_Input_W()
 *
 ********************************************/

int EslMktVolCalibration(
    MKTVOL_DATA*        mktvol_data,  /* (O) Volatility data */
    char const*         Index,        /* (I) Index to calibrate */
    T_CURVE const*      t_curve,      /* (I) Term structure data */
    BASEVOL_DATA const* baseVol_data, /* (I) Base vol */
    SWAPVOL_DATA const* swapVol_data) /* (I) Swaption matrix */
{
    return EslMktVolCalibrationAndRecording(
                mktvol_data,
                NULL,
                NULL,
                Index,
                t_curve,
                baseVol_data,
                swapVol_data);
}

/********************************************
*
*  Populate MKTVOL_DATA structure
*
*  Takes model parameters and overwrites
*  from the deal file and populates
*  MKTVOL_DATA, FX_DATA, TREE_DATA
*  see also: MktVol_Input_W()
*
*  Records dates of selected base or swaption
*  vols in selected_baseVols or selected_swapVols,
*  respectively, if these params are not NULL.
*
********************************************/
int EslMktVolCalibrationAndRecording(
    MKTVOL_DATA*        mktvol_data,  /* (O) Volatility data */
    BASEVOL_EXPOSURE_DATA* selected_baseVols, /* (O) Selected base vol info */
    SWAPVOL_EXPOSURE_DATA* selected_swapVols, /* (O) Selected swap vol info */
    char const*         Index,        /* (I) Index to calibrate */
    T_CURVE const*      t_curve,      /* (I) Term structure data */
    BASEVOL_DATA const* baseVol_data, /* (I) Base vol */
    SWAPVOL_DATA const* swapVol_data) /* (I) Swaption matrix */
{

    char     IndexL[MAXINDEX];   /* Local copy of index */
    char    *StrIdx = NULL;
    long     IdxMat;             /* Maturity of the index (final or forward) */
    long     Mat;
    int      NbRows = 0;         /* Swaption volatility matrix */
    int      NbCol = 0;
    long const* Expiry = NULL;
    long const* FwdMat = NULL;
    double   VolMatrix[NBSWAPTION][NBSWAPTION];
    int      i, j, tenorIndex;
    int      status = FAILURE;   /* Error status = FAILURE initially         */

    if (selected_baseVols != NULL)
        selected_baseVols->NbPoints = 0;

    if (selected_swapVols != NULL)
        selected_swapVols->NbPoints = 0;

    strcpy(IndexL,Index);

    /* No calibration case */
    if (strstr (IndexL, "nil") != NULL)
    {
        mktvol_data->CalibFlag = FALSE;

        /* Initialize unused variables */
        mktvol_data->BaseDate = 0;
        mktvol_data->NbVol = 0;
        mktvol_data->Freq = 'z';
        mktvol_data->DCC = 'z';
        mktvol_data->SkipFlag = FALSE;

        return (SUCCESS);
    }


    mktvol_data->CalibFlag = TRUE;


    /* Search for * character in calibration index name */
    StrIdx = strchr(IndexL, '*');

    if (StrIdx == NULL)
    {
        mktvol_data->SkipFlag = FALSE;  /* vol points skipping not allowed */
    }
    else
    {
        mktvol_data->SkipFlag = TRUE;   /* vol points skipping allowed      */
        *StrIdx = '\0';                 /* terminated index name before '*' */
    }


    /* Cms indices */
    StrIdx = strstr (IndexL, "yCms");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL) * 12;

        if (swapVol_data == NULL)
        {
            DR_Error ("swapVol_data not provided.");
            goto RETURN;
        }

        /* Read conventions from t_curve */
        mktvol_data->BaseDate = t_curve->Today;
        mktvol_data->Freq     = t_curve->SwapFreq;

        if (!strcmp(t_curve->SwapDCC, "360"))
        {
            mktvol_data->DCC = '0';
        }
        else if (!strcmp(t_curve->SwapDCC, "365"))
        {
            mktvol_data->DCC = '5';
        }
        else
        {
            mktvol_data->DCC = '3';
        }


        NbRows = swapVol_data->NbSwaptionExpiries;
        NbCol = swapVol_data->NbSwapTenors;
        Expiry = swapVol_data->SwaptionExpiries;
        FwdMat = swapVol_data->SwapTenors;

        for (i = 0; i < NbRows; i++)
        {
            for (j = 0; j < NbCol; j++)
            {
                VolMatrix[i][j] = swapVol_data->VolMatrix[i][j];
            }
        }

        /* Find required Cms column */
        j = 0;
        while ((j < NbCol-1) && (IdxMat > FwdMat[j]))
            j++;

        if (IdxMat != FwdMat[j])
        {
            DR_Error("Cms index is not in swaption matrix (MktVol_Input_W)!");
            goto RETURN;
        }


        mktvol_data->NbVol = NbRows;

        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            /*
             *  VolDate and SwapSt are identical: i.e. we don't calibrate
             *  options on forward starting swaps (e.g. mid curve options)
             */
            if (Expiry[i]==0) {
            	DR_Error("Expiry shorter than a month");
            	goto RETURN;
            }
            mktvol_data->VolDate[i] = Nxtmth(mktvol_data->BaseDate,Expiry[i],1L);
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i],IdxMat,1L);
            mktvol_data->Vol[i]     = VolMatrix[i][j];
            mktvol_data->VolUsed[i] = TRUE;

            if (selected_swapVols != NULL)
            {
                selected_swapVols->SwaptionExpiryDates[
                    selected_swapVols->NbPoints] = mktvol_data->VolDate[i];
                selected_swapVols->SwaptionExpiryIndices[
                    selected_swapVols->NbPoints] = i;
                selected_swapVols->SwapMatDates[
                    selected_swapVols->NbPoints] = mktvol_data->SwapMat[i];
                selected_swapVols->SwapTenorIndices[
                    selected_swapVols->NbPoints++] = j;
            }
        }


        if (MktVol_Check_W (mktvol_data)  == FAILURE)
        {
            goto RETURN;
        }

        goto DONE;
    }


    /* Final maturity indices */
    StrIdx = strstr (IndexL, "yFix");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL) * 12;

        if (swapVol_data == NULL)
        {
            DR_Error ("swapVol_data not provided.");
            goto RETURN;
        }

        /* Read conventions from t_curve */
        mktvol_data->BaseDate = t_curve->Today;
        mktvol_data->Freq     = t_curve->SwapFreq;

        if (!strcmp(t_curve->SwapDCC, "360"))
        {
            mktvol_data->DCC = '0';
        }
        else if (!strcmp(t_curve->SwapDCC, "365"))
        {
            mktvol_data->DCC = '5';
        }
        else
        {
            mktvol_data->DCC = '3';
        }


        NbRows = swapVol_data->NbSwaptionExpiries;
        NbCol = swapVol_data->NbSwapTenors;
        Expiry = swapVol_data->SwaptionExpiries;
        FwdMat = swapVol_data->SwapTenors;

        for (i = 0; i < NbRows; i++)
        {
            for (j = 0; j < NbCol; j++)
            {
                VolMatrix[i][j] = swapVol_data->VolMatrix[i][j];
            }
        }


        /* Process the final maturity index */
        for (i = 0; i < NbRows; i++)
        {
            Mat = IdxMat - Expiry[i];

            /* We got out of the swaption matrix */
            if (Mat < FwdMat[0])
                break;

            j = 0;
            while ((j < NbCol-1) && (Mat >= FwdMat[j]))
                j++;

            /* Use higher end of bracket: we don't interpolate to avoid stub */
            /* This includes the case FwdMat[i] > FwdMat[NbCol-1] so that    */
            /* we use a flat volatility after the last forward maturity.     */
            if (2 * Mat >= FwdMat[j-1] + FwdMat[j])
            {
                Mat = FwdMat[j];
                mktvol_data->Vol[i] = VolMatrix[i][j];
                tenorIndex = j;
            }
            else
            {
                Mat = FwdMat[j-1];
                mktvol_data->Vol[i] = VolMatrix[i][j-1];
                tenorIndex = j - 1;
            }

            mktvol_data->VolDate[i] = Nxtmth(mktvol_data->BaseDate,Expiry[i],1L);
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], Mat, 1L);
            mktvol_data->VolUsed[i] = TRUE;

            if (selected_swapVols != NULL)
            {
                selected_swapVols->SwaptionExpiryDates[
                    selected_swapVols->NbPoints] = mktvol_data->VolDate[i];
                selected_swapVols->SwaptionExpiryIndices[
                    selected_swapVols->NbPoints] = i;
                selected_swapVols->SwapMatDates[
                    selected_swapVols->NbPoints] = mktvol_data->SwapMat[i];
                selected_swapVols->SwapTenorIndices[
                    selected_swapVols->NbPoints++] = tenorIndex;
            }

        }  /* for i */

        if (i == 0)
        {
            DR_Error ("nyFix calibration falls outside swaption matrix "
                        "(MktVol_Input_W)!");
            goto RETURN;
        }

        mktvol_data->NbVol = i;


        if (MktVol_Check_W (mktvol_data)  == FAILURE)
        {
            goto RETURN;
        }

        goto DONE;
    }


    /* Base vol indices */
    StrIdx = strchr(IndexL, 'm');

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL);

        if (baseVol_data == NULL)
        {
            DR_Error ("baseVol_data not provided.");
            goto RETURN;
        }

        /* Read conventions from t_curve */
        mktvol_data->BaseDate = t_curve->Today;

        if (t_curve->MMB == 365)
        {
            mktvol_data->DCC = '5';
        }
        else
        {
            mktvol_data->DCC = '0';
        }

        /* Populate mktvol_data structure with base vols - this */
        /* replicates the work of BaseVol_Input_W(...) function */

        mktvol_data->Freq = baseVol_data->Frequency;
        mktvol_data->NbVol = baseVol_data->NbVols;
        if (mktvol_data->NbVol > MAXNBDATE)
        {
            DR_Error("Nb of base vols exceeds maximum of %d! "
                     "(BaseVol_Input)", MAXNBDATE);
        }

        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            mktvol_data->VolDate[i] = baseVol_data->VolDates[i];
            mktvol_data->Vol[i] = baseVol_data->Vols[i];
        }

        /* Eliminate dates falling before base date */
        j = 0;
        while (mktvol_data->VolDate[j] <= mktvol_data->BaseDate)
            j++;

        mktvol_data->NbVol -= j;

        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            mktvol_data->VolDate[i] = mktvol_data->VolDate[i+j];
            mktvol_data->Vol[i] = mktvol_data->Vol[i+j];
        }

        /* end of BaseVol_Input_W function logic */



        if (12 / Conv_Freq (mktvol_data->Freq) != IdxMat)
        {
            DR_Error("Base vol curve frequency different from calibration "
                        "index (MktVol_Input_W)!");
            goto RETURN;
        }


        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth(mktvol_data->VolDate[i],IdxMat,1L);
            mktvol_data->VolUsed[i] = TRUE;

            if (selected_baseVols != NULL)
            {
                selected_baseVols->VolDates[selected_baseVols->NbPoints] =
                    mktvol_data->VolDate[i];
                selected_baseVols->VolIndices[selected_baseVols->NbPoints++] = i;
            }
        }


        if (MktVol_Check_W (mktvol_data)  == FAILURE)
        {
            goto RETURN;
        }

        return (SUCCESS);
    }


    if (StrIdx == NULL)
    {
        DR_Error ("Incorrect calibration index (MktVol_Input_W)!");
        goto RETURN;
    }


    DONE:

    status = SUCCESS;

    RETURN:

    return (status);

}


/* function based on BaseVol_Input_W in stdinput.c */
static int EslReadBaseVolW(
    BASEVOL_DATA*   baseVol,
    char const*     FileName)     /* (I) File name including extension */
{

    int     i;
    int     readerror;          /* Reading error status             */
    const char *routine = "EslReadBaseVolW: ";
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;


    /* Open the base volatility curve data file (see termodel.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        DR_Error("%s Could not open file %s! ", 
                 routine,
                 FileName);
        
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "base vol frequency", routine, 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(baseVol->Frequency));
    if (readerror != 1)
    {        
        DR_Error("%s Could not read base vol frequency in file %s! ",
                 routine,
                 FileName);
        
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "number of base vols", routine, 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(baseVol->NbVols));
    if (readerror != 1)
    {        
        DR_Error("%s Could not read number of base vols in file %s! ",
                 routine,
                 FileName);
        
        goto RETURN;
    }
    
    if (baseVol->NbVols > MAXNBDATE)
    {        
        DR_Error("%s Nb of base vols in file %s exceeds maximum of %d! ",
                 routine,
                 FileName, 
                 MAXNBDATE);
        
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "base vol dates and rates", routine,
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < baseVol->NbVols; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf\n", 
                            &(baseVol->VolDates[i]),
                            &(baseVol->Vols[i]));

        baseVol->VolDates[i] = IRDateFromYMDDate(baseVol->VolDates[i]);
        baseVol->Vols[i] /= 100.0;

        if (readerror != 2)
        {        
            DR_Error("%s Could not read base vol date and rate #%d in file %s! ",
                     routine,
                     i+1, 
                     FileName);
            
            goto RETURN;
        }
    }

    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }
        
    return (status);

}

int EslReadVolsW(
    BASEVOL_DATA*   baseVol,
    SWAPVOL_DATA*   swapVol,
    char const*     baseVolFilename,
    char const*     swapVolFilename)
{
    int  status = FAILURE;
    const char *routine = "EslReadVolsW: ";

    int      i,j;
    int      NbRows = 0;         /* Swaption volatility matrix */
    int      NbCol = 0;
    long    *expiry = NULL;
    long    *fwdMat = NULL;
   	char    *DoMoY = NULL;
    double  **volMatrix = NULL;

    /* read swapvol matrix using existing esl functionality and
       copy results into swapVol structure */
    if (SwapVol_Input_W(&NbRows,
                        &NbCol,
                        &expiry,
                        &fwdMat,
                        &volMatrix,
                        &DoMoY,
                        swapVolFilename) == FAILURE)
    {
        goto RETURN;
    }

    if (NbRows > NBSWAPTION ||
        NbCol > NBSWAPTION)
    {
        DR_Error("%s Number of swaption vols %d > max size of internal structure %d",
                 routine, MAX(NbRows, NbCol), NBSWAPTION);
        goto RETURN;
    }

    swapVol->NbSwaptionExpiries = NbRows;
    swapVol->NbSwapTenors = NbCol;
    for (i = 0; i < NbRows; i++)
    {
        swapVol->SwaptionExpiries[i] = expiry[i];
        for (j = 0; j < NbCol; j++)
        {
            /* ??? hack - store vols as % when read from file - have to
                   convert them back as SwapVol_Input function /100 */
            swapVol->VolMatrix[i][j] = volMatrix[i][j];
        }
    }
    for (i = 0; i < NbCol; i++)
    {
        swapVol->SwapTenors[i] = fwdMat[i];
    }

    /* read base vols */
    if (EslReadBaseVolW(baseVol,
                         baseVolFilename) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    Free_DR_Array  (expiry    , LONG,   0, NbRows);
    Free_DR_Array  (fwdMat    , LONG,   0, NbCol);
    Free_DR_Array  (DoMoY     , CHAR,   0, NbRows);
    Free_DR_Matrix (volMatrix , DOUBLE, 0, NbRows, 0, NbCol);

    if (status == FAILURE)
        DR_Error("%s failed", routine);

    return status;
}
