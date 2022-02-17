/****************************************************************************/
/*      Standard input output for yield and volatility curves.              */
/****************************************************************************/
/*      STDINPUT.c                                                          */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/stdinput.c,v 1.21 2005/02/04 20:01:48 skuzniar Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include "fix123head.h"



/*****	FindAndSkipComLine  *************************************************/
/*
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine (FILE*       stream,
                            char const* Label,
                            char const* Routine,
                            char const* FileName)
{

    char    ErrorMsg[MAXBUFF];
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */


    if (fgets (string, MAXBUFF, stream) == NULL)
    {
        sprintf (ErrorMsg, "%s: %s comment line expected and not found in file"
                    "%s!", Routine, Label, FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    else
    {
        if (string[0] != '#')
        {
            sprintf (ErrorMsg, "%s: expected '#' in %s comment line in file"
                        "%s!", Routine, Label, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }  /* if then else */

    
    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* FindAndSkipComLine */

/*****	FindAndSkipComLine1  *************************************************/
/*
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine1 (FILE*       stream,
                            char const* Label,
                            char const* Routine,
                            char const* FileName)
{

    char    ErrorMsg[MAXBUFF];
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */


    if (fgets (string, MAXBUFF, stream) == NULL)
    {
        sprintf (ErrorMsg, "%s: %s comment line expected and not found in file  "
                    "%s!", Routine, Label, FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    else
    {
        if (string[0] != '#')
        {
            sprintf (ErrorMsg, "%s: expected '#' in %s comment line in file"
                        "%s!", Routine, Label, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }  /* if then else */

    
    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* FindAndSkipComLine */



/*****	Term_Input_W  *******************************************************/
/*
*       Read term structure input for DR Wrapper and check validity of input.
*/
int  Term_Input_W (T_CURVE*    t_curve,  /* (O) Structure of zero curve data  */
                   char const* FileName) /* (I) File name including extension */
{

    long    LastDate;           /* Last date in zero curve */
    int     i;
    int     readerror;          /* Reading error status */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;


    /* Open the yield curve data file (see termodel.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "value date", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", &(t_curve->Today));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read value date in %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Today's date and spot days are not available */
    t_curve->ValueDate = t_curve->Today;
    t_curve->SpotDays  = 0;
    
    if (FindAndSkipComLine (stream, "money market basis", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(t_curve->MMB));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read money market basis in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }                                
        
    if (FindAndSkipComLine (stream, "yield curve frequency", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(t_curve->SwapFreq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read benchmark swap frequency in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "swap DCC", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%s \n", t_curve->SwapDCC);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read swap DCC in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "number of zeros", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(t_curve->NbZero));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of zeros in file %s! (Term_Input_W)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (t_curve->NbZero > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of rates in file %s exceeds maximum of %d! (Term_Input_W)", FileName, MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "zero dates and rates", "Term_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < t_curve->NbZero; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \n",
                            &(t_curve->ZeroDate[i]),
                            &(t_curve->Zero[i]));
                        
        t_curve->Zero[i] /= 100.;
                
        if (readerror != 2)
        {        
            sprintf (ErrorMsg, "Could not read zero date and rate #%d in file %s! (Term_Input_W)", i + 1, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }
        

    LastDate = Nxtmth (t_curve->ValueDate, 1200L, 1L);                          /* We require at least 100 years of zero curve for swaption vol bootstrapping */

        
    if (t_curve->ZeroDate[t_curve->NbZero-1] < LastDate)                        /* If the zero curve does not extend up to last date we add an extra point */
    {
        t_curve->NbZero += 1;                       
        t_curve->ZeroDate[t_curve->NbZero-1] = LastDate;
        t_curve->Zero[t_curve->NbZero-1] = t_curve->Zero[t_curve->NbZero-2];    /* Flat zero curve */
    }
    
                
    if (Term_Check_W (t_curve) == FAILURE)                    
    {        
        goto RETURN;
    }
        

    status = SUCCESS;
        
    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

}  /* Term_Input_W */



/*****  Term_Check_W  *******************************************************/
/*
*       Check validity of DR Wrapper term structure inputs.
*/
int     Term_Check_W (T_CURVE  *t_curve) /* (I) Structure of zero curve data */
{

    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */

    char    ErrorMsg[MAXBUFF];
    

    if (Dateok(t_curve->ValueDate))
    {
        DR_Error("Incorrect format for value date!");
        goto RETURN;
    }

    if ((t_curve->MMB != 360) && (t_curve->MMB != 365))
    {
        DR_Error("Money Market basis should be 360 or 365!");
        goto RETURN;        
    }

    if ((t_curve->SwapFreq != 'A') && (t_curve->SwapFreq != 'S') && 
		(t_curve->SwapFreq != 'Q') && (t_curve->SwapFreq != 'M'))
    
	{
        DR_Error("Incorrect benchmark swap frequency (please use A, S, Q or M)!");
        goto RETURN;
    }

    if (       strcmp(t_curve->SwapDCC,"ACT") 
            && strcmp(t_curve->SwapDCC,"365") 
            && strcmp(t_curve->SwapDCC,"360"))
    {
        DR_Error("Specify 'ACT' for Actual, '365' for 365 Fixed or '360' for swap DCC!");
        goto RETURN;
    }

    if (t_curve->NbZero > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of zero rates exceeds maximum of %d!", MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    for (i = 1; i < t_curve->NbZero; i++)
        if (t_curve->ZeroDate[i] <= t_curve->ZeroDate[i-1])
        {
            DR_Error("Zero dates must be entered in ascending order!");
            goto RETURN;
        }

    for (i = 0; i < t_curve->NbZero; i++)
    {
        if ((t_curve->Zero[i] < -100.) || (t_curve->Zero[i] > 100.))
        {
            DR_Error("Zero rates out of range!");
            goto RETURN;
        }
    }


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Term_Check_W */



/*****	BaseVol_Input_W  ****************************************************/
/*
*       Read base volatility input and check validity of input.
*/
int     BaseVol_Input_W (
            MKTVOL_DATA  *mktvol_data,  /* (O) Volatility data               */
            char const*   FileName)     /* (I) File name including extension */
{

    int     i, j;
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;


    /* Open the base volatility curve data file (see termodel.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (BaseVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "base vol frequency", "BaseVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(mktvol_data->Freq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read base vol frequency in file %s! (BaseVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "number of base vols", "BaseVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(mktvol_data->NbVol));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of base vols in file %s! (BaseVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (mktvol_data->NbVol > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of base vols in file %s exceeds maximum of %d! (BaseVol_Input)", FileName, MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "base vol dates and rates", "BaseVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf\n", 
                            &(mktvol_data->VolDate[i]),
                            &(mktvol_data->Vol[i]));

        if (readerror != 2)
        {        
            sprintf (ErrorMsg, "Could not read base vol date and rate #%d in file %s! (BaseVol_Input)", i+1, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        mktvol_data->Vol[i] /= 100.;
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


    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }
        
    return (status);

}  /* BaseVol_Input_W */



/*****	SwapVol_Input_W  ****************************************************/
/*
*       Read swaption volatility input and check validity of input.
*/
int     SwapVol_Input_W (int          *NbRows,      /* (O) Volatility matrix */
                         int          *NbCol,
                         long         **Expiry,
                         long         **FwdMat,
                         double      ***VolMatrix,
                         char const    *FileName)     /* (I) File name */
{
    int     i, j;
    int     readerror;          /* Reading error status         */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;


    /* Open the swaption volatility curve data file (see fix123.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (SwapVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;        
    }
        
    if (FindAndSkipComLine (stream, "number of rows", "SwapVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbRows);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of rows in file %s! (SwapVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (*NbRows > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of expiry in file %s exceeds maximum of %d! (SwapVol_Input)", FileName, MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "number of columns", "SwapVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbCol);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of columns in file %s! (SwapVol_Input)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "swaption matrix", "SwapVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    *Expiry    = (long *)    DR_Array  (LONG,    0, *NbRows);
    *FwdMat    = (long *)    DR_Array  (LONG,    0, *NbCol);
    *VolMatrix = (double **) DR_Matrix (DOUBLE, 0, *NbRows, 0, *NbCol);

    if    ((*Expiry    == NULL)
        || (*FwdMat    == NULL)
        || (*VolMatrix == NULL))
    {
            DR_Error("Could not allocate memory (SwapVol_Input)!");
            goto RETURN;
    }

    for (j = 0; j < *NbCol; j++)
    {
        readerror = fscanf (stream, "\t%ld", &((*FwdMat)[j]));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read forward maturity in file %s! (SwapVol_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        (*FwdMat)[j] *= 12;    /* Convert to months */
    }
        
    for (i = 0; i < *NbRows; i++)
    {
        readerror = fscanf (stream, "%ld", &((*Expiry)[i]));
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read expiry in file %s! (SwapVol_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
            
        for (j = 0; j < *NbCol; j++)
        {
            readerror = fscanf (stream, "\t%lf", &((*VolMatrix)[i][j]));
            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read volatility in file %s! (SwapVol_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            (*VolMatrix)[i][j] /= 100.;
        }
    }  /* for i */  
        

    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

}  /* SwapVol_Input_W */

/*****  MktVol_Input_W  *********************************************************/
/*
*       Utility routine converting an index name to a series of option expi-
*       ries, index maturities and volatilities taken from either the base
*       volatility curve or the swaption matrix.
*/
int     MktVol_Input_W (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                        char            *Index,       /* (I) Index to calibrate   */
                        char            *Index2,      /* (I) Other index          */
                        T_CURVE         *t_curve,     /* (I) Term structure data  */
                        char const*      BaseVolFile, /* (I) Base vol curve file  */
                        char const*      SwapVolFile) /* (I) Swaption matrix file */
{

    char    *StrIdx = NULL, *StrIdx2 = NULL;     
    long    IdxMat;             /* Maturity of the index (final or forward) */
    long    IdxMat2;            /* Maturity of the index (final or forward) */
    long    Mat;
    int     NbRows=0;           /* Swaption volatility matrix */
    int     NbCol=0;
    long    *Expiry = NULL;
    long    *FwdMat = NULL;
    double  **VolMatrix = NULL;
    int     i, j, j2;
    int     status = FAILURE;   /* Error status = FAILURE initially         */


    /* No calibration case */
    if (strstr (Index, "nil") != NULL)
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


    /* Search for * character in calibr ation index name */
    StrIdx = strchr(Index, '*');

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
    StrIdx = strstr (Index, "yCms");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index) * 12;


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


        /* Read full swaption matrix */
        if (SwapVol_Input_W (   &NbRows,
                                &NbCol,
                                &Expiry,
                                &FwdMat,
                                &VolMatrix,
                                SwapVolFile) == FAILURE)
        {
            goto RETURN;
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

        if (strcmp (Index, Index2))
        {
            StrIdx2 = strstr (Index2, "yCms");
            if (StrIdx2 == NULL)
            {
                DR_Error ("Calibration indices can be different only "
                          "if both are of XyCms type");
                goto RETURN;
            }
            *StrIdx2 = '\0';
            IdxMat2 = atoi(Index2) * 12;

            /* Check other Cms column */
            j2 = 0; 
            while ((j2 < NbCol-1) && (IdxMat2 > FwdMat[j2]))
                j2++;

            if (IdxMat2 != FwdMat[j2])
            {        
                DR_Error("Cms index2 is not in swaption matrix (MktVol_Input_W)!");
                goto RETURN;
            }
        }

        mktvol_data->NbVol = NbRows;

        for (i = 0; i < mktvol_data->NbVol; i++)
        {                                                                       
            /* 
             *  VolDate and SwapSt are identical: i.e. we don't calibrate 
             *  options on forward starting swaps (e.g. mid curve options) 
             */
            mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], IdxMat, 1L);
            mktvol_data->Vol[i]     = VolMatrix[i][j];
            mktvol_data->VolUsed[i] = TRUE;

            /* if (strcmp (Index, Index2))
            { */
                mktvol_data->SwapMat2[i] = 
                    Nxtmth (mktvol_data->VolDate[i], IdxMat2, 1L);
            /* } */
        }


        if (MktVol_Check_W (mktvol_data)  == FAILURE)
        {
            goto RETURN;
        }

        goto DONE;
    }


    /* Final maturity indices */
    StrIdx = strstr (Index, "yFix");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index) * 12;


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


        /* Read full swaption matrix */
        if (SwapVol_Input_W (   &NbRows,
                                &NbCol,
                                &Expiry,
                                &FwdMat,
                                &VolMatrix,
                                SwapVolFile) == FAILURE)
        {
            goto RETURN;
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
            }
            else                            
            {
                Mat = FwdMat[j-1];
                mktvol_data->Vol[i] = VolMatrix[i][j-1];
            }

            mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], Mat, 1L);
            mktvol_data->VolUsed[i] = TRUE;

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
    StrIdx = strchr(Index, 'm');

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(Index);


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


        /* Read vol curve */
        if (BaseVol_Input_W (   mktvol_data,
                                BaseVolFile) == FAILURE)
        {
            goto RETURN;
        }


        if (12 / Conv_Freq (mktvol_data->Freq) != IdxMat)
        {
            DR_Error("Base vol curve frequency different from calibration "
                        "index (MktVol_Input_W)!");
            goto RETURN;
        }


        for (i = 0; i < mktvol_data->NbVol; i++)
        {
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], IdxMat, 1L);
            mktvol_data->VolUsed[i] = TRUE;
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

    Free_DR_Array  (Expiry    , LONG,   0, NbRows);
    Free_DR_Array  (FwdMat    , LONG,   0, NbCol);
    Free_DR_Matrix (VolMatrix , DOUBLE, 0, NbRows, 0, NbCol);

    return (status);

}  /* MktVol_Input_W */


/*****	MktVol_Check_W  *****************************************************/
/*
*       Check validity of base volatility inputs.
*/
int     MktVol_Check_W (MKTVOL_DATA     *mktvol_data) /* (O) Volatility data */
{

    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */

    char    ErrorMsg[MAXBUFF];
    

    if ( (mktvol_data->FilterSpotVolFlag != TRUE) &&
         (mktvol_data->FilterSpotVolFlag != FALSE) )
    {
        DR_Error("FilterSpotVol flag must be TRUE or FALSE (MktVol_Check_W))");
        goto RETURN;
    }

    if ( (mktvol_data->SmoothingFlag != 'Y') &&
         (mktvol_data->SmoothingFlag != 'N') )
    {
        DR_Error("Cet smoothing flag must be 'Y' or 'N' (MktVol_Check_W))");
        goto RETURN;
    }

    if ( (mktvol_data->FilterSpotVolFlag != TRUE) &&
         (mktvol_data->FilterSpotVolFlag != FALSE) )
    {
        DR_Error("FilterSpotVol flag must be TRUE or FALSE (MktVol_Check_W))");
        goto RETURN;
    }

    if ( (mktvol_data->SkipFlag != TRUE) &&
         (mktvol_data->SkipFlag != FALSE) )
    {
        DR_Error("SkipFlag must be TRUE or FALSE (MktVol_Check_W))");
        goto RETURN;
    }

    if ( (mktvol_data->CalibFlag != TRUE) &&
         (mktvol_data->CalibFlag != FALSE) )
    {
        DR_Error("CalibFlag must be TRUE or FALSE (MktVol_Check_W))");
        goto RETURN;
    }

    /* Nothing else to check */
    if (mktvol_data->CalibFlag == FALSE)
    {
        return (SUCCESS);
    }

    if (mktvol_data->NbVol > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of vols exceeds maximum of %d! (MktVol_Check_W)", MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
        if (Dateok(mktvol_data->VolDate[i]))
        {
            DR_Error("Incorrect format for volatility date (MktVol_Check_W)!");
            goto RETURN;
        }

    for (i = 1; i < mktvol_data->NbVol; i++)
        if (mktvol_data->VolDate[i] <= mktvol_data->VolDate[i-1])
        {
            DR_Error("Volatilty dates must be entered in ascending order "
                        "(MktVol_Check_W)!");
            goto RETURN;
        }

    if (Dateok(mktvol_data->BaseDate))
    {
        DR_Error("Incorrect format for base date (MktVol_Check_W)!");
        goto RETURN;
    }

    if (mktvol_data->BaseDate >= mktvol_data->VolDate[0])
    {
        DR_Error("Base date must be before vol date (MktVol_Check_W)!");
        goto RETURN;
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
        if ((mktvol_data->Vol[i] < .0001) || (mktvol_data->Vol[i] > 999.))
        {
            DR_Error("Volatilities out of range (MktVol_Check_W)!");
            goto RETURN;
        }

    for (i = 0; i < mktvol_data->NbVol; i++)
        if (mktvol_data->VolDate[i] > mktvol_data->SwapSt[i])
        {
            DR_Error("Benchmark swap starts before expiry date (MktVol_Check_W)!");
            goto RETURN;
        }

    for (i = 0; i < mktvol_data->NbVol; i++)
        if (mktvol_data->SwapSt[i] >= mktvol_data->SwapMat[i])
        {
            DR_Error("Benchmark swap matures before swap start (MktVol_Check_W)!");
            goto RETURN;
        }

    if (   (mktvol_data->Freq != 'A')
        && (mktvol_data->Freq != 'S')
        && (mktvol_data->Freq != 'Q')
        && (mktvol_data->Freq != 'M'))
    {
        DR_Error("Incorrect vol frequency (MktVol_Check_W)!");
        goto RETURN;
    }

    if ((   mktvol_data->DCC != '3') 
        && (mktvol_data->DCC != '5') 
        && (mktvol_data->DCC != '0'))
    {
        DR_Error("Incorrect day count convention (MktVol_Check_W)!");
        goto RETURN;
    }

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* MktVol_Check_W */



/*****  Param_Input  ********************************************************/
/*
*  	Read model parameters and check validity of input.
*/
int     Param_Input (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             TREE_DATA     *tree_data,                  /* (O) Tree data structure           */
             int           NbFactor,                    /* (I) Number of factors             */
             char          OverWriteString[6][MAXBUFF], /* (I) Overwrite strings             */
             char const*   FileName)                    /* (I) File name including extension */
{

    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    *getserror;         /* Reading error for fgets          */
    char    string[MAXBUFF];
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    stream = fopen (FileName, "r");

    /*
     *  If there is no parameter file use the overwrite strings.
     */

    if (stream == NULL)
    {
        readerror = sscanf (OverWriteString[0],
                            "%d \n", 
                            &(tree_data->Ppy));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not find file %s: Ppy overwrite required! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (strstr(OverWriteString[1], "N") != NULL)
        {
            mktvol_data->QLeft    = 0.;
            mktvol_data->QRight   = 0.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }
        else if (strstr(OverWriteString[1], "L") != NULL)
        {
            mktvol_data->QLeft    = 1.;
            mktvol_data->QRight   = 1.;
            mktvol_data->FwdShift = 0.;
            mktvol_data->CetNbIter = 0;
        }
        else
        {
            readerror = sscanf (OverWriteString[1],
                                "%lf %lf %lf %d\n", 
                                &(mktvol_data->QLeft),
                                &(mktvol_data->QRight),
                                &(mktvol_data->FwdShift),
                                &(mktvol_data->CetNbIter));

            if (readerror != 4)
            {      
                sprintf (ErrorMsg, "Could not find file %s: Q overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* At input level q=0 means log-normal whereas */
            /* internally q=0 means normal: we switch here */
            mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
            mktvol_data->QRight = 1. - mktvol_data->QRight;
        }

        if (NbFactor == 1)
        {
            readerror = sscanf (OverWriteString[2], 
                                "%lf \n", 
                                &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* Fill in the unused parameters with N/A values */

            mktvol_data->Alpha[1] = -999.;
            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[1]  = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[0]   = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;
        }
        else if (NbFactor == 2)
        {
            readerror = sscanf (OverWriteString[2], 
                                "%lf \t%lf \n", 
                                &(mktvol_data->Alpha[0]),
                                &(mktvol_data->Alpha[1]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;
        }
        else if (NbFactor == 3)
        {
            readerror = sscanf (OverWriteString[2], 
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Alpha[0]),
                                &(mktvol_data->Alpha[1]),
                                &(mktvol_data->Alpha[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]),
                                &(mktvol_data->Beta[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Rho[0]),
                                &(mktvol_data->Rho[1]),
                                &(mktvol_data->Rho[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }  /* if NbFactor == */

        readerror = sscanf (OverWriteString[5],
                            "%lf \n", 
                            &(mktvol_data->Bbq));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find file %s: Bbq overwrite required! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* At input level q=0 means log-normal whereas */
        /* internally q=0 means normal: we switch here */
        mktvol_data->Bbq  = 1. - mktvol_data->Bbq;
    
    }
    else
    {        
        /* 
         *  Parameter file exists: use the overwrite strings if they are not 
         *  "nil", otherwise use the values in the parameter file.
         */
        if (NbFactor == 1)
        {
            /*
             *  Read mean reversion in parameter file.
             */

            if (FindAndSkipComLine (stream, "one factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* 
             *  If overwrite exists, use it.
             */

            if (strstr (OverWriteString[3], "nil") == NULL)                     /* Need to use strstr() here, strcmp won't do */
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \n", 
                                    &(mktvol_data->Beta[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \n", 
                                    &(mktvol_data->Alpha[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor ppy", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%d \n", 
                                &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            /* Skip two and three factor parameters in the file */
            for (i = 0; i < 12; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find two factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find three factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }


            mktvol_data->Alpha[1] = -999.;
            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[1]  = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[0]   = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;

        }
        else if (NbFactor == 2)
        {
            /* Skip one factor parameters in the file */
            for (i = 0; i < 6; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find one factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]));

                if (readerror != 2)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]));

                if (readerror != 2)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor correlation", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \n", 
                                    &(mktvol_data->Rho[0]));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for correlation! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor ppy", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find three factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;

        }
        else if (NbFactor == 3)
        {
            /* Skip one and two factor parameters in the file */
            for (i = 0; i < 6; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find one factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 12; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find two factor parameters in file %s! (Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]),
                                    &(mktvol_data->Beta[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]),
                                    &(mktvol_data->Alpha[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Rho[0]),
                                    &(mktvol_data->Rho[1]),
                                    &(mktvol_data->Rho[2]));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for correlation! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor ppy", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

        }  /* if NbFactor == */

        if (FindAndSkipComLine (stream, "QLeft", "Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->QLeft));
        
        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find QLeft in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "QRight", "Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->QRight));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find QRight in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "FwdShift", "Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->FwdShift));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find FwdShift in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "CetNbIter", "Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%d \n", 
                            &(mktvol_data->CetNbIter));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find CetNbIter in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
        mktvol_data->QRight = 1. - mktvol_data->QRight;

        if (strstr (OverWriteString[1], "nil") == NULL)
        {
            if (strstr(OverWriteString[1], "N") != NULL)
            {
                mktvol_data->QLeft    = 0.;
                mktvol_data->QRight   = 0.;
                mktvol_data->FwdShift = 0.;
                mktvol_data->CetNbIter = 0;
            }       
            else if (strstr(OverWriteString[1], "L") != NULL)
            {
                mktvol_data->QLeft    = 1.;
                mktvol_data->QRight   = 1.;
                mktvol_data->FwdShift = 0.;
                mktvol_data->CetNbIter = 0;
            }
            else
            {

                readerror = sscanf (OverWriteString[1],
                                    "%lf %lf %lf %d\n", 
                                    &(mktvol_data->QLeft),
                                    &(mktvol_data->QRight),
                                    &(mktvol_data->FwdShift),
                                    &(mktvol_data->CetNbIter));

                if (readerror != 4)
                {      
                    sprintf (ErrorMsg, "Could not read overwrite string for Qweight! (Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }       
                
                mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                mktvol_data->QRight = 1. - mktvol_data->QRight;
            }
        }  /* if then else */
        
        if (strstr (OverWriteString[5], "nil") != NULL)
        {
            if (FindAndSkipComLine (stream, "Bbq", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                "%lf \n", 
                                &(mktvol_data->Bbq));

            if (readerror != 1)
            {      
                sprintf (ErrorMsg, "Could not find Bbq in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            mktvol_data->Bbq = 1. - mktvol_data->Bbq;

        }
        else
        {
            readerror = sscanf (OverWriteString[5],
                                "%lf \n", 
                                &(mktvol_data->Bbq));

            if (readerror != 1)
            {      
                sprintf (ErrorMsg, "Could not read overwrite string for Bbq! (Param_Input)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }       
            
            mktvol_data->Bbq  = 1. - mktvol_data->Bbq;

        }  /* if then else */

    }  /* if then else */

    /* Check validity of input */
    if (Param_Check (   NbFactor,
                        mktvol_data,
                        tree_data) == FAILURE)              
    {        
        goto RETURN;
    }
        
    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

}  /* Param_Input */



/*****  Param_Check  ********************************************************/
/*
*  	Read term structure input and check validity of input.
*/
int     Param_Check (  	int             NbFactor,                               /* (I) Number of factors                     */
                        MKTVOL_DATA     *mktvol_data,                           /* (I) Structure of swaption volatility data */
                        TREE_DATA       *tree_data)                             /* (I) Tree data structure                   */
{
    int
        status = FAILURE;                                                       /* Error status = FAILURE initially */
       
    char
        ErrorMsg[MAXBUFF];

    double  norm;
    int     i;

    /* Nb of factors */
    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    /* Set total vol constants */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) norm += mktvol_data->Alpha[i] * mktvol_data->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
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

    /* Nb of st dev check */
    if (tree_data->NbSigmaMax < 3)
    {
        DR_Error("Can't cut the tree at less than three std devs!");
        goto RETURN;
    }              

    /* Ppy check */
    if ((tree_data->Ppy < 1) || (tree_data->Ppy > 365))
    {
        DR_Error("Ppy has to be between 1 and 365!");
        goto RETURN;
    }

    if ((mktvol_data->QLeft  < -10.) || (mktvol_data->QLeft  > 10.) ||
        (mktvol_data->QRight < -10.) || (mktvol_data->QRight > 10.))
    {
        DR_Error("Q out of range!");
        goto RETURN;
    }

    if (IS_EQUAL(1. + mktvol_data->FwdShift,0))
    {
        DR_Error("Fwd shift is singular!");
        goto RETURN;
    }

    if (NbFactor == 1)
    {
        if (mktvol_data->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[0] < 0.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }
    }
    else if (NbFactor == 2)
    {
        if (mktvol_data->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if (mktvol_data->Alpha[1] < 0.0001)
        {
            DR_Error("Weight #2 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[0] < 0.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[1] < 0.) || (mktvol_data->Beta[1] > 10.))
        {
            DR_Error("Beta #2 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Rho[0] < -.95) || (mktvol_data->Rho[0] > .95))
        {
            DR_Error("Correlation out of range!");
            goto RETURN;
        }
    }
    else if (NbFactor == 3)
    {
        if (mktvol_data->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if (mktvol_data->Alpha[1] < 0.0001)
        {
            DR_Error("Weight #2 out of range!");
            goto RETURN;
        }

        if (mktvol_data->Alpha[2] < 0.0001)
        {
            DR_Error("Weight #3 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[0] < 0.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[1] < 0.) || (mktvol_data->Beta[1] > 10.))
        {
            DR_Error("Beta #2 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[2] < 0.) || (mktvol_data->Beta[2] > 10.))
        {
            DR_Error("Beta #3 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Rho[0] < -.95) || (mktvol_data->Rho[0] > .95))
        {
            DR_Error("Correlation #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Rho[1] < -.95) || (mktvol_data->Rho[1] > .95))
        {
            DR_Error("Correlation #2 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Rho[2] < -.95) || (mktvol_data->Rho[2] > .95))
        {
            DR_Error("Correlation #3 out of range!");
            goto RETURN;
        }
    }  /* if then else */

    /* Nb of iterations */
    if (mktvol_data->CetNbIter > MAX_ITERATIONS)
    {
        sprintf(ErrorMsg,"Maximum allowed number of iterations is %d!",
                MAX_ITERATIONS);

        DR_Error(ErrorMsg);
        goto RETURN;
    }   

    /* Bbq parameters */
    if ((mktvol_data->Bbq  < -10.) || (mktvol_data->Bbq  > 10.))
    {
        DR_Error("Bbq out of range!");
        goto RETURN;
    }

    if ((mktvol_data->VolNorm < 0.00001) && (mktvol_data->VolLogn < 0.0001))
    {
        DR_Error("Bbq volatility out of range!");
        goto RETURN;
    }


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Param_Check */

/*****  MR_Input_W  *********************************************************/
/*
* Read time-dependent mean-reversion input      
*/
int     MR_Input_W (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                    char            *FileName) 
{
    double MrVNFM;
    int NbMr;
    int CalibFlag; 
    int NbMrCet;
    long currDate, prevDate;
    double tmp;
    int t;

    FILE     *stream = NULL;
    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int status = FAILURE;

    stream = fopen (FileName, "r");

    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* VNFM mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "VNFM mean-reversion", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf\n", &MrVNFM);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->MrVNFM = MrVNFM;

    /* MR calibration flag */
    if (FindAndSkipComLine (stream, 
                            "MR calibration flag", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &CalibFlag);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->MrCalibFlag = CalibFlag;

    /* Max number of MR CET iterations */
    if (FindAndSkipComLine (stream, 
                            "Max number of MR CET iterations", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbMrCet);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMrCet < 0)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of MR CET in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
        
    mktvol_data->NbMrCet = NbMrCet;
    
    /* Number of mean-reversion dates */
    if (FindAndSkipComLine (stream, 
                            "Number of mean-reversion dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbMr);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMr < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbMr = NbMr;

    /* Term structure of mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "Term structure of mean-reversion", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevDate = -1;
    for ( t = 0; t < NbMr; t++)
    {
        readerror = fscanf (stream, "%ld ", &currDate);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read mr date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currDate <= prevDate)
        {
            sprintf(ErrorMsg, 
                "Invalid input for date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrDate[t] = currDate;
        prevDate = currDate;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read mean-reversion input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    /* Vol ratio dates and values */
    if (FindAndSkipComLine (stream, 
                            "Vol ratio dates and values", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevDate = -1;
    for ( t = 0; t < NbMr; t++)
    {
        readerror = fscanf (stream, "%ld ", &currDate);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol ratio date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currDate <= prevDate)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol ratio date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->RatioDate[t] = currDate;
        prevDate = currDate;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol ratio input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->VolRatioInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    status = SUCCESS;
        
    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }
        
    return (status);

} /* MR_Input_W */


/*****  MR_Input_W1  *********************************************************/
/*
* Read time-dependent mean-reversion input      
*/
int     MR_Input_W1 (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                    char            *FileName) 
{
    int NbVolShift;
    double MrShift;
    int NbMr;
    int CalibFlag; 
    int NbMrCet;
    long currExp, prevExp;
    double tmp;
    int t;

    FILE     *stream = NULL;
    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int status = FAILURE;

    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open file %s!", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Volatility shift for CET */
    /* Number of vol shift dates */
    if (FindAndSkipComLine (stream, 
                            "Number of vol shift dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbVolShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbVolShift < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbCetVolShift = NbVolShift;

    /* Term structure of vol shifts */
    if (FindAndSkipComLine (stream, 
                            "Term structure of vol shifts", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbVolShift; t++)
    {
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol shift expiry %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol shift expiry %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftExp[t] = currExp;
        mktvol_data->CetVolShiftDate[t] = Nxtmth (mktvol_data->BaseDate, currExp, 1L);
	
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol shift input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftVal[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */


    /* Mean-reversion shift */
    if (FindAndSkipComLine (stream, 
                            "Mean-reversion shift", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf\n", &MrShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read mean-reversion shift in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->CetMrShift = MrShift;

    /* MR calibration flag */
    if (FindAndSkipComLine (stream, 
                            "MR calibration flag", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &CalibFlag);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->MrCalibFlag = CalibFlag;

    /* Number of MR CET iterations */
    if (FindAndSkipComLine (stream, 
                            "Number of MR CET iterations", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbMrCet);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMrCet < 0)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of MR CET in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
        
    mktvol_data->NbMrCet = NbMrCet;
    
    /* Number of mean-reversion dates */
    if (FindAndSkipComLine (stream, 
                            "Number of mean-reversion dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbMr);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMr < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbMr = NbMr;

    /* Term structure of mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "Term structure of mean-reversion", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbMr; t++)
    {
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read mr date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrExp[t] = currExp;
        mktvol_data->MrDate[t] = Nxtmth (mktvol_data->BaseDate, currExp, 1L);
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read mean-reversion input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    /* Vol ratio dates and values */
    if (FindAndSkipComLine (stream, 
                            "Vol ratio dates and values", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbMr; t++)
    {
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol ratio date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol ratio date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->RatioExp[t] = currExp;
        mktvol_data->RatioDate[t] = Nxtmth (mktvol_data->BaseDate, currExp, 1L);
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol ratio input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->VolRatioInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    status = SUCCESS;
        
    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }
        
    return (status);

} /* MR_Input_W1 */

/*****  MR_Input_W2  *********************************************************/
/*
* Read time-dependent mean-reversion input      
*/
int     MR_Input_W2 (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                    char            *FileName) 
{
    int NbVolShift;
    double MrShift;
    int NbMr;
    int CalibFlag; 
    int NbMrCet;
    long currExp, prevExp;
    double tmp;
    int t;

    FILE     *stream = NULL;
    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int status = FAILURE;

	FILE     *stream1 = NULL;
	FILE     *stream2 = NULL;
	char     myDataString[MAXBUFF];

    stream = fopen (FileName, "r");
	stream1 = fopen ("C:/koswapm_t/mr_copy.dat", "w");
	stream2 = fopen (FileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open file %s!", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

	if (stream2 == NULL)
    {
        sprintf (ErrorMsg,
                 "stream2: Could not open file %s!", 
                 FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

	if (stream1 == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open file mr_copy.dat!");
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* Volatility shift for CET */
    /* Number of vol shift dates */
    if (FindAndSkipComLine (stream, 
                            "Number of vol shift dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbVolShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbVolShift < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbCetVolShift = NbVolShift;

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* Term structure of vol shifts */
    if (FindAndSkipComLine (stream, 
                            "Term structure of vol shifts", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbVolShift; t++)
    {
		fgets (myDataString, MAXBUFF, stream2);
	    fprintf (stream1, "%s", myDataString);
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol shift expiry %ld in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol shift expiry %ld in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftExp[t] = currExp;
	
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol shift input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftVal[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */


	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* Mean-reversion shift */
    if (FindAndSkipComLine (stream, 
                            "Mean-reversion shift", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%lf\n", &MrShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read mean-reversion shift in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->CetMrShift = MrShift;

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* MR calibration flag */
    if (FindAndSkipComLine (stream, 
                            "MR calibration flag", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &CalibFlag);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->MrCalibFlag = CalibFlag;

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* Number of MR CET iterations */
    if (FindAndSkipComLine (stream, 
                            "Number of MR CET iterations", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbMrCet);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMrCet < 0)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of MR CET in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
        
    mktvol_data->NbMrCet = NbMrCet;

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);
    
    /* Number of mean-reversion dates */
    if (FindAndSkipComLine (stream, 
                            "Number of mean-reversion dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbMr);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMr < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbMr = NbMr;

	fgets (myDataString, MAXBUFF, stream2);
	fprintf (stream1, "%s", myDataString);

    /* Term structure of mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "Term structure of mean-reversion", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbMr; t++)
    {
		fgets (myDataString, MAXBUFF, stream2);
	    fprintf (stream1, "%s", myDataString);

        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read mr date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrExp[t] = currExp;
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read mean-reversion input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    fprintf (stream1, "***");
	fscanf (stream2, "%s", myDataString);
	fprintf (stream1, "%s", myDataString);

    /* Vol ratio dates and values */
    if (FindAndSkipComLine (stream, 
                            "Vol ratio dates and values", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbMr; t++)
    {
		fgets (myDataString, MAXBUFF, stream2);
	    fprintf (stream1, "%s", myDataString);

        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol ratio date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol ratio date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->RatioExp[t] = currExp;
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol ratio input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->VolRatioInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    status = SUCCESS;
        
    RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }

	if (stream1 != NULL)
    {
        fclose (stream1);
    }

	if (stream2 != NULL)
    {
        fclose (stream1);
    }
        
    return (status);

} /* MR_Input_W2 */

/*****  MR_Input_W3  *********************************************************/
/*
* Read time-dependent mean-reversion input from file other than mr.dat      
*/
int     MR_Input_W3 (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                     FILE            *stream,
					 char            *FileName) 
{
    int NbVolShift;
    double MrShift;
    int NbMr;
    int CalibFlag; 
    int NbMrCet;
    long currExp, prevExp;
    double tmp;
    int t;

    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int status = FAILURE;

//	FILE     *stream1 = NULL;
//	FILE     *stream2 = NULL;
//	char     myDataString[MAXBUFF];

//	stream1 = fopen ("C:/koswapm_t/mr_copy.dat", "w");
//	stream2 = stream;

    if (stream == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open data file");
        DR_Error (ErrorMsg);
        goto RETURN;
    }

/* 	if (stream2 == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open data file");
        DR_Error (ErrorMsg);
        goto RETURN;
    }

	if (stream1 == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open file mr_copy.dat!");
        DR_Error (ErrorMsg);
        goto RETURN;
    } */

//  fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* Volatility shift for CET */
    /* Number of vol shift dates */
    if (FindAndSkipComLine (stream, 
                            "Number of vol shift dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//  fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbVolShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbVolShift < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of vol shift dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbCetVolShift = NbVolShift;

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* Term structure of vol shifts */
    if (FindAndSkipComLine (stream, 
                            "Term structure of vol shifts", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbVolShift; t++)
    {
//		fgets (myDataString, MAXBUFF, stream2);
//	    fprintf (stream1, "%s", myDataString);
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol shift expiry %ld in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol shift expiry %ld in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftExp[t] = currExp;
	
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol shift input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->CetVolShiftVal[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */


//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* Mean-reversion shift */
    if (FindAndSkipComLine (stream, 
                            "Mean-reversion shift", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%lf\n", &MrShift);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read mean-reversion shift in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->CetMrShift = MrShift;

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* MR calibration flag */
    if (FindAndSkipComLine (stream, 
                            "MR calibration flag", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &CalibFlag);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->MrCalibFlag = CalibFlag;

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* Number of MR CET iterations */
    if (FindAndSkipComLine (stream, 
                            "Number of MR CET iterations", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbMrCet);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read VNFM mean-reversion in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMrCet < 0)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of MR CET in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
        
    mktvol_data->NbMrCet = NbMrCet;

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);
    
    /* Number of mean-reversion dates */
    if (FindAndSkipComLine (stream, 
                            "Number of mean-reversion dates", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    readerror = fscanf (stream, "%d\n", &NbMr);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMr < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of mr dates in file %s! "
                "(MR_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbMr = NbMr;

//	fgets (myDataString, MAXBUFF, stream2);
//	fprintf (stream1, "%s", myDataString);

    /* Term structure of mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "Term structure of mean-reversion", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = 0;
    for ( t = 0; t < NbMr; t++)
    {
//		fgets (myDataString, MAXBUFF, stream2);
//	    fprintf (stream1, "%s", myDataString);

        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read mr date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrExp[t] = currExp;
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read mean-reversion input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

//    fprintf (stream1, "***");
//	fscanf (stream2, "%s", myDataString);
//	fprintf (stream1, "%s", myDataString);

    /* Vol ratio dates and values */
    if (FindAndSkipComLine (stream, 
                            "Vol ratio dates and values", 
                            "MR_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevExp = -1;
    for ( t = 0; t < NbMr; t++)
    {
//		fgets (myDataString, MAXBUFF, stream2);
//	    fprintf (stream1, "%s", myDataString);

        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read vol ratio date %d in file %s! "
                "(MR_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currExp <= prevExp)
        {
            sprintf(ErrorMsg, 
                "Invalid input for vol ratio date %d in file %s! "
                "(MR_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->RatioExp[t] = currExp;
        prevExp = currExp;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read vol ratio input at date %d in file %s! "
                    "(MR_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->VolRatioInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    status = SUCCESS;
        
RETURN:

/* 	if (stream1 != NULL)
    {
        fclose (stream1);
    } */
        
    return (status);

} /* MR_Input_W3 */

int MrExpToDates(int             BaseVolFlag,
		         MKTVOL_DATA    *mktvol_data)
{
    int i;
	int status = FAILURE;

	if (BaseVolFlag == TRUE)
    /* calibration to base vols, pass dates */
	{   
		if (mktvol_data->CetVolShiftExp[0] < 19000000)
		{
			DR_Error("Invalid dates for vol shift (MrExpToDates)");
			goto RETURN;
		}

		if (mktvol_data->MrExp[0] < 19000000)
		{
			DR_Error("Invalid dates for mean-reversion (MrExpToDates)");
			goto RETURN;
		}

		if (mktvol_data->RatioExp[0] < 19000000)
		{
			DR_Error("Invalid dates for vol ratios (MrExpToDates)");
			goto RETURN;
		}

		for (i = 0; i < mktvol_data->NbCetVolShift; i++)
		{
			mktvol_data->CetVolShiftDate[i] = mktvol_data->CetVolShiftExp[i];
		}

		for (i = 0; i < mktvol_data->NbMr; i++)
		{
			mktvol_data->MrDate[i] = mktvol_data->MrExp[i];
			mktvol_data->RatioDate[i] = mktvol_data->RatioExp[i];
		}
	}
	else
	/* calibration to swap vols, generate dates from exp in months */
	{
		for (i = 0; i < mktvol_data->NbCetVolShift; i++)
		{
			mktvol_data->CetVolShiftDate[i] = Nxtmth (mktvol_data->BaseDate, mktvol_data->CetVolShiftExp[i], 1L);
		}

		for (i = 0; i < mktvol_data->NbMr; i++)
		{
			mktvol_data->MrDate[i] = Nxtmth (mktvol_data->BaseDate, mktvol_data->MrExp[i], 1L);
			mktvol_data->RatioDate[i] = Nxtmth (mktvol_data->BaseDate, mktvol_data->RatioExp[i], 1L);
		}
	}

	status = SUCCESS;

RETURN:

	return(status);
}

/*****  MRTS_Input_W  *********************************************************/
/*
* Read mean-reversion term structure      
*/
int     MRTS_Input_W (MKTVOL_DATA     *mktvol_data, /* (O) Volatility data      */
                     FILE            *stream,
					 char            *FileName) 
{
    int NbMr;
    long currDate, prevDate;
    double tmp;
    int t;

    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int status = FAILURE;

    if (stream == NULL)
    {
        sprintf (ErrorMsg,
                 "Could not open data file");
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line", 
                            "MRTS_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Number of mean-reversion dates */
    if (FindAndSkipComLine (stream, 
                            "Number of mean-reversion dates", 
                            "MRTS_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbMr);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of mr dates in file %s! "
                "(MRTS_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    if (NbMr < 1)
    {
        sprintf(ErrorMsg, 
                "Invalid input for number of mr dates in file %s! "
                "(MRTS_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    mktvol_data->NbMr = NbMr;

    /* Term structure of mean-reversion */
    if (FindAndSkipComLine (stream, 
                            "Term structure of mean-reversion", 
                            "MRTS_Input_W", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    prevDate = 0;
    for ( t = 0; t < NbMr; t++)
    {
        readerror = fscanf (stream, "%ld ", &currDate);

        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                "Could not read mr date %d in file %s! "
                "(MRTS_Input_W)", 
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        if (currDate <= prevDate)
        {
            sprintf(ErrorMsg, 
                "Invalid input for date %d in file %s! "
                "(MRTS_Input_W)",
                t,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrDate[t] = currDate;
        prevDate = currDate;

        readerror = fscanf (stream, "%lf ", &tmp);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, 
                    "Could not read mean-reversion input at date %d in file %s! "
                    "(MRTS_Input_W)", 
                    t,
                    FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }

        mktvol_data->MrInput[t] = tmp;

        fscanf (stream, "\n");
    } /* for t */

    status = SUCCESS;
        
RETURN:
        
    return (status);

} /* MRTS_Input */

