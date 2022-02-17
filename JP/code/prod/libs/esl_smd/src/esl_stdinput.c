/****************************************************************************/
/*      Standard input output for yield and volatility curves.              */
/****************************************************************************/
/*      STDINPUT.c                                                          */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "esl_stdinput.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"
#include "esl_util.h"

/*****	FindAndSkipComLine  *************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine (FILE         *stream,
                            char const   *Label,
                            char const   *Routine,
                            char const   *FileName)
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

/*****	Term_Input_W  *******************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int  Term_Input_W (T_CURVE  *t_curve    /** (O) Structure of zero curve data  */
                  ,char const* FileName /** (I) File name including extension */
		  )
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
/**
*       Check validity of DR Wrapper term structure inputs.
*/
int     Term_Check_W (
		T_CURVE  *t_curve  /** (I) Structure of zero curve data */
		)
{

    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */

    char    ErrorMsg[MAXBUFF];
    

    if (Dateok(t_curve->ValueDate))
    {
        DR_Error("Incorrect format for value date [%d]", t_curve->ValueDate);
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
        if ((t_curve->Zero[i] < -100.) || (t_curve->Zero[i] < .000001) || (t_curve->Zero[i] > 100.))
        {
            DR_Error("Zero rates out of range!");
            printf("my test i=%d, nbZero=%d, zero point=%lf\n",i,t_curve->NbZero, t_curve->Zero[i]);
            goto RETURN;
        }
    }

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Term_Check_W */



/*****	BaseVol_Input_W  ****************************************************/
/**
*       Read base volatility input and check validity of input.
*/
int     BaseVol_Input_W (
            MKTVOL_DATA  *mktvol_data   /** (O) Volatility data               */
           ,char   const* FileName      /** (I) File name including extension */
	    )
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
/**
*       Read swaption volatility input and check validity of input.
*/
int     SwapVol_Input_W (int          *NbRows,      /** (O) Volatility matrix */
                         int          *NbCol,
                         long         **Expiry,
                         long         **FwdMat,
                         double       ***VolMatrix,
                         char         **DoMoY,     
                         char     const *FileName  /** (I) File name */
		)
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
    *VolMatrix = (double **) DR_Matrix (DOUBLE,  0, *NbRows, 0, *NbCol);
    *DoMoY     = (char *)    DR_Array  (CHAR,    0, *NbRows);

    if    ((*Expiry    == NULL)
        || (*FwdMat    == NULL)
        || (*VolMatrix == NULL)
        || (*DoMoY     == NULL))
	
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
	   	char token[64];

		readerror = fscanf (stream, "%s ", token);
        if (readerror != 1)
		{        
			sprintf (ErrorMsg, "Could not read expiry in file %s! (SwapVol_Input)", FileName);
			DR_Error(ErrorMsg);
			goto RETURN;
		}

      
		(*Expiry)[i] = atol(token);

		if (strchr(token, 'D') != NULL)
		{       
			(*DoMoY)[i] = 'D';
		}
		else if (strchr(token, 'Y') != NULL)
		{
			(*DoMoY)[i] = 'Y';
		}
		else
        {
			(*DoMoY)[i] = 'M';
        }   
	
        for (j = 0; j < *NbCol; j++)
        {
			readerror = fscanf (stream, "%s", token);
            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read volatility in file %s! (SwapVol_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            (*VolMatrix)[i][j] = atof(token) / 100.;
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


/*****  MktVol_Input_W  ******************************************************/
/**
*       Utility routine converting an index name to a series of option expi-
*       ries, index maturities and volatilities taken from either the base
*       volatility curve or the swaption matrix.
*/
int     MktVol_Input_W (
		MKTVOL_DATA*   mktvol_data, /**< (O) Volatility data      */
                char const*    Index,       /**< (I) Index to calibrate   */
                T_CURVE const* t_curve,     /**< (I) Term structure data  */
                char const*    BaseVolFile, /**< (I) Base vol curve file  */
                char const*    SwapVolFile) /**< (I) Swaption matrix file */
{

    char     IndexL[MAXINDEX];   /* Local copy of index */
    char    *StrIdx = NULL;     
    long    IdxMat;             /* Maturity of the index (final or forward) */
    long    Mat;
    int     NbRows=0;           /* Swaption volatility matrix */
    int     NbCol=0;
    long    *Expiry = NULL;
    long    *FwdMat = NULL;
    double  **VolMatrix = NULL;
    char    *DoMoY = NULL;
    int     i, j;
    int     status = FAILURE;   /* Error status = FAILURE initially         */

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


    /* Search for * character in calibr ation index name */
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
                                &DoMoY,
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


        mktvol_data->NbVol = NbRows;

        for (i = 0; i < mktvol_data->NbVol; i++)
        {                                                                       
            /* 
             *  VolDate and SwapSt are identical: i.e. we don't calibrate 
             *  options on forward starting swaps (e.g. mid curve options) 
             */
			if( DoMoY[i] == 'M')
                mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
			else if (DoMoY[i] == 'Y')
				mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, 12 * Expiry[i], 1L);
		    else if (DoMoY[i] == 'D')
			    mktvol_data->VolDate[i] = Nxtday (mktvol_data->BaseDate, Expiry[i]);

            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], IdxMat, 1L);
            mktvol_data->Vol[i]     = VolMatrix[i][j];
            mktvol_data->VolUsed[i] = TRUE;
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
                                &DoMoY,
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

            if( DoMoY[i] == 'M')
                mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
			else if (DoMoY[i] == 'Y')
				mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, 12 * Expiry[i], 1L);
		    else if (DoMoY[i] == 'D')
			    mktvol_data->VolDate[i] = Nxtday (mktvol_data->BaseDate, Expiry[i]);

           
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
    StrIdx = strchr(IndexL, 'm');

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL);


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
    Free_DR_Array  (DoMoY     , CHAR,   0, NbRows);
    Free_DR_Array  (FwdMat    , LONG,   0, NbCol);
    Free_DR_Matrix (VolMatrix , DOUBLE, 0, NbRows, 0, NbCol);

    return (status);

}  /* MktVol_Input_W */




/*****	VolDiag_Input_W  ****************************************************/
/**
 *       Read volatility input
 */
int     VolDiag_Input_W (
            long    *BaseDate    /** (O) Volatility data */
           ,int     *VolUnit
           ,int     *NbVol
           ,long    *VolDate
           ,long    *SwapSt
           ,long    *SwapMat
           ,double  *Vol
           ,char    *VolType
           ,char    *FileName    /** (I) File name including extension */
	   )
{
    int     i;
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;

    if ((BaseDate == NULL) ||
        (VolUnit  == NULL) ||
        (NbVol    == NULL) ||
        (VolDate  == NULL) ||
        (SwapSt   == NULL) ||
        (SwapMat  == NULL) ||
        (Vol      == NULL) ||
        (VolType  == NULL) ||
        (FileName == NULL)) goto RETURN;

    /* Open and read the volatility data file */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        DR_Error("Could not open file %s! (VolDiag_Input_W)", FileName);
        
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "start date", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", BaseDate);
    if (readerror != 1)
    {        
        DR_Error("Could not read start date in file %s! (VolDiag_Input_W)", FileName);
        
        goto RETURN;
    }
    
    if (FindAndSkipComLine (stream, "volatility unit", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", VolUnit);
    if (readerror != 1)
    {        
        DR_Error("Could not read vol unit in file %s! (VolDiag_Input_W)", FileName);
        
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "Nb of vol points", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbVol);
    if (readerror != 1)
    {        
        DR_Error("Could not read number of vols in file %s! (VolDiag_Input_W)", FileName);
        
        goto RETURN;
    }


    if ((*NbVol > MAXNBDATE) ||
        (*NbVol < 0))
    {        
        DR_Error("Nb of vols in file %s must be >= 0 and <= %d! (VolDiag_Input_W)", FileName, MAXNBDATE);
        
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "vol dates and rates", "VolDiag_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < *NbVol; i++)
    {
        readerror = fscanf (stream, "%ld %ld %ld %lf %c\n", 
                            &(VolDate[i]),
                            &(SwapSt[i]),
                            &(SwapMat[i]),
                            &(Vol[i]),
                            &(VolType[i]));

        if (readerror != 5)
        {        
            DR_Error("Could not read vol dates/rate/type #%d in file %s! (VolDiag_Input_W)", i+1, FileName);
            
            goto RETURN;
        }

        Vol[i] /= 100.;
    }

    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status == FAILURE)
    {
        DR_Error("VolDiag_Input_W: failed!");
    }
    

    return (status);

}  /* VolDiag_Input_W */


/*****  MktVol_Input_Plus_W  ************************************************/
/**
        Utility routine converting an index name to a series of option expi-
        ries, index maturities and volatilities taken from either the base
        volatility curve or the swaption matrix.
 
        Same as MktVol_Input_W except this function allows the reading 
        from ir_voldiag_i.dat
        The reading logic is as follows:
        - treats nil calibration first
        - checks if it's basevol caibration and read from basevol.dat
        - checks availability of d(f)swapvol.dat in the folder
        - if not available, then reads from ir_voldiag_i.dat
        - otherwise, always reads from d(f)swapvol.dat
*/
int MktVol_Input_Plus_W (
	MKTVOL_DATA *mktvol_data     /** (O) Volatility data               */
       ,char       *Index           /** (I) Index to calibrate             */
       ,T_CURVE    *t_curve         /** (I) Term structure data            */
       ,char       *BaseVolFile     /** (I) Base vol curve file            */
       ,char       *SwapVolFile     /** (I) Swaption matrix file           */
       ,char       *SwapVolFileSRM3 /** (I) Swaption matrix file from SRM3 */
		)
{
    char     IndexL[MAXINDEX];   /* Local copy of index */
    char    *StrIdx = NULL;     
    long     IdxMat;             /* Maturity of the index (final or forward) */
    long     Mat;
    int      NbRows = 0;         /* Swaption volatility matrix */
    int      NbCol = 0;
    long    *Expiry = NULL;
    long    *FwdMat = NULL;
    FILE    *stream = NULL;
    double  **VolMatrix = NULL;
    char     *DoMoY = NULL;
    int      i, j;
	int      SRM3Ind;            /* (0 or 1). 1 implies we need to read from ir_vol_diag_i.dat */
    int      status = FAILURE;   /* Error status = FAILURE initially         */

    int      is_yFix, is_yCms, is_BaseVol;

	char    VolType[MAXNBDATE];     /* Caplet or swaption vol          */
	int     VolUnit;                /* Volatility unit (Form SRM3)     */

    
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

    is_yFix    = (strstr(IndexL, "yFix") != NULL);
    is_yCms    = (strstr(IndexL, "yCms") != NULL);
    is_BaseVol = ((strstr(IndexL, "m") != NULL) && (!is_yCms));

    if ((!is_yFix) && (!is_yCms) && (!is_BaseVol))
    {
        DR_Error("Invalid calibration index.");
        goto RETURN;
    }

	/*Base vol calibration, still read old files*/
	/* Base vol indices */
    if (is_BaseVol)
	{
        StrIdx = strchr(IndexL, 'm');
		*StrIdx = '\0';
		IdxMat = atoi(IndexL);

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
		             "index (MktVol_Input_Plus_W)!");
		 goto RETURN;
		}


		for (i = 0; i < mktvol_data->NbVol; i++)
		{
		 mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
		 mktvol_data->SwapMat[i] = Nxtmth(mktvol_data->VolDate[i],IdxMat,1L);
		 mktvol_data->VolUsed[i] = TRUE;	
		}


		if (MktVol_Check_W (mktvol_data)  == FAILURE)
		{
		 goto RETURN;
		}
		return (SUCCESS);
	}


	/* Other cases Fix and CMS, try to open HYB3 files first*/
	/* Now try to open Hyb3 swapVol file */
    stream = fopen (SwapVolFile, "r");

    if (stream == NULL)
    {
		DR_Error("%s not available.", SwapVolFile);
		SRM3Ind = 1;    /* =1 implies we need the ir_voldaig_0.dat */
    }
	else
	{
		SRM3Ind = 0;    /* =0 implies we need read from d(f)swapvol.dat */
		fclose(stream);
	}

	if(SRM3Ind == 1)
	{

	    /**********************/
	    /*  Read Vol diagonal */
		/**********************/
		if (VolDiag_Input_W (
                        &(mktvol_data->BaseDate),
                        &VolUnit,
                        &(mktvol_data->NbVol),
                        mktvol_data->VolDate,
                        mktvol_data->SwapSt,
                        mktvol_data->SwapMat,
                        mktvol_data->Vol,
                        VolType,
                        SwapVolFileSRM3) == FAILURE) goto RETURN;

		/* Fill vector volUsed (HYB3 vector)*/
		for(i = 0; i < mktvol_data->NbVol; i++)
		{
			mktvol_data->VolUsed[i] = TRUE;

		}

		/* Calib is off either CalibFlag from deal = 'N' or NbVol from env = 0 */
		if (mktvol_data->NbVol <= 0)
		{
		    DR_Error("Nb Vol has to be > 0!");
            goto RETURN;
		}
		else 
		{
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
			/* Set Freq and DCC for vols */
			/* only need to check the first 'Vol Type' entry because 
			in irInputCheck we will check that all Vol Types are equal */
    		for(i = 0; i < mktvol_data->NbVol; i++)
	    	{
                if (VolType[i] != 'S')
	    		{
		    		 DR_Error("Some VolType entry is not S in ir_vol_diag_i.dat.");
			    	 goto RETURN; 
			    }
		    }
		}
	}

	if(SRM3Ind == 0)
	{
		/* Cms indices */
		if (is_yCms)
		{
            StrIdx = strstr (IndexL, "yCms");
		    *StrIdx = '\0';
		    IdxMat = atoi(IndexL) * 12;


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
			if (SwapVol_Input_W ( &NbRows,
			                      &NbCol,
			                      &Expiry,
			                      &FwdMat,
			                      &VolMatrix,
                                              &DoMoY,
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
			    DR_Error("Cms index is not in swaption matrix (MktVol_Input_Plus_W)!");
			    goto RETURN;
			}


			mktvol_data->NbVol = NbRows;

			for (i = 0; i < mktvol_data->NbVol; i++)
			{                                                                       
				/* 
				*  VolDate and SwapSt are identical: i.e. we don't calibrate 
				*  options on forward starting swaps (e.g. mid curve options) 
				*/
                if( DoMoY[i] == 'M')
					mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
			    else if (DoMoY[i] == 'Y')
					mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, 12 * Expiry[i], 1L);
		        else if (DoMoY[i] == 'D')
					mktvol_data->VolDate[i] = Nxtday (mktvol_data->BaseDate, Expiry[i]);

			
				mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
				mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i],IdxMat,1L);
				mktvol_data->Vol[i]     = VolMatrix[i][j];
				mktvol_data->VolUsed[i] = TRUE;
			}


			if (MktVol_Check_W (mktvol_data)  == FAILURE)
			{
			    goto RETURN;
			}

			goto DONE;
		}


		/* Final maturity indices */
		if (is_yFix)
		{
         StrIdx = strstr (IndexL, "yFix");
		 *StrIdx = '\0';
		 IdxMat = atoi(IndexL) * 12;


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
                                         &DoMoY,
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

            if( DoMoY[i] == 'M')
                mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, Expiry[i], 1L);
			else if (DoMoY[i] == 'Y')
				mktvol_data->VolDate[i] = Nxtmth (mktvol_data->BaseDate, 12 * Expiry[i], 1L);
		    else if (DoMoY[i] == 'D')
			    mktvol_data->VolDate[i] = Nxtday (mktvol_data->BaseDate, Expiry[i]);
         
            mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];
            mktvol_data->SwapMat[i] = Nxtmth (mktvol_data->VolDate[i], Mat, 1L);
            mktvol_data->VolUsed[i] = TRUE;

        }  /* for i */
        
        if (i == 0)
        {        
            DR_Error ("nyFix calibration falls outside swaption matrix "
                        "(MktVol_Input_Plus_W)!");
            goto RETURN;
        }

        mktvol_data->NbVol = i;


        if (MktVol_Check_W (mktvol_data)  == FAILURE)
        {
            goto RETURN;
        }

		goto DONE;
		}

        if (is_BaseVol)
        {
            DR_Error("Incorrect logic. Should not get here.");
            goto RETURN;
        }

    }


    DONE:

    status = SUCCESS;

    RETURN:

    Free_DR_Array  (Expiry    , LONG,   0, NbRows);
    Free_DR_Array  (DoMoY     , CHAR,   0, NbRows);
    Free_DR_Array  (FwdMat    , LONG,   0, NbCol);
    Free_DR_Matrix (VolMatrix , DOUBLE, 0, NbRows, 0, NbCol);

    return (status);

}  /* MktVol_Input_Plus_W */



/*****	MktVol_Check_W  *****************************************************/
/**
*       Check validity of base volatility inputs.
*/
int     MktVol_Check_W (
		MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
		)
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
            DR_Error("Incorrect format for volatility date [%d] (MktVol_Check_W)!", mktvol_data->VolDate[i]);
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
        DR_Error("Incorrect format for base date [%d] (MktVol_Check_W)!", mktvol_data->BaseDate);
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
            DR_Error("Benchmark swap matures %d before swap start %d (MktVol_Check_W)!", 
		mktvol_data->SwapSt[i], mktvol_data->SwapMat[i]);
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

/** 
*   Print probability in output file for wrapper.
*/
void    printProbabilityInFile(
		OPT_OUT_DATA *opt_out_data)
{
    FILE*    fproba;
    char     ErrorMsg[MAXBUFF];
    int      i,j;

    fproba = fopen ("probability.dat", "w");
    
    if (fproba == NULL)
    {
        sprintf (ErrorMsg, "Could not open file probability! (main)");
        DR_Error (ErrorMsg);
        return;
    }

    fprintf (fproba, "Option price   \t%16.4f\n", opt_out_data->Option);

    for (i=0; i<ESL_NB_EVENT ; i++)
    {
        if (opt_out_data->prob_calc[i] == TRUE)
        {
            if (i == EXER_EVENT) 
            {
                fprintf (fproba, "Exer Prob      \t%12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].TotalEventProb);
                fprintf (fproba, "Exer Time      \t%12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].ExpEventTime);
                fprintf (fproba, "Exer Time s.d. \t%12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].StdEventTime);
                fprintf (fproba, "Fugit          \t%12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].Fugit);
                fprintf (fproba, "Exercise probability schedule \t \n");
            }
	    else if (i == KO_EVENT)
            {
                fprintf (fproba, "KO Prob      \t%12.6f\n", opt_out_data->prob_out_data[KO_EVENT].TotalEventProb);
                fprintf (fproba, "KO Time      \t%12.6f\n", opt_out_data->prob_out_data[KO_EVENT].ExpEventTime);
                fprintf (fproba, "KO Time s.d. \t%12.6f\n", opt_out_data->prob_out_data[KO_EVENT].StdEventTime);
                fprintf (fproba, "Fugit        \t%12.6f\n", opt_out_data->prob_out_data[KO_EVENT].Fugit);
                fprintf (fproba, "KO probability schedule \t \n");
            }

            fprintf (fproba, "Date           \tProbability\n");

            for (j=0; j<opt_out_data->prob_out_data[i].Count; j++)
                fprintf (fproba, "%ld\t%8.6f\n", opt_out_data->prob_out_data[i].EventDate[j], 
                                opt_out_data->prob_out_data[i].EventProb[j]);
        }
    }

    fclose (fproba);

    return;
}

