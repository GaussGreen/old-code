/****************************************************************************/
/*      Standard input output for yield and volatility curves.              */
/****************************************************************************/
/*      STDINPUT.c                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "esl_stdinput.h"
#include "esl_error.h"
#include "esl_date.h"
#include "esl_alloc.h"
#include "esl_util.h"
#include "esl_types.h"
#include "irx/irxflow.h"
#include "irx/zerocurve.h"
#include "irx/strutils.h"
#include "esl_zeros.h"

#define BIG_STRING 50000
#define SMALL_STRING 100

/*****  FindAndSkipSectionLine  *************************************************/
/*
*       This function does:
*       1) Mode 0: Find and skip mode
*           Skip everything until it finds line start with '###' (ignore spaces)
*           Flag an error message if it can't find
*       2) Mode 1: Skip mode
*           Skip all the empty lines and spaces
*           Read the first non-empty character and check it is '#' followed by '##' 
*           Flag an error message if it isn't.
*       3) Mode 2: Find and Skip mode with no error message printed
*           Same as Mode 0, but no error message is printed if failed to find
*           next section line.
*           Note the status is still FAILURE in this case
*
*       The routine returns with the file pointer at the start of the next line
*/
int     FindAndSkipSectionLine (int      Mode,
                                FILE    *stream,
                                char const  *SectionLabel,
                                char const  *Routine,
                                char const  *FileName
                                )
{
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i;
    int     length;

    if (stream == NULL) goto RETURN;

    if ((Mode != 0) && (Mode != 1) && (Mode != 2))
    {
        DR_Error("FindAndSkipSectionLine: Mode has to be 0 or 1 or 2!");
        goto RETURN;
    }

    while (fgets (string, MAXBUFF, stream) != NULL)
    {
        length = strlen(string);

        for (i = 0; i < length; i++)
        {
            if (string[i] == '#')
            {
                if ((i > length - 3) || 
                    (string[i+1] != '#') || 
                    (string[i+2] != '#')) 
                {
                    if (Mode == 0) break;
                    else goto RETURN;
                }
                else
                {
                    status = SUCCESS;
                    goto RETURN;
                }

            }
            else if ((string[i] != ' ' ) && 
                     (string[i] != '\n') && 
                     (string[i] != '\t')) 
            {
                if ((Mode == 0) || (Mode == 2)) break;
                else goto RETURN;
            }
        }
    }  

    RETURN:


    if (status == FAILURE)
    {
        if (Mode != 2)
        {
            DR_Error ("%s: %s section line expected and not found in file"
                      " %s!", Routine, SectionLabel, FileName);
        }
    }
        
    return (status);

}  /* FindAndSkipSectionLine */

/*****  FindAndSkipComLine  *************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine (FILE         *stream,
                            char const   *Label,
                            char const   *Routine,
                            char const   *FileName)
{
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */

    if (fgets (string, MAXBUFF, stream) == NULL)
    {
        DR_Error ("%s: %s comment line expected and not found in file %s!", 
                  Routine, Label, FileName);
        goto RETURN;
    }
    else
    {
        if (string[0] != '#')
        {
            DR_Error ("%s: expected '#' in %s comment line in file %s!", 
                       Routine, Label, FileName);
            goto RETURN;
        }
    }  /* if then else */

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* FindAndSkipComLine */



/*****  FindAndSkipComLine_2  *************************************************/
/*
*       This function does:
*       1) Mode 0: Find and skip mode
*           Skip everything until it finds line start with '#' (ignore spaces)
*           Flag an error message if it can't find
*       2) Mode 1: Skip mode
*           Skip all the empty lines and spaces
*           Read the first non-empty character and check it is '#' 
*           Flag an error message if it isn't.
*       3) Mode 2: Find and Skip mode with no error message printed
*           Same as Mode 0, but no error message is printed if failed to find
*           next section line.
*           Note the status is still FAILURE in this case
*
*       The routine returns with the file pointer at the start of the next line
*/
int     FindAndSkipComLine_2(int      Mode,
                            FILE    *stream,
                            char const *SectionLabel,
                            char const *Routine,
                            char const *FileName)
{
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    int     i;
    int     length;

    if (stream == NULL) goto RETURN;

    if ((Mode != 0) && (Mode != 1) && (Mode != 2))
    {
        DR_Error("FindAndSkipComLine_2: Mode has to be 0, 1 or 2!");
        goto RETURN;
    }

    while (fgets (string, MAXBUFF, stream) != NULL)
    {
        length = strlen(string);

        for (i = 0; i < length; i++)
        {
            if (string[i] == '#')
            {
                status = SUCCESS;
                goto RETURN;
            }
            else if ((string[i] != ' ' ) && 
                     (string[i] != '\n') && 
                     (string[i] != '\t')) 
            {
                if (Mode == 0) break;
                else goto RETURN;
            }
        }
    }  

    RETURN:


    if (status == FAILURE)
    {
        if (Mode != 2)
        {
            DR_Error("%s: %s section line expected and not found in file"
                    " %s!", Routine, SectionLabel, FileName);
        }
    }
        
    return (status);

}  /* FindAndSkipComLine_2 */




/*****  Term_Input_W  *******************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int  Term_Input_W (T_CURVE*   crv   /** (O) Structure of zero curve data  */
                  ,char const*      fileName /** (I) File name including extension */)
{

    long    LastDate;           /* Last date in zero curve */
    int     i, cnt;
    int     readerror;          /* Reading error status */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;

    long    dates[MAXNBDATE];
    double  zeros[MAXNBDATE];


    /* Open the yield curve data file (see termodel.h) */
    stream = fopen (fileName, "r");

    if (stream == NULL)
    {
        DR_Error ("Could not open file %s! (Term_Input_W)", fileName);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "value date", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", &crv->Today); 
    if (readerror != 1)
    {        
        DR_Error ("Could not read value date in %s! (Term_Input_W)", fileName);
        goto RETURN;
    }
    crv->Today = IRDateFromYMDDate(crv->Today);

    /* Today's date and spot days are not available */
    crv->ValueDate = crv->Today;
    crv->SpotDays  = 0;


    if (FindAndSkipComLine (stream, "money market basis", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &crv->MMB);
    if (readerror != 1)
    {        
        DR_Error ("Could not read money market basis in file %s! (Term_Input_W)", fileName);
        goto RETURN;
    }     
   
    if (FindAndSkipComLine (stream, "yield curve frequency", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &crv->SwapFreq);
    if (readerror != 1)
    {        
        DR_Error ("Could not read benchmark swap frequency in file %s! (Term_Input_W)", fileName);
        goto RETURN;
    }
    
    if (FindAndSkipComLine (stream, "swap DCC", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }


    readerror = fscanf (stream, "%s \n", crv->SwapDCC);
    if (readerror != 1)
    {        
        DR_Error ("Could not read swap DCC in file %s! (Term_Input_W)", fileName);
        goto RETURN;
    }


    if (FindAndSkipComLine (stream, "number of zeros", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &cnt);
    if (readerror != 1)
    {        
        DR_Error ("Could not read number of zeros in file %s! (Term_Input_W)", fileName);
        goto RETURN;
    }
    
    if (cnt > MAXNBDATE)
    {        
        DR_Error ("Nb of rates in file %s exceeds maximum of %d! (Term_Input_W)", fileName, MAXNBDATE);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "zero dates and rates", "Term_Input_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < cnt; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \n", &dates[i], &zeros[i]); 
        if (readerror != 2)
        {        
            DR_Error ("Could not read zero date and rate #%d in file %s! (Term_Input_W)", i + 1, fileName);
            goto RETURN;
        }
        dates[i] = IRDateFromYMDDate(dates[i]);
        zeros[i] /= 100.;
    }
    

    LastDate = Nxtmth (crv->ValueDate, 1200L, 1L);                          /* We require at least 100 years of zero curve for swaption vol bootstrapping */

    
    if (dates[cnt-1] < LastDate)                        /* If the zero curve does not extend up to last date we add an extra point */
    {
        cnt += 1;                       
        dates[cnt-1] = LastDate;
        zeros[cnt-1] = zeros[cnt-2];    /* Flat zero curve */
    }
            
    /* copy local curve object to the output */
#ifdef ESL_NEW_CURVE
    if (irxZeroCurveConstructFromRates(
                crv,
                crv->ValueDate,
                cnt,
                dates, 
                zeros,
                IRX_ANNUAL_RATE,
                IRX_ACT_365F) != SUCCESS)
        goto RETURN;


#else
    for (i = 0; i < cnt; i++)
    {
        crv->ZeroDate[i] = dates[i];
        crv->Zero[i] = zeros[i];
    }
    crv->NbZero = cnt;
    crv->InterpType =  EslGetZeroInterpolation();  /* should check srm3 is using linear */
  
#endif

    if (Term_Check_W (crv) == FAILURE)                    
        goto RETURN;

    status = SUCCESS;

    
    RETURN:

    if (stream != NULL)
        fclose (stream);

    return (status);

}  /* Term_Input_W */



/*****  Term_Input_New_W  *******************************************************/
/**
*       Read term structure input from ir_curve and check validity of input.
*/
int  Term_Input_New_W (
             T_CURVE*   crv,            /** (O) Structure of zero curve data  */
             char const*      fileName) /** (I) File name including extension */
          
{

    long    LastDate;           /* Last date in zero curve */
    int     i, cnt;
    int     readerror;          /* Reading error status */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;
    char    ErrorMsg[MAXBUFF];


    char    CurveDCC[MAXBUFF], CurveFreq;

    long    dates[MAXNBDATE];
    double  zeros[MAXNBDATE];

    /* Open the yield curve data file (see termodel.h) */
    stream = fopen (fileName, "r");
   

    if (stream == NULL)
    {
        DR_Error ("Could not open file %s! (Term_Input_New_W)", fileName);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "value date", "Term_Input_New_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", &crv->ValueDate); 
    if (readerror != 1)
    {        
        DR_Error ("Could not read value date in %s! (Term_Input_New_W)", fileName);
        goto RETURN;
    }
    crv->ValueDate = IRDateFromYMDDate(crv->ValueDate);

    /* Today's date and spot days are not available */
    crv->Today = crv->ValueDate;
    crv->SpotDays  = 0;
    
    /* day count convention for zero rates */
    if (FindAndSkipComLine (stream, "dcc for zero rates", "Term_Input_New_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%s \n", CurveDCC);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read dcc for zero rates in file %s! (Term_Input_New_W)", fileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    /* only ACT/365 allowed */
    if (strstr(CurveDCC, "ACT/365F") == NULL)
    {
        sprintf (ErrorMsg, "dcc in file %s has to be ACT/365F' ! (Term_Input_New_W)", fileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* compounding frequency for zero rates */    
    if (FindAndSkipComLine (stream, "compounding frequency", "Term_Input_New_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(CurveFreq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read compounding frequency in file %s! (Term_Input_New_W)", fileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    /* only annual frequency allowed */
    if (CurveFreq != 'A')
    {
        sprintf (ErrorMsg, "Only curves of annual frequency sre supported! (Term_Input_New_W)");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }
        
    
    if (FindAndSkipComLine (stream, "number of zeros", "Term_Input_New_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &cnt);
    if (readerror != 1)
    {        
        DR_Error ("Could not read number of zeros in file %s! (Term_Input_New_W)", fileName);
        goto RETURN;
    }
        
    if (cnt > MAXNBDATE)
    {        
        DR_Error ("Nb of rates in file %s exceeds maximum of %d! (Term_Input_New_W)", fileName, MAXNBDATE);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "zero dates and rates", "Term_Input_New_W", fileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < cnt; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \n", &dates[i], &zeros[i]); 
        if (readerror != 2)
        {        
            DR_Error ("Could not read zero date and rate #%d in file %s! (Term_Input_New_W)", i + 1, fileName);
            goto RETURN;
        }
        dates[i] = IRDateFromYMDDate(dates[i]);
        zeros[i] /= 100.;
    }
        

    LastDate = Nxtmth (crv->ValueDate, 1200L, 1L);  
    /* We require at least 100 years of zero curve for swaption vol bootstrapping */

        
    if (dates[cnt-1] < LastDate)     /* If the zero curve does not extend up to last date we add an extra point */
    {
        cnt += 1;                       
        dates[cnt-1] = LastDate;
        zeros[cnt-1] = zeros[cnt-2];    /* Flat zero curve */
    }
                
    /* copy local curve object to the output */
#ifdef ESL_NEW_CURVE
    if (irxZeroCurveConstructFromRates(
                crv,
                crv->ValueDate,
                cnt,
                dates, 
                zeros,
                IRX_ANNUAL_RATE,
                IRX_ACT_365F) != SUCCESS)
        goto RETURN;
#else
    for (i = 0; i < cnt; i++)
    {
        crv->ZeroDate[i] = dates[i];
        crv->Zero[i] = zeros[i];
    }
    crv->NbZero = cnt;
    crv->InterpType = EslGetZeroInterpolation();  /* should check interp type for srm3 */
   
#endif

   if (Term_Check_W (crv) == FAILURE)                    
        goto RETURN;


    status = SUCCESS;
        
    RETURN:

    if (stream != NULL)
        fclose (stream);

    return (status);

}  /* Term_Input_New_W */



/*****  Term_Check_W  *******************************************************/
/**
*       Check validity of DR Wrapper term structure inputs.
*/
int     Term_Check_W (T_CURVE const*  t_curve  /** (I) Structure of zero curve data */
        )
{

    int     status = FAILURE;   /* Error status = FAILURE initially */

#ifndef ESL_NEW_CURVE
    int   i;
#endif

    if (Dateok(t_curve->ValueDate))
    {
        DR_Error("Incorrect format for value date [%ld]", t_curve->ValueDate);
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

    if (   strcmp(t_curve->SwapDCC,"ACT") 
        && strcmp(t_curve->SwapDCC,"365") 
        && strcmp(t_curve->SwapDCC,"360"))
    {
        DR_Error("Specify 'ACT' for Actual, '365' for 365 Fixed or '360' for swap DCC!");
        goto RETURN;
    }

#ifndef ESL_NEW_CURVE

    if (t_curve->NbZero > MAXNBDATE)
    {        
        DR_Error ("Nb of zero rates exceeds maximum of %d!", MAXNBDATE);
        goto RETURN;
    }

    for (i = 1; i < t_curve->NbZero; i++) {
        if (t_curve->ZeroDate[i] <= t_curve->ZeroDate[i-1])
        {
            DR_Error("Zero dates must be entered in ascending order!");
            goto RETURN;
        }
    }

    for (i = 0; i < t_curve->NbZero; i++)
    {
        /*if ((t_curve->Zero[i] < -100.) || (t_curve->Zero[i] < .000001) || (t_curve->Zero[i] > 100.))*/
        if ((t_curve->Zero[i] < -100.) || (t_curve->Zero[i] > 100.))
        {
            DR_Error("Zero rates out of range!");
            printf("my test i=%d, nbZero=%d, zero point=%lf\n",i,t_curve->NbZero, t_curve->Zero[i]);
            goto RETURN;
        }
    }

    if (t_curve->InterpType != ESL_INTERP_FLATFWD && t_curve->InterpType != ESL_INTERP_LINEAR)
    {
        DR_Error("InterpType is not set");
        goto RETURN;
    }

#endif

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Term_Check_W */



/*****  BaseVol_Input_W  ****************************************************/
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
    FILE    *stream = NULL;


    /* Open the base volatility curve data file (see termodel.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        DR_Error ("Could not open file %s! (BaseVol_Input)", FileName);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "base vol frequency", "BaseVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(mktvol_data->Freq));
    if (readerror != 1)
    {        
        DR_Error ("Could not read base vol frequency in file %s! (BaseVol_Input)", FileName);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "number of base vols", "BaseVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(mktvol_data->NbVol));
    if (readerror != 1)
    {        
        DR_Error ("Could not read number of base vols in file %s! (BaseVol_Input)", FileName);
        goto RETURN;
    }
    
    if (mktvol_data->NbVol > MAXNBDATE)
    {        
        DR_Error ("Nb of base vols in file %s exceeds maximum of %d! (BaseVol_Input)", FileName, MAXNBDATE);
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
            DR_Error ("Could not read base vol date and rate #%d in file %s! (BaseVol_Input)", i+1, FileName);
            goto RETURN;
        }
        mktvol_data->VolDate[i] = IRDateFromYMDDate(mktvol_data->VolDate[i]);
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


/*****  SwapVol_Input_W  ****************************************************/
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
    FILE    *stream = NULL;

    /* Open the swaption volatility curve data file (see fix123.h) */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        DR_Error ("Could not open file %s! (SwapVol_Input)", FileName);
        goto RETURN;        
    }
        
    if (FindAndSkipComLine (stream, "number of rows", "SwapVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbRows);
    if (readerror != 1)
    {        
        DR_Error ("Could not read number of rows in file %s! (SwapVol_Input)", FileName);
        goto RETURN;
    }
        
    if (*NbRows > MAXNBDATE)
    {        
        DR_Error ("Nb of expiry in file %s exceeds maximum of %d! (SwapVol_Input)", FileName, MAXNBDATE);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "number of columns", "SwapVol_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbCol);
    if (readerror != 1)
    {        
        DR_Error ("Could not read number of columns in file %s! (SwapVol_Input)", FileName);
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
            DR_Error ("Could not read forward maturity in file %s! (SwapVol_Input)", FileName);
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
            DR_Error ("Could not read expiry in file %s! (SwapVol_Input)", FileName);
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
                DR_Error ("Could not read volatility in file %s! (SwapVol_Input)", FileName);
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
        MKTVOL_DATA*            mktvol_data, /**< (O) Volatility data      */
        char const*             Index,       /**< (I) Index to calibrate   */
        T_CURVE const*          t_curve,     /**< (I) Term structure data  */
        char const*             BaseVolFile, /**< (I) Base vol curve file  */
        char const*             SwapVolFile) /**< (I) Swaption matrix file */
{

    
    char    IndexL[MAXINDEX];   /* Local copy of index */
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
        mktvol_data->BaseDate = t_curve->Today;
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
            mktvol_data->VolType[i] = 'S'; /* set it to swaption */
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
            mktvol_data->VolType[i] = 'S'; /* set it to swaption */

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
            mktvol_data->VolType[i] = 'C'; /* set it to cap */
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


/*****  VolDiag_Input_New_W  ****************************************************/
/**
 *       Read volatility input
 */
int     VolDiag_Input_New_W (
            long    *   BaseDate    /** (O) Volatility data */
           ,int     *   VolUnit
           ,int     *   NbVol
           ,long    *   VolDate
           ,long    *   SwapSt
           ,long    *   SwapMat
           ,double  *   Vol
           ,char    *   VolType
           ,int     *   VolSkipFlag
           ,char const* FileName    /** (I) File name including extension */
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
        (VolSkipFlag == NULL) ||
        (FileName == NULL)) goto RETURN;

    /* Open and read the volatility data file */
    stream = fopen (FileName, "r");

    if (stream == NULL)
    {
        DR_Error("Could not open file %s! (VolDiag_Input_New_W)", FileName);
        
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "start date", "VolDiag_Input_New_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld \n", BaseDate);
    if (readerror != 1)
    {        
        DR_Error("Could not read start date in file %s! (VolDiag_Input_New_W)", FileName);
        
        goto RETURN;
    }
    *BaseDate = IRDateFromYMDDate(*BaseDate);
    
    if (FindAndSkipComLine (stream, "volatility unit", "VolDiag_Input_New_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", VolUnit);
    if (readerror != 1)
    {        
        DR_Error("Could not read vol unit in file %s! (VolDiag_Input_New_W)", FileName);
        
        goto RETURN;
    }

    *VolUnit = 1 - *VolUnit;

    if (FindAndSkipComLine (stream, "Nb of vol points", "VolDiag_Input_New_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", NbVol);
    if (readerror != 1)
    {        
        DR_Error("Could not read number of vols in file %s! (VolDiag_Input_New_W)", FileName);
        
        goto RETURN;
    }

    if ((*NbVol > MAXNBDATE) ||
        (*NbVol < 0))
    {        
        DR_Error("Nb of vols in file %s must be >= 0 and <= %d! (VolDiag_Input_New_W)", FileName, MAXNBDATE);
        
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "vol dates and rates", "VolDiag_Input_New_W", FileName) == FAILURE)
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
            DR_Error("Could not read vol dates/rate/type #%d in file %s! (VolDiag_Input_New_W)", i+1, FileName);
            goto RETURN;
        }

        VolDate[i] = IRDateFromYMDDate(VolDate[i]);
        SwapSt[i] = IRDateFromYMDDate(SwapSt[i]);
        SwapMat[i] = IRDateFromYMDDate(SwapMat[i]);

        Vol[i] /= 100.;
    }

    if (FindAndSkipComLine_2 (2, stream, "vol points skip flag", "VolDiag_Input_New_W", (char*)FileName) != FAILURE)
    {
        char VolSkipFlagYN;
        readerror = fscanf (stream, "%c\n", &VolSkipFlagYN);
        if (readerror != 1)
        {
            DR_Error("Could not read vol points skip flag in file %s! (VolDiag_Input_New_W)", FileName);
            goto RETURN;
        }
        *VolSkipFlag = toupper(VolSkipFlagYN) == 'Y' ? TRUE : FALSE;
    }
    else
        *VolSkipFlag = FALSE;

    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status == FAILURE)
    {
        DR_Error("VolDiag_Input_New_W: failed!");
    }
    

    return (status);

}  /* VolDiag_Input_New_W */

/*****  VolDiag_Input_W  ****************************************************/
/**
 *       Read volatility input
 */
int     VolDiag_Input_W (
            long    *   BaseDate    /** (O) Volatility data */
           ,int     *   VolUnit
           ,int     *   NbVol
           ,long    *   VolDate
           ,long    *   SwapSt
           ,long    *   SwapMat
           ,double  *   Vol
           ,char    *   VolType
           ,char const* FileName    /** (I) File name including extension */
       )
{
    int VolSkipFlag;       /* Dummy */

    return VolDiag_Input_New_W(BaseDate, VolUnit, NbVol, VolDate, SwapSt, SwapMat, Vol, VolType, &VolSkipFlag, FileName);
}

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
    MKTVOL_DATA*            mktvol_data     /** (O) Volatility data               */
    ,char const*            Index           /** (I) Index to calibrate             */
    ,T_CURVE const*   t_curve         /** (I) Term structure data            */
    ,char const*            BaseVolFile     /** (I) Base vol curve file            */
    ,char const*            SwapVolFile     /** (I) Swaption matrix file           */
    ,char const*            SwapVolFileSRM3 /** (I) Swaption matrix file from SRM3 */
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
        mktvol_data->BaseDate = t_curve->Today;
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
            mktvol_data->VolType[i] = 'C';
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
        DR_Error("%s not available.\n", SwapVolFile);
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
                mktvol_data->VolType[i] = 'S';
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
            mktvol_data->VolType[i] = 'S';

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




/*****  MktVol_Check_W  *****************************************************/
/**
*       Check validity of base volatility inputs.
*/
int     MktVol_Check_W (
        MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
        )
{

    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */


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
        DR_Error ("Nb of vols exceeds maximum of %d! (MktVol_Check_W)", MAXNBDATE);
        goto RETURN;
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
        if (Dateok(mktvol_data->VolDate[i]))
        {
            DR_Error("Incorrect format for volatility date [%ld] (MktVol_Check_W)!", mktvol_data->VolDate[i]);
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
        DR_Error("Incorrect format for base date [%ld] (MktVol_Check_W)!", mktvol_data->BaseDate);
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
            DR_Error("Benchmark swap matures %ld before swap start %ld (MktVol_Check_W)!", 
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


/*****  MktVol_PrintStructure  **********************************************/
/**
*       Print out all fields of the structure to provided stream 
*       to assist debugging
*/
int     MktVol_PrintStructure (
        FILE*        stream,
        MKTVOL_DATA  *mktvol_data
        )
{
    int i;
    fprintf(stream, "VolatilityCurve\n");
    fprintf(stream, "BaseDate %ld\n", mktvol_data->BaseDate);
    fprintf(stream, "NbVol %d\n", mktvol_data->NbVol);
    fprintf(stream, "NbCetVol %d\n", mktvol_data->NbCetVol);
    fprintf(stream, "VolDate      Vol       VolUsed\n");
    for (i = 0; i < mktvol_data->NbVol; i++)
    {
        fprintf(stream, "%ld   %lf   %d\n", 
                mktvol_data->VolDate[i],
                mktvol_data->Vol[i],
                mktvol_data->VolUsed[i]);
    }
    fprintf(stream, "CetNbIter %d\n", mktvol_data->CetNbIter);
    fprintf(stream, "CetVegaError %lf\n", mktvol_data->CetVegaError);

    fprintf(stream, "ForwardSwap\n");
    fprintf(stream, "Freq %c\n", mktvol_data->Freq);
    fprintf(stream, "DCC %c\n", mktvol_data->DCC);
    fprintf(stream, "SwapSt    SwapMat\n");
    for (i = 0; i < mktvol_data->NbVol; i++)
    {
        fprintf(stream, "%ld  %ld\n", 
                mktvol_data->SwapSt[i], 
                mktvol_data->SwapMat[i]);
    }

    fprintf(stream, "Model Parameters\n");
    fprintf(stream, "QLeft %lf\n", mktvol_data->QLeft);
    fprintf(stream, "QRight %lf\n", mktvol_data->QRight);
    fprintf(stream, "FwdShift %lf\n", mktvol_data->FwdShift);
    fprintf(stream, "Alpha   Beta   Rho\n");
    for (i = 0; i < 3; i++)
    {
        fprintf(stream, "%lf  %lf  %lf\n", 
                mktvol_data->Alpha[i],
                mktvol_data->Beta[i],
                mktvol_data->Rho[i]);
    }

    fprintf(stream, "Backbone Parameters\n");
    fprintf(stream, "Bbq %lf\n", mktvol_data->Bbq);
    fprintf(stream, "VolNorm %lf\n", mktvol_data->VolNorm);
    fprintf(stream, "VolLogn %lf\n", mktvol_data->VolLogn);

    fprintf(stream, "Calibration Flags\n");
    fprintf(stream, "SkipFlag %d\n", mktvol_data->SkipFlag);
    fprintf(stream, "CalibFlag %d\n", mktvol_data->CalibFlag);
    fprintf(stream, "FilterSpotVolFlag %d\n", mktvol_data->FilterSpotVolFlag);
    fprintf(stream, "SmoothingFlag %c\n", mktvol_data->SmoothingFlag);
    fprintf(stream, "TraceFlag %c\n", mktvol_data->TraceFlag);

    return SUCCESS;
}


/** 
*   Print probability in output file for wrapper.
*/
void    printProbabilityInFile(
        OPT_OUT_DATA *opt_out_data)
{
    FILE*    fproba;
    int      i,j;

    fproba = fopen ("probability.dat", "w");
    
    if (fproba == NULL)
    {
        DR_Error ("Could not open file probability! (main)");
        return;
    }

   
    for (i=0; i<ESL_NB_EVENT ; i++)
    {
        if (opt_out_data->prob_calc[i] == TRUE)
        {
            if (i == EXER_EVENT) 
            {
                /* Title */
                fprintf (fproba, "#\"ExerciseStatistics\"\n");
                fprintf (fproba, "# Exer Prob\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].TotalEventProb);
                fprintf (fproba, "# Exer Time \n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].ExpEventTime);
                fprintf (fproba, "# Exer Time s.d.\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].StdEventTime);
                fprintf (fproba, "# Fugit\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[EXER_EVENT].Fugit);
          

            }
	        else if (i == KO_EVENT)
            {
                /* Title */
                fprintf (fproba, "#\"KOStatistics\"\n");
                fprintf (fproba, "# KO Prob\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[KO_EVENT].TotalEventProb);
                fprintf (fproba, "# KO Time\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[KO_EVENT].ExpEventTime);
                fprintf (fproba, "# KO Time s.d.\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[KO_EVENT].StdEventTime);
                fprintf (fproba, "# Fugit\n");
                fprintf (fproba, "%-12.6f\n", opt_out_data->prob_out_data[KO_EVENT].Fugit);
               
            }

            
            fprintf (fproba, "# Cardinality probability dates \n");
            fprintf (fproba, "%d\n", opt_out_data->prob_out_data[i].Count);

            fprintf (fproba, "# Date \tProbability\n");

            for (j=0; j<opt_out_data->prob_out_data[i].Count; j++)
                fprintf (fproba, "%ld\t%8.6f\n", opt_out_data->prob_out_data[i].EventDate[j], 
                                opt_out_data->prob_out_data[i].EventProb[j]);
        }   
    }

    fclose (fproba);

    return;
}

#ifdef ESL_NEW_CURVE
void EslPrintZeroCurveAndDiscFactors(T_CURVE const* crv,
                                     int includeDiscFactors,
                                     FILE*          file)
{  
    int i;

 /* REVIEW */
    fprintf(file, "# Start date\n%ld\n", YMDDateFromIRDate(crv->baseDate));
    fprintf(file, "# Money Market basis (360 or 365)\n%d\n", crv->MMB);
    fprintf(file, "# Annual or semi-annual curve (A or S)\n%c\n", crv->SwapFreq);
    fprintf(file, "# Year basis for benchmark swaps (ACT, 365 or 360)\n%s\n", crv->SwapDCC);

    fprintf(file, "# Number of entries\n%d\n", crv->numItems - 1);
    if (includeDiscFactors)
        fprintf(file, "# Zero Dates (yyyymmdd), Rates (ACT/365F annual), and Disc Factors\n");
    else
        fprintf(file, "# Zero Dates (yyyymmdd) and Rates (ACT/365F annual)\n");
    for (i=0; i<crv->numItems; ++i)
    {
       double rate;

    /* REVIEW */
       if (crv->startDates[i] == crv->baseDate)
           continue;

       if (irxZeroRate(crv,
                       crv->startDates[i],
                       /*strcmp(crv->SwapDCC, "360") == 0 ? IRX_ACT_360 : strcmp(crv->SwapDCC, "ACT") == 0 ? IRX_ACT_ACT : IRX_ACT_365F,*/
                       IRX_ACT_365F,
                       IRX_ANNUAL_RATE,
                       &rate) != SUCCESS)
           break;

       if (includeDiscFactors)
            fprintf(file, "%ld %16.12f %.17lf\n", YMDDateFromIRDate(crv->startDates[i]), rate * 100.0, GetZeroPrice(crv->startDates[i], crv));
       else
            fprintf(file, "%ld %16.12f\n", YMDDateFromIRDate(crv->startDates[i]), rate * 100.0);
    }
}  

/* MAW 2.0 Compliant version */
void EslPrintToStringZeroCurve2(T_CURVE const* crv, char** curveString)
{
    static char routine[] = "EslPrintToStringZeroCurve2";
    int i;
    char* tempString=0;

    *curveString = (char*)calloc(sizeof(char), BIG_STRING);
    if (curveString == 0) {
        DR_Error("%s: could not allocate memory for string!", routine);
        goto done;
    }

    tempString = (char*)calloc(sizeof(char), SMALL_STRING);
    if (tempString == 0) {
        DR_Error("%s: could not allocate memory for string!", routine);
        goto done;
    }

 /* REVIEW */
    sprintf(*curveString, "# Start date\n%ld\n", YMDDateFromIRDate(crv->baseDate));
    sprintf(tempString, "# DCC\n%s\n", "ACT/365F");
    strcat(*curveString, tempString);
    sprintf(tempString, "# Freq\n%c\n", 'A');
    strcat(*curveString, tempString);
    sprintf(tempString, "# Number of entries\n%d\n", crv->numItems - 1);
    strcat(*curveString, tempString);

    sprintf(tempString, "# Zero Dates (yyyymmdd) and Rates (ACT/365F annual)\n");
    strcat(*curveString, tempString);
    for (i=0; i<crv->numItems; ++i)
    {
       double rate;

    /* REVIEW */
       if (crv->startDates[i] == crv->baseDate)
           continue;

       if (irxZeroRate(crv,
                       crv->startDates[i],
                       /*strcmp(crv->SwapDCC, "360") == 0 ? IRX_ACT_360 : strcmp(crv->SwapDCC, "ACT") == 0 ? IRX_ACT_ACT : IRX_ACT_365F,*/
                       IRX_ACT_365F,
                       IRX_ANNUAL_RATE,
                       &rate) != SUCCESS)
           break;

           sprintf(tempString, "%ld %16.12f\n", YMDDateFromIRDate(crv->startDates[i]), rate * 100.0);
           strcat(*curveString, tempString);
    }
done:
    if (tempString) {
        free(tempString);
    }
}

void EslPrintZeroCurve2(T_CURVE const* crv, FILE* file)
{  
    static char routine[] = "EslPrintZeroCurve2";
    char* curveString;

    curveString = 0;

    EslPrintToStringZeroCurve2(crv, &curveString);

    fprintf(file, curveString);

done:
    if (curveString) {
        free(curveString);
    }
}
#else
void EslPrintZeroCurveAndDiscFactors(T_CURVE const* crv,  
                                     int includeDiscFactors,
                                     FILE*          file)
{  
   int i;

   fprintf(file, "# Start date\n%ld\n", YMDDateFromIRDate(crv->ValueDate));
   fprintf(file, "# Money Market basis (360 or 365)\n%d\n", crv->MMB);
   fprintf(file, "# Annual or semi-annual curve (A or S)\n%c\n", crv->SwapFreq);
   fprintf(file, "# Year basis for benchmark swaps (ACT, 365 or 360)\n%s\n", crv->SwapDCC);

   fprintf(file, "# Number of entries\n%d\n", crv->NbZero);
    if (includeDiscFactors)
        fprintf(file, "# Zero Dates (yyyymmdd), Rates (ACT/365F annual), and Disc Factors\n");
    else
        fprintf(file, "# Zero Dates (yyyymmdd) and Rates (ACT/365F annual)\n");
   for (i=0; i<crv->NbZero; ++i)
   {
        if (includeDiscFactors)
            fprintf(file, "%ld %16.12f %.17lf\n", YMDDateFromIRDate(crv->ZeroDate[i]),crv->Zero[i] * 100.0, GetZeroPrice(crv->ZeroDate[i], crv));
        else
            fprintf(file, "%ld %16.12f\n", YMDDateFromIRDate(crv->ZeroDate[i]),crv->Zero[i] * 100.0);
   }
}  

/* MAW 2.0 Compliant Version */
void EslPrintToStringZeroCurve2(T_CURVE const* crv, char** curveString)
{
    static char routine[] = "EslPrintToStringZeroCurve2";
    int i = 0;
    char* tempString=0;

    *curveString = (char *)calloc(sizeof(char), BIG_STRING);
    if (*curveString == 0) {
        DR_Error("%s: could not allocate memory for string!", routine);
        goto done;
    }

    tempString = (char *)calloc(sizeof(char), SMALL_STRING);
    if (tempString == 0) {
        DR_Error("%s: could not allocate memory for string!", routine);
        goto done;
    }

    sprintf(*curveString, 
            "# Start date\n%ld\n", 
            YMDDateFromIRDate(crv->ValueDate));
    sprintf(tempString, "# DCC\n%s\n", "ACT/365F");
    strcat(*curveString, tempString);
    sprintf(tempString, "# Freq\n%c\n", 'A');
    strcat(*curveString, tempString);
    sprintf(tempString, "# Number of entries\n%d\n", crv->NbZero);
    strcat(*curveString, tempString);
    sprintf(tempString, "# Zero Dates (yyyymmdd) and Rates (ACT/365F annual)\n");
    strcat(*curveString, tempString);

    for (i=0; i<crv->NbZero; ++i)
    {
        sprintf(tempString, "%ld %16.12f\n", YMDDateFromIRDate(crv->ZeroDate[i]), crv->Zero[i] * 100.0);
        strcat(*curveString, tempString);
    }

done:
    if (tempString) {
        free(tempString);
    }
}

void EslPrintZeroCurve2(T_CURVE const* crv, FILE* file)
{  
    char* curveString;

    curveString = 0;

    EslPrintToStringZeroCurve2(crv, &curveString);

    fprintf(file, curveString);

/* done: */
    if (curveString) {
        free(curveString);
    }
}  
#endif /* ifdef ESL_NEW_CURVE else */

void EslPrintZeroCurve(T_CURVE const* crv,  
                       FILE*          file)
{
    EslPrintZeroCurveAndDiscFactors(crv, FALSE, file);
}

/*****  MktVol_Input_New_W  ************************************************/
/**
        Utility routine to read a series of option expiries, index maturities 
        and volatilities from ir_voldiag_i.dat
*/
int MktVol_Input_New_W (
    MKTVOL_DATA*       mktvol_data,    /** (O) Volatility data          */
    T_CURVE *          t_curve,        /** (I) Term structure data      */
    char const*        VolFileSRM3)    /** (I) Swaption matrix file     */
    
{
  /*char     IndexL[MAXINDEX]; */   /* Local copy of index */
    double   temp1, temp2;
    int     nbMthPerPeriod;
 
    int      i;
    int      status = FAILURE;   /* Error status = FAILURE initially         */

    /**********************/
    /*  Read Vol diagonal */
    /**********************/
    if (VolDiag_Input_New_W (
                        &(mktvol_data->BaseDate),
                        &(mktvol_data->VolUnit),
                        &(mktvol_data->NbVol),
                        mktvol_data->VolDate,
                        mktvol_data->SwapSt,
                        mktvol_data->SwapMat,
                        mktvol_data->Vol,
                        mktvol_data->VolType,
                        &mktvol_data->SkipFlag,
                        VolFileSRM3) == FAILURE) goto RETURN;
    /* Calib is off either CalibFlag from deal = 'N' or NbVol from env = 0 */ 
    if (mktvol_data->NbVol < 0) 
    {   
        DR_Error("Nb Vol has to be >= 0!");
        goto RETURN;    
    }
    else if (mktvol_data->NbVol == 0)
    {
        mktvol_data->BaseDate = t_curve[0].Today;
        mktvol_data->CalibFlag = FALSE;
        mktvol_data->SkipFlag = FALSE;
        mktvol_data->VolUnit = -99999;
        mktvol_data->NbVol = -99999;
        mktvol_data->Freq = 'z';
        mktvol_data->DCC  = 'z';
    }
    else
    {
        mktvol_data->CalibFlag = TRUE;
        /* Fill vector volUsed */
        for(i = 0; i < mktvol_data->NbVol; i++)
        {   
            mktvol_data->VolUsed[i] = TRUE; 
        }

        /* if vol type is cash use MM conventions */
        if (mktvol_data->VolType[0] == 'C')
        {
            
            temp1 = (double) Daysact(mktvol_data->SwapSt[0], mktvol_data->SwapMat[0]);
            temp2 = temp1/30.0;
            nbMthPerPeriod = (int) ceil(temp2 - 0.5);

            if( (nbMthPerPeriod == 1) )
            {
                mktvol_data->Freq = 'M';
            }
            else if (nbMthPerPeriod == 3)
            {
                mktvol_data->Freq = 'Q';
            }
            else if (nbMthPerPeriod == 6)
            {
                mktvol_data->Freq = 'S';
            }
         
            if (t_curve[0].MMB == 360)
            {
                mktvol_data->DCC= '0';
            }
            else if (t_curve[0].MMB == 365)
            {
                mktvol_data->DCC = '5';
            }
        }
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
       mktvol_data->SwapSt[i]  = mktvol_data->VolDate[i];

    
       
    /* check validity of vol inputs */
    if (MktVol_Check_New_W (mktvol_data) == FAILURE)
    {
        goto RETURN;
    }
      

    status = SUCCESS;

    RETURN:

    return (status);

}  /* MktVol_Input_New_W */


/*****  MktVol_Check_New_W  *****************************************************/
/**
*       Check validity of volatility inputs.
*/
int     MktVol_Check_New_W (
        MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
        )

{
    int     status = FAILURE;
    char    ErrorMsg[MAXBUFF];
    long    i;

    if (mktvol_data == NULL) goto RETURN;

    /* check filter spot vol flag */
    if ( (mktvol_data->FilterSpotVolFlag != TRUE) &&
         (mktvol_data->FilterSpotVolFlag != FALSE) )
    {
        DR_Error("FilterSpotVol flag must be TRUE or FALSE (MktVol_Check_W))");
        goto RETURN;
    }

    /* check smoothing flag */
    if ( (mktvol_data->SmoothingFlag != 'Y') &&
         (mktvol_data->SmoothingFlag != 'N') )
    {
        DR_Error("Cet smoothing flag must be 'Y' or 'N' (MktVol_Check_W))");
        goto RETURN;
    }
    /* check Calib Flag */
    if ((mktvol_data->CalibFlag != TRUE) &&
        (mktvol_data->CalibFlag != FALSE))
    {
        DR_Error("invalid Calib Flag! (MktVolCheck_New)");
        goto RETURN;
    }

     /* Nothing else to check */
    if (mktvol_data->CalibFlag == FALSE)
    {
        return (SUCCESS);
    }
    
    if ((mktvol_data->DCC != '0') &&
        (mktvol_data->DCC != '3') &&
        (mktvol_data->DCC != '5'))
    {
        sprintf (ErrorMsg, "dcc has to be '0' or"
                " '3' or '5'! (MktVolCheck_New)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (   (mktvol_data->Freq != 'A')
        && (mktvol_data->Freq != 'S')
        && (mktvol_data->Freq != 'Q')
        && (mktvol_data->Freq != 'M'))
    {
        sprintf (ErrorMsg, "swap freq has to be 'A' or"
                " 'S' or 'Q' or 'M'! (MktVolCheck_New)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Check vol inputs */
    if (mktvol_data->CalibFlag == TRUE)
    {
         if (VolDiag_Check_W (mktvol_data->BaseDate,        /* Volatility data */
                              mktvol_data->VolUnit,
                              mktvol_data->NbVol,
                              mktvol_data->VolDate,
                              mktvol_data->SwapSt,
                              mktvol_data->SwapMat,
                              mktvol_data->Vol,
                              mktvol_data->VolType) == FAILURE) goto RETURN;

         /* NbVol = 0 is not allowed for IR */
         if (mktvol_data->NbVol == 0)
         {
            sprintf (ErrorMsg, "Nb vol can't be zero! (MktVolCheck_New)");
            DR_Error(ErrorMsg);
            goto RETURN;
         }
         
         /* Check VolType[i] are all 'C' or all 'S' */
         /* Not done inside VolDiag_Check_W because other assets may
            allow combination of 'C' and 'S' */
         if ((mktvol_data->VolType[0] != 'C') &&
             (mktvol_data->VolType[0] != 'S'))
         {
            sprintf (ErrorMsg, "vol type has to be 'C' or 'S'! (MktVolCheck_New)");
            DR_Error(ErrorMsg);
            goto RETURN;
         }

         for (i=1; i<mktvol_data->NbVol; i++)
         {
             if (mktvol_data->VolType[i] != mktvol_data->VolType[0])
             {
                sprintf (ErrorMsg, "vol types have to be equal! (MktVolCheck_New)");
                DR_Error(ErrorMsg);
                goto RETURN;
             }
         }

    }
        

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
    {
        DR_Error("MktVolCheck_New: failed!");
    }
    return (status);
}


/*****  VolDiag_Check_W  *****************************************************/
/*
*       Check validity of volatility inputs.
*/

int     VolDiag_Check_W (long    BaseDate,        /* Volatility data */
                         int     VolUnit,
                         int     NbVol,
                         long    *VolDate,
                         long    *SwapSt,
                         long    *SwapMat,
                         double  *Vol,
                         char    *VolType)    
{
    int     i;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];

    /* Basic Check */
    if ((VolDate == NULL) ||
        (SwapSt  == NULL) ||
        (SwapMat == NULL) ||
        (Vol     == NULL) ||
        (VolType == NULL)) goto RETURN;


    if (Dateok(BaseDate))
    {
        sprintf (ErrorMsg, "Incorrect format for start date! (VolDiag_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if ((VolUnit != 0) &&
        (VolUnit != 1))
    {
        sprintf (ErrorMsg, "Volatility unit must be 0 or 1! (VolDiag_Check_W)");
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    if ((NbVol > MAXNBDATE) ||
        (NbVol < 0))
    {        
        sprintf (ErrorMsg, "Nb of vols must be >= 0 and <= %d! (VolDiag_Check_W)", MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    for (i=0; i<NbVol; i++)
    {
        if ((Dateok(VolDate[i]))||
            (Dateok(SwapSt [i]))||
            (Dateok(SwapMat[i]))) 
        {
            sprintf (ErrorMsg, "incorrect date format #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (SwapMat[i] <= SwapSt[i])
        {
            sprintf (ErrorMsg, "rate end <= rate start #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (VolDate[i] >  SwapSt[i])
        {
            sprintf (ErrorMsg, "expiration date > rate start #%d! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (Vol[i] < 0.0)
        {
            sprintf (ErrorMsg, "vol #%d < 0! (VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if ((VolType[i] != 'C') &&
            (VolType[i] != 'S'))
        {
            sprintf (ErrorMsg, "Volatility type #%d has to be 'C' or 'S' ! "
                               "(VolDiag_Check_W)", i+1);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }


    for (i = 1; i < NbVol; i++)
    {
        if (VolDate[i] <= VolDate[i-1])
        {
            sprintf (ErrorMsg, "Volatilty dates must be entered in ascending order! (VolDiag_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    if (NbVol != 0)
    {
        if (BaseDate >= VolDate[0])
        {
            sprintf (ErrorMsg, "Base date must be before vol date! (VolDiag_Check_W)");
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }

    status = SUCCESS;
        
    RETURN:

    if (status == FAILURE)
    {
        DR_Error("VolDiag_Check_W: failed!");
    }
        
    return (status);

}  /* VolDiag_Check_W */



/*****  MktInfo  *****************************************************/
/*
*       Read inputs from ir_info
*/
int MktInfo(
         MKTVOL_DATA*       mktvol_data,       /** (O) Volatility data         */
         int*               diffusionCurveIdx, /** (0) Index of diffused curve */
         int*               discountCurveIdx,  /** (0) Index of discount curve */
         T_CURVE *          t_curve,           /** (I) Term structure data     */
         char const*        IRInfoFileName)    /** (I) IR info file            */
{   
    char    SwapDCC[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */
    FILE    *stream = NULL;
    char    ErrorMsg[MAXNBDATE];
    int     MMBinput, readerror;

    /* Open and read the IR info file */
    stream = fopen (IRInfoFileName, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (MktInfo)", IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* read Money Market Basis */
    if (FindAndSkipComLine (stream, "money market basis", "MktInfo", 
        IRInfoFileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(MMBinput));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read money market basis in %s!(MktInfo)",
           IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    
    /* define money market basis same for all curves */
    if ((MMBinput == 360) || (MMBinput == 365))
    {
        t_curve[0].MMB = MMBinput;
        t_curve[1].MMB = MMBinput;
        t_curve[2].MMB = MMBinput;
    }
    else
    {
        sprintf (ErrorMsg, "MMB in %s has to be 360 or 365! (MktInfo)", IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
     /* read Swap Day Count Convention */
    if (FindAndSkipComLine (stream, "swap day count convention", "MktInfo", IRInfoFileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%s \n", SwapDCC);
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read swap dcc in file %s! (MktInfo)", IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }                                
    
    if (strstr(SwapDCC, "30/360") != NULL)
    {
        mktvol_data->DCC = '3'; 
        strcpy (t_curve[0].SwapDCC, "ACT");
        strcpy (t_curve[1].SwapDCC, "ACT");
        strcpy (t_curve[2].SwapDCC, "ACT");

    }
    else if (strstr(SwapDCC, "ACT/365F") != NULL)
    {
        mktvol_data->DCC = '5';
        strcpy (t_curve[0].SwapDCC, "365");
        strcpy (t_curve[1].SwapDCC, "365");
        strcpy (t_curve[2].SwapDCC, "365");
    }
    else if (strstr(SwapDCC, "ACT/360") != NULL)
    {
        mktvol_data->DCC = '0';
        strcpy (t_curve[0].SwapDCC, "360");
        strcpy (t_curve[1].SwapDCC, "360");
        strcpy (t_curve[2].SwapDCC, "360");
    }
    else
    {
        sprintf (ErrorMsg, "swap dcc in file %s has to be '30/360' or"
                " 'ACT/365F' or 'ACT/360'! (MktInfo)", IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* read swap frequency */
    if (FindAndSkipComLine (stream, "swap frequency", "MktInfo", IRInfoFileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%c \n", &(mktvol_data->Freq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read swap frequency in %s! (MktInfo)", IRInfoFileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    /* define zero curves benchmark swap frequencies */
    t_curve[0].SwapFreq = mktvol_data->Freq;
    t_curve[1].SwapFreq = mktvol_data->Freq;
    t_curve[2].SwapFreq = mktvol_data->Freq;

    if (FindAndSkipComLine_2(FIND_SKIP_SILENT,
                             stream,
                             "Diffusion curve index",
                             "MktInfo",
                             (char*)IRInfoFileName) != FAILURE)
    {
        readerror = fscanf (stream, "%d\n", diffusionCurveIdx);
        if (readerror != 1)
        {
            DR_Error("Could not find diffusion curve index in file %s! (MktInfo).\n", IRInfoFileName);
            goto RETURN;
        }

        if (FindAndSkipComLine(stream, "Discount curve index", "MktInfo", (char*)IRInfoFileName) == FAILURE)
            goto RETURN;

        readerror = fscanf (stream, "%d\n", discountCurveIdx);
        if (readerror != 1)
        {
            DR_Error("Could not find discount curve index in file %s! (MktInfo).\n", IRInfoFileName);
            goto RETURN;
        }
    }
    else
    {
        *diffusionCurveIdx = 0;
        *discountCurveIdx = 1;
    }

    status = SUCCESS;
    
  RETURN:

    if (stream != NULL)
        fclose (stream);
 
    return (status);

} /* MktInfo */    



/*****  SummaryInfo  *****************************************************/
/*
*       Read inputs from summary.dat
*/
int SummaryInfo(
         T_CURVE *          t_curve,         /** (I) Term structure data    */
         char const*        FileName)     /** (I) Swaption matrix file   */
{
    char ErrorMsg[MAXBUFF];
    long Today, NbIr;
    int  readerror;
    FILE     *stream = NULL;
    int status = FAILURE;
    


    stream = fopen (FileName, "r");

    /* if file does not exist nothing to do */
    if (stream == NULL)
        return (SUCCESS);

    if (FindAndSkipSectionLine (
                             1, /* skip mode */
                             stream, 
                            "Environment section", 
                            "SummaryInfo", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipComLine(
                      stream, "today", "SummaryInfo", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream,
                        "%ld \n", 
                        &(Today));
    if (readerror != 1)
    {      
        sprintf (ErrorMsg, "Could not today in file %s! (SummaryInfo)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    Today = IRDateFromYMDDate(Today);

   
    if (Dateok(Today))
    {
        sprintf(ErrorMsg, "Incorrect format for today date in file %s! (SummaryInfo)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }


    /* set today for zero curves to the input value */
    t_curve[0].Today = Today;
    t_curve[1].Today = Today;
    t_curve[2].Today = Today;

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End section", 
                            "SummaryInfo", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    if (FindAndSkipSectionLine (
                             1, /* skip mode */
                             stream, 
                            "IR Section", 
                            "Summary Info", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipComLine(
                      stream, "NbIR", "Summary Info", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream,
                        "%ld \n", 
                        &NbIr);

    if (readerror != 1)
    {      
        sprintf (ErrorMsg, "Could not find NbIR in file %s! (SummaryInfo)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End section", 
                            "SummaryInfo", 
                            FileName) == FAILURE)
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

} /* SummaryInfo */

int ModelInfo(const char* fileName, const char* expectedEngineName, char** modelChoiceString)
{
    int status = FAILURE;
    const char* routine = "ModelInfo";

    int     readerror;
    char    ModelChoiceString[MAXBUFF];
    char    EngineChoiceString[MAXBUFF];

    FILE*   fp = NULL;

    if (modelChoiceString == NULL)
    {
        DR_Error("NULL modelChoiceString parameter for routine: %s.\n", routine);
        goto RETURN;
    }
    *modelChoiceString = NULL;

    /* Read model choice */
    fp = fopen (fileName, "r");
    if (fp != NULL)
    {
         if (FindAndSkipComLine (fp, "Model File", routine, fileName) == FAILURE)
             goto RETURN;

        if (FindAndSkipComLine(fp, "Engine Choice", routine, fileName) == FAILURE)
            goto RETURN;

        readerror = fscanf (fp, "%s\n", EngineChoiceString);
        if (readerror != 1)
        {
            DR_Error ("Could not read engine choice in file %s.\n", fileName);
            goto RETURN;
        }

        if (expectedEngineName != NULL && irxStrcmpi(EngineChoiceString, expectedEngineName) != 0)
        {
            DR_Error ("Engine choice  should be %s, not %s", expectedEngineName, EngineChoiceString);
            goto RETURN;
        }

        if (FindAndSkipComLine (fp, "Model Choice", routine, fileName) == FAILURE)
        {
            goto RETURN;
        }

        readerror = fscanf (fp, "%s\n", ModelChoiceString);
        if (readerror != 1)
        {
            DR_Error ("Could not read model choice in file %s.", fileName);
            goto RETURN;
        }

        *modelChoiceString = irxStringToUpper(ModelChoiceString);
    }
    else /* backward compatibility -- model information file not found */
    {
        *modelChoiceString = strdup("ORIGINAL");
    }

    if (*modelChoiceString != NULL)
        status = SUCCESS;

RETURN:
    if (fp != NULL)
        fclose(fp);

    return status;
}
