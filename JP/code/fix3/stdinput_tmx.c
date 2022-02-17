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
#include "fix123head.h"




/*****  BaseSmile_Input_W  ***************************************************/
/*
*       Read base volatility smile input and check validity of input.
*       Read whole file as smile parameters between idxLo and idxHi
*/ 
int     BaseSmile_Input_W (
            MKTSMILE_DATA  *ms,             /* (O) Volatility data           */
            int            idxLo,           /* (I) Lower index               */
            int            idxHi,           /* (I) Higher index              */
            long           BaseDate,        /* (I) Base date                 */  
            char           *File)           /* (I) File name                 */
{

    char    Code[MAXBUFF] = "BaseSmile_Input";
    int     i, j, s;
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;


    /* Open the base volatility smile data file */
    stream = fopen (File, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! %s", File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "frequency", Code, File) == FAILURE) 
        goto RETURN;
    
    readerror = fscanf (stream, "%c \n", &(ms->Freq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read frequency in %s! %s", 
                 File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (FindAndSkipComLine (stream, "nb of expiries", Code, File) == FAILURE)
        goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(ms->NbExpiry));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of expiries in %s! %s", 
                 File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (ms->NbExpiry > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of vols in %s exceeds maximum of %d! %s", 
                 File, MAXNBDATE, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    ms->NbFwdMat = 1;
    switch (ms->Freq)
    {
        case ('A'): ms->FwdMat[0] = 12; break;
        case ('S'): ms->FwdMat[0] = 6;  break;
        case ('Q'): ms->FwdMat[0] = 3;  break;
        case ('M'): ms->FwdMat[0] = 1;  break;
        default:    ms->FwdMat[0] = -999;
    }

    for (s = 0; s <= idxHi-idxLo; s++) 
    {
        if (FindAndSkipComLine (stream, "vol info", Code, File) == FAILURE)
            goto RETURN;

        for (i = 0; i < ms->NbExpiry; i++)
        {
            readerror = fscanf (stream, "%ld \t%lf\n", 
                                &(ms->VolDate[i]),
                                &(ms->Vol[s+idxLo][i][0]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not read vol date & rate #%d "
                                   "in %s! %s", i+1, File, Code);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }

    }

    /* basevol file case */
    if (idxLo == 0)
    {
        /* Preprocess vol formats -- vol is 0th parameter! */
        for (i = 0; i < ms->NbExpiry; i++)
        {
            ms->Vol[0][i][0] /= 100.;
        }
    }

    /* Eliminate dates falling before base date */
    j = 0;
    while (ms->VolDate[j] <= BaseDate)
        j++;

    ms->NbExpiry -= j;

    for (i = 0; i < ms->NbExpiry; i++)
    {
        ms->VolDate[i] = ms->VolDate[i+j];
        for (s = idxLo; s <= idxHi ; s++) ms->Vol[s][i][0] = ms->Vol[s][i+j][0];
    }

    status = SUCCESS;
        
    RETURN:
        
    if (status != SUCCESS)
    {        
        DR_Error ("%s failed!", Code);
    }

    if (stream != NULL)
    {
        fclose (stream);
    }
        
    return (status);

}  /* BaseSmile_Input_W */


/*****  SwapSmile_Input_W  ***************************************************/
/*
*       Read swaption volatility input and check validity of input.
*       Read whole file as smile parameters between idxLo and idxHi 
*/
int     SwapSmile_Input_W (
            MKTSMILE_DATA  *ms,             /* (O) Volatility data           */
            int            idxLo,           /* (I) Lower index               */
            int            idxHi,           /* (I) Higher index              */
            long           BaseDate,        /* (I) Base date                 */  
            int            isTMXStylFile,   /* (I) Is the file TMX format    */
            char           *File)           /* (I) File name                 */
{

    char    Code[MAXBUFF] = "SwapSmile_Input";
    int     i, j, s, k;
    int     readerror;          /* Reading error status         */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    /* Open the swaption volatility smile data file */
    stream = fopen (File, "r");

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! %s", File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;        
    }
        
    if (FindAndSkipComLine (stream, "nb of expiries", Code, File) == FAILURE)
        goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(ms->NbExpiry));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of expiry in  %s! %s", 
                 File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
        
    if (ms->NbExpiry > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of expiry in %s exceeds maximum of %d! %s",
                 File, MAXNBDATE, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "number of tenors", Code, File) == FAILURE)
        goto RETURN;

    readerror = fscanf (stream, "%ld \n", &(ms->NbFwdMat));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of swap tenors in %s! %s", 
                 File, Code);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
      
    for (s = 0; s <= idxHi-idxLo; s++)
    {
        if (FindAndSkipComLine (stream, "vol matrix", Code, File) == FAILURE)
        goto RETURN;

        for (j = 0; j < ms->NbFwdMat; j++)
        {
            readerror = fscanf (stream, "\t%ld", &(ms->FwdMat[j]));
            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not read mat in %s! %s", File, Code);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        
            ms->FwdMat[j] *= 12;    /* Convert to months */
        }
            
        for (i = 0; i < ms->NbExpiry; i++)
        {
            readerror = fscanf (stream, "%ld", &(ms->Expiry[i]));
            if (readerror != 1)
            {        
                sprintf (ErrorMsg,"Could not read expiry in %s! %s",File,Code);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if ((ms->Expiry[i] < 10000)&&(isTMXStylFile == TRUE))
            {
                sprintf(ErrorMsg,"Wrong date format in the expiry in file %s! %s",File,Code);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            for (j = 0; j < ms->NbFwdMat; j++)
            {
                readerror = fscanf (stream, "\t%lf", &(ms->Vol[s+idxLo][i][j]));
                if (readerror != 1)
                {        
                    sprintf (ErrorMsg,"Could not read vol in %s! %s",File,Code);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }
            readerror = fscanf (stream, "\n");
        }
    }
        
    /* swapvol file case  */
    if (idxLo == 0)
    {
        /* Preprocess vol formats -- vol is 0th parameter! */
        for (i = 0; i < ms->NbExpiry; i++)
            for (j = 0; j < ms->NbFwdMat; j++)
                ms->Vol[0][i][j] /= 100.;

    } /* if idxLo == 0 */

    /* Create vol dates */
    for (i = 0; i < ms->NbExpiry; i++)
    {
        if (isTMXStylFile == TRUE)
        {
            ms->VolDate[i] =  ms->Expiry[i];
        }
        else
    {
            ms->VolDate[i] = Nxtmth (BaseDate, ms->Expiry[i], 1L);
        }
    }

    /* Eliminate dates falling before base date */
    j = 0;

    while (ms->VolDate[j] <= BaseDate)
        j++;

    ms->NbExpiry -= j;

    for (i = 0; i < ms->NbExpiry; i++)
    {
        ms->VolDate[i] = ms->VolDate[i+j];
        for (s = idxLo; s <= idxHi ; s++) 
            for (k = 0; k < ms->NbFwdMat; k++)
                ms->Vol[s][i][k] = ms->Vol[s][i+j][k];
    }

    status = SUCCESS;
        
    RETURN:
           
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error ("%s failed!", Code);
    }

    return (status);

}  /* SwapSmile_Input_W */




/*****  smileinterp *********************************************************/
/*      
*       Smile parameter interpolation; size and dates of smile inputs do not
*       have to agree with those of the vol data structure. Assumes that the 
*       benchmark swaps in MKTVOL_DATA have already been selected.
*/
int smileinterp (
               MKTVOL_DATA      *mvd,       /* (I/O) Vol data   */
               MKTSMILE_DATA    *msSmile)   /* (I)   Smile data */
{
    double  expiry, expFrac, mat, matFrac;
    double  expiry1, expiry2;
    int     nbRows = msSmile->NbExpiry;
    int     nbCols = msSmile->NbFwdMat;
    int     row    = 0;
    int     col    = 0;
    int     i,s;
    int     status = FAILURE;

    
    /* special case: only one expiration dats in msSmile */
    if (nbRows == 1)
    {
        for (i = 0; i < mvd->NbVol; i++)
        {
            if (nbCols == 1)
            {
                for (s = 1; s < NBVOLPARS; s++) 
                {
                    mvd->Smile[s][i] = msSmile->Vol[s][0][0];
                }
            }
            else
            {
                col = 0;
                while ((col < nbCols - 2) &&
                       (mvd->SwapMatMos[i] > msSmile->FwdMat[col+1]))
                {
                    col++;
                }
                mat = (double) COLLAR (mvd->SwapMatMos[i],
                               msSmile->FwdMat[nbCols-1], msSmile->FwdMat[0]);
                matFrac = (mat - msSmile->FwdMat[col]) /
                          (msSmile->FwdMat[col+1] - msSmile->FwdMat[col]);

                for (s = 1; s < NBVOLPARS; s++) 
                {
                    mvd->Smile[s][i] = (1.-matFrac) * msSmile->Vol[s][0][col] +
                                       matFrac      * msSmile->Vol[s][0][col+1];
                }
            }
        }
    }
    else /* interpolate over expiry as well */
    {
    
        for (i = 0; i < mvd->NbVol; i++)
        {
            row = 0;
            while ((row < nbRows-2) &&
                   (mvd->VolDate[i] > msSmile->VolDate[row+1]))
            {
                row++;
            }

            /* extension is flat outside the msSmile vol dates */
            expiry = (double) Daysact (mvd->BaseDate, COLLAR (mvd->VolDate[i],
                         msSmile->VolDate[nbRows-1],msSmile->VolDate[0]))/365.;
            expiry1 = (double) Daysact (mvd->BaseDate, msSmile->VolDate[row])
                            / 365.;
            expiry2 = (double) Daysact (mvd->BaseDate, msSmile->VolDate[row+1])
                            / 365.;
            expFrac = (expiry - expiry1) / (expiry2 - expiry1);

            if (nbCols == 1)
            {
                for (s = 1; s < NBVOLPARS; s++) 
                {
                    mvd->Smile[s][i] = (1.-expFrac) * msSmile->Vol[s][row][0] +
                                       expFrac      * msSmile->Vol[s][row+1][0];
                }
            }
            else
            {
                col = 0;
                while ((col < nbCols - 2) &&
                       (mvd->SwapMatMos[i] > msSmile->FwdMat[col+1]))
                {
                    col++;
                }
                mat = (double) COLLAR (mvd->SwapMatMos[i],
                               msSmile->FwdMat[nbCols-1], msSmile->FwdMat[0]);
                matFrac = (mat - msSmile->FwdMat[col]) /
                          (msSmile->FwdMat[col+1] - msSmile->FwdMat[col]);

                for (s = 1; s < NBVOLPARS; s++) 
                {
                    mvd->Smile[s][i] = 
                        (1.-expFrac)*(1.-matFrac)*msSmile->Vol[s][row][col]   +
                        (1.-expFrac)*matFrac     *msSmile->Vol[s][row][col+1] +
                        expFrac     *(1.-matFrac)*msSmile->Vol[s][row+1][col] +
                        expFrac     *matFrac     *msSmile->Vol[s][row+1][col+1];
                }
            }
        }
    }

    status = SUCCESS;

    return (status);
}


/*****  MktSmile_Input_W  ***************************************************/
/*
*       Read flat smile input.
*/
int     MktSmile_Input_W (
            MKTSMILE_DATA  *ms,             /* (O) Volatility data           */
            char           OWS[MAXBUFF],    /* (I) Overwrite strings         */
            long           BaseDate,        /* (I) Base Date                 */
            char           BoS,             /* (I) Base or swap smile        */
            int            isTMXStylFile,   /* (I) Tmx style file            */
            char           *File)           /* (I) File name                 */
{

    char    Code[MAXBUFF] = "MktSmile_Input";
    double  OWVol[NBVOLPARS];   /* Owevright vol smile              */
    int     i, j, s;
    int     readerror;          /* Reading error status             */
    int     validOWS = FALSE;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;

    static  double LognSmile[] = {-999,0.,0.,0.,0.,0.25,0.,0.75,0.};
    static  double NormSmile[] = {-999,1.,0.,0.,0.,0.25,1.,0.75,1.};


    /* Open the swaption volatility smile data file */
    stream = fopen (File, "r");

    if (stream == NULL)
    {
        /* File does not exist */
        sprintf (ErrorMsg, "Could not find file %s", File);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    else
    {
        if (BoS == 'B')
        {
            if (BaseSmile_Input_W (ms, 1, NBVOLPARS-1,
                                   BaseDate,
                                   File) == FAILURE)
            {
                DR_Error ("Incorrect format for base smile file %s", File);
                goto RETURN;
            }
        }
        else if (BoS == 'S')
        {
            if (SwapSmile_Input_W (ms, 1, NBVOLPARS-1,
                                   BaseDate,
                                   isTMXStylFile,
                                   File) == FAILURE)
            {
                DR_Error ("Incorrect format for swap smile file %s", File);
                goto RETURN;
            }
        }
        else 
        {
            DR_Error ("Smile input type must be 'B'ase or 'S'wap");
            goto RETURN;
        }

        /* if there is overwrite, use it */
        if (strchr (OWS, 'L') != NULL) 
        {
            for (s=1;s<NBVOLPARS;s++) OWVol[s] = LognSmile[s];
            validOWS = TRUE;
        }
        else if (strchr (OWS, 'N') != NULL)
        {
            for (s=1;s<NBVOLPARS;s++) OWVol[s] = NormSmile[s];
            validOWS = TRUE;
        }
        else if (strstr(OWS, "nil") == NULL)
        {
            readerror = sscanf (OWS,
                            "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
                            &(OWVol[1]),
                            &(OWVol[2]),
                            &(OWVol[3]),
                            &(OWVol[4]),
                            &(OWVol[5]),
                            &(OWVol[6]),
                            &(OWVol[7]),
                            &(OWVol[8]),
                            &(OWVol[9]));

            if (readerror != NBVOLPARS-1)
            {      
                sprintf (ErrorMsg, "Could not read smile owrt! %s",
                         Code);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
            validOWS = TRUE;
        }

        if (validOWS)
        {
            /* Copy to ms structure */
            for (i = 0; i < ms->NbExpiry; i++)
            {
                for(j = 0; j < ms->NbFwdMat; j++)
                {
                    ms->Vol[0][i][j] = -999;

                    for (s=1;s<NBVOLPARS;s++) ms->Vol[s][i][j] = OWVol[s];
                }
            }
        }

    } /* if file exists */

    status = SUCCESS;
        
    RETURN:
           
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (status != SUCCESS)
    {        
        DR_Error ("%s failed!", Code);
    }

    return (status);

}  /* MktSmile_Input_W */



#define NILNBVOL 15
/*****  Fix3_MktVol_Input_W_Tmx ***********************************************/
/*
*       Utility routine converting an index name to a series of option expi-
*       ries, index maturities and volatilities taken from either the base
*       volatility curve or the swaption matrix.
*/
int     Fix3_MktVol_Input_W_Tmx (
            MKTVOL_DATA     *mvd,             /* (O) Volatility data          */
            char            *Index,           /* (I) Index to calibrate       */
            T_CURVE const   *t_curve,         /* (I) Term structure data      */
            char            OWS[MAXBUFF],     /* (I) Overwrite strings        */
            char            *BaseVolFile,     /* (I) Base vol curve file      */
            char            *BaseSmlFile,     /* (I) Base smile file          */
            char            *SwapVolFile,     /* (I) Swaption matrix file     */
            char            *TMXSwapVolFile,  /* (I) TMX Swaption matrix file */
            char            *SwapSmlFile)     /* (I) Swaption smile file      */
{

    char    *IndexL = NULL;       
    char    *StrIdx = NULL;
    long    IdxMat;             /* Maturity of the index (final or forward) */
    long    Mat;
    long    SmlExp;
    int     NbRows=0;           /* Swaption volatility matrix size */
    int     NbCols=0; 
    int     ILiqBmk;
    int     SmlLiqFound = FALSE;
    int     isTMXStylFile;      /* test if swapvolx.dat exist      */
    MKTSMILE_DATA msVol;        /* For ATM vol parameters */
    MKTSMILE_DATA msSmile;      /* For smile parameters   */
    char    *LocalSwapVolFile;  /* used swaption matrix file     */

    int     i, j, matIdx, res;
    int     status = FAILURE;   /* Error status = FAILURE initially         */
    struct stat    testFile;

    /* Save index */
    IndexL = Index;
 
    /* test if swapvolx.dat exist */
    res = stat (TMXSwapVolFile, &testFile);

    if (res == 0)
    {
         isTMXStylFile = TRUE;
         LocalSwapVolFile = TMXSwapVolFile;
    }
    else
    {
         isTMXStylFile = FALSE;
         LocalSwapVolFile = SwapVolFile;
    }

    /* Clean Up */
    msVol.NbExpiry = 0;
    msVol.NbFwdMat = 0;
    msVol.Freq     = 'z';

    /* No calibration case */
    if (strstr (IndexL, "nil") != NULL)
    {
        DR_Error ("MktVol_Input_W: A Calibration Index must be specified!");
        goto RETURN;
    }    

    mvd->CalibFlag      = TRUE;
    mvd->CalibSmileFlag = TRUE;
    mvd->SkipFlag       = FALSE;


    /* Search for * character in calibration index name */
    StrIdx = strchr(IndexL, '*');

    if (StrIdx != NULL)
    {
        mvd->SkipFlag = TRUE;   
        *StrIdx = '\0';                   /* terminated index name before '*' */
    }


    /* Search for % character in calibration index name */
    StrIdx = strchr(IndexL, '%');

    if (StrIdx != NULL)
    {
        mvd->CalibSmileFlag = FALSE;   
        *StrIdx = '\0';                   /* terminated index name before '*' */
    }


    /* Cms indices */
    StrIdx = strstr (IndexL, "yCms");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL) * 12;


        /* Read conventions from t_curve */
        mvd->BaseDate = t_curve->Today;
        mvd->Freq     = t_curve->SwapFreq;
        msVol.Freq    = mvd->Freq;

        if (!strcmp(t_curve->SwapDCC, "360"))
        {
            mvd->DCC = '0';
        }
        else if (!strcmp(t_curve->SwapDCC, "365"))
        {
            mvd->DCC = '5';
        }
        else
        {
            mvd->DCC = '3';
        }

        /* Read full swaption vol matrix */
        if ( SwapSmile_Input_W (&msVol,
                                0,
                                0,
                                mvd->BaseDate,
                                isTMXStylFile,
                                LocalSwapVolFile) == FAILURE )
        {
            goto RETURN;
        }

        NbCols = msVol.NbFwdMat;
        NbRows = msVol.NbExpiry;

        /* Read full swaption smile matrix */
        if (MktSmile_Input_W (&msSmile,
                              OWS,
                              mvd->BaseDate,
                              'S',
                              isTMXStylFile,
                              SwapSmlFile) == FAILURE)
        {
            goto RETURN;
        }

        /* Find required Cms column */
        matIdx = 0; 
        while ((matIdx < NbCols-1) && (IdxMat > msVol.FwdMat[matIdx]))
            matIdx++;

        if (IdxMat != msVol.FwdMat[matIdx])
        {        
            DR_Error("Cms index is not in swaption matrix (MktVol_Input_W)!");
            goto RETURN;
        }

        for (i = 0; i < NbRows; i++)
        {                                                                       
            /* 
             * Generate benchmark swap start and end dates
             * VolDate and SwapSt are identical: i.e. we don't calibrate 
             * options on forward starting swaps (e.g. mid curve options) 
             */
            mvd->VolDate[i]    = msVol.VolDate[i];
            mvd->SwapSt[i]     = msVol.VolDate[i];
            mvd->SwapMat[i]    = Nxtmth (msVol.VolDate[i], IdxMat, 1L);
            mvd->SwapMatMos[i] = IdxMat;
            mvd->Vol[i]        = msVol.Vol[0][i][matIdx];
            mvd->VolUsed[i]    = TRUE;
            mvd->SmlLiqDate[i] = 0;        /* initial value */
        }

        mvd->NbVol = NbRows;

        /* Mark Liquid Expiries */
        NbRows = msSmile.NbExpiry;

        for (i = 0; i < NbRows; i++)
        {
            if (isTMXStylFile == TRUE)
            {
                SmlExp  = msSmile.Expiry[i];
            }
            else
            {
                SmlExp  = Nxtmth (mvd->BaseDate, msSmile.Expiry[i], 1L);
            }

            if (SmlExp <= mvd->SwapSt[mvd->NbVol-1])
            {
                ILiqBmk = GetDLOffset (mvd->NbVol,
                                       mvd->SwapSt,
                                       SmlExp,
                                       CbkEXACT);

                if (ILiqBmk >= 0)
                {
                    mvd->SmlLiqDate[ILiqBmk] = 1;
                    SmlLiqFound = TRUE;
                }
                else
                {
                    DR_Error("Liquid Expiry is not in swaption matrix (MktVol_Input_W)!");
                    goto RETURN;
                }
            }
        }

        /* check that at least one liquid smile benchmark has been found */
        if (!SmlLiqFound)
        {
            DR_Error ("There must be at least one smile liquid date");
            goto RETURN;
        }

        /* interpolate smile parameters */
        if (smileinterp (mvd, &msSmile) == FAILURE) goto RETURN;

        if (Convert_VoV (mvd, t_curve) == FAILURE) goto RETURN;

        if (MktVol_Check_W (mvd)  == FAILURE)
        {
            goto RETURN;
        }

        return (SUCCESS);

    } /* Cms case */


    /* Final maturity indices */
    StrIdx = strstr (IndexL, "yFix");

    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL) * 12;


        /* Read conventions from t_curve */
        mvd->BaseDate = t_curve->Today;
        mvd->Freq     = t_curve->SwapFreq;
        msVol.Freq    = mvd->Freq;

        if (!strcmp(t_curve->SwapDCC, "360"))
        {
            mvd->DCC = '0';
        }
        else if (!strcmp(t_curve->SwapDCC, "365"))
        {
            mvd->DCC = '5';
        }
        else
        {
            mvd->DCC = '3';
        }

        /* Read full swaption vol matrix */
        if ( SwapSmile_Input_W (&msVol,
                                0,
                                0,
                                mvd->BaseDate,
                                isTMXStylFile,
                                LocalSwapVolFile) == FAILURE )
        {
            goto RETURN;
        }

        NbCols = msVol.NbFwdMat;
        NbRows = msVol.NbExpiry;

        /* Read full swaption smile matrix */
        if (MktSmile_Input_W (&msSmile,
                              OWS,
                              mvd->BaseDate,
                              'S',
                              isTMXStylFile,
                              SwapSmlFile) == FAILURE)
        {
            goto RETURN;
        }


        /* Process the final maturity index */
        for (i = 0; i < NbRows; i++)
        {
            if (isTMXStylFile == TRUE)
            {
                Mat = IdxMat - Months360(mvd->BaseDate, msVol.Expiry[i]);
            }
        else
            {
                Mat = IdxMat - msVol.Expiry[i];
            }

            if (Mat < msVol.FwdMat[0])
                break;

            j = 0; 
            while ((j < NbCols-1) && (Mat >= msVol.FwdMat[j]))
                j++;

            /* Use higher end of bracket: we don't interpolate to avoid stub */
            /* This includes the case FwdMat[i] > FwdMat[NbCol-1] so that    */
            /* we use a flat volatility after the last forward maturity.     */
            if (2 * Mat >= msVol.FwdMat[j-1] + msVol.FwdMat[j])
            {                                                                       
                Mat = msVol.FwdMat[j];                                          
                mvd->Vol[i]=msVol.Vol[0][i][j];
            }
            else                            
            {
                Mat = msVol.FwdMat[j-1];
                mvd->Vol[i]=msVol.Vol[0][i][j-1];
            }
            
            /* Generate benchmark swap start and end dates */
            mvd->VolDate[i]    = msVol.VolDate[i];
            mvd->SwapSt[i]     = msVol.VolDate[i];
            mvd->SwapMat[i]    = Nxtmth (msVol.VolDate[i], Mat, 1L);
            mvd->SwapMatMos[i] = Mat;
            mvd->VolUsed[i]    = TRUE;
            mvd->SmlLiqDate[i] = 0;
        }

        mvd->NbVol = i;

        if (i == 0)
        {
            DR_Error ("nyFix calibration falls outside swaption matrix "
                      "(MktVol_Input_W)!");
            goto RETURN;
        }
         
        /* Mark Liquid Expiries */
        NbRows = msSmile.NbExpiry;

        for (i = 0; i < NbRows; i++)
        {
            if (isTMXStylFile == TRUE)
            {
                SmlExp  = msSmile.Expiry[i];
            }
            else
            {
                SmlExp  = Nxtmth (mvd->BaseDate, msSmile.Expiry[i], 1L);
            }

            if (SmlExp <= mvd->SwapSt[mvd->NbVol-1])
            {
                ILiqBmk = GetDLOffset (mvd->NbVol,
                                       mvd->SwapSt,
                                       SmlExp,
                                       CbkEXACT);

                if (ILiqBmk >= 0)
                {
                    mvd->SmlLiqDate[ILiqBmk] = 1;
                    SmlLiqFound = TRUE;
                }
                else
                {
                    DR_Error("Liquid Expiry is not in swaption matrix (MktVol_Input_W)!");
                    goto RETURN;
                }
            }
        }

        /* check that at least one liquid smile benchmark has been found */
        if (!SmlLiqFound)
        {
            DR_Error ("There must be at least one smile liquid date");
            goto RETURN;
        }

        /* interpolate smile parameters */
        if (smileinterp (mvd, &msSmile) == FAILURE) goto RETURN;

        if (Convert_VoV (mvd, t_curve) == FAILURE) goto RETURN;

        if (MktVol_Check_W (mvd)  == FAILURE)
        {
            goto RETURN;
        }

        return (SUCCESS);

    } /* Fix case */


    /* Base vol indices */
    StrIdx = strchr(IndexL, 'm');

    /* Check that we are not reprocessing the "yCms" case */
    if (StrIdx != NULL)
    {
        *StrIdx = '\0';
        IdxMat = atoi(IndexL);


        /* Read conventions from t_curve */
        mvd->BaseDate = t_curve->Today;

        if (t_curve->MMB == 365)
        {
            mvd->DCC = '5';
        }
        else
        {
            mvd->DCC = '0';
        }


        /* Read base smile curves */
        if ( BaseSmile_Input_W (&msVol,
                                0,
                                0,
                                mvd->BaseDate,
                                BaseVolFile) == FAILURE )
        {
            goto RETURN;
        }

        NbCols = msVol.NbFwdMat;
        NbRows = msVol.NbExpiry;

        /* Read base smile */
        if (MktSmile_Input_W (&msSmile,
                              OWS,
                              mvd->BaseDate,
                              'B',
                              isTMXStylFile,
                              BaseSmlFile) == FAILURE)
        {
            goto RETURN;
        }


        if (12 / Conv_Freq (msVol.Freq) != IdxMat)
        {
            DR_Error("Base vol curve frequency different from calibration "
                        "index (MktVol_Input_W)!");
            goto RETURN;
        }

        /* Copy to mvd->Frequency (to be used as frequency nmr dates) */
        mvd->Freq = msVol.Freq;

        for (i = 0; i < NbRows; i++)
        {
            /* Generate benchmark swap start and end dates */
            mvd->VolDate[i]    = msVol.VolDate[i];
            mvd->SwapSt[i]     = msVol.VolDate[i];
            mvd->SwapMat[i]    = Nxtmth (msVol.VolDate[i], IdxMat, 1L);
            mvd->Vol[i]     = msVol.Vol[0][i][0];
            mvd->SwapMatMos[i] = IdxMat;
            mvd->VolUsed[i]    = TRUE;
            mvd->SmlLiqDate[i] = 0;
        }

        mvd->NbVol = NbRows;

        /* Mark Liquid Expiries */
        NbRows = msSmile.NbExpiry;

        for (i = 0; i < NbRows; i++)
        {
            SmlExp  = msSmile.VolDate[i];

            if (SmlExp <= mvd->SwapSt[mvd->NbVol-1])
            {
                ILiqBmk = GetDLOffset (mvd->NbVol,
                                       mvd->SwapSt,
                                       SmlExp,
                                       CbkEXACT);

                if (ILiqBmk >= 0)
                {
                    mvd->SmlLiqDate[ILiqBmk] = 1;
                    SmlLiqFound = TRUE;
                }
                else
                {
                    DR_Error("Liquid Expiry is not in swaption matrix (MktVol_Input_W)!");
                    goto RETURN;
                }
            }
        }

        /* check that at least one liquid smile benchmark has been found */
        if (!SmlLiqFound)
        {
            DR_Error ("There must be at least one smile liquid date");
            goto RETURN;
        }

        /* interpolate smile parameters */
        if (smileinterp (mvd, &msSmile) == FAILURE) goto RETURN;

        if (Convert_VoV (mvd, t_curve) == FAILURE) goto RETURN;

        if (MktVol_Check_W (mvd)  == FAILURE)
        {
            goto RETURN;
        }

        return (SUCCESS);

    } /* Basevol case */

    if (StrIdx == NULL)
    {
        DR_Error ("Incorrect calibration index (MktVol_Input_W)!");
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    return (status);

}  /* MktVol_Input_W */



/*****  Fix3_Param_Input_Tmx_Old  ************************************************/
/*
*   Read model parameters and check validity of input.
*/
int     Fix3_Param_Input_Tmx_Old (   
             MKTVOL_DATA    *mvd,               /* (O) Volatility data        */
             FIX3_TREE_DATA *tree_data,         /* (O) Tree data structure    */
             int            NbFactor,           /* (I) Number of factors      */
             char           OWS[6][MAXBUFF],    /* (I) Overwrite strings      */
             char const*    FileName,           /* (I) File name for OU data  */
             char const*    FileNameXSMILE)     /* (I) File name for XSMILE   */
{

    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    *getserror;         /* Reading error for fgets          */
    char    string[MAXBUFF];
    double  Fix3Stuff[MAXBUFF]; /* Unused Fix3 stuff (smile, fwd shift) */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream       = NULL;
    FILE    *streamXSMILE = NULL;
    double  norm;

    /* TMX3: only allow one factor for now */
    if ((NbFactor != 1) /*&& 
        (NbFactor != 2) &&
        (NbFactor != 3) */)
    {
        DR_Error("Nb of factors must be 1"/*, 2 or 3*/"!\n");
        goto RETURN;
    }

    stream = fopen (FileName, "r");

    /*
     *  If there is no parameter file use the overwrite strings.
     */

    if (stream == NULL)
    {
        readerror = sscanf (OWS[0],
                            "%d \n", 
                            &(tree_data->Ppy));

        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not find file %s: Ppy overwrite required! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (NbFactor == 1)
        {
            readerror = sscanf (OWS[2], 
                                "%lf \n", 
                                &(mvd->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OWS[3],
                                "%lf \n", 
                                &(mvd->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* Fill in the unused parameters with N/A values */

            mvd->Alpha[1] = -999.;
            mvd->Alpha[2] = -999.;
            mvd->Beta[1]  = -999.;
            mvd->Beta[2]  = -999.;
            mvd->Rho[0]   = -999.;
            mvd->Rho[1]   = -999.;
            mvd->Rho[2]   = -999.;
        }
        else if (NbFactor == 2)
        {
            readerror = sscanf (OWS[2], 
                                "%lf \t%lf \n", 
                                &(mvd->Alpha[0]),
                                &(mvd->Alpha[1]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OWS[3],
                                "%lf \t%lf \n", 
                                &(mvd->Beta[0]),
                                &(mvd->Beta[1]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OWS[4],
                                "%lf \n", 
                                &(mvd->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            mvd->Alpha[2] = -999.;
            mvd->Beta[2]  = -999.;
            mvd->Rho[1]   = -999.;
            mvd->Rho[2]   = -999.;
        }
        else if (NbFactor == 3)
        {
            readerror = sscanf (OWS[2], 
                                "%lf \t%lf \t%lf \n", 
                                &(mvd->Alpha[0]),
                                &(mvd->Alpha[1]),
                                &(mvd->Alpha[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OWS[3],
                                "%lf \t%lf \t%lf \n", 
                                &(mvd->Beta[0]),
                                &(mvd->Beta[1]),
                                &(mvd->Beta[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OWS[4],
                                "%lf \t%lf \t%lf \n", 
                                &(mvd->Rho[0]),
                                &(mvd->Rho[1]),
                                &(mvd->Rho[2]));

            if (readerror != 3)
            {        
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }  /* if NbFactor == */

        readerror = sscanf (OWS[5],
                            "%lf \n", 
                            &(mvd->Bbq));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find file %s: Bbq overwrite required! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        /* At input level q=0 means log-normal whereas */
        /* internally q=0 means normal: we switch here */
        mvd->Bbq  = 1. - mvd->Bbq;
    
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
                                &(mvd->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* 
             *  If overwrite exists, use it.
             */

            if (strstr (OWS[3], "nil") == NULL)  /* Need to use strstr() here, strcmp won't do */
            {
                readerror = sscanf (OWS[3],
                                    "%lf \n", 
                                    &(mvd->Beta[0]));

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
                                &(mvd->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[2], "nil") == NULL)
            {
                readerror = sscanf (OWS[2], 
                                    "%lf \n", 
                                    &(mvd->Alpha[0]));

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

            if (strstr (OWS[0], "nil") == NULL)
            {
                readerror = sscanf (OWS[0],
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


            mvd->Alpha[1] = -999.;
            mvd->Alpha[2] = -999.;
            mvd->Beta[1]  = -999.;
            mvd->Beta[2]  = -999.;
            mvd->Rho[0]   = -999.;
            mvd->Rho[1]   = -999.;
            mvd->Rho[2]   = -999.;

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Beta[0]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Beta[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[3], "nil") == NULL)
            {
                readerror = sscanf (OWS[3],
                                    "%lf \t%lf \n", 
                                    &(mvd->Beta[0]),
                                    &(mvd->Beta[1]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Alpha[0]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Alpha[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[2], "nil") == NULL)
            {
                readerror = sscanf (OWS[2], 
                                    "%lf \t%lf \n", 
                                    &(mvd->Alpha[0]),
                                    &(mvd->Alpha[1]));

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
                                &(mvd->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[4], "nil") == NULL)
            {
                readerror = sscanf (OWS[4],
                                    "%lf \n", 
                                    &(mvd->Rho[0]));

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

            if (strstr (OWS[0], "nil") == NULL)
            {
                readerror = sscanf (OWS[0],
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

            mvd->Alpha[2] = -999.;
            mvd->Beta[2]  = -999.;
            mvd->Rho[1]   = -999.;
            mvd->Rho[2]   = -999.;

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Beta[0]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Beta[1]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Beta[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[3], "nil") == NULL)
            {
                readerror = sscanf (OWS[3],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mvd->Beta[0]),
                                    &(mvd->Beta[1]),
                                    &(mvd->Beta[2]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Alpha[0]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Alpha[1]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Alpha[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[2], "nil") == NULL)
            {
                readerror = sscanf (OWS[2], 
                                    "%lf \t%lf \t%lf \n", 
                                    &(mvd->Alpha[0]),
                                    &(mvd->Alpha[1]),
                                    &(mvd->Alpha[2]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Rho[0]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Rho[1]));

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

            readerror = fscanf (stream, "%lf \n", &(mvd->Rho[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OWS[4], "nil") == NULL)
            {
                readerror = sscanf (OWS[4],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mvd->Rho[0]),
                                    &(mvd->Rho[1]),
                                    &(mvd->Rho[2]));

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

            if (strstr (OWS[0], "nil") == NULL)
            {
                readerror = sscanf (OWS[0],
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
                            Fix3Stuff);
        
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
                            Fix3Stuff);

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
                            Fix3Stuff);

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
                            &(mvd->CetNbIter));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find CetNbIter in file %s! (Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
        if (strstr (OWS[5], "nil") != NULL)
        {
            if (FindAndSkipComLine (stream, "Bbq", "Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                "%lf \n", 
                                &(mvd->Bbq));

            if (readerror != 1)
            {      
                sprintf (ErrorMsg, "Could not find Bbq in file %s! (Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            mvd->Bbq = 1. - mvd->Bbq;

        }
        else
        {
            readerror = sscanf (OWS[5],
                                "%lf \n", 
                                &(mvd->Bbq));

            if (readerror != 1)
            {      
                sprintf (ErrorMsg, "Could not read overwrite string for Bbq! (Param_Input)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }       
            
            mvd->Bbq  = 1. - mvd->Bbq;

        }  /* if then else */

    }  /* if then else */

    /* Set total vol constants */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) norm += mvd->Alpha[i] * mvd->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
        goto RETURN;
    }       
    if (IS_EQUAL(mvd->Bbq,1))
    {
        ;
    }
    else
    {
        DR_Error("Bbq parameter must be 0 (lognormal Bbq)!");
        goto RETURN;
    }      
    
    /*-------------------------------------------*/
    /* New in TMX: tail definitions and numerics */
    /*-------------------------------------------*/
    streamXSMILE = fopen (FileNameXSMILE, "r");

    if (streamXSMILE == NULL)
    {
        DR_Error ("Could not open parameter file %s\n", FileNameXSMILE);
        goto RETURN;
    }

    /* First line is always a comment                         */
    if (FindAndSkipComLine (streamXSMILE, "Header/MultiQ Nb Std", "Param_Input",
                            FileNameXSMILE) == FAILURE)
    {        
        goto RETURN;    
    }

    /* Second line is either a comment or the Normal Cutoff parameter: */    
    /* find out which                                                  */
    if (fgets (string, MAXBUFF, streamXSMILE) == NULL)
    {
        DR_Error ("Could not find 2nd line - comment or value - in file %s! (Param_Input)", 
                  FileNameXSMILE);
        goto RETURN;    
    }

    if (string[0] == '#')
    {    
        /* Ignore the line as it is a comment                      */      
        if (fgets (string, MAXBUFF, streamXSMILE) == NULL)
        {
            DR_Error ("Could not find MultiQ Nb Std in file %s! (Param_Input)", 
                      FileNameXSMILE);       
            goto RETURN;
        }
    }

    readerror = sscanf (string, "%lf \n", &(mvd->NbSigmaMQ));

    if (readerror != 1)
    {
        DR_Error ("Could not find MultiQ Nb Std in file %s! (Param_Input)", 
                  FileNameXSMILE);
        goto RETURN;    
    }

    if (FindAndSkipComLine (streamXSMILE, "MultiQ Nb Nck", "Param_Input", 
                            FileNameXSMILE) == FAILURE)        
    {        
        goto RETURN;
    }

    readerror = fscanf (streamXSMILE, "%lf \n", &(mvd->NckMQ));
        
    if (readerror != 1)
    {      
        DR_Error ("Could not find MultiQ Nck in file %s! (Param_Input)", 
                  FileNameXSMILE);
        goto RETURN;
    }

    /* Check validity of input */
    if (Fix3_Param_Check_Tmx (NbFactor,
                              mvd,
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
    if (streamXSMILE != NULL)
    {
        fclose (streamXSMILE);
    }

    return (status);

}  /* Fix3_Param_Input_Tmx_Old */



/*****  Fix3_Param_Input_Tmx ************************************************/
/*
*   Read model parameters and check validity of input.
*/
int     Fix3_Param_Input_Tmx (   
             MKTVOL_DATA    *mktvol_data,       /* (O) Volatility data        */
             FIX3_TREE_DATA *tree_data,         /* (O) Tree data structure    */
             char const*   FileNameTreeNum,     /* (I) File name tree num par */
             char const*   FileNameModelPar)    /* (I) File name model par    */
             
{
    int     status = FAILURE;   /* Error status = FAILURE initially */

    /* Read model parameters file */
    if (Fix3_Model_Input_Tmx (mktvol_data, 
                              tree_data,
                              FileNameModelPar) == FAILURE)
    {
        goto RETURN;
    }

    if (tree_data->NbFactor != 1)
    {
        DR_Error ("TMX model: number of factors must be equal to 1.\n");
        goto RETURN;
    }
            
    /* Read numerical parameters file */
    if (Fix3_Num_Input_New(mktvol_data, 
                           tree_data,
                           FileNameTreeNum) == FAILURE)
    {
        goto RETURN;
    }

    /* Check validity of input */
    if (Fix3_Param_Check_Tmx (tree_data->NbFactor,
                              mktvol_data,
                              tree_data) == FAILURE)              
    {        
        goto RETURN;
    }
        
    status = SUCCESS;
        
  RETURN:
        
    return status;

}  /* Fix3_Param_Input_Tmx */


/*****  Fix3_Param_Check_Tmx  ************************************************/
/*
*   Read term structure input and check validity of input.
*/
int     Fix3_Param_Check_Tmx (
                     int             NbFactor,    /* (I) Number of factors   */
                     MKTVOL_DATA     *mvd,        /* (I) Market vol data     */
                     FIX3_TREE_DATA  *tree_data)  /* (I) Tree data structure */
{
    int  status = FAILURE;     /* Error status = FAILURE initially */
    int  i;
    double norm;

    /* Nb of factors */
    if (mvd->NbFactor != NbFactor || tree_data->NbFactor != NbFactor)
    {
        DR_Error("Inconsistent number of factors! ");
        goto RETURN;
    }

    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    /* Set total vol constants */
    norm = 0.;
    for (i = 0; i < NbFactor; i++) norm += mvd->Alpha[i] * mvd->Alpha[i]; 
    norm = sqrt(norm);
    if (fabs(norm) < ERROR)
    {      
        DR_Error("Total alpha is too small !");
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

    if (NbFactor == 1)
    {
        if (mvd->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[0] < -10.) || (mvd->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }
    }
    else if (NbFactor == 2)
    {
        if (mvd->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if (mvd->Alpha[1] < 0.0001)
        {
            DR_Error("Weight #2 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[0] < -10.) || (mvd->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[1] < -10.) || (mvd->Beta[1] > 10.))
        {
            DR_Error("Beta #2 out of range!");
            goto RETURN;
        }

        if ((mvd->Rho[0] < -.95) || (mvd->Rho[0] > .95))
        {
            DR_Error("Correlation out of range!");
            goto RETURN;
        }
    }
    else if (NbFactor == 3)
    {
        if (mvd->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 out of range!");
            goto RETURN;
        }

        if (mvd->Alpha[1] < 0.0001)
        {
            DR_Error("Weight #2 out of range!");
            goto RETURN;
        }

        if (mvd->Alpha[2] < 0.0001)
        {
            DR_Error("Weight #3 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[0] < -10.) || (mvd->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[1] < -10.) || (mvd->Beta[1] > 10.))
        {
            DR_Error("Beta #2 out of range!");
            goto RETURN;
        }

        if ((mvd->Beta[2] < -10.) || (mvd->Beta[2] > 10.))
        {
            DR_Error("Beta #3 out of range!");
            goto RETURN;
        }

        if ((mvd->Rho[0] < -.95) || (mvd->Rho[0] > .95))
        {
            DR_Error("Correlation #1 out of range!");
            goto RETURN;
        }

        if ((mvd->Rho[1] < -.95) || (mvd->Rho[1] > .95))
        {
            DR_Error("Correlation #2 out of range!");
            goto RETURN;
        }

        if ((mvd->Rho[2] < -.95) || (mvd->Rho[2] > .95))
        {
            DR_Error("Correlation #3 out of range!");
            goto RETURN;
        }
    }  /* if then else */

    /* Check Cet Nb Iter */
    if (mvd->CetNbIter < 0)
    {
        DR_Error("The Nb of Cet Iters must be greater than 0!");
        goto RETURN;
    }

    /* Bbq parameters */
    if ((mvd->Bbq  < -10.) || (mvd->Bbq  > 10.))
    {
        DR_Error("Bbq out of range!");
        goto RETURN;
    }


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Fix3_Param_Check_Tmx */


/*****  Convert_VoV  *******************************************************/
/*
*       Returns % VoV for Vol Bmks and resets BB's to 0..
*/
int   Convert_VoV (MKTVOL_DATA   *mvd,        /* (I/O) Volatility data     */
                   T_CURVE const *t_curve)    /* (I) Term structure data   */
{
    int     i, NbVol;
    double  ParYield0, Annuity0;

    int     status = FAILURE;   /* Error status = FAILURE initially        */

    /* Initialisations                                                     */
    NbVol = mvd->NbVol;

    /* Conversion                                                          */
    for (i = 0; i < NbVol; i++)
    {
        if (ParYieldFromDates (&ParYield0,
                               &Annuity0,
                               mvd->SwapSt[i],
                               mvd->SwapMat[i],
                               mvd->DCC,
                               mvd->Freq,
                               'F',
                               t_curve) == FAILURE)
        {
            goto RETURN;
        }

        mvd->Smile[2][i] = mvd->Smile[2][i]    * 
            pow(mvd->Vol[i], mvd->Smile[3][i]) *
            pow(ParYield0, mvd->Smile[4][i]);
        mvd->Smile[3][i] = 0.;
        mvd->Smile[4][i] = 0.;
    }
    status = SUCCESS;

RETURN:

    return (status);

}  /* Convert_VoV */




/*****Fix3_Model_Input_Tmx  *********************************************************/
/*
* Read model parameters
*/
int Fix3_Model_Input_Tmx (MKTVOL_DATA    *mktvol_data, /* (O) Volatility data               */
                          FIX3_TREE_DATA *tree_data,   /* (0) Tree data                     */
                          char const     *FileName)    /* (I) File name including extension */   
{
    int     NbTDInp, NbFactor,  NbCorr;
    long    idxLo, idxHi, s;
    int     NbRows=0;           /* Swaption volatility matrix size */
    int     ILiqBmk;
    int     i, j, k;
    long    SmlExp;
    int     SmlLiqFound = FALSE;

    MKTSMILE_DATA ms;           /*temporary smile structure         */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;
    int     readerror;
    int     status = FAILURE;

    stream = fopen (FileName, "r");
    if (stream == NULL)
    {
        DR_Error("Could not open model parameters file: %s.\n", FileName);
        goto RETURN;
    }
  
    idxLo = 1; idxHi = 8;


    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line",
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Section:VNFM parameters", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
   

     /* Number of factors */
    if (FindAndSkipComLine (stream, 
                            "Number of factors", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbFactor);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of factors in file %s! "
                "(Model_Input_W)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
   
    tree_data->NbFactor = mktvol_data->NbFactor = NbFactor;

    /*compute number of correlations */
    NbCorr = NbFactor * (NbFactor - 1) / 2;

    /* Number of time dependent input  dates */
    if (FindAndSkipComLine (stream, 
                            "Number of benchmark dates", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbTDInp);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of time dependent inputs dates in file %s! "
                "(Fix3_Model_Input_Tmx)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /* if classic version of the tree one and only one time dep input date */
    if (NbTDInp != 1)
    {
        sprintf(ErrorMsg, 
                "Only one time dependent input date should be specified for TMX version in file %s! "
                "(Fix3_Model_Input_Tmx)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
  
   
    mktvol_data->NbTDInp = NbTDInp;
   

    if (NbFactor == 1)
    {
        mktvol_data->Alpha[1] = -999.;
        mktvol_data->Alpha[2] = -999.;
        mktvol_data->Beta[1]  = -999.;
        mktvol_data->Beta[2]  = -999.;
        mktvol_data->Rho[0]  = -999.;
        mktvol_data->Rho[1]   = -999.;
        mktvol_data->Rho[2]   = -999.;

    }
    else if (NbFactor == 2)
    {
        mktvol_data->Alpha[2] = -999.;
        mktvol_data->Beta[2]  = -999.;
        mktvol_data->Rho[1]  = -999.;
        mktvol_data->Rho[2]   = -999.;

    }

    /* Term structure of weights, mean-reversion, correlation */
    if (FindAndSkipComLine (stream, 
                            "Term structure of , weights,mean-reversion, correlation", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld ", &(mktvol_data->TDInpDate[0]));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
            "Could not read td input date in file %s! "
            "(Fix3_Model_Input_Tmx)", 
            FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    } 
    mktvol_data->TDInpDate[0] = IRDateFromYMDDate(mktvol_data->TDInpDate[0]);
  
	for (k = 0; k < NbFactor; k++)
	{
        readerror = fscanf (stream, "%lf ", &mktvol_data->Beta[k]);
        if (readerror != 1)
		{        
            sprintf(ErrorMsg, 
                "Could not read mean-reversion input for factor %d in file %s! "
                "(Fix3_Model_Input_Tmx)", 
				k,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
	}

    for (k = 0; k < NbFactor; k++)
	{
        readerror = fscanf (stream, "%lf ", &mktvol_data->Alpha[k]);
        if (readerror != 1)
		{        
            sprintf(ErrorMsg, 
                "Could not read factor weight input for factor %d  in file %s! "
                "(Fix3_Model_Input_Tmx)", 
				k,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
		}
	}

    for (k = 0; k < NbCorr; k++)
	{
        readerror = fscanf (stream, "%lf ", &mktvol_data->Rho[k]);
        if (readerror != 1)
		{        
            sprintf(ErrorMsg, 
                "Could not read correlation input for factor %d in file %s! "
                "(Fix3_Model_Input_Tmx)", 
				k,
                FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
		}
	}
    fscanf (stream, "\n");
    
    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Backbone section */
    if (FindAndSkipSectionLine (
                             1, /* skip mode */
                             stream, 
                            "Section:Backbone", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipComLine(
                      stream, "Backbone", "Fix_Model_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream,
                        "%lf \n", 
                        &(mktvol_data->Bbq));

    if (readerror != 1)
    {      
        sprintf (ErrorMsg, "Could not find Bbq in file %s! (Fix3_Model_Input_Tmx)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    mktvol_data->Bbq = 1. - mktvol_data->Bbq;

     if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    /* SMILE SECTION */
    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Smile", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Smile number of entries */
    if (FindAndSkipComLine (stream, 
                            "Number of entries", 
                            "Fix3_Model_Input_Tmx", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

  

    readerror = fscanf (stream, "%ld \n", &(ms.NbExpiry));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of expiries in %s!", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    
    if (ms.NbExpiry > MAXNBDATE)
    {        
        sprintf (ErrorMsg, "Nb of vols in %s exceeds maximum of %d! ", FileName, MAXNBDATE);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
    ms.NbFwdMat = 1;

    /* set frequency to mktvol_data freq */
    ms.Freq = mktvol_data->Freq;
    switch (ms.Freq)
    {
        case ('A'): ms.FwdMat[0] = 12; break;
        case ('S'): ms.FwdMat[0] = 6;  break;
        case ('Q'): ms.FwdMat[0] = 3;  break;
        case ('M'): ms.FwdMat[0] = 1;  break;
        default:    ms.FwdMat[0] = -999;
    }

    
    if (FindAndSkipComLine (stream, 
                            "Dates SmileParams", 
                            "Fix3_Model_Input_Tmx",
                            FileName) == FAILURE)
    {
        goto RETURN;
    }

    for (i = 0; i < ms.NbExpiry; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf \t%lf\n", 
                            &(ms.VolDate[i]),
                            &(ms.Vol[1][i][0]),&(ms.Vol[2][i][0]),
                            &(ms.Vol[3][i][0]),&(ms.Vol[4][i][0]),
                            &(ms.Vol[5][i][0]),&(ms.Vol[6][i][0]),
                            &(ms.Vol[7][i][0]),&(ms.Vol[8][i][0]));

        if (readerror != 9)
        {        
            sprintf (ErrorMsg, "Could not read vol date & smile params #%d "
                               "in %s! ", i+1, FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        ms.VolDate[i] = IRDateFromYMDDate (ms.VolDate[i]);
    }

    /* Eliminate dates falling before base date */
    j = 0;
    while (ms.VolDate[j] <= mktvol_data->BaseDate)
        j++;

    ms.NbExpiry -= j;

    for (i = 0; i < ms.NbExpiry; i++)
    {
        ms.VolDate[i] = ms.VolDate[i+j];
        for (s = idxLo; s <= idxHi ; s++) ms.Vol[s][i][0] = ms.Vol[s][i+j][0];
    }

    for (i = 0; i < mktvol_data->NbVol; i++)
    {
        mktvol_data->SmlLiqDate[i] = 0;
    }

    NbRows = ms.NbExpiry;
    for (i = 0; i < NbRows; i++)
    {
        SmlExp  = ms.VolDate[i];
       
        
        if (SmlExp <= mktvol_data->SwapSt[mktvol_data->NbVol-1])
        {
            ILiqBmk = GetDLOffset (mktvol_data->NbVol,
                                   mktvol_data->SwapSt,
                                   SmlExp,
                                   CbkEXACT);

            if (ILiqBmk >= 0)
            {
                mktvol_data->SmlLiqDate[ILiqBmk] = 1;
                SmlLiqFound = TRUE;
            }
            else
            {
                DR_Error("Liquid Expiry is not in swaption matrix (MktVol_Input_W)!");
                goto RETURN;
            }
        }
    }

    /* check that at least one liquid smile benchmark has been found */
    if (!SmlLiqFound)
    {
        DR_Error ("There must be at least one smile liquid date");
        goto RETURN;
    }

    /* smile parameter interpolation */
    if (smileinterp (mktvol_data,
                     &ms) == FAILURE)
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

} /* Fix3_Model_Input_Tmx */
