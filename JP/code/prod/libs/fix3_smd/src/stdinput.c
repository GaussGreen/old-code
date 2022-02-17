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
#include "fix123head.h"


/*****  Fix3_Param_Input  ********************************************************/
/*
*   Read model parameters and check validity of input.
*/
int     Fix3_Param_Input (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             FIX3_TREE_DATA     *tree_data,                  /* (O) Tree data structure           */
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
            sprintf (ErrorMsg, "Could not find file %s: Ppy overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: Q overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]));

            if (readerror != 2)
            {        
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: factor weight overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: mean reversion overwrite required! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not find file %s: correlation overwrite required! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        }  /* if NbFactor == */

        readerror = sscanf (OverWriteString[5],
                            "%lf \n", 
                            &(mktvol_data->Bbq));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find file %s: Bbq overwrite required! (Fix3_Param_Input)", FileName);
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

            if (FindAndSkipComLine (stream, "one factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor ppy", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%d \n", 
                                &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Fix3_Param_Input)");
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
                    sprintf (ErrorMsg, "Could not find two factor parameters in file %s! (Fix3_Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find three factor parameters in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not find one factor parameters in file %s! (Fix3_Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor correlation", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for correlation! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor ppy", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find three factor parameters in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not find one factor parameters in file %s! (Fix3_Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            for (i = 0; i < 12; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    sprintf (ErrorMsg, "Could not find two factor parameters in file %s! (Fix3_Param_Input)", FileName);
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find mean reversion in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for mean reversion! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find factor weight in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for factor weight! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[1]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Fix3_Param_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[2]));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find correlation in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for correlation! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor ppy", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Ppy in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for Ppy! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

        }  /* if NbFactor == */

        if (FindAndSkipComLine (stream, "QLeft", "Fix3_Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->QLeft));
        
        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find QLeft in file %s! (Fix3_Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "QRight", "Fix3_Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->QRight));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find QRight in file %s! (Fix3_Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "FwdShift", "Fix3_Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%lf \n", 
                            &(mktvol_data->FwdShift));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find FwdShift in file %s! (Fix3_Param_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

        if (FindAndSkipComLine (stream, "CetNbIter", "Fix3_Param_Input", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream,
                            "%d \n", 
                            &(mktvol_data->CetNbIter));

        if (readerror != 1)
        {      
            sprintf (ErrorMsg, "Could not find CetNbIter in file %s! (Fix3_Param_Input)", FileName);
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
                    sprintf (ErrorMsg, "Could not read overwrite string for Qweight! (Fix3_Param_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }       
                
                mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                mktvol_data->QRight = 1. - mktvol_data->QRight;
            }
        }  /* if then else */
        
        if (strstr (OverWriteString[5], "nil") != NULL)
        {
            if (FindAndSkipComLine (stream, "Bbq", "Fix3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                "%lf \n", 
                                &(mktvol_data->Bbq));

            if (readerror != 1)
            {      
                sprintf (ErrorMsg, "Could not find Bbq in file %s! (Fix3_Param_Input)", FileName);
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
                sprintf (ErrorMsg, "Could not read overwrite string for Bbq! (Fix3_Param_Input)");
                DR_Error(ErrorMsg);
                goto RETURN;
            }       
            
            mktvol_data->Bbq  = 1. - mktvol_data->Bbq;

        }  /* if then else */

    }  /* if then else */

    /* Check validity of input */
    if (Fix3_Param_Check (   NbFactor,
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

}  /* Fix3_Param_Input */



/*****  Fix3_Param_Check  ********************************************************/
/*
*   Read term structure input and check validity of input.
*/
int     Fix3_Param_Check (int         NbFactor,     /* (I) Number of factors                     */
                     MKTVOL_DATA *mktvol_data, /* (I) Structure of swaption volatility data */
                     FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure                   */
{
    int
        status = FAILURE;                      /* Error status = FAILURE initially */
       
    char
        ErrorMsg[MAXBUFF];

    double  norm;

    /* Nb of factors */
    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }
    norm = 1.0;


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

}  /* Fix3_Param_Check */

/*****  Smd_Corr_Input  ********************************************************/
/*
*   Read new SMD model parameters.
*/
int     Smd_Corr_Input (   
             MKTVOL_DATA   *mktvol_data,             /* (O) Volatility data               */
             char          OverWriteString[MAXBUFF], /* (I) Overwrite strings             */
             char const*   FileName)                 /* (I) File name including extension */
{

    int     readerror;          /* Reading error status             */
    double  dummy;
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    ErrorMsg[MAXBUFF];
    FILE    *stream = NULL;


    stream = fopen (FileName, "r");

    /*
     *  If there is no parameter file use the overwrite strings.
     */

    if (stream == NULL)
    {
        readerror = sscanf (OverWriteString,
                            "%lf %lf %lf %lf\n", 
                             &(mktvol_data->Afac),
                             &(mktvol_data->Bfac),
                             &(mktvol_data->Cfac),
                             &(mktvol_data->Dfac));
        
        if (readerror != 4)
        {
            sprintf (ErrorMsg, "Could not find file %s: Q overwrite required! (Smd_Corr_Input)", FileName);
            DR_Error(ErrorMsg);
            goto RETURN;
        }

           
    }
    else
    {
        if (FindAndSkipComLine (stream, "A factor", "Smd_Corr_Input", FileName) == FAILURE)
        {
            goto RETURN;
        }

        /* 
         *  Parameter file exists: use the overwrite strings if they are not 
         *  "nil", otherwise use the values in the parameter file.
         */
         
         /*
          *  Read A factor in parameter file.
          */

            if (FindAndSkipComLine (stream, "A factor", "Smd_Corr_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Afac));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find Afactor in file %s! (Smd_Corr_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            /* 
             *  If overwrite exists, use it.
             */

            if (strstr (OverWriteString, "nil") == NULL)                     /* Need to use strstr() here, strcmp won't do */
            {
                readerror = sscanf (OverWriteString,
                                    "%lf \n", 
                                    &(mktvol_data->Afac));

                if (readerror != 1)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for A factor! (Smd_Corr_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }


            /* B Factor */
            if (FindAndSkipComLine (stream, "B Factor", "Smd_Corr_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Bfac));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find B Factor in file %s! (Smd_Corr_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString, "nil") == NULL)
            {
                readerror = sscanf (OverWriteString, 
                                    "%lf %lf \n",
                                    &dummy,
                                    &(mktvol_data->Bfac));

                if (readerror != 2)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for B Factor! (Smd_Corr_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }




            /* C Factor */
            if (FindAndSkipComLine (stream, "C Factor", "Smd_Corr_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Cfac));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find C Factor in file %s! (Smd_Corr_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString, "nil") == NULL)
            {
                readerror = sscanf (OverWriteString, 
                                    "%lf %lf %lf \n",
                                    &dummy, &dummy,
                                    &(mktvol_data->Cfac));

                if (readerror != 3)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for C Factor! (Smd_Corr_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }






            /* D Factor */
            if (FindAndSkipComLine (stream, "D Factor", "Smd_Corr_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Dfac));

            if (readerror != 1)
            {        
                sprintf (ErrorMsg, "Could not find D Factor in file %s! (Smd_Corr_Input)", FileName);
                DR_Error(ErrorMsg);
                goto RETURN;
            }

            if (strstr (OverWriteString, "nil") == NULL)
            {
                readerror = sscanf (OverWriteString, 
                                    "%lf %lf %lf %lf \n",
                                    &dummy, &dummy, &dummy,
                                    &(mktvol_data->Dfac));

                if (readerror != 4)
                {        
                    sprintf (ErrorMsg, "Could not read overwrite string for D Factor! (Smd_Corr_Input)");
                    DR_Error(ErrorMsg);
                    goto RETURN;
                }
            }

    }  /* if then else */

        
    status = SUCCESS;
        
    RETURN:
        
    if (stream != NULL)
    {
        fclose (stream);
    }

    return (status);

}  /* Smd_Corr_Input */

