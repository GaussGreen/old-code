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
#include <ctype.h>
#include <malloc.h>
#include <math.h>
#include "cupslib.h"


/*****  FindAndSkipComLine  *************************************************/
/**
*       Try to read a comment line and skip it is present. If no comment line
*       is present, make sure it was the last line in code
*/
int     Hyb3_FindAndSkipComLineOptional (FILE         *stream,
                            char const   *Label,
                            char const   *Routine,
                            char const   *FileName)
{
    char    ErrorMsg[MAXBUFF];
    char    string[MAXBUFF];
    int     status = FAILURE;   /* Error status = FAILURE initially */

    if (fgets (string, MAXBUFF, stream) == NULL) 
    {
        goto RETURN; /* return FAILURE, but print no error message */
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
    }  /* if */


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* FindAndSkipComLine */


/*****  Hyb3_Param_Input  ********************************************************/
/**
*   Read model parameters and check validity of input.
*/
int     Hyb3_Param_Input (   
             MKTVOL_DATA    *mktvol_data,                /**< (O) Volatility data               */
             char            calibIndex[MAXINDEX],       /**< (O) IR Calibration index          */
             HYB3_TREE_DATA *tree_data,                  /**< (O) Tree data structure           */
             int             NbFactor,                   /**< (I) Number of factors             */
             char            OverWriteString[5][MAXBUFF],/**< (I) Overwrite strings             */
             char const     *FileName)                   /**< (I) File name including extension */
{

    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    *getserror;         /* Reading error for fgets          */
    char    string[MAXBUFF];
    FILE    *stream = NULL;
    int     backboneFlagGiven = TRUE; /* needed for IR-Index overwrite */

    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    mktvol_data->NbFactor = NbFactor;

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
            DR_Error("Could not find file %s: Ppy overwrite required! (Hyb3_Param_Input)", FileName);
            
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
                DR_Error("Could not find file %s: Q overwrite required! (Hyb3_Param_Input)", FileName);
                
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
                DR_Error("Could not find file %s: factor weight overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Hyb3_Param_Input)", FileName);
                
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
                DR_Error("Could not find file %s: factor weight overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]));

            if (readerror != 2)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find file %s: correlation overwrite required! (Hyb3_Param_Input)", FileName);
                
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
                DR_Error("Could not find file %s: factor weight overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]),
                                &(mktvol_data->Beta[2]));

            if (readerror != 3)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Rho[0]),
                                &(mktvol_data->Rho[1]),
                                &(mktvol_data->Rho[2]));

            if (readerror != 3)
            {        
                DR_Error("Could not find file %s: correlation overwrite required! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }
        }  /* if then else */
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

            if (FindAndSkipComLine (stream, "one factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
                
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
                    DR_Error("Could not read overwrite string for mean reversion! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \n", 
                                    &(mktvol_data->Alpha[0]));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for factor weight! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "one factor ppy", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%d \n", 
                                &(tree_data->Ppy));

            if (readerror != 1)
            {        
                DR_Error("Could not find Ppy in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for Ppy! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            /* No overwrite for Qweight: look for it in the parameter file */
            /* Skip two and three factor parameters in the file */
            for (i = 0; i < 12; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("Could not find two factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
                    goto RETURN;
                }
            }

            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("Could not find three factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
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
                    DR_Error("Could not find one factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
                
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
                    DR_Error("Could not read overwrite string for mean reversion! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                
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
                    DR_Error("Could not read overwrite string for factor weight! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor correlation", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find correlation in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \n", 
                                    &(mktvol_data->Rho[0]));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for correlation! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "two factor ppy", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                DR_Error("Could not find Ppy in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for Ppy! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            /* skip the 3 factor section in the file */
            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("Could not find three factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
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
                    DR_Error("Could not find one factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
                    goto RETURN;
                }
            }

            for (i = 0; i < 12; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("Could not find two factor parameters in file %s! (Hyb3_Param_Input)", FileName);
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[1]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor mean reversion", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Beta[2]));

            if (readerror != 1)
            {        
                DR_Error("Could not find mean reversion in file %s! (Hyb3_Param_Input)", FileName);
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
                    DR_Error("Could not read overwrite string for mean reversion! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[1]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor weight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Alpha[2]));

            if (readerror != 1)
            {        
                DR_Error("Could not find factor weight in file %s! (Hyb3_Param_Input)", FileName);
                
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
                    DR_Error("Could not read overwrite string for factor weight! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find correlation in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[1]));

            if (readerror != 1)
            {        
                DR_Error("Could not find correlation in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "three factor correlation", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%lf \n", &(mktvol_data->Rho[2]));

            if (readerror != 1)
            {        
                DR_Error("Could not find correlation in file %s! (Hyb3_Param_Input)", FileName);
                
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
                    DR_Error("Could not read overwrite string for correlation! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "three factor ppy", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));

            if (readerror != 1)
            {        
                DR_Error("Could not find Ppy in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }

            if (strstr (OverWriteString[0], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[0],
                                    "%d \n", 
                                    &(tree_data->Ppy));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for Ppy! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }
            }

        }  /* if then else */


        /* read environment driven Calibration index (for smile families) */
        /* this part is the same for any IR dimension                     */
        if (strstr (OverWriteString[1], "nil") != NULL)
        {
            
            if (FindAndSkipComLine (stream, "QLeft", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }
            
            readerror = fscanf (stream,
                "%lf \n", 
                &(mktvol_data->QLeft));
            
            if (readerror != 1)
            {      
                DR_Error("Could not find QLeft in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }
            
            if (FindAndSkipComLine (stream, "QRight", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }
            
            readerror = fscanf (stream,
                "%lf \n", 
                &(mktvol_data->QRight));
            
            if (readerror != 1)
            {      
                DR_Error("Could not find QRight in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }
            
            if (FindAndSkipComLine (stream, "FwdShift", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }
            
            readerror = fscanf (stream,
                "%lf \n", 
                &(mktvol_data->FwdShift));
            
            if (readerror != 1)
            {      
                DR_Error("Could not find FwdShift in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }
            
            if (FindAndSkipComLine (stream, "CetNbIter", "Hyb3_Param_Input", FileName) == FAILURE)
            {        
                goto RETURN;
            }
            
            readerror = fscanf (stream,
                "%d \n", 
                &(mktvol_data->CetNbIter));
            
            if (readerror != 1)
            {      
                DR_Error("Could not find CetNbIter in file %s! (Hyb3_Param_Input)", FileName);
                
                goto RETURN;
            }
            
            mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
            mktvol_data->QRight = 1. - mktvol_data->QRight;
            
        }
        else
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
                    DR_Error("Could not read overwrite string for Qs! (Hyb3_Param_Input)");
                    
                    goto RETURN;
                }       
                
                mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                mktvol_data->QRight = 1. - mktvol_data->QRight;
            }

            /* move file pointer to the end of the QLeft,Qright overwrite position */
            for (i = 0; i < 8; i++)
            {
                getserror = fgets (string, MAXBUFF, stream); /* no error if non-existent */
            }

        }  /* if then else */


        /* the next inputs are optional inputs                  */
        /* skip the next line (backbone flag: not used in Hyb3) */

        for (i = 0; i < 2; i++)
        {
            getserror = fgets (string, MAXBUFF, stream);
            if (getserror  == NULL)
            {      
                backboneFlagGiven = FALSE;
            }
        }

        /* overwrite parameter for the calibration index (if not set by trade) */
        if ( backboneFlagGiven && ( strstr (calibIndex, "nil") != NULL )  ) 
        {
            /* check if the overwrite parameter is given (otherwise ignore) */
            getserror = fgets (string, MAXBUFF, stream); /* read comment line if necessary */

            /* skip the comment flag if existent) */
            if ( (getserror != NULL ) && (string[0] == '#') )
            {        
                readerror = fscanf (stream, "%s\n", calibIndex );          
                if (readerror != 1)
                {        
                    DR_Error ("Could not read volatility index in "
                                "file %s! (Hyb3_Param_Input)", FileName);
                    goto RETURN;
                }

            }
        } /* backboneFlagGiven and not set by trade*/


    }  /* if then else */

    
    /* Check validity of input */
    if (Hyb3_Param_Check (   NbFactor,
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

}  /* Hyb3_Param_Input */



/*****  Hyb3_Param_Check  ********************************************************/
/**
*   Read term structure input and check validity of input.
*/
int  Hyb3_Param_Check ( 
                   int          NbFactor,           /**< (I) Number of factors                     */
                   MKTVOL_DATA *mktvol_data,        /**< (I) Structure of swaption volatility data */
                   HYB3_TREE_DATA   *tree_data)     /**< (I) Tree data structure                   */
{
    int
        status = FAILURE;                      /* Error status = FAILURE initially */

    /* Nb of factors */
    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    /* Nb of st dev check */
    if (tree_data->NbSigmaMax < 3)
    {
        DR_Error("Can't cut the tree at less than three std devs, "
                 "supplied value = %d", tree_data->NbSigmaMax);
        goto RETURN;
    }              

    /* Ppy check */
    if ((tree_data->Ppy < 1) || (tree_data->Ppy > 365))
    {
        DR_Error("Ppy has to be between 1 and 365! Supplied value = %d",
                 tree_data->Ppy);
        goto RETURN;
    }

    if ((mktvol_data->QLeft  < -10.) || (mktvol_data->QLeft  > 10.) ||
        (mktvol_data->QRight < -10.) || (mktvol_data->QRight > 10.))
    {
        DR_Error("Q (QLeft = %lf, QRight = %lf) out of range!",
                 mktvol_data->QLeft, mktvol_data->QRight);
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
            DR_Error("Weight #1 (%lf) out of range!", mktvol_data->Alpha[0]);
            goto RETURN;
        }

        if ((mktvol_data->Beta[0] < 0.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 (%lf) out of range!", mktvol_data->Beta[0]);
            goto RETURN;
        }
    }
    else if (NbFactor == 2)
    {
        if (mktvol_data->Alpha[0] < 0.0001)
        {
            DR_Error("Weight #1 (%lf) out of range!", mktvol_data->Alpha[0]);
            goto RETURN;
        }

        if (mktvol_data->Alpha[1] < 0.0001)
        {
            DR_Error("Weight #2 (%lf) out of range!", mktvol_data->Alpha[1]);
            goto RETURN;
        }

        if ((mktvol_data->Beta[0] < 0.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 (%lf) out of range!", mktvol_data->Beta[0]);
            goto RETURN;
        }

        if ((mktvol_data->Beta[1] < 0.) || (mktvol_data->Beta[1] > 10.))
        {
            DR_Error("Beta #2 (%lf) out of range!", mktvol_data->Beta[1]);
            goto RETURN;
        }

        if ((mktvol_data->Rho[0] < -.95) || (mktvol_data->Rho[0] > .95))
        {
            DR_Error("Correlation (%lf) out of range!", mktvol_data->Rho[0]);
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
        DR_Error("Maximum allowed number of iterations is %d, "
                 "number provided = %d", MAX_ITERATIONS, mktvol_data->CetNbIter);
        goto RETURN;
    }   

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Hyb3_Param_Check */


/*****  Hyb3_Eq_Check_W  *********************************************************/
/*
*	Read term structure input and check validity of input.
*/
int    Hyb3_Eq_Check_W(EQ_DATA *eq_data)  /* (I) Structure of eq data */
{



     int
     i,
     status = FAILURE;     /* Error status = FAILURE initially */
                
                
     if (Dateok(eq_data->ValueDate))
     {
          DR_Error ("Incorrect format for eq value date!");
          goto RETURN;
     }  /* if */
     
     if (eq_data->NbVol < 0) 
     {
         DR_Error("Nb of composite vols must be >= 0\n");
         goto RETURN;
     }

     if (eq_data->NbVol > MAXNBDATE)
     {
          DR_Error ("Nb of vols exceeds max limit\n");
          goto RETURN;
     }

     if (eq_data->NbInpSpotVol < 0)
     {
         DR_Error(" Nb of input spot vols must be >= 0 \n");
         goto RETURN;
     }

     if (eq_data->NbInpSpotVol > MAXNBDATE)
     {
          DR_Error ("Nb of input spot vols exceeds max limit\n");
          goto RETURN;
     }
     if((eq_data->NbInpSpotVol == 0) &&
         (eq_data->NbVol == 0))
     {
         DR_Error("Must have at least one equity vol point (composite or spot)\n");
         goto RETURN;
     }

    
     for (i = 1; i <= eq_data->NbVol; i++)
     {    
          if (Dateok(eq_data->VolDate[i]))
          {
              DR_Error("Invalid format for fx composite vol dates\n");
              goto RETURN;
          }
     }
     
     for (i = 1; i <= eq_data->NbVol; i++)
     {              
         if (eq_data->VolDate[i] <= eq_data->VolDate[i-1])
         {
              DR_Error ("EQ volatility points must be "
                        " entered in ascending order!");
              goto RETURN;
         }  /* if */
     }

     if (eq_data->NbVol > 0)
     {

        /****
        
        VEZI We need to fix up the handling of different value dates

        if (eq_data->VolDate[0] != eq_data->ValueDate)
        {
             DR_Error("CompositeVolDate[0] must be same as eq.value date\n");
            goto RETURN;
        }
        
        ****/
     
        if (fabs(eq_data->Vol[0] - eq_data->Vol[1]) > TINY)
        {
             DR_Error("The internal compositevol[0] must be the same"
                        " as compositevol[1] \n");
            goto RETURN;             
        }
     }

     for (i = 0 ; i < eq_data->NbInpSpotVol ; i++)
     {
         if(Dateok(eq_data->InpSpotVolDate[i]))
         {
             DR_Error("incorrect format for Input spot vol dates\n");
             goto RETURN;
         }
     }

     for (i = 0 ; i < (eq_data->NbInpSpotVol - 1) ; i++)
     {
         if(eq_data->InpSpotVolDate[i] >= eq_data->InpSpotVolDate[i+1])
         {
             DR_Error("Spot Vol dates must be entered in ascending order\n");
             goto RETURN;
         }
     }
    
     if (eq_data->NbInpSpotVol > 0) 
     {
         if(eq_data->InpSpotVolDate[0] <= eq_data->ValueDate)
         {
            DR_Error("First SpotVol date must be > eq Value date\n");
            goto RETURN;
         }
     }

     for (i = 0 ; i < eq_data->NbInpSpotVol ; i++)
     {
         if((eq_data->InpSpotVol[i] < 0.000001) || 
            (eq_data->InpSpotVol[i] > 9.99))
         {
             DR_Error("Equity spot vol entry %d (value = %lf) not entered correctly. "
                      "Must be between 0.000001 and 9.99 (decimal value)\n",
                      i, eq_data->InpSpotVol[i]);
             goto RETURN;
         }
     }

     if ((eq_data->NbVol > 0) && (eq_data->NbInpSpotVol > 0))
     {
         if (eq_data->InpSpotVolDate[0] <= eq_data->VolDate[eq_data->NbVol])
         {
             DR_Error("Last composite vol date must be strictly"
                      " < first spot vol date \n");
             goto RETURN;
         }
     }

     for (i = 1; i <= eq_data->NbVol; i++)
          if ((eq_data->Vol[i] < .000001) || (eq_data->Vol[i] > 9.99))
          {
               DR_Error ("Equity comp vol entry %d (value = %lf) not entered correctly, "
                         "Must be between 0.000001 and 9.99 (decimal value)\n",
                         i, eq_data->Vol[i]);
               goto RETURN;
          }  /* if */

     if (eq_data->EqCutOffLevel < 0.0)
     {
         DR_Error("EQ cut off level must be positive\n");
         goto RETURN;
     }
     
     if ((eq_data->EqCutOffFlag != FALSE) &&
         (eq_data->EqCutOffFlag != TRUE))
     {
         DR_Error("EqCutOffFlag should be internally initialised to "
                  "TRUE of FALSE\n");
         goto RETURN;
     }

     if((eq_data->EqCutOffLast != TRUE) &&
         (eq_data->EqCutOffLast != FALSE))
     {
         DR_Error("eqCutOffLast should be internally initialised to "
                    " TRUE or FALSE\n");
         goto RETURN;
     }
     
        
     if( (eq_data->EqBootStrapMode != 0L) &&
         (eq_data->EqBootStrapMode != 1L)  &&
         (eq_data->EqBootStrapMode != 2L))
     {

         DR_Error("The internal initialisation of EqBootstrapMode is"
                    " invalid\n");
         goto RETURN; 
     }

        /* check FX smile parameters */
        if (eq_data->NbSmilePt < 1)
        {
            DR_Error("eq_Input_W: Nb of EQ smile param pts < 1!");
            goto RETURN;
        }
        
        for (i=0; i<eq_data->NbSmilePt; i++)
        {
            if (Dateok(eq_data->SmileDate[i]))
            {
                DR_Error ("Eq_Input_W: incorrect eq smile date format!");
                goto RETURN;
            }
        }
        
        for (i=1; i<eq_data->NbSmilePt; i++)
        {
            if (eq_data->SmileDate[i] <= eq_data->SmileDate[i-1])
            {
                DR_Error ("Eq_Input_W: "
                          "EQ smile dates are not in ascending order!");
                goto RETURN;
            }
        }
        
        if (eq_data->SmileDate[0] <= eq_data->ValueDate)
        {
            DR_Error ("Eq_Input_W: some EQ smile dates <= today!");
            goto RETURN;
        }

    status = SUCCESS;
          
RETURN:
          
    return (status);

}  /* Hyb3_Eq_Check_W */

/*****  Eq_Input_W_WithSmile  **********************************************
*
*	Read eq input for wrappers and check validity of input.
*  Composite vols are input with a 1 OFFSET, wheras SpotVols are input
*  with 0 OFFSET                                                           
****************************************************************************/

int  Hyb3_Eq_Input_W_WithSmile(EQ_DATA   *eq_data,                 /* (O) Fx data            */
                               long      Today,                    /* (I) Zeroth tree node   */
                               char      OverWriteString[MAXBUFF], /* (I) Owrite str         */
                               char      *VolFileName,             /* (I) FXVols file name   */
                               char      *SmileFileName,           /* (I) FX smile file name */
                               char      NbEQSmileOWS[MAXBUFF],    /* (I)                    */
                               char      EQSmileParamOWS[MAXNBDATE][MAXBUFF]   /* (I) */)
{


     long
          Year[2],
          Month[2],
          Day[2];
     int
          i,
          readerror,                 /* Reading error status             */
          status = FAILURE;          /* Error status = FAILURE initially */
     char
          string[MAXBUFF];
     FILE
         *stream        = NULL,
         *streamEQSmile = NULL; 
     


     stream = fopen (VolFileName, "r");           /* Open the fx data file  */

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Eq_Input_W)", VolFileName);
          
          goto RETURN;
          
     }  /* if */
          
     fgets  (string, 80, stream);    /* Read the comment line */
     readerror = fscanf (stream, "%ld \n", &(eq_data->ValueDate));
     eq_data->ValueDate = IRDateFromYMDDate(eq_data->ValueDate);
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Eq_Input_W)", VolFileName);
          
          goto RETURN;
     }
          
     fgets  (string, 80, stream);
     readerror = fscanf (stream, "%lf \n", &(eq_data->Spot));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Eq_Input_W)", VolFileName);
          
          goto RETURN;
     }
     
     fgets ( string, 80, stream);  /* Ignore this input */
     fgets ( string, 80, stream);
     
     fgets ( string, 80, stream);
     readerror = fscanf (stream, "%d \n", &(eq_data->NbVol));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Eq_Input_W)", VolFileName);
          
          goto RETURN;
     }

     if (eq_data->NbVol > MAXNBDATE)
     {
          DR_Error("Nb of vols exceeds max limit of %d in file %s! "
                             "(Eq_Input_W)", MAXNBDATE, VolFileName);
          
          goto RETURN;
     }
     
     Dsplit(eq_data->ValueDate, /* Split value date into month, day and year */
              &(Month[0]), 
              &(Day[0]), 
              &(Year[0]));

     fgets ( string, 80, stream);
     for (i = 1; i <= eq_data->NbVol; i++)  /* Composite vols -> 1 offset*/
     {
          readerror = fscanf (stream, "%ld \t%lf \n", &(eq_data->VolDate[i]), &(eq_data->Vol[i]));
          eq_data->VolDate[i] = IRDateFromYMDDate(eq_data->VolDate[i]);
          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read file %s! "
                   "(Eq_Input_W)", VolFileName);
               
               goto RETURN;
          }
          
          eq_data->Vol[i] /= 100.0;
     }  /* for i */

     /* Repeat the first values for interpolation of pseudocomposite vols */
     /* in Get_TreeSpotVols                                               */

     if (eq_data->NbVol > 0)
     {
        eq_data->VolDate[0] = Today;
        eq_data->Vol[0]     = eq_data->Vol[1];
     }

    /* Input Spot Vol section */
    fgets  (string, 80, stream);    /* comment line */
    readerror = fscanf (stream, "%d \n", &(eq_data->NbInpSpotVol));
    if (readerror != 1) 
    {          
         DR_Error("Could not read Nb of Input Spot Vol points "
                   "in file %s! (Fx_Input_W)", VolFileName);
         
         goto RETURN;
    }
     
    if (eq_data->NbInpSpotVol > MAXNBDATE)
    {
         DR_Error("Nb of Inp Spot Vol exceeds max limit of %d "
                   "in file %s! (Eq_Input_W)", MAXNBDATE, VolFileName);
         
         goto RETURN;
    }

    if (eq_data->NbInpSpotVol > 0)
    {
       fgets  (string, 80, stream);    /* comment line */
       for (i = 0; i < eq_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
       {
             readerror = fscanf (stream, "%ld \t%lf \n",
                                 &(eq_data->InpSpotVolDate[i]),
                                 &(eq_data->InpSpotVol[i]));
             eq_data->InpSpotVolDate[i] = IRDateFromYMDDate(eq_data->InpSpotVolDate[i]);

           if (readerror != 2)
           {          
                  DR_Error("Could not read input Spot Vol"
                      " dates/values in file %s! (Fx_Input_W)", VolFileName);
               
               goto RETURN;
           }
          
           eq_data->InpSpotVol[i] /= 100.0;
       }  /* for i */
    }


    /* deal with possible overlap between composite and spot vols */
    if ((eq_data->NbInpSpotVol > 0) && (eq_data->NbVol > 0))
    {
        long Idx;
        Idx = GetDLOffset(eq_data->NbVol,
                          &(eq_data->VolDate[1]),
                          eq_data->InpSpotVolDate[0],
                          CbkHIGHER);
        /*******************************************************/
        /* if all composite dates are < spotvoldate[0]         */
        /* then Idx = -999, in which case NbVol does not need  */
        /* to be changed.                                      */
        /*******************************************************/
        if (Idx >= 0L)
        {
            eq_data->NbVol = Idx;
        }        
    }
    
    if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
     {

        eq_data->EqCutOffFlag = FALSE;
        eq_data->EqCutOffLast = FALSE;
        eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED; 

     }
     else if (strcmp(OverWriteString,"last")== 0) /* 2.Cut off at last level */
     {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        eq_data->EqCutOffFlag = TRUE;
        eq_data->EqCutOffLast = TRUE;
        eq_data->EqCutOffLevel = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
        eq_data->EqBootStrapMode = EQ_USE_LAST_LEVEL ;
     } 
     else                                         /* 3.Cut off at user level */
     {
          /* Cut off is allowed and there */
          /* is volatility cutoff data    */
          eq_data->EqCutOffFlag = TRUE;
          eq_data->EqCutOffLast = FALSE;
          eq_data->EqCutOffLevel = atof(OverWriteString);
          if (eq_data->EqCutOffLevel < TINY) /* The FAILURE return */
          {
              DR_Error("Unable to convert EQ cut off value!\n");
              goto RETURN;
          }
          eq_data->EqCutOffLevel /= 100.;
          eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;
          
     }
     

    /* Read Smile data file */
    

    if (strstr(NbEQSmileOWS, "nil") != NULL)
    {
        if (strstr(EQSmileParamOWS[0], "nil") == NULL)
        {
            DR_Error("Eq_Input_W: invalid EQ Smile Params input\n");
            goto RETURN;
        }

        /* read from environment file.  If environment smile parameters do not exist,
           the model will default to smile off, by setting a1=a2=a3=0 at the valuedate */
        streamEQSmile = fopen (SmileFileName, "r");
        if (streamEQSmile == NULL)
        {
            eq_data->NbSmilePt = 1;

            /* add arbitary number of days after the value date to ensure the smile date is
               > today as in some cases valueDate = Today. - chose to use 5 days but could have
               use any number */
            eq_data->SmileDate[0] = Nxtday(eq_data->ValueDate, 5);
            eq_data->a1[0] = 0.0;
            eq_data->a2[0] = 0.0;
            eq_data->a3[0] = 0.0;
        }
        else
        {
            /* read Nb FX Smile Param lines */
            if (FindAndSkipComLine (streamEQSmile, 
                                    "Nb EQ Param Lines",
                                    "Eq_Input_W", 
                                    SmileFileName) == FAILURE) goto RETURN;

            readerror = fscanf (streamEQSmile, "%d \n", &(eq_data->NbSmilePt));
            if (readerror != 1)
            {
                DR_Error("Eq_Input_W: Cannot read Nb FX smile param lines");
                goto RETURN;
            }

            if ((eq_data->NbSmilePt < 0) ||
                (eq_data->NbSmilePt > MAXNBDATE))
            {
                DR_Error("Eq_Input_W: Nb EQ smile param lines out of range!");
                goto RETURN;
            }

            if (FindAndSkipComLine (streamEQSmile, 
                                   "EQ Smile Params",
                                   "EQ_Input_W", 
                                   SmileFileName) == FAILURE) goto RETURN;

            for (i=0; i < eq_data->NbSmilePt; i++)
            {
                readerror = fscanf(streamEQSmile, "%ld\t%lf\t%lf\t%lf \n",
                                   &(eq_data->SmileDate[i]),
                                   &(eq_data->a1[i]),
                                   &(eq_data->a2[i]),
                                   &(eq_data->a3[i]));
                eq_data->SmileDate[i] = IRDateFromYMDDate(eq_data->SmileDate[i]);
                if (readerror != 4)
                {          
                    DR_Error("Eq_Input_W: Cannot read EQ smile params");
                    goto RETURN;
                }
            }  
        }
    }
    else
    {
        /* read the overwrite strings */
        readerror = sscanf (NbEQSmileOWS, 
                            "%d \n", 
                            &(eq_data->NbSmilePt));

        if (readerror != 1)
        {        
            DR_Error("Eq_Input_W: Cannot read Nb EQ Smile Param OWS");
            goto RETURN;
        }

        if ((eq_data->NbSmilePt < 0) ||
            (eq_data->NbSmilePt > MAXNBDATE))
        {
            DR_Error("Eq_Input_W: Nb FX smile param lines OWS out of range!");
            goto RETURN;
        }

        for (i = 0; i < eq_data->NbSmilePt; i++)
        {
            readerror = sscanf (EQSmileParamOWS[i], 
                                "%ld\t%lf\t%lf\t%lf \n",
                                &(eq_data->SmileDate[i]),
                                &(eq_data->a1[i]),
                                &(eq_data->a2[i]),
                                &(eq_data->a3[i]));
            eq_data->SmileDate[i] = IRDateFromYMDDate(eq_data->SmileDate[i]);
            if (readerror != 4)
            {          
                DR_Error("Eq_Input_W: Cannot read FX smile params OWS");
                goto RETURN;
            }
        }
    }

    /* Check validity of input */ 
    if (Hyb3_Eq_Check_W (eq_data) == FAILURE)
    {
         goto RETURN;
    }

    status = SUCCESS;
         
RETURN:
       
    if (stream != NULL) fclose(stream);

    return (status);
}  /* Eq_Input_W_WithSmile */

/*****  Hyb3_Eq_Input_Dyn_WType6	**************************************************/
/**
 *     Read dynamic equity input  for DR Wrappers of type 6 (i.e. equity is the 
 *     third asset and the EQ_DATA structure contains three correlation curves).
 *
 *     Since  Type6 DR  Wrapper  include a correlation.dat  file, the overwrite 
 *     strings are dealt with by the correlation input function.
 *
 */
int Hyb3_Eq_Input_Dyn_WType6
        (EQ_DATA    *eq_data,               /**< (O) Equity data              */
         long        ValueDate,             /**< (I) Value date               */
         char       *OverWriteString,       /**< (I) Owrite string            */
         char        DomForFlag,            /**< (I) D for domestic eq        */
         char       *FileName,              /**< (I) Including extension      */
         char        NbEQSmileOWS[MAXBUFF],
         char        EQSmileParamOWS[MAXNBDATE][MAXBUFF])
{
          char
                     string[MAXBUFF];
          int
                     i,
                     readerror,             /* Reading error status         */
                     status = FAILURE;      /* Error status=FAILURE         */
          double
                     Correlation;
          FILE
                     *stream = NULL;

          /* overwrites */

          if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
          {

              eq_data->EqCutOffFlag = FALSE;
              eq_data->EqCutOffLast = FALSE;
              eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
              eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED;  

          }
          else if (strcmp(OverWriteString,"last")== 0) /* 2.Cut off at last level */
          {
              /* Cut off is allowed and there */
              /* is volatility cutoff data    */
              eq_data->EqCutOffFlag = TRUE;
              eq_data->EqCutOffLast = TRUE;
              eq_data->EqCutOffLevel = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
              eq_data->EqBootStrapMode = EQ_USE_LAST_LEVEL ;
          } 
          else                                         /* 3.Cut off at user level */
          {
              /* Cut off is allowed and there */
              /* is volatility cutoff data    */
              eq_data->EqCutOffFlag = TRUE;
              eq_data->EqCutOffLast = FALSE;
              eq_data->EqCutOffLevel = atof(OverWriteString);
              if (eq_data->EqCutOffLevel < TINY) /* The FAILURE return */
              {
                  DR_Error("Unable to convert EQ cut off value!\n");
                  goto RETURN;
              }
              eq_data->EqCutOffLevel /= 100.;
              eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;
          
          }

          /* this function is only used for environments lacking smile information */

          eq_data->NbSmilePt = 1;

          eq_data->SmileDate[0] = Nxtday(ValueDate, 5);

          eq_data->a1[0] = 0.0;
          eq_data->a2[0] = 0.0;
          eq_data->a3[0] = 0.0;

          /* but there may be overwrites */

          if (strstr(NbEQSmileOWS, "nil") != NULL)
          {
              if (strstr(EQSmileParamOWS[0], "nil") == NULL)
              {
                  DR_Error("Eq_Input_Dyn_WType6: invalid EQ Smile Params overwrites.\n");
                  goto RETURN;
              }
          }
          else
          {
              /* read the overwrite strings */
              readerror = sscanf (NbEQSmileOWS, 
                                  "%d \n", 
                                  &(eq_data->NbSmilePt));

              if (readerror != 1)
              {        
                  DR_Error("Eq_Input_Dyn_WType6: Cannot read Nb EQ Smile Param OWS");
                  goto RETURN;
              }

              if ((eq_data->NbSmilePt < 0) ||
                  (eq_data->NbSmilePt > MAXNBDATE))
              {
                  DR_Error("Eq_Input_Dyn_WType6: Nb EQ smile param lines OWS out of range!");
                  goto RETURN;
              }

              for (i = 0; i < eq_data->NbSmilePt; i++)
              {
                  readerror = sscanf (EQSmileParamOWS[i], 
                                      "%ld\t%lf\t%lf\t%lf \n",
                                      &(eq_data->SmileDate[i]),
                                      &(eq_data->a1[i]),
                                      &(eq_data->a2[i]),
                                      &(eq_data->a3[i]));
                  eq_data->SmileDate[i] = IRDateFromYMDDate(eq_data->SmileDate[i]);
                  if (readerror != 4)
                  {          
                      DR_Error("Eq_Input_Dyn_WType6: Cannot read EQ smile params OWS");
                      goto RETURN;
                  }
              }
          }

          /* Open the eq data file */

          stream = fopen (FileName, "r");   

          if (stream == NULL)
          {
               DR_Error("Could not open file %s! (Hyb3_Eq_Input_Dyn_WType6)",
                          FileName);
               
               goto RETURN;
          
          }  /* if */                                                                
          

          /* Value date of vol curve */
          if (FindAndSkipComLine(stream,"value date","Hyb3_Eq_Input_Dyn_WType6", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          fgets (string, MAXBUFF, stream);  /* Ignore this input */


          /* Equity spot price */
          if (FindAndSkipComLine(stream, "spot level", "Hyb3_Eq_Input_Dyn_WType6", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%lf \n", &(eq_data->Spot));
          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read spot level file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType6)", FileName);
               goto RETURN;
          }
          
          /* Correlation between equity and IR in which eq is denominated */
          if (FindAndSkipComLine(stream,"correlation","Hyb3_Eq_Input_Dyn_WType6", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%lf \n", &(Correlation));
          if ((readerror == 0) || (readerror == EOF))
          {          
                     DR_Error("Could not read correlation in file %s! "
                                 "(Hyb3_Eq_Input_Dyn_WType6)", FileName);
                     
                     goto RETURN;
          }

          /* Number of points in vol curve */
          if (FindAndSkipComLine(stream,"number of vols","Hyb3_Eq_Input_Dyn_WType6",
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%d \n", &(eq_data->NbVol));
          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read number of vols in file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType6)", FileName);
               
               goto RETURN;
          }
          
          if (eq_data->NbVol > MAXNBDATE)
          {
               DR_Error("Number of vols exceeds maximum limit of %d in file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType6)", MAXNBDATE, FileName);
               
               goto RETURN;
          }

          /* Vol curve */
          if (FindAndSkipComLine(stream,"vol dates and rates",
                                 "Hyb3_Eq_Input_Dyn_WType6",FileName) == FAILURE)
          {          
               goto RETURN;
          }

          for (i = 1; i <= eq_data->NbVol; i++)
          {
               readerror = fscanf (stream, "%ld \t%lf \n",
                                          &(eq_data->VolDate[i]),
                                          &(eq_data->Vol[i]));

               eq_data->VolDate[i] = IRDateFromYMDDate(eq_data->VolDate[i]);
               if ((readerror == 0) || (readerror == EOF))
                {          
                     DR_Error("Could not read vol date and rate #%d in "
                                 "file %s! (Hyb3_Eq_Input_Dyn_WType6)", i+1, FileName);
                     
                     goto RETURN;
                }
               eq_data->Vol[i] /= 100.0;

          }  /* for i */

          /* Repeat the first values for interpolation */

          /* eq_data->VolDate[0] = Nxtday (eq_data->VolDate[1], (long) -1); */
          eq_data->VolDate[0] = ValueDate;
          eq_data->Vol[0]     = eq_data->Vol[1];

          /********************************************************/
          /* Construct flat correlation curve                     */
          /* although the data structure supports a term structure*/
          /* of correlation, we only use the first value          */
          /********************************************************/
          
          
          if(DomForFlag == 'D')
          {    /* Domestic equity */
                     
               eq_data->Rho[1][0] = Correlation;
          }
          else
          {    /* Foreign cups and foreign composite equity */
                     
               eq_data->Rho[0][0] = Correlation;
          }

          /* Check validity of input */
          if (Hyb3_Eq_Check_Dyn_WType6(eq_data,
                                          ValueDate,
                                          DomForFlag) == FAILURE)
          {          
                     goto RETURN;
          }

          status = SUCCESS;
RETURN:
          
          if (stream != NULL) fclose (stream);

          return (status);

}  /* Eq_Input_Dyn_W */


/*****  Hyb3_Eq_Check_Dyn_WType6	*************************************************/
/**
*   Check validity of static equity input.
*/
int Hyb3_Eq_Check_Dyn_WType6
          (EQ_DATA     *eq_data,       /**< (I) Structure of equity data      */
           long         ValueDate,     /**< (I) Value date                    */
           char         DomForFlag)    /**< (I) Equity denomination flag      */
{
          long
                     LastDate;
          int
                     CorrIdx,          /* Index of corr Eq vs. IrEquity     */
                     i,
                     status = FAILURE; /* Error status = FAILURE initially  */
                          
          if (eq_data->NbVol > MAXNBDATE)
          {
             DR_Error (  "Nb of equity vols exceeds maximum limit!");
             goto RETURN;
          }  /* if */

          if (eq_data->NbVol < 1)
          {
                     DR_Error (  "At least one equity volatility is required!");
                     goto RETURN;
          
          }  /* if */

          LastDate = Nxtmth (ValueDate, MAXMAT * 12, 1L);

          if (eq_data->VolDate[eq_data->NbVol] > LastDate)
          {
                     DR_Error ("Equity volatility point must be less than 60 years!");
                     goto RETURN;
          
          }  /* if */
          
          for (i = 1; i <= eq_data->NbVol; i++)
                     if (eq_data->VolDate[i] <= eq_data->VolDate[i-1])
                     {
                               DR_Error ("Equity volatility points must be entered in ascending order!");
                               goto RETURN;
          
                     }  /* if */                                                                
          
          for (i = 1; i <= eq_data->NbVol; i++)
                     if ((eq_data->Vol[i] < .000001) || (eq_data->Vol[i] > 9.99))
                     {
                               DR_Error ("Equity volatilities out of range!");
                               goto RETURN;
          
                     }  /* if */                                                                


          CorrIdx = (DomForFlag == 'D') ? 1 : 0;

          
          if ((eq_data->Rho[CorrIdx][0] < - MAXCORRELATION) || 
               (eq_data->Rho[CorrIdx][0] > MAXCORRELATION))
          {

              DR_Error ("Correlation out of range!");
              goto RETURN;
          
          }  /* if */


          status = SUCCESS;
          
          RETURN:
          
          return (status);

}  /* Hyb3_Eq_Check_Dyn_WType6 */


/*****  Hyb3_Eq_Input_Dyn_WType4	**************************************************/
/**
 *     Read dynamic equity input  for DR Wrappers of type 4 (i.e. equity is the 
 *     second asset and the EQ_DATA structure contains one correlation curve).
 *
 *     Correlation overwrite strings are dealt with here.
 *
 */
int Hyb3_Eq_Input_Dyn_WType4
                        (EQ_DATA   *eq_data,         /**< (O) Equity data         */
                         long       ValueDate,       /**< (I) Value date          */
                         char      *OwriteCorr,      /**< (I) String              */
                         char      *OwriteCutoff,    /**< (I) String              */
                         char      *FileName,        /**< (I) Including extension */
                         char       NbEQSmileOWS[MAXBUFF],
                         char       EQSmileParamOWS[MAXNBDATE][MAXBUFF])
{
          char
                     string[MAXBUFF];
          int
                     i,
                     readerror,                 /* Reading error status             */
                     status = FAILURE;          /* Error status = FAILURE initially */
          double
                     Correlation;
          FILE
                     *stream = NULL;


          if (strcmp (OwriteCutoff, "nil") == 0)    /* 1.Disallow cut off      */
          {

              eq_data->EqCutOffFlag = FALSE;
              eq_data->EqCutOffLast = FALSE;
              eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
              eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED;  

          }
          else if (strcmp(OwriteCutoff,"last")== 0) /* 2.Cut off at last level */
          {
              /* Cut off is allowed and there */
              /* is volatility cutoff data    */
              eq_data->EqCutOffFlag = TRUE;
              eq_data->EqCutOffLast = TRUE;
              eq_data->EqCutOffLevel = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
              eq_data->EqBootStrapMode = EQ_USE_LAST_LEVEL ;
          } 
          else                                         /* 3.Cut off at user level */
          {
              /* Cut off is allowed and there */
              /* is volatility cutoff data    */
              eq_data->EqCutOffFlag = TRUE;
              eq_data->EqCutOffLast = FALSE;
              eq_data->EqCutOffLevel = atof(OwriteCutoff);
              if (eq_data->EqCutOffLevel < TINY) /* The FAILURE return */
              {
                  DR_Error("Unable to convert EQ cut off value!\n");
                  goto RETURN;
              }
              eq_data->EqCutOffLevel /= 100.;
              eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;
          
          }

          /* this function is only used for environments lacking smile information */

          eq_data->NbSmilePt = 1;

          eq_data->SmileDate[0] = Nxtday(ValueDate, 5);

          eq_data->a1[0] = 0.0;
          eq_data->a2[0] = 0.0;
          eq_data->a3[0] = 0.0;

          /* this function is only used for environments lacking smile information */

          eq_data->NbSmilePt = 1;

          eq_data->SmileDate[0] = Nxtday(ValueDate, 5);

          eq_data->a1[0] = 0.0;
          eq_data->a2[0] = 0.0;
          eq_data->a3[0] = 0.0;

          /* but there may be overwrites */

          if (strstr(NbEQSmileOWS, "nil") != NULL)
          {
              if (strstr(EQSmileParamOWS[0], "nil") == NULL)
              {
                  DR_Error("Eq_Input_Dyn_WType4: invalid EQ Smile Params overwrites.\n");
                  goto RETURN;
              }
          }
          else
          {
              /* read the overwrite strings */
              readerror = sscanf (NbEQSmileOWS, 
                                  "%d \n", 
                                  &(eq_data->NbSmilePt));

              if (readerror != 1)
              {        
                  DR_Error("Eq_Input_Dyn_WType4: Cannot read Nb EQ Smile Param OWS");
                  goto RETURN;
              }

              if ((eq_data->NbSmilePt < 0) ||
                  (eq_data->NbSmilePt > MAXNBDATE))
              {
                  DR_Error("Eq_Input_Dyn_WType4: Nb EQ smile param lines OWS out of range!");
                  goto RETURN;
              }

              for (i = 0; i < eq_data->NbSmilePt; i++)
              {
                  readerror = sscanf (EQSmileParamOWS[i], 
                                      "%ld\t%lf\t%lf\t%lf \n",
                                      &(eq_data->SmileDate[i]),
                                      &(eq_data->a1[i]),
                                      &(eq_data->a2[i]),
                                      &(eq_data->a3[i]));

                  eq_data->SmileDate[i] = IRDateFromYMDDate(eq_data->SmileDate[i]);
                  if (readerror != 4)
                  {          
                      DR_Error("Eq_Input_Dyn_WType4: Cannot read EQ smile params OWS");
                      goto RETURN;
                  }
              }
          }

          /* Open the eq data file */

          stream = fopen (FileName, "r");

          if (stream == NULL)
          {
               DR_Error("Could not open file %s! (Hyb3_Eq_Input_Dyn_WType4)",
                         FileName);
               
               goto RETURN;
          
          }  /* if */                                                                
          

          /* Value date of vol curve */
          if (FindAndSkipComLine(stream,"value date","Hyb3_Eq_Input_Dyn_WType4", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          fgets (string, MAXBUFF, stream);  /* Ignore this input */


          /* Equity spot price */
          if (FindAndSkipComLine(stream, "spot level", "Hyb3_Eq_Input_Dyn_WType4", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%lf \n", &(eq_data->Spot));
          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read spot level file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType4)", FileName);
               
               goto RETURN;
          }
          
          /* Correlation between equity and IR in which eq is denominated */
          if (FindAndSkipComLine(stream,"correlation","Hyb3_Eq_Input_Dyn_WType4", 
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%lf \n", &(Correlation));
          if ((readerror == 0) || (readerror == EOF))
          {          
                     DR_Error("Could not read correlation in file %s! "
                                 "(Hyb3_Eq_Input_Dyn_WType4)", FileName);
                     
                     goto RETURN;
          }

          /* Number of points in vol curve */
          if (FindAndSkipComLine(stream,"number of vols","Hyb3_Eq_Input_Dyn_WType4",
                                         FileName) == FAILURE)
          {          
               goto RETURN;
          }
          readerror = fscanf (stream, "%d \n", &(eq_data->NbVol));
          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read number of vols in file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType4)", FileName);
               
               goto RETURN;
          }

          if (eq_data->NbVol > MAXNBDATE)
          {
               DR_Error("Number of vols exceeds maximum limit of %d in file %s! "
                            "(Hyb3_Eq_Input_Dyn_WType4)", MAXNBDATE, FileName);
               
               goto RETURN;
          }
          
          /* Vol curve */
          if (FindAndSkipComLine(stream,"vol dates and rates",
                                         "Hyb3_Eq_Input_Dyn_WType4",FileName) == FAILURE)
          {          
               goto RETURN;
          }

          for (i = 1; i <= eq_data->NbVol; i++)
          {
               readerror = fscanf (stream, "%ld \t%lf \n",
                                          &(eq_data->VolDate[i]),
                                          &(eq_data->Vol[i]));

               eq_data->VolDate[i] = IRDateFromYMDDate(eq_data->VolDate[i]);
               if ((readerror == 0) || (readerror == EOF))
                {          
                     DR_Error("Could not read vol date and rate #%d in "
                                 "file %s! (Hyb3_Eq_Input_Dyn_WType4)", i+1, FileName);
                     
                     goto RETURN;
                }
               eq_data->Vol[i] /= 100.0;

          }  /* for i */

          /* eq_data->VolDate[0] = Nxtday (eq_data->VolDate[1], (long) -1); */

          eq_data->VolDate[0] = ValueDate;
          eq_data->Vol[0]     = eq_data->Vol[1];
          
          /* Construct flat correlation curve    */
          /* although the data structure supports*/
          /* a term-strucutre of correlation     */
          /* we only use the firts elelment      */
          
          
          eq_data->Rho[0][0] = Correlation;

          /* Overwrite it if necessary */     
          if (strcmp (OwriteCorr, "nil"))
          {                            
              eq_data->Rho[0][0] = atof(OwriteCorr);
          }

          /* Check validity of input */
          if (Hyb3_Eq_Check_Dyn_WType4(eq_data,
                                  ValueDate) == FAILURE)
          {          
              goto RETURN;
          }
          

          status = SUCCESS;

RETURN:
          
          if (stream != NULL) fclose (stream);

          return (status);

}  /* Hyb3_Eq_Input_Dyn_WType4  */





/*****  Hyb3_Eq_Check_Dyn_WType4	*************************************************/
/**
*   Check validity of static equity input.
*/
int	Hyb3_Eq_Check_Dyn_WType4
             (EQ_DATA      *eq_data,       /**< (I) Structure of equity data  */
              long          ValueDate)     /**< (I) Value date                */
{
        long
            LastDate;
        int
            i,
            status = FAILURE;         /* Error status = FAILURE initially */
                        
                        
        if (eq_data->NbVol > MAXNBDATE)
        {
            DR_Error ("Nb of equity vols exceeds maximum limit!");
            goto RETURN;
        }

        if (eq_data->NbVol < 1)
        {
            DR_Error (  "At least one equity volatility is required!");
            goto RETURN;
        
        }  /* if */

        LastDate = Nxtmth (ValueDate, MAXMAT * 12, 1L);

        if (eq_data->VolDate[eq_data->NbVol] > LastDate)
        {
            DR_Error ("Equity volatility point must be less than 60 years!");
            goto RETURN;
        
        }  /* if */
        
        for (i = 1; i <= eq_data->NbVol; i++)
            if (eq_data->VolDate[i] <= eq_data->VolDate[i-1])
            {
                DR_Error ("Equity volatility points must be entered in ascending order!");
                goto RETURN;
        
            }  /* if */                                                                
        
        for (i = 1; i <= eq_data->NbVol; i++)
            if ((eq_data->Vol[i] < .000001) || (eq_data->Vol[i] > 9.99))
            {
                DR_Error ("Equity volatilities out of range!");
                goto RETURN;
        
            }  /* if */                                                                


        
        
        if ((eq_data->Rho[0][0] < - MAXCORRELATION) || 
             (eq_data->Rho[0][0] > MAXCORRELATION))
        {
             DR_Error ("Correlation out of range!");
             goto RETURN;
        
        }  /* if */

        status = SUCCESS;
        
        RETURN:
        
        return (status);

}  /* Hyb3_Eq_Check_Dyn_WType4 */



/*****  Hyb3_Eq_Input_Sta_W  *****************************************************/
/**
*   Read static equity input for wrappers and check validity of input.
*/
int  Hyb3_Eq_Input_Sta_W(
                    EQ_DATA *eq_data,       /**< (O) Structure of equity data */
                    long      ValueDate,    /**< (I) Value date               */
                    char     *FileName)     /**< (I) File name including ext  */
{
          long
                     bDate;
          int
                     i,
                     numPeriods,
                     readerror,               /* Reading error status             */
                     status = FAILURE;        /* Error status = FAILURE initially */
          char
                     delimiter[] = " ,%\t\n", /* file delimiters                  */
                     period,
                     *token,                  /*                                  */
                     string[MAXBUFF];
          FILE
                     *stream = NULL;


          stream = fopen (FileName, "r");       /* Open the eq data file */

          if (stream == NULL)
          {
                     DR_Error("Could not open file %s! (Eq_Input_Sta_W)", FileName);  
                     goto RETURN;

          }

          if (FindAndSkipComLine(   stream,
                                    "number of borrowing curve points",
                                    "Eq_Input_Sta_W", FileName) == FAILURE)
          {
                     goto RETURN;
          }

          if (FindAndSkipComLine(   stream,
                                    "number of borrowing curve points",
                                    "Eq_Input_Sta_W", FileName) == FAILURE)
          {
                     goto RETURN;
          }

          readerror = fscanf (stream, "%d \n", &(eq_data->NbBorrow));

          if ((readerror == 0) || (readerror == EOF))
          {
             DR_Error(  "Could not read number of borrowing curve points in file %s! "
                        "(Eq_Input_Sta_W)", FileName);
             goto RETURN;
          }

          if (eq_data->NbBorrow > MAXNBDATE)
          {
             DR_Error(  "Nb of borrow curve points exceeds maximum limit of %d in file %s! "
                        "(Eq_Input_Sta_W)", 
                        MAXNBDATE,   
                        FileName);
             goto RETURN;
          }

          for (i = 0; i < eq_data->NbBorrow; i++)
          {
                readerror = fscanf (stream, "%s \t%lf \n",
                                    string, &(eq_data->Borrow[i]));

                if (readerror != 2)
                {
                    DR_Error(   "Could not read label and borrowing rate #%d in file %s! "
                                "(Eq_Input_Sta_W)", i+1, FileName);
                    goto RETURN;
                }

                eq_data->Borrow[i] /= 100;

                period     = string[strlen(string)-1];
                bDate      = atol(string);
                numPeriods = atoi(string);

                if (!isalpha(period))
                {
                    if (Dateok(bDate))
                    {
                        DR_Error(   "Could not read label #%d in file %s! "
                                    "(Eq_Input_Sta_W)", i+1, FileName);
                        goto RETURN;
                    }
                    else
                    {
                        eq_data->BorrowDate[i] = bDate;
                    }
                }
                else
                {
                    if (period == 'I') /* others are checked in DrDateFwdAny */
                    {
                        DR_Error(   "Could not read label #%d in file %s! "
                                    "(Eq_Input_Sta_W)", i+1, FileName);
                        goto RETURN;
                    }

                    if (period == 'Y')
                    {
                        period = 'A';
                    }

                    if(DrDateFwdAny(ValueDate,
                                    numPeriods,
                                    period,
                                    'F',
                                    &(eq_data->BorrowDate[i])) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
          }  /* for i */

          if (FindAndSkipComLine (stream, "number of dividends",
                                  "Eq_Input_Sta_W", FileName) == FAILURE)
          {
                     goto RETURN;
          }

          readerror = fscanf (stream, "%d \n", &(eq_data->NbFwd));
          if ((readerror == 0) || (readerror == EOF))
          {          
                     DR_Error("Could not read number of dividends in file %s! "
                              "(Eq_Input_Sta_W)", FileName);
                     
                     goto RETURN;
          }

          if (eq_data->NbFwd > MAXNBEQDIV)
          {          
                     DR_Error("Given number of divs (%d) exceeds maximum allowed (%d)! "
                              "(Eq_Input_Sta_W)", eq_data->NbFwd, MAXNBEQDIV);
                     
                     goto RETURN;
          }
          
          for (i = 0; i < eq_data->NbFwd; i++)
          {
                     fgets (string, MAXBUFF, stream); /* Can't use scanf but have to use strtok */

                     token = strtok (string, delimiter);

                     if (token == NULL)
                     {
                               DR_Error("Could not read dividend date #%d in file %s! "
                                        "(Eq_Input_Sta_W)", i+1, FileName);
                               goto RETURN;
                                                                
                     }  /* if */
                                                     
                     eq_data->FwdDate[i] = atol (token);

                     token = strtok (NULL, delimiter);
                                                                
                     if (token == NULL)
                     {
                               DR_Error("Could not read dividend rate #%d in file %s! "
                                        "(Eq_Input_Sta_W)", i+1, FileName);
                               
                               goto RETURN;
                                                                
                     }  /* if */
                                                     
                     eq_data->Fwd[i] = atof (token);

                     token = strtok (NULL, delimiter);
                                                                
                     if (token == NULL)
                     {
                               DR_Error("Could not read dividend type #%d in file %s! "
                                        "(Eq_Input_Sta_W)", i+1, FileName);
                               
                               goto RETURN;
                                                                
                     }  /* if */
                                                     
                     eq_data->FwdType[i] = token[0];
                     if ((eq_data->FwdType[i] == 'Y') ||
                         (eq_data->FwdType[i] == 'C'))
                     {
                          eq_data->Fwd[i] /= 100.0;
                     }

          }  /* for i */
          
          if (FindAndSkipComLine (stream, "number of settlement points",
                                  "Eq_Input_Sta_W", FileName) == FAILURE)
          {          
                     goto RETURN;
          }

          readerror = fscanf (stream, "%d \n", &(eq_data->NbSettle));
          if ((readerror == 0) || (readerror == EOF))
          {          
                     DR_Error("Could not read number of settlement points in file %s! "
                              "(Eq_Input_Sta_W)", FileName);
                     goto RETURN;
          }

          if (eq_data->NbSettle > MAXNBDATE)
          {
             DR_Error("number of settlement points exceeds limit of %d in file %s! "
                                "(Eq_Input_Sta_W)",MAXNBDATE, FileName);
             
             goto RETURN;
          }

          if ((eq_data->NbSettle == 0)||(eq_data->NbSettle == -1))  /* No settlement date list: rolling settlement */
          {          
                     eq_data->SettleType = 'R';

                     readerror = fscanf (stream, "%d \n", &(eq_data->NbSettle));
                     if ((readerror == 0) || (readerror == EOF))
                     {          
                               DR_Error("Could not read rolling settlement days in file %s! "
                                        "(Hyb3_Eq_Input_Sta_W)", FileName);
                               
                               goto RETURN;
                     }
          }
          else /* Fixed settlement */
          {
                     eq_data->SettleType = 'F';

                     for (i = 0; i < eq_data->NbSettle; i++)
                     {
                               readerror = fscanf (stream, "%ld \t%ld \n",
                                          &(eq_data->LastTrading[i]),
                                          &(eq_data->SettleDate[i]));

                               eq_data->LastTrading[i] = IRDateFromYMDDate(eq_data->LastTrading[i]);
                               eq_data->SettleDate[i] = IRDateFromYMDDate(eq_data->SettleDate[i]);
                               if ((readerror == 0) || (readerror == EOF))
                                {          
                                          DR_Error("Could not read last trading and settlement dates #%d in file %s! (Hyb3_Eq_Input_Sta_W)", i+1, FileName);
                                          
                                          goto RETURN;
                                }
                     }  /* for i */
          }  /* if then else */
          

          if (Hyb3_Eq_Check_Sta_W (eq_data) == FAILURE) /* Check validity of input */
          {          
                     goto RETURN;
          }
          

          status = SUCCESS;
          
RETURN:
          
          if (stream != NULL) fclose (stream);

          return (status);

}  /* Hyb3_Eq_Input_Sta_W */


/*****  Hyb3_Eq_Check_Sta_W  *****************************************************/
/**
*   Check validity of dynamic equity input.
*/
int Hyb3_Eq_Check_Sta_W (EQ_DATA *eq_data)  /**< (I) Structure of equity data	*/
{
          int
                     i,
                     status = FAILURE;            /* Error status = FAILURE initially */
                          
                          
          if (eq_data->NbBorrow > MAXNBDATE)
          {
             DR_Error ("Nb of borrowing curve points exceeds max limit!");
             goto RETURN;
          }

          if (eq_data->NbFwd > MAXNBEQDIV)
          {          
             DR_Error ("number of divs exceeds max limit!");
             goto RETURN;
          }

          if (eq_data->NbSettle > MAXNBDATE)
          {
             DR_Error ("number of settlement points exceeds limit!");
             goto RETURN;
          }

          if (eq_data->NbFwd < 1)
          {
                     DR_Error ("At least one dividend point is required!");
                     goto RETURN;
          
          }  /* if */

          for (i = 1; i < eq_data->NbFwd; i++)
                     if (eq_data->FwdDate[i] <= eq_data->FwdDate[i-1])
                     {
                               DR_Error ("Dividend dates must be entered in ascending order!");
                               goto RETURN;
          
                     }  /* if */ 
          
          for (i = 0; i < eq_data->NbFwd; i++)
                     if (     (eq_data->FwdType[i] != 'D')
                          &&  (eq_data->FwdType[i] != 'Y')
                          &&  (eq_data->FwdType[i] != 'C'))
                     {
                               DR_Error ("Dividend type must be 'D', 'Y' or 'C'!");
                               goto RETURN;
          
                     }  /* if */ 

          i = 0;           /* Find first continuous dividend */
          while ((i < eq_data->NbFwd) && (eq_data->FwdType[i] != 'C'))
                     i++;

          for (; i < eq_data->NbFwd; i++)
                     if (eq_data->FwdType[i] != 'C')
                     {
                               DR_Error ("Can't have discrete dividends after a continuous one!");
                               goto RETURN;
          
                     }  /* if */ 

          if (eq_data->SettleType == 'R')
          {
                     if (eq_data->NbSettle < 0)
                     {
                               DR_Error ("Rolling settlement period should be positive!");
                               goto RETURN;
          
                     }  /* if */
          }
          else
          {
                     if (eq_data->NbSettle < 1)
                     {
                               DR_Error ("At least one settlement date is required!");
                               goto RETURN;
          
                     }  /* if */

                     for (i = 1; i < eq_data->NbSettle; i++)
                     {
                               if (eq_data->LastTrading[i] <= eq_data->LastTrading[i-1])
                                {
                                          DR_Error ("Last trading dates must be entered in ascending order!");
                                          goto RETURN;
          
                                }  /* if */                                                                
                                
                               if (eq_data->SettleDate[i] <= eq_data->SettleDate[i-1])
                                {
                                          DR_Error ("Settlement dates must be entered in ascending order!");
                                          goto RETURN;
          
                                }  /* if */                                                                
                     }  /* for i */
                                
                     for (i = 0; i < eq_data->NbSettle; i++)
                               if (eq_data->LastTrading[i] > eq_data->SettleDate[i])
                                {
                                          DR_Error ("Last trading date must be before settlement date!");
                                          goto RETURN;
          
                                }  /* if */                                                                
          }  /* if then else */
          
          
          status = SUCCESS;
          
          RETURN:
          
          return (status);

}  /* Hyb3_Eq_Check_Sta_W */




/*****  Hyb3_Fx_Input_W_WithSmile  *********************************************************/
/**
    Read fx input for wrappers and check validity of input.
    Composite vols are input with a 1 OFFSET, wheras SpotVols are input
    with 0 OFFSET                                                           
 ****************************************************************************/
int  Hyb3_Fx_Input_W_WithSmile(
                          FX_DATA*    fx_data,                 /**< (O) Fx data            */
                          char const  OverWriteString[MAXBUFF], /**< (I) Owrite str         */
                          char const* VolFileName,             /**< (I) FXVols file name   */
                          char const* SmileFileName,           /**< (I) FX smile file name */
                          char const  NbFXSmileOWS[MAXBUFF],    /**< (I)                    */
                          char const  FXSmileParamOWS[MAXNBDATE][MAXBUFF]   /**< (I) */)
{

     long
          Year[2],
          Month[2],
          Day[2];
     int
          i,
          readerror,                 /* Reading error status             */
          status = FAILURE;          /* Error status = FAILURE initially */
     char
          string[MAXBUFF];
     FILE
         *stream        = NULL,
         *streamFXSmile = NULL; 
     


     stream = fopen (VolFileName, "r");           /* Open the fx data file  */

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Hyb3_Fx_Input_W)", VolFileName);
          
          goto RETURN;
          
     }  /* if */
          
     fgets  (string, 80, stream);    /* Read the comment line */
     readerror = fscanf (stream, "%ld \n", &(fx_data->Today));
     fx_data->Today = IRDateFromYMDDate(fx_data->Today);
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", VolFileName);
          
          goto RETURN;
     }
          
     fx_data->ValueDate = fx_data->Today;
     fx_data->SpotDays  = 0;
     
     fgets  (string, 80, stream);
     readerror = fscanf (stream, "%lf \n", &(fx_data->Spot));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", VolFileName);
          
          goto RETURN;
     }
     
     fgets ( string, 80, stream);  /* Ignore this input */
     fgets ( string, 80, stream);
     
     fgets ( string, 80, stream);
     readerror = fscanf (stream, "%d \n", &(fx_data->NbVol));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", VolFileName);
          
          goto RETURN;
     }

     if (fx_data->NbVol > MAXNBDATE)
     {
          DR_Error("Nb of vols exceeds max limit of %d in file %s! "
                             "(Hyb3_Fx_Input_W)", MAXNBDATE, VolFileName);
          
          goto RETURN;
     }
     
     Dsplit(fx_data->ValueDate, /* Split value date into month, day and year */
              &(Month[0]), 
              &(Day[0]), 
              &(Year[0]));

     fgets ( string, 80, stream);
     for (i = 1; i <= fx_data->NbVol; i++)  /* Composite vols -> 1 offset*/
     {
          readerror = fscanf (stream, "%ld \t%lf \n",
          &(fx_data->VolDate[i]),
          &(fx_data->FxVol[i]));
          fx_data->VolDate[i] = IRDateFromYMDDate(fx_data->VolDate[i]);

          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read file %s! "
                   "(Hyb3_Fx_Input_W)", VolFileName);
               
               goto RETURN;
          }
          
          fx_data->FxVol[i] /= 100.0;
     }  /* for i */

     /* Repeat the first values for interpolation of pseudocomposite vols */
     /* in Get_TreeSpotVols                                               */

     if (fx_data->NbVol > 0)
     {
        fx_data->VolDate[0] = fx_data->ValueDate;
        fx_data->FxVol[0]   = fx_data->FxVol[1];
     }
     

    /* Input Spot Vol section */
    fgets  (string, 80, stream);    /* comment line */
    readerror = fscanf (stream, "%d \n", &(fx_data->NbInpSpotVol));
    if (readerror != 1) 
    {          
         DR_Error("Could not read Nb of Input Spot Vol points "
                   "in file %s! (Hyb3_Fx_Input_W)", VolFileName);
         
         goto RETURN;
    }
     
    if (fx_data->NbInpSpotVol > MAXNBDATE)
    {
         DR_Error("Nb of Inp Spot Vol exceeds max limit of %d "
                   "in file %s! (Hyb3_Fx_Input_W)", MAXNBDATE, VolFileName);
         
         goto RETURN;
    }

    if (fx_data->NbInpSpotVol > 0)
    {
       fgets  (string, 80, stream);    /* comment line */
       for (i = 0; i < fx_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
       {
             readerror = fscanf (stream, "%ld \t%lf \n",
                                 &(fx_data->InpSpotVolDate[i]),
                                 &(fx_data->InpSpotVol[i]));
             fx_data->InpSpotVolDate[i] = IRDateFromYMDDate(fx_data->InpSpotVolDate[i]);

           if (readerror != 2)
           {          
                  DR_Error("Could not read input Spot Vol"
                      " dates/values in file %s! (Hyb3_Fx_Input_W)", VolFileName);
               
               goto RETURN;
           }
          
           fx_data->InpSpotVol[i] /= 100.0;
       }  /* for i */
    }


    /* deal with possible overlap between composite and spot vols */
    if ((fx_data->NbInpSpotVol > 0) && (fx_data->NbVol > 0))
    {
        long Idx;
        Idx = GetDLOffset(fx_data->NbVol,
                          &(fx_data->VolDate[1]),
                          fx_data->InpSpotVolDate[0],
                          CbkHIGHER);
        /*******************************************************/
        /* if all composite dates are < spotvoldate[0]         */
        /* then Idx = -999, in which case NbVol does not need  */
        /* to be changed.                                      */
        /*******************************************************/
        if (Idx >= 0L)
        {
            fx_data->NbVol = Idx;
        }        
    }
    
    if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
     {

        fx_data->FxCutOffFlag = FALSE;
        fx_data->FxCutOffLast = FALSE;
        fx_data->FxCutOffLevel   = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        fx_data->FxBootStrapMode = FX_NO_FAILURE_ALLOWED; 

     }
     else if (strcmp(OverWriteString,"last")== 0) /* 2.Cut off at last level */
     {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        fx_data->FxCutOffFlag = TRUE;
        fx_data->FxCutOffLast = TRUE;
        fx_data->FxCutOffLevel = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
        fx_data->FxBootStrapMode = FX_USE_LAST_LEVEL ;
     } 
     else                                         /* 3.Cut off at user level */
     {
          /* Cut off is allowed and there */
          /* is volatility cutoff data    */
          fx_data->FxCutOffFlag = TRUE;
          fx_data->FxCutOffLast = FALSE;
          fx_data->FxCutOffLevel = atof(OverWriteString);
          if (fabs (fx_data->FxCutOffLevel) < TINY ) /* The FAILURE return */
          {
              DR_Error("Unable to convert FX cut off value!\n");
              goto RETURN;
          }
          fx_data->FxCutOffLevel /= 100.;
          fx_data->FxBootStrapMode = FX_CONSTANT_SPOT_VOL;
          
     }
     

    /* Read Smile data file */
    

    if (strstr(NbFXSmileOWS, "nil") != NULL)
    {
        if (strstr(FXSmileParamOWS[0], "nil") == NULL)
        {
            DR_Error("Hyb3_Fx_Input_W: invalid FX Smile Params input\n");
            goto RETURN;
        }

        /* read from environment file.  If environment smile parameters do not exist,
           the model will default to smile off, by setting a1=a2=a3=0 at the valuedate */
        streamFXSmile = fopen (SmileFileName, "r");
        if (streamFXSmile == NULL)
        {
            fx_data->NbSmilePt = 1;

            /* add arbitary number of days after the value date to ensure the smile date is
               > today as in some cases valueDate = Today. - chose to use 5 days but could have
               use any number */
            fx_data->SmileDate[0] = Nxtday(fx_data->ValueDate, 5);
            fx_data->a1[0] = 0.0;
            fx_data->a2[0] = 0.0;
            fx_data->a3[0] = 0.0;
        }
        else
        {
            /* read Nb FX Smile Param lines */
            if (FindAndSkipComLine (streamFXSmile, 
                                    "Nb FX Param Lines",
                                    "Hyb3_Fx_Input_W", 
                                    SmileFileName) == FAILURE) goto RETURN;

            readerror = fscanf (streamFXSmile, "%d \n", &(fx_data->NbSmilePt));
            if (readerror != 1)
            {
                DR_Error("Hyb3_Fx_Input_W: Cannot read Nb FX smile param lines");
                goto RETURN;
            }

            if ((fx_data->NbSmilePt < 0) ||
                (fx_data->NbSmilePt > MAXNBDATE))
            {
                DR_Error("Hyb3_Fx_Input_W: Nb FX smile param lines out of range!");
                goto RETURN;
            }

            if (FindAndSkipComLine (streamFXSmile, 
                                   "FX Smile Params",
                                   "Hyb3_Fx_Input_W", 
                                   SmileFileName) == FAILURE) goto RETURN;

            for (i=0; i < fx_data->NbSmilePt; i++)
            {
                readerror = fscanf(streamFXSmile, "%ld\t%lf\t%lf\t%lf \n",
                                   &(fx_data->SmileDate[i]),
                                   &(fx_data->a1[i]),
                                   &(fx_data->a2[i]),
                                   &(fx_data->a3[i]));

                fx_data->SmileDate[i] = IRDateFromYMDDate(fx_data->SmileDate[i]);
                if (readerror != 4)
                {          
                    DR_Error("Hyb3_Fx_Input_W: Cannot read FX smile params");
                    goto RETURN;
                }
            } 
            fclose( streamFXSmile );
        }
    }
    else
    {
        /* read the overwrite strings */
        readerror = sscanf (NbFXSmileOWS, 
                            "%d \n", 
                            &(fx_data->NbSmilePt));

        if (readerror != 1)
        {        
            DR_Error("Hyb3_Fx_Input_W: Cannot read Nb FX Smile Param OWS");
            goto RETURN;
        }

        if ((fx_data->NbSmilePt < 0) ||
            (fx_data->NbSmilePt > MAXNBDATE))
        {
            DR_Error("Hyb3_Fx_Input_W: Nb FX smile param lines OWS out of range!");
            goto RETURN;
        }

        for (i = 0; i < fx_data->NbSmilePt; i++)
        {
            readerror = sscanf (FXSmileParamOWS[i], 
                                "%ld\t%lf\t%lf\t%lf \n",
                                &(fx_data->SmileDate[i]),
                                &(fx_data->a1[i]),
                                &(fx_data->a2[i]),
                                &(fx_data->a3[i]));

            fx_data->SmileDate[i] = IRDateFromYMDDate(fx_data->SmileDate[i]);
            if (readerror != 4)
            {          
                DR_Error("Hyb3_Fx_Input_W: Cannot read FX smile params OWS");
                goto RETURN;
            }
        }
    }





     /* Check validity of input */ 
     if (Hyb3_Fx_Check_W (fx_data) == FAILURE)
     {
          goto RETURN;
     }

     status = SUCCESS;
          
     RETURN:

     if (stream != NULL) fclose (stream);
          
     return (status);

}  /* Hyb3_Fx_Input_W_WithSmile */




/*****  Hyb3_Fx_Input_W*********************************************************/
/**
    Read fx input for wrappers and check validity of input.
    Composite vols are input with a 1 OFFSET, wheras SpotVols are input
    with 0 OFFSET                                                             
 ****************************************************************************/
int  Hyb3_Fx_Input_W(
                FX_DATA   *fx_data,                 /**< (O) Fx data            */
                char      OverWriteString[MAXBUFF], /**< (I) Owrite str         */
                char      *FileName)             /**< (I) FXVols file name   */
{
     long
          Year[2],
          Month[2],
          Day[2];
     int
          i,
          readerror,                 /* Reading error status             */
          status = FAILURE;          /* Error status = FAILURE initially */
     char
          string[MAXBUFF];
     FILE
         *stream        = NULL;

     


     stream = fopen (FileName, "r");           /* Open the fx data file  */

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Hyb3_Fx_Input_W)", FileName);
          
          goto RETURN;
          
     }  /* if */
          
     fgets  (string, 80, stream);    /* Read the comment line */
     readerror = fscanf (stream, "%ld \n", &(fx_data->Today));
     fx_data->Today = IRDateFromYMDDate(fx_data->Today);
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", FileName);
          
          goto RETURN;
     }
          
     fx_data->ValueDate = fx_data->Today;
     fx_data->SpotDays  = 0;
     
     fgets  (string, 80, stream);
     readerror = fscanf (stream, "%lf \n", &(fx_data->Spot));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", FileName);
          
          goto RETURN;
     }
     
     fgets ( string, 80, stream);  /* Ignore this input */
     fgets ( string, 80, stream);
     
     fgets ( string, 80, stream);
     readerror = fscanf (stream, "%d \n", &(fx_data->NbVol));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Hyb3_Fx_Input_W)", FileName);
          
          goto RETURN;
     }

     if (fx_data->NbVol > MAXNBDATE)
     {
          DR_Error("Nb of vols exceeds max limit of %d in file %s! "
                             "(Hyb3_Fx_Input_W)", MAXNBDATE, FileName);
          
          goto RETURN;
     }
     
     Dsplit(fx_data->ValueDate, /* Split value date into month, day and year */
              &(Month[0]), 
              &(Day[0]), 
              &(Year[0]));

     fgets ( string, 80, stream);
     for (i = 1; i <= fx_data->NbVol; i++)  /* Composite vols -> 1 offset*/
     {
          readerror = fscanf (stream, "%ld \t%lf \n",
          &(fx_data->VolDate[i]),
          &(fx_data->FxVol[i]));
          fx_data->VolDate[i] = IRDateFromYMDDate(fx_data->VolDate[i]);

          if ((readerror == 0) || (readerror == EOF))
          {          
               DR_Error("Could not read file %s! "
                   "(Hyb3_Fx_Input_W)", FileName);
               
               goto RETURN;
          }
          
          fx_data->FxVol[i] /= 100.0;
     }  /* for i */

     /* Repeat the first values for interpolation of pseudocomposite vols */
     /* in Get_TreeSpotVols                                               */

     if (fx_data->NbVol > 0)
     {
        fx_data->VolDate[0] = fx_data->ValueDate;
        fx_data->FxVol[0]   = fx_data->FxVol[1];
     }
     

    /* Input Spot Vol section */
    fgets  (string, 80, stream);    /* comment line */
    readerror = fscanf (stream, "%d \n", &(fx_data->NbInpSpotVol));
    if (readerror != 1) 
    {          
         DR_Error("Could not read Nb of Input Spot Vol points "
                   "in file %s! (Hyb3_Fx_Input_W)", FileName);
         
         goto RETURN;
    }
     
    if (fx_data->NbInpSpotVol > MAXNBDATE)
    {
         DR_Error("Nb of Inp Spot Vol exceeds max limit of %d "
                   "in file %s! (Hyb3_Fx_Input_W)", MAXNBDATE, FileName);
         
         goto RETURN;
    }

    if (fx_data->NbInpSpotVol > 0)
    {
       fgets  (string, 80, stream);    /* comment line */
       for (i = 0; i < fx_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
       {
             readerror = fscanf (stream, "%ld \t%lf \n",
                                 &(fx_data->InpSpotVolDate[i]),
                                 &(fx_data->InpSpotVol[i]));

             fx_data->InpSpotVolDate[i] = IRDateFromYMDDate(fx_data->InpSpotVolDate[i]);
           if (readerror != 2)
           {          
                  DR_Error("Could not read input Spot Vol"
                      " dates/values in file %s! (Hyb3_Fx_Input_W)", FileName);
               
               goto RETURN;
           }
          
           fx_data->InpSpotVol[i] /= 100.0;
       }  /* for i */
    }
    fclose (stream);


    /* deal with possible overlap between composite and spot vols */
    if ((fx_data->NbInpSpotVol > 0) && (fx_data->NbVol > 0))
    {
        long Idx;
        Idx = GetDLOffset(fx_data->NbVol,
                          &(fx_data->VolDate[1]),
                          fx_data->InpSpotVolDate[0],
                          CbkHIGHER);
        /*******************************************************/
        /* if all composite dates are < spotvoldate[0]         */
        /* then Idx = -999, in which case NbVol does not need  */
        /* to be changed.                                      */
        /*******************************************************/
        if (Idx >= 0L)
        {
            fx_data->NbVol = Idx;
        }        
    }
    
    if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
     {

        fx_data->FxCutOffFlag = FALSE;
        fx_data->FxCutOffLast = FALSE;
        fx_data->FxCutOffLevel   = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        fx_data->FxBootStrapMode = FX_NO_FAILURE_ALLOWED; 

     }
     else if (strcmp(OverWriteString,"last")== 0) /* 2.Cut off at last level */
     {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        fx_data->FxCutOffFlag = TRUE;
        fx_data->FxCutOffLast = TRUE;
        fx_data->FxCutOffLevel = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
        fx_data->FxBootStrapMode = FX_USE_LAST_LEVEL ;
     } 
     else                                         /* 3.Cut off at user level */
     {
          /* Cut off is allowed and there */
          /* is volatility cutoff data    */
          fx_data->FxCutOffFlag = TRUE;
          fx_data->FxCutOffLast = FALSE;
          fx_data->FxCutOffLevel = atof(OverWriteString);
          if (fabs(fx_data->FxCutOffLevel) < TINY) /* The FAILURE return */
          {
              DR_Error("Unable to convert FX cut off value!\n");
              goto RETURN;
          }
          fx_data->FxCutOffLevel /= 100.;
          fx_data->FxBootStrapMode = FX_CONSTANT_SPOT_VOL;
          
     }
     

    /* Hard-Code Smile parameters for back-ward compatibility */
    fx_data->NbSmilePt = 1L;

    if(DrDateFwdAny(fx_data->ValueDate,
                    1,
                    'D',
                    'F',
                    &(fx_data->SmileDate[0])) == FAILURE)
    {
        goto RETURN;
    }
    fx_data->a1[0] = 0.0;
    fx_data->a2[0] = 0.0;
    fx_data->a3[0] = 0.0;


     /* Check validity of input */ 
     if (Hyb3_Fx_Check_W (fx_data) == FAILURE)
     {
          goto RETURN;
     }

     status = SUCCESS;
          
     RETURN:
          
     return (status);

}  /* Hyb3_Fx_Input_W */


/*****  Hyb3_Fx_Check_W  *********************************************************/
/**
*	Read term structure input and check validity of input.
*/
int    Hyb3_Fx_Check_W(FX_DATA     *fx_data)  /**< (I) Structure of fx data */
{
     int
     i,
     status = FAILURE;     /* Error status = FAILURE initially */
                
                
     if (Dateok(fx_data->ValueDate))
     {
          DR_Error ("Incorrect format for fx value date!");
          goto RETURN;
     }  /* if */
     
     if (Dateok(fx_data->Today))
     {
         DR_Error("incorrect format for fx.Today date\n");
         goto RETURN;
     }

     if (fx_data->NbVol < 0) 
     {
         DR_Error("Nb of composite vols must be >= 0\n");
         goto RETURN;
     }

     if (fx_data->NbVol > MAXNBDATE)
     {
          DR_Error ("Nb of vols exceeds max limit\n");
          goto RETURN;
     }

     if (fx_data->NbInpSpotVol < 0)
     {
         DR_Error(" Nb of input spot vols must be >= 0 \n");
         goto RETURN;
     }

     if (fx_data->NbInpSpotVol > MAXNBDATE)
     {
          DR_Error ("Nb of input spot vols exceeds max limit\n");
          goto RETURN;
     }
     if((fx_data->NbInpSpotVol == 0) &&
         (fx_data->NbVol == 0))
     {
         DR_Error("Must have at least one fx vol point (composite or spot)\n");
         goto RETURN;
     }

    
     for (i = 1; i <= fx_data->NbVol; i++)
     {    
          if (Dateok(fx_data->VolDate[i]))
          {
              DR_Error("Invalid format for fx composite vol dates\n");
              goto RETURN;
          }
     }
     
     for (i = 1; i <= fx_data->NbVol; i++)
     {              
         if (fx_data->VolDate[i] <= fx_data->VolDate[i-1])
         {
              DR_Error ("FX volatility points must be "
                        " entered in ascending order!");
              goto RETURN;
         }  /* if */
     }
     
     if (fx_data->NbVol > 0)
     {
        if (fx_data->VolDate[0] != fx_data->ValueDate)
        {
             DR_Error("CompositeVolDate[0] must be same as fx.value date\n");
            goto RETURN;
        }

     
        if (fabs(fx_data->FxVol[0] - fx_data->FxVol[1]) > TINY)
        {
             DR_Error("The internal compositevol[0] must be the same"
                        " as compositevol[1] \n");
            goto RETURN;             
        }
     }

     for (i = 0 ; i < fx_data->NbInpSpotVol ; i++)
     {
         if(Dateok(fx_data->InpSpotVolDate[i]))
         {
             DR_Error("incorrect format for Input spot vol dates\n");
             goto RETURN;
         }
     }

     for (i = 0 ; i < (fx_data->NbInpSpotVol - 1) ; i++)
     {
         if(fx_data->InpSpotVolDate[i] >= fx_data->InpSpotVolDate[i+1])
         {
             DR_Error("Spot Vol dates must be entered in ascending order\n");
             goto RETURN;
         }
     }
    
     if (fx_data->NbInpSpotVol > 0) 
     {
         if(fx_data->InpSpotVolDate[0] <= fx_data->ValueDate)
         {
            DR_Error("First Hyb3_SpotVol date must be > fx Value date\n");
            goto RETURN;
         }
     }

     for (i = 0 ; i < fx_data->NbInpSpotVol ; i++)
     {
         if((fx_data->InpSpotVol[i] < 0.000001) || 
            (fx_data->InpSpotVol[i] > 9.99))
         {
             DR_Error("Spot Vols not entered correctly \n");
             goto RETURN;
         }
     }

     if ((fx_data->NbVol > 0) && (fx_data->NbInpSpotVol > 0))
     {
         if (fx_data->InpSpotVolDate[0] <= fx_data->VolDate[fx_data->NbVol])
         {
             DR_Error("Last composite vol date must be strictly"
                      " < first spot vol date \n");
             goto RETURN;
         }
     }


     for (i = 1; i <= fx_data->NbVol; i++)
          if ((fx_data->FxVol[i] < .000001) || (fx_data->FxVol[i] > 9.99))
          {
               DR_Error ("FX volatilities not entered correctly!");
               goto RETURN;
          }  /* if */

     if (fx_data->FxCutOffLevel < 0.0)
     {
         DR_Error("FX cut off level must be positive!\n");
         goto RETURN;
     }
     
     if ((fx_data->FxCutOffFlag != FALSE) &&
         (fx_data->FxCutOffFlag != TRUE))
     {
         DR_Error("FxCutOffFlag should be internally initialised to "
                  "TRUE of FALSE\n");
         goto RETURN;
     }

     if((fx_data->FxCutOffLast != TRUE) &&
         (fx_data->FxCutOffLast != FALSE))
     {
         DR_Error("fxCutOffLast should be internally initialised to "
                    " TRUE or FALSE\n");
         goto RETURN;
     }
     
        
     if( (fx_data->FxBootStrapMode != 0L) &&
         (fx_data->FxBootStrapMode != 1L)  &&
         (fx_data->FxBootStrapMode != 2L))
     {

         DR_Error("The internal initialisation of FxBootstrapMode is"
                    " invalid\n");
         goto RETURN; 
     }



        /* check FX smile parameters */
        if (fx_data->NbSmilePt < 1)
        {
            DR_Error("Hyb3_Fx_Input_W: Nb of FX smile param pts < 1!");
            goto RETURN;
        }
        
        for (i=0; i<fx_data->NbSmilePt; i++)
        {
            if (Dateok(fx_data->SmileDate[i]))
            {
                DR_Error ("Hyb3_Fx_Input_W: incorrect fx smile date format!");
                goto RETURN;
            }
        }
        
        for (i=1; i<fx_data->NbSmilePt; i++)
        {
            if (fx_data->SmileDate[i] <= fx_data->SmileDate[i-1])
            {
                DR_Error ("Hyb3_Fx_Input_W: "
                          "FX smile dates are not in ascending order!");
                goto RETURN;
            }
        }
        
        if (fx_data->SmileDate[0] <= fx_data->ValueDate)
        {
            DR_Error ("Hyb3_Fx_Input_W: some FX smile dates <= today!");
            goto RETURN;
        }


     
     status = SUCCESS;
          
     RETURN:
          
     return (status);

}  /* Hyb3_Fx_Check_W */


/*****  Hyb3_Correl_Input_WType3  *************************************************/
/**
 *  Read fx input for DR Wrappers of Type 3 and check validity of input. Note
 *  that Type3 Wrapper implies 6 overwrite strings.
 *
 */
int    Hyb3_Correl_Input_WType3
            (FX_DATA    *fx_data,                    /**< (O) Fx data          */
             char       OverWriteString[3][MAXBUFF], /**< (I) Overwrite strings*/
             char const *FileName)                   /**< (I) Filename incl ext*/
{
     double
          temp[3];

     int         
          readerror,                      /* Reading error status                 */
          status = FAILURE;             /* Error status = FAILURE initially */
     FILE
          *stream = NULL;


     stream = fopen (FileName, "r");             /* Open the correlation data file */

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Correl_Input_W3)", FileName);
          
          goto RETURN;
          
     }  /* if */                                                                
          
     if (FindAndSkipComLine (stream, "corr dom vs. foreign IR", 
                                     "Hyb3_Correl_Input_WType3", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[0]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W3)", FileName);
          
          goto RETURN;
     }


     if (FindAndSkipComLine (stream, "corr domestic IR vs. FX", 
                                     "Hyb3_Correl_Input_WType3", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[1]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W3)", FileName);
          
          goto RETURN;
     }

     if (FindAndSkipComLine (stream, "corr foreign IR vs. FX", 
                                     "Hyb3_Correl_Input_WType3", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[2]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W3)", FileName);
          
          goto RETURN;
     }
     /*********************************************/
     /* Although the fx_data structure supports a */
     /* term-structure of correlations, we only   */
     /* use the first element in Build Tree       */
     /*********************************************/
     fx_data->Rho[0][0] = temp[0];   /* Flat correlation curves */
     fx_data->Rho[1][0] = temp[2];
     fx_data->Rho[2][0] = temp[1];   /* Reverse here: Rho[1]=temp[2] and Rho[2]=temp[1] */
                


          

     if (strcmp (OverWriteString[0], "nil"))                      /* Correlation overwrite */
     {                    
          fx_data->Rho[0][0] = atof (OverWriteString[0]);/* same comment as above*/
                               
     }  /* if */
          
     if (strcmp (OverWriteString[1], "nil"))
     {
         fx_data->Rho[1][0] = atof (OverWriteString[1]);          
     }  /* if */
          
     if (strcmp (OverWriteString[2], "nil"))
     {                    
         fx_data->Rho[2][0] = atof (OverWriteString[2]);                               
     }  /* if */
          

     if (Hyb3_Correl_Check_WType3 (fx_data) == FAILURE)                  /* Check validity of input */
     {          
          goto RETURN;
     }
          

     status = SUCCESS;
          
     RETURN:

     if (stream != NULL) fclose (stream);
          
     return (status);

}  /* Hyb3_Correl_Input_WType3 */


/*****  Hyb3_Correl_Check_WType3  *************************************************/
/**
*	Check validity of correlation inputs.
*/
int    Hyb3_Correl_Check_WType3(FX_DATA     *fx_data)  /**< (I) Structure of fx data */
{
     int
          k,
          status = FAILURE;                     /* Error status = FAILURE initially */
                
    
    for (k = 0; k < 3; k++)

        if ((fx_data->Rho[k][0] < -MAXCORRELATION) || 
                     (fx_data->Rho[k][0] >  MAXCORRELATION))
        {
            DR_Error("Correlation outside allowable range (%6.4f < x < %6.4f)!",
                      -MAXCORRELATION, MAXCORRELATION);
            
            goto RETURN;
          
        }  /* if */                                                                


     status = SUCCESS;
          
     RETURN:
          
     return (status);

}  /* Hyb3_Correl_Check_WType3 */


/*****  Hyb3_Correl_Input_WType4  *************************************************/
/**
 *  Read eq input for DR Wrappers of Type 4 and check validity of input.
 *
 */
int    Hyb3_Correl_Input_WType4
            (EQ_DATA    *eq_data,                     /**< (O) Eq data           */
             char        OwriteCorrel[MAXBUFF],       /**< (I) Overwrite string  */
             char const *FileName)                    /**< (I) Filename incl ext */
{
     int         
          readerror,
          status = FAILURE;
     FILE
          *stream = NULL;


     stream = fopen (FileName, "r");

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Hyb3_Correl_Input_WType4)", FileName);
          
          goto RETURN;
          
     }                                                                
          
     if (FindAndSkipComLine (stream, "corr dom i.r. vs. equity", 
                                     "Hyb3_Correl_Input_WType3", FileName) == FAILURE)
     {          
          goto RETURN;
     }

     readerror = fscanf (stream, "%lf \n", &(eq_data->Rho[0][0]));

     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W3)", FileName);
          
          goto RETURN;
     }


     if (strcmp (OwriteCorrel, "nil")) 
     {                    
          eq_data->Rho[0][0] = atof(OwriteCorrel);
     }
          
     if (Hyb3_Correl_Check_WType4 (eq_data) == FAILURE)
     {          
          goto RETURN;
     }

     status = SUCCESS;
          
RETURN:
     if (stream != NULL) fclose(stream);
          
     return (status);

}  /* Hyb3_Correl_Input_WType4 */

/*****  Hyb3_Correl_Check_WType4  *************************************************/
/**
*	Check validity of correlation inputs.
*/
int    Hyb3_Correl_Check_WType4(EQ_DATA *eq_data)  /**< (I) Structure of eq data */
{
    int status = FAILURE;
    
    if ((eq_data->Rho[0][0] < -MAXCORRELATION) || 
        (eq_data->Rho[0][0] >  MAXCORRELATION))
    {
        DR_Error("Correlation outside allowable range (%6.4f < x < %6.4f)!",
                  -MAXCORRELATION, MAXCORRELATION);
        
        goto RETURN;
    } 

     status = SUCCESS;
          
RETURN:
          
     return (status);

}  /* Hyb3_Correl_Check_WType4 */

/*****  Hyb3_Correl_Input_WType6  *************************************************/
/**
 *  Read fx input for DR Wrappers of Type 6 and check validity of input. Note
 *  that Type6 Wrapper implies 9 overwrite strings. The first 3 are used  for
 *  the FX spot vol  cut-off facility and the  final 6 are correlation values
 *  used below.
 *
 *  Note that both  EQ_DATA  and  FX_DATA structures  contain  3  correlation 
 *  curves and the order is IrDom vs. IrFor, IrFor vs. FX and IrDom vs FX for
 *  FX _DATA and Eq vs. IrFor, Eq vs. IrDom and Eq vs. FX for EQ_DATA.
 *
 *  Dr Wrapper Type 6 implies a 3-factor  tree running with an equity as  the
 *  third asset.
 *
 *
 */
int    Hyb3_Correl_Input_WType6
           (FX_DATA   *fx_data,                     /**< (O) Fx data           */
            EQ_DATA   *eq_data,                     /**< (O) Equity data       */
            char       DomForFlag,                  /**< (I) Flag for eq denom */
            char       OverWriteString[6][MAXBUFF], /**< (I) Overwrite strings */
            char const *FileName)                   /**< (I) Filename incl ext */
{


     double
          temp[6];         /* Temporary vars holding corr values in corr.dat */
     int
          valueInFile[6];  /* Flag to indicate if a -999.0 value is found    */

     int
          i,
          readerror,                     /* Reading error status             */
          status = FAILURE;              /* Error status = FAILURE initially */

     FILE
          *stream = NULL;


     stream = fopen (FileName, "r"); /* Open the correlation data file   */

     if (stream == NULL)
     {
          DR_Error("Could not open file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
          
     }  /* if */                                                                
         
     /* IRFor vs. IRDom */ 
     if (FindAndSkipComLine (stream, "correl IRDom vs. IRFor", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[0]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }

     /* IRDom vs. FX */
     if (FindAndSkipComLine (stream, "correl IRDom vs. FX", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[1]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }

     /* IRFor vs. FX */
     if (FindAndSkipComLine (stream, "correl IRFor vs. FX", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[2]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }

     /* IRDom vs. Equity */
     if (FindAndSkipComLine (stream, "correl IRDom vs. Equity", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[3]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }

     /* IRFor vs. Equity */
     if (FindAndSkipComLine (stream, "correl IRFor vs. Equity", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[4]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }

     /* FX vs. Equity */
     if (FindAndSkipComLine (stream, "correl FX vs. Equity", 
                                     "Hyb3_Correl_Input_WType6", FileName) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(temp[5]));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("Could not read file %s! (Correl_Input_W6)", FileName);
          
          goto RETURN;
     }


     /* Kapital must place a -999.0 value if correl N/A */
     for (i=0; i<6; i++)
     {
          /* if temp[i] == CORR_NOT_AVAILABLE */
          if(ABS(temp[i] - CORR_NOT_AVAIL) < TINY)
          {
               valueInFile[i] = FALSE;
          }
          else
          {
               valueInFile[i] = TRUE;
          }
     }

     
          
     /* First correl overwrite */
     if (strcmp (OverWriteString[0], "nil"))                                                                                                                                                                                                                                                        	
    {
       
       fx_data->Rho[0][0] = atof (OverWriteString[0]);
                
    } 
    else
    {
        if (valueInFile[0] == FALSE)
        {
            DR_Error("Correl IRFor x IRDom not in corr.dat and not provided "
                    "as an overwrite!\n");
            goto RETURN;
        }
        
        fx_data->Rho[0][0] = temp[0];
    }
    

    /* Second correl overwrite */
    if (strcmp (OverWriteString[1], "nil"))	
    {
        
        fx_data->Rho[1][0] = atof (OverWriteString[1]);
                
    } 
    else
    {
        if (valueInFile[2] == FALSE)
        {
            DR_Error("Correl IRFor x FX not in corr.dat and not provided "
                    "as an overwrite!\n");
            goto RETURN;
        }
        
        fx_data->Rho[1][0] = temp[2];
    }

    /* Third correl overwrite */
    if (strcmp (OverWriteString[2], "nil"))	
    {
        
        fx_data->Rho[2][0] = atof (OverWriteString[2]);
                
    } 
    else
    {
        if (valueInFile[1] == FALSE)
        {
            DR_Error("Correl IRDom x FX not in corr.dat and not provided "
                    "as an overwrite!\n");
            goto RETURN;
        }
        
        fx_data->Rho[2][0] = temp[1];
    }

    /* 4th/5th o'writes: may not be necessary as corr may be in equity.dyn */
    /* Fourth correl overwrite (Eq x IRFor) */
    if (strcmp (OverWriteString[3], "nil"))	
    {
        for (i = 0; i <= eq_data->NbVol; i++)
        {
            eq_data->Rho[0][i] = atof (OverWriteString[3]);
        }	    
    } 
    else
    {  
        /* if IRFor vs Eq not in file and eq is D, then data missing */
        if (valueInFile[4]) 
        {
            for (i = 0; i <= eq_data->NbVol; i++)
            { 
                eq_data->Rho[0][i] = temp[4];
            }  
        }
        else
        {
            if (DomForFlag == 'D')
            { 
                DR_Error("Correl IRFor x Eq not in corr.dat and not provided "
                        "as an overwrite!\n");
                goto RETURN;
            }
        }
             
    }



    /* Fifth correl overwrite (Eq x IRDom) */
    if (strcmp (OverWriteString[4], "nil"))	
    {
                
        eq_data->Rho[1][0] = atof (OverWriteString[4]);        	    
    } 
    else
    {   
        /* if Eq x IRDom not in file and eq is F or C then data missing */
        if (valueInFile[3])
        {                        
            eq_data->Rho[1][0] = temp[3];            
        }
        else
        {
            if ((DomForFlag == 'F') || (DomForFlag == 'C'))
            {
                DR_Error("Correl IRDom x Eq not in corr.dat and not provided "
                        "as an overwrite!\n");
                goto RETURN;
            }
        }
    }




    /* Sixth correl overwrite (Eq x FX) */
    if (strcmp (OverWriteString[5], "nil"))	
    {                
        eq_data->Rho[2][0] = atof (OverWriteString[5]);        	    
    } 
    else
    {
        if (valueInFile[5])
        {                        
            eq_data->Rho[2][0] = temp[5];            
        }
        else
        {
            DR_Error("Correl FX x Eq not in corr.dat and not provided "
                    "as an overwrite!\n");
            goto RETURN;
        }            
    }


        

    if (Hyb3_Correl_Check_WType6(fx_data,
                            eq_data) == FAILURE)
    {        
        goto RETURN;
    }

            

    status = SUCCESS;
        
    RETURN:

    if (stream != NULL) fclose(stream);

    return (status);

}  /* Hyb3_Correl_Input_WType6 */




/*****  Hyb3_Correl_Check_WType6  *************************************************/
/**
*   Check validity of correlation inputs.
*/
int    Hyb3_Correl_Check_WType6(FX_DATA    *fx_data,  /**< (I) Structure of fx data */
                           EQ_DATA    *eq_data)
{
    int
        k,
        status = FAILURE;                /* Error status = FAILURE initially */
            
    
    for (k = 0; k < 3; k++)
        if ((fx_data->Rho[k][0] < -MAXCORRELATION) || 
            (fx_data->Rho[k][0] >  MAXCORRELATION))
        {
                DR_Error("Correlation outside allowable range (%6.4f < x < %6.4f)!",
                         -MAXCORRELATION, MAXCORRELATION);
                
                goto RETURN;
        
        }  /* if */     
            
            
    
    for (k = 0; k < 3; k++)
        if ((eq_data->Rho[k][0] < -MAXCORRELATION) || 
                (eq_data->Rho[k][0] >  MAXCORRELATION))
        {
                DR_Error("Correlation outside allowable range (%6.4f < x < %6.4f)!",
                         -MAXCORRELATION, MAXCORRELATION);
                
                goto RETURN;
        
        }  /* if */                                               


    /* One further check for Type 6 Dr Wrapper: correlation matrix must be */
    /* positive definite.  We know  it is symmetric, so  we check that the */
    /* determinant is above a certain minimum.                             */
    {
        double det;
        double r1, r2, r3;

        r1 = fx_data->Rho[0][0]; /* IR vs IR  */
        r2 = eq_data->Rho[0][0]; /* IRf vs Eq */
        r3 = eq_data->Rho[1][0]; /* IRd vs Eq */


        det = 1.0 - SQUARE(r1) - SQUARE(r2) - SQUARE(r3) + 2.0*r1*r2*r3;

        if (det < (1.0 - SQUARE(MAXCORRELATION)))
        {
            DR_Error("Correlation matrix ill-defined!\n");
            goto RETURN;
        }
    }


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Hyb3_Correl_Check_WType6 */


/* read overwrites for the various model parameters */

int   Hyb3_ModelOverwrites_Input_W(
         FILE* stream,
         char* FileName,
         char OwriteFxSpot[MAXBUFF],
         char NbFXSmileOWS[MAXBUFF],
         char FXSmileParamOWS[MAXNBDATE][MAXBUFF],
         char OwriteCorrel[3][MAXBUFF],    /**< O'writes related to correlation.dat */
         char OwriteMParamF[5][MAXBUFF],   /**< O'writes related to modelParams.dat */
         char OwriteMParamD[5][MAXBUFF],   /**< O'writes related to modelParams.dat */
         HYB3_TREE_DATA *tree_data)
{

    int i;
    int readerror;
    int status = FAILURE;
    char ErrorMsg[MAXBUFF];
    char* getserror;
    
    char   momMatchFlag;             /* (Y)es or (N)o: for moment matching method */

    /* FX volatility calibration choice  */
    if (FindAndSkipComLine (stream, "FX vol calibration overwrite",         
                        "Hyb3_ModelOverwrites_Input_W", FileName) == FAILURE)
    {                                                            
        goto RETURN;                                             
    }                                                            
    readerror = fscanf (stream, "%s \n",OwriteFxSpot);        
    if (readerror != 1)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read FX vol calibration overwrite in " 
                "file %s! (Hyb3_ModelOverwrites_Input_W)", FileName);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }                                                            

    /* FX smile parameters */
    if (FindAndSkipComLine (stream,
                            "Nb of FX smile parameter sets OWS",
                            "Hyb3_ModelOverwrites_Input_W",
                            FileName) == FAILURE) goto RETURN;

    readerror = fscanf (stream, "%s \n", NbFXSmileOWS);
    if (readerror != 1)
    {
        DR_Error("Unable to read nb FX smile params OWS! " 
                 "(Hyb3_ModelOverwrites_Input_W)");
        goto RETURN;
    }

    if (FindAndSkipComLine (stream,
                            "FX smile parameters OWS",
                            "Hyb3_ModelOverwrites_Input_W",
                            FileName) == FAILURE) goto RETURN;

    if (strstr(NbFXSmileOWS, "nil") == NULL)
    {
        long NbFXSmile = atoi(NbFXSmileOWS);

        if ((NbFXSmile < 0) ||
            (NbFXSmile > MAXNBDATE))
        {
            DR_Error("Nb FX smile params OWS out of range! " 
                     "(Hyb3_ModelOverwrites_Input_W)");
            goto RETURN;
        }

        for (i=0; i<NbFXSmile; i++)
        {
            getserror = fgets (FXSmileParamOWS[i], MAXBUFF, stream);
            if (getserror == NULL)
            {
                DR_Error("Can't read FX Smile Params!" 
                         "(Hyb3_ModelOverwrites_Input_W)");
                goto RETURN;
            }
        }
    }
    else
    {
        /* FXsmileParamOWS[0] has to be set to "nil" for internal usage.  For the wrapper input file,
           if the number of fxsmile parameters = "nil", then the overwrite parameters should be
           a blank line - as if NbFXSmileOWS = 0 */
        strcpy(FXSmileParamOWS[0], "nil");
    }

    /* Correlation overwrites */
    for (i=0; i<3; i++)
    {
        if (FindAndSkipComLine (stream, "corr overwrites", 
                            "Hyb3_ModelOverwrites_Input_W", FileName) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, "%s \n",OwriteCorrel[i]);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, "Could not read corr overwrites in "
                    "file %s! (ModelOverwrite_Input_W)", FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }
    
    
    /* Foreign interest rate model specification (normal, lognormal etc) */
    if (FindAndSkipComLine (stream, "Foreign IR model type", 
                            "ModelOverwrite_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    getserror = fgets (OwriteMParamF[1], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, "Could not read foreign model type in file %s! "
                        "(ModelOverwrite_Input_W)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    

    /* Domestic interest rate model specification (normal, lognormal etc) */
    if (FindAndSkipComLine (stream, "Domestic IR model", 
                            "ModelOverwrite_Input_W", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    getserror = fgets (OwriteMParamD[1], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, "Could not read domestic model type in file %s! "
                        "(ModelOverwrite_Input_W)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
   
    /* Foreign mean reversion */
    if (FindAndSkipComLine (stream, "m.r. overwrite",         
                        "ModelOverwrite_Input_W", FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    readerror = fscanf (stream, "%s \n",OwriteMParamF[3]);        
    if (readerror != 1)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read m.r. overwrite in " 
                "file %s! (ModelOverwrite_Input_W)", FileName);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }    
     
    /* Domestic mean reversion */
    if (FindAndSkipComLine (stream, "m.r. overwrite",         
                        "ModelOverwrite_Input_W", FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    readerror = fscanf (stream, "%s \n",OwriteMParamD[3]);        
    if (readerror != 1)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read m.r. overwrite in " 
                "file %s! (ModelOverwrite_Input_W)", FileName);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }

    /* check if Moment Matching on/off (last argument thus we can treat it as */
    /* purely optional) : default: switch on                                  */
    tree_data->FxMomentMatching = TRUE;
    if (Hyb3_FindAndSkipComLineOptional (stream, "moment matching flag",
                        "Fxseries_Manager", FileName) == FAILURE)
    {
      /*  printf("Using default value FALSE for moment matching flag!\n");*/
    }
    else
    {
        readerror = fscanf (stream, "%c\n", &momMatchFlag );
        if (readerror != 1)
        {
            sprintf(ErrorMsg, "Could not read moment matching flag "
                "file %s! (Fxseries_Manager)", FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }  
        tree_data->FxMomentMatching = ( (char)toupper(momMatchFlag) == 'Y' ); 
    }

    

    /* Initialise the remaining strings with the default values */
    sprintf(OwriteMParamF[0],"%3d",tree_data->Ppy);
    sprintf(OwriteMParamD[0],"%3d",tree_data->Ppy);
    sprintf(OwriteMParamF[2],"nil");  /* Weights  */
    sprintf(OwriteMParamD[2],"nil");
    sprintf(OwriteMParamF[4],"nil");  /* Factor correlations */
    sprintf(OwriteMParamD[4],"nil");
 
    status = SUCCESS;

    RETURN:

    return status;

}

/* read overwrites for the various model parameters */
/* in the Hyb2+1 CUPS mode with multi factor IR     */
int   Hyb3_ModelOverwritesMultFactCups_Input_W(
         FILE* stream,
         char* FileName,
         char OwriteFxSpot[MAXBUFF],
         char NbFXSmileOWS[MAXBUFF],
         char FXSmileParamOWS[MAXNBDATE][MAXBUFF],
         char OwriteCorrel[3][MAXBUFF],    /**< O'writes related to correlation.dat */
         char OwriteMParamF[5][MAXBUFF],   /**< O'writes related to modelParams.dat */
         char OwriteMParamD[5][MAXBUFF],   /**< O'writes related to modelParams.dat */
         HYB3_TREE_DATA *tree_data)
{

    int i;
    int readerror;
    int status = FAILURE;
    char ErrorMsg[MAXBUFF];
    char* getserror;
    static char routine[] = "ModelOverwritesMultiFactorCups_Input_W";


    /* FX volatility calibration choice  */
    if (FindAndSkipComLine (stream, "FX vol calibration overwrite",         
                        routine, FileName) == FAILURE)
    {                                                            
        goto RETURN;                                             
    }                                                            
    readerror = fscanf (stream, "%s \n",OwriteFxSpot);        
    if (readerror != 1)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read FX vol calibration overwrite in " 
                "file %s! (%s)",FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }                                                            



    /* Correlation overwrites */
    for (i=0; i<3; i++)
    {
        if (FindAndSkipComLine (stream, "corr overwrites", 
                            routine, FileName) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, "%s \n",OwriteCorrel[i]);
        if (readerror != 1)
        {        
            sprintf(ErrorMsg, "Could not read corr overwrites in "
                    "file %s! (%s)", FileName, routine);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }
    
    
    /* Foreign interest rate model specification (normal, lognormal etc) */
    if (FindAndSkipComLine (stream, "Foreign IR model type", 
                            routine, FileName) == FAILURE)
    {        
        goto RETURN;
    }

    getserror = fgets (OwriteMParamF[1], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, "Could not read foreign model type in file %s! "
                        "(%s)", FileName, routine);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    

    /* Domestic interest rate model specification (normal, lognormal etc) */
    if (FindAndSkipComLine (stream, "Domestic IR model", 
                            routine, FileName) == FAILURE)
    {        
        goto RETURN;
    }

    getserror = fgets (OwriteMParamD[1], MAXBUFF, stream);
    if (getserror  == NULL)
    {
        sprintf (ErrorMsg, "Could not read domestic model type in file %s! "
                        "(%s)", FileName, routine);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
   
    /*Foreign Factor weights */
    if (FindAndSkipComLine (stream, "factor weight overwrite",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamF[2], MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read foreign factor weight overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }

    /*Domestic Factor weights */
    if (FindAndSkipComLine (stream, "factor weight overwrite",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamD[2], MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read domestic factor weight overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }

    /* Foreign mean reversion */
    if (FindAndSkipComLine (stream, "m.r. overwrite",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamF[3], MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read m.r. overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }    
     
    /* Domestic mean reversion */
    if (FindAndSkipComLine (stream, "m.r. overwrite",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamD[3],MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read m.r. overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }   

    /* Foreign factor correlation */
    if (FindAndSkipComLine (stream, "factor correlation",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamF[4], MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read factor correlation overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }    

    /* Domestic factor correlation */
    if (FindAndSkipComLine (stream, "factor correlation",         
                        routine, FileName) == FAILURE)    
    {                                                            
        goto RETURN;                                             
    }                                                            
    getserror = fgets (OwriteMParamD[4], MAXBUFF, stream);        
    if (getserror == NULL)                                          
    {                                                            
        sprintf(ErrorMsg, "Could not read factor correlation overwrite in " 
                "file %s! (%s)", FileName, routine);            
        DR_Error (ErrorMsg);                                     
        goto RETURN;                                             
    }    


    

    /* Initialise the remaining strings with the default values */
    sprintf(OwriteMParamF[0],"%3d",tree_data->Ppy);
    sprintf(OwriteMParamD[0],"%3d",tree_data->Ppy);
    strcpy(NbFXSmileOWS, "nil");
    strcpy(FXSmileParamOWS[0], "nil");



 
    status = SUCCESS;

    RETURN:

    return status;

}
