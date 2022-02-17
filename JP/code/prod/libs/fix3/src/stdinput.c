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
#include <irx/irxutilio.h>

/*****  Fix3_MktVolAndModel_Input  **********************************************/
/*
*  	Read volatility and model parameters  and zero curves
*/
int Fix3_MktAndModel_Input (
        MKTVOL_DATA    *mktvol_data,                  /* (O) Volatility data     */
        FIX3_TREE_DATA *tree_data,                    /* (O) Tree data           */
        T_CURVE        *t_curve,                      /* (I) Zero curve          */
        char           Index[2][8],                   /* (I) Calibration indices */ 
        char           OverWriteString[6][MAXBUFF] )  /* (I) Overwrite strings   */
{
    int status = FAILURE;

     /* read zero curve info */
    if (Fix3_ZeroCurve_Input (t_curve,
                              mktvol_data,
                              tree_data) == FAILURE)
    {
        goto RETURN;
    }

    /* read volatility and model data */
    if (Fix3_MktVolAndModel_Input (mktvol_data,
                                   tree_data,
                                   t_curve,
                                   Index,
                                   OverWriteString) == FAILURE)
    {
        goto RETURN;
    }
    
    status = SUCCESS;
    
RETURN:

    return(status);
} /* Fix3_MktAndModel_Input*/




/*****  Fix3_ZeroCurve_Input **********************************************/
/*
*  	Read zero curves
*/
int Fix3_ZeroCurve_Input(
        T_CURVE        *t_curve,         /* (O) Zero curve            */
        MKTVOL_DATA    *mktvol_data,     /* (I/0) Volatility data     */
        FIX3_TREE_DATA *tree_data)       /* (I/0) Tree data           */
       
{
    FILE *streamRisk = NULL;
    int status = FAILURE;
    char* modelChoiceString = NULL;

    if (ModelInfo("model.dat", "FIX3", &modelChoiceString) != SUCCESS)
        goto RETURN;
  
    /* Identify model string */
    if (!strcmp(modelChoiceString,"ORIGINAL"))
        mktvol_data->ModelChoice = FIX3_ORIGINAL;    
    else if (!strcmp(modelChoiceString,"CLASSIC"))
        mktvol_data->ModelChoice = FIX3_CLASSIC;
    else if (!strcmp(modelChoiceString,"TIMEDEP"))
        mktvol_data->ModelChoice = FIX3_TIMEDEP;
    else if (!strcmp(modelChoiceString,"2Q"))
        mktvol_data->ModelChoice = FIX3_TIMEDEP;
    else if (!strcmp(modelChoiceString,"SMD"))
        mktvol_data->ModelChoice = FIX3_SMD;
    else if (!strcmp(modelChoiceString,"TMX"))
        mktvol_data->ModelChoice = FIX3_TMX;
    else if (!strcmp(modelChoiceString,"E2Q"))
        mktvol_data->ModelChoice = FIX3_E2Q;
    else 
    {
        DR_Error("Unrecognized model choice: %s.\n", modelChoiceString);
        goto RETURN;
    }

    /* Is this a numeraire model? */
    mktvol_data->IsNmrModel = (mktvol_data->ModelChoice == FIX3_TMX);

    /* Initialize model interface */
    if (Fix3_Model_Interface_Init(mktvol_data->ModelChoice) == FAILURE)
    {
        goto RETURN;
    }
    
    switch (mktvol_data->ModelChoice)
    {
    case (FIX3_ORIGINAL):

        /* Read index curve from zero.dat */ 
        if (Term_Input_W (&(t_curve[0]),
                          "zero.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from disczero.dat */
        if (Term_Input_W (&(t_curve[1]),
                          "disczero.dat") == FAILURE)
        {
            goto RETURN;
        }

      
        streamRisk = fopen ("riskzero.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_W ( &(t_curve[2]),
                              "disczero.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_W ( &(t_curve[2]),
                              "riskzero.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        break;

    case (FIX3_CLASSIC):

        /* read info frum "summary.dat" */
        if (SummaryInfo(t_curve,        
                        "summary.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* read ir info */
        if( MktInfo(mktvol_data, 
                    &tree_data->CvDiff,
                    &tree_data->CvDisc,
                    t_curve,         
                    "ir_info_0.dat")  == FAILURE)
        {        
            goto RETURN;
        }

        /* Read index curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[0]),
                          "ir_curve0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[1]),
                          "ir_curve1_0.dat") == FAILURE)
        {
            goto RETURN;
        }

         /* Read riskzero curve, if failed use discount */
        streamRisk = fopen ("ir_curve2_0.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve1_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve2_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }

        break;

    case (FIX3_TIMEDEP):
        /* read info from summary.dat */
        if (SummaryInfo(t_curve,        
                        "summary.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* read ir info */
        if( MktInfo(mktvol_data, 
                    &tree_data->CvDiff,
                    &tree_data->CvDisc,
                    t_curve,         
                    "ir_info_0.dat")  == FAILURE)
        {        
            goto RETURN;
        }

        /* Read index curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[0]),
                          "ir_curve0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from MAW style file */
        if (Term_Input_New_W (&(t_curve[1]),
                          "ir_curve1_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read riskzero curve, if failed use discount */
        streamRisk = fopen ("ir_curve2_0.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve1_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve2_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }

        break;

    case (FIX3_TMX):
        
        /* read info from summary.dat */
        if (SummaryInfo(t_curve,        
                        "summary.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* read ir info */
        if( MktInfo(mktvol_data, 
                    &tree_data->CvDiff,
                    &tree_data->CvDisc,
                    t_curve,         
                    "ir_info_0.dat")  == FAILURE)
        {        
            goto RETURN;
        }

        /* Read index curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[0]),
                          "ir_curve0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from MAW style file */
        if (Term_Input_New_W (&(t_curve[1]),
                          "ir_curve1_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read riskzero curve, if failed use discount */
        streamRisk = fopen ("ir_curve2_0.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve1_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve2_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }

        break;
  
    case(FIX3_SMD):

         /* read info from summary.dat */
        if (SummaryInfo(t_curve,        
                        "summary.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* read ir info */
        if( MktInfo(mktvol_data, 
                    &tree_data->CvDiff,
                    &tree_data->CvDisc,
                    t_curve,         
                    "ir_info_0.dat")  == FAILURE)
        {        
            goto RETURN;
        }

        /* Read index curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[0]),
                          "ir_curve0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from MAW style file */
        if (Term_Input_New_W (&(t_curve[1]),
                          "ir_curve1_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read riskzero curve, if failed use discount */
        streamRisk = fopen ("ir_curve2_0.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve1_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve2_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
  
        break;


     case(FIX3_E2Q):

         /* read info frum "summary.dat" */
        if (SummaryInfo(t_curve,        
                        "summary.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* read ir info */
        if( MktInfo(mktvol_data, 
                    &tree_data->CvDiff,
                    &tree_data->CvDisc,
                    t_curve,         
                    "ir_info_0.dat")  == FAILURE)
        {        
            goto RETURN;
        }

        /* Read index curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[0]),
                          "ir_curve0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read discount curve from MAW style file*/
        if (Term_Input_New_W (&(t_curve[1]),
                          "ir_curve1_0.dat") == FAILURE)
        {
            goto RETURN;
        }

         /* Read riskzero curve, if failed use discount */
        streamRisk = fopen ("ir_curve2_0.dat", "r");
        if (streamRisk == NULL)
        {
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve1_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }
        else
        {
            fclose(streamRisk);
            if (Term_Input_New_W (&(t_curve[2]),
                              "ir_curve2_0.dat") == FAILURE)
            {
                goto RETURN;
            }
        }

        
        break;


    default:
        DR_Error ("Unsupported model choice %s \n", modelChoiceString);
        goto RETURN;
   
    }
       
     
     status = SUCCESS;

    RETURN:

    if (modelChoiceString != NULL)
        free(modelChoiceString);

    return(status);

} /* Fix3_ZeroCurve_Input */


/*****  Fix3_MktVolAndModel_Input  **********************************************/
/*
*  	Read volatility and model parameters and check validity of input.
*
*  NOTE:  Calibration indices and overwrite strings are only relevant when operating
*         in legacy (FIX3_ORIGINAL) mode.
*/
int Fix3_MktVolAndModel_Input (
        MKTVOL_DATA    *mktvol_data,                  /* (O) Volatility data     */
        FIX3_TREE_DATA *tree_data,                    /* (O) Tree data           */
        T_CURVE        *t_curve,                      /* (I) Zero curve          */
        char           Index[2][8],                   /* (I) Calibration indices */ 
        char           OverWriteString[6][MAXBUFF] )  /* (I) Overwrite strings   */
{

    int     status = FAILURE;   /* Error status = FAILURE initially */
   

    switch (mktvol_data->ModelChoice)
    {

    case (FIX3_ORIGINAL):
        if (strcmp (Index[0], Index[1]))
        {
            DR_Error ("Volatility calibration indices should be identical");
            goto RETURN;
        }

        /* Read volatility curve */
        if (MktVol_Input_W (mktvol_data,
                            Index[0],
                            &(t_curve[tree_data->CvDiff]),
                            "basevol.dat",
                            "swapvol.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read Vnfm model parameters */
        if (Fix3_Param_Input (mktvol_data,
                              tree_data,
                              tree_data->NbFactor,
                              OverWriteString,
                              "modelParameters.dat") == FAILURE)
        {
            goto RETURN;
        }

        break;

    case (FIX3_CLASSIC):

        if (MktVol_Input_New_W (mktvol_data,
                            t_curve,
                            "ir_voldiag0_0.dat") == FAILURE)
        {
            goto RETURN;
        }


        /* Read Vnfm model parameters */
        if (Fix3_Param_Input_Classic (mktvol_data,
                                      tree_data,
                                      "numerics.dat",
                                      "ir_mpar_0.dat") == FAILURE)
        {
            goto RETURN;
        }

       
        break;

    case (FIX3_TIMEDEP):

        if (MktVol_Input_New_W (mktvol_data,
                            t_curve,
                            "ir_voldiag0_0.dat") == FAILURE)
        {
            goto RETURN;
        }


        /* Read Vnfm model parameters */
        if (Fix3_Param_Input_TimeDep (mktvol_data,
                                      tree_data,
                                      "numerics.dat",
                                      "ir_mpar_0.dat") == FAILURE)
        {
            goto RETURN;
        }

      
        break;

    case (FIX3_SMD):


        if (MktVol_Input_New_W (mktvol_data,
                            t_curve,
                            "ir_voldiag0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Read Vnfm model parameters */
        if (Fix3_Param_Input_Smd (mktvol_data,
                                  tree_data,
                                  "numerics.dat",
                                  "ir_mpar_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        break;

     case (FIX3_E2Q):

        if (MktVol_Input_New_W (mktvol_data,
                            t_curve,
                            "ir_voldiag0_0.dat") == FAILURE)
        {
            goto RETURN;
        }


        /* Read Vnfm model parameters */
        if (Fix3_Param_Input_E2Q (mktvol_data,
                                  tree_data,
                                  "numerics.dat",
                                  "ir_mpar_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        break;

    case (FIX3_TMX):
        if (MktVol_Input_New_W (mktvol_data,
                            t_curve,
                            "ir_voldiag0_0.dat") == FAILURE)
        {
            goto RETURN;
        }

        /* Only lognormal vols allowed -- temporary restriction */
        if (mktvol_data->VolUnit == 0)
        {
            DR_Error ("Only lognormal vols allowed in TMX model");
            goto RETURN;
        }

        /* Read Vnfm model parameters */
        if (Fix3_Param_Input_Tmx (mktvol_data,
                                  tree_data,
                                  "numerics.dat",
                                  "ir_mpar_0.dat") == FAILURE)
        {
            goto RETURN;
        }


        break;
    }
   
    
    status = SUCCESS;

    RETURN:

    return (status);
} /* Fix3MktVolAndModel_Input */




/*****  Fix3_Param_Input  ********************************************************/
/*
*  	Read model parameters and check validity of input.
*/
int     Fix3_Param_Input (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             FIX3_TREE_DATA     *tree_data,             /* (O) Tree data structure           */
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

    /* make a copy of tree_data->NbFactor in mktvol_data */
    mktvol_data->NbFactor = NbFactor;

    /* Initialize model interface -- only for backward compatibility;
     * No longer needed when model I/O is done through Fix3_MktVolAndModel_Input */
    if (Fix3_Model_Interface_Init(mktvol_data->ModelChoice) == FAILURE)
    {
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

    
     mktvol_data->NbFactor = NbFactor;
        
    /* Check validity of input */
    if (Fix3_Param_Check (NbFactor,
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


/*****  Fix3_Param_Check_Original***********************************************/
/*
*  	Read term structure input and check validity of input.
*/
int    Fix3_Param_Check_Original (  
                int              NbFactor,     /* (I) Number of factors                     */
                MKTVOL_DATA      *mktvol_data, /* (I) Structure of swaption volatility data */
                FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure                   */
{
    int
        status = FAILURE;                      /* Error status = FAILURE initially */
       
    char
        ErrorMsg[MAXBUFF];

    double  norm;
    int     i;

    mktvol_data->NbFactor = tree_data->NbFactor;

    /* Legacy DLL and XLL support -- set Model Interface here */
    if (Fix3_Model_Interface_Init (mktvol_data->ModelChoice) == FAILURE)
    {
        goto RETURN;
    }

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

        if ((mktvol_data->Beta[0] < -10.) || (mktvol_data->Beta[0] > 10.))
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

        if ((mktvol_data->Beta[0] < -10.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[1] < -10.) || (mktvol_data->Beta[1] > 10.))
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

        if ((mktvol_data->Beta[0] < -10.) || (mktvol_data->Beta[0] > 10.))
        {
            DR_Error("Beta #1 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[1] < -10.) || (mktvol_data->Beta[1] > 10.))
        {
            DR_Error("Beta #2 out of range!");
            goto RETURN;
        }

        if ((mktvol_data->Beta[2] < -10.) || (mktvol_data->Beta[2] > 10.))
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

}  /* Fix3_Param_Check*/


/*****  Fix3_Param_Input_Classic  ********************************************************/
/*
*  	Read model parameters and check validity of input.
*/
int     Fix3_Param_Input_Classic (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             FIX3_TREE_DATA*tree_data,                  /* (O) Tree data structure           */
             char const*   FileNameTreeNum,             /* (I) File name tree numerical pars */
             char const*   FileNameModelPar)            /* (I) File name model parameters    */
{
    int     status = FAILURE;   /* Error status = FAILURE initially */

    /* Read model parameters file */
    if (Fix3_Model_Input_Classic (mktvol_data, 
                                  tree_data,
                                  FileNameModelPar) == FAILURE)
    {
        goto RETURN;
    }

    /* Read numerical parameters file */
    if (Fix3_Num_Input_New(mktvol_data, 
                           tree_data,
                           FileNameTreeNum) == FAILURE)
    {
        goto RETURN;
    }

    /* check validity of input */
    if (Fix3_Param_Check_Classic (mktvol_data->NbFactor, mktvol_data, tree_data) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;
        
  RETURN:
        
    return status;
}


/*****  Fix3_Param_Check_Classic  ********************************************************/
/*
*  	Read time dependent term structure input and check validity of input.
*/
int     Fix3_Param_Check_Classic(int NbFactor,                  /* (I) NbFactor for interface consistency    */
                                 MKTVOL_DATA *mktvol_data,      /* (I) Structure of swaption volatility data */
                                 FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure                   */
{
    int status = FAILURE;                      /* Error status = FAILURE initially */ 
    int     k;     
    int     NbTDInp;
    int     NbCorr;

    /* Nb of factors */
    if (mktvol_data->NbFactor != NbFactor || tree_data->NbFactor != NbFactor)
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

    
    if (IS_EQUAL(mktvol_data->Bbq,1))
    {
        mktvol_data->VolNorm = 0.;
        mktvol_data->VolLogn = 1.;
    }
    else
    if (IS_EQUAL(mktvol_data->Bbq,0))
    {
        mktvol_data->VolNorm = 1.;
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

    
    /* Nb of iterations */
    if (mktvol_data->CetNbIter > MAX_ITERATIONS)
    {
        DR_Error("Maximum allowed number of iterations is %d!", MAX_ITERATIONS);
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


    NbCorr = NbFactor * (NbFactor - 1) / 2;

    NbTDInp = mktvol_data->NbTDInp;

    /* check number of TD inp*/
    if (mktvol_data->NbTDInp != 1)
    {
        DR_Error ("Number of time dependent inputs should be  1 - number supplied "
                  " = %d", mktvol_data->NbTDInp);
        goto RETURN;
    }
  
    /* check TD input date format */
    if (Dateok(mktvol_data->TDInpDate[0]))
    {
        DR_Error("Incorrect format for time dependent input date [%ld] !", 
            mktvol_data->TDInpDate[0]);
        goto RETURN;
    }
    

    for (k = 0; k < NbFactor; k++)
	{
        /* check factor weights */
        if (mktvol_data->Alpha[k] < 0.0001)
        {
            DR_Error("Weight nb %d out of range!",  k );
            goto RETURN;
        }

        if ((mktvol_data->Beta[k] < -10.) || (mktvol_data->Beta[k] > 10.))
        {
            DR_Error("Beta nb %d  out of range!", k);
            goto RETURN;
        }
        
    }
    /* check correlation*/
    for (k = 0; k < NbCorr; k++)
	{
        if ((mktvol_data->Rho[k]< -.95) || (mktvol_data->Rho[k] > .95))
        {
            DR_Error("Correlation nb %d  out of range!", k );
            goto RETURN;
        }
    }


    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Fix3_Param_Check_Classic */



/*****Fix3_Model_Input_Classic  *********************************************************/
/*
* Read time-dependent mean-reversion, correlation, weight input      
*/
int Fix3_Model_Input_Classic (MKTVOL_DATA    *mktvol_data, /* (O) Volatility data                   */
                              FIX3_TREE_DATA *tree_data,   /* (0) For NbFactors & curve assignments */
                              char const     *FileName)    /* (I) File name including extension     */           
{
    int     NbTDInp, NbFactor, SmileNbEntries;
    long    SmileDate;
    int     k;
  
    int     NbCorr; /* nb correlation needed */
    char    ErrorMsg[MAXBUFF];
    int     readerror;
    FILE     *stream = NULL;
    int status = FAILURE;
    

    stream = fopen (FileName, "r");
    if (stream == NULL)
    {
        DR_Error("Could not open model parameters file: %s.\n", FileName);
        goto RETURN;
    }

    /* 
     *  Parameter file exists: and some of the overwrites are nil
     *  use the values in the parameter file.
     */
    /*  Title */
	if (FindAndSkipComLine (stream, 
                            "Title line",
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Section:VNFM parameters", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
   

     /* Number of factors */
    if (FindAndSkipComLine (stream, 
                            "Number of factors", 
                            "Fix3_Model_Input_Classic", 
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
    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }
    tree_data->NbFactor = mktvol_data->NbFactor = NbFactor;

    /*compute number of correlations */
    NbCorr = NbFactor * (NbFactor - 1) / 2;

    /* Number of time dependent input  dates */
    if (FindAndSkipComLine (stream, 
                            "Number of benchmark dates", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbTDInp);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of time dependent inputs dates in file %s! "
                "(Fix3_Model_Input_Classic)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* if classic version of the tree one and only one time dep input date */
    if (NbTDInp != 1)
    {
        sprintf(ErrorMsg, 
                "Only one time dependent input date should be specified for classic version in file %s! "
                "(Fix3_Model_Input_Classic)", 
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
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld ", &(mktvol_data->TDInpDate[0]));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
            "Could not read td input date in file %s! "
            "(Fix3_Model_Input_Classic)", 
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
                "(Fix3_Model_Input_Classic)", 
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
                "(Fix3_Model_Input_Classic)", 
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
                "(Fix3_Model_Input_Classic)", 
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
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Backbone section */
    if (FindAndSkipSectionLine (
                             1, /* skip mode */
                             stream, 
                            "Section:Backbone", 
                            "Fix3_Model_Input_Classic", 
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
        sprintf (ErrorMsg, "Could not find Bbq in file %s! (Fix3_Model_Input_Classic)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    mktvol_data->Bbq = 1. - mktvol_data->Bbq;

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    /* SMILE SECTION */
    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Smile", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Smile number of entries */
    if (FindAndSkipComLine (stream, 
                            "Number of entries", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream,
                        "%d \n", 
                         &SmileNbEntries);

    if (readerror != 1)
    {      
        sprintf (ErrorMsg, "Could not find Smile Nb entries in file %s! (Fix3_Model_Input_Classic)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if(SmileNbEntries != 1)
    {
        sprintf (ErrorMsg, "Number of entries in Smile section has to be 1 in file %s! (Fix3_Model_Input_Classic)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Read smile parameters */
    if (FindAndSkipComLine (stream, 
                            "Date, QLeft, QRight, FwdShift", 
                            "Fix3_Model_Input_Classic", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%ld \t%lf \t%lf \t%lf \n", 
                        &SmileDate, 
                        &(mktvol_data->QLeft),
                        &(mktvol_data->QRight),
                        &(mktvol_data->FwdShift));

    if (readerror != 4)
    {        
        sprintf (ErrorMsg, "Could not read smile info in file %s! (Fix3_Model_Input_Classic)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    SmileDate = IRDateFromYMDDate(SmileDate);
    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
    mktvol_data->QRight = 1. - mktvol_data->QRight;

    if (Dateok(SmileDate))
    {
        sprintf(ErrorMsg, "Incorrect format for smile date in file %s (Fix3_Model_Input_TimeDep)!", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_Classic", 
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

} /* Fix3_Model_Input_Classic */


/*****  Fix3_Num_Input_Common  *********************************************************/
/*
* Read numerical parameters: ppy, nb stdev, nb cet iterations from numerics.dat     
*/
static int Fix3_Num_Input_Common (MKTVOL_DATA     *mktvol_data,   /* (O) Volatility data                 */
                                  FIX3_TREE_DATA  *tree_data,     /* (O) Tree data                       */
                                  char const*     FileName,       /* (I) File name including extension   */ 
                                  int fileMustExist,              /* (I) The file must exist?            */
                                  int populateTreeNbSigmaMax)     /* (I) Populate tree_data->NbSigmaMax? */
{
    int status = FAILURE;           /* Error status = FAILURE initially */
    int     readerror;              /* Reading error status             */
    FILE    *stream = NULL;
    int     NbSigmaMax;

    stream = fopen (FileName, "r");
    if (stream == NULL)
    {
        if (fileMustExist)
        {
            DR_Error("Could not open tree numerics file: %s.\n", FileName);
            goto RETURN;
        }
        else
        {
            /* Legacy mode:  If file does not exist, we assume that params in OverWriteString
               (processed by caller). */ 
            return SUCCESS;
        }
    }
    
    /*number of Ppy's*/
    if (FindAndSkipComLine (stream, "Number of periods per year", "Fix3_Num_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                        &(tree_data->Ppy));

    if (readerror != 1)
    {        
        DR_Error("Could not find Ppy in file %s! (Fix3_Num_Input)", FileName);
        goto RETURN;
    }

    /* number of st dev; we will use the one from the deal file so read this one into a dummy variable*/
    if (FindAndSkipComLine (stream, "Number of standard deviations", "Fix3_Num_Input",FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                        &NbSigmaMax);

    if (readerror != 1)
    {        
        DR_Error("Could not find number of std in file %s! (Fix3_Num_Input)", FileName);
        goto RETURN;
    }

    if (populateTreeNbSigmaMax)
        tree_data->NbSigmaMax = NbSigmaMax;

    /* number of CET iterations */
    if (FindAndSkipComLine (stream, "Number of cet iterations", "Fix3_Num_Input",FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                        &(mktvol_data->CetNbIter));

    if (readerror != 1)
    {        
        DR_Error("Could not find number of cet iterations in file %s! (Fix3_Num_Input)", FileName);
        goto RETURN;
    }

    /* number of st dev for state variables*/
    if (FindAndSkipComLine (stream, "Number of standard deviations for state variables", 
                            "Fix3_Num_Input",FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                        &(tree_data->NbStdDevStates));

    if (readerror != 1)
    {        
        DR_Error("Could not find number of std for state var in file %s! (Fix3_Num_Input)", FileName);
        goto RETURN;
    }

     /* number of state variables */
    if (FindAndSkipComLine (stream, "Number of state variables", 
                            "Fix3_Num_Input",FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, 
                        "%d \n", 
                        &(tree_data->NbStates));

    if (readerror != 1)
    {        
        DR_Error("Could not find number of states in file %s! (Fix3_Num_Input)", FileName);
        goto RETURN;
    }
    
    if (mktvol_data->ModelChoice == FIX3_TMX)
    {
        /* number of multiq stdev iterations */
        if (FindAndSkipComLine (stream, "Number of MultiQ stdev", "Fix3_Num_Input",FileName) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, 
                        "%lf \n", 
                        &(mktvol_data->NbSigmaMQ));

        if (readerror != 1)
        {        
            DR_Error("Could not find MultiQ NbStd in file %s! (Fix3_NumInput)", FileName);
            goto RETURN;
        }

        /* number of multiq stdev iterations */
        if (FindAndSkipComLine (stream, "Number of MultiQ Nck", "Fix3_Num_Input",FileName) == FAILURE)
        {        
            goto RETURN;
        }
        readerror = fscanf (stream, 
                        "%lf \n", 
                        &(mktvol_data->NckMQ));

        if (readerror != 1)
        {        
            DR_Error("Could not find MultiQ Nck in file %s! (Fix3_NumInput)", FileName);
            goto RETURN;
        }
    }

    status = SUCCESS;

RETURN:

    if (stream != NULL)
    {
        fclose (stream);
    }
   
    return status;      
} /* Fix3_Num_Input_Common */


/*****  Fix3_Num_Input  *********************************************************/
/*
* See Fix3_Num_Input_Common.  Legacy (ORIGINAL mode) version that allows for deal file
* OverWriteString.
*/
int Fix3_Num_Input (MKTVOL_DATA    *mktvol_data,                 /* (O) Volatility data               */
                    FIX3_TREE_DATA *tree_data,                   /* (O) Tree data                     */
                    char            OverWriteString[6][MAXBUFF], /* (I) Overwrite strings             */
                    char const     *FileName)                    /* (I) File name including extension */ 
{
    int status = FAILURE;           /* Error status = FAILURE initially */
    int   readerror;                /* Reading error status             */

    if (Fix3_Num_Input_Common(mktvol_data,
                              tree_data,
                              FileName,
                              FALSE, /* fileMustExist */
                              FALSE  /* populateTreeNbSigmaMax */) != SUCCESS)
        goto RETURN;

    /* if overwrite exists in the deal file use it */
    if (strstr (OverWriteString[0], "nil") == NULL)
    {
        readerror = sscanf (OverWriteString[0],
                            "%d \n", 
                            &(tree_data->Ppy));

        if (readerror != 1)
        {        
            DR_Error("Could not read overwrite string for Ppy in file %s!", FileName);
            goto RETURN;
        }
    }

    /* if overwrites exist for smile parameters use the CET value */
    if ((strstr (OverWriteString[1], "nil") == NULL) &&  (mktvol_data->ModelChoice != FIX3_TMX))
    {
        if (strstr(OverWriteString[1], "N") != NULL)
        {
            mktvol_data->CetNbIter = 0;
        }       
        else if (strstr(OverWriteString[1], "L") != NULL)
        {
            mktvol_data->CetNbIter = 0;
        }
        else
        {
            double QLeft, QRight, FwdShift; /* Dummy local variables for smile parameters */

            readerror = sscanf (OverWriteString[1],
                                "%lf %lf %lf %d\n", 
                                &QLeft,
                                &QRight,
                                &FwdShift,
                                &(mktvol_data->CetNbIter));

            if (readerror != 4)
            {      
                DR_Error("Could not read overwrite string for Qweight! (Fix3_Num_Input)");
                goto RETURN;
            }       
        }
    }  /* if then else */

    status = SUCCESS;

RETURN:
   
    return (status);      
} /* Fix3_Num_Input */


/*****  Fix3_Num_Input_New  *********************************************************/
/*
* See Fix3_Num_Input_Common
*/
int Fix3_Num_Input_New (MKTVOL_DATA    *mktvol_data,                 /* (O) Volatility data               */
                        FIX3_TREE_DATA *tree_data,                   /* (O) Tree data                     */
                        char const     *FileName)                    /* (I) File name including extension */ 
{
    return Fix3_Num_Input_Common(mktvol_data, tree_data, FileName, TRUE, TRUE);
}


int Fix3_Env_Logging (
          char             Index[2][8],        /* (I) Calib index     */
          MKTVOL_DATA      *mktvol_data,       /* (I) Vol  data       */
          FIX3_TREE_DATA   *tree_data,         /* (I) Tree data       */
          T_CURVE          *t_curve)           /* (I) Zero curve      */
{
    char IndexL[MAXINDEX];   /* Local copy of index */
    int  status = FAILURE;
    int  i;
    char SoC;
    FILE *stream = NULL;
    int  is_yFix, is_yCms, is_BaseVol;

    strcpy(IndexL,Index[0]);

    is_yFix    = (strstr(IndexL, "yFix") != NULL);
    is_yCms    = (strstr(IndexL, "yCms") != NULL);
    is_BaseVol = ((strstr(IndexL, "m") != NULL) && (!is_yCms));

    stream = fopen ("numerics.dat","w");
    fprintf(stream, "# Number of periods per year\n");
    fprintf(stream, "%d\n", tree_data->Ppy);

    fprintf(stream, "# Number of standard deviations\n");
    fprintf(stream, "%d\n", tree_data->NbSigmaMax);

    fprintf(stream, "# Number of CET iter\n");
    fprintf(stream, "%d\n", mktvol_data->CetNbIter);

    fprintf(stream, "# Number of std dev for state variable\n");
    fprintf(stream, "%d\n", tree_data->NbStdDevStates);

    fprintf(stream, "# Number of states\n");
    fprintf(stream, "%d\n", tree_data->NbStates); 

    fclose(stream);
   
   /* create model parameters file */
    stream = fopen("ir_mpar_0.dat","w");
     /* Title */
    fprintf(stream, "####### FIX3 model parameters file ############\n");

    /* Number of factors */
    fprintf(stream, "### Section: VNFM parameters \n");
    fprintf(stream, "# Number of factors \n");
    fprintf(stream, "%d \n", (mktvol_data->NbFactor));

    /* Number of benchmark dates */
    fprintf(stream, "# Number of benchmark dates\n");
    fprintf(stream, "%d\n", 1);

    /* Number of benchmark dates */
    fprintf(stream, "# Dates, beta, alpha, rho\n");
    if (mktvol_data->NbFactor == 1)
        fprintf(stream, "%ld %lf %lf\n", YMDDateFromIRDate(mktvol_data->BaseDate), mktvol_data->Beta[0], mktvol_data->Alpha[0]);
    else if (mktvol_data->NbFactor == 2)
        fprintf(stream, "%ld %lf %lf %lf %lf %lf\n", YMDDateFromIRDate(mktvol_data->BaseDate), mktvol_data->Beta[0],
        mktvol_data->Beta[1],  mktvol_data->Alpha[0], mktvol_data->Alpha[1], mktvol_data->Rho[0]);
    else if (mktvol_data->NbFactor == 3)
        fprintf(stream, "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf \n", YMDDateFromIRDate(mktvol_data->BaseDate),
        mktvol_data->Beta[0], mktvol_data->Beta[1], mktvol_data->Beta[2], mktvol_data->Alpha[0], mktvol_data->Alpha[1],
        mktvol_data->Alpha[2], mktvol_data->Rho[0], mktvol_data->Rho[1], mktvol_data->Rho[2]);
    fprintf(stream, "###End Section\n");

    fprintf(stream," \n");

    /* Backbone */ 
    fprintf(stream, "### Section:Backbone\n");
    fprintf(stream, "#Backbone \n");
    fprintf(stream, "%lf\n",1- mktvol_data->Bbq);
    fprintf(stream, "###EndSection\n");

     fprintf(stream," \n");
    /* Smile */
    fprintf(stream, "###Section:Smile\n");
    fprintf(stream, "#Number of entries \n");
    fprintf(stream, "%d\n", 1);
    fprintf(stream, "#Date, qLeft, qRight, fwdShift \n");
    fprintf(stream, "%ld %lf %lf %lf\n", YMDDateFromIRDate(mktvol_data->BaseDate), 1-mktvol_data->QLeft, 1-mktvol_data->QRight,
                    mktvol_data->FwdShift);
    fprintf(stream, "###EndSection\n");

    fclose(stream);

   /* create model parameters file */
    if (is_BaseVol)
        SoC = 'C';
    else SoC = 'S';
    
    stream = fopen("ir_voldiag0_0.dat","w");

    /* Start date */
    fprintf(stream, "# Start date \n");
    fprintf(stream, "%ld \n", YMDDateFromIRDate(mktvol_data->BaseDate));

    fprintf(stream, "# Units \n");
    fprintf(stream, "%d \n", 0);

    fprintf(stream, "# NbEntries \n");
    fprintf(stream, "%d \n", mktvol_data->NbVol);

    fprintf(stream, "# Opt exp date, swap st, mat, vol \n");
    for (i = 0; i < mktvol_data->NbVol; i++)
        fprintf(stream, "%ld %ld %ld %12.8lf %c\n", YMDDateFromIRDate(mktvol_data->SwapSt[i]),
                YMDDateFromIRDate(mktvol_data->SwapSt[i]), YMDDateFromIRDate(mktvol_data->SwapMat[i]),
                mktvol_data->Vol[i]*100, SoC);

    fprintf(stream, "# Skip vols that fail to calibrate \n");
    fprintf(stream, "%s \n", mktvol_data->SkipFlag ? "Y" : "N");

    fclose(stream);

    stream = fopen ("ir_info_0.dat","w");
    fprintf(stream, "# Money market basis\n");
    fprintf(stream, "%d\n", t_curve[0].MMB);
    fprintf(stream, "# SwapDCC\n");
    if (strstr(t_curve[0].SwapDCC , "ACT") != NULL)
    fprintf(stream, "%s\n", "30/360");
    if (strstr(t_curve[0].SwapDCC, "365") != NULL)
    fprintf(stream, "%s\n", "ACT/365F");
    if (strstr(t_curve[0].SwapDCC, "360") != NULL)
    fprintf(stream, "%s\n", "ACT/360");
    fprintf(stream, "# SwapFreq\n");
    fprintf(stream, "%c\n", t_curve[0].SwapFreq);
    fprintf(stream, "# Diffusion curve index \n");
    fprintf(stream, "%d\n", tree_data->CvDiff);
    fprintf(stream, "# Discount curve index \n");
    fprintf(stream, "%d\n", tree_data->CvDisc);
    fclose(stream);

    stream = fopen ("ir_curve0_0.dat","w");
    EslPrintZeroCurve2(&t_curve[0], stream);
    fclose(stream);

    stream = fopen ("ir_curve1_0.dat","w");
    EslPrintZeroCurve2(&t_curve[1], stream);
    fclose(stream); 

    stream = fopen ("ir_curve2_0.dat","w");
    EslPrintZeroCurve2(&t_curve[2], stream);
    fclose(stream); 

    stream = fopen("summary.dat", "w");
    fprintf(stream, "### Environment Section\n");
    fprintf(stream, "# today's date\n");
    fprintf(stream, "%ld\n", mktvol_data->BaseDate);
    fprintf(stream, "### End Section\n\n");
    fprintf(stream, "### IR Section\n");
    fprintf(stream, "# Nb IR\n");
    fprintf(stream, "%d\n", 1);
    fprintf(stream, "### End Section\n");
    fclose(stream);

    stream = fopen("model.dat", "w");
    fprintf(stream, "# Model File\n");
    fprintf(stream, "# Pricing engine\n");
    fprintf(stream, "fix3\n");
    fprintf(stream, "# Model Choices: CLASSIC TIMEDEP SMD TMX\n");
    switch (mktvol_data->ModelChoice)
    {
    case FIX3_ORIGINAL:
        fprintf(stream, "ORIGINAL\n");
        break;
    case FIX3_CLASSIC:
        fprintf(stream, "CLASSIC\n");
        break;
    case FIX3_TIMEDEP:
        fprintf(stream, "TIMEDEP\n");
        break;
    case FIX3_SMD:
        fprintf(stream, "SMD\n");
        break;
    case FIX3_TMX:
        fprintf(stream, "TMX\n");
        break;
    case FIX3_E2Q:
        fprintf(stream, "E2Q\n");
        break;
    default:
        fprintf(stream, "\n");
    }
    fclose(stream);
    
    status = SUCCESS;
    return(status);

} /* Fix3_Env_Logging */
