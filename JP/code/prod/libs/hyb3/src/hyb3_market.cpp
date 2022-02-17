#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <hyb3_market.h>
#include <esl_market.h>
#include <string.h>
#include <stdlib.h>


/*************************************************
 *
 *  Read Raw Market Data
 *
 *  These functions read the "raw" wrapper type
 *  market environment data (<files>.dat) and 
 *  store in market data structures, independently
 *  of deal specifications
 *
 *  Interface functions read different market data
 *  for different wrapper types is....
 *
 *  HYB3Read<Tree mode>MarketW
 *
 * 
 **************************************************/


int HYB3ReadType3MarketW(
    Hyb3Market   *market)
{    
    int  status = FAILURE;
    const char *routine = "HYB3ReadType3MarketW:";

    FILE *streamRisk;

    /* read zero curves */
    market->mNbZeroCurve[FOR] = 3;
    market->mNbZeroCurve[DOM] = 3;

    if (Term_Input_W(&(market->mZeroCurve[FOR][0]),
                      "fzero.dat") == FAILURE)
    {        
        goto RETURN;
    }
    if (Term_Input_W(&(market->mZeroCurve[FOR][1]),
                      "fdisczero.dat") == FAILURE)
    {        
        goto RETURN;
    }
    /* if riskzero.dat doesn't exist use disczero.dat */
    streamRisk = fopen ("friskzero.dat", "r");

    if (streamRisk == NULL)
    {
        if (Term_Input_W(&(market->mZeroCurve[FOR][2]),
                      "fdisczero.dat") == FAILURE)
        {        
            goto RETURN;
        }
    }
    else
    {
        if (Term_Input_W(&(market->mZeroCurve[FOR][2]),
                      "friskzero.dat") == FAILURE)
        {        
            goto RETURN;
        }
        fclose(streamRisk);
    }
    
    if (Term_Input_W(&(market->mZeroCurve[DOM][0]),
                      "dzero.dat") == FAILURE)
    {        
        goto RETURN;
    }
    if (Term_Input_W(&(market->mZeroCurve[DOM][1]),
                      "ddisczero.dat") == FAILURE)
    {        
        goto RETURN;
    }
    /* if riskzero.dat doesn't exist use disczero.dat */
    streamRisk = fopen ("driskzero.dat", "r");

    if (streamRisk == NULL)
    {
        if (Term_Input_W(&(market->mZeroCurve[DOM][2]),
                      "ddisczero.dat") == FAILURE)
        {        
            goto RETURN;
        }
    }
    else
    {
        if (Term_Input_W(&(market->mZeroCurve[DOM][2]),
                      "driskzero.dat") == FAILURE)
        {        
            goto RETURN;
        }
        fclose(streamRisk);
    }

    /* read base/swap vols */
    if (EslReadVolsW(&(market->mBaseVol[FOR]),
                     &(market->mSwapVol[FOR]),
                     "fbasevol.dat",
                     "fswapvol.dat") == FAILURE)
    {
        goto RETURN;
    }
    if (EslReadVolsW(&(market->mBaseVol[DOM]),
                     &(market->mSwapVol[DOM]),
                     "dbasevol.dat",
                     "dswapvol.dat") == FAILURE)
    {
        goto RETURN;
    }

    /* read FX vols and smile */
    if (HYB3ReadFXEnvW(&(market->mFXVolatility),
                       &(market->mFXSmile),
                       "FXVolatility.dat",
                       "fxsmile_0.dat") == FAILURE)
    {
        goto RETURN;
    }

    /* read correlation */
    if (HYB3ReadCorrelation(&(market->mCorrelation),
                            "correlation.dat") == FAILURE)
    {
        goto RETURN;
    }

    /* read model parameters */
    if (HYB3ReadModelParams(&(market->mModelParameter[FOR]),
                            1,
                            &(market->mModelParamsSupplied[FOR]),
                            FALSE,  /* don't read Q parameters from market file */
                            "fmodelParameters.dat") == FAILURE)
    {
        goto RETURN;
    }
    if (HYB3ReadModelParams(&(market->mModelParameter[DOM]),
                            1,
                            &(market->mModelParamsSupplied[DOM]),
                            FALSE,  /* don't read Q parameters from market file */
                            "dmodelParameters.dat") == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
        DR_Error("%s failed", routine);        

    return status;

}

int HYB3ReadFXEnvW(
    FXVOLATILITY_DATA  *FXVol,
    FXSMILE_DATA       *FXSmile,
    char const         *FXVolFilename,
    char const         *FXSmileFilename)
{
    int  i;
    int  status = FAILURE;
    const char *routine = "HYB3ReadFXEnvW: ";
    char string[MAXBUFF];
    int  readerror;
          
    FILE
        *stream        = NULL,
        *streamFXSmile = NULL; 


    stream = fopen (FXVolFilename, "r");

    if (stream == NULL)
    {
        DR_Error("%s Could not open file %s! ", 
                 routine,
                 FXVolFilename);
          
        goto RETURN;
          
    }  /* if */
          
    fgets  (string, 80, stream);    /* Read the comment line */
    readerror = fscanf (stream, "%ld \n", &(FXVol->ValueDate));
    if ((readerror == 0) || (readerror == EOF))
    {          
        DR_Error("%s Could not read file %s! ",
                  routine,
                  FXVolFilename);

        goto RETURN;
    }
     
    fgets  (string, 80, stream);
    readerror = fscanf (stream, "%lf \n", &(FXVol->FXSpotRate));
    if ((readerror == 0) || (readerror == EOF))
    {          
        DR_Error("%s Could not read FX spot in file %s! ",
                 routine,
                 FXVolFilename);
          
        goto RETURN;
    }
     
    fgets ( string, 80, stream);  /* Ignore this input */
    fgets ( string, 80, stream);

    fgets ( string, 80, stream);
    readerror = fscanf (stream, "%d \n", &(FXVol->NbBaseVols));
    if ((readerror == 0) || (readerror == EOF))
    {          
        DR_Error("%s Could not read Nb base vols in file %s! ",
                 routine,
                 FXVolFilename);
          
        goto RETURN;
    }

    if (FXVol->NbBaseVols > MAXNBDATE)
    {
        DR_Error("%s Nb of vols exceeds max limit of %d in file %s! ",
                 routine,
                 MAXNBDATE, 
                 FXVolFilename);

        goto RETURN;
    }

    fgets ( string, 80, stream);
    for (i = 0; i < FXVol->NbBaseVols; i++)
    {
        readerror = fscanf (stream, "%ld \t%lf \n",
                            &(FXVol->BaseVolDates[i]),
                            &(FXVol->BaseVols[i]));

        if ((readerror == 0) || (readerror == EOF))
        {          
            DR_Error("%s Could not read file %s! ",
                     routine,
                     FXVolFilename);

            goto RETURN;
        }
    }

    /* Input Spot Vol section */
    fgets  (string, 80, stream);    /* comment line */
    readerror = fscanf (stream, "%d \n", &(FXVol->NbSpotVols));
    if (readerror != 1) 
    {          
        DR_Error("%s Could not read Nb of Input Spot Vol points "
                 "in file %s! ", 
                 routine, 
                 FXVolFilename);
        goto RETURN;
    }

    if (FXVol->NbSpotVols > MAXNBDATE)
    {
         DR_Error("%s Nb of Inp Spot Vol exceeds max limit of %d "
                   "in file %s! ", 
                   routine, 
                   MAXNBDATE, 
                   FXVolFilename);
         
         goto RETURN;
    }

    if (FXVol->NbSpotVols > 0)
    {
        fgets  (string, 80, stream);
        for (i = 0; i < FXVol->NbSpotVols; i++)
        {
            readerror = fscanf (stream, "%ld \t%lf \n",
                                &(FXVol->SpotVolDates[i]),
                                &(FXVol->SpotVols[i]));

            if (readerror != 2)
            {          
                DR_Error("%s Could not read input Spot Vol"
                         " dates/values in file %s! ",
                         routine,
                         FXVolFilename);
               goto RETURN;
            }
        }
    }
    fclose (stream);

    /* Read Smile data file */
    if (strlen(FXSmileFilename) == 0)
        streamFXSmile = NULL;
    else
        streamFXSmile = fopen (FXSmileFilename, "r");

    /* if fx smiel file doesn't exist (it is optional), enter default values for smile */
    if (streamFXSmile == NULL)
    {
        FXSmile->NbSmileDates = 1;

        /* add arbitary number of days after the value date to ensure the smile date is
           > today as in some cases valueDate = Today. - chose to use 5 days but could have
           use any number */
        FXSmile->SmileDates[0] = Nxtday(FXVol->ValueDate, 5);
        FXSmile->A1[0] = 0.0;
        FXSmile->A2[0] = 0.0;
        FXSmile->A3[0] = 0.0;
    }
    else
    {
        /* read Nb FX Smile Param lines */
        if (FindAndSkipComLine (streamFXSmile, 
                                "Nb FX Param Lines",
                                routine,
                                FXSmileFilename) == FAILURE) goto RETURN;

        readerror = fscanf (streamFXSmile, "%d \n", &(FXSmile->NbSmileDates));
        if (readerror != 1)
        {
            DR_Error("%s Cannot read Nb FX smile param lines in file %s! ",
                     routine,
                     FXSmileFilename);
            goto RETURN;
        }

        if ((FXSmile->NbSmileDates < 0) ||
            (FXSmile->NbSmileDates > MAXNBDATE))
        {
            DR_Error("%s Nb FX smile param lines in file %s out of range!",
                     routine,
                     FXSmileFilename);
            goto RETURN;
        }

        if (FindAndSkipComLine (streamFXSmile, 
                                "FX Smile Params",
                                routine,
                                FXSmileFilename) == FAILURE) goto RETURN;

        for (i=0; i < FXSmile->NbSmileDates; i++)
        {
            readerror = fscanf(streamFXSmile, "%ld\t%lf\t%lf\t%lf \n",
                               &(FXSmile->SmileDates[i]),
                               &(FXSmile->A1[i]),
                               &(FXSmile->A2[i]),
                               &(FXSmile->A3[i]));

            if (readerror != 4)
            {          
                DR_Error("%s Cannot read FX smile params in file %s! ",
                         routine,
                         FXSmileFilename);
                goto RETURN;
            }
        }
        fclose(streamFXSmile);
    }

    status = SUCCESS;
          
    RETURN:
          
    return (status);

}


/* ???? double check that these correlation are in the correct order */
int HYB3ReadCorrelation(
    CORRELATION_DATA   *correl,
    char const         *filename)
{          
    int  status = FAILURE;
    const char *routine = "HYB3ReadCorrelation:";
    int   readerror;
    
    FILE *stream = NULL;

    stream = fopen (filename, "r");             /* Open the correlation data file */

    if (stream == NULL)
    {
         DR_Error("%s Could not open file %s! ", 
                  routine,
                  filename);
         goto RETURN;
    }
          
     if (FindAndSkipComLine (stream, "corr dom vs. foreign IR", 
                                      routine, filename) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(correl->CorrIR));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("%s Could not read file %s! ",
                   routine,
                   filename);
          
          goto RETURN;
     }


     if (FindAndSkipComLine (stream, "corr domestic IR vs. FX", 
                                     routine, filename) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(correl->CorrDomIRFX));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("%s Could not read file %s! ",
                   routine,
                   filename);
          
          goto RETURN;
     }

     if (FindAndSkipComLine (stream, "corr foreign IR vs. FX", 
                                     routine, filename) == FAILURE)
     {          
          goto RETURN;
     }
     readerror = fscanf (stream, "%lf \n", &(correl->CorrForIRFX));
     if ((readerror == 0) || (readerror == EOF))
     {          
          DR_Error("%s Could not read file %s! ",
                   routine,
                   filename);
          
          goto RETURN;
     }

     fclose (stream);

     status = SUCCESS;
          
     RETURN:
          
     return (status);

}


int HYB3ReadModelParams(
    MODELPARAMETERS_DATA *modelParams,
    int                  NbFactor,        // always 1 factor IR for hyb3, except for hyb2+1 mode
    int                  *paramsSupplied, // (O) true/false
    int                  readSmileParams, // (I) true/false
    char const           *filename)
{
    int     i; 
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    const char    *routine = "HYB3ReadModelParams:";
    char    *getserror;         /* Reading error for fgets          */
    char    string[MAXBUFF];
    FILE    *stream = NULL;

    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    stream = fopen (filename, "r");

    /* model parameters file is optional, so if it is not supplied, simply
       set the paramsSupplied variable to indicate this and return */
    if (stream == NULL)
    {
        *paramsSupplied = FALSE;
        return SUCCESS;
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

            if (FindAndSkipComLine (stream, "one factor mean reversion", 
                                    routine,
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->OneFactorMR));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find mean reversion in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "one factor weight", 
                                    routine, 
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->OneFactorVol));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find factor weight in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "one factor ppy", 
                                    routine,
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%d \n", 
                                &(modelParams->OneFactorPPY));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find Ppy in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            /* No overwrite for Qweight: we will need to skip two and three factor param.
             * in the file to access Qweight data at the end of the file
             */
            if (readSmileParams)
            {
                /* Skip two factor parameters in the file */
                for (i = 0; i < 12; i++)
                {
                    getserror = fgets (string, MAXBUFF, stream);
                    if (getserror  == NULL)
                    {        
                        DR_Error("%s Could not find two factor parameters in file %s! ",
                                 routine,
                                 filename);
                        
                        goto RETURN;
                    }
                }
            }
        }
        else if (NbFactor == 2)
        {
            /* Skip one factor parameters in the file */
            for (i = 0; i < 6; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("%s Could not find one factor parameters in file %s! ",
                             routine,
                             filename);
                    
                    goto RETURN;
                }
             }

            /*
             *  Read mean reversion in parameter file.
             */

            if (FindAndSkipComLine (stream, "two factor mean reversion", 
                                    routine,
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->TwoFactorMR1));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find mean reversion in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }
            if (FindAndSkipComLine (stream, "two factor mean reversion", 
                                    routine,
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->TwoFactorMR2));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find mean reversion in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }


            if (FindAndSkipComLine (stream, "two factor weight", 
                                    routine, 
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->TwoFactorVol1));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find factor weight in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor weight", 
                                    routine, 
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->TwoFactorVol2));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find factor weight in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "two factor correlation", 
                                    routine, 
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%lf \n", 
                                &(modelParams->TwoFactorCorr));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find factor correlation in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }
            if (FindAndSkipComLine (stream, "two factor ppy", 
                                    routine,
                                    filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream, 
                                "%d \n", 
                                &(modelParams->TwoFactorPPY));

            if (readerror != 1)
            {        
                DR_Error("%s Could not find Ppy in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }
        }
        else 
        {
            DR_Error("%s Hyb3 only supports 1 and 2 factor interest rates in file %s! ",
                     routine,
                     filename);
            goto RETURN;
        }

        /* No overwrite for Qweight: look for it in the parameter file */
        if (readSmileParams)
        {
             /* Skip three factor parameters in the file */
            for (i = 0; i < 20; i++)
            {
                getserror = fgets (string, MAXBUFF, stream);
                if (getserror  == NULL)
                {        
                    DR_Error("%s Could not find three factor parameters in file %s! ",
                             routine,
                             filename);
                    
                    goto RETURN;
                }
            }

            if (FindAndSkipComLine (stream, "QLeft", routine, filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                    "%lf \n", 
                                    &(modelParams->QLeft));

            if (readerror != 1)
            {      
                DR_Error("%s Could not find QLeft in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "QRight", routine, filename) == FAILURE)
            {        
                goto RETURN;
            }

                readerror = fscanf (stream,
                                    "%lf \n", 
                                    &(modelParams->QRight));

            if (readerror != 1)
            {      
                DR_Error("%s Could not find QRight in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "FwdShift", routine, filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                    "%lf \n", 
                                    &(modelParams->FwdShift));

            if (readerror != 1)
            {      
                DR_Error("%s Could not find FwdShift in file %s! ", 
                         routine,
                         filename);
                
                goto RETURN;
            }

            if (FindAndSkipComLine (stream, "CetNbIter", routine, filename) == FAILURE)
            {        
                goto RETURN;
            }

            readerror = fscanf (stream,
                                    "%d \n", 
                                    &(modelParams->CetNbIter));

            if (readerror != 1)
            {      
                DR_Error("%s Could not find CetNbIter in file %s! ",
                         routine,
                         filename);
                
                goto RETURN;
            }

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

/********************************************
 *
 *  Populate underlying HYB3 tree structures
 *
 *  Takes model parameters and overwrites
 *  from the deal file and populates
 *  MKTVOL_DATA, FX_DATA, TREE_DATA
 *
 ********************************************/

int HYB3ParamOverwrites(
    MKTVOL_DATA           *mktvol_data,              // (O) Volatility data
    HYB3_TREE_DATA        *tree_data,                // (O) Tree data structure
    int                   NbFactor,                  // (I) Number of factors
    MODELPARAMETERS_DATA* modelParameters_data,
    int                   modelParametersSupplied,
    char                  OverWriteString[5][MAXBUFF])
{

    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */
    char    FileName[MAXBUFF];
    
    FileName[0] = '\0';


    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    //
    //  If there is no parameter file use the overwrite strings.
    //

    if (!modelParametersSupplied)
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
                DR_Error("Could not find file %s: Q overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            // At input level q=0 means log-normal whereas
            // internally q=0 means normal: we switch here
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
                DR_Error("Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \n", 
                                &(mktvol_data->Beta[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            // Fill in the unused parameters with N/A values

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
                DR_Error("Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]));

            if (readerror != 2)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \n", 
                                &(mktvol_data->Rho[0]));

            if (readerror != 1)
            {        
                DR_Error("Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                
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
                DR_Error("Could not find file %s: factor weight overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[3],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Beta[0]),
                                &(mktvol_data->Beta[1]),
                                &(mktvol_data->Beta[2]));

            if (readerror != 3)
            {        
                DR_Error("Could not find file %s: mean reversion overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }

            readerror = sscanf (OverWriteString[4],
                                "%lf \t%lf \t%lf \n", 
                                &(mktvol_data->Rho[0]),
                                &(mktvol_data->Rho[1]),
                                &(mktvol_data->Rho[2]));

            if (readerror != 3)
            {        
                DR_Error("Could not find file %s: correlation overwrite required! (Param_Input)", FileName);
                
                goto RETURN;
            }
        } 
    }
    else
    {        
        //
        //  Parameter file exists: use the overwrite strings if they are not 
        //  "nil", otherwise use the values in the parameter file.
        //
        if (NbFactor == 1)
        {
            // mean reversion
            mktvol_data->Beta[0] = modelParameters_data->OneFactorMR;

            // overwrite if requried
            if (strstr (OverWriteString[3], "nil") == NULL)                     /* Need to use strstr() here, strcmp won't do */
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \n", 
                                    &(mktvol_data->Beta[0]));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for mean reversion! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            // factor weight
            mktvol_data->Alpha[0] = modelParameters_data->OneFactorVol;

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \n", 
                                    &(mktvol_data->Alpha[0]));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for factor weight! (Param_Input)");
                    
                    goto RETURN;
                }
            }


            /* No overwrite for Qweight: look for it in the parameter file */
            if (strstr (OverWriteString[1], "nil") != NULL)
            {
                mktvol_data->QLeft    = 1.0 - modelParameters_data->QLeft;
				mktvol_data->QRight   = 1.0 - modelParameters_data->QRight;
                mktvol_data->FwdShift = modelParameters_data->FwdShift;
                mktvol_data->CetNbIter= modelParameters_data->CetNbIter;
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
                        DR_Error("Could not read overwrite string for Qs ! (Param_Input)");
                        
                        goto RETURN;
                    }       

                    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                    mktvol_data->QRight = 1. - mktvol_data->QRight;
                }
            }  /* if then else */
 
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
            mktvol_data->Beta[0] = modelParameters_data->TwoFactorMR1;
            mktvol_data->Beta[1] = modelParameters_data->TwoFactorMR2;

            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]));

                if (readerror != 2)
                {        
                    DR_Error("Could not read overwrite string for mean reversion! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            mktvol_data->Alpha[0] = modelParameters_data->TwoFactorVol1;
            mktvol_data->Alpha[0] = modelParameters_data->TwoFactorVol2;

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]));

                if (readerror != 2)
                {        
                    DR_Error("Could not read overwrite string for factor weight! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            mktvol_data->Rho[0] = modelParameters_data->TwoFactorCorr;
            
            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \n", 
                                    &(mktvol_data->Rho[0]));

                if (readerror != 1)
                {        
                    DR_Error("Could not read overwrite string for correlation! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (strstr (OverWriteString[1], "nil") != NULL)
            {
                DR_Error("No smile parameters in market model parameters");
                goto RETURN;
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
                        DR_Error("Could not read overwrite string for Qs! (Param_Input)");
                        
                        goto RETURN;
                    }       

                    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                    mktvol_data->QRight = 1. - mktvol_data->QRight;
                }
            }  /* if then else */

            mktvol_data->Alpha[2] = -999.;
            mktvol_data->Beta[2]  = -999.;
            mktvol_data->Rho[1]   = -999.;
            mktvol_data->Rho[2]   = -999.;
        }
        else if (NbFactor == 3)
        {
            
            mktvol_data->Beta[0] = modelParameters_data->ThreeFactorMR1;
            mktvol_data->Beta[1] = modelParameters_data->ThreeFactorMR2;
            mktvol_data->Beta[2] = modelParameters_data->ThreeFactorMR3;

            if (strstr (OverWriteString[3], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[3],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Beta[0]),
                                    &(mktvol_data->Beta[1]),
                                    &(mktvol_data->Beta[2]));

                if (readerror != 3)
                {        
                    DR_Error("Could not read overwrite string for mean reversion! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            mktvol_data->Alpha[0] = modelParameters_data->ThreeFactorVol1;
            mktvol_data->Alpha[1] = modelParameters_data->ThreeFactorVol2;
            mktvol_data->Alpha[2] = modelParameters_data->ThreeFactorVol3;

            if (strstr (OverWriteString[2], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[2], 
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Alpha[0]),
                                    &(mktvol_data->Alpha[1]),
                                    &(mktvol_data->Alpha[2]));

                if (readerror != 3)
                {        
                    DR_Error("Could not read overwrite string for factor weight! (Param_Input)");
                    
                    goto RETURN;
                }
            }
            
            mktvol_data->Rho[0] = modelParameters_data->ThreeFactorCorr12;
            mktvol_data->Rho[1] = modelParameters_data->ThreeFactorCorr13;
            mktvol_data->Rho[2] = modelParameters_data->ThreeFactorCorr23;

            
            if (strstr (OverWriteString[4], "nil") == NULL)
            {
                readerror = sscanf (OverWriteString[4],
                                    "%lf \t%lf \t%lf \n", 
                                    &(mktvol_data->Rho[0]),
                                    &(mktvol_data->Rho[1]),
                                    &(mktvol_data->Rho[2]));

                if (readerror != 3)
                {        
                    DR_Error("Could not read overwrite string for correlation! (Param_Input)");
                    
                    goto RETURN;
                }
            }

            if (strstr (OverWriteString[1], "nil") != NULL)
            {
                DR_Error("No smile parameters in market model parameters");
                goto RETURN;
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
                        DR_Error("Could not read overwrite string for Qweight! (Param_Input)");
                        
                        goto RETURN;
                    }       
        
                    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
                    mktvol_data->QRight = 1. - mktvol_data->QRight;
                }
            }  /* if then else */
        }  /* if then else */
    }  /* if then else */
    
    /* Check validity of input */
    if (Hyb3_Param_Check (NbFactor,
                          mktvol_data,
                          tree_data) == FAILURE)              
    {        
        goto RETURN;
    }
        
    status = SUCCESS;
        
    RETURN:

    return (status);

}  /* Param_Input */



int HYB3CorrelOverwrites(
    FX_DATA            *fx_data,                    // (O) Fx data
    char               OverWriteString[3][MAXBUFF], // (I) Overwrite strings
    CORRELATION_DATA   *correlation_data)
{
    double temp[3];

    int status = FAILURE;             /* Error status = FAILURE initially */

              
    temp[0] = correlation_data->CorrIR;
    temp[1] = correlation_data->CorrDomIRFX;
    temp[2] = correlation_data->CorrForIRFX;


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
    }
          
    if (strcmp (OverWriteString[1], "nil"))
    {
        fx_data->Rho[1][0] = atof (OverWriteString[1]);          
    }
          
    if (strcmp (OverWriteString[2], "nil"))
    {                    
        fx_data->Rho[2][0] = atof (OverWriteString[2]);                               
    }
          

    if (Hyb3_Correl_Check_WType3 (fx_data) == FAILURE)                  /* Check validity of input */
    {          
         goto RETURN;
    }
          

    status = SUCCESS;
          
    RETURN:
          
    return (status);

}


int HYB3FXSmileOverwrites(
    FX_DATA            *fx_data,                 // (O) Fx data
    char               OverWriteString[MAXBUFF], // (I) FX Owrite str
    FXVOLATILITY_DATA  *Volatility_data,         // (I) FXVols file name
    FXSMILE_DATA       *Smile_data,              // (I) FX smile file name
    char               NbFXSmileOWS[MAXBUFF],    // (I)
    char               FXSmileParamOWS[MAXNBDATE][MAXBUFF])
{
     long
          Year[2],
          Month[2],
          Day[2];
     int
          i,
          readerror,                 /* Reading error status             */
          status = FAILURE;          /* Error status = FAILURE initially */


     fx_data->Today = Volatility_data->ValueDate;
          
     fx_data->ValueDate = fx_data->Today;
     fx_data->SpotDays  = 0;

     fx_data->Spot = Volatility_data->FXSpotRate;
     // ignore BaseVolFreq entry

     fx_data->NbVol = Volatility_data->NbBaseVols;
     if (fx_data->NbVol > MAXNBDATE)
     {
          DR_Error("Nb of vols exceeds max limit of %d in "
                             "(Fx_Input_W_WithSmile_DRI)", MAXNBDATE);
          
          goto RETURN;
     }
     
     Dsplit(fx_data->ValueDate, /* Split value date into month, day and year */
              &(Month[0]), 
              &(Day[0]), 
              &(Year[0]));

     for (i = 1; i <= fx_data->NbVol; i++)  /* Composite vols -> 1 offset*/
     {
        fx_data->VolDate[i] = Volatility_data->BaseVolDates[i-1];
        fx_data->FxVol[i] = Volatility_data->BaseVols[i-1]/100.0;

     }  /* for i */

     /* Repeat the first values for interpolation of pseudocomposite vols */
     /* in Get_TreeSpotVols                                               */

    if (fx_data->NbVol > 0)
    {
        fx_data->VolDate[0] = fx_data->ValueDate;
        fx_data->FxVol[0]   = fx_data->FxVol[1];
    }
     

    fx_data->NbInpSpotVol = Volatility_data->NbSpotVols;
    if (fx_data->NbInpSpotVol > MAXNBDATE)
    {
         DR_Error("Nb of Inp Spot Vol exceeds max limit of %d "
                   "in ! (Fx_Input_W)", MAXNBDATE);
         
         goto RETURN;
    }

    if (fx_data->NbInpSpotVol > 0)
    {
       for (i = 0; i < fx_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
       {
             fx_data->InpSpotVolDate[i] = Volatility_data->SpotVolDates[i];
             fx_data->InpSpotVol[i] = Volatility_data->SpotVols[i]/100.0;
       }
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
     else if (strcmp(OverWriteString,"last") == 0) /* 2.Cut off at last level */
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
          if (ABS(fx_data->FxCutOffLevel) < TINY) /* The FAILURE return */
          {
              DR_Error("Unable to convert FX cut off value!\n");
              goto RETURN;
          }
          fx_data->FxCutOffLevel /= 100.;
          fx_data->FxBootStrapMode = FX_CONSTANT_SPOT_VOL;
          
     }
     

    // Read FX Smile data
    

    if (strstr(NbFXSmileOWS, "nil") != NULL)
    {
        if (strstr(FXSmileParamOWS[0], "nil") == NULL)
        {
            DR_Error("Fx_Input_W: invalid FX Smile Params input\n");
            goto RETURN;
        }

        // if smile data structure is NULL, turn off the FX smile by setting all
        // the smile parameters to 0.0
        if (Smile_data == NULL)
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
            fx_data->NbSmilePt = Smile_data->NbSmileDates;
            if ((fx_data->NbSmilePt < 0) ||
                (fx_data->NbSmilePt > MAXNBDATE))
            {
                DR_Error("Fx_Input_W: Nb FX smile param lines out of range!");
                goto RETURN;
            }

            for (i=0; i < fx_data->NbSmilePt; i++)
            {
                fx_data->SmileDate[i] = Smile_data->SmileDates[i];
                fx_data->a1[i] = Smile_data->A1[i];
                fx_data->a2[i] = Smile_data->A2[i];
                fx_data->a3[i] = Smile_data->A3[i];
            }  
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
            DR_Error("Fx_Input_W: Cannot read Nb FX Smile Param OWS");
            goto RETURN;
        }

        if ((fx_data->NbSmilePt < 0) ||
            (fx_data->NbSmilePt > MAXNBDATE))
        {
            DR_Error("Fx_Input_W: Nb FX smile param lines OWS out of range!");
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

            if (readerror != 4)
            {          
                DR_Error("Fx_Input_W: Cannot read FX smile params OWS");
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
          
     return (status);

}
