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

/*****  Fix3_Param_Input_E2Q ********************************************************/
/*
*  	Read model parameters and check validity of input.
*/
int     Fix3_Param_Input_E2Q (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             FIX3_TREE_DATA*tree_data,                  /* (O) Tree data structure           */
             char const*   FileNameTreeNum,             /* (I) File name tree numerical pars */
             char const*   FileNameModelPar)            /* (I) File name model parameters    */
{
    int     readerror;          /* Reading error status             */
    int     status = FAILURE;   /* Error status = FAILURE initially */

    /* Read model parameters file */
    if (Fix3_Model_Input_E2Q (mktvol_data, 
                              tree_data,
                              FileNameModelPar) == FAILURE)
    {
        goto RETURN;
    }
            
    /* Read numerical parameters file */
    if (Fix3_Num_Input_New (mktvol_data, 
                            tree_data,
                            FileNameTreeNum) == FAILURE)
    {
        goto RETURN;
    }
  
    /* check validity of input */
    if (Fix3_Param_Check_E2Q (tree_data->NbFactor,
                              mktvol_data,
                              tree_data) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;
        
  RETURN:

    return status;
}




/*****  Fix3_Param_Check_E2Q ***********************************************/
/*
*  	Read term structure input and check validity of input.
*/
int     Fix3_Param_Check_E2Q (  
                int              NbFactor,     /* (I) Number of factors   */
                MKTVOL_DATA      *mktvol_data, /* (I) Volatility data     */
                FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure */
{
    int
        status = FAILURE;       /* Error status = FAILURE initially */
       
    char
        ErrorMsg[MAXBUFF];

    double  norm;

    /* Nb of factors */
    if (mktvol_data->NbFactor != NbFactor || tree_data->NbFactor != NbFactor)
    {
        DR_Error("Inconsistent number of factors! ");
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
    norm = 1.;
 
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

    if (!(mktvol_data->VolUnit == 1) && !(mktvol_data->VolUnit == 0))
        {
            DR_Error ("Incorrect value of VolUnit");
            goto RETURN;
        }
      

    status = SUCCESS;
        
    RETURN:
        
    return (status);

}  /* Fix3_Param_Check_E2Q */


/*****Fix3_Model_Input_E2Q  *********************************************************/
/*
* Read time-dependent mean-reversion, correlation, weight input      
*/
int Fix3_Model_Input_E2Q (MKTVOL_DATA    *mktvol_data, /* (O) Volatility data               */
                          FIX3_TREE_DATA *tree_data,   /* (0) Tree data                     */ 
                          char const     *FileName)    /* (I) File name including extension */ 
{
    int status = FAILURE;
    int     NbTDInp, NbFactor, SmileNbEntries;
    long    SmileDate;
    int     k;
   
    int     NbCorr; /* nb correlation needed */
    char    ErrorMsg[MAXBUFF];
    int     readerror;
    int     NbCet;
    FILE     *stream = NULL;

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
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Section:VNFM parameters", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
   

     /* Number of factors */
    if (FindAndSkipComLine (stream, 
                            "Number of factors", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    readerror = fscanf (stream, "%d\n", &NbFactor);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of time dependent inputs dates in file %s! "
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
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    readerror = fscanf (stream, "%d\n", &NbTDInp);

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
                "Could not read number of time dependent inputs dates in file %s! "
                "(Fix3_Model_Input_E2Q)", 
                FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /* if classic version of the tree one and only one time dep input date */
    if (NbTDInp != 1)
    {
        sprintf(ErrorMsg, 
                "Only one time dependent input date should be specified for classic version in file %s! "
                "(Fix3_Model_Input_E2Q)", 
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
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%ld ", &(mktvol_data->TDInpDate[0]));

    if (readerror != 1)
    {        
        sprintf(ErrorMsg, 
            "Could not read td input date in file %s! "
            "(Fix3_Model_Input_E2Q)", 
            FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    } 

	for (k = 0; k < NbFactor; k++)
	{
        readerror = fscanf (stream, "%lf ", &mktvol_data->Beta[k]);
        if (readerror != 1)
		{        
            sprintf(ErrorMsg, 
                "Could not read mean-reversion input for factor %d in file %s! "
                "(Fix3_Model_Input_E2Q)", 
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
                "(Fix3_Model_Input_E2Q)", 
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
                "(Fix3_Model_Input_E2Q)", 
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
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    

    /* Backbone section */
    if (FindAndSkipSectionLine (
                             1, /* skip mode */
                             stream, 
                            "Section:Backbone", 
                            "Fix3_Model_Input_E2Q", 
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
        sprintf (ErrorMsg, "Could not find Bbq in file %s! (Fix3_Model_Input_E2Q)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    mktvol_data->Bbq = 1. - mktvol_data->Bbq;


     if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


    /* SMILE SECTION */
    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Smile", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    

    /* Smile number of entries */
    if (FindAndSkipComLine (stream, 
                            "Number of entries", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream,
                        "%d \n", 
                         &SmileNbEntries);

    if (readerror != 1)
    {      
        sprintf (ErrorMsg, "Could not find Smile Nb entries in file %s! (Fix3_Model_Input_E2Q)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if(SmileNbEntries != 1)
    {
        sprintf (ErrorMsg, "Number of entries in Smile section has to be 1 in file %s! (Fix3_Model_Input_E2Q)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Read smile parameters */
    if (FindAndSkipComLine (stream, 
                            "Date, QLeft, QRight, FwdShift", 
                            "Fix3_Model_Input_E2Q", 
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
        sprintf (ErrorMsg, "Could not read smile info in file %s! (Fix3_Model_Input_E2Q)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (Dateok(SmileDate))
    {
        sprintf(ErrorMsg, "Incorrect format for smile date in file %s (Fix3_Model_Input_E2Q)!", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* At input level q=0 means log-normal whereas */
    /* internally q=0 means normal: we switch here */
    mktvol_data->QLeft  = 1. - mktvol_data->QLeft;
    mktvol_data->QRight = 1. - mktvol_data->QRight;

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
        {        
            goto RETURN;
        }



     /* SMILE SECTION */
   if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "E2Q parameters", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }



    /* Smile number of entries */
    if (FindAndSkipComLine (stream, 
                            "Parameters", 
                            "Fix3_Model_Input_E2Q", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }


     readerror = fscanf (stream, 
                    "%lf \t%lf  \n", 
                    &(mktvol_data->Amap), 
                    &(mktvol_data->Bmap));

    if (readerror != 2)
    {        
        sprintf (ErrorMsg, "Could not read E2Q info in file %s! (Fix3_Model_Input_E2Q)", FileName);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                         0, /* skip mode */
                         stream, 
                        "End Section", 
                        "Fix3_Model_Input_E2Q", 
                        FileName) == FAILURE)
    {        
        goto RETURN;
    }

	status = SUCCESS;
        
  RETURN:

    if (stream != NULL)
        fclose (stream);
        
    return status;

} /* Fix3_Model_Input_E2Q */
