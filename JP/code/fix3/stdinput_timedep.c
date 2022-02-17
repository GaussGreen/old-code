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



/*****  Fix3_Param_Input_TimeDep  ********************************************************/
/*
*  	Read model parameters and check validity of input.
*/
int     Fix3_Param_Input_TimeDep (   
             MKTVOL_DATA   *mktvol_data,                /* (O) Volatility data               */
             FIX3_TREE_DATA*tree_data,                  /* (O) Tree data structure           */
             char const*   FileNameTreeNum,             /* (I) File name tree numerical pars */
             char const*   FileNameModelPar)            /* (I) File name model parameters    */
{

    int     status = FAILURE;   /* Error status = FAILURE initially */
  
    /* Read model parameters file */
    if (Fix3_Model_Input_TimeDep (mktvol_data, 
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
    if (Fix3_Param_Check_TimeDep (tree_data->NbFactor,
                     mktvol_data,
                     tree_data) == FAILURE)
    {
        goto RETURN;
    }

    status = SUCCESS;
        
  RETURN:
        
    return status;

}  /* Fix3_Param_Input_TimeDep */



/*****Fix3_Model_Input_TimeDep  *********************************************************/
/*
* Read time-dependent mean-reversion, correlation, weight input      
*/
int Fix3_Model_Input_TimeDep (MKTVOL_DATA    *mktvol_data,   /* (O) Volatility data               */
                              FIX3_TREE_DATA *tree_data,     /* (0) Tree data                     */
                              char const     *FileName)      /* (I) File name including extension */           
{
    int     NbTDInp, NbFactor;
    long    currExp, prevExp;
    int     t, k;

    int     NbCorr; /* nb correlation needed */
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
                             "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Section : VNFM parameters", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
   
     /* Number of factors */
    if (FindAndSkipComLine (stream, 
                            "Number of factors", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbFactor);
    if (readerror != 1)
    {        
        DR_Error("Could not read number of time dependent inputs dates in file %s! "
                 "(Fix3_Model_Input_TimeDep)", 
                 FileName);
        goto RETURN;
    }

    tree_data->NbFactor = mktvol_data->NbFactor = NbFactor;
    
    /*compute number of correlations */
    NbCorr = NbFactor * (NbFactor - 1) / 2;

    /* Number of time dependent input  dates */
    if (FindAndSkipComLine (stream, 
                            "Number of benchmark dates", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d\n", &NbTDInp);
    if (readerror != 1)
    {        
        DR_Error("Could not read number of time dependent inputs dates in file %s! "
                 "(Fix3_Model_Input_TimeDep)", 
                 FileName);
        goto RETURN;
    }

    /* at least one time dependent input date needs to be specified */
    if (NbTDInp < 1)
    {
        DR_Error("Invalid input for number of time dependent input dates in file %s! "
                 "(Fix3_Model_Input_TimeDep)", 
                 FileName);
        goto RETURN;
    }

    if (NbTDInp > MAXNBTD)
    {
         DR_Error("Maximum number of time dependent input dates in file %s! surpassed"
                 "(Fix3_Model_Input_TimeDep)", 
                 FileName);
        goto RETURN;
    }


    mktvol_data->NbTDInp = NbTDInp;
   
    if (NbFactor == 1)
    {
        mktvol_data->AlphaTD[1][0] = -999.;
        mktvol_data->AlphaTD[2][0] = -999.;
        mktvol_data->BetaTD[1][0]  = -999.;
        mktvol_data->BetaTD[2][0]  = -999.;
        mktvol_data->RhoTD[0][0]   = -999.;
        mktvol_data->RhoTD[1][0]   = -999.;
        mktvol_data->RhoTD[2][0]   = -999.;
    }
    else if (NbFactor == 2)
    {
        mktvol_data->AlphaTD[2][0] = -999.;
        mktvol_data->BetaTD[2][0]  = -999.;
        mktvol_data->RhoTD[1][0]   = -999.;
        mktvol_data->RhoTD[2][0]   = -999.;
    }

    /* Term structure of weights, mean-reversion, correlation */
    if (FindAndSkipComLine (stream, 
                            "Term structure of , weights,mean-reversion, correlation", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    prevExp = 0;

    for ( t = 0; t < NbTDInp; t++)
    {
        readerror = fscanf (stream, "%ld ", &currExp);

        if (readerror != 1)
        {        
            DR_Error("Could not read td input date %d in file %s! "
                     "(Fix3_Model_Input_TimeDep)", 
                     t,
                     FileName);
            goto RETURN;
        }

        if (t > 0 && currExp <= prevExp)
        {
            DR_Error("Invalid input for date %d in file %s! "
                     "(Fix3_Model_Input_TimeDep)",
                     t,
                     FileName);
            goto RETURN;
        }

        mktvol_data->TDInpDate[t] = currExp;
        prevExp = currExp;

		for (k = 0; k < NbFactor; k++)
		{
            readerror = fscanf (stream, "%lf ", &mktvol_data->BetaTD[k][t]);
            if (readerror != 1)
			{        
                DR_Error("Could not read mean-reversion input for factor %d at date %d in file %s! "
                         "(Fix3_Model_Input_TimeDep)", 
					     k,
                         t,
                         FileName);
                goto RETURN;
            }
		}

        for (k = 0; k < NbFactor; k++)
		{
            readerror = fscanf (stream, "%lf ", &mktvol_data->AlphaTD[k][t]);
            if (readerror != 1)
			{        
                DR_Error("Could not read factor weight input for factor %d at date %d in file %s! "
                         "(Fix3_Model_Input_TimeDep)", 
					     k,
                         t,
                         FileName);
                goto RETURN;
			}
		}

       	for (k = 0; k < NbCorr; k++)
		{
            readerror = fscanf (stream, "%lf ", &mktvol_data->RhoTD[k][t]);
            if (readerror != 1)
			{        
                DR_Error("Could not read correlation input for factor %d at date %d in file %s! "
                         "(Fix3_Model_Input_TimeDep)", 
					     k,
                         t,
                         FileName);
                goto RETURN;
			}
		}
        fscanf (stream, "\n");
    } /* for t */

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Backbone section */
    if (FindAndSkipSectionLine (
                             1,/* skip mode */
                             stream, 
                            "Section:Backbone", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    if (FindAndSkipComLine (stream, "Backbone", "Fix_Model_Input", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%lf \n", &(mktvol_data->Bbq));
    if (readerror != 1)
    {      
        DR_Error("Could not find Bbq in file %s! (Fix3_Model_Input_TimeDep)", FileName);
        goto RETURN;
    }

    mktvol_data->Bbq = 1. - mktvol_data->Bbq;

    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    
    /* SMILE SECTION */
    if (FindAndSkipSectionLine (
                             1,
                             stream, 
                            "Section:Smile", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    /* Smile number of entries */
    if (FindAndSkipComLine (stream, 
                            "Number of entries", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(mktvol_data->NbSmileDates)); 
    if (readerror != 1)
    {      
        DR_Error("Could not find Smile Nb entries in file %s! (Fix3_Model_Input_TimeDep)", FileName);
        goto RETURN;
    }

    if(mktvol_data->NbSmileDates < 1)
    {
        DR_Error("Number of entries in Smile section has to be at least 1 in file %s! (Fix3_Model_Input_TimeDep)", FileName);
        goto RETURN;
    }

    if (mktvol_data->NbSmileDates > MAXNBTD)
    {
         DR_Error("Maximum number of smile dependent input dates in file %s! surpassed"
                 "(Fix3_Model_Input_TimeDep)", 
                 FileName);
        goto RETURN;
    }

    /* Read smile parameters */
    if (FindAndSkipComLine (stream, 
                            "Date, QLeft, QRight, FwdShift", 
                            "Fix3_Model_Input_TimeDep", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (t = 0; t < mktvol_data->NbSmileDates; t++)
    {
        readerror = fscanf (stream, 
                            "%ld \t%lf \t%lf \t%lf \n", 
                            &(mktvol_data->SmileDate[t]), 
                            &(mktvol_data->QLeftTD[t]),
                            &(mktvol_data->QRightTD[t]),
                            &(mktvol_data->FwdShiftTD[t]));

        if (readerror != 4)
        {        
            DR_Error("Could not read smile info in file %s! (Fix3_Model_Input_TimeDep)", FileName);
            goto RETURN;
        }

        mktvol_data->QLeftTD[t]  = 1. - mktvol_data->QLeftTD[t];
        mktvol_data->QRightTD[t] = 1. - mktvol_data->QRightTD[t];

        if (Dateok(mktvol_data->SmileDate[t]))
        {
            DR_Error("Incorrect format for smile date in file %s (Fix3_Model_Input_TimeDep)!", FileName);
            goto RETURN;
        }
    }
   
    if (FindAndSkipSectionLine (
                             0, /* skip mode */
                             stream, 
                            "End Section", 
                            "Fix3_Model_Input_TimeDep", 
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

} /* Fix3_Model_Input_TimeDep */





/*****  Fix3_Param_Check_TimeDep  ********************************************************/
/*
*  	Read time dependent term structure input and check validity of input.
*/
int     Fix3_Param_Check_TimeDep(int         NbFactor, /* (I) Number of factors                     */
                     MKTVOL_DATA *mktvol_data,      /* (I) Structure of swaption volatility data */
                     FIX3_TREE_DATA   *tree_data)   /* (I) Tree data structure                   */
{
    int status = FAILURE;                      /* Error status = FAILURE initially */
       
    int     i, t, k;     
    int     NbTDInp;
    int     NbCorr;


    /* Nb of factors */
    if ((NbFactor != 1) && 
        (NbFactor != 2) &&
        (NbFactor != 3))
    {
        DR_Error("Nb of factors must be 1, 2 or 3!\n");
        goto RETURN;
    }

    /* Set total vol constants */
   
    if (mktvol_data->NbFactor != NbFactor || tree_data->NbFactor != NbFactor)
    {
        DR_Error("Inconsistent number of factors ! ");
        goto RETURN;
    }

    
    if (IS_EQUAL(mktvol_data->Bbq,1))
    {
        mktvol_data->VolNorm = 0.;
        mktvol_data->VolLogn = 1;
    }
    else
    if (IS_EQUAL(mktvol_data->Bbq,0))
    {
        mktvol_data->VolNorm = 1;
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

    for( i = 0; i < mktvol_data->NbSmileDates; i++)
    {
        if ((mktvol_data->QLeftTD[i]  < -10.) || (mktvol_data->QLeftTD[i]  > 10.) ||
            (mktvol_data->QRightTD[i] < -10.) || (mktvol_data->QRightTD[i] > 10.))
        {
            DR_Error("Q out of range!");
            goto RETURN;
        }

        if (IS_EQUAL(1. + mktvol_data->FwdShiftTD[i],0))
        {
            DR_Error("Fwd shift is singular!");
            goto RETURN;
        }
    }

    
    /* Nb of iterations */
    if (mktvol_data->CetNbIter > MAX_ITERATIONS)
    {
        DR_Error("Maximum allowed number of iterations is %d!",
                 MAX_ITERATIONS);
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
    if (mktvol_data->NbTDInp < 1)
    {
        DR_Error ("Number of time dependent inputs should be at least 1");
        goto RETURN;
    }
  
    /* check TD input date format */
    for (i = 0; i < NbTDInp; i++)
    {
        if (Dateok(mktvol_data->TDInpDate[i]))
        {
            DR_Error("Incorrect format for volatility date [%ld] (TDInput_Check_W)!", 
                mktvol_data->TDInpDate[i]);
            goto RETURN;
        }
    }

    /* check that dates are in ascending order */
    for (i = 1; i < NbTDInp; i++)
    {
        if (mktvol_data->TDInpDate[i] <= mktvol_data->TDInpDate[i-1])
        {
            DR_Error("TD Input dates must be entered in ascending order "
                        "(TDInput_Check_W)!");
            goto RETURN;
        }
    }


    for ( t = 0; t < NbTDInp; t++)
    {
        for (k = 0; k < NbFactor; k++)
		{
            /* check factor weights */
            if (mktvol_data->AlphaTD[k][t] < 0.0001)
            {
                DR_Error("VNMF Factor weight for time dependent entry nb %d, "
                         "factor nb %d is out of range (supplied value = %lf, "
                         "must be > 0.0001)", 
                         t, k, mktvol_data->AlphaTD[k][t]);
                goto RETURN;
            }

            if ((mktvol_data->BetaTD[k][t] < -10.) || (mktvol_data->BetaTD[k][t] > 10.))
            {
                DR_Error("VNFM mean reversion for time dependent entry nb %d "
                         "factor nb %d is out of range (supplied value = %lf, "
                         "must be between -10.0 and +10.0)",
                         t , k, mktvol_data->BetaTD[k][t]);
                goto RETURN;
            }
            
        }
        /* check correlation*/
        for (k = 0; k < NbCorr; k++)
		{
            if ((mktvol_data->RhoTD[k][t] < -.95) || (mktvol_data->RhoTD[k][t] > .95))
            {
                DR_Error("VNFM factor correlation for time dependent entry nb %d "
                         "column nb %d is out of range (supplied value = %lf,"
                         "must be between -0.95 and +0.95)",
                         t, k, mktvol_data->RhoTD[k][t]);
                goto RETURN;
            }
        }

    }

    /* verify TD input dates subset of benchmark swaps */
    k = 0;
    if ((mktvol_data->CalibFlag == TRUE) && (mktvol_data->NbTDInp > 1))
    {
        for (i = 0; i < mktvol_data->NbTDInp; i++)
        {
            while (( k < mktvol_data->NbVol) && (mktvol_data->TDInpDate[i]) != (mktvol_data->SwapSt[k]))
            {
                 k++;
            }
        
            if ( k == mktvol_data->NbVol)
            {
                DR_Error ("TD Inp Date no. %d (%ld) is not a subset of benchmark vol dates !", 
                          i, mktvol_data->TDInpDate[i]);
                goto RETURN;
            }
        }
    }

    /* verify smile dates are subset of benchmark swaps if calib flag is TRUE */
    k = 0;
    if ((mktvol_data->CalibFlag == TRUE) && (mktvol_data->NbSmileDates > 1))
    {
        for (i = 0; i < mktvol_data->NbSmileDates; i++)
        {
            while (( k < mktvol_data->NbVol) && (mktvol_data->SmileDate[i]) != (mktvol_data->SwapSt[k]))
            {
                 k++;
            }
    
            if ( k == mktvol_data->NbVol)
            {
                DR_Error ("Smile date no. %d (%ld) is not a subset of benchmark vol dates !", 
                          i, mktvol_data->SmileDate[i]);
                goto RETURN;
            }
        }
    }



  

    status = SUCCESS;
        
  RETURN:
        
    return (status);

}  /* Fix3_Param_Check_TimeDep */

