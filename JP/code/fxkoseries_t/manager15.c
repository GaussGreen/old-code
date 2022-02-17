/****************************************************************************/
/*      General I/O routine for FX option series data.                      */
/****************************************************************************/
/*      MANAGER15.C                                                         */
/****************************************************************************/

/*
$Header$
*/


#include     <stdio.h>
#include     <string.h>
#include     <stdlib.h>
#include     <math.h>
#include     "template15.h"





/*****  Fxkoseries_Manager  **************************************************/
/*
 *       Manage input & output, reading from ascii files into data structures.
 *
 */

int    Fxkoseries_Manager
    (T_CURVE          t_curve[2][3],      /* (O) Zero curves (C, I, R)       */
     MKTVOL_DATA      *mktvol_data,       /* (O) Base volatility data        */
     FX_DATA          *fx_data,           /* (O) FX data                     */
     FXKOSERIES_DATA  *fxkoseries_data,   /* (O) FX option series            */
     HYB3_TREE_DATA        *tree_data)         /* (O) Structure of tree data      */
{


    int
            i, k,
            readerror,                    /* Reading error status            */
            status = FAILURE;             /* Error status=FAILURE initially  */

    char
            FileName[MAXBUFF],
            ErrorMsg[MAXBUFF],
            *ErrorString[2];

    char   OwriteFxSpot[MAXBUFF];    /* O'writes related to FXVolatility.dat*/
    char   NbFXSmileOWS[MAXBUFF];    /* Number of entered FX Smile Parameters */
    char   FXSmileParamOWS[MAXNBDATE][MAXBUFF];  /* FX Smile Parameters */
    char   OwriteCorrel[3][MAXBUFF]; /* O'writes related to correlation.dat */
    char   OwriteMParamF[5][MAXBUFF];/* O'writes related to modelParams.dat */
    char   OwriteMParamD[5][MAXBUFF];/* O'writes related to modelParams.dat */

    FILE
            *stream = NULL;


    /* Assign strings for error messages as product contains for/dom data */
    ErrorString[0] = "foreign";
    ErrorString[1] = "domestic";



    /* Read deal file */
    strcpy (FileName, "fxkoseries_t.dat");
    stream = fopen (FileName, "r");  /* Open the deal data file (see fxkoseries_t.h) */
                
    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        return (status);              /* Don't proceed any further */
        
    }  /* if */ 
    

    /* Title line */
    if (FindAndSkipComLine (stream, "title line", 
                             "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    
    /* Long/short flag */
    if (FindAndSkipComLine(stream, "long/short flag", 
                           "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->LoS));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read long/short flag in %s!"
                           " (Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->LoS = (char)toupper(fxkoseries_data->LoS);                                               

    /* Knock out frequency specification */
    if (FindAndSkipComLine(stream, 
                           "ko frequency", 
                           "Fxkoseries_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->KoFreq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read ko frequency in file %s! "
                 "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->KoFreq = (char)toupper(fxkoseries_data->KoFreq);


    /* KO data */
    if (FindAndSkipComLine(stream, 
                           "number of ko dates", 
                           "Fxkoseries_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(fxkoseries_data->NbKoDates));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, 
                 "Could not read number of exercise dates in file %s! "
                 "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    if (FindAndSkipComLine (stream, 
                            "ko dates and barrier levels", 
                            "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    for (i = 0; i < fxkoseries_data->NbKoDates; i++)
    {
        readerror = fscanf(stream, "%ld %lf %lf %lf\n",
                           &(fxkoseries_data->KoDates[i]),
                           &(fxkoseries_data->LoBarrier[i]),
                           &(fxkoseries_data->HiBarrier[i]),
                           &(fxkoseries_data->Rebate[i]));
        if (readerror != 4)
        {        
            sprintf (ErrorMsg, "Could not read exercise data line "
                     "#%d in file %s! (Fxkoseries_Manager)", 
                     i+1, 
                     FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }	


    /* Knock in or knock out */
    if (FindAndSkipComLine(stream, 
                           "KI/KO flag", 
                           "Fxkoseries_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->KnockIoO));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read knock in/out flag in file %s! "
                 "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->KnockIoO = (char)toupper(fxkoseries_data->KnockIoO);


    /* Knock out inside or outside levels */
    if (FindAndSkipComLine(stream, 
                           "inside/outside KO flag", 
                           "Fxkoseries_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->IoO));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read inside/outside ko flag in file %s! "
                 "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->IoO = (char)toupper(fxkoseries_data->IoO);


   
    /* SmoothFlag flag */
    if (FindAndSkipComLine(stream, 
                           "smoothing flag", 
                           "Fxkoseries_Manager", 
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->SmoothFlag));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read smoothing flag in file %s! "
                 "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->SmoothFlag = (char)toupper(fxkoseries_data->SmoothFlag);


    /* Notional in domestic currency */
    if (FindAndSkipComLine (stream, "notional", "Fxkoseries_Manager", 
                            FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%lf \n", &(fxkoseries_data->Notional));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read notional in file %s! "
                           "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }


    /* Exercise frequency (E, N, M, Q, S, A) */
    if (FindAndSkipComLine(stream, "exercise freq", "Fxkoseries_Manager",
                           FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->ExerFreq));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,"Could not read exercise freq in %s! "
                          "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->ExerFreq = (char)toupper(fxkoseries_data->ExerFreq);


    /* Number of exercise dates */
    if (FindAndSkipComLine(stream, "number of exercise", 
                           "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(fxkoseries_data->NbExer));
    if (readerror != 1)
    {        
        sprintf (ErrorMsg,"Could not read number of exercise in file %s! "
                          "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Exercise dates and strikes */
    if (FindAndSkipComLine(stream, "exercise dates and strikes",
                           "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    for (i = 0; i < fxkoseries_data->NbExer; i++)
    {
        readerror = fscanf(stream, "%ld \t%lf \n",
                           &(fxkoseries_data->Exer[i]),
                           &(fxkoseries_data->Strike[i]));
        if (readerror != 2)
        {        
            sprintf (ErrorMsg,"Could not read exercise date and strike #%d "
                              "in file %s! (Fxkoseries_Manager)", i+1, FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }	

    /* Call or put */
    if (FindAndSkipComLine (stream, "call or put", "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%c \n", &(fxkoseries_data->CoP));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg,"Could not read call or put in file %s! "
                         "(Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }
    fxkoseries_data->CoP = (char)toupper(fxkoseries_data->CoP);
 

    /*     MODEL DATA                                                    */
    /* Index to be used for volatility calibration (foreign and domestic)*/
    /* The index names are placed directly into tree_data                */
    for (k=0; k<2; k++)
    {
        if (FindAndSkipComLine (stream, "volatility index", 
                                "Fxkoseries_Manager", FileName) == FAILURE)
        {        
            goto RETURN;
        }

        readerror = fscanf (stream, "%s \n", tree_data->Index[k]);	        
        if (readerror != 1)
        {        
            sprintf (ErrorMsg, "Could not read %s volatility index in "
                     "file %s! (Fxkoseries_Manager)",ErrorString[k],FileName);
            DR_Error (ErrorMsg);
            goto RETURN;
        }
    }

    /* Number of standard deviations at which to cut the tree */
    if (FindAndSkipComLine (stream, "number of standard deviations", 
                            "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }

    readerror = fscanf (stream, "%d \n", &(tree_data->NbSigmaMax));        
    if (readerror != 1)
    {        
        sprintf (ErrorMsg, "Could not read number of standard "
                 "deviations in file %s! (Fxkoseries_Manager)",FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /*  Number of periods per year in the tree */
    if (FindAndSkipComLine(stream, "periods per year", 
                           "Fxkoseries_Manager", FileName) == FAILURE)
    {        
        goto RETURN;
    }
    readerror = fscanf (stream, "%d \n", &(tree_data->Ppy));
    if (readerror != 1)
    {        
        sprintf(ErrorMsg, "Could not read periods per year in "
                "file %s! (Fxkoseries_Manager)", FileName);
        DR_Error (ErrorMsg);
        goto RETURN;
    }

    /* Finally read overwrite strings */

    /* This function reads all of the overwrite entries in the deal file
       referenced by the FILE stream pointer and populates the overwrite strings */
    if (Hyb3_ModelOverwrites_Input_W(stream,
                                FileName,   
                                OwriteFxSpot,
                                NbFXSmileOWS,
                                FXSmileParamOWS,
                                OwriteCorrel,
                                OwriteMParamF,
                                OwriteMParamD,
                                tree_data) == FAILURE)
    {
        goto RETURN;
    }


    /* ------------------------------------------------------------------*/
    /*    STANDARD ENVIRONMENT ELEMENTS (i.e. zero.dat, basevol.dat etc.)*/
    /*-------------------------------------------------------------------*/

    for (k=0; k<2; k++)
    {

        /* Hardcoded  assigment of zero  curves for the  engine */
        tree_data->CvDiff[k] = 0;  /* zero.dat is index curve   */
        tree_data->CvIdx1[k] = 1;
        tree_data->CvIdx2[k] = 2;
        tree_data->CvDisc[k] = 1;  

	MktVol_Init(&mktvol_data[k]);

        if (Term_Input_W (&(t_curve[k][0]),     /* Read index curve */
                          k ? "dzero.dat" 
                            : "fzero.dat") == FAILURE)
        {        
            goto RETURN;
        }
                    
        if (Term_Input_W (&(t_curve[k][1]),     /* Read discount curve */
                          k ? "ddisczero.dat"
                          :   "fdisczero.dat") == FAILURE)
        {        
            goto RETURN;
        }
                    
        if (Term_Input_W (&(t_curve[k][2]),    /* Read credit curve:       */
                          k ? "ddisczero.dat"  /* same as discount for now */
                          :   "fdisczero.dat") == FAILURE)
        {        
            goto RETURN;
        }
                    
        if (Hyb3_Param_Input (&(mktvol_data[k]),  /* Read  mean rev parameter */
                         tree_data->Index[k],     /* Calib Index can be overwritten in env */
                         tree_data,
                         1,                      /* Always one-factor for IR */
                         k ? OwriteMParamD
                           : OwriteMParamF,
                         k ? "dmodelParameters.dat"
                           : "fmodelParameters.dat") == FAILURE)
        {        
            goto RETURN;
        }
        
        if (MktVol_Input_W (&(mktvol_data[k]),   /* Read vol data */
                            tree_data->Index[k],
                            &(t_curve[k][tree_data->CvDiff[k]]),
                            k ? "dbasevol.dat"
                              : "fbasevol.dat",
                            k ? "dswapvol.dat" 
                             : "fswapvol.dat"  ) == FAILURE)
        {        
            goto RETURN;
        }


    } /* for k */


    /*  FX data  */
    if (Hyb3_Fx_Input_W_WithSmile(fx_data,
                   OwriteFxSpot,
                   "FXVolatility.dat",
                   "fxsmile_0.dat",
                    NbFXSmileOWS,
                    FXSmileParamOWS) ==  FAILURE)
    {
        goto RETURN;      
    }  


    /* Read in DR WRapper Type3 correlation.dat file and */
    /* process the overwrite strings at the same time    */
    if(Hyb3_Correl_Input_WType3(fx_data,
                           OwriteCorrel, /* (O) O'write strings*/
                           "correlation.dat") == FAILURE)
    {
        goto RETURN;
    }
                      
        
    if (Fxkoseries_Check(fxkoseries_data,
                       fx_data,
                       tree_data) == FAILURE)
    {
        goto RETURN;      
    } 




    status = SUCCESS;
        
    RETURN:

    fclose (stream);

    return (status);



}  /* End of Fxkoseries_Manager */



/*****  Fxkoseries_Check  ****************************************************/
/*
*       Check the inputs.
*/
int     Fxkoseries_Check 
            (FXKOSERIES_DATA  *fxkoseries_data,
             FX_DATA          *fx_data,      /* (I)                          */
             HYB3_TREE_DATA        *tree_data)    /* (I) Structure of tree data   */
{

    char
        ErrorMsg[MAXBUFF],
        *ErrorString[2];

    int
        i,
        k,
        status = FAILURE;        /* Error status = FAILURE initially */

 
    ErrorString[0] = "foreign";
    ErrorString[1] = "domestic";


    /* Check all flags */
    if ((fxkoseries_data->LoS != 'L') && 
        (fxkoseries_data->LoS != 'S') ) 
    {
        DR_Error("Long/short flag must be L or S.\n");
        goto RETURN;
    }

    if ((fxkoseries_data->KnockIoO != 'I') && 
        (fxkoseries_data->KnockIoO != 'O') ) 
    {
        DR_Error("Knock in/knock out flag must be I or O.\n");
        goto RETURN;
    }

    if ((fxkoseries_data->IoO != 'I') && 
        (fxkoseries_data->IoO != 'O') ) 
    {
        DR_Error("Inside/outside barrier range flag must be I or O.\n");
        goto RETURN;
    }

    if ((fxkoseries_data->SmoothFlag != 'Y') && 
        (fxkoseries_data->SmoothFlag != 'N') )   
    {
        DR_Error("Invalid entry for smoothing flag ('Y' or 'N')!");
        goto RETURN;
    } 

    if ((fxkoseries_data->CoP != 'C') && 
        (fxkoseries_data->CoP != 'P') ) 
    {
        DR_Error("Call or put must be C or P.\n");
        goto RETURN;
    } 

    if ((fxkoseries_data->ExerFreq != 'E') && /* European  */
        (fxkoseries_data->ExerFreq != 'M') && /* Monthly   */
        (fxkoseries_data->ExerFreq != 'Q') && /* Quarterly */
        (fxkoseries_data->ExerFreq != 'S') && /* Semi-a    */
        (fxkoseries_data->ExerFreq != 'A') )  /* Annual    */
    {
        DR_Error("Exercise frequency must be I, M, Q, S or A.\n");
        goto RETURN;
    }

    if ((fxkoseries_data->KoFreq != 'I') && /* Input     */
        (fxkoseries_data->KoFreq != 'N') && /* americaN  */
        (fxkoseries_data->KoFreq != 'D') && /* Daily     */
        (fxkoseries_data->KoFreq != 'W') && /* Weekly    */
        (fxkoseries_data->KoFreq != 'M') && /* Monthly   */
        (fxkoseries_data->KoFreq != 'Q') && /* Quarterly */ 
        (fxkoseries_data->KoFreq != 'S') && /* Semi-a    */
        (fxkoseries_data->KoFreq != 'A') )  /* Annual    */ 
    {
        DR_Error ("Invalid entry for monitoring style!");
        goto RETURN;
    } 

    /* Nb of ko dates can only be 1 for European monitoring */
    if (fxkoseries_data->NbKoDates < 1)
    {
        DR_Error ("Specify at least one ko date !");
        goto RETURN;
        
    }  /* if */      
    if (   (fxkoseries_data->NbKoDates == 1) 
        && (fxkoseries_data->KoFreq != 'I')   )            
    {
        DR_Error ("A minimum of two ko dates is required, "
                 "except for European ko monitoring!");
        goto RETURN;
    }   

    /* Ko date format and ascending order */
    for (i = 0; i < fxkoseries_data->NbKoDates; i++)
        if (Dateok(fxkoseries_data->KoDates[i]))
        {
            DR_Error ("Incorrect format for ko date!");
            goto RETURN;
        
        }  /* if */                                                

    for (i = 1; i < fxkoseries_data->NbKoDates; i++)
        if (fxkoseries_data->KoDates[i] <= fxkoseries_data->KoDates[i-1])
        {
            DR_Error ("Ko dates must be entered in ascending order!");
            goto RETURN;
        
        }  /* if */                                                

    /* At least one ko date after value date */    
    if (fxkoseries_data->KoDates[fxkoseries_data->NbKoDates-1] < fx_data->ValueDate)
    {
        DR_Error ("Specify at least one ko date after value date!");
        goto RETURN;
        
    }  /* if */      
        
    /* Low barrier smaller than high barrier */
    for (i = 0; i < fxkoseries_data->NbKoDates; i++)
        if (fxkoseries_data->LoBarrier[i] > fxkoseries_data->HiBarrier[i])
        {
            DR_Error ("Low barrier has to be smaller than high barrier !");
            goto RETURN;
        
        }  /* if */                      
                          
    /* Ko dates must match its frequency */
    if ((fxkoseries_data->KoFreq != 'I') && (fxkoseries_data->KoFreq != 'N'))
    {
        if (DrDatesInSchedule(fxkoseries_data->NbKoDates,
                              fxkoseries_data->KoDates,
                              fxkoseries_data->KoDates[0],
                              fxkoseries_data->KoDates[fxkoseries_data->NbKoDates-1],
                              fxkoseries_data->KoFreq,
                              'N') == FAILURE)
        {
            DR_Error ("Ko date does not agree with monitoring frequency!");
            goto RETURN;
        }
    }

    /* No rebate for knock-in option */
    if (fxkoseries_data->KnockIoO == 'I')
        for (i = 0; i < fxkoseries_data->NbKoDates; i++)
            if (fabs(fxkoseries_data->Rebate[i]) > ERROR)
            {
                DR_Error ("Rebate feature is not supported for knock-in options!");
                goto RETURN;
        
            }  /* if */                      
                          
    /* Check notional */
    if (fxkoseries_data->Notional < 0.0)
    {
        DR_Error("Notional must be positive (use long or short flag).\n");
        goto RETURN;
    }


    /* At least one exercise date if European, two otherwise */
    if (fxkoseries_data->ExerFreq == 'E')
    {
        if (fxkoseries_data->NbExer < 1)
        {
            DR_Error("There must be at least one exercise date.\n");
            goto RETURN;
        }
    }
    else
    {
        if (fxkoseries_data->NbExer < 2)
        {
            DR_Error("When an exercise frequency is specified, there must be at least\n"
                    "two exercise dates.\n");
            goto RETURN;
        }
    }
                       
    /* Exercise date format and ascending order */
    for (i = 0; i < fxkoseries_data->NbExer; i++)
        if (Dateok(fxkoseries_data->Exer[i]))
        {
            DR_Error ("Incorrect format for exercise date!");
            goto RETURN;
        
        }  /* if */                                                

    for (i=1; i<fxkoseries_data->NbExer; i++)
    {
        if (fxkoseries_data->Exer[i] <= fxkoseries_data->Exer[i-1])
        {
            sprintf(ErrorMsg,"Exercise date #%d <= previous date.\n", i);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
    }


    /* At least one of the exercise dates must be > value date */
    if (fxkoseries_data->Exer[fxkoseries_data->NbExer-1] < fx_data->ValueDate)
    {
        sprintf(ErrorMsg,"At least on exercise date must be on or "
                "after fx value date.\n");
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    /* Ko window within exercise window */
    if (fxkoseries_data->KoDates[fxkoseries_data->NbKoDates-1] > fxkoseries_data->Exer[fxkoseries_data->NbExer-1])
    {
        DR_Error ("Ko window must end before exercise window !");
        goto RETURN;
        
    }  /* if */      
        
    /* Unless exercise on given dates only, dates must match frequency */
    if (fxkoseries_data->ExerFreq != 'E')
    {
        if (DrDatesInSchedule(fxkoseries_data->NbExer,
                              fxkoseries_data->Exer,
                              fxkoseries_data->Exer[0],
                              fxkoseries_data->Exer[fxkoseries_data->NbExer-1],
                              fxkoseries_data->ExerFreq,
                              'N') == FAILURE)
        {
            DR_Error("Specified exercise dates are not subset of dates\n"
                    "obtained from exercise frequency.\n");
            goto RETURN;
        }
    }


    /* Checks over foreign/domestic data */
    for (k=0; k<2; k++)
    {
               
        
        /* Diffused curve number */
        if ((tree_data->CvDiff[k] != 0) &&
            (tree_data->CvDiff[k] != 1) &&
            (tree_data->CvDiff[k] != 2) )
        {
            DR_Error("Diffused curve number must be 0, 1, 2 !");
            goto RETURN;
        } 
        
        /* Index curve number */
        if ((tree_data->CvIdx1[k] != 0) &&
            (tree_data->CvIdx1[k] != 1) &&
            (tree_data->CvIdx1[k] != 2) )
        {
            DR_Error("Index1 curve number must be 0, 1, 2 !");
            goto RETURN;
        } 
        if ((tree_data->CvIdx2[k] != 0) &&
            (tree_data->CvIdx2[k] != 1) &&
            (tree_data->CvIdx2[k] != 2) )
        {
            DR_Error("Index2 curve number must be 0, 1, 2 !");
            goto RETURN;
        } 
        
        /* All three curves have different indexes in tree */
        if ((tree_data->CvDiff[k] == tree_data->CvIdx1[k]) || 
            (tree_data->CvDiff[k] == tree_data->CvIdx2[k]) || 
            (tree_data->CvIdx1[k] == tree_data->CvIdx2[k]))
        {
            DR_Error("Curve indexing is wrong in the tree !\n"
                     "The same curve index cannot be used more\n"
                     "than once (3 different curves are modelled).\n");
            goto RETURN;
        } 
        
        /* Discount curve number */
        if ((tree_data->CvDisc[k] != 0) &&
            (tree_data->CvDisc[k] != 1) &&
            (tree_data->CvDisc[k] != 2) )
        {
            DR_Error("Discount curve number must be 0, 1, 2 !");
            goto RETURN;
        } 

        
            

    } /* for k (=0 for foreign, =1 for domestic) */


    /* Model parameters entered through deal file */
    if (tree_data->NbSigmaMax < 3)
    {
        DR_Error("Can't cut the tree at less than three std devs!");
        goto RETURN;
     
    }  /* if */                                                


    if (tree_data->Ppy < 1)
    {
        DR_Error("Can't build a tree with less than 1 period per year!");
        goto RETURN;
        
    }  /* if */  


    if (tree_data->Ppy > 365)
    {
        DR_Error("Can't build a tree with more than 365 periods per year!");
        goto RETURN;
        
    }  /* if */  


    status = SUCCESS;

    RETURN:

    return (status);

}  /* End of Fxseries_Check */





/*****  Print_Fxkoseries  ****************************************************/
/*
*       Print debug information in an ascii file.
*/
int     Print_Fxkoseries
            (T_CURVE       t_curve[2][3],/* (I) Structure of zero curve data */
             HYB3_TREE_DATA     *tree_data)   /* (I) Structure of tree data       */
{
    int
            i, k, l,
            status = FAILURE; /* Error status = FAILURE initially        */

    double          
            days,      /* Number of days from today to the current node  */
            Forward,   /* Forward rate for the current period            */
            discount;

    char
            *ErrorString[2];

    FILE
            *stream;


    ErrorString[0] = "Foreign";
    ErrorString[1] = "Domestic";



    stream = fopen ("TERM.prn", "w");

    if (stream == NULL)
    {
        DR_Error ("Could not open TERM.prn! (Print_Fxseries)");
        return (status);              /* Don't proceed any further */

    }  /* if */


    /* Schedule of this particular product */
    fprintf (stream, "\nPRODUCT SCHEDULE:\n");
    fprintf (stream,"\nNode   Date        Ko             Low            High          "
                    "Rebate    Exercise      Strike \n");

    for (i = 0; i <= tree_data->NbTP; i++)
    {
        fprintf (stream,
                "[%3d]%10ld    (%d)   %14.4f %14.4f %14.4f    (%d)   %14.4f \n",
                i,
                tree_data->TPDate[i],
                tree_data->TPType[0][i],
                tree_data->CritDate[0][i].Value[0],
                tree_data->CritDate[0][i].Value[1],
                tree_data->CritDate[0][i].Value[2],
                tree_data->TPType[1][i],
                tree_data->CritDate[1][i].Value[0]);
        
    }  /* for i */


    fprintf (stream, "\n\n");



    /* Slice sizes used for allocation */
    fprintf(stream, "   W1   HW1    W2   HW2    W3   HW3\n");
    fprintf(stream, "%5d %5d %5d %5d %5d %5d\n\n\n",
                    tree_data->Width[0], tree_data->HalfWidth[0],
                    tree_data->Width[1], tree_data->HalfWidth[1],
                    tree_data->Width[2], tree_data->HalfWidth[2]);




    for (l=0; l<2; l++)
    {
        fprintf (stream,"\nTIMELINE INFORMATION (%s IR)\n",ErrorString[l]);

        fprintf (stream,
                 "Node      Date     Days  Max  Forward   "
                 "Zero0   Discount0   Zero1   Discount1   "
                 "Zero2   Discount2    IrMidNode    SpotVol \n");
        
        for (i = 0; i <= tree_data->NbTP; i++)
        {
                
            /* Conversion: 1.+FwdRate[i]=(1.+Forward)^(Length[i]/365.) */
            /*Forward = 100.*(pow (1.+tree_data->FwdRate[l][0][i],   
                                1./tree_data->Length[i])-1.);*/
            Forward = 100. * tree_data->FwdRate[l][0][i]/tree_data->Length[i];

            days = Daysact (tree_data->TPDate[0], tree_data->TPDate[i]);

            fprintf(stream,
                    "[%3d] \t%8ld  %5.0f  %3d  %7.4f  %7.4f   %8.6f  "
                    "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                    i,
                    tree_data->TPDate[i],
                    days,
                    l ? tree_data->Top2[i][0] : tree_data->Top1[i],
                    Forward,
                    tree_data->ZeroRate[l][0][i] * 100.,
                    tree_data->ZeroCoupon[l][0][i],
                    tree_data->ZeroRate[l][1][i] * 100.,
                    tree_data->ZeroCoupon[l][1][i],
                    tree_data->ZeroRate[l][2][i] * 100.,
                    tree_data->ZeroCoupon[l][2][i],
                    exp(tree_data->IrZCenter[l][i]),
                    tree_data->SpotVol[l][i] * 100.);

        }  /* for i */
    } /* For l */



    fprintf (stream, "\n\n");

    fprintf (stream,"\nTIMELINE INFORMATION (FX)\n");
        
    fprintf (stream, "Node \t  Date \t\tCWidth   Forward    FxVol   "
                     " Spotvol   Rho's (Irf-Ird, Irf-FX, Ird-FX)\n");
                                        
    for (i = 0; i <= tree_data->NbTP; i++)
    {
        fprintf (stream, "[%3d]  %8ld  %6d  %10.4f    "
                         "%5.2f    %5.2f       %4.2f       %4.2f      %4.2f\n",
                         i, 
                         tree_data->TPDate[i],
                         tree_data->Top3[i][0][0],
                         tree_data->FwdFx[i],
                         tree_data->FxVol[i] * 100,
                         tree_data->SpotFxVol[i] * 100.,
                         tree_data->Rho[0][i],
                         tree_data->Rho[1][i],
                         tree_data->Rho[2][i]);

    }  /* for i */


    /* Print out input zero curves */
    for (l=0; l<2; l++)
    {
        fprintf (stream, "\n\n");
        fprintf (stream, "%s Currency Curves (Index, COF, Risk)\n\n",ErrorString[l]);
        for (k = 0; k < 3; k++)
        {
                
            fprintf (stream, "Maturity      Zero       Discount \n");

            for (i = 0; i < (t_curve[l][k]).NbZero; i++)
            {

                days = Daysact ((t_curve[l][k]).ValueDate, 
                                (t_curve[l][k]).ZeroDate[i]);
                /* Discount factor up to the current date */
                discount = pow(1. + (t_curve[l][k]).Zero[i], -days/365.);          

                fprintf (stream,
                         "%ld   %9.6f     %8.6f \n",
                         (t_curve[l][k]).ZeroDate[i],
                         (t_curve[l][k]).Zero[i] * 100.,                            
                         discount);
            }  /* for i */

            fprintf (stream, "\n");

        }  /* for k */

    } /* For l */

    fclose (stream);
        
        
    status = SUCCESS;

    return (status);

}  /* End of Print_Fxkoseries */
