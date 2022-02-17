/****************************************************************************/
/*      Calibration Enhancement Tool                                        */
/****************************************************************************/
/*      CET.C                                                               */
/****************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fix123head.h"
#include "q3_tmx.h"





#ifndef  CET_ERROR_TOLERANCE
#define  CET_ERROR_TOLERANCE    0.50           /* in vegas     */
#endif

#ifndef  CET_WARNING_TOLERANCE
#define  CET_WARNING_TOLERANCE  0.05           /* in vegas     */
#endif

/* Delta points for Smile definition & calibration */
#define  MIDSTRIKE         2
extern const double DSTRIKE[];



/*****  Fix3_Cet_Tmx  ********************************************************/
/*
* Determine the Libor SV-MQ distribution needed so that the tree-based Yield
* - calculated via Libor - fits the input smile.
*
*/
int Fix3_Cet_Tmx (
              int            CetOutputFlag,    /* (I) On TRUE write to file */
              T_CURVE const  *t_curve,         /* (I) Zero curve            */
              MKTVOL_DATA    *mvd,             /* (I/O) Volatility data     */
              FIX3_TREE_DATA const* tree_dataProd)   /* (I) Tree data       */
{


    FIX3_TREE_DATA tree_data;          /* Local Tree                        */
    long           ValueDate;

    /* Indices                 */
    int            i, NbCetIter, IterIdx,CvDiff;
    char           SaveTrace;

    /* Local array for printing    */
    double         *VolPrev = NULL,
                    VolDiff;

    int            isFirstWarning = TRUE;  /* Prints warning header info    */
    int            status = FAILURE;

    
    
    /* Initialisations         */
    ValueDate   = t_curve->ValueDate;
    CvDiff      = tree_dataProd->CvDiff;
    NbCetIter   = mvd->CetNbIter;
    SaveTrace   = mvd->TraceFlag;
    mvd->TraceFlag = 'N';   /* No printing to NMR.prn during Cet */

    VolPrev     = (double *) DR_Array(DOUBLE,0,mvd->NbVol);

    if (VolPrev == NULL)
    {
        DR_Error("Cet_Main: Could not allocate memory.\n");
        goto FREE_MEM_AND_RETURN;
    }

    /* Initialize tree */
    Fix3_Tree_Init (&tree_data);

    /*  Determine the Target Swap Atm at each Vol Bmk & the Target Swap    */
    /*  Sml at each Sml Bmk.                                               */
    if (Cet_NmrVanl(ValueDate,
                    t_curve,             
                    mvd,
                    CvDiff) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    /* Initialise Local Tree                                               */
    if (Cet_Manager(tree_dataProd,
                    &tree_data) == FAILURE)
    {                                      
        goto FREE_MEM_AND_RETURN;
    }

    /* Perform the Cet iteration                                           */
    for (IterIdx = 1; IterIdx <= NbCetIter; IterIdx++)
    {
        /* Build CET Time Line for Local Tree (use Nmr Dates)              */
        if (Cet_Schedule(ValueDate,
                         mvd,                              
                         &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;   
        }

        /* Build Local Tree (inc. Nmr Calc., i.e. Bmk Pricing on the Tree)  */
        if (Fix3_Build_Tree (t_curve,
                        mvd,
                        &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;   
        }


        /* Calculate the Smile for each Bmk Swap                            */
        if (Cet_Calc (t_curve,
                      mvd,
                      &tree_data) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }
        
        /* If it has not converged, update the Atm Vol & Smile at Vol/Sml   */
        /* Bmk Dates & then interpolate to all Nmr Dates                    */
        /* If it has not converged, update the Atm Vol & Smile at Vol/Sml   */
        /* Bmk Dates & then interpolate to all Nmr Dates                    */
        for (i = 0; i < mvd->NbVol; i++)
        {
            VolPrev[i] = mvd->Vol[i];
        }

        if(IterIdx < NbCetIter)
        {
            if (Cet_UpdSwapSml(ValueDate,
                               mvd) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }

        /* Print current Cet iteration to CET.prn & SMILE.prn               */
        if (CetOutputFlag)
        {
            if (Cet_Print (IterIdx,
                           VolPrev,
                           mvd) == FAILURE) 
            {
                goto FREE_MEM_AND_RETURN;
            }
        }
        {
            int NbFactor = tree_data.NbFactor;
            int NbTP = tree_data.NbTP;
            int NbEDevDates = tree_data.NbEDevDates;
            int JumpPpy = tree_data.JumpPpy;
            int NbDailyPts = tree_data.NbDailyPts;
            int NbNmr = tree_data.NbNmr;
            Fix3_Tree_Free(&tree_data);
            tree_data.NbFactor = NbFactor;
            tree_data.NbTP = NbTP;
            tree_data.NbEDevDates = NbEDevDates;
            tree_data.JumpPpy = JumpPpy;
            tree_data.NbDailyPts = NbDailyPts;
            tree_data.NbNmr = NbNmr;
        }
        /* re-initialize pointers back to null */
        Fix3_Cet_Tree_Reset(&tree_data);

    } /* for IterIdx */

    /* Convergence Report:Check final CET vols are tolerable & report       */
    /* error/warning as appropriate. Do this AFTER call to Cet_Print for    */
    /* accesss to full debug info via CET.prn                               */ 
    for (i = 0; i < mvd->NbVol; i++)
    {
        VolDiff = mvd->SwapVnVol[MIDSTRIKE][i] - mvd->SwapTreeVol[MIDSTRIKE][i];
        
        if (fabs(VolDiff) > CET_ERROR_TOLERANCE)
        {
            /* print error to stdout and exit */
            DR_Error("Cet_Main Error: Final CET vol has not converged to within %4.2lf vega!\n"
                    "Expiry = %ld SwapSt = %ld SwapMat = %ld VegaDiff = %4.2lf vega\n",
                     CET_ERROR_TOLERANCE, 
                     mvd->SwapSt[i],
                     mvd->SwapSt[i], 
                     mvd->SwapMat[i],
                     VolDiff);
            
            goto FREE_MEM_AND_RETURN;
        }
        else
        {
            if (fabs(VolDiff) > CET_WARNING_TOLERANCE)
            {
                if (isFirstWarning)
                {
                    /* print header info */
                    DR_Warning(TRUE, /* open clean file */
                            "Cet_Main Warning: Final CET vol has not converged to within %6.4lf vega!\n"
                            "Expiry  \t SwapSt  \t SwapMat \t VegaDiff", 
                            CET_WARNING_TOLERANCE);
                    
                    isFirstWarning = FALSE;
                }
                
                /* print warning to stdout and warning.log but do not exit */
                DR_Warning(FALSE, /* append to file */
                        "%ld\t %ld\t %ld\t %6.4lf", 
                        mvd->SwapSt[i],
                        mvd->SwapSt[i], 
                        mvd->SwapMat[i],
                        VolDiff);
            }
        }

    } /* for i < NbVol */

    /* Restore Trace Flag before exit */
    mvd->TraceFlag = SaveTrace;

    status = SUCCESS;

FREE_MEM_AND_RETURN:

    Free_DR_Array (VolPrev, DOUBLE, 0, mvd->NbVol);
    
    Fix3_Tree_Free (&tree_data);
                                               
    return (status);

}  /* Fix3_Cet_Tmx */




/*****  Cet_Calc  **********************************************************/
/*
*      Price the Smile for each Bmk Swap
*/
int Cet_Calc (T_CURVE const        *t_curve,       /* (I) Zero curve       */
              MKTVOL_DATA          *mvd,           /* (I) Volatility data  */
              FIX3_TREE_DATA       *tree_data)     /* (O) Tree data        */
{
    FIX3_DEV_DATA    dev_data;         /* Dev data structure                    */

    /* Slices */
    double      *Annuity[MAX_INST]    ; /* Annuity slice                      */
    double      *LastZero[MAX_INST]   ; /* Longest maturity zero              */
    double      *BMarkPrice[MAX_INST] ; /* Benchmark Swaption Price           */
    double      *Yield       = NULL   ; /* Swap Yield                         */       

    long        ExerFlag[MAX_INST];     /* Indicates an exercise date         */
    long        CpnFlag[MAX_INST];      /* Indicates a fixed payment date     */
    long        MatFlag[MAX_INST];      /* Indicates a maturity date          */
    long        NmrFlag;
    double      Dccfrac[MAX_INST];
    
    /* Numeraire date info */
    int         NmrIdx = 0;
    int         Smoothing;              /* Use smoothing if TRUE               */


    /* Numerical variables */
    long        CurrentDate;
    int         CvDiff;           /* Diffuse curve number                  */
    int         i;                /* Node indices                          */
    int         t;                /* Current time point                    */   
    int         T;                /* Last time point                       */
    int         status = FAILURE; /* Error status = FAILURE initially      */
    
    /* Dates and indices */
    T            = tree_data->NbTP;
    CvDiff       = tree_data->CvDiff;

     /* Smoothing flag */
    Smoothing = (mvd->SmoothingFlag == 'Y');

    /* Initialize and allocate DEV_DATA structure */
    Fix3_Dev_Init (&dev_data);
    if (Fix3_Dev_Alloc (&dev_data, tree_data) == FAILURE) goto RETURN;

    for (i = 0; i < MAX_INST; i++)
    {
        Annuity[i]  = NULL;
        LastZero[i] = NULL;
    }

    /* Numeraire must be interpolated */
    tree_data->NmrInterpOn = TRUE;


    /* Allocate variables for yield calculation */
    for (i = 0; i < mvd->NbVol; i++)
    {
        /* Smile liquid date */
        if (mvd->SmlLiqDate[i])
        {
            Annuity[i]  = Fix3_Alloc_Slice (tree_data);
            LastZero[i] = Fix3_Alloc_Slice (tree_data);
            
            if ( (Annuity[i]  == NULL) ||
                 (LastZero[i] == NULL) )    
            {
                DR_Error("Cet_Calc: Could not allocate memory.\n");
                goto RETURN;
            }
        }
        else  /* Benchmark vol date */
        {
            BMarkPrice[i] = Fix3_Alloc_Slice (tree_data);

            if (BMarkPrice[i] == NULL)
            {
                DR_Error("Cet_Calc: Could not allocate memory.\n");
                goto RETURN;
            }
        }
        
    }

    /* Allocate memory - Extra slices */
    Yield  = Fix3_Alloc_Slice(tree_data);
    
    if (Yield    == NULL)
    {
        DR_Error("Cet_Calc: could not allocate memory for b'mark swaptions!");
        goto RETURN;        
    }

    /* Step backward through the tree */
    for (t = T; t >= 0; t--)
    {
        /* Current Date */
        CurrentDate = tree_data->TPDate[t];
        
        /* Set Flags for Coupon Pmts, Maturities and Exercises */
        NmrFlag = FALSE;

        for (i = 0; i < mvd->NbVol; i++)
        {
            CpnFlag[i] = tree_data->TPtype[i][t]; 

            if (CpnFlag[i])
            {
                Dccfrac[i] = (tree_data->CritDate[i][t]).Value[0]; 
            }
        
            ExerFlag[i] = FALSE;
            
            if (CurrentDate == mvd->SwapSt[i])
            {
                ExerFlag[i] = TRUE;
            }

            MatFlag[i] = FALSE;
            
            if (CurrentDate == mvd->SwapMat[i])
            {
                MatFlag[i] = TRUE;
            }

            NmrFlag = NmrFlag || CpnFlag[i] || ExerFlag[i] || MatFlag[i];            
        }

        dev_data.NmrToCcy = NmrFlag;
         
        /* 'Update' tree */
        if (Fix3_Lattice_Tmx (
                     &dev_data,
                     t,
                     T,
                     mvd,
                     tree_data) == FAILURE) goto RETURN;


        for (i = 0; i < mvd->NbVol; i++)
        {
            /* If necessary, Ev Annuity & Last Zero */
            if ( (CurrentDate >= mvd->SwapSt[i]) && (CurrentDate <= mvd->SwapMat[i]))
            {
                /* Smile Liquid expiry */
                if (mvd->SmlLiqDate[i])
                {
                    if (Fix3_Dev_Tmx (
                             Annuity[i],
                             t,
                             T,
                             CvDiff,
                             &dev_data,
                             tree_data) == FAILURE) goto RETURN;
                    
                    if (Fix3_Dev_Tmx (
                             LastZero[i],
                             t,
                             T,
                             CvDiff,
                             &dev_data,
                             tree_data) == FAILURE) goto RETURN;
                }

                /* Vol Benchmark */
                else
                {
                    if (Fix3_Dev_Tmx (
                             BMarkPrice[i],
                             t,
                             T,
                             CvDiff,
                             &dev_data,
                             tree_data) == FAILURE) goto RETURN;
                }
            }
            
            /* At coupon dates, add a new zero to the annuity */
            if (CpnFlag[i])
            {
                /* Smile liquid expiry */
                if (mvd->SmlLiqDate[i])
                {
                    if (Fix3_AddScalar (Annuity[i],
                                        Dccfrac[i],
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                }

                 /* Vol Benchmark */
                else
                {
                    if (Fix3_AddScalar (BMarkPrice[i],
                                        Dccfrac[i] * mvd->ParYield0[i],
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
            }

            /* At Swap Maturity, initialise Last Zero to 1 */
            if (MatFlag[i])
            {
                /* Smile liquid expiry */
                if (mvd->SmlLiqDate[i])
                {
                    if (Fix3_Set_Slice (LastZero[i],
                                        1, 
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
                else
                {
                    if (Fix3_AddScalar (BMarkPrice[i],
                                        1.0,
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
            }

            /* At Exercise, compute the swap yield & determine the swap vols */
            if (ExerFlag[i])
            {
                /* Find corresponding numeraire dates */
                NmrIdx = (int)(tree_data->CritDate[NMREVENT][t]).Value[0];
                
                if (mvd->NmrDate[NmrIdx] != mvd->SwapSt[i])
                {
                    DR_Error("Cet_Calc : inconsistency in the numeraire timeline !");
                    goto RETURN;
                }

                /* Smile liquid expiry */
                if (mvd->SmlLiqDate[i])
                {
                    
                    /* Set Swap Yield Slice */
                    if (Fix3_Set_Slice (Yield,
                                        1.,
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                    
                    /* Add Last Zero */
                    if (Fix3_LCombTwoSlices (Yield,
                                             Yield,
                                             1.0,
                                             LastZero[i],
                                             -1.0,
                                             t,
                                             tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                    
                    /* Divide by annuity */
                    if (Fix3_DivideTwoSlices (Yield,
                                              Yield,
                                              Annuity[i],
                                              t,
                                              tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }
                    
                    
                    /* Deduce Swaption Implied Vols */
                    if (Cet_SwapVols(i,  /* Benchmark instrument index */
                                     Smoothing, 
                                     Annuity[i],
                                     dev_data.NmrInv[CvDiff],
                                     Yield,
                                     NULL, 
                                     mvd, 
                                     tree_data,  
                                     &dev_data,
                                     t) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
                
                /* Vol benchmark expiry */
                else
                {
                    /* extra extra coupon (funding) */
                    if (Fix3_AddScalar (BMarkPrice[i],
                                        -1.0,
                                        t,
                                        tree_data) == FAILURE)
                    {
                        goto RETURN;
                    }

                    /* Deduce Swaption Implied Vols */
                    if (Cet_SwapVols(i,  /* Benchmark instrument index */
                                     Smoothing, 
                                     NULL,
                                     dev_data.NmrInv[CvDiff],
                                     NULL,
                                     BMarkPrice[i], 
                                     mvd, 
                                     tree_data,  
                                     &dev_data,
                                    t) == FAILURE)
                    {
                        goto RETURN;
                    }
                }
            } /* if ExerFlag */

        } /* for i */


        /* Record conversion to Ccy so that slices can be converted 
         * back to numeraire at next time point */
         dev_data.CcyToNmr = dev_data.NmrToCcy;

    } /* for t */

    /* From now on it is possible to intepolate numeraire */
  
    status = SUCCESS;

    RETURN:

    for (i = 0; i < mvd->NbVol; i++)
    {
        if (mvd->SmlLiqDate[i])
        {
            Fix3_Free_Slice(Annuity[i], tree_data);
            Fix3_Free_Slice(LastZero[i],tree_data);
        }
        else
        {
            Fix3_Free_Slice(BMarkPrice[i], tree_data);
        }
    }    

    Fix3_Free_Slice (Yield,     tree_data); 
    Fix3_Dev_Free   (&dev_data, tree_data);

    if (status != SUCCESS)
    {        
        DR_Error("Cet_Calc: Failed!");
    }

    return (status);

}  /* Cet_Calc */



/*****  Cet_Manager  ********************************************************/
/*
*   Initialise Local Tree
*/
int     Cet_Manager 
            (FIX3_TREE_DATA const*tree_dataProd,    /* (I) Tree data        */
             FIX3_TREE_DATA      *tree_data)        /* (O) Tree data        */
{
    int     status = FAILURE;               /* Status = FAILURE initially   */

    /*  Model specific data  */
    tree_data->NbSigmaMax = tree_dataProd->NbSigmaMax;
    tree_data->NbFactor   = tree_dataProd->NbFactor;

    /*  Standard environment elements */
    tree_data->CvDiff      = tree_dataProd->CvDiff;
    tree_data->CvIdx1      = tree_dataProd->CvIdx1;
    tree_data->CvIdx2      = tree_dataProd->CvIdx2;
    tree_data->CvDisc      = tree_dataProd->CvDisc;
    tree_data->JumpPpy     = tree_dataProd->Ppy;
    tree_data->NbDailyPts  = tree_dataProd->NbDailyPts;
    tree_data->Ppy         = tree_dataProd->Ppy;
    
    status = SUCCESS;
        
    return (status);

}  /* Cet_Manager */


/***** Cet_NmrVanl ********************************************************/
/*
*   Determine the Target Swap Atm at each Vol Bmk & the Target Swap Sml 
*   at each Sml Bmk.                  
*
*   NB: Cet will use Nmr Dates only; other Bmk Dates will be ignored.
*
*/
int Cet_NmrVanl (long          ValueDate,  /* (I) Value date              */
                 T_CURVE const *t_curve,    /* (I) Zero curve data        */
                 MKTVOL_DATA   *mvd,        /* (I/O) Mkt vol data         */
                 int           CvDiff)     /* (I) Diffuse curve number    */
{
    MQDATA  NmrMQ;               /* Local MultiQ structure for vanilla    */

    int     i, s;    
    double  Expiry;
    double  ParYield0, SwapVnPr;
    double  TotVol, NmrSwapAtm, Volguess;
    double  NmrSmile[NBVOLPARS+4];     /* Extra space for NbSig,Nck,dN,tN */
    double  Delta[5] = {0.05, 0.25, 0.50, 0.75, 0.95};

    int     status = FAILURE;    /* Error status = FAILURE initially      */  
    
    /* Initialisations */
    for (i = 0; i < mvd->NbVol; i++)
    {
        for (s = 0; s < NBSTRIKE; s++)
        {
            mvd->SwapVnVol[s][i]   = 0.0;
            mvd->SwapTreeVol[s][i] = 0.0;
        }
    }
    
    /* Compute Fwd Yield, Fwd Libor & Zero to Nmr Date from curve         */
    for (i = 0; i < mvd->NbVol; i++ )
    {
        if (ParYieldFromDates (&(mvd->ParYield0[i]),
                               &(mvd->Annuity0[i]),
                               mvd->SwapSt[i],
                               mvd->SwapMat[i],
                               mvd->DCC,
                               mvd->Freq,
                               'F',
                               &(t_curve[CvDiff])) == FAILURE)
        {
            goto RETURN;
        }
        
        
        ParYield0  = mvd->ParYield0[i];
        NmrSwapAtm = mvd->Vol[i];
        Expiry     = Daysact(ValueDate, mvd->SwapSt[i]) / 365.;


        /* If smile benchmark date */
        if (mvd->SmlLiqDate[i])
        {
            TotVol   = NmrSwapAtm * sqrt(Expiry);
            
            /* Overwrite 2nd & 4th pts with DeltaL & DeltaR, resp.    */
            Delta[1] = mvd->Smile[5][i];  
            Delta[0] = MIN(Delta[0],Delta[1]);
            
            Delta[3] = mvd->Smile[7][i];
            Delta[4] = MAX(Delta[4], Delta[3]);
            
            
            /* MultiQ: fill numerical tree parameters and load all    */
            /* into smile vector                                      */
            NmrSmile[0] = mvd->Vol[i];
            for (s=1; s<NBVOLPARS; s++) NmrSmile[s] = mvd->Smile[s][i];
            
            NmrSmile[NBVOLPARS]   = mvd->NbSigmaMQ;
            NmrSmile[NBVOLPARS+1] = mvd->NckMQ;
            NmrSmile[NBVOLPARS+2] = Q3TMX_6Q_DELTAN_DEFAULT;           
            NmrSmile[NBVOLPARS+3] = Q3TMX_6Q_TAUN_DEFAULT;
            
            if (Q3TMXSVToMQ (ParYield0,
                          NmrSwapAtm,
                          Expiry, 
                          NmrSmile+1,
                          &(NmrMQ)) == FAILURE)
            {
                DR_Error ("Calibration of MultiQ parameters failed "
                          "at date %ld.\n", mvd->SwapSt[i]);
                goto RETURN;
            }
            
            /* Calculate Vanilla Prices & Vols at each Delta for Cet  */
            for (s = 0; s < NBSTRIKE; s++)
            {
                mvd->SwapStrk[s][i] = exp (Normal_InvH(Delta[s]) * TotVol) *
                                              ParYield0;
                
                if (Q3TMXMQPricer (&(NmrMQ),
                                Q3TMX_CALL,
                                mvd->SwapStrk[s][i],
                                &(SwapVnPr)) == FAILURE)
                {
                    goto RETURN;
                }

                Volguess = NmrSwapAtm;

                if ((s != MIDSTRIKE) && (i > 0) && (mvd->SwapVnVol[s][i-1] > TINY))
                {
                    Volguess = mvd->SwapVnVol[s][i-1];
                }
                
                mvd->SwapVnVol[s][i] = ImpVol_BSQ (ParYield0,
                                                   mvd->SwapStrk[s][i],
                                                   Expiry,
                                                   SwapVnPr,
                                                   'C',
                                                   1., /*  1, 0, */
                                                   Volguess);

                if ((s >= MIDSTRIKE-1) && (s <= MIDSTRIKE+1) && (mvd->SwapVnVol[s][i] < TINY))              
                {                
                    DR_Error ("Cet_NmrVanl: Unable to calculate Vanilla Smile "
                              "at date %ld. Check MQ parameters for that expiry!\n",
                              mvd->SwapSt[i]);              
                    goto RETURN;
                }

            } /* for s < NBSTRIKE */
        }
        else
        {
            /* Only Atm Vol needed at non-Smile Bmk's */
            mvd->SwapStrk[MIDSTRIKE][i]  = ParYield0;
            mvd->SwapVnVol[MIDSTRIKE][i] = NmrSwapAtm;   
            
        } /*If smile benchmark date*/

        
    } /* for i < NbVol */

    status = SUCCESS;

RETURN:


    if (status == FAILURE)
    {
        DR_Error("Cet_NmrVanl: Failed!");
    }

    return (status);

}  /* Cet_NmrVanl */


/*****  Cet_Print **********************************************************/
/*
*  Print Residuals & Stoch Vol Parameters at Nmr Dates.
*/
int Cet_Print (int          CetIter,        /* (I) Current cet iter        */
               double      *VolPrev,        /* (I) Atm before updating     */
               MKTVOL_DATA *mvd)            /* (I) Mkt vol data            */
{
    int     status = FAILURE; /* Error status = FAILURE initially */  
    char    ErrorMsg[MAXBUFF];
    int     i, s, n;
    double  Delta[5] = {0.05, 0.25, 0.50, 0.75, 0.95};
    double  VanlVol[NBSTRIKE];
    double  TreeVol[NBSTRIKE];
    double  Expiry, ErrVol, ErrAtm;
    char    FileName[MAXBUFF], 
            FileName2[MAXBUFF];

    FILE   *stream = NULL,
           *stream2=NULL;

  


    strcpy (FileName,  "SMILE.prn");
    strcpy (FileName2, "CET.prn");

    if (CetIter == 1)
    {
        stream  = fopen (FileName, "w");
        stream2 = fopen (FileName2, "w");
    }
    else
    {
        stream  = fopen (FileName, "a");
        stream2 = fopen (FileName2, "a");
    }

    if (stream == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Cet_Print)",FileName);                 
        DR_Error (ErrorMsg);
        goto FREE_MEM_AND_RETURN;
    }     

    if (stream2 == NULL)
    {
        sprintf (ErrorMsg, "Could not open file %s! (Cet_Print)",FileName2);                 
        DR_Error (ErrorMsg);
        goto FREE_MEM_AND_RETURN;
    }     

    /* Printing SMILE                                  */
    fprintf (stream, "######################################################\n");
    fprintf (stream, "Cet Iter No:    %6d \n", CetIter);
    fprintf (stream, "                                 atm          skew");
    fprintf (stream, "           vov          tauL          tauR         ErrSmile       ErrAtm\n");

    for (i = 0; i < mvd->NbVol; i++)
    {

        /* Find corresponding numeraire index */
        for (n=1; n < mvd->NbNmr-1; n++)
        {
            if (mvd->NmrDate[n] == mvd->SwapSt[i])
                break;
        }



        if (mvd->SmlLiqDate[i])
        {
            fprintf (stream, "[%4d] %8ld Swap MQ:   %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n",
                                 i,
                                 mvd->SwapSt[i],
                                 mvd->Vol[i],
                                 mvd->Smile[1][i],
                                 mvd->Smile[2][i],
                                 mvd->Smile[6][i],
                                 mvd->Smile[8][i]
                    );

            fprintf (stream, "       %8ld Libor MQ:  %12.8f  %12.8f  %12.8f  %12.8f  %12.8f\n",
                                 mvd->SwapMat[i],
                                 mvd->NmrLibVol[0][n],
                                 mvd->NmrLibVol[1][n],
                                 mvd->NmrLibVol[2][n],
                                 mvd->NmrLibVol[6][n],
                                 mvd->NmrLibVol[8][n]
                    );

            fprintf(stream, "\n");

            Expiry = Daysact(mvd->NmrDate[0], mvd->NmrDate[n]) / 365.;
            ErrVol = 0.;
            ErrAtm = 0.;

            for (s=0; s<NBSTRIKE; s++)
            {
                VanlVol[s] = mvd->SwapVnVol[s][i];
                TreeVol[s] = mvd->SwapTreeVol[s][i];
                ErrVol     = MAX(ErrVol,fabs(TreeVol[s]-VanlVol[s]));
            }

            ErrAtm = MAX (ErrAtm, fabs (TreeVol[2]-VanlVol[2]));
            
            fprintf (stream, "                Delta:   ");

            Delta[1] = mvd->Smile[5][i];
            Delta[0] = MIN(Delta[0],Delta[1]);
            Delta[3] = mvd->Smile[7][i];
            Delta[4] = MAX(Delta[3],Delta[4]);

            for (s=0;s<NBSTRIKE;s++) 
                fprintf(stream,"%7.0f   ", Delta[s]*100);
                
            fprintf(stream, "\n");
            
            fprintf (stream, "                Strike:   ");
            
            for (s=0;s<NBSTRIKE;s++) 
                fprintf(stream,"%8.4f  ", mvd->SwapStrk[s][i]*100);

            fprintf(stream, "                        %8.2f      %8.2f", 
                    ErrVol*100, ErrAtm*100);

            fprintf(stream, "\n");
        
            fprintf (stream, "                Target:   ");

            for (s=0;s<NBSTRIKE;s++)
                fprintf(stream,"%8.2f  ",VanlVol[s]*100);

            fprintf(stream, "\n");
            
            fprintf (stream, "                Tr. Price:");

            for (s=0;s<NBSTRIKE;s++)
                fprintf(stream,"%8.2f  ",TreeVol[s]*100);

            fprintf(stream, "\n");
            
            fprintf (stream, "                Diff:     ");

            for (s=0;s<NBSTRIKE;s++)    
                fprintf(stream,"%8.2f  ",(VanlVol[s]-TreeVol[s])*100);
                
            fprintf(stream, "\n\n");

        } /* if (mvd->SmlLiqDate[i]) */

    } /* for i */

    fprintf (stream, "\n");

    /* Printing CET                   */
    fprintf (stream2, "###################################################################################\n");
    fprintf (stream2, "ITERATION:  %d\n\n",CetIter);

    fprintf (stream2, "       Vol Bmk              Tree       Vanl   VolDiff   Vol     NxtVol \n");
    fprintf (stream2, "           #                (%%)        (%%)     (%%)      (%%)      (%%) \n");


    for (i = 0; i < mvd->NbVol; i++)
    {
        fprintf (stream2,
                 "%2d (%8ld/%8ld)   %8.2f  %8.2f  %6.2f  %8.4f  (%8.4f)\n",
                 i,
                 mvd->SwapSt[i],
                 mvd->SwapMat[i],
                 100.0 * mvd->SwapTreeVol[MIDSTRIKE][i],
                 100.0 * mvd->SwapVnVol[MIDSTRIKE][i],
                 100.0 * (mvd->SwapVnVol[MIDSTRIKE][i] - mvd->SwapTreeVol[MIDSTRIKE][i]),
                 100.0 * VolPrev[i],
                 100.0 * mvd->Vol[i]
                 );
        
    } 

    fprintf (stream2, "\n\n");

    status = SUCCESS;

FREE_MEM_AND_RETURN:
 
    if (stream != NULL)
    {
        fclose (stream);
    }

    if (stream2 != NULL)
    {
        fclose (stream2);
    }

    return (status);

}  /* Cet_Print */


/*****  Cet_SwapVols  **********************************************************/
/*
*  Compute Tree Vols for Benchmark Swaps at given Deltas
*/

int   Cet_SwapVols (int            VolIndex,        /* (I) Numeraire index    */
                    int            Smoothing,       /* (I) Smoothing          */
                    double         *Annuity,        /* (I) Annuity slice      */
                    double         *InvNmr,         /* (I) Inverse Numeraire  */
                    double         *Yield,          /* (I) Yield slice        */
                    double         *Swap,           /* (I) Swap Slice         */
                    MKTVOL_DATA    *mvd,            /* (I/O) Volatility data  */
                    FIX3_TREE_DATA *tree_data,      /* (I) Tree data          */
                    FIX3_DEV_DATA  *dev_data,       /* (I) Dev data           */
                    int         t)                  /* (I) Current time point */
{
    double *AuxSlice    = NULL;
    double *Swaption    = NULL;
    double *BMarkPrice  = NULL;
    
    double  AAtoOUFact;                     /* Change of measure constant  */
    double  Expiry;
    double  Aux, Volguess;

    int     EDevIdx;                        /* Express DEV index           */
    int     s;
    int     status = FAILURE;               /* Error status */

    AAtoOUFact = mvd->TermZero0 / mvd->Annuity0[VolIndex];
    Expiry     = (double ) (Daysact(mvd->NmrDate[0], mvd->SwapSt[VolIndex])) / 365.;

    /* Get state price slice index */
    EDevIdx = GetDLOffset (tree_data->NbEDevDates,
                           tree_data->EDevDate,
                           tree_data->TPDate[t],
                           CbkEXACT);

    if (EDevIdx < 0)
    {
        DR_Error ("State prices not available for current numeraire date %ld\n",
                  tree_data->TPDate[t]);
        return (FAILURE);
    }


    /* Allocate memory */
    AuxSlice   = Fix3_Alloc_Slice(tree_data);
    Swaption   = Fix3_Alloc_Slice(tree_data);
    BMarkPrice = Fix3_Alloc_Slice(tree_data);

    if ( (AuxSlice    == NULL) ||
         (Swaption    == NULL) ||
         (BMarkPrice  == NULL) )
    {
        DR_Error("Cet_SwapVols: Cannot allocate memory!");
        goto RETURN;
    }

    /* Check memory has been properly allocated */
    if (mvd->SmlLiqDate[VolIndex])
    {
        if ( (Yield   == NULL) ||
             (Annuity == NULL) )
        {
            DR_Error("Cet_SwapVols: Can not perform Smile calculation!");
            goto RETURN;
        }
    }
    else 
    {
        if (Swap == NULL)
        {
            DR_Error("Cet_SwapVols: Can not perform Vol calculation!");
            goto RETURN;
        }
    }



    /* Smile liquid date */
    if (mvd->SmlLiqDate[VolIndex])
    {
        /* Smile liquid expiry */
        for (s = 0; s < NBSTRIKE ; s++) 
        {
            
            /* clean option slice */
            if (Fix3_Set_Slice (Swaption,
                                0.,
                                t,
                                tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            /* Set option price */
            if (Fix3_Copy_Slice (BMarkPrice,
                                 Yield,
                                 t,
                                 tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            if (Fix3_AddScalar (BMarkPrice,
                                -mvd->SwapStrk[s][VolIndex],
                                t,
                                tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            /* Multiply by annuity */
            if (Fix3_MultiplyTwoSlices (BMarkPrice,
                                        BMarkPrice,
                                        Annuity,
                                        t,
                                        tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            /* evaluate option claim */
            if (Fix3_OptionPlus_t (Swaption,
                             NULL,          /* Statistics not needed */
                             NULL,          /* Statistics not needed */
                             NULL,          /* Statistics not needed */
                             BMarkPrice,
                             1.,            /* Notional   */
                             0.,            /* Strike     */
                             TRUE,          /* ExerFlag   */
                             1,             /* always price call = 1 */
                             Smoothing,
                             AuxSlice,
                             t,
                             t,             /* no discounting */
                             0,        /* => not used    */
                             dev_data,
                             tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            /* Divide Annuity by Numeraire for Smile calculation */
            if (Fix3_MultiplyTwoSlices (Swaption,
                                        Swaption,
                                        InvNmr,
                                        t,
                                        tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            if (Fix3_MultTwoSlicesAddAll (&Aux,
                                          Swaption,
                                          tree_data->EDevStPrice[EDevIdx],
                                          t,
                                          tree_data) == FAILURE)
            {
                goto RETURN;
            }
            
            Aux  *= AAtoOUFact;
            Volguess = mvd->Vol[VolIndex];

            if (mvd->SwapTreeVol[s][VolIndex] > TINY)               
            {           
                Volguess = mvd->SwapTreeVol[s][VolIndex];               
            }
            else if (mvd->SwapVnVol[s][VolIndex] > TINY)
            {           
                Volguess = mvd->SwapVnVol[s][VolIndex];             
            }

            mvd->SwapTreeVol[s][VolIndex] = ImpVol_BSQ (mvd->ParYield0[VolIndex],
                                                        mvd->SwapStrk[s][VolIndex],
                                                        Expiry,
                                                        Aux,
                                                        'C',
                                                        1,
                                                        Volguess);

            if ((s >= MIDSTRIKE-1) && (s <= MIDSTRIKE+1) && (mvd->SwapTreeVol[s][VolIndex] < TINY))         
            {                           
                DR_Error ("Cet_SwapVols: Unable to calculate the Smile on the Tree "
                          "at date %ld. Check MQ parameters for that expiry!\n",                     
                          mvd->SwapSt[VolIndex]);               
                goto RETURN;                
            }
        }
    }
    else
    {
        /* clean option slice */
        if (Fix3_Set_Slice (Swaption,
                            0.,
                            t,
                            tree_data) == FAILURE)
        {
            goto RETURN;
        }

        /* Set option price */
        if (Fix3_Copy_Slice (BMarkPrice,
                             Swap,
                             t,
                             tree_data) == FAILURE)
        {
            goto RETURN;
        }

        
        /* evaluate option claim */
        if (Fix3_OptionPlus_t (Swaption,
                         NULL,          /* Statistics not needed */
                         NULL,          /* Statistics not needed */
                         NULL,          /* Statistics not needed */
                         BMarkPrice,
                         1.,            /* Notional   */
                         0.,            /* Strike     */
                         TRUE,          /* ExerFlag   */
                         1,             /* always price call = 1 */
                         Smoothing,
                         AuxSlice,
                         t,
                         t,             /* no discounting */
                         0,        /* => not used    */
                         dev_data,
                         tree_data) == FAILURE)
        {
            goto RETURN;
        }
        
        
        /* Divide Annuity by Numeraire for Smile calculation */
        if (Fix3_MultiplyTwoSlices (Swaption,
                                    Swaption,
                                    InvNmr,
                                    t,
                                    tree_data) == FAILURE)
        {
            goto RETURN;
        }
        
        if (Fix3_MultTwoSlicesAddAll (&Aux,
                                      Swaption,
                                      tree_data->EDevStPrice[EDevIdx],
                                      t,
                                      tree_data) == FAILURE)
        {
            goto RETURN;
        }


        Aux  *= AAtoOUFact;        
        Volguess = mvd->Vol[VolIndex];

        if (mvd->SwapTreeVol[MIDSTRIKE][VolIndex] > TINY)                       
        {                   
            Volguess = mvd->SwapTreeVol[MIDSTRIKE][VolIndex];                           
        }
        else if (mvd->SwapVnVol[MIDSTRIKE][VolIndex] > TINY)        
        {                   
            Volguess = mvd->SwapVnVol[MIDSTRIKE][VolIndex];             
        }

        mvd->SwapTreeVol[MIDSTRIKE][VolIndex] = ImpVol_BSQ (mvd->ParYield0[VolIndex],
                                                            mvd->ParYield0[VolIndex],
                                                            Expiry,
                                                            Aux,
                                                            'C',
                                                            1,
                                                            Volguess);

        if (mvd->SwapTreeVol[MIDSTRIKE][VolIndex] < TINY)                   
        {                                   
            DR_Error ("Cet_SwapVols: Unable to calculate the ATM Vol on the Tree "
                      "at date %ld. Check MQ parameters for that expiry!\n",                                         
                      mvd->SwapSt[VolIndex]);                           
            goto RETURN;                            
        }
    }

    status = SUCCESS;

 

RETURN:

    Fix3_Free_Slice (AuxSlice, tree_data);
    Fix3_Free_Slice (Swaption, tree_data);
    Fix3_Free_Slice (BMarkPrice, tree_data);


    if (status != SUCCESS)
    {        
        DR_Error("Cet_SwapVols: Failed!");
    }

    return (status);

} /* Cet_SwapVols */


/*****  Cet_Schedule ******************************************************/
/*
*   Sets up the time line using Nmr Dates. 
*
*   This routine deals with the detailed date generation from input
*   parameters and the subsequent addition of dates to the timeline
*   critical date list.
*
*/
int     Cet_Schedule 
            (long           ValueDate,     /* (I) Value date               */
             MKTVOL_DATA    *mktvol_data,  /* (I) Market Vol data for CET  */
             FIX3_TREE_DATA *tree_data)    /* (O) Tree data for CET        */
{
    int         NbCritDate = 0;
    CRIT_DATE   *CritDate = NULL;        /* Critical date list             */

    EVENT_LIST  *PmtEventListFix[NBCRITDATE]; /* Event list for fixed pmts */

    long        CurrCouponDate;          /* Date of prossesed coupon       */
    double      CurrCouponDayCount;      /* Day count of prossesed coupon  */

    
    int         i, k;                    /* Index for convenience          */
    int         status       = FAILURE;  /* Status = FAILURE initially     */
    int         MemAllocCrit = FALSE;
    int         MemAllocCrit2[NBCRITDATE];
    
    
    for (i = 0; i < NBCRITDATE; i++) MemAllocCrit2[i] = FALSE;


    /* Check to ensure that there are enough critical dates */
    /*       for all chosen benchmark instruments           */
    if (mktvol_data->NbVol > NBCRITDATE)
    {
        DR_Error("The nb of swaptions to calibrate exceeds the pre-set\n"
                 "memory allocated for it. (Cet_Schedule)\n");
        goto FREE_MEM_AND_RETURN;
    }

    
    /* Initialise all event lists to NULL */
    for (k=0; k< mktvol_data->NbVol; k++) 
    {
        PmtEventListFix[k] = NULL;
    }

    /* CRITICAL DATE LIST */
    /*  Allocate empty critical date list and then add dates successively.*/
    CritDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
    
    if (CritDate == NULL)
    {
        DR_Error("Cet_Schedule: unable to allocate memory for critical "
                    "dates array !");
        goto FREE_MEM_AND_RETURN;
    }

    MemAllocCrit = TRUE;

    /************************************/
    /*                                  */
    /*  TYPES 0 TO NBCRITDATE:  N/A     */
    /*                                  */ 
    /************************************/                       
    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->CritType[i] = 'D';
        tree_data->NbZeros[i]  =  0;
    }

    /********************************/
    /*                              */
    /*      FIXED COUPON DATES      */
    /*    (Final princ included     */
    /********************************/
    for (k = 0; k < mktvol_data->NbVol; k++)
    {
        /* Process the swap inputs to generate the appropriate fixed cp  */
        /* payment dates which will be placed in a temporary EVENT_LIST. */
        long TempDates[2];  /* Only accrual start and final mat */
        
        TempDates[0] = mktvol_data->SwapSt[k];
        TempDates[1] = mktvol_data->SwapMat[k];
        
        PmtEventListFix[k] = DrNewEventListFromFreq 
                                  (2,
                                   TempDates,
                                   mktvol_data->Freq,
                                   'F',    /* Forward stub                  */
                                   'N',    /* 'Dates in' check not required */
                                   NULL, NULL, NULL, NULL, NULL);
        
        if (PmtEventListFix[k] == NULL)
        {
            goto FREE_MEM_AND_RETURN;
        }

        MemAllocCrit2[k] = TRUE;

        for (i = 1; i < PmtEventListFix[k]->NbEntries; i++)
        {
            CurrCouponDate = PmtEventListFix[k]->Dates[i];

            if (DrDayCountFraction (PmtEventListFix[k]->Dates[i-1],
                                    PmtEventListFix[k]->Dates[i],
                                    mktvol_data->DCC,
                                    &CurrCouponDayCount) == FAILURE)
            {
                goto FREE_MEM_AND_RETURN;
            }
            
            if (Add_To_DateList (&NbCritDate,
                                 &CritDate,
                                 PmtEventListFix[k]->Dates[i],
                                 k,
                                 CurrCouponDayCount,
                                 0, 0, 0, 0,
                                 0, 0, 0) == FAILURE)     
            {
                goto FREE_MEM_AND_RETURN;
            }

        }  /* for i */

    } /* For k */

    /********************************/
    /*  TMX: Make numeraire data    */ 
    /********************************/              

    /* Add Nmr dates to Critical List  */
    for (i = 0; i <mktvol_data->NbNmr; i++)                
    {
        if (Add_To_DateList (&NbCritDate,
                             &CritDate,
                             mktvol_data->NmrDate[i],
                             NMREVENT,
                             i, 0, 0, 0, 0, 
                             0, 0, 0) == FAILURE)
        {
            goto FREE_MEM_AND_RETURN;
        }

    }  /* for i */

    /* 
    *   Finally construct the time line using 'J' time steps
    *   i.e. increasing time steps with modified jump size
    */
    if (Fix3_Time_Line (ValueDate,                      
                        NbCritDate,
                        CritDate,
                        'J',                            
                        tree_data) == FAILURE)
    {
        goto FREE_MEM_AND_RETURN;
    }

    status = SUCCESS;

FREE_MEM_AND_RETURN:                            

    if (MemAllocCrit)
    {
        Free_DR_Array (CritDate, CRITDATE, 0, NbCritDate-1);
    }

    for (k=0; k<mktvol_data->NbVol; k++)
    {
        if (MemAllocCrit2[k])
        {
            DrFreeEventList(PmtEventListFix[k]);
        }
    }

    if (status == FAILURE)
    {
        DR_Error("Cet_Schedule: Failed.\n");
    }

    return (status);

}  /* Cet_Schedule */



/*****  Cet_UpdSwapSml **********************************************************/
/*
*    Given the current Swap SV-MQ Smile - implied by the current Libor SV-MQ
*    distribution -, calculate the error wrt Target & adjust it accordingly:
*    firstly, adjust the AtmVol, then the Skew & VoV.
*    (NB:  TauL, TauR, DeltaL & DeltaR remain unchanged at their intial values).
*/
int Cet_UpdSwapSml (long        ValueDate,    /* (I) Value date                 */
                    MKTVOL_DATA *mvd)         /* (I/O) Mkt vol data             */
{
    int     i;            
    double  Expiry,RRTree, RRVanl, StrTree, StrVanl, diff, der, deltaStrk, totVol;
    int     status = FAILURE;               /* Error status = FAILURE initially */  
 

    /*  Update target diffs and improve Atm Vols on Bmk Dates           */
    for (i = 0; i < mvd->NbVol; i++)
    {

        /* Expiry */
        Expiry   = Daysact(ValueDate, mvd->SwapSt[i]) / 365.;


        /* 
         * FIRST UPDATE THE ATM VOL FOR THE SWAPS 
         */
        mvd->Vol[i] += mvd->SwapVnVol[MIDSTRIKE][i] - mvd->SwapTreeVol[MIDSTRIKE][i];



        /* Smile liquid date */
        if (mvd->SmlLiqDate[i] && mvd->CalibSmileFlag)
        {
            totVol = mvd->Vol[i] * sqrt(Expiry);

            /* 
             * NOW UPDATE SKEW
             */
            
            RRTree    = mvd->SwapTreeVol[MIDSTRIKE-1][i] - mvd->SwapTreeVol[MIDSTRIKE+1][i];
            RRVanl    = mvd->SwapVnVol[MIDSTRIKE-1][i]   - mvd->SwapVnVol[MIDSTRIKE+1][i];
            deltaStrk = mvd->SwapStrk[MIDSTRIKE-1][i]    - mvd->SwapStrk[MIDSTRIKE+1][i];
            diff   = RRTree - RRVanl; 
            
            /* dRR/dSkew = ATM * dK/(2Y0) */
            der    = - mvd->SwapTreeVol[MIDSTRIKE][i] / (2.0 * mvd->ParYield0[i]) * deltaStrk;

            if (fabs(der) > TINY && fabs(diff) > CET_WARNING_TOLERANCE / 100.0)
            {
                mvd->Smile[1][i] -= diff/der;
            }

            /* Apply MultiQ limits to skew */
            mvd->Smile[1][i] = MAX(mvd->Smile[1][i], 
                                 MAX(1. - Q3TMX_MAX_SKEW/totVol,(-1.0)*TMX_MAX_SKEW));
            mvd->Smile[1][i] = MIN(mvd->Smile[1][i], 
                                 MIN(1. + Q3TMX_MAX_SKEW/totVol,TMX_MAX_SKEW));



            /* 
             * NOW UPDATE VOV
             */
            StrTree = mvd->SwapTreeVol[MIDSTRIKE-1][i] + mvd->SwapTreeVol[MIDSTRIKE+1][i] 
                            -  2.0 * mvd->SwapTreeVol[MIDSTRIKE][i];
            StrVanl = mvd->SwapVnVol[MIDSTRIKE-1][i]   + mvd->SwapVnVol[MIDSTRIKE+1][i] 
                            -  2.0 * mvd->SwapVnVol[MIDSTRIKE][i];          
            diff    = StrTree - StrVanl;
            
            if ((fabs(mvd->Smile[2][i]) > 1.0e-04) && 
                (fabs(diff) > CET_WARNING_TOLERANCE / 100.0) && 
                (fabs(StrVanl) > TINY))
            {
                mvd->Smile[2][i] -= diff * mvd->Smile[2][i] / fabs(StrVanl);     
            }

            /* Apply MultiQ limits to vvol */
            mvd->Smile[2][i] = MAX(mvd->Smile[2][i], TINY);
            mvd->Smile[2][i] = MIN(mvd->Smile[2][i], 
                             MIN(Q3TMX_MAX_VOL_VOL/(sqrt(Expiry/3.)),TMX_MAX_VOL_VOL));



        } /* if Smile liquid date */
    }

    status = SUCCESS;

    if (status == FAILURE)
    {
        DR_Error("Cet_UpdSwapSml: Failed!");
    }

    return (status);

}  /* Cet_UpdSwapSml */




