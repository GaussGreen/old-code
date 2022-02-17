/****************************************************************************/
/*      CALLTURBO.C                                                         */
/****************************************************************************/
/*                                                                          */
/*      Calculation of cancellable turbo floater in the lattice.            */
/*                                                                          */
/****************************************************************************/



/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


/*****  Hyb3_Callturbo_t  *********************************************************/
/*                                                                           */
/*       Cancellable turbo. Strike amount is interpreted as a penalty paid   */
/*       upon cancellation.                                                  */
/*                                                                           */
/*       The implied strike variable is offered in order to enable the user  */
/*       to call  the function with a 3-D correction to the strike  (e.g if  */
/*       notional exchanges must be accounted for).                          */
/*                                                                           */
/*       Assumes dimensions as follows:                                      */
/*           FX rates :              3-Dim                                   */
/*           Turbo pmt :             3-Dim                                   */
/*           Funding pmt :           3-Dim                                   */
/*           Implied strike :        3-Dim (for full generality)             */
/*                                                                           */
int    Hyb3_Callturbo_t(TSLICE      Callturbo,     /* (I/O) Callable swap         */
                   TSLICE      TurboPmt,      /* (I) Turbo payment (if reset)*/
                   TSLICE      FundPmt,       /* (I) Funding payment ( " )   */
                   int         TurboRFlag,    /* (I) TRUE if reset of turbo  */
                   int         FundRFlag,     /* (I) TRUE if reset of funding*/
                   double      Strike,        /* (I) Strk amount in dom ccy  */
                   int         ExerFlag,      /* (I) TRUE if exercise        */
                   TSLICE      ImpliedStrike, /* (I) Auxiliary strike amount */
                   char        TurboArrears,  /* (I) 'Y' for arrears         */
                   char        FundArrears,   /* (I) 'Y' for arrears         */

                   int         DPrincipalFlag,/* (I) Flag for dom principal  */
                   double      DPrincipal,    /* (I) Domestic princ pmt      */
                   int         FPrincipalFlag,/* (I) Flag for foreign princ  */
                   double      FPrincipal,    /* (I) Foreign princ pmt       */
                   TSLICE      FxSpotPtr,     /* (I) FX rates at time point  */     

                   int         t,             /* (I) Current time point      */
                   int         T,             /* (I) Last time point         */
                   int         DCurve,        /* (I) Discount curve          */
          
                   HYB3_DEV_DATA    *dev_data,     /* (I) Hyb3_Dev data structure      */
                   HYB3_TREE_DATA   *tree_data)    /* (I) Tree data structure     */
{





    double  *CallturboL;               /* Local slice pointers */
    
    double  *TurboL;
    double  *FundL;

    double  *ImpStrikeL;
    double  *FXRate;
 
        
    int       Top1,   Bottom1;          /* Tree limits (1rst dim) */
    int      *Top2,  *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;          /* Tree limits (3rd dim)  */

    int     i, j, k;                    /* Node indices           */
    int     offset;                     /* Node offset            */
    int     status = FAILURE;           /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];


    /* Hyb3_Dev the cancellable floater  */
    if (Hyb3_Dev(Callturbo,
            t,
            T,
            DCurve,
            DISC_3D_CUPS,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;                    
    }

  

    /*  Deal with principal exchanges first     */
    /* NB: The  full use of the strike variable */
    /* will  enable the  caller to  account for */
    /* exchanges of principal when dealing with */
    /* the exercise value.                      */
    if (FPrincipalFlag)  /* Foreign principal   */
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                CallturboL = (double *)Callturbo + offset;
                FXRate     = (double *)FxSpotPtr + offset;

                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {
                    CallturboL[k] += FPrincipal * FXRate[k];

                }  /* for k */
            }
        }	
    }  /* if */

    if (DPrincipalFlag) /* Domestic principal */
    {
        for (i = Bottom1; i <= Top1; i ++)
        {
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                CallturboL= (double *)Callturbo + offset;
                
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)	            
                {
                    CallturboL[k] += DPrincipal;

                }  /* for k */
            }
        }
    } 


    /* Now deal with addition of accruals in conjunction  */
    /* with the cancellation exercise                     */
    if (TurboArrears == 'Y')
    {

        if (FundRFlag && (FundArrears == 'N'))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
        
                    CallturboL = Callturbo + offset;
                    FundL      = FundPmt   + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k] += FundL[k]; /* $$$ funding is ADDED */
                    }
                } /* for j */
        } /* if funding in advance */

        if (ExerFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                    CallturboL = Callturbo     + offset;
                    ImpStrikeL = ImpliedStrike + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        /* Call pay-off including notional exchanges */
                        CallturboL[k] = MAX(CallturboL[k],
                                              (ImpStrikeL[k]-Strike));
                        
                    }
                }
        } /* if exercise */

        if (FundRFlag && (FundArrears == 'Y'))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                    CallturboL = Callturbo + offset;
                    FundL      = FundPmt   + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k] += FundL[k]; /* $$$ ADDED */
                    }
                }
        } /* if funding in arrears */

        if (TurboRFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                    CallturboL = Callturbo + offset;
                    TurboL     = TurboPmt  + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k] += TurboL[k];
                    }
                }
        } /* if interest */
    }

    else /* Turbo reset in advance */
    {
        if (FundRFlag && (FundArrears == 'N'))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
        
                    CallturboL = Callturbo + offset;
                    FundL      = FundPmt   + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k] += FundL[k];
                    }
                } /* for j */
        } /* if funding in advance */

        if (TurboRFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                    CallturboL = Callturbo + offset;
                    TurboL     = TurboPmt  + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k]   += TurboL[k];
                    }
                }
        } /* if interest */

        if (ExerFlag)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                    CallturboL = Callturbo + offset;
                    ImpStrikeL = ImpliedStrike + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        /* Call pay-off including notional exchanges */
                        CallturboL[k] = MAX(CallturboL[k],
                                              (ImpStrikeL[k]-Strike));
                    }
                }
        } /* if exercise */

        if (FundRFlag && (FundArrears == 'Y'))
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                    CallturboL = Callturbo + offset;
                    FundL      = FundPmt   + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                    {
                        CallturboL[k] += FundL[k];
                    }
                }
        } /* if funding in advance and interest day ! */
    }
    

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Callturbo_t */



