/****************************************************************************/
/*      Calculation of callable zero swap price in the lattice.             */
/****************************************************************************/
/*      CALLZERO.C                                                          */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"


/*****  Hyb3_Callzero_t  *********************************************************/
/*
*       Callable zero swap.
*/
int     Hyb3_Callzero_t(TSLICE      Callzero,     /* (I/O) Callable swap         */
                   TSLICE      Index,        /* (I) Compounding index       */
                   TSLICE      Interest,     /* (I) Interest                */
                   TSLICE      Funding,      /* (I) Funding                 */
                   TSLICE      TotFunding,   /* (I/O) Total funding (adv)   */
                   int         IntrFlag,     /* (I) TRUE if refix           */
                   int         FundFlag,     /* (I) TRUE if refix           */
                   int         ExerFlag,     /* (I) TRUE if exercise        */
                   double      Strike,       /* (I) Strike for exercise     */
                   double      DayCntFtn,    /* (I) Day count fraction      */
                   double      Spread,       /* (I) Spread                  */
                   char        CoS,          /* (I) 'C'mp or 'S'imple       */
                   char        IntrArrears,  /* (I) 'Y' for arrears         */
                   char        FundArrears,  /* (I) 'Y' for arrears         */
                   int         t,            /* (I) Current time point      */
                   int         T,            /* (I) Last time point         */
                   int         DCurve,       /* (I) Discount curve          */
                   int         DMode,
                   HYB3_DEV_DATA    *dev_data,    /* (I) Hyb3_Dev data structure      */
                   HYB3_TREE_DATA   *tree_data)   /* (I) Tree data structure     */
{

    double  *CallzeroL;               /* Local slice pointers */
    double  *IndexL;
    double  *InterestL;
    double  *FundingL;
    double  *TotFundingL;
        
    int     Top1, Bottom1;              /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;            /* Tree limits (2nd dim)  */
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


    if (Hyb3_Dev(Callzero,
            t,
            T,
            DCurve,
            DMode,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;                    
    }

    if (IntrArrears == 'N')
    {
        if (Hyb3_Dev(TotFunding,
                t,
                T,
                DCurve,
                DMode,
                dev_data,
                tree_data) == FAILURE)
        {
            goto RETURN;                    
        }
    }

    if (DMode == DISC_1D_NOCUPS)
    {
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        CallzeroL    = Callzero   + offset;
        IndexL       = Index      + offset;
        InterestL    = Interest   + offset;
        FundingL     = Funding    + offset;
        TotFundingL  = TotFunding + offset;
    
        if (IntrArrears == 'Y')
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] -= FundingL[i];
                }
            } /* if funding in advance */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] = MAX(CallzeroL[i], -Strike);
                }
            } /* if exercise */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CallzeroL[i] *= (1. + (IndexL[i] + Spread) * DayCntFtn);
                    }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CallzeroL[i] *= pow(1. + IndexL[i] + Spread, DayCntFtn);
                    }
                }
            } /* if compounding */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] -= FundingL[i];
                }
            } /* if funding in arrears */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] += InterestL[i];
                }
            } /* if interest */
        }
        else
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    TotFundingL[i] += FundingL[i];
                }
            } /* if funding in advance */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CallzeroL[i] *= (1. + (IndexL[i] + Spread) * DayCntFtn);
                    }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        CallzeroL[i] *= pow(1. + IndexL[i] + Spread, DayCntFtn);
                    }
                }
            } /* if compounding */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i]   -= TotFundingL[i];
                    TotFundingL[i]  = 0.;
                }
            } /* if funding in advance and interest day ! */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    TotFundingL[i] += FundingL[i];
                }
            } /* if funding in advance */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] += InterestL[i];
                }
            } /* if interest */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    CallzeroL[i] = MAX(CallzeroL[i], -Strike);
                }
            } /* if exercise */
        }
    }
    else if (DMode == DISC_2D_CUPS || DMode == DISC_2D_NOCUPS)
    {
        if (IntrArrears == 'Y')
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;
                    FundingL  = Funding  + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] -= FundingL[j];
                    }
                }
            } /* if funding in advance */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] = MAX(CallzeroL[j], -Strike);
                    }
                }
            } /* if exercise */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= (1. + (IndexL[j] + Spread) * DayCntFtn);
                        }
                    }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= pow(1. + IndexL[j] + Spread, DayCntFtn);
                        }
                    }
                }
            } /* if compounding */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;
                    FundingL  = Funding  + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] -= FundingL[j];
                    }
                }
            } /* if funding in arrears */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;
                    InterestL = Interest + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] += InterestL[j];
                    }
                }
            } /* if interest */
        }
        else
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    TotFundingL = TotFunding + offset;
                    FundingL    = Funding    + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        TotFundingL[j] += FundingL[j];
                    }
                }
            } /* if funding in advance */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= (1. + (IndexL[j] + Spread) * DayCntFtn);
                        }
                    }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                    {
                        offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= pow(1. + IndexL[j] + Spread, DayCntFtn);
                        }
                    }
                }
            } /* if compounding */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL    = Callzero    + offset;
                    TotFundingL  = TotFunding  + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j]   -= TotFundingL[j];
                        TotFundingL[j]  = 0.;
                    }
                }
            } /* if funding in advance and interest day ! */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    TotFundingL = TotFunding + offset;
                    FundingL    = Funding    + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        TotFundingL[j] += FundingL[j];
                    }
                }
            } /* if funding in advance */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;
                    InterestL = Interest + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] += InterestL[j];
                    }
                }
            } /* if interest */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] = MAX(CallzeroL[j], -Strike);
                    }
                }
            } /* if exercise */
        }
    }
    else if (DMode == DISC_3D_CUPS)
    {
        if (IntrArrears == 'Y')
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
            
                        CallzeroL = Callzero + offset;
                        FundingL  = Funding  + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] -= FundingL[k];
                        }
                    } /* for j */
            } /* if funding in advance */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                        CallzeroL = Callzero + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] = MAX(CallzeroL[k], -Strike);
                        }
                    }
            } /* if exercise */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
        
                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= (1. + (IndexL[k] + Spread) * DayCntFtn);
                            }
                        }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= pow(1. + IndexL[k] + Spread, DayCntFtn);
                            }
                        }
                }
            } /* if compounding */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                        CallzeroL = Callzero + offset;
                        FundingL  = Funding  + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] -= FundingL[k];
                        }
                    }
            } /* if funding in arrears */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CallzeroL = Callzero + offset;
                        InterestL = Interest + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] += InterestL[k];
                        }
                    }
            } /* if interest */
        }
        else
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
            
                        TotFundingL = TotFunding + offset;
                        FundingL    = Funding    + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            TotFundingL[k] += FundingL[k];
                        }
                    } /* for j */
            } /* if funding in advance */

            if (IntrFlag)
            {
                if (CoS == 'S')
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
        
                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= (1. + (IndexL[k] + Spread) * DayCntFtn);
                            }
                        }
                }
                else
                {
                    for (i = Bottom1; i <= Top1; i ++)
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= pow(1. + IndexL[k] + Spread, DayCntFtn);
                            }
                        }
                }
            } /* if compounding */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CallzeroL   = Callzero   + offset;
                        TotFundingL = TotFunding + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k]   -= TotFundingL[k];
                            TotFundingL[k]  = 0.;
                        }
                    }
            } /* if interest */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                        TotFundingL = TotFunding + offset;
                        FundingL    = Funding    + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            TotFundingL[k] += FundingL[k];
                        }
                    }
            } /* if funding in advance and interest day ! */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);

                        CallzeroL = Callzero + offset;
                        InterestL = Interest + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] += InterestL[k];
                        }
                    }
            } /* if interest */

            if (ExerFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Hyb3_Node_Offset(3, i, j, t, tree_data);
    
                        CallzeroL = Callzero + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] = MAX(CallzeroL[k], -Strike);
                        }
                    }
            } /* if exercise */
        }
    }  /* if DMode ... */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Callzero_t */



