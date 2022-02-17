/****************************************************************************/
/*      Calculation of callable zero swap price in the lattice.             */
/****************************************************************************/
/*      CALLZERO.c                                                          */
/****************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/callzero.c,v 1.1 1998/04/16 10:59:49 plewicki Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"


/*****  Callzero_t  *********************************************************/
/*
*       Callable zero swap.
*/
int     Callzero_t (double      *Callzero,     /* (I/O) Callable swap     */
                    double      *Index,        /* (I) Compounding index   */
                    double      *Interest,     /* (I) Interest            */
                    double      *Funding,      /* (I) Funding             */
                    double      *TotFunding,   /* (I) Total funding (adv) */
                    long        IntrFlag,      /* (I) TRUE if refix       */
                    long        FundFlag,      /* (I) TRUE if refix       */
                    long        ExerFlag,      /* (I) TRUE if exercise    */
                    double      Strike,        /* (I) Strike for exercise */
                    double      DayCntFtn,     /* (I) Day count fraction  */
                    double      Spread,        /* (I) Spread              */
                    char        CoS,           /* (I) 'C'mp or 'S'imple   */
                    char        IntrArrears,   /* (I) 'Y' for arrears     */
                    char        FundArrears,   /* (I) 'Y' for arrears     */
                    int         t,             /* (I) Current time point  */
                    int         T,             /* (I) Last time point     */
                    int         DCurve,        /* (I) Discount curve      */
                    DEV_DATA    *dev_data,     /* (I) Dev data structure  */
                    TREE_DATA   *tree_data)    /* (I) Tree data structure */
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


    if (Dev (   Callzero,
                t,
                T,
                DCurve,
                dev_data,
                tree_data) == FAILURE)
    {
        goto RETURN;                    
    }

    if (IntrArrears == 'N')
    {
        if (Dev (   TotFunding,
                    t,
                    T,
                    DCurve,
                    dev_data,
                    tree_data) == FAILURE)
        {
            goto RETURN;                    
        }
    }

    if (tree_data->NbFactor == 1)
    {
        offset = Node_Offset(1, 0, 0, t, tree_data);

        CallzeroL    = Callzero   + offset;
        IndexL       = Index   + offset;
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
            } /* if compunding */

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
            } /* if compunding */

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
    else if (tree_data->NbFactor == 2)
    {
        if (IntrArrears == 'Y')
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                        offset = Node_Offset(2, i, 0, t, tree_data);
        
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
                        offset = Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= pow(1. + IndexL[j] + Spread, DayCntFtn);
                        }
                    }
                }
            } /* if compunding */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                        offset = Node_Offset(2, i, 0, t, tree_data);
        
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
                        offset = Node_Offset(2, i, 0, t, tree_data);
        
                        CallzeroL = Callzero + offset;
                        IndexL    = Index    + offset;
            
                        for (j = Bottom2[i]; j <= Top2[i]; j++)
                        {
                            CallzeroL[j] *= pow(1. + IndexL[j] + Spread, DayCntFtn);
                        }
                    }
                }
            } /* if compunding */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                {
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

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
                    offset = Node_Offset(2, i, 0, t, tree_data);

                    CallzeroL = Callzero + offset;

                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        CallzeroL[j] = MAX(CallzeroL[j], -Strike);
                    }
                }
            } /* if exercise */
        }
    }
    else if (tree_data->NbFactor == 3)
    {
        if (IntrArrears == 'Y')
        {
            if (FundFlag && (FundArrears == 'N'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
            
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
                        offset = Node_Offset(3, i, j, t, tree_data);
    
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
                            offset = Node_Offset(3, i, j, t, tree_data);
        
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
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= pow(1. + IndexL[k] + Spread, DayCntFtn);
                            }
                        }
                }
            } /* if compunding */

            if (FundFlag && (FundArrears == 'Y'))
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);
    
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
                        offset = Node_Offset(3, i, j, t, tree_data);

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
                        offset = Node_Offset(3, i, j, t, tree_data);
            
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
                            offset = Node_Offset(3, i, j, t, tree_data);
        
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
                            offset = Node_Offset(3, i, j, t, tree_data);

                            CallzeroL = Callzero + offset;
                            IndexL    = Index    + offset;
            
                            for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                            {
                                CallzeroL[k] *= pow(1. + IndexL[k] + Spread, DayCntFtn);
                            }
                        }
                }
            } /* if compunding */

            if (IntrFlag)
            {
                for (i = Bottom1; i <= Top1; i ++)
                    for (j = Bottom2[i]; j <= Top2[i]; j++)
                    {
                        offset = Node_Offset(3, i, j, t, tree_data);

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
                        offset = Node_Offset(3, i, j, t, tree_data);
    
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
                        offset = Node_Offset(3, i, j, t, tree_data);

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
                        offset = Node_Offset(3, i, j, t, tree_data);
    
                        CallzeroL = Callzero + offset;

                        for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)               
                        {
                            CallzeroL[k] = MAX(CallzeroL[k], -Strike);
                        }
                    }
            } /* if exercise */
        }
    }  /* if NbFactor ... */

    
    status = SUCCESS;

    RETURN:

    return (status);

}  /* Callzero_t */



