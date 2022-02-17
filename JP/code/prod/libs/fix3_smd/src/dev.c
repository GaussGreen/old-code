/****************************************************************************/
/*      Calculate discounted expected value in the tree.                    */
/****************************************************************************/
/*      DEV.c                                                               */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"



/*****  Fix3_Dev  ****************************************************************/
/*
*       Discounted expected value function.
*/
int     Fix3_Dev (   double     *Price,       /* (I/O) Prices to be discounted */
                int         t,           /* (I) Current time point        */
                int         T,           /* (I) Last time point           */
                int         DCurve,      /* (I) Discount curve            */
                FIX3_DEV_DATA   const* dev_data,	 /* (I) Fix3_Dev data structure        */
                FIX3_TREE_DATA  const* tree_data)   /* (I) Tree data structure       */
{

    double  *Discount, *DiscountL;          /* One period discount factors */
    double  *PriceL, *NewPriceL;            /* Local slice pointers        */
    double  *Price0, *Price1, *Price2;
    double  *Price00, *Price01, *Price02;
    double  *Price10, *Price11, *Price12;
    double  *Price20, *Price21, *Price22;

    double  *pu, *p0, *pd;
    double  *quu, *qu0, *qud, *q0u, *q00, *q0d, *qdu, *qd0, *qdd;
    double  Quu, Qu0, Qud, Q0u, Q00, Q0d, Qdu, Qd0, Qdd;
    double  *ru, *r0, *rd;

    double  PFlat;

    int     *Shift1, *Shift2, *Shift3;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */
    int     NxtTop1, NxtBottom1;            /* idem for time t+1      */
    int     *NxtTop2, *NxtBottom2;
    int     **NxtTop3, **NxtBottom3;
    int     OutTop1, OutBottom1;            /* Outer ellipsoid at time t+1 */
    int     *OutTop2, *OutBottom2;
    int     **OutTop3, **OutBottom3;

    int     i;                              /* Node indices */
    int     j, j1, j2, j3;
    int     k, k1, k2, k3;
    int     offset;                         /* Node offset */
    int     l, m;                           /* Node branching indices */
    int     status = FAILURE;               /* Error status */


    /* Nothing to do at the back of the tree */

    if (t == T)                                                        
    {
        return (SUCCESS);
    }

        
    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    NxtTop1    = tree_data->Top1[t+1];
    NxtTop2    = tree_data->Top2[t+1];
    NxtTop3    = tree_data->Top3[t+1];
    NxtBottom1 = tree_data->Bottom1[t+1];
    NxtBottom2 = tree_data->Bottom2[t+1];
    NxtBottom3 = tree_data->Bottom3[t+1];

    OutTop1    = tree_data->OutTop1[t+1];
    OutTop2    = tree_data->OutTop2[t+1];
    OutTop3    = tree_data->OutTop3[t+1];
    OutBottom1 = tree_data->OutBottom1[t+1];
    OutBottom2 = tree_data->OutBottom2[t+1];
    OutBottom3 = tree_data->OutBottom3[t+1];


    /* 
     *  Choose the discount curve: index, cost of fund or risky.
     */
    if ((DCurve == 0) || (DCurve == 1) || (DCurve == 2))
    {
        Discount = dev_data->Discount[DCurve];
    }
    else
    {
        DR_Error ("Fix3_Dev: incorrect specification for discount curve!");
        goto RETURN;
    }


    if (tree_data->NbFactor == 1)
    {
        /*
         *  Flat values for nodes between the inner ellipsoid 
         *  where values have been calculated and the outer 
         *  ellipsoid defined by OutBottom and OutTop.
         */
        PriceL = Price + Fix3_Node_Offset(1, 0, 0, t+1, tree_data);

        for (i = OutBottom1; i < NxtBottom1; i++)
        {
            PriceL[i] = PriceL[NxtBottom1];
        }

        for (i = NxtTop1+1; i <= OutTop1; i++)
        {
            PriceL[i] = PriceL[NxtTop1];
        }


        /*
         *  Discounted expected value
         */
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        Shift1 = dev_data->Shift1 + offset;

        pu = dev_data->pu + offset;
        p0 = dev_data->p0 + offset;
        pd = dev_data->pd + offset;

        DiscountL = Discount           + offset;
        NewPriceL = dev_data->NewPrice + offset;

        for (i = Bottom1; i <= Top1; i ++)                  
        {
            l = i + Shift1[i];

            NewPriceL[i] = pu[i] * PriceL[l+1] + p0[i] * PriceL[l] + pd[i] * PriceL[l-1];

            NewPriceL[i] *= DiscountL[i];
        }


        /* Put the prices back in the original array */
        PriceL = Price + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            PriceL[i] = NewPriceL[i];
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        PriceL = Price + Fix3_Node_Offset (2, NxtBottom1, 0, t+1, tree_data);

        PFlat = PriceL[NxtTop2[NxtBottom1]];

        for (i = OutBottom1; i < NxtBottom1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */

        for (i = NxtBottom1; i <= NxtTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            PFlat = PriceL[NxtBottom2[i]];

            for (j = OutBottom2[i]; j < NxtBottom2[i]; j++)
            {
                PriceL[j] = PFlat;
            }

            PFlat = PriceL[NxtTop2[i]];

            for (j = NxtTop2[i]+1; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */


        PriceL = Price + Fix3_Node_Offset (2, NxtTop1, 0, t+1, tree_data);

        PFlat = PriceL[NxtTop2[NxtTop1]];

        for (i = NxtTop1+1; i <= OutTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */


        Shift1 = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        for (i = Bottom1; i <= Top1; i++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quu = dev_data->quu + offset;
            qu0 = dev_data->qu0 + offset;
            qud = dev_data->qud + offset;
            q0u = dev_data->q0u + offset;
            q00 = dev_data->q00 + offset;
            q0d = dev_data->q0d + offset;
            qdu = dev_data->qdu + offset;
            qd0 = dev_data->qd0 + offset;
            qdd = dev_data->qdd + offset;

            Shift2 = dev_data->Shift2 + offset;

            DiscountL = Discount           + offset;
            NewPriceL = dev_data->NewPrice + offset;
                        
            l = i + Shift1[i];
                
            Price0 = Price + Fix3_Node_Offset (2, l+1, 0, t+1, tree_data);
            Price1 = Price + Fix3_Node_Offset (2, l  , 0, t+1, tree_data);
            Price2 = Price + Fix3_Node_Offset (2, l-1, 0, t+1, tree_data);
                        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {        
                j2 = j  + Shift2[j];
                j1 = j2 + 1;
                j3 = j2 - 1;
                
                NewPriceL[j] = quu[j] * Price0[j1] + qu0[j] * Price0[j2] + qud[j] * Price0[j3]
                             + q0u[j] * Price1[j1] + q00[j] * Price1[j2] + q0d[j] * Price1[j3]
                             + qdu[j] * Price2[j1] + qd0[j] * Price2[j2] + qdd[j] * Price2[j3];

                NewPriceL[j] *= DiscountL[j];
            }
        }  /* for i */


        for (i = Bottom1; i <= Top1; i++)
        {        
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            PriceL    = Price              + offset;
            NewPriceL = dev_data->NewPrice + offset;
                        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {	
                PriceL[j] = NewPriceL[j];
            }
        }  /* for i */        
    }
    else if (tree_data->NbFactor == 3)
    {
        PriceL = Price + Fix3_Node_Offset (3, NxtBottom1, NxtTop2[NxtBottom1], t+1, tree_data);

        PFlat = PriceL[NxtTop3[NxtBottom1][NxtTop2[NxtBottom1]]];

        for (i = OutBottom1; i < NxtBottom1; i++)
            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

        for (i = NxtBottom1; i <= NxtTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (3, i, NxtBottom2[i], t+1, tree_data);

            PFlat = PriceL[NxtTop3[i][NxtBottom2[i]]];

            for (j = OutBottom2[i]; j < NxtBottom2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

            for (j = NxtBottom2[i]; j <= NxtTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                PFlat = PriceL[NxtBottom3[i][j]];

                for (k = OutBottom3[i][j]; k < NxtBottom3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }

                PFlat = PriceL[NxtTop3[i][j]];

                for (k = NxtTop3[i][j]+1; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

            PriceL = Price + Fix3_Node_Offset (3, i, NxtTop2[i], t+1, tree_data);

            PFlat = PriceL[NxtTop3[i][NxtTop2[i]]];

            for (j = NxtTop2[i]+1; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */
        }  /* for i */


        PriceL = Price + Fix3_Node_Offset (3, NxtTop1, NxtTop2[NxtTop1], t+1, tree_data);

        PFlat = PriceL[NxtTop3[NxtTop1][NxtTop2[NxtTop1]]];

        for (i = NxtTop1+1; i <= OutTop1; i++)
            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */


        Shift1 = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        for (i = Bottom1; i <= Top1; i++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quu = dev_data->quu + offset;
            qu0 = dev_data->qu0 + offset;
            qud = dev_data->qud + offset;
            q0u = dev_data->q0u + offset;
            q00 = dev_data->q00 + offset;
            q0d = dev_data->q0d + offset;
            qdu = dev_data->qdu + offset;
            qd0 = dev_data->qd0 + offset;
            qdd = dev_data->qdd + offset;

            Shift2 = dev_data->Shift2 + offset;

            l = i + Shift1[i];
                
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
                Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];

                ru = dev_data->ru + offset;
                r0 = dev_data->r0 + offset;
                rd = dev_data->rd + offset;

                Shift3 = dev_data->Shift3 + offset;

                DiscountL = Discount           + offset;
                NewPriceL = dev_data->NewPrice + offset;
                        
                m = j + Shift2[j];
                
                Price00 = Price + Fix3_Node_Offset (3, l+1, m+1, t+1, tree_data);
                Price01 = Price + Fix3_Node_Offset (3, l+1, m  , t+1, tree_data);
                Price02 = Price + Fix3_Node_Offset (3, l+1, m-1, t+1, tree_data);
                Price10 = Price + Fix3_Node_Offset (3, l  , m+1, t+1, tree_data);
                Price11 = Price + Fix3_Node_Offset (3, l  , m  , t+1, tree_data);
                Price12 = Price + Fix3_Node_Offset (3, l  , m-1, t+1, tree_data);
                Price20 = Price + Fix3_Node_Offset (3, l-1, m+1, t+1, tree_data);
                Price21 = Price + Fix3_Node_Offset (3, l-1, m  , t+1, tree_data);
                Price22 = Price + Fix3_Node_Offset (3, l-1, m-1, t+1, tree_data);
                        
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    k2 = k  + Shift3[k];
                    k1 = k2 + 1;
                    k3 = k2 - 1;
                
                    NewPriceL[k] = Quu * (ru[k] * Price00[k1] + r0[k] * Price00[k2] + rd[k] * Price00[k3])
                                 + Qu0 * (ru[k] * Price01[k1] + r0[k] * Price01[k2] + rd[k] * Price01[k3])
                                 + Qud * (ru[k] * Price02[k1] + r0[k] * Price02[k2] + rd[k] * Price02[k3])
                                 + Q0u * (ru[k] * Price10[k1] + r0[k] * Price10[k2] + rd[k] * Price10[k3])
                                 + Q00 * (ru[k] * Price11[k1] + r0[k] * Price11[k2] + rd[k] * Price11[k3])
                                 + Q0d * (ru[k] * Price12[k1] + r0[k] * Price12[k2] + rd[k] * Price12[k3])
                                 + Qdu * (ru[k] * Price20[k1] + r0[k] * Price20[k2] + rd[k] * Price20[k3])
                                 + Qd0 * (ru[k] * Price21[k1] + r0[k] * Price21[k2] + rd[k] * Price21[k3])
                                 + Qdd * (ru[k] * Price22[k1] + r0[k] * Price22[k2] + rd[k] * Price22[k3]);

                    NewPriceL[k] *= DiscountL[k];
                }
            }  /* for j */
        }  /* for i */


        for (i = Bottom1; i <= Top1; i++)
        {        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {	
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                PriceL    = Price              + offset;
                NewPriceL = dev_data->NewPrice + offset;
                        
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    PriceL[k] = NewPriceL[k];
                }
            }  /* for j */                        
        }  /* for i */        
    }  /* if then else */


    status = SUCCESS;

    RETURN:

    return (status);

}  /* Fix3_Dev */


/*****  Fix3_Ev  ****************************************************************/
/*
 *       Expected value function.(no discounting)
 *       Used for stats on payoff
 */

int     Fix3_Ev (double     *Price,       /* (I/O) Values to be EV'd   */
            int         t,           /* (I)   Current time point        */
            int         T,           /* (I)   Last time point           */
            FIX3_DEV_DATA   const* dev_data,	 /* (I)   Fix3_Dev data structure        */
            FIX3_TREE_DATA  const* tree_data)   /* (I)   Tree data structure       */
{
    double  *PriceL,  *NewPriceL;            /* Local slice pointers    */
    double  *Price0,  *Price1,  *Price2;
    double  *Price00, *Price01, *Price02;
    double  *Price10, *Price11, *Price12;
    double  *Price20, *Price21, *Price22;

    double  *pu, *p0, *pd;
    double  *quu, *qu0, *qud, *q0u, *q00, *q0d, *qdu, *qd0, *qdd;
    double  Quu, Qu0, Qud, Q0u, Q00, Q0d, Qdu, Qd0, Qdd;
    double  *ru, *r0, *rd;

    double  PFlat;

    int     *Shift1, *Shift2, *Shift3;

    int     Top1, Bottom1;                  /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;                /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;              /* Tree limits (3rd dim)  */
    int     NxtTop1, NxtBottom1;            /* idem for time t+1      */
    int     *NxtTop2, *NxtBottom2;
    int     **NxtTop3, **NxtBottom3;
    int     OutTop1, OutBottom1;            /* Outer ellipsoid at time t+1 */
    int     *OutTop2, *OutBottom2;
    int     **OutTop3, **OutBottom3;

    int     i;                              /* Node indices */
    int     j, j1, j2, j3;
    int     k, k1, k2, k3;
    int     offset;                         /* Node offset */
    int     l, m;                           /* Node branching indices */
    int     status = FAILURE;               /* Error status */

    /* Nothing to do at the back of the tree */

    if (t == T)                                                        
    {
        return (SUCCESS);
    }

    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    NxtTop1    = tree_data->Top1[t+1];
    NxtTop2    = tree_data->Top2[t+1];
    NxtTop3    = tree_data->Top3[t+1];
    NxtBottom1 = tree_data->Bottom1[t+1];
    NxtBottom2 = tree_data->Bottom2[t+1];
    NxtBottom3 = tree_data->Bottom3[t+1];

    OutTop1    = tree_data->OutTop1[t+1];
    OutTop2    = tree_data->OutTop2[t+1];
    OutTop3    = tree_data->OutTop3[t+1];
    OutBottom1 = tree_data->OutBottom1[t+1];
    OutBottom2 = tree_data->OutBottom2[t+1];
    OutBottom3 = tree_data->OutBottom3[t+1];

    if (tree_data->NbFactor == 1)
    {
        /*
         *  Flat values for nodes between the inner ellipsoid 
         *  where values have been calculated and the outer 
         *  ellipsoid defined by OutBottom and OutTop.
         */
        PriceL = Price + Fix3_Node_Offset(1, 0, 0, t+1, tree_data);

        for (i = OutBottom1; i < NxtBottom1; i++)
        {
            PriceL[i] = PriceL[NxtBottom1];
        }

        for (i = NxtTop1+1; i <= OutTop1; i++)
        {
            PriceL[i] = PriceL[NxtTop1];
        }

        /*
         *  Do the expected value
         */
        offset = Fix3_Node_Offset(1, 0, 0, t, tree_data);

        Shift1 = dev_data->Shift1 + offset;

        pu = dev_data->pu + offset;
        p0 = dev_data->p0 + offset;
        pd = dev_data->pd + offset;

        NewPriceL = dev_data->NewPrice + offset;

        for (i = Bottom1; i <= Top1; i ++)                  
        {
            l = i + Shift1[i];

            NewPriceL[i] = pu[i] * PriceL[l+1] + 
                           p0[i] * PriceL[l] + 
                           pd[i] * PriceL[l-1];
        }


        /* Put the prices back in the original array */
        PriceL = Price + offset;

        for (i = Bottom1; i <= Top1; i ++)
        {
            PriceL[i] = NewPriceL[i];
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        PriceL = Price + Fix3_Node_Offset (2, NxtBottom1, 0, t+1, tree_data);

        PFlat = PriceL[NxtTop2[NxtBottom1]];

        for (i = OutBottom1; i < NxtBottom1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */

        for (i = NxtBottom1; i <= NxtTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            PFlat = PriceL[NxtBottom2[i]];

            for (j = OutBottom2[i]; j < NxtBottom2[i]; j++)
            {
                PriceL[j] = PFlat;
            }

            PFlat = PriceL[NxtTop2[i]];

            for (j = NxtTop2[i]+1; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */


        PriceL = Price + Fix3_Node_Offset (2, NxtTop1, 0, t+1, tree_data);

        PFlat = PriceL[NxtTop2[NxtTop1]];

        for (i = NxtTop1+1; i <= OutTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (2, i, 0, t+1, tree_data);

            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL[j] = PFlat;
            }
        }  /* for i */


        Shift1 = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        for (i = Bottom1; i <= Top1; i++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quu = dev_data->quu + offset;
            qu0 = dev_data->qu0 + offset;
            qud = dev_data->qud + offset;
            q0u = dev_data->q0u + offset;
            q00 = dev_data->q00 + offset;
            q0d = dev_data->q0d + offset;
            qdu = dev_data->qdu + offset;
            qd0 = dev_data->qd0 + offset;
            qdd = dev_data->qdd + offset;

            Shift2 = dev_data->Shift2 + offset;

            NewPriceL = dev_data->NewPrice + offset;
                        
            l = i + Shift1[i];
                
            Price0 = Price + Fix3_Node_Offset (2, l+1, 0, t+1, tree_data);
            Price1 = Price + Fix3_Node_Offset (2, l  , 0, t+1, tree_data);
            Price2 = Price + Fix3_Node_Offset (2, l-1, 0, t+1, tree_data);
                        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {        
                j2 = j  + Shift2[j];
                j1 = j2 + 1;
                j3 = j2 - 1;
                
                NewPriceL[j] = quu[j] * Price0[j1] + 
                               qu0[j] * Price0[j2] + 
                               qud[j] * Price0[j3] + 
                               q0u[j] * Price1[j1] + 
                               q00[j] * Price1[j2] + 
                               q0d[j] * Price1[j3] + 
                               qdu[j] * Price2[j1] + 
                               qd0[j] * Price2[j2] + 
                               qdd[j] * Price2[j3];
            }
        }  /* for i */

        for (i = Bottom1; i <= Top1; i++)
        {        
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            PriceL    = Price              + offset;
            NewPriceL = dev_data->NewPrice + offset;
                        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {	
                PriceL[j] = NewPriceL[j];
            }
        }  /* for i */        
    }
    else if (tree_data->NbFactor == 3)
    {
        PriceL = Price + 
                 Fix3_Node_Offset (3, NxtBottom1, NxtTop2[NxtBottom1], 
                              t+1, tree_data);

        PFlat = PriceL[NxtTop3[NxtBottom1][NxtTop2[NxtBottom1]]];

        for (i = OutBottom1; i < NxtBottom1; i++)
            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

        for (i = NxtBottom1; i <= NxtTop1; i++)
        {
            PriceL = Price + Fix3_Node_Offset (3, i, NxtBottom2[i], t+1, tree_data);

            PFlat = PriceL[NxtTop3[i][NxtBottom2[i]]];

            for (j = OutBottom2[i]; j < NxtBottom2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

            for (j = NxtBottom2[i]; j <= NxtTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                PFlat = PriceL[NxtBottom3[i][j]];

                for (k = OutBottom3[i][j]; k < NxtBottom3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }

                PFlat = PriceL[NxtTop3[i][j]];

                for (k = NxtTop3[i][j]+1; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */

            PriceL = Price + Fix3_Node_Offset (3, i, NxtTop2[i], t+1, tree_data);

            PFlat = PriceL[NxtTop3[i][NxtTop2[i]]];

            for (j = NxtTop2[i]+1; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */
        }  /* for i */


        PriceL = Price + 
                 Fix3_Node_Offset (3, NxtTop1, NxtTop2[NxtTop1], t+1, tree_data);

        PFlat = PriceL[NxtTop3[NxtTop1][NxtTop2[NxtTop1]]];

        for (i = NxtTop1+1; i <= OutTop1; i++)
            for (j = OutBottom2[i]; j <= OutTop2[i]; j++)
            {
                PriceL = Price + Fix3_Node_Offset (3, i, j, t+1, tree_data);

                for (k = OutBottom3[i][j]; k <= OutTop3[i][j]; k++)
                {
                    PriceL[k] = PFlat;
                }
            }  /* for j */


        Shift1 = dev_data->Shift1 + Fix3_Node_Offset(1, 0, 0, t, tree_data);

        for (i = Bottom1; i <= Top1; i++)
        {
            offset = Fix3_Node_Offset(2, i, 0, t, tree_data);

            quu = dev_data->quu + offset;
            qu0 = dev_data->qu0 + offset;
            qud = dev_data->qud + offset;
            q0u = dev_data->q0u + offset;
            q00 = dev_data->q00 + offset;
            q0d = dev_data->q0d + offset;
            qdu = dev_data->qdu + offset;
            qd0 = dev_data->qd0 + offset;
            qdd = dev_data->qdd + offset;

            Shift2 = dev_data->Shift2 + offset;

            l = i + Shift1[i];
                
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
                Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];

                ru = dev_data->ru + offset;
                r0 = dev_data->r0 + offset;
                rd = dev_data->rd + offset;

                Shift3 = dev_data->Shift3 + offset;

                NewPriceL = dev_data->NewPrice + offset;
                        
                m = j + Shift2[j];
                
                Price00 = Price + Fix3_Node_Offset (3, l+1, m+1, t+1, tree_data);
                Price01 = Price + Fix3_Node_Offset (3, l+1, m  , t+1, tree_data);
                Price02 = Price + Fix3_Node_Offset (3, l+1, m-1, t+1, tree_data);
                Price10 = Price + Fix3_Node_Offset (3, l  , m+1, t+1, tree_data);
                Price11 = Price + Fix3_Node_Offset (3, l  , m  , t+1, tree_data);
                Price12 = Price + Fix3_Node_Offset (3, l  , m-1, t+1, tree_data);
                Price20 = Price + Fix3_Node_Offset (3, l-1, m+1, t+1, tree_data);
                Price21 = Price + Fix3_Node_Offset (3, l-1, m  , t+1, tree_data);
                Price22 = Price + Fix3_Node_Offset (3, l-1, m-1, t+1, tree_data);
                        
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    k2 = k  + Shift3[k];
                    k1 = k2 + 1;
                    k3 = k2 - 1;
                
                    NewPriceL[k] = Quu * (ru[k] * Price00[k1] + 
                                          r0[k] * Price00[k2] + 
                                          rd[k] * Price00[k3]) 
                                 + Qu0 * (ru[k] * Price01[k1] + 
                                          r0[k] * Price01[k2] + 
                                          rd[k] * Price01[k3])
                                 + Qud * (ru[k] * Price02[k1] + 
                                          r0[k] * Price02[k2] + 
                                          rd[k] * Price02[k3])
                                 + Q0u * (ru[k] * Price10[k1] + 
                                          r0[k] * Price10[k2] + 
                                          rd[k] * Price10[k3])
                                 + Q00 * (ru[k] * Price11[k1] + 
                                          r0[k] * Price11[k2] + 
                                          rd[k] * Price11[k3])
                                 + Q0d * (ru[k] * Price12[k1] + 
                                          r0[k] * Price12[k2] + 
                                          rd[k] * Price12[k3])
                                 + Qdu * (ru[k] * Price20[k1] + 
                                          r0[k] * Price20[k2] + 
                                          rd[k] * Price20[k3])
                                 + Qd0 * (ru[k] * Price21[k1] + 
                                          r0[k] * Price21[k2] + 
                                          rd[k] * Price21[k3])
                                 + Qdd * (ru[k] * Price22[k1] + 
                                          r0[k] * Price22[k2] + 
                                          rd[k] * Price22[k3]);
                }
            }  /* for j */
        }  /* for i */


        for (i = Bottom1; i <= Top1; i++)
        {        
            for (j = Bottom2[i]; j <= Top2[i]; j++)
            {	
                offset = Fix3_Node_Offset(3, i, j, t, tree_data);

                PriceL    = Price              + offset;
                NewPriceL = dev_data->NewPrice + offset;
                        
                for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                {        
                    PriceL[k] = NewPriceL[k];
                }
            }  /* for j */                        
        }  /* for i */        
    }  /* if then else */


    status = SUCCESS;

    return (status);

}  /* Fix3_Ev */

