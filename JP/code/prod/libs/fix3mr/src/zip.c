/**************************************************************************/
/*      Payoff function for the zip                                       */
/**************************************************************************/
/*      zip.c                                                             */
/**************************************************************************/


/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/src/zip.c,v 1.1 1999/01/05 14:52:12 plewicki Exp $
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fix123head.h"


/*****  Zip_t  ***********************************************************/
/*
*   Calculates the zip price
*/
int  Zip_t
             (double     *Zip,            /* (O) Price                   */
              long        AmortFlag,      /* (I) Princ pmt at timept     */
              long        CmpFlag,        /* (I) Compounding at timept   */
              double      Notional,       /* (I) Notinal of deal         */
              double     *AmortIndex,     /* (I) Amort index             */
              double      FixCpn,         /* (I) Fix cpn amt             */
              double      AmortBase,      /* (I) Base rate               */
              double      AmortFactor,    /* (I) Amort factor            */
              double      *ATRate,        /* (I) Stoch amort rate values */
              double      *ATSprd,        /* (I) Stoch amort rate spreads*/
              int         ATDim,          /* (I) Stoch amort rate dim    */
              int         t,              /* (I) Current time point      */
              int         T,              /* (I) Last time point         */
              int         DCurve,         /* (I) Discount curve          */
              DEV_DATA    *dev_data,      /* (I) Dev data structure      */
              TREE_DATA   *tree_data)     /* (I) Tree data structure     */
{
    /* Local slice pointers */

    double  *ZipL;
    double  *AmortIndexL;

    /* payoff variables */

    double  Rate;

    /* Tree variables */
                                                                        
    int     Top1, Bottom1;            /* Tree limits (1rst dim) */
    int     *Top2, *Bottom2;          /* Tree limits (2nd dim)  */
    int     **Top3, **Bottom3;        /* Tree limits (3rd dim)  */

    int     i, j, k;                  /* Node indices           */
    int     offset;                   /* Node offset            */
    int     status = FAILURE;         /* Error status           */


    Top1    = tree_data->Top1[t];
    Top2    = tree_data->Top2[t];
    Top3    = tree_data->Top3[t];
    Bottom1 = tree_data->Bottom1[t];
    Bottom2 = tree_data->Bottom2[t];
    Bottom3 = tree_data->Bottom3[t];

    /* Discount all the state variables */
    if (Dev(Zip,
            t,
            T,
            DCurve,
            dev_data,
            tree_data) == FAILURE)
    {
        goto RETURN;
    }

    /* Compound */

    if (CmpFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);
            ZipL   = Zip + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {   
                ZipL[i] *= (1. + FixCpn);

            } /* for i */
        } 
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);
                ZipL   = Zip + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    ZipL[j] *= (1. + FixCpn);

                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);
                    ZipL   = Zip + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        ZipL[k] *= (1. + FixCpn);

                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if CmpFlag */

    /* Amortize */

    if (AmortFlag)
    {
        if (tree_data->NbFactor == 1)
        {
            offset = Node_Offset(1, 0, 0, t, tree_data);
            AmortIndexL = AmortIndex  + offset;
            ZipL = Zip + offset;

            for (i = Bottom1; i <= Top1; i ++)
            {                
                /* find stochastic amortization rate */
                tableinterp(AmortIndexL[i]-AmortBase,
                            &Rate,
                            ATSprd,
                            ATRate,
                            ATDim);
                Rate *= AmortFactor;

                ZipL[i] = ZipL[i] * (1. - Rate) + Rate * Notional;

            } /* for i */

        } /* if NbFactor == 1 */
        else if (tree_data->NbFactor == 2)
        {
            for (i = Bottom1; i <= Top1; i ++)
            {
                offset = Node_Offset(2, i, 0, t, tree_data);
                AmortIndexL = AmortIndex  + offset;
                ZipL = Zip + offset;

                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    /* find stochastic amortization rate */
                    tableinterp(AmortIndexL[j]-AmortBase,
                                &Rate,
                                ATSprd,
                                ATRate,
                                ATDim);
                    Rate *= AmortFactor;

                    ZipL[j] = ZipL[j] * (1. - Rate) + Rate * Notional;
                }  /* for j */
            }  /* for i */
        }
        else if (tree_data->NbFactor == 3)
        {
            for (i = Bottom1; i <= Top1; i ++)
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    offset = Node_Offset(3, i, j, t, tree_data);
                    AmortIndexL = AmortIndex  + offset;
                    ZipL = Zip + offset;

                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        /* find stochastic amortization rate */
                        tableinterp(AmortIndexL[k]-AmortBase,
                                    &Rate,
                                    ATSprd,
                                    ATRate,
                                    ATDim);
                        Rate *= AmortFactor;
        
                        ZipL[k] = ZipL[k] * (1. - Rate) + Rate * Notional;
                    }  /* for k */
                }  /* for j */
        }  /* if then else */
    }  /* if AmortFlag */


    status = SUCCESS;
    
RETURN:

    return (status);

}  /* ZipSwap_t */
