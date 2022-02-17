/****************************************************************************/
/*        Memory allocation.	                                            */
/****************************************************************************/
/*        ALLOC.c                                                           */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "fix123head.h"




/*****  Fix3_Tree_Init  **********************************************************/
/*
*       Initialize tree pointers to NULL.
*/
void    Fix3_Tree_Init (FIX3_TREE_DATA    *tree_data) /* Tree building data structure */
{
    int
            i;


    tree_data->NbFactor = 0;
    tree_data->NbTP = 0;

    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->CritDate[i] = NULL;
        tree_data->TPtype[i]   = NULL;
        tree_data->CritType[i] = 'D';
    }

    tree_data->Top1    = NULL;
    tree_data->Bottom1 = NULL;
    tree_data->Top2    = NULL;
    tree_data->Bottom2 = NULL;
    tree_data->Top3    = NULL;
    tree_data->Bottom3 = NULL;

    tree_data->OutTop1    = NULL;
    tree_data->OutBottom1 = NULL;
    tree_data->OutTop2    = NULL;
    tree_data->OutBottom2 = NULL;
    tree_data->OutTop3    = NULL;
    tree_data->OutBottom3 = NULL;


    tree_data->QLeft    = NULL;
    tree_data->QRight   = NULL;
    tree_data->FwdShift = NULL;
    for (i = 0; i < 3; i++)
    {
        tree_data->ZeroCoupon[i] = NULL;
        tree_data->ZeroRate[i]   = NULL;
        tree_data->FwdRate[i]    = NULL;
        tree_data->BetaTD[i]     = NULL;
        tree_data->TermZero[i]   = NULL;
        tree_data->Width[i]      = 0;
    }

    for (i = 0; i < 6; i++)
    {
        tree_data->Aweight[i] = NULL;
    }

    tree_data->TPDate  = NULL;
    tree_data->ZCenter = NULL;
    tree_data->Length  = NULL;
    tree_data->LengthJ = NULL;

    tree_data->NbEDevDates = 0;
    tree_data->EDevDate    = NULL;
    tree_data->EDevStPrice = NULL;

    tree_data->JumpPpy = 0;
    tree_data->NbDailyPts = Nb_Daily_Pts;

    tree_data->PpyCet[0] = 24;
    tree_data->PpyCet[1] = 4;
    tree_data->PpyCet[2] = 4;

    tree_data->NbNmr    = 0;
    tree_data->NmrInv   = NULL;
    tree_data->LastZero = NULL;
    tree_data->Libor    = NULL;

    return;

}  /* Fix3_Tree_Init */



/*****  Fix3_Tree_Alloc  *********************************************************/
/*
*       Allocation of memory for the arrays constituting the tree.
*       The limits of the tree are allocated directly in Fix3_Tree_Limits.
*/
int     Fix3_Tree_Alloc (FIX3_TREE_DATA   *tree_data) /* Tree building data structure */
{
    int     i;
    int     NbTP;               /* Number of time points in the tree */
    int     status = FAILURE;   /* Error status                      */


    NbTP = tree_data->NbTP;

    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->TPtype[i]   = (int *)       DR_Array(INT,      -1, NbTP+1);
        tree_data->CritDate[i] = (CRIT_DATE *) DR_Array(CRITDATE, -1, NbTP+1);

        if (  (tree_data->TPtype[i] == NULL)
           || (tree_data->CritDate[i] == NULL))
        {
            goto RETURN;
        }
    }

    tree_data->Top1    = (int *)   DR_Array (INT,       -1, NbTP+1);
    tree_data->Bottom1 = (int *)   DR_Array (INT,       -1, NbTP+1);
    tree_data->Top2    = (int **)  DR_Array (INT_PTR,   -1, NbTP+1);
    tree_data->Bottom2 = (int **)  DR_Array (INT_PTR,   -1, NbTP+1);
    tree_data->Top3    = (int ***) DR_Array (INT_D_PTR, -1, NbTP+1);
    tree_data->Bottom3 = (int ***) DR_Array (INT_D_PTR, -1, NbTP+1);

    if (  (tree_data->Top1    == NULL)
       || (tree_data->Bottom1 == NULL)
       || (tree_data->Top2    == NULL)
       || (tree_data->Bottom2 == NULL)
       || (tree_data->Top3    == NULL)
       || (tree_data->Bottom3 == NULL))
    {
        goto RETURN;
    }

    tree_data->OutTop1    = (int *)   DR_Array (INT,       -1, NbTP+1);
    tree_data->OutBottom1 = (int *)   DR_Array (INT,       -1, NbTP+1);
    tree_data->OutTop2    = (int **)  DR_Array (INT_PTR,   -1, NbTP+1);
    tree_data->OutBottom2 = (int **)  DR_Array (INT_PTR,   -1, NbTP+1);
    tree_data->OutTop3    = (int ***) DR_Array (INT_D_PTR, -1, NbTP+1);
    tree_data->OutBottom3 = (int ***) DR_Array (INT_D_PTR, -1, NbTP+1);

    if (  (tree_data->OutTop1    == NULL)
       || (tree_data->OutBottom1 == NULL)
       || (tree_data->OutTop2    == NULL)
       || (tree_data->OutBottom2 == NULL)
       || (tree_data->OutTop3    == NULL)
       || (tree_data->OutBottom3 == NULL))
    {
        goto RETURN;
    }
 
    tree_data->QLeft    = (double *) DR_Array (DOUBLE, -1, NbTP+1);
    tree_data->QRight   = (double *) DR_Array (DOUBLE, -1, NbTP+1);
    tree_data->FwdShift = (double *) DR_Array (DOUBLE, -1, NbTP+1);

    for (i = 0; i < 3; i++)
    {
        tree_data->ZeroCoupon[i] = (double *) DR_Array (DOUBLE, -1, NbTP+1);
        tree_data->ZeroRate[i]   = (double *) DR_Array (DOUBLE, -1, NbTP+1);
        tree_data->FwdRate[i]    = (double *) DR_Array (DOUBLE, -1, NbTP+1);
        tree_data->TermZero[i]   = (double *) DR_Array (DOUBLE, -1, NbTP+1);

        if (  (tree_data->ZeroCoupon[i] == NULL)
           || (tree_data->ZeroRate[i]   == NULL)
           || (tree_data->FwdRate[i]    == NULL)
           || (tree_data->TermZero[i]   == NULL))
        {
            goto RETURN;
        }
    }

    for (i = 0; i < 6; i++)
    {
        tree_data->Aweight[i] = (double *) DR_Array (DOUBLE, -1, NbTP+1);

        if (tree_data->Aweight[i] == NULL)
        {
            goto RETURN;
        }
    }


    for (i = 0; i < 3; i++)
    {
        tree_data->BetaTD[i]  = (double *) DR_Array (DOUBLE, 0, NbTP );
        if (tree_data->BetaTD[i] == NULL)
        {
            goto RETURN;
        }
        
    }
    tree_data->TPDate  = (long *)   DR_Array (LONG,   -1, NbTP+1);
    tree_data->ZCenter = (double *) DR_Array (DOUBLE, -1, NbTP+1);
    tree_data->Length  = (double *) DR_Array (DOUBLE, -1, NbTP+1);
    tree_data->LengthJ = (double *) DR_Array (DOUBLE, -1, NbTP+1);

    if (  (tree_data->TPDate  == NULL)
       || (tree_data->ZCenter == NULL)
       || (tree_data->Length  == NULL)
       || (tree_data->LengthJ == NULL))
    {
        goto RETURN;
    }



    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Fix3_Tree_Alloc: could not allocate memory for the tree!");
    }

    return (status);

}  /* Fix3_Tree_Alloc */



/*****  Fix3_Tree_Free  **********************************************************/
/*
*       Free the arrays constituting the tree.
*/
int     Fix3_Tree_Free ( FIX3_TREE_DATA   *tree_data) /* Tree building data structure */
{
    int     i;
    int     t;
    int     NbTP;               /* Number of time points in the tree */
    int     Area = 0;


    NbTP = tree_data->NbTP;

    if (tree_data->NbFactor == 1)
    {
        Area = tree_data->Width[0];
    }
    if (tree_data->NbFactor == 2)
    {
        Area = tree_data->Width[0] * tree_data->Width[1];
    }
    if (tree_data->NbFactor == 3)
    {
        Area = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];
    }


    /* First of all, free slices living in FIX3_TREE_DATA */
    if (tree_data->EDevStPrice != NULL)
    {
        for (i = 0; i < tree_data->NbEDevDates; i++)
        {
            Fix3_Free_Slice(tree_data->EDevStPrice[i], tree_data);
        }
    }
    Free_DR_Array(tree_data->EDevStPrice,DOUBLE_PTR,0,tree_data->NbEDevDates-1);
    Free_DR_Array(tree_data->EDevDate, LONG, 0, tree_data->NbEDevDates - 1);

    if (tree_data->NmrInv != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Fix3_Free_Slice(tree_data->NmrInv[i],tree_data);
        }
        Free_DR_Array(tree_data->NmrInv,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    if (tree_data->LastZero != NULL)
    {
        Fix3_Free_Slice(tree_data->LastZero, tree_data);
    }
    if (tree_data->Libor != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_DR_Array (tree_data->Libor[i],DOUBLE,-1,Area);
        }
        Free_DR_Array(tree_data->Libor,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }

    Free_DR_Array (tree_data->QLeft,    DOUBLE, -1, NbTP+1);
    Free_DR_Array (tree_data->QRight,   DOUBLE, -1, NbTP+1);
    Free_DR_Array (tree_data->FwdShift, DOUBLE, -1, NbTP+1);

    Free_DR_Array (tree_data->LengthJ, DOUBLE, -1, NbTP+1);
    Free_DR_Array (tree_data->Length,  DOUBLE, -1, NbTP+1);
    Free_DR_Array (tree_data->ZCenter, DOUBLE, -1, NbTP+1);
    Free_DR_Array (tree_data->TPDate,  LONG,   -1, NbTP+1);

    for (i = 0; i < 6; i++)
    {
        Free_DR_Array(tree_data->Aweight[i], DOUBLE, -1, NbTP+1);
    }

    
    for (i = 0; i < 3; i++)
    {
        Free_DR_Array(tree_data->BetaTD[i], DOUBLE, 0, NbTP );
    }


   
    for (i = 0; i < 3; i++)
    {
        Free_DR_Array(tree_data->FwdRate[i],    DOUBLE, -1, NbTP+1);
        Free_DR_Array(tree_data->ZeroRate[i],   DOUBLE, -1, NbTP+1);
        Free_DR_Array(tree_data->ZeroCoupon[i], DOUBLE, -1, NbTP+1);
        Free_DR_Array(tree_data->TermZero[i],   DOUBLE, -1, NbTP+1);
    }
    /* 
     *  Free the limits of the two ellipses.
     */

    if (   (tree_data->NbFactor >= 3) 
        && (tree_data->Top1 != NULL) && (tree_data->Bottom1 != NULL)
        && (tree_data->Top2 != NULL) && (tree_data->Bottom2 != NULL)
        && (tree_data->Top3 != NULL) && (tree_data->Bottom3 != NULL))
    {
        for (t = 0; t <= NbTP; t++)
        {    
            if (   (tree_data->Top2[t]    != NULL) 
                && (tree_data->Bottom2[t] != NULL)
                && (tree_data->Top3[t]    != NULL)
                && (tree_data->Bottom3[t] != NULL))
            {
                for (i = tree_data->Bottom1[t]; i <= tree_data->Top1[t]; i++)
                {
                    Free_DR_Array ( tree_data->Top3[t][i],
                                        INT,
                                        tree_data->Bottom2[t][i],
                                        tree_data->Top2[t][i]);

                    Free_DR_Array ( tree_data->Bottom3[t][i],
                                        INT,
                                        tree_data->Bottom2[t][i],
                                        tree_data->Top2[t][i]);
                }  /* for i */
            }  /* if */

            Free_DR_Array ( tree_data->Top3[t],
                            INT_PTR,
                            tree_data->Bottom1[t],
                            tree_data->Top1[t]);

            Free_DR_Array ( tree_data->Bottom3[t],
                            INT_PTR,
                            tree_data->Bottom1[t],
                            tree_data->Top1[t]);
        }  /* for t */
    }  /* if */

    if (   (tree_data->NbFactor >= 2) 
        && (tree_data->Top1 != NULL) && (tree_data->Bottom1 != NULL)
        && (tree_data->Top2 != NULL) && (tree_data->Bottom2 != NULL))
    {
        for (t = 0; t <= NbTP; t++)
        {    
            Free_DR_Array ( tree_data->Top2[t],
                            INT,
                            tree_data->Bottom1[t],
                            tree_data->Top1[t]);

            Free_DR_Array ( tree_data->Bottom2[t],
                            INT,
                            tree_data->Bottom1[t],
                            tree_data->Top1[t]);
        }  /* for t */
    }  /* if */

    Free_DR_Array (tree_data->Top1,    INT,       -1, NbTP+1);
    Free_DR_Array (tree_data->Bottom1, INT,       -1, NbTP+1);
    Free_DR_Array (tree_data->Top2,    INT_PTR,   -1, NbTP+1);
    Free_DR_Array (tree_data->Bottom2, INT_PTR,   -1, NbTP+1);
    Free_DR_Array (tree_data->Top3,    INT_D_PTR, -1, NbTP+1);
    Free_DR_Array (tree_data->Bottom3, INT_D_PTR, -1, NbTP+1);

    if (   (tree_data->NbFactor >= 3) 
        && (tree_data->OutTop1 != NULL) && (tree_data->OutBottom1 != NULL)
        && (tree_data->OutTop2 != NULL) && (tree_data->OutBottom2 != NULL)
        && (tree_data->OutTop3 != NULL) && (tree_data->OutBottom3 != NULL))
    {
        for (t = 0; t <= NbTP; t++)
        {    
            if (   (tree_data->OutTop2[t]    != NULL) 
                && (tree_data->OutBottom2[t] != NULL)
                && (tree_data->OutTop3[t]    != NULL)
                && (tree_data->OutBottom3[t] != NULL))
            {
                for (i=tree_data->OutBottom1[t]; i<=tree_data->OutTop1[t]; i++)
                {
                    Free_DR_Array ( tree_data->OutTop3[t][i],
                                    INT,
                                    tree_data->OutBottom2[t][i],
                                    tree_data->OutTop2[t][i]);

                    Free_DR_Array ( tree_data->OutBottom3[t][i],
                                    INT,
                                    tree_data->OutBottom2[t][i],
                                    tree_data->OutTop2[t][i]);
                }  /* for i */
            }  /* if */

            Free_DR_Array ( tree_data->OutTop3[t],
                            INT_PTR,
                            tree_data->OutBottom1[t],
                            tree_data->OutTop1[t]);

            Free_DR_Array ( tree_data->OutBottom3[t],
                            INT_PTR,
                            tree_data->OutBottom1[t],
                            tree_data->OutTop1[t]);
        }  /* for t */
    }  /* if */

    if (   (tree_data->NbFactor >= 2) 
        && (tree_data->OutTop1 != NULL) && (tree_data->OutBottom1 != NULL)
        && (tree_data->OutTop2 != NULL) && (tree_data->OutBottom2 != NULL))
    {
        for (t = 0; t <= NbTP; t++)
        {    
            Free_DR_Array ( tree_data->OutTop2[t],
                            INT,
                            tree_data->OutBottom1[t],
                            tree_data->OutTop1[t]);

            Free_DR_Array ( tree_data->OutBottom2[t],
                            INT,
                            tree_data->OutBottom1[t],
                            tree_data->OutTop1[t]);
        }  /* for t */
    }  /* if */

    Free_DR_Array (tree_data->OutTop1,    INT,       -1, NbTP+1);
    Free_DR_Array (tree_data->OutBottom1, INT,       -1, NbTP+1);
    Free_DR_Array (tree_data->OutTop2,    INT_PTR,   -1, NbTP+1);
    Free_DR_Array (tree_data->OutBottom2, INT_PTR,   -1, NbTP+1);
    Free_DR_Array (tree_data->OutTop3,    INT_D_PTR, -1, NbTP+1);
    Free_DR_Array(tree_data->OutBottom3,  INT_D_PTR, -1, NbTP+1);

    for (i = 0; i < NBCRITDATE; i++)
    {
        Free_DR_Array (tree_data->CritDate[i], CRITDATE, -1, NbTP+1);
        Free_DR_Array (tree_data->TPtype[i],   INT,      -1, NbTP+1);
    }
    

    Fix3_Tree_Init(tree_data);
    return (SUCCESS);

}  /* Fix3_Tree_Free */



/*****  Fix3_Dev_Init  **********************************************************/
/*
*       Initialize FIX3_DEV_DATA pointers to NULL.
*/
void    Fix3_Dev_Init (FIX3_DEV_DATA    *dev_data)
{
    int   i;

    dev_data->Shift1 = dev_data->Shift2 = dev_data->Shift3 = NULL;

    for (i=0; i<3; i++) dev_data->Discount[i] = NULL;

    dev_data->pu  = dev_data->p0  = dev_data->pd  = NULL;
    dev_data->quu = dev_data->qu0 = dev_data->qud = NULL;
    dev_data->q0u = dev_data->q00 = dev_data->q0d = NULL;
    dev_data->qdu = dev_data->qd0 = dev_data->qdd = NULL;
    dev_data->ru  = dev_data->r0  = dev_data->rd  = NULL;

    dev_data->NewPrice = NULL;

    dev_data->NmrToCcy = FALSE;
    dev_data->CcyToNmr = FALSE;

    return;

}  /* Fix3_Dev_Init */


/*****  Fix3_Dev_Alloc  **********************************************************/
/*
*       Allocation of memory for the dev structure. We allocate arrays and
*       the addressing will be done linearly.
*/
int     Fix3_Dev_Alloc ( FIX3_DEV_DATA   *dev_data,  /* (I) Fix3_Dev data structure  */
                    FIX3_TREE_DATA const* tree_data) /* (I) Tree data structure */
{
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;
    int     status = FAILURE;   /* Error status */


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];

    if (tree_data->NbFactor == 1)
    {
        dev_data->NewPrice = (double *) DR_Array (DOUBLE, 0, Area1);

        dev_data->Discount[0] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->Discount[1] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->Discount[2] = (double *) DR_Array (DOUBLE, 0, Area1);

        dev_data->NmrInv[0] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->NmrInv[1] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->NmrInv[2] = (double *) DR_Array (DOUBLE, 0, Area1);

        dev_data->NmrInvLag[0] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->NmrInvLag[1] = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->NmrInvLag[2] = (double *) DR_Array (DOUBLE, 0, Area1);

        dev_data->Shift1 = (int *) DR_Array (INT, 0, Area1);

        dev_data->pu = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->p0 = (double *) DR_Array (DOUBLE, 0, Area1);
        dev_data->pd = (double *) DR_Array (DOUBLE, 0, Area1);

        if (  (dev_data->NewPrice    == NULL)
           || (dev_data->Discount[0] == NULL)
           || (dev_data->Discount[1] == NULL)
           || (dev_data->Discount[2] == NULL)
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
           || (dev_data->NmrInvLag[0]== NULL)
           || (dev_data->NmrInvLag[1]== NULL)
           || (dev_data->NmrInvLag[2]== NULL)
           || (dev_data->Shift1      == NULL)
           || (dev_data->pu          == NULL)
           || (dev_data->p0          == NULL)
           || (dev_data->pd          == NULL))
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 2)
    {
        dev_data->NewPrice = (double *) DR_Array (DOUBLE, 0, Area2);

        dev_data->Discount[0] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->Discount[1] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->Discount[2] = (double *) DR_Array (DOUBLE, 0, Area2);

        dev_data->NmrInv[0] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->NmrInv[1] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->NmrInv[2] = (double *) DR_Array (DOUBLE, 0, Area2);

        dev_data->NmrInvLag[0] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->NmrInvLag[1] = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->NmrInvLag[2] = (double *) DR_Array (DOUBLE, 0, Area2);

        dev_data->Shift1 = (int *) DR_Array (INT, 0, Area1);
        dev_data->Shift2 = (int *) DR_Array (INT, 0, Area2);

        dev_data->quu = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qu0 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qud = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q0u = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q00 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q0d = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qdu = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qd0 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qdd = (double *) DR_Array (DOUBLE, 0, Area2);

        if (  (dev_data->NewPrice    == NULL)
           || (dev_data->Discount[0] == NULL)
           || (dev_data->Discount[1] == NULL)
           || (dev_data->Discount[2] == NULL)
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
           || (dev_data->NmrInvLag[0]== NULL)
           || (dev_data->NmrInvLag[1]== NULL)
           || (dev_data->NmrInvLag[2]== NULL)
           || (dev_data->Shift1      == NULL)
           || (dev_data->Shift2      == NULL)
           || (dev_data->quu         == NULL)
           || (dev_data->qu0         == NULL)
           || (dev_data->qud         == NULL)
           || (dev_data->q0u         == NULL)
           || (dev_data->q00         == NULL)
           || (dev_data->q0d         == NULL)
           || (dev_data->qdu         == NULL)
           || (dev_data->qd0         == NULL)
           || (dev_data->qdd         == NULL))
        {
            goto RETURN;
        }
    }
    else if (tree_data->NbFactor == 3)
    {
        dev_data->NewPrice = (double *) DR_Array (DOUBLE, 0, Area3);

        dev_data->Discount[0] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->Discount[1] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->Discount[2] = (double *) DR_Array (DOUBLE, 0, Area3);

        dev_data->NmrInv[0] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->NmrInv[1] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->NmrInv[2] = (double *) DR_Array (DOUBLE, 0, Area3);

        dev_data->NmrInvLag[0] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->NmrInvLag[1] = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->NmrInvLag[2] = (double *) DR_Array (DOUBLE, 0, Area3);

        dev_data->Shift1 = (int *) DR_Array (INT, 0, Area1);
        dev_data->Shift2 = (int *) DR_Array (INT, 0, Area2);
        dev_data->Shift3 = (int *) DR_Array (INT, 0, Area3);

        dev_data->quu = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qu0 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qud = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q0u = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q00 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->q0d = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qdu = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qd0 = (double *) DR_Array (DOUBLE, 0, Area2);
        dev_data->qdd = (double *) DR_Array (DOUBLE, 0, Area2);

        dev_data->ru = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->r0 = (double *) DR_Array (DOUBLE, 0, Area3);
        dev_data->rd = (double *) DR_Array (DOUBLE, 0, Area3);

        if (  (dev_data->NewPrice    == NULL)
           || (dev_data->Discount[0] == NULL)
           || (dev_data->Discount[1] == NULL)
           || (dev_data->Discount[2] == NULL)
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
           || (dev_data->NmrInvLag[0]== NULL)
           || (dev_data->NmrInvLag[1]== NULL)
           || (dev_data->NmrInvLag[2]== NULL)
           || (dev_data->Shift1      == NULL)
           || (dev_data->Shift2      == NULL)
           || (dev_data->Shift3      == NULL)
           || (dev_data->quu         == NULL)
           || (dev_data->qu0         == NULL)
           || (dev_data->qud         == NULL)
           || (dev_data->q0u         == NULL)
           || (dev_data->q00         == NULL)
           || (dev_data->q0d         == NULL)
           || (dev_data->qdu         == NULL)
           || (dev_data->qd0         == NULL)
           || (dev_data->qdd         == NULL)
           || (dev_data->ru          == NULL)
           || (dev_data->r0          == NULL)
           || (dev_data->rd          == NULL))
        {
            goto RETURN;
        }
    }  /* if then else */

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
    {
        DR_Error ("Fix3_Dev_Alloc: could not allocate memory for dev structure!");
    }

    return (status);

}  /* Fix3_Dev_Alloc */



/*****  Fix3_Dev_Free  ***********************************************************/
/*
*       Free memory for the dev structure.
*/
int     Fix3_Dev_Free (	FIX3_DEV_DATA    *dev_data,  /* (I) Fix3_Dev data structure  */
                    FIX3_TREE_DATA const* tree_data) /* (I) Tree data structure */
{
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];

    if (tree_data->NbFactor == 1)
    {
        Free_DR_Array (dev_data->NewPrice,    DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->Discount[0], DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->Discount[1], DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->Discount[2], DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInv[0],   DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInv[1],   DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInv[2],   DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInvLag[0],DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInvLag[1],DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->NmrInvLag[2],DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->Shift1,      INT,    0, Area1);
        Free_DR_Array (dev_data->pu,          DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->p0,          DOUBLE, 0, Area1);
        Free_DR_Array (dev_data->pd,          DOUBLE, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
        Free_DR_Array (dev_data->NewPrice,    DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->Discount[0], DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->Discount[1], DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->Discount[2], DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInv[0],   DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInv[1],   DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInv[2],   DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInvLag[0],DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInvLag[1],DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->NmrInvLag[2],DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->Shift1,      INT,    0, Area1);
        Free_DR_Array (dev_data->Shift2,      INT,    0, Area2);
        Free_DR_Array (dev_data->quu,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qu0,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qud,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q0u,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q00,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q0d,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qdu,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qd0,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qdd,         DOUBLE, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
        Free_DR_Array (dev_data->NewPrice,    DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->Discount[0], DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->Discount[1], DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->Discount[2], DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInv[0],   DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInv[1],   DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInv[2],   DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInvLag[0],DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInvLag[1],DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->NmrInvLag[2],DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->Shift1,      INT,    0, Area1);
        Free_DR_Array (dev_data->Shift2,      INT,    0, Area2);
        Free_DR_Array (dev_data->Shift3,      INT,    0, Area3);
        Free_DR_Array (dev_data->quu,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qu0,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qud,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q0u,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q00,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->q0d,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qdu,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qd0,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->qdd,         DOUBLE, 0, Area2);
        Free_DR_Array (dev_data->ru,          DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->r0,          DOUBLE, 0, Area3);
        Free_DR_Array (dev_data->rd,          DOUBLE, 0, Area3);

    }  /* if then else */

    Fix3_Dev_Init(dev_data);
    return (SUCCESS);

}  /* Fix3_Dev_Free */



/*****  Fix3_Alloc_Slice  ********************************************************/
/*
*       Allocate memory for one time slice.
*/
double  *Fix3_Alloc_Slice (FIX3_TREE_DATA   const* tree_data) /* (I) Tree data structure */
{
    void    *Slice;
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];


    Slice = (double *) NULL;

    if (tree_data->NbFactor == 1)
    {
         Slice = (double *) DR_Array (DOUBLE, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
         Slice = (double *) DR_Array (DOUBLE, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
         Slice = (double *) DR_Array (DOUBLE, 0, Area3);

    }  /* if then else */


    return (Slice);


}  /* Fix3_Alloc_Slice */



/*****  Fix3_Free_Slice  *********************************************************/
/*
*       Free memory for time slice.
*/
int     Fix3_Free_Slice (double      *Slice,     /* (I) Time slice          */
                    FIX3_TREE_DATA const* tree_data)	/* (I) Tree data structure */
{
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];

    if (tree_data->NbFactor == 1)
    {
        Free_DR_Array (Slice, DOUBLE, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
        Free_DR_Array (Slice, DOUBLE, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
        Free_DR_Array (Slice, DOUBLE, 0, Area3);

    }  /* if then else */

    return (SUCCESS);

}  /* Fix3_Free_Slice */


/*****  Fix3_Alloc_Slice_Int  *****************************************************/
/*
*       Allocate memory for one time slice.
*/
int     *Fix3_Alloc_Slice_Int (FIX3_TREE_DATA   const* tree_data) /* (I) Tree data structure */
{
    void    *Slice;
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];


    Slice = (int *) NULL;

    if (tree_data->NbFactor == 1)
    {
         Slice = (double *) DR_Array (INT, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
         Slice = (double *) DR_Array (INT, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
         Slice = (double *) DR_Array (INT, 0, Area3);

    }  /* if then else */


    return (Slice);


}  /* Fix3_Alloc_Slice_Int */



/*****  Fix3_Free_Slice_Int  ******************************************************/
/*
*       Free memory for time slice.
*/
int     Fix3_Free_Slice_Int (int        *Slice,     /* (I) Time slice          */
                    FIX3_TREE_DATA const* tree_data)/* (I) Tree data structure */
{
    int     Area1;              /* Surface of the ellipse */
    int     Area2;
    int     Area3;


    Area1 = tree_data->Width[0];
    Area2 = tree_data->Width[0] * tree_data->Width[1];
    Area3 = tree_data->Width[0] * tree_data->Width[1] * tree_data->Width[2];

    if (tree_data->NbFactor == 1)
    {
        Free_DR_Array (Slice, INT, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
        Free_DR_Array (Slice, INT, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
        Free_DR_Array (Slice, INT, 0, Area3);

    }  /* if then else */

    return (SUCCESS);

}  /* Fix3_Free_Slice_Int */

