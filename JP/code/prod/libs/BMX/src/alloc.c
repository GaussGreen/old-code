/****************************************************************************/
/*        Memory allocation.                                                */
/****************************************************************************/
/*        ALLOC.c                                                           */
/****************************************************************************/

/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "bmx123head.h"


/*****  Tree_Init  **********************************************************/
/*
*       Initialize tree pointers to NULL.
*/
void    Tree_Init (TREE_DATA    *tree_data) /* Tree building data structure */
{
    int i;

    tree_data->NbTP = 0;

    for (i = 0; i < NBCRITDATE; i++)
    {
        tree_data->CritDate[i] = NULL;
        tree_data->TPtype[i]   = NULL;
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

    for (i = 0; i < 3; i++)
    {
        tree_data->ZeroCoupon[i] = NULL;
        tree_data->ZeroRate[i]   = NULL;
        tree_data->TermZero[i]   = NULL;
        tree_data->FwdRate[i]    = NULL;
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

    tree_data->NbNmr   = 0;
    tree_data->NmrDate = NULL;
    tree_data->NmrInv  = NULL;

    tree_data->Ann1  = NULL;
    tree_data->Ann0  = NULL;
    tree_data->Zero1 = NULL;
    tree_data->Zero0 = NULL;

    tree_data->YLMin = NULL;
    tree_data->YLMax = NULL;
    tree_data->Yield = NULL;
    tree_data->ZAInv = NULL;
    tree_data->YCumP = NULL;

    tree_data->JumpPpy = 0;

    tree_data->XValues = NULL;
    tree_data->YValues = NULL;

    return;

}  /* Tree_Init */



/*****  Tree_Alloc  *********************************************************/
/*
*       Allocation of memory for the arrays constituting the tree.
*       The limits of the tree are allocated directly in Tree_Limits.
*/
int     Tree_Alloc (TREE_DATA   *tree_data) /* Tree building data structure */
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
        DR_Error ("Tree_Alloc: could not allocate memory for the tree!");
    }

    return (status);

}  /* Tree_Alloc */



/*****  Tree_Free  **********************************************************/
/*
*       Free the arrays constituting the tree.
*/
int     Tree_Free ( TREE_DATA   *tree_data) /* Tree building data structure */
{
    int     i;
    int     t;
    long    Area = 0;
    int     NbTP;               /* Number of time points in the tree */
    long    EDevIdx, NbNodes;


    /* Slices size */
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


    NbTP = tree_data->NbTP;


    if (tree_data->Yield != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_DR_Array(tree_data->Yield[i],DOUBLE,-1,Area);
        }

        Free_DR_Array(tree_data->Yield,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }


    if (tree_data->ZAInv != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_DR_Array(tree_data->ZAInv[i],DOUBLE,-1,Area);
        }

        Free_DR_Array(tree_data->ZAInv,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    if (tree_data->YCumP != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_DR_Array(tree_data->YCumP[i],DOUBLE,-1,Area);
        }

        Free_DR_Array(tree_data->YCumP,DOUBLE_PTR, 0,tree_data->NbNmr-1);
    }
    


    /* First of all, free slices living in TREE_DATA */
    if (tree_data->EDevDate != NULL)
    {
        Free_DR_Array (tree_data->EDevDate,LONG,0, 
                       tree_data->NbEDevDates-1);
    }

    if (tree_data->EDevStPrice != NULL)
    {
        for (i = 0; i < tree_data->NbEDevDates; i++)
        {
            Free_Slice(tree_data->EDevStPrice[i],tree_data);
        }

        Free_DR_Array(tree_data->EDevStPrice,DOUBLE_PTR,0,
                      tree_data->NbEDevDates-1);
    }

    if (tree_data->XValues != NULL)
    {
        for (i = 0; i < tree_data->NbEDevDates; i++)
        {
            Free_Slice(tree_data->XValues[i],tree_data);
        }

        Free_DR_Array(tree_data->XValues,DOUBLE_PTR,0,
                      tree_data->NbEDevDates-1);
    }

    if (tree_data->YValues != NULL)
    {
        for (i = 0; i < tree_data->NbEDevDates; i++)
        {
            Free_Slice(tree_data->YValues[i],tree_data);
        }

        Free_DR_Array(tree_data->YValues,DOUBLE_PTR,0,
                      tree_data->NbEDevDates-1);
    } 



    if (tree_data->NmrDate != NULL)
    {
        Free_DR_Array (tree_data->NmrDate,LONG,0,tree_data->NbNmr-1);
    }

    if (tree_data->NmrInv != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_Slice(tree_data->NmrInv[i],tree_data);
        }

        Free_DR_Array (tree_data->NmrInv,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    if (tree_data->Ann1 != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_Slice(tree_data->Ann1[i], tree_data);
        }

        Free_DR_Array(tree_data->Ann1,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    if (tree_data->Ann0 != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_Slice(tree_data->Ann0[i],tree_data);
        }

        Free_DR_Array(tree_data->Ann0,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    if (tree_data->Zero1 != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_Slice(tree_data->Zero1[i], tree_data);
        }

        Free_DR_Array(tree_data->Zero1,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    if (tree_data->Zero0 != NULL)
    {
        for (i = 0; i < tree_data->NbNmr; i++)
        {
            Free_Slice(tree_data->Zero0[i], tree_data);
        }

        Free_DR_Array(tree_data->Zero0,DOUBLE_PTR,0,tree_data->NbNmr-1);
    }
    

    
    if (tree_data->YLMin != NULL)
    {
        Free_DR_Array (tree_data->YLMin,INT,0,tree_data->NbNmr-1);
    }

    if (tree_data->YLMax != NULL)
    {
        Free_DR_Array (tree_data->YLMax,INT,0,tree_data->NbNmr-1);
    }


    if (tree_data->MappSize != NULL)
    {
        Free_DR_Array (tree_data->MappSize,LONG,0,tree_data->NbNmr-1);
    }

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

    
    

    return (SUCCESS);

}  /* Tree_Free */



/*****  Dev_Init  **********************************************************/
/*
*       Initialize DEV_DATA pointers to NULL.
*/
void    Dev_Init (DEV_DATA    *dev_data)
{
    int   i;

    dev_data->Shift1 = dev_data->Shift2 = dev_data->Shift3 = NULL;

    for (i=0; i<3; i++) dev_data->NmrInv[i]    = NULL;
    for (i=0; i<3; i++) dev_data->NmrInvLag[i] = NULL;

    dev_data->pu  = dev_data->p0  = dev_data->pd  = NULL;
    dev_data->quu = dev_data->qu0 = dev_data->qud = NULL;
    dev_data->q0u = dev_data->q00 = dev_data->q0d = NULL;
    dev_data->qdu = dev_data->qd0 = dev_data->qdd = NULL;
    dev_data->ru  = dev_data->r0  = dev_data->rd  = NULL;

    dev_data->NewPrice = NULL;

    /* Initialize numeraire flags */
    dev_data->NmrToCcy = FALSE;
    dev_data->CcyToNmr = FALSE;

    return;

}  /* Dev_Init */


/*****  Dev_Alloc  **********************************************************/
/*
*       Allocation of memory for the dev structure. We allocate arrays and
*       the addressing will be done linearly.
*/
int     Dev_Alloc ( DEV_DATA    *dev_data,  /* (I) Dev data structure  */
                    TREE_DATA   *tree_data) /* (I) Tree data structure */
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
        dev_data->NewPrice  = (double *) DR_Array (DOUBLE, 0, Area1);

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
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
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
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
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
           || (dev_data->NmrInv[0]   == NULL)
           || (dev_data->NmrInv[1]   == NULL)
           || (dev_data->NmrInv[2]   == NULL)
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
        DR_Error ("Dev_Alloc: could not allocate memory for dev structure!");
    }

    return (status);

}  /* Dev_Alloc */



/*****  Dev_Free  ***********************************************************/
/*
*       Free memory for the dev structure.
*/
int     Dev_Free (  DEV_DATA    *dev_data,  /* (I) Dev data structure  */
                    TREE_DATA   *tree_data) /* (I) Tree data structure */
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

    return (SUCCESS);

}  /* Dev_Free */



/*****  MktVol_Init  **********************************************************/
/*
*       Initialize MKTVOL_DATA structure.
*/
void    MktVol_Init (MKTVOL_DATA    *mktvol_data)
{
    mktvol_data->FilterSpotVolFlag = FALSE;
    mktvol_data->SmoothingFlag     = 'N';
    mktvol_data->NmrStatFlag       = TRUE;
    mktvol_data->Trace             = 1;
    return;

}  /* MktVol_Init */



/*****  Opt_Out_Data_Init  **************************************************/
/*
 *      Initialize the OPT_OUT_DATA structure
 */
void     Opt_Out_Data_Init(OPT_OUT_DATA  *ood)
{
    int  i;

    if (ood == NULL) return;

    ood->Option = 0.0;
    for (i=0; i<10; i++) ood->Price[i] = 0.0;

    return;

} /* Opt_Out_Data_Init */



/*****  Alloc_Slice  ********************************************************/
/*
*       Allocate memory for one time slice.
*/
double  *Alloc_Slice (  TREE_DATA   *tree_data) /* (I) Tree data structure */
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


}  /* Alloc_Slice */



/*****  Free_Slice  *********************************************************/
/*
*       Free memory for time slice.
*/
int     Free_Slice (double      *Slice,     /* (I) Time slice          */
                    TREE_DATA   *tree_data) /* (I) Tree data structure */
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

}  /* Free_Slice */



/*****  DR_Array  ***********************************************************/
/*
*       Allocation of memory for an array of arbitrary type. Returns void*.
*/
void    *DR_Array ( int     type,   /* (I) Type         */
                    int     nl,     /* (I) Lower bound  */
                    int     nh)     /* (I) Higher bound */
{

    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v;

            v = (CRIT_DATE *) calloc ((unsigned) (nh-nl+1), sizeof (CRIT_DATE));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            /* Shift pointer and cast to void* */

            return ((void *) (v-nl));

        }                          
        case INT:
        {
            int *v;

            v = (int *) calloc ((unsigned) (nh-nl+1), sizeof (int));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_PTR:
        {
            int **v;

            v = (int **) calloc ((unsigned) (nh-nl+1), sizeof (int *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case INT_D_PTR:
        {
            int ***v;

            v = (int ***) calloc ((unsigned) (nh-nl+1), sizeof (int **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG:
        {
            long *v;

            v = (long *) calloc ((unsigned) (nh-nl+1), sizeof (long));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_PTR:
        {
            long **v;

            v = (long **) calloc ((unsigned) (nh-nl+1), sizeof (long *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case LONG_D_PTR:
        {
            long ***v;

            v = (long ***) calloc ((unsigned) (nh-nl+1), sizeof (long **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE:
        {
            double *v;

            v = (double *) calloc ((unsigned) (nh-nl+1), sizeof (double));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_PTR:
        {
            double **v;

            v = (double **) calloc ((unsigned) (nh-nl+1), sizeof (double *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }                          
        case DOUBLE_D_PTR:
        {
            double ***v;

            v = (double ***) calloc ((unsigned) (nh-nl+1), sizeof (double **));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }
        case CHAR:
        {
            char *v;

            v = (char *) calloc ((unsigned) (nh-nl+1), sizeof (char));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        } 
        case CHAR_PTR:
        {
            char **v;

            v = (char **) calloc ((unsigned) (nh-nl+1), sizeof (char *));

            if (v == NULL)
            {
                DR_Error ("DR_Array: allocation failure!");
                return (NULL);
            }

            return ((void *) (v-nl));

        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in DR_Array!");
            return (NULL);
        }                          
    }  /* switch */

}  /* DR_Array */


/*****  DR_Matrix  **********************************************************/
/*
 *      Allocation of memory for a matrix of arbitrary type. Returns void *.
 */
void  *DR_Matrix (int     type,   /* (I) Type                */
                  int     nrl,    /* (I) Lower row index     */
                  int     nrh,    /* (I) Higher row index    */
                  int     ncl,    /* (I) Lower column index  */
                  int     nch)    /* (I) Higher column index */
{
    void  **Matx = NULL;
    int     i;
    int     type_ptr;

    switch (type)
    {
    case INT:
        type_ptr = INT_PTR;
        break;
    case LONG:
        type_ptr = LONG_PTR;
        break;
    case DOUBLE:
        type_ptr = DOUBLE_PTR;
        break;
    default:
        DR_Error ("DR_Matrix: unsupported type.");
        goto RETURN;
    }  /* switch */

    Matx = (void **) DR_Array (type_ptr, nrl, nrh);

    if (Matx == NULL)
    {
        DR_Error ("DR_Matrix: allocation failure #1!");
        goto RETURN;
    }

    for (i = nrl; i <= nrh; i++)
    {
        Matx[i] = DR_Array (type, ncl, nch);

        if (Matx[i] == NULL)
        {
            DR_Error ("DR_Matrix: allocation failure #2!");
            goto RETURN;
        }
    }

    return ((void *) Matx);

RETURN:

    /* allocation failed, release memory alloc'd so far */
    Free_DR_Matrix (Matx, type, nrl, nrh, ncl, nch);

    return (NULL);

}  /* DR_Matrix */


/*****  DR_Cube  ****************************************************************/
/*
 *       Allocation of memory for a 3D array of arbitrary type. Returns void *.
 */
void    *DR_Cube (  int     type,   /* (I) Type                */
                    int     nl,     /* (I) Lower bound         */
                    int     nh,     /* (I) Higher bound        */
                    int     nrl,    /* (I) Lower row index     */
                    int     nrh,    /* (I) Higher row index    */
                    int     ncl,    /* (I) Lower column index  */
                    int     nch)    /* (I) Higher column index */
{
    void ***Cube = NULL;
    int     i, j;
    int     type_ptr, type_dbl_ptr;

    switch (type)
    {
    case INT:
        type_dbl_ptr = INT_D_PTR;
        type_ptr     = INT_PTR;
        break;
    case LONG:
        type_dbl_ptr = LONG_D_PTR;
        type_ptr     = LONG_PTR;
        break;
    case DOUBLE:
        type_dbl_ptr = DOUBLE_D_PTR;
        type_ptr     = DOUBLE_PTR;
        break;
    default:
        DR_Error ("DR_Cube: unsupported type.");
        goto RETURN;
    }  /* switch */

    Cube = (void ***) DR_Array (type_dbl_ptr, nl, nh);

    if (Cube == NULL)
    {
        DR_Error ("DR_Cube: allocation failure #1!");
        goto RETURN;
    }

    for (i = nl; i <= nh; i++)
    {
        Cube[i] = (void **) DR_Array (type_ptr, nrl, nrh);
        if (Cube[i] == NULL)
        {
            DR_Error ("DR_Cube: allocation failure #2!");
            goto RETURN;
        }
        
        for (j = nrl; j <= nrh; j++)
        {
            Cube[i][j] = DR_Array (type, ncl, nch);
            
            if (Cube[i][j] == NULL)
            {
                DR_Error ("DR_Cube: allocation failure #3!");
                goto RETURN;
            }
        } /* j */
    } /* i */

    return ((void *) Cube);

RETURN:

    /* allocation failed, release memory alloc'd so far */
    Free_DR_Cube (Cube, type, nl, nh, nrl, nrh, ncl, nch);

    return (NULL);

}  /* DR_Cube */


/*****  Free_DR_Array  ******************************************************/
/*
*       Free DR array.
*/
int     Free_DR_Array ( void   *Array,  /* (I) Array        */
                        int     type,   /* (I) Type         */
                        int     nl,     /* (I) Lower bound  */
                        int     nh)     /* (I) Higher bound */
{
    /* To avoid warning message */
    nh += 0;

    switch (type)
    {
        case CRITDATE:
        {
            CRIT_DATE *v = (CRIT_DATE *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT:
        {
            int *v = (int *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_PTR:
        {
            int **v = (int **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case INT_D_PTR:
        {
            int ***v = (int ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG:
        {
            long *v = (long *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_PTR:
        {
            long **v = (long **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case LONG_D_PTR:
        {
            long ***v = (long ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double *v = (double *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_PTR:
        {
            double **v = (double **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case DOUBLE_D_PTR:
        {
            double ***v = (double ***) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }
        case CHAR:
        {
            char *v = (char *) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }                          
        case CHAR_PTR:
        {
            char **v = (char **) Array;

            if (v != NULL)
            {
                free ((char *) (v+nl));
            }

            break;
        }  
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Array!");
            return (FAILURE);
        }                          
    }  /* switch */

    Array = NULL;

    return (SUCCESS);

}  /* Free_DR_Array */



/*****  Free_DR_Matrix  *****************************************************/
/*
*       Free DR matrix.
*/
int     Free_DR_Matrix (void   *Matrix,     /* (I) Matrix              */
                        int     type,       /* (I) Type                */
                        int     nrl,        /* (I) Lower row index     */
                        int     nrh,        /* (I) Higher row index    */
                        int     ncl,        /* (I) Lower column index  */
                        int     nch)        /* (I) Higher column index */
{
    int i;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int **m = (int **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case LONG:
        {
            long **m = (long **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        case DOUBLE:
        {
            double **m = (double **) Matrix;

            if (m != NULL)
            {
                for (i = nrh; i >= nrl; i--)
                {
                    if (m[i] != NULL)
                    {
                        free ((char *) (m[i] + ncl));
                    }
                }

                free ((char *) (m+nrl));
            }

            break;
        }                          
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Matrix!");
            return (FAILURE);
        }                          
    }  /* switch */

    Matrix = NULL;

    return (SUCCESS);

}  /* Free_DR_Matrix */



/*****  Free_DR_Cube  *******************************************************/
/*
*       Free 3D DR array.
*/
int     Free_DR_Cube (  void   *Cube,   /* (I) Cube                */
                        int     type,   /* (I) Type                */
                        int     nl,     /* (I) Lower bound         */
                        int     nh,     /* (I) Higher bound        */
                        int     nrl,    /* (I) Lower row index     */
                        int     nrh,    /* (I) Higher row index    */
                        int     ncl,    /* (I) Lower column index  */
                        int     nch)    /* (I) Higher column index */
{
    int i;
    int j;

    /* To avoid warning message */
    nch += 0;

    switch (type)
    {
        case INT:
        {
            int ***c = (int ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case LONG:
        {
            long ***c = (long ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                          
        case DOUBLE:
        {
            double ***c = (double ***) Cube;

            if (c != NULL)
            {
                for (i = nh; i >= nl; i--)
                {
                    if (c[i] != NULL)
                    {
                        for (j = nrh; j >= nrl; j--)
                        {
                            if (c[i][j] != NULL)
                            {
                                free ((char *) (c[i][j] + ncl));
                            }
                        }

                        free ((char *) (c[i] + nrl));
                    }
                }

                free ((char *) (c+nl));
            }

            break;

        }                        
        default:
        {
            DR_Error ("Program bug: unsupported type in Free_DR_Cube!");
            return (FAILURE);
        }                          
    }  /* switch */

    Cube = NULL;

    return (SUCCESS);

}  /* Free_DR_Cube */




/*****  Alloc_Slice_Int  *****************************************************/
/*
*       Allocate integer memory for one time slice.
*/
int  *Alloc_Slice_Int (  TREE_DATA   *tree_data) /* (I)Tree data structure*/
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
         Slice = (int *) DR_Array (INT, 0, Area1);
    }
    else if (tree_data->NbFactor == 2)
    {
         Slice = (int *) DR_Array (INT, 0, Area2);
    }
    else if (tree_data->NbFactor == 3)
    {
         Slice = (int *) DR_Array (INT, 0, Area3);

    }  /* if then else */

    return (Slice);

}  /* Alloc_Slice */







/*****  Free_Slice_Int  *****************************************************/
/*
*       Free integer memory for time slice.
*/
int     Free_Slice_Int (int         *Slice,     /* (I) Time slice          */
                        TREE_DATA   *tree_data) /* (I) Tree data structure */

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

}  /* Free_Slice */

