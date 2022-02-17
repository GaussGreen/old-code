/****************************************************************************/
/*        Memory allocation.                                                */
/****************************************************************************/
/*        HYB4_ALLOC.c                                                      */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include "hyb4_lib.h"

void Hyb4_Tree_Init(HYB4_TREE_DATA *tree_data)
{
    int i;

    tree_data->NodeOffset0 = NULL;
    tree_data->NodeOffset1 = NULL;
    tree_data->NodeOffset2 = NULL;
    tree_data->NodeOffset3 = NULL;

    tree_data->NbAsset =  0;
    tree_data->NbTP    = -1;

    tree_data->LatticeTime = 0;

    for (i=0; i<4; ++i)
    {
        tree_data->AssetType  [i] = -1;

        tree_data->asset      [i] = NULL;

        tree_data->Width      [i] = 0;
        tree_data->HalfWidth  [i] = 0;

        tree_data->xWidth     [i] = 0;
        tree_data->xHalfWidth [i] = 0;

        tree_data->MaxIndex   [i] = 0;
    }

    tree_data->xDate = 0;
    tree_data->xT    = 0;

    tree_data->NbThread = 1;

    tree_data->iMin = NULL;
    tree_data->iMax = NULL;
    tree_data->jMin = NULL;
    tree_data->jMax = NULL;
    tree_data->kMin = NULL;
    tree_data->kMax = NULL;
    tree_data->LMin = NULL;
    tree_data->LMax = NULL;

    tree_data->iOutMin = NULL;
    tree_data->iOutMax = NULL;
    tree_data->jOutMin = NULL;
    tree_data->jOutMax = NULL;
    tree_data->kOutMin = NULL;
    tree_data->kOutMax = NULL;
    tree_data->LOutMin = NULL;
    tree_data->LOutMax = NULL;

    tree_data->Drift0 = NULL;
    tree_data->DriftF = NULL;
    tree_data->Drift1 = NULL;
    tree_data->Drift2 = NULL;
    tree_data->Drift3 = NULL;

    tree_data->Length    = NULL;
    tree_data->LengthJ   = NULL;
    tree_data->JumpCoeff = NULL;

    tree_data->Time      = NULL;

    tree_data->J00 = NULL;
    tree_data->J10 = NULL;
    tree_data->J11 = NULL;
    tree_data->J20 = NULL;
    tree_data->J21 = NULL;
    tree_data->J22 = NULL;
    tree_data->J30 = NULL;
    tree_data->J31 = NULL;
    tree_data->J32 = NULL;
    tree_data->J33 = NULL;

    tree_data->R10 = NULL;
    tree_data->R20 = NULL;
    tree_data->R21 = NULL;
    tree_data->R30 = NULL;
    tree_data->R31 = NULL;
    tree_data->R32 = NULL;

    tree_data->M0 = NULL;
    tree_data->M1 = NULL;
    tree_data->M2 = NULL;
    tree_data->M3 = NULL;

    tree_data->V00 = NULL;
    tree_data->V10 = NULL;
    tree_data->V11 = NULL;
    tree_data->V20 = NULL;
    tree_data->V21 = NULL;
    tree_data->V22 = NULL;
    tree_data->V30 = NULL;
    tree_data->V31 = NULL;
    tree_data->V32 = NULL;
    tree_data->V33 = NULL;

    tree_data->StateProb = NULL;

    tree_data->LimLoMin = NULL;
    tree_data->LimLoMax = NULL;
    tree_data->LimHiMin = NULL;
    tree_data->LimHiMax = NULL;
}

int Hyb4_Tree_Limits_Assign(HYB4_TREE_DATA *tree_data_new,
                            HYB3_TREE_DATA *tree_data_old)
{
    int NbTP = tree_data_new->NbTP;

    int NbAssetOn = tree_data_new->NbAssetOn;

    int t, xT;

    int status = FAILURE;

    tree_data_new->TPDate = tree_data_old->TPDate;
    tree_data_new->xDate  = tree_data_old->xDate;
    tree_data_new->xT     = tree_data_old->xT;

    xT = tree_data_new->xT;

    for (t=0; t<=NbTP; ++t)
    {
        if (NbAssetOn > 0)
        {
            tree_data_new->iMin[t] = tree_data_old->Bottom1[t];
            tree_data_new->iMax[t] = tree_data_old->Top1[t];
        }

        if (NbAssetOn > 1)
        {
            tree_data_new->jMin[t] = tree_data_old->Bottom2[t];
            tree_data_new->jMax[t] = tree_data_old->Top2[t];
        }

        if (NbAssetOn > 2)
        {
            tree_data_new->kMin[t] = tree_data_old->Bottom3[t];
            tree_data_new->kMax[t] = tree_data_old->Top3[t];
        }

        if ((NbAssetOn > 3) && (t <= xT))
        {
            tree_data_new->LMin[t] = tree_data_old->Bottom4[t];
            tree_data_new->LMax[t] = tree_data_old->Top4[t];
        }
    }

    status = SUCCESS;

    return status;
}

int Hyb4_Tree_Limits_Unassign(HYB4_TREE_DATA *tree_data_new)
{
    int NbTP = tree_data_new->NbTP;

    int NbAssetOn = tree_data_new->NbAssetOn;

    int t;
    int xT = tree_data_new->xT;

    int status = FAILURE;

    for (t=0; t<=NbTP; ++t)
    {
        if (NbAssetOn > 1)
        {
            tree_data_new->jMin[t] = NULL;
            tree_data_new->jMax[t] = NULL;
        }

        if (NbAssetOn > 2)
        {
            tree_data_new->kMin[t] = NULL;
            tree_data_new->kMax[t] = NULL;
        }

        if ((NbAssetOn > 3) && (t <= xT))
        {
            tree_data_new->LMin[t] = NULL;
            tree_data_new->LMax[t] = NULL;
        }
    }

    tree_data_new->TPDate = NULL;

    status = SUCCESS;

    return status;
}

int Hyb4_Offset_Alloc(HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    int NbAssetOn = tree_data->NbAssetOn;

    int i, j;
    int t;
    int xT = tree_data->xT;

    int iMin = 0;
    int iMax = 0;

    int *jMin = NULL;
    int *jMax = NULL;

    int jMini = 0;
    int jMaxi = 0;

    int **kMin = NULL;
    int **kMax = NULL;

    int *kMini = NULL;
    int *kMaxi = NULL;
    
    int kMinij = 0;
    int kMaxij = 0;

    int status = FAILURE;

    for (t=0; t<=NbTP; ++t)
    {
        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        if (NbAssetOn > 1)
        {
            jMin = tree_data->jMin[t];
            jMax = tree_data->jMax[t];

            tree_data->NodeOffset1[t] = (int *)DR_Array(INT, iMin, iMax);

            if (NbAssetOn > 2)
            {
                kMin = tree_data->kMin[t];
                kMax = tree_data->kMax[t];

                tree_data->NodeOffset2[t] = (int **)DR_Array(INT_PTR, iMin, iMax);

                if ((NbAssetOn > 3) && (t <= xT))
                {
                    tree_data->NodeOffset3[t] = (int ***)DR_Array(INT_D_PTR, iMin, iMax);
                }
            }

            for (i=iMin; i<=iMax; ++i)
            {
                jMini = jMin[i];
                jMaxi = jMax[i];

                if (NbAssetOn > 2)
                {
                    kMini = kMin[i];
                    kMaxi = kMax[i];

                    tree_data->NodeOffset2[t][i] = (int *)DR_Array(INT, jMini, jMaxi);

                    if ((NbAssetOn > 3) && (t <= xT))
                    {
                        tree_data->NodeOffset3[t][i] = (int **)DR_Array(INT_PTR, jMini, jMaxi);
                    }

                    for (j=jMini; j<=jMaxi; ++j)
                    {
                        kMinij = kMini[j];
                        kMaxij = kMaxi[j];

                        if ((NbAssetOn > 3) && (t <= xT))
                        {
                            tree_data->NodeOffset3[t][i][j] = (int *)DR_Array(INT, kMinij, kMaxij);
                        }
                    }
                }
            }
        }
    }

    status = SUCCESS;

    return status;
}
    
int Hyb4_Offset_Free(HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    int NbAssetOn = tree_data->NbAssetOn;

    int i, j;
    int t;
    int xT = tree_data->xT;

    int iMin = 0;
    int iMax = 0;

    int *jMin = NULL;
    int *jMax = NULL;

    int jMini = 0;
    int jMaxi = 0;

    int **kMin = NULL;
    int **kMax = NULL;

    int *kMini = NULL;
    int *kMaxi = NULL;
    
    int kMinij = 0;
    int kMaxij = 0;

    int status = FAILURE;

    for (t=0; t<=NbTP; ++t)
    {
        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        if (NbAssetOn > 1)
        {
            jMin = tree_data->jMin[t];
            jMax = tree_data->jMax[t];

            if (NbAssetOn > 2)
            {
                kMin = tree_data->kMin[t];
                kMax = tree_data->kMax[t];
            }

            for (i=iMin; i<=iMax; ++i)
            {
                jMini = jMin[i];
                jMaxi = jMax[i];

                if (NbAssetOn > 2)
                {
                    kMini = kMin[i];
                    kMaxi = kMax[i];

                    for (j=jMini; j<=jMaxi; ++j)
                    {
                        kMinij = kMini[j];
                        kMaxij = kMaxi[j];

                        if ((NbAssetOn > 3) && (t <= xT))
                        {
                            Free_DR_Array(tree_data->NodeOffset3[t][i][j], INT, kMinij, kMaxij);
                        }
                    }

                    if ((NbAssetOn > 3) && (t <= xT))
                    {
                        Free_DR_Array(tree_data->NodeOffset3[t][i], INT_PTR, jMini, jMaxi);
                    }

                    Free_DR_Array(tree_data->NodeOffset2[t][i], INT, jMini, jMaxi);
                }
            }

            if (NbAssetOn > 2)
            {
                if ((NbAssetOn > 3) && (t <= xT))
                {
                    Free_DR_Array(tree_data->NodeOffset3[t], INT_D_PTR, iMin, iMax);
                }

                Free_DR_Array(tree_data->NodeOffset2[t], INT_PTR, iMin, iMax);
            }

            Free_DR_Array(tree_data->NodeOffset1[t], INT, iMin, iMax);
        }
    }

    status = SUCCESS;

    return status;
}
        
int Hyb4_Tree_Alloc(HYB4_TREE_DATA *tree_data_new,
                    HYB3_TREE_DATA *tree_data_old)
{
    int NbTP = tree_data_old->NbTP;

    int h0   = tree_data_old->HalfWidth[0];
    int h1   = tree_data_old->HalfWidth[1];
    int h2   = tree_data_old->HalfWidth[2];
    int h3   = tree_data_old->HalfWidth[3];

    int w0   = tree_data_old->Width[0];
    int w1   = tree_data_old->Width[1];
    int w2   = tree_data_old->Width[2];
    int w3   = tree_data_old->Width[3];

    int xh0  = tree_data_old->xHalfWidth[0];
    int xh1  = tree_data_old->xHalfWidth[1];
    int xh2  = tree_data_old->xHalfWidth[2];
    int xh3  = tree_data_old->xHalfWidth[3];

    int xw0  = tree_data_old->xWidth[0];
    int xw1  = tree_data_old->xWidth[1];
    int xw2  = tree_data_old->xWidth[2];
    int xw3  = tree_data_old->xWidth[3];

    int status = FAILURE;

    tree_data_new->NbTP = NbTP;

    tree_data_new->HalfWidth[0] = h0;
    tree_data_new->HalfWidth[1] = h1;
    tree_data_new->HalfWidth[2] = h2;
    tree_data_new->HalfWidth[3] = h3;

    tree_data_new->Width[0] = w0;
    tree_data_new->Width[1] = w1;
    tree_data_new->Width[2] = w2;
    tree_data_new->Width[3] = w3;

    tree_data_new->xHalfWidth[0] = xh0;
    tree_data_new->xHalfWidth[1] = xh1;
    tree_data_new->xHalfWidth[2] = xh2;
    tree_data_new->xHalfWidth[3] = xh3;

    tree_data_new->xWidth[0] = xw0;
    tree_data_new->xWidth[1] = xw1;
    tree_data_new->xWidth[2] = xw2;
    tree_data_new->xWidth[3] = xw3;

    tree_data_new->NodeOffset0 = (int    *) DR_Array(INT,       -1, NbTP+1);
    tree_data_new->NodeOffset1 = (int   **) DR_Array(INT_PTR,   -1, NbTP+1);
    tree_data_new->NodeOffset2 = (int  ***) DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data_new->NodeOffset3 = (int ****) DR_Array(INT_T_PTR, -1, NbTP+1);

    tree_data_new->iMin      = (int *)   DR_Array(INT,       -1, NbTP+1);
    tree_data_new->iMax      = (int *)   DR_Array(INT,       -1, NbTP+1);
        
    tree_data_new->jMin      = (int **)  DR_Array(INT_PTR,   -1, NbTP+1);
    tree_data_new->jMax      = (int **)  DR_Array(INT_PTR,   -1, NbTP+1);
        
    tree_data_new->kMin      = (int ***) DR_Array(INT_D_PTR, -1, NbTP+1);
    tree_data_new->kMax      = (int ***) DR_Array(INT_D_PTR, -1, NbTP+1);
        
    tree_data_new->LMin      = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);
    tree_data_new->LMax      = (int ****)DR_Array(INT_T_PTR, -1, NbTP+1);
        
    tree_data_new->iOutMin   = tree_data_new->iMin;
    tree_data_new->iOutMax   = tree_data_new->iMax;

    tree_data_new->jOutMin   = tree_data_new->jMin;
    tree_data_new->jOutMax   = tree_data_new->jMax;
        
    tree_data_new->kOutMin   = tree_data_new->kMin;
    tree_data_new->kOutMax   = tree_data_new->kMax;
        
    tree_data_new->LOutMin   = tree_data_new->LMin;
    tree_data_new->LOutMax   = tree_data_new->LMax;
    
    tree_data_new->Drift0    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->DriftF    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->Drift1    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->Drift2    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->Drift3    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->Length    = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->LengthJ   = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->JumpCoeff = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->Time      = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->J00       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J10       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J11       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J20       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J21       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J22       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J30       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J31       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J32       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->J33       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->R10       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->R20       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->R21       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->R30       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->R31       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->R32       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->M0        = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->M1        = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->M2        = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->M3        = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->V00       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V10       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V11       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V20       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V21       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V22       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V30       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V31       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V32       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->V33       = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->StateProb = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    tree_data_new->LimLoMin  = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->LimLoMax  = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->LimHiMin  = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);
    tree_data_new->LimHiMax  = (double *)   DR_Array(DOUBLE,    -1, NbTP+1);

    status = SUCCESS;

    return status;
}

int Hyb4_Tree_Free(HYB4_TREE_DATA *tree_data)
{
    int    NbTP = tree_data->NbTP;

    int    status = FAILURE;

    Free_DR_Array(tree_data->NodeOffset0, INT,       -1, NbTP+1);
    Free_DR_Array(tree_data->NodeOffset1, INT_PTR,   -1, NbTP+1);
    Free_DR_Array(tree_data->NodeOffset2, INT_D_PTR, -1, NbTP+1);
    Free_DR_Array(tree_data->NodeOffset3, INT_T_PTR, -1, NbTP+1);

    tree_data->iOutMin = NULL;
    tree_data->iOutMax = NULL;

    tree_data->jOutMin = NULL;
    tree_data->jOutMax = NULL;

    tree_data->kOutMin = NULL;
    tree_data->kOutMax = NULL;

    tree_data->LOutMin = NULL;
    tree_data->LOutMax = NULL;

    Free_DR_Array(tree_data->iMin, INT,       -1, NbTP+1);
    Free_DR_Array(tree_data->iMax, INT,       -1, NbTP+1);
        
    Free_DR_Array(tree_data->jMin, INT_PTR,   -1, NbTP+1);
    Free_DR_Array(tree_data->jMax, INT_PTR,   -1, NbTP+1);
        
    Free_DR_Array(tree_data->kMin, INT_D_PTR, -1, NbTP+1);
    Free_DR_Array(tree_data->kMax, INT_D_PTR, -1, NbTP+1);
        
    Free_DR_Array(tree_data->LMin, INT_T_PTR, -1, NbTP+1);
    Free_DR_Array(tree_data->LMax, INT_T_PTR, -1, NbTP+1);
        
    Free_DR_Array(tree_data->Drift0, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->DriftF, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->Drift1, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->Drift2, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->Drift3, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->Length,    DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->LengthJ,   DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->JumpCoeff, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->Time,      DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->J00, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J10, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J11, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J20, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J21, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J22, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J30, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J31, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J32, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->J33, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->R10, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->R20, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->R21, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->R30, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->R31, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->R32, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->M0,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->M1,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->M2,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->M3,  DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->V00, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V10, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V11, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V20, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V21, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V22, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V30, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V31, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V32, DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->V33, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->StateProb, DOUBLE,    -1, NbTP+1);

    Free_DR_Array(tree_data->LimLoMin,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->LimLoMax,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->LimHiMin,  DOUBLE,    -1, NbTP+1);
    Free_DR_Array(tree_data->LimHiMax,  DOUBLE,    -1, NbTP+1);

    status = SUCCESS;

    return status;
}

void Hyb4_Dev_Init(HYB4_DEV_DATA *dev_data)
{
    int i;
    
    for (i=0; i<2; ++i)
    {
        dev_data->AssetValue0[i] = NULL;
        dev_data->AssetValue1[i] = NULL;
        dev_data->AssetValue2[i] = NULL;
        dev_data->AssetValue3[i] = NULL;
    }

    dev_data->Aux0 = NULL;
    dev_data->Aux1 = NULL;
    dev_data->Aux2 = NULL;
    dev_data->Aux3 = NULL;

    dev_data->Shift0 = NULL;
    dev_data->ShiftF = NULL;
    dev_data->Shift1 = NULL;
    dev_data->Shift2 = NULL;
    dev_data->Shift3 = NULL;

    dev_data->p0 = NULL;
    dev_data->pF = NULL;
    dev_data->p1 = NULL;
    dev_data->p2 = NULL;
    dev_data->p3 = NULL;

    dev_data->r1 = NULL;
    dev_data->r2 = NULL;
}
    
int Hyb4_Dev_Alloc(HYB4_DEV_DATA *dev_data, HYB4_TREE_DATA *tree_data)
{
    int w0 = tree_data->Width[0];
    int w1 = tree_data->Width[1];
    int w2 = tree_data->Width[2];

    int w3 = tree_data->xWidth[3];

    int m0 = tree_data->MaxIndex[0];
    int m1 = tree_data->MaxIndex[1];
    int m2 = tree_data->MaxIndex[2];
    int m3 = tree_data->MaxIndex[3];

    int NbAssetOn = tree_data->NbAssetOn;

    if (NbAssetOn > 0)
    {
        dev_data->AssetValue0[0] = (double *) DR_Array(DOUBLE, 0, m0);

        if ((tree_data->AssetType[0] == IRD) ||
            (tree_data->AssetType[0] == IRF))
        {
            dev_data->AssetValue0[1] = (double *) DR_Array(DOUBLE, 0, m0);
        }
        else
        {
            dev_data->AssetValue0[1] = dev_data->AssetValue0[0];
        }

        dev_data->Aux0 = (double *) DR_Array(DOUBLE, 0, m0);

        dev_data->Shift0 = (int *) DR_Array(INT, 0, m0);
        dev_data->ShiftF = (int *) DR_Array(INT, 0, m0);

        dev_data->p0 = (TPROB_0 *) DR_Array(TYPE_TPROB_0, 0, m0);
        dev_data->pF = (TPROB_0 *) DR_Array(TYPE_TPROB_0, 0, m0);
    }

    if (NbAssetOn > 1)
    {
        dev_data->AssetValue1[0] = (double *) DR_Array(DOUBLE, 0, m1);
    
        if ((tree_data->AssetType[1] == IRD) ||
            (tree_data->AssetType[1] == IRF))
        {
            dev_data->AssetValue1[1] = (double *) DR_Array(DOUBLE, 0, m1);
        }
        else
        {
            dev_data->AssetValue1[1] = dev_data->AssetValue1[0];
        }

        dev_data->Aux1 = (double *) DR_Array(DOUBLE, 0, m1);

        dev_data->Shift1 = (int *) DR_Array(INT, 0, m1);

        dev_data->p1 = (TPROB_0 *) DR_Array(TYPE_TPROB_0, 0, m1);

        dev_data->t1 = (double *) DR_Array(DOUBLE, 0, w1);
    }

    if (NbAssetOn > 2)
    {
        dev_data->AssetValue2[0] = (double *) DR_Array(DOUBLE, 0, m2);
    
        if ((tree_data->AssetType[2] == IRD) ||
            (tree_data->AssetType[2] == IRF))
        {
            dev_data->AssetValue2[1] = (double *) DR_Array(DOUBLE, 0, m2);
        }
        else
        {
            dev_data->AssetValue2[1] = dev_data->AssetValue2[0];
        }

        dev_data->Aux2 = (double *) DR_Array(DOUBLE, 0, m2);

        dev_data->Shift2 = (int *) DR_Array(INT, 0, m2);

        dev_data->p2 = (TPROB_0 *) DR_Array(TYPE_TPROB_0, 0, m2);

        dev_data->r1 = (TPROB_1 *) DR_Array(TYPE_TPROB_1, 0, m1);

        dev_data->t2 = (double *) DR_Array(DOUBLE, 0, w2);
    }

    if (NbAssetOn > 3)
    {
        dev_data->AssetValue3[0] = (double *) DR_Array(DOUBLE, 0, m3);

        if ((tree_data->AssetType[3] == IRD) ||
            (tree_data->AssetType[3] == IRF))
        {
            dev_data->AssetValue3[1] = (double *) DR_Array(DOUBLE, 0, m3);
        }
        else
        {
            dev_data->AssetValue3[1] = dev_data->AssetValue3[0];
        }

        dev_data->Aux3 = (double *) DR_Array(DOUBLE, 0, m3);

        dev_data->Shift3 = (int *) DR_Array(INT, 0, m3);

        dev_data->p3 = (TPROB_0 *) DR_Array(TYPE_TPROB_0, 0, m3);

        dev_data->r2 = (TPROB_2 *) DR_Array(TYPE_TPROB_2, 0, m2);

        dev_data->t3 = (double *) DR_Array(DOUBLE, 0, w3);
    }

    return SUCCESS;
}    

int Hyb4_Dev_Free(HYB4_DEV_DATA *dev_data, HYB4_TREE_DATA *tree_data)
{
    int w0 = tree_data->Width[0];
    int w1 = tree_data->Width[1];
    int w2 = tree_data->Width[2];

    int w3 = tree_data->xWidth[3];

    int m0 = tree_data->MaxIndex[0];
    int m1 = tree_data->MaxIndex[1];
    int m2 = tree_data->MaxIndex[2];
    int m3 = tree_data->MaxIndex[3];

    int NbAssetOn = tree_data->NbAssetOn;

    if (NbAssetOn > 0)
    {
        Free_DR_Array(dev_data->AssetValue0[0], DOUBLE, 0, m0);

        if ((tree_data->AssetType[0] == IRD) ||
            (tree_data->AssetType[0] == IRF))
        {
            Free_DR_Array(dev_data->AssetValue0[1], DOUBLE, 0, m0);
        }

        Free_DR_Array(dev_data->Aux0, DOUBLE, 0, m0);

        Free_DR_Array(dev_data->Shift0, INT, 0, m0);
        Free_DR_Array(dev_data->ShiftF, INT, 0, m0);

        Free_DR_Array(dev_data->p0, TYPE_TPROB_0, 0, m0);
        Free_DR_Array(dev_data->pF, TYPE_TPROB_0, 0, m0);
    }

    if (NbAssetOn > 1)
    {
        Free_DR_Array(dev_data->AssetValue1[0], DOUBLE, 0, m1);

        if ((tree_data->AssetType[1] == IRD) ||
            (tree_data->AssetType[1] == IRF))
        {
            Free_DR_Array(dev_data->AssetValue1[1], DOUBLE, 0, m1);
        }

        Free_DR_Array(dev_data->Aux1, DOUBLE, 0, m1);

        Free_DR_Array(dev_data->Shift1, INT, 0, m1);

        Free_DR_Array(dev_data->p1, TYPE_TPROB_0, 0, m1);

        Free_DR_Array(dev_data->t1, DOUBLE, 0, w1);
    }

    if (NbAssetOn > 2)
    {
        Free_DR_Array(dev_data->AssetValue2[0], DOUBLE, 0, m2);

        if ((tree_data->AssetType[2] == IRD) ||
            (tree_data->AssetType[2] == IRF))
        {
            Free_DR_Array(dev_data->AssetValue2[1], DOUBLE, 0, m2);
        }

        Free_DR_Array(dev_data->Aux2, DOUBLE, 0, m2);

        Free_DR_Array(dev_data->Shift2, INT, 0, m2);

        Free_DR_Array(dev_data->p2, TYPE_TPROB_0, 0, m2);

        Free_DR_Array(dev_data->r1, TYPE_TPROB_1, 0, m1);

        Free_DR_Array(dev_data->t2, DOUBLE, 0, w2);
    }

    if (NbAssetOn > 3)
    {
        Free_DR_Array(dev_data->AssetValue3[0], DOUBLE, 0, m3);

        if ((tree_data->AssetType[3] == IRD) ||
            (tree_data->AssetType[3] == IRF))
        {
            Free_DR_Array(dev_data->AssetValue3[1], DOUBLE, 0, m3);
        }

        Free_DR_Array(dev_data->Aux3, DOUBLE, 0, m3);

        Free_DR_Array(dev_data->Shift3, INT, 0, m3);

        Free_DR_Array(dev_data->p3, TYPE_TPROB_0, 0, m3);

        Free_DR_Array(dev_data->r2, TYPE_TPROB_2, 0, m2);

        Free_DR_Array(dev_data->t3, DOUBLE, 0, w3);
    }

    return SUCCESS;
}    

void Hyb4_Asset_IR_Init(ASSET_IR *asset)
{
    asset->MR     = 0;
    asset->QLeft  = 0;
    asset->QRight = 0;
    asset->VolBbq = 0;

    asset->FwdRateA      = NULL;
    asset->MLeft         = NULL;
    asset->MRight        = NULL;
    asset->SLeft         = NULL;
    asset->SRight        = NULL;
    asset->ZFRatio1      = NULL;
    asset->ZFRatio2      = NULL;
    asset->ZCenter       = NULL;

    asset->SpotVol       = NULL;
}

int Hyb4_Asset_IR_Alloc(ASSET_IR *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    asset->FwdRateA      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->MLeft         = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->MRight        = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->SLeft         = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->SRight        = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->ZFRatio1      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->ZFRatio2      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->ZCenter       = (double *)DR_Array(DOUBLE, -1, NbTP+1);

    asset->SpotVol       = (double *)DR_Array(DOUBLE, -1, NbTP+1);

    return SUCCESS;
}

int Hyb4_Asset_IR_Free(ASSET_IR *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    Free_DR_Array(asset->FwdRateA     , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->MLeft        , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->MRight       , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->SLeft        , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->SRight       , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->ZFRatio1     , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->ZFRatio2     , DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->ZCenter      , DOUBLE, -1, NbTP+1);

    Free_DR_Array(asset->SpotVol      , DOUBLE, -1, NbTP+1);

    return SUCCESS;
}

void Hyb4_Asset_EQ_Init(ASSET_EQ *asset)
{
    asset->Fwd      = NULL;
    asset->SpotVol  = NULL;
    asset->Center   = NULL;
}

int Hyb4_Asset_EQ_Alloc(ASSET_EQ *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    asset->Fwd      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->SpotVol  = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->Center   = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    
    return SUCCESS;
}

int Hyb4_Asset_EQ_Free(ASSET_EQ *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    Free_DR_Array(asset->Fwd,      DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->SpotVol,  DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->Center,   DOUBLE, -1, NbTP+1);

    return SUCCESS;
}

void Hyb4_Asset_FX_Init(ASSET_FX *asset)
{
    asset->Fwd      = NULL;
    asset->SpotVol  = NULL;
    asset->Center   = NULL;
}

int Hyb4_Asset_FX_Alloc(ASSET_FX *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    asset->Fwd      = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->SpotVol  = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    asset->Center   = (double *)DR_Array(DOUBLE, -1, NbTP+1);
    
    return SUCCESS;
}

int Hyb4_Asset_FX_Free(ASSET_FX *asset, HYB4_TREE_DATA *tree_data)
{
    int NbTP = tree_data->NbTP;

    Free_DR_Array(asset->Fwd,      DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->SpotVol,  DOUBLE, -1, NbTP+1);
    Free_DR_Array(asset->Center,   DOUBLE, -1, NbTP+1);

    return SUCCESS;
}

/*****  Hyb4_Alloc_Slice  ********************************************************/
/**
*         Allocate memory for one time slice.
*/
TSLICE    Hyb4_Alloc_Slice (HYB4_TREE_DATA  *tree_data, /**< (I) Tree data structure  */
                            int              dimension) /**< (I) Dim of desired slice */

{
    TSLICE Slice = NULL;

    int MaxIndex;

    if ((dimension < 1) ||
        (dimension > 4))
    {
        DR_Error("Invalid slice dimension! (Hyb4_Alloc_Slice)\n");
    }
    else
    {
        MaxIndex = tree_data->MaxIndex[dimension-1];

        Slice = (TSLICE) DR_Array(DOUBLE, 0, MaxIndex);
    }

     return Slice;

}  /* Hyb4_Alloc_Slice */


/*****  Hyb4_Free_Slice  *********************************************************/
/**
*         Free memory for time slice.
*/
void     Hyb4_Free_Slice(TSLICE            Slice,        /**< (I) Time slice            */
                         HYB4_TREE_DATA   *tree_data,    /**< (I) Tree data structure   */
                         int               dimension)
{
    int MaxIndex;

    if ((dimension < 1) ||
        (dimension > 4))
    {
        DR_Error("Invalid slice dimension! (Hyb4_Free_Slice)\n");
    }
    else
    {
        MaxIndex = tree_data->MaxIndex[dimension-1];

        Free_DR_Array(Slice, DOUBLE, 0, MaxIndex);
    }

     return;

}  /* Hyb4_Free_Slice */
