#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hyb4_lib.h"

/******************************************/
/*                                        */
/*   Low-level esl functions support the  */
/*   calculation of moments in all        */
/*   generic cases. (For example,         */
/*   quanto fx.) See esl_moments.c        */
/*                                        */
/*   Currently supported configurations   */
/*   at the highest level:                */
/*                                        */
/*   IRF*                                 */
/*   IRD                                  */
/*   X*                                   */
/*   Y*                                   */
/*   Z*                                   */
/*                                        */
/*   where "*" indicates an optional      */
/*   asset. Each of X, Y, Z can be        */
/*   either EQD, EQF, or FX. Further      */
/*   restrictions:                        */
/*                                        */
/*   1) If FX or EQF is present, then     */
/*      so must IRF be present.           */
/*   2) Maximum of 1 FX process.          */
/*   3) FX, EQF, and EQD are lognormal.   */
/*      (But rates are 2q.)               */
/*   4) Maximum of 4 assets.              */
/*                                        */
/******************************************/

int Hyb4_Moments_Analytic(HYB4_TREE_DATA *tree_data)
{
    int     status = FAILURE;

    int     NbAsset = tree_data->NbAsset;

    double  One  [MAXNBDATE];

    double *M[4];
    double *V[4][4];
    double *R[4][4];

    double *Time = tree_data->Time;

    double  Fw[4][MAXNBDATE];
    double  Mr[4];
    double *Sv[4];

    double *Forward;

    double dt;

    int i0, i1;
    int j, k;
    int t;

    int NbTP = tree_data->NbTP;
    int xT   = tree_data->xT;

    int Type;

    int IrDom [4];
    int IrFor [4];
    int Qx    [4];

    int IsIr  [4];
    int IsFx  [4];
    int IsEq  [4];

    M[0] = tree_data->M0;
    M[1] = tree_data->M1;
    M[2] = tree_data->M2;
    M[3] = tree_data->M3;

    V[0][0] = tree_data->V00;
    V[1][0] = tree_data->V10;
    V[1][1] = tree_data->V11;
    V[2][0] = tree_data->V20;
    V[2][1] = tree_data->V21;
    V[2][2] = tree_data->V22;
    V[3][0] = tree_data->V30;
    V[3][1] = tree_data->V31;
    V[3][2] = tree_data->V32;
    V[3][3] = tree_data->V33;

    V[0][1] = V[1][0];
    V[0][2] = V[2][0];
    V[1][2] = V[2][1];
    V[0][3] = V[3][0];
    V[1][3] = V[3][1];
    V[2][3] = V[3][2];

    for (t=0; t<MAXNBDATE; ++t)
    {
        One  [t] = 1;
    }

    R[0][0] = One;
    R[1][0] = tree_data->R10;
    R[1][1] = One;
    R[2][0] = tree_data->R20;
    R[2][1] = tree_data->R21;
    R[2][2] = One;
    R[3][0] = tree_data->R30;
    R[3][1] = tree_data->R31;
    R[3][2] = tree_data->R32;
    R[3][3] = One;

    R[0][1] = R[1][0];
    R[0][2] = R[2][0];
    R[1][2] = R[2][1];
    R[0][3] = R[3][0];
    R[1][3] = R[3][1];
    R[2][3] = R[3][2];

    /* set up spot vols, forward rates, and mean reversion */

    for (j=0; j<NbAsset; ++j)
    {
        IsIr[j] = FALSE;
        IsFx[j] = FALSE;
        IsEq[j] = FALSE;

        if ((tree_data->AssetType[j] == IRD) ||
            (tree_data->AssetType[j] == IRF))
        {
            IsIr[j] = TRUE;

            Forward = ((ASSET_IR *)(tree_data->asset[j]))->FwdRateA;

            for (t=0; t<=NbTP; ++t)
            {
                dt = tree_data->Length[t];

                Fw[j][t] = Forward[t] / dt;
            }

            Mr[j] = ((ASSET_IR *)(tree_data->asset[j]))->MR;

            Sv[j] = ((ASSET_IR *)(tree_data->asset[j]))->SpotVol;
        }
        else if (tree_data->AssetType[j] == FX)
        {
            IsFx[j] = TRUE;

            for (t=0; t<=NbTP; ++t)
            {
                Fw[j][t] = 0;
            }

            Mr[j] = 0;

            Sv[j] = ((ASSET_FX *)(tree_data->asset[j]))->SpotVol;
        }
        else
        {
            IsEq[j] = TRUE;

            for (t=0; t<=NbTP; ++t)
            {
                Fw[j][t] = 0;
            }

            Mr[j] = 0;

            Sv[j] = ((ASSET_EQ *)(tree_data->asset[j]))->SpotVol;
        }
    }

    /* set up domestic / foreign / qx */

    for (j=0; j<NbAsset; ++j)
    {
        if (IsIr[j])
        {
            Qx[j] = ((ASSET_IR *)(tree_data->asset[j]))->Fx;
        }

        if (IsFx[j])
        {
            IrDom[j] = ((ASSET_FX *)(tree_data->asset[j]))->IrDom;
            IrFor[j] = ((ASSET_FX *)(tree_data->asset[j]))->IrFor;

            Qx[j] = ((ASSET_IR *)(tree_data->asset[IrDom[j]]))->Fx;
        }

        if (IsEq[j])
        {
            IrDom[j] = ((ASSET_EQ *)(tree_data->asset[j]))->IrDom;

            Qx[j] = ((ASSET_IR *)(tree_data->asset[IrDom[j]]))->Fx;
        }
    }

    /* moments */

    for (j=0; j<NbAsset; ++j)
    {
        Type = tree_data->AssetType[j];

        if (Type == IRF)
        {
            TreeMean_IrExt(M[j],
                           xT,
                           Time,
                           Mr[j],
                           Sv[j],
                           Sv[Qx[j]],
                           R[j][Qx[j]]);
        }
        else if (Type == IRD)
        {
            TreeMean_IrNum(M[j],
                           xT);
        }
        else if (Type == FX)
        {
            TreeMean_FxNum(M[j],
                           xT,
                           Time,
                           Fw[IrDom[j]],
                           Mr[IrDom[j]],
                           Sv[IrDom[j]],
                           Fw[IrFor[j]],
                           Mr[IrFor[j]],
                           Sv[IrFor[j]],
                           Sv[j],
                           R[IrFor[j]][j]);
        }
        else if (Type == EQD)
        {
            TreeMean_EqNum(M[j],
                           xT,
                           Time,
                           Fw[IrDom[j]],
                           Mr[IrDom[j]],
                           Sv[IrDom[j]]);
        }
        else if (Type == EQF)
        {
            TreeMean_EqExt(M[j],
                           xT,
                           Time,
                           Fw[IrDom[j]],
                           Mr[IrDom[j]],
                           Sv[IrDom[j]],
                           Sv[Qx[j]],
                           R[IrDom[j]][Qx[j]]);
        }
        else
        {
            DR_Error("Unknown asset type.\n");
            goto RETURN;
        }

        for (k=j; k<NbAsset; ++k)
        {
            if (IsIr[j] && IsIr[k])
            {
                TreeCovariance_IrIr(V[j][k],
                                    xT,
                                    Time,
                                    Mr[j],
                                    Sv[j],
                                    Mr[k],
                                    Sv[k],
                                    R[j][k]);
            }
            else if ((IsIr[j] && IsFx[k]) ||
                     (IsIr[k] && IsFx[j]))
            {
                i0 = IsIr[j] ? j : k;
                i1 = IsFx[k] ? k : j;

                TreeCovariance_IrFx(V[j][k],
                                    xT,
                                    Time,
                                    Mr[i0],
                                    Sv[i0],
                                    Fw[IrFor[i1]],
                                    Mr[IrFor[i1]],
                                    Sv[IrFor[i1]],
                                    Fw[IrDom[i1]],
                                    Mr[IrDom[i1]],
                                    Sv[IrDom[i1]],
                                    Sv[i1],
                                    R[i0][IrFor[i1]],
                                    R[i0][IrDom[i1]],
                                    R[i0][i1]);
            }
            else if ((IsIr[j] && IsEq[k]) ||
                     (IsIr[k] && IsEq[j]))
            {
                i0 = IsIr[j] ? j : k;
                i1 = IsEq[k] ? k : j;

                TreeCovariance_IrEq(V[j][k],
                                    xT,
                                    Time,
                                    Mr[i0],
                                    Sv[i0],
                                    Fw[IrDom[i1]],
                                    Mr[IrDom[i1]],
                                    Sv[IrDom[i1]],
                                    Sv[i1],
                                    R[i0][IrDom[i1]],
                                    R[i0][i1]);
            }
            else if (IsFx[j] && IsFx[k])
            {
                TreeCovariance_FxFx(V[j][k],
                                    xT,
                                    Time,
                                    Fw[IrFor[j]],
                                    Mr[IrFor[j]],
                                    Sv[IrFor[j]],
                                    Fw[IrDom[j]],
                                    Mr[IrDom[j]],
                                    Sv[IrDom[j]],
                                    Fw[IrFor[k]],
                                    Mr[IrFor[k]],
                                    Sv[IrFor[k]],
                                    Fw[IrDom[k]],
                                    Mr[IrDom[k]],
                                    Sv[IrDom[k]],
                                    Sv[j],
                                    Sv[k],
                                    R[IrFor[j]][IrFor[k]],
                                    R[IrFor[j]][IrDom[k]],
                                    R[IrFor[j]][k],
                                    R[IrDom[j]][IrFor[k]],
                                    R[IrDom[j]][IrDom[k]],
                                    R[IrDom[j]][k],
                                    R[j][IrFor[k]],
                                    R[j][IrDom[k]],
                                    R[j][k]);
            }
            else if ((IsFx[j] && IsEq[k]) ||
                     (IsFx[k] && IsEq[j]))
            {
                i0 = IsFx[j] ? j : k;
                i1 = IsEq[k] ? k : j;

                TreeCovariance_FxEq(V[j][k],
                                    xT,
                                    Time,
                                    Fw[IrFor[i0]],
                                    Mr[IrFor[i0]],
                                    Sv[IrFor[i0]],
                                    Fw[IrDom[i0]],
                                    Mr[IrDom[i0]],
                                    Sv[IrDom[i0]],
                                    Fw[IrDom[i1]],
                                    Mr[IrDom[i1]],
                                    Sv[IrDom[i1]],
                                    Sv[i0],
                                    Sv[i1],
                                    R[IrFor[i0]][IrDom[i1]],
                                    R[IrFor[i0]][i1],
                                    R[IrDom[i0]][IrDom[i1]],
                                    R[IrDom[i0]][i1],
                                    R[i0][IrDom[i1]],
                                    R[i0][i1]);
            }
            else if (IsEq[j] && IsEq[k])
            {
                TreeCovariance_EqEq(V[j][k],
                                    xT,
                                    Time,
                                    Fw[IrDom[j]],
                                    Mr[IrDom[j]],
                                    Sv[IrDom[j]],
                                    Fw[IrDom[k]],
                                    Mr[IrDom[k]],
                                    Sv[IrDom[k]],
                                    Sv[j],
                                    Sv[k],
                                    R[IrDom[j]][IrDom[k]],
                                    R[IrDom[j]][k],
                                    R[j][IrDom[k]],
                                    R[j][k]);
            }
            else
            {
                DR_Error("Unknown tree configuration.\n");
                goto RETURN;
            }
        }
    }

    status = SUCCESS;

RETURN:

    return status;
}

int Hyb4_Cholesky4d(double R[4][4],
                    double N[4][4],
                    double G[4][4])
{
    int status = FAILURE;

    int u, v;

    double det;

    /* check for positivity */

    det = 1 -       R[1][0] * R[1][0]
            -       R[2][0] * R[2][0]
            -       R[2][1] * R[2][1]
            -       R[3][0] * R[3][0]
            -       R[3][1] * R[3][1]
            -       R[3][2] * R[3][2]

            + 2 * ( R[2][1] * R[3][1] * R[3][2]
            +       R[2][0] * R[3][0] * R[3][2]
            +       R[1][0] * R[3][0] * R[3][1]
            +       R[1][0] * R[2][0] * R[2][1] )

            +       R[1][0] * R[1][0] * R[3][2] * R[3][2]
            +       R[2][0] * R[2][0] * R[3][1] * R[3][1]
            +       R[2][1] * R[2][1] * R[3][0] * R[3][0]

            - 2 * ( R[1][0] * R[2][1] * R[3][2] * R[3][0]
            +       R[1][0] * R[2][0] * R[3][1] * R[3][2]
            +       R[2][0] * R[2][1] * R[3][0] * R[3][1] );

    if (det < TINY)
    {
        DR_Error("Resulting composite correl mtx is not positive\n");
        goto RETURN;
    }

    /* R = N * N^T */

    N[0][0] = 1;

    N[1][0] = R[1][0] / N[0][0];

    N[1][1] = sqrt(1 - SQUARE(N[1][0]));

    N[2][0] = R[2][0] / N[0][0];

    N[2][1] = (R[2][1] - N[2][0] * N[1][0]) / N[1][1];

    N[2][2] = sqrt(1 - SQUARE(N[2][0]) - SQUARE(N[2][1]));

    N[3][0] = R[3][0] / N[0][0];

    N[3][1] = (R[3][1] - N[3][0] * N[1][0]) / N[1][1];

    N[3][2] = (R[3][2] - N[3][0] * N[2][0] - N[3][1] * N[2][1]) / N[2][2];

    N[3][3] = sqrt(1 - SQUARE(N[3][0]) - SQUARE(N[3][1]) - SQUARE(N[3][2]));

    /* G = inverse of N */

    G[0][0] = 1;

    G[1][1] = 1 / N[1][1];

    G[1][0] = - G[1][1] * N[1][0] / N[0][0];

    G[2][2] = 1 / N[2][2];

    G[2][1] = - G[2][2] * N[2][1] / N[1][1];

    G[2][0] = - (G[2][1] * N[1][0] + G[2][2] * N[2][0]) / N[0][0];

    G[3][3] = 1 / N[3][3];

    G[3][2] = - G[3][3] * N[3][2] / N[2][2];

    G[3][1] = - (G[3][2] * N[2][1] + G[3][3] * N[3][1]) / N[1][1];

    G[3][0] = - (G[3][1] * N[1][0] + G[3][2] * N[2][0] + G[3][3] * N[3][0]) / N[0][0];

    for (u=0; u<4; ++u)
    {
        for (v=u+1; v<4; ++v)
        {
            N[u][v] = 0;
            G[u][v] = 0;
        }
    }

    status = SUCCESS;

RETURN:

    return SUCCESS;
}

int Hyb4_Moments_Empirical(HYB4_TREE_DATA *tree_data, LATTICE_PROG *LP)
{
    int    NbAsset = tree_data->NbAsset;

    double *M[4];
    double *V[4][4];

    double  Sigma[4];
    double  R[4][4];
    double  N[4][4];
    double  G[4][4];

    double X[4];

    int    status = FAILURE;

    int    t, i, j, k, L, u, v;

    double PJ00;

    double PJ10;
    double PJ11;

    double PJ20;
    double PJ21;
    double PJ22;

    double PJ30;
    double PJ31;
    double PJ32;
    double PJ33;

    double X0i;

    double X1i;
    double X1ij;

    double X2i;
    double X2ij;
    double X2ijk;

    double X3i;
    double X3ij;
    double X3ijk;
    double X3ijkL;

    double X3Lo;
    double X3Hi;

    double YLo0;
    double YLo1;
    double YLo2;
    double YLo3;

    double YHi0;
    double YHi1;
    double YHi2;
    double YHi3;

    double YLo;
    double YHi;

    double YLoMin = 0;
    double YLoMax = 0;
    double YHiMin = 0;
    double YHiMax = 0;

    int    iMin;
    int    iMax;

    int   *jMin;
    int   *jMax;

    int  **kMin;
    int  **kMax;

    int ***LMin;
    int ***LMax;

    int    jMini;
    int    jMaxi;
    int   *kMini;
    int   *kMaxi;
    int  **LMini;
    int  **LMaxi;

    int    kMinij;
    int    kMaxij;
    int   *LMinij;
    int   *LMaxij;

    int    LMinijk;
    int    LMaxijk;

    int    IsFirstIJK;

    double SP;
    double SPsum;

    TSLICE StatePr  = NULL;
    TSLICE StatePr0 = NULL;
    TSLICE StatePr1 = NULL;
    TSLICE Temp     = NULL;

    int xT   = tree_data->xT;

    HYB4_DEV_DATA dev_data;

    Hyb4_Dev_Init(&dev_data);

    if (Hyb4_Dev_Alloc(&dev_data, tree_data) == FAILURE)
    {
        goto RETURN;
    }

    StatePr0 = Hyb4_Alloc_Slice(tree_data, 4);
    StatePr1 = Hyb4_Alloc_Slice(tree_data, 4);

    if ((StatePr0 == NULL) ||
        (StatePr1 == NULL))
    {
        goto RETURN;
    }

    M[0] = tree_data->M0;
    M[1] = tree_data->M1;
    M[2] = tree_data->M2;
    M[3] = tree_data->M3;

    V[0][0] = tree_data->V00;
    V[1][0] = tree_data->V10;
    V[1][1] = tree_data->V11;
    V[2][0] = tree_data->V20;
    V[2][1] = tree_data->V21;
    V[2][2] = tree_data->V22;
    V[3][0] = tree_data->V30;
    V[3][1] = tree_data->V31;
    V[3][2] = tree_data->V32;
    V[3][3] = tree_data->V33;

    V[0][1] = V[1][0];
    V[0][2] = V[2][0];
    V[1][2] = V[2][1];
    V[0][3] = V[3][0];
    V[1][3] = V[3][1];
    V[2][3] = V[3][2];

    if (NbAsset != 4)
    {
        DR_Error("Only 4-d mode is supported by Hyb4_Moments_Empirical.\n");
        goto RETURN;
    }

    for (u=0; u<4; ++u)
    {
        M[u][0] = 0;
        
        for (v=0; v<=u; ++v)
        {
            V[u][v][0] = 0;
        }
    }

    tree_data->StateProb[0] = 1;

    if (Hyb4_UpdateStatePrices4D(-1,
                                 0,
                                 tree_data,
                                 &dev_data,
                                 LP,
                                 StatePr0,
                                 StatePr1) == FAILURE)
    {
        goto RETURN;
    }

    Temp     = StatePr0;
    StatePr0 = StatePr1;
    StatePr1 = Temp;

    /* moments and state probabilities */

    for (t=1; t<=xT; ++t)
    {
        if (Hyb4_UpdateStatePrices4D(-1,
                                     t,
                                     tree_data,
                                     &dev_data,
                                     LP,
                                     StatePr0,
                                     StatePr1) == FAILURE)
        {
            goto RETURN;
        }

        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        jMin = tree_data->jMin[t];
        jMax = tree_data->jMax[t];

        kMin = tree_data->kMin[t];
        kMax = tree_data->kMax[t];

        LMin = tree_data->LMin[t];
        LMax = tree_data->LMax[t];

        PJ00 = tree_data->J00[t-1];

        PJ10 = tree_data->J10[t-1];
        PJ11 = tree_data->J11[t-1];

        PJ20 = tree_data->J20[t-1];
        PJ21 = tree_data->J21[t-1];
        PJ22 = tree_data->J22[t-1];

        PJ30 = tree_data->J30[t-1];
        PJ31 = tree_data->J31[t-1];
        PJ32 = tree_data->J32[t-1];
        PJ33 = tree_data->J33[t-1];

        for (u=0; u<4; ++u)
        {
            M[u][t] = 0;
            
            for (v=0; v<=u; ++v)
            {
                V[u][v][t] = 0;
            }
        }

        SPsum = 0;

        for (i=iMin; i<=iMax; ++i)
        {
            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            LMini = LMin[i];
            LMaxi = LMax[i];

            X0i = PJ00 * i;
            X1i = PJ10 * i;
            X2i = PJ20 * i;
            X3i = PJ30 * i;

            X[0] = X0i;

            for (j=jMini; j<=jMaxi; ++j)
            {
                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                LMinij = LMini[j];
                LMaxij = LMaxi[j];

                X1ij = X1i + PJ11 * j;
                X2ij = X2i + PJ21 * j;
                X3ij = X3i + PJ31 * j;

                X[1] = X1ij;

                for (k=kMinij; k<=kMaxij; ++k)
                {
                    LMinijk = LMinij[k];
                    LMaxijk = LMaxij[k];

                    X2ijk = X2ij + PJ22 * k;
                    X3ijk = X3ij + PJ32 * k;

                    X[2] = X2ijk;

                    StatePr = StatePr0 + tree_data->NodeOffset3[t][i][j][k];

                    for (L=LMinijk; L<=LMaxijk; ++L)
                    {
                        X3ijkL = X3ijk + PJ33 * L;

                        X[3] = X3ijkL;

                        SP = StatePr[L];

                        /* moments */

                        for (u=0; u<4; ++u)
                        {
                            M[u][t] += X[u] * SP;

                            for (v=0; v<=u; ++v)
                            {
                                V[u][v][t] += X[u] * X[v] * SP;
                            }
                        }

                        SPsum += SP;
                    }
                }
            }
        }

        for (u=0; u<4; ++u)
        {
            for (v=0; v<=u; ++v)
            {
                V[u][v][t] -= M[u][t] * M[v][t];
            }
        }

        tree_data->StateProb[t] = SPsum;

        Temp     = StatePr0;
        StatePr0 = StatePr1;
        StatePr1 = Temp;
    }

    /* limits */

    for (t=1; t<=xT; ++t)
    {
        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        jMin = tree_data->jMin[t];
        jMax = tree_data->jMax[t];

        kMin = tree_data->kMin[t];
        kMax = tree_data->kMax[t];

        LMin = tree_data->LMin[t];
        LMax = tree_data->LMax[t];

        PJ00 = tree_data->J00[t-1];

        PJ10 = tree_data->J10[t-1];
        PJ11 = tree_data->J11[t-1];

        PJ20 = tree_data->J20[t-1];
        PJ21 = tree_data->J21[t-1];
        PJ22 = tree_data->J22[t-1];

        PJ30 = tree_data->J30[t-1];
        PJ31 = tree_data->J31[t-1];
        PJ32 = tree_data->J32[t-1];
        PJ33 = tree_data->J33[t-1];

        IsFirstIJK = 1;

        tree_data->LimLoMin[t] = 0;
        tree_data->LimLoMax[t] = 0;
        tree_data->LimHiMin[t] = 0;
        tree_data->LimHiMax[t] = 0;

        for (u=0; u<4; ++u)
        {
            Sigma[u] = sqrt(V[u][u][t]);
        }

        if ((fabs(Sigma[0]) < TINY) ||
            (fabs(Sigma[1]) < TINY) ||
            (fabs(Sigma[2]) < TINY) ||
            (fabs(Sigma[3]) < TINY))
        {
            continue;
        }

        for (u=0; u<4; ++u)
        {
            for (v=0; v<4; ++v)
            {
                R[u][v] = V[u][v][t] / (Sigma[u] * Sigma[v]);
            }
        }

        if (Hyb4_Cholesky4d(R, N, G) == FAILURE)
        {
            continue;
        }

        for (i=iMin; i<=iMax; ++i)
        {
            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            LMini = LMin[i];
            LMaxi = LMax[i];

            X0i = PJ00 * i;
            X1i = PJ10 * i;
            X2i = PJ20 * i;
            X3i = PJ30 * i;

            X[0] = X0i;

            for (j=jMini; j<=jMaxi; ++j)
            {
                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                LMinij = LMini[j];
                LMaxij = LMaxi[j];

                X1ij = X1i + PJ11 * j;
                X2ij = X2i + PJ21 * j;
                X3ij = X3i + PJ31 * j;

                X[1] = X1ij;

                for (k=kMinij; k<=kMaxij; ++k)
                {
                    LMinijk = LMinij[k];
                    LMaxijk = LMaxij[k];

                    X2ijk = X2ij + PJ22 * k;
                    X3ijk = X3ij + PJ32 * k;

                    X[2] = X2ijk;

                    X3Lo = X3ijk + PJ33 * LMinijk;
                    X3Hi = X3ijk + PJ33 * LMaxijk;

                    YLo0 = (G[0][0] / Sigma[0]) * X[0];

                    YLo1 = (G[1][0] / Sigma[0]) * X[0] +
                           (G[1][1] / Sigma[1]) * X[1];

                    YLo2 = (G[2][0] / Sigma[0]) * X[0] +
                           (G[2][1] / Sigma[1]) * X[1] +
                           (G[2][2] / Sigma[2]) * X[2];

                    YLo3 = (G[3][0] / Sigma[0]) * X[0] +
                           (G[3][1] / Sigma[1]) * X[1] +
                           (G[3][2] / Sigma[2]) * X[2] +
                           (G[3][3] / Sigma[3]) * X3Lo;

                    YLo = sqrt(YLo0 * YLo0 +
                               YLo1 * YLo1 +
                               YLo2 * YLo2 +
                               YLo3 * YLo3);

                    YHi0 = YLo0;
                    YHi1 = YLo1;
                    YHi2 = YLo2;

                    YHi3 = (G[3][0] / Sigma[0]) * X[0] +
                           (G[3][1] / Sigma[1]) * X[1] +
                           (G[3][2] / Sigma[2]) * X[2] +
                           (G[3][3] / Sigma[3]) * X3Hi;

                    YHi = sqrt(YHi0 * YHi0 +
                               YHi1 * YHi1 +
                               YHi2 * YHi2 +
                               YHi3 * YHi3);

                    if (IsFirstIJK)
                    {
                        YLoMin = YLo;
                        YLoMax = YLo;
                        YHiMin = YHi;
                        YHiMax = YHi;
                    }
                    else
                    {
                        YLoMin = MIN(YLoMin, YLo);
                        YLoMax = MAX(YLoMax, YLo);
                        YHiMin = MIN(YHiMin, YHi);
                        YHiMax = MAX(YHiMax, YHi);
                    }

                    IsFirstIJK = 0;
                }
            }
        }

        tree_data->LimLoMin[t] = YLoMin;
        tree_data->LimLoMax[t] = YLoMax;
        tree_data->LimHiMin[t] = YHiMin;
        tree_data->LimHiMax[t] = YHiMax;
    }

    status = SUCCESS;

RETURN:

    Hyb4_Dev_Free(&dev_data, tree_data);

    Hyb4_Free_Slice(StatePr0, tree_data, 4);
    Hyb4_Free_Slice(StatePr1, tree_data, 4);

    return status;
}

int Hyb4_Moments_Print(HYB4_TREE_DATA *tree_data,
                       char           *FileName0,
                       char           *FileName1)
{
    int   status = FAILURE;
    int   t      = 0;
    int   xT     = tree_data->xT;
    FILE *fp     = NULL;

    fp = fopen(FileName0, "w");

    if (fp == NULL)
    {
        DR_Error("Could not open file for printing means.\n");
        goto RETURN;
    }

    for (t=0; t<=xT; ++t)
    {
        if (t<100)
        {
            fprintf(fp, " ");
        }

        if (t<10)
        {
            fprintf(fp, " ");
        }

        fprintf(fp, "%d ", t);

        fprintf(fp, "%16.8lf %16.8lf %16.8lf %16.8lf\n",
                    tree_data->M0[t],
                    tree_data->M1[t],
                    tree_data->M2[t],
                    tree_data->M3[t]);
    }

    fclose(fp);

    fp = NULL;

    fp = fopen(FileName1, "w");

    if (fp == NULL)
    {
        DR_Error("Could not open file for printing covariances.\n");
        goto RETURN;
    }

    for (t=0; t<=xT; ++t)
    {
        if (t<100)
        {
            fprintf(fp, " ");
        }

        if (t<10)
        {
            fprintf(fp, " ");
        }

        fprintf(fp, "%d ", t);

        fprintf(fp, "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf"
                    "%16.8lf\n",
                    tree_data->V00[t],
                    tree_data->V10[t],
                    tree_data->V11[t],
                    tree_data->V20[t],
                    tree_data->V21[t],
                    tree_data->V22[t],
                    tree_data->V30[t],
                    tree_data->V31[t],
                    tree_data->V32[t],
                    tree_data->V33[t]);
    }

    fclose(fp);

    status = SUCCESS;

RETURN:

    return status;
}

int Hyb4_StateProbs_Print(HYB4_TREE_DATA *tree_data,
                          char           *FileName)
{
    int   status = FAILURE;
    int   t      = 0;
    int   xT     = tree_data->xT;
    FILE *fp     = NULL;

    fp = fopen(FileName, "w");

    if (fp == NULL)
    {
        DR_Error("Could not open file for printing state probs.\n");
        goto RETURN;
    }

    for (t=0; t<=xT; ++t)
    {
        if (t<100)
        {
            fprintf(fp, " ");
        }

        if (t<10)
        {
            fprintf(fp, " ");
        }

        fprintf(fp, "%d ", t);

        fprintf(fp, "%16.8lf\n",
                    tree_data->StateProb[t]);
    }

    fclose(fp);

    status = SUCCESS;

RETURN:

    return status;
}

int Hyb4_Limits_Print(HYB4_TREE_DATA *tree_data,
                      char           *FileName)
{
    int   status = FAILURE;
    int   t      = 0;
    int   xT     = tree_data->xT;
    FILE *fp     = NULL;

    fp = fopen(FileName, "w");

    if (fp == NULL)
    {
        DR_Error("Could not open file for printing empirical limits.\n");
        goto RETURN;
    }

    for (t=0; t<=xT; ++t)
    {
        if (t<100)
        {
            fprintf(fp, " ");
        }

        if (t<10)
        {
            fprintf(fp, " ");
        }

        fprintf(fp, "%d ", t);

        fprintf(fp, "%16.8lf "
                    "%16.8lf "
                    "%16.8lf "
                    "%16.8lf\n",
                    tree_data->LimLoMin[t],
                    tree_data->LimLoMax[t],
                    tree_data->LimHiMin[t],
                    tree_data->LimHiMax[t]);
    }

    fclose(fp);

    status = SUCCESS;

RETURN:

    return status;
}

int Hyb4_AssetPayoff_Empirical(HYB4_TREE_DATA *tree_data,
                               LATTICE_PROG   *LP,
                               int             DfDim,
                               int             AssetDim,
                               double         *X)
{
    int    NbAssetOn = tree_data->NbAssetOn;

    int    status = FAILURE;

    int    t, i, j, k, L;

    double val = 0;

    int    iMin;
    int    iMax;

    int   *jMin;
    int   *jMax;

    int  **kMin;
    int  **kMax;

    int ***LMin;
    int ***LMax;

    int    jMini;
    int    jMaxi;
    int   *kMini;
    int   *kMaxi;
    int  **LMini;
    int  **LMaxi;

    int    kMinij;
    int    kMaxij;
    int   *LMinij;
    int   *LMaxij;

    int    LMinijk;
    int    LMaxijk;

    int    offset0;
    int    offset1;
    int    offset2;
    int    offset3;

    TSLICE StatePr  = NULL;
    TSLICE StatePr0 = NULL;
    TSLICE StatePr1 = NULL;
    TSLICE Temp     = NULL;

    TSLICE AV0;
    TSLICE AV1;
    TSLICE AV2;
    TSLICE AV3;

    int xT   = tree_data->xT;

    HYB4_DEV_DATA dev_data;

    Hyb4_Dev_Init(&dev_data);

    if (Hyb4_Dev_Alloc(&dev_data, tree_data) == FAILURE)
    {
        goto RETURN;
    }

    StatePr0 = Hyb4_Alloc_Slice(tree_data, 4);
    StatePr1 = Hyb4_Alloc_Slice(tree_data, 4);

    if ((StatePr0 == NULL) ||
        (StatePr1 == NULL))
    {
        goto RETURN;
    }

    if (AssetDim >= tree_data->NbAssetOn)
    {
        DR_Error("Can't price asset number %d; only %d assets are on.\n",
                 AssetDim,
                 tree_data->NbAssetOn);

        goto RETURN;
    }

    if (DfDim >= tree_data->NbAssetOn)
    {
        DR_Error("Can't discount with asset number %d; only %d assets are on.\n",
                 DfDim,
                 tree_data->NbAssetOn);

        goto RETURN;
    }

    if (NbAssetOn != 4)
    {
        DR_Error("Only 4d mode is currently supported "
                 "(Hyb4_AssetPayoff_Empirical).\n");

        goto RETURN;
    }

    if (Hyb4_UpdateStatePrices4D(DfDim,
                                 0,
                                 tree_data,
                                 &dev_data,
                                 LP,
                                 StatePr0,
                                 StatePr1) == FAILURE)
    {
        goto RETURN;
    }

    Temp     = StatePr0;
    StatePr0 = StatePr1;
    StatePr1 = Temp;

    for (t=1; t<=xT; ++t)
    {
        if (Hyb4_UpdateStatePrices4D(DfDim,
                                     t,
                                     tree_data,
                                     &dev_data,
                                     LP,
                                     StatePr0,
                                     StatePr1) == FAILURE)
        {
            goto RETURN;
        }

        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        jMin = tree_data->jMin[t];
        jMax = tree_data->jMax[t];

        kMin = tree_data->kMin[t];
        kMax = tree_data->kMax[t];

        LMin = tree_data->LMin[t];
        LMax = tree_data->LMax[t];

        X[t] = 0;

        offset0 = tree_data->NodeOffset0[t];

        AV0     = dev_data.AssetValue0[0] + offset0;

        for (i=iMin; i<=iMax; ++i)
        {
            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            LMini = LMin[i];
            LMaxi = LMax[i];

            if (AssetDim == 0)
            {
                val = AV0[i];
            }

            offset1 = tree_data->NodeOffset1[t][i];

            AV1     = dev_data.AssetValue1[0] + offset1;

            for (j=jMini; j<=jMaxi; ++j)
            {
                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                LMinij = LMini[j];
                LMaxij = LMaxi[j];

                if (AssetDim == 1)
                {
                    val = AV1[j];
                }

                offset2 = tree_data->NodeOffset2[t][i][j];

                AV2     = dev_data.AssetValue2[0] + offset2;

                for (k=kMinij; k<=kMaxij; ++k)
                {
                    LMinijk = LMinij[k];
                    LMaxijk = LMaxij[k];

                    if (AssetDim == 2)
                    {
                        val = AV2[k];
                    }

                    offset3 = tree_data->NodeOffset3[t][i][j][k];

                    AV3     = dev_data.AssetValue3[0] + offset3;
                    StatePr = StatePr0                + offset3;

                    for (L=LMinijk; L<=LMaxijk; ++L)
                    {
                        if (AssetDim == 3)
                        {
                            val = AV3[L];
                        }

                        X[t] += val * StatePr[L];
                    }
                }
            }
        }

        Temp     = StatePr0;
        StatePr0 = StatePr1;
        StatePr1 = Temp;
    }

    status = SUCCESS;

RETURN:

    Hyb4_Dev_Free(&dev_data, tree_data);

    Hyb4_Free_Slice(StatePr0, tree_data, 4);
    Hyb4_Free_Slice(StatePr1, tree_data, 4);

    return status;
}

