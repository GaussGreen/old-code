/*****************************************************************************/
/*        Copy structures from hyb3 to hyb4.                                 */
/*****************************************************************************/
/*        HYB4_COPY.C                                                        */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

int Hyb4_CopyTreeCups(HYB4_TREE_DATA *tn,
                      HYB3_TREE_DATA *tt)
{
    int NbTP = tt->NbTP;

    double du;

    int t;

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;
        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->R10[t] = tt->Rho[0][t];
        tn->R20[t] = tt->Rho[1][t];
        tn->R21[t] = tt->Rho[2][t];

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;
        tn->Drift2[t] = 0;
        tn->Drift3[t] = 0;
        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

int Hyb4_CopyTreeEq(HYB4_TREE_DATA *tn,
                    HYB3_TREE_DATA *tt,
                    int             RateDim,
                    int             EqDim)
{
    int NbTP = tt->NbTP;
    double du;
    double FwdDrift;

    int t;

    if (EqDim != 1)
    {
        DR_Error("Hyb4_CopyTreeEq only supports EqDim = 1.\n");
        return FAILURE;
    }

    if (RateDim != 0)
    {
        DR_Error("Hyb4_CopyTreeEq only supports RateDim = 0.\n");
        return FAILURE;
    }

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;
        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->R10[t] = tt->Rho[3][t];

        FwdDrift = (log(tt->ZeroCoupon[RateDim][1][t+1]) -
                    log(tt->ZeroCoupon[RateDim][1][t  ]));

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;
        tn->Drift2[t] = 0;
        tn->Drift3[t] = 0;

        FwdDrift /= tn->J11[t];
        tn->Drift1[t] = FwdDrift;

        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

int Hyb4_CopyTreeEq3d(HYB4_TREE_DATA *tn,
                      HYB3_TREE_DATA *tt,
                      int             RateDim,
                      int             EqDim)
{
    int NbTP = tt->NbTP;
    double du;
    double FwdDrift;

    int t;

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;

        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->J20[t] = tt->Aweight[3][t] * du;
        tn->J21[t] = tt->Aweight[4][t] * du;
        tn->J22[t] = tt->Aweight[5][t] * du;

        tn->R10[t] = tt->Rho[0][t];
        tn->R20[t] = tt->Rho[3][t];
        tn->R21[t] = tt->Rho[4][t];

        FwdDrift = (log(tt->ZeroCoupon[RateDim][1][t+1]) -
                    log(tt->ZeroCoupon[RateDim][1][t  ]));

        FwdDrift /= tn->J22[t];

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;
        tn->Drift2[t] = FwdDrift;
        tn->Drift3[t] = 0;

        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

int Hyb4_CopyTreeFx(HYB4_TREE_DATA *tn,
                    HYB3_TREE_DATA *tt,
                    int             FxDim)
{
    int NbTP = tt->NbTP;

    double du;
    double FwdDrift;

    int t;

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;

        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->J20[t] = tt->Aweight[3][t] * du;
        tn->J21[t] = tt->Aweight[4][t] * du;
        tn->J22[t] = tt->Aweight[5][t] * du;

        tn->R10[t] = tt->Rho[0][t];
        tn->R20[t] = tt->Rho[1][t];
        tn->R21[t] = tt->Rho[2][t];

        FwdDrift = (log(tt->ZeroCoupon[1][1][t+1]) -
                    log(tt->ZeroCoupon[1][1][t  ]) -
                    log(tt->ZeroCoupon[0][1][t+1]) +
                    log(tt->ZeroCoupon[0][1][t  ]) );

        FwdDrift /= tn->J22[t];

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;
        tn->Drift2[t] = FwdDrift;
        tn->Drift3[t] = 0;

        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

/* dummy fourth dimension */
int Hyb4_CopyTreeFx2FxEq(HYB4_TREE_DATA *tn,
                         HYB3_TREE_DATA *tt,
                         int             FxDim,
                         double          EqSpotVol)
{
    int NbTP = tt->NbTP;

    double du;
    double FwdDrift;

    int t;

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;

        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->J20[t] = tt->Aweight[3][t] * du;
        tn->J21[t] = tt->Aweight[4][t] * du;
        tn->J22[t] = tt->Aweight[5][t] * du;

        tn->J30[t] = 0;
        tn->J31[t] = 0;
        tn->J32[t] = 0;
        tn->J33[t] = EqSpotVol * du;

        tn->R10[t] = tt->Rho[0][t];
        tn->R20[t] = tt->Rho[1][t];
        tn->R21[t] = tt->Rho[2][t];

        tn->R30[t] = 0;
        tn->R31[t] = 0;
        tn->R32[t] = 0;

        FwdDrift = (log(tt->ZeroCoupon[1][1][t+1]) -
                    log(tt->ZeroCoupon[1][1][t  ]) -
                    log(tt->ZeroCoupon[0][1][t+1]) +
                    log(tt->ZeroCoupon[0][1][t  ]) );

        FwdDrift /= tn->J22[t];

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;
        tn->Drift2[t] = FwdDrift;

        FwdDrift = (log(tt->ZeroCoupon[1][1][t+1]) -
                    log(tt->ZeroCoupon[1][1][t  ]));

        FwdDrift /= tn->J33[t];

        tn->Drift3[t] = FwdDrift;

        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

int Hyb4_CopyTreeFxEq(HYB4_TREE_DATA *tn,
                      HYB3_TREE_DATA *tt,
                      int             EqRateDim)
{
    int NbTP = tt->NbTP;

    double du;
    double FwdDrift;

    int t;

    for (t=0;t<=NbTP; ++t)
    {
        tn->Time[t] = (double)Daysact(tt->TPDate[0],
                                      tt->TPDate[t]) / 365;

        du = sqrt(JUMPCOEFF * tt->LengthJ[t]);

        tn->J00[t] = tt->Aweight[0][t] * du;

        tn->J10[t] = tt->Aweight[1][t] * du;
        tn->J11[t] = tt->Aweight[2][t] * du;

        tn->J20[t] = tt->Aweight[3][t] * du;
        tn->J21[t] = tt->Aweight[4][t] * du;
        tn->J22[t] = tt->Aweight[5][t] * du;

        tn->J30[t] = tt->Aweight[6][t] * du;
        tn->J31[t] = tt->Aweight[7][t] * du;
        tn->J32[t] = tt->Aweight[8][t] * du;
        tn->J33[t] = tt->Aweight[9][t] * du;

        tn->R10[t] = tt->Rho[0][t];
        tn->R20[t] = tt->Rho[1][t];
        tn->R21[t] = tt->Rho[2][t];
        tn->R30[t] = tt->Rho[3][t];
        tn->R31[t] = tt->Rho[4][t];
        tn->R32[t] = tt->Rho[5][t];

        tn->Drift0[t] = tt->DriftCUPS[0][t]/tn->J00[t];
        tn->DriftF[t] = 0;
        tn->Drift1[t] = 0;

        FwdDrift = (log(tt->ZeroCoupon[1][1][t+1]) -
                    log(tt->ZeroCoupon[1][1][t  ]) -
                    log(tt->ZeroCoupon[0][1][t+1]) +
                    log(tt->ZeroCoupon[0][1][t  ]) );

        FwdDrift /= tn->J22[t];

        tn->Drift2[t] = FwdDrift;

        FwdDrift = (log(tt->ZeroCoupon[EqRateDim][1][t+1]) -
                    log(tt->ZeroCoupon[EqRateDim][1][t  ]));

        FwdDrift /= tn->J33[t];

        tn->Drift3[t] = FwdDrift;

        tn->Length[t] = tt->Length[t];
        tn->LengthJ[t] = tt->LengthJ[t];
        tn->JumpCoeff[t] = tt->Length[t]/(tt->LengthJ[t] * JUMPCOEFF);
    }

    return SUCCESS;
}

int Hyb4_CopyAssetIR(HYB4_TREE_DATA  *tn,
                     HYB3_TREE_DATA       *tt,
                     ASSET_IR             *asset,
                     MKTVOL_DATA          *m,
                     int                   e)
{
    int NbTP = tt->NbTP;

    double FwdShift;
    double QMid;
    double FwdRateA;

    int t;

    int CvDiffF;
    int CvIdx1F;
    int CvIdx2F;
    int CvDiscF;

    asset->MR = m->Beta[0];
    asset->QLeft = m->QLeft;
    asset->QRight = m->QRight;

    for (t=0; t<=NbTP; ++t)
    {
        CvDiffF = tt->CvDiff[e];    /* Internal assigments of zero curves  */
        CvIdx1F = tt->CvIdx1[e];    
        CvIdx2F = tt->CvIdx2[e];    
        CvDiscF = tt->CvDisc[e];

        asset->ZFRatio1[t] = (1. + tt->FwdRate[e][CvDiffF][t]) 
                           / (1. + tt->FwdRate[e][CvIdx1F][t]);
        asset->ZFRatio2[t] = (1. + tt->FwdRate[e][CvDiffF][t]) 
                           / (1. + tt->FwdRate[e][CvIdx2F][t]);

        FwdShift = m->FwdShift;
        QMid = (asset->QLeft + asset->QRight)/2;

        asset->VolBbq = (1. + FwdShift) / (1. + QMid * FwdShift);

        FwdRateA = tt->FwdRate[e][tt->CvDiff[e]][t] / (1. + FwdShift);

        asset->FwdRateA[t] = FwdRateA;

        asset->MLeft[t] = asset->MRight[t] = FwdRateA;
        asset->SLeft[t] = asset->SRight[t] = 1 + FwdRateA;

        if (fabs(asset->QLeft) > QCUTOFF)
        {
            asset->MLeft[t] /= asset->QLeft;
            asset->SLeft[t] -= FwdRateA / asset->QLeft;      
        }

        if (fabs(asset->QRight) > QCUTOFF)
        {
            asset->MRight[t] /= asset->QRight;
            asset->SRight[t] -= FwdRateA / asset->QRight;      
        }

        asset->ZCenter[t] = tt->IrZCenter[e][t];

        asset->SpotVol[t] = tt->SpotVol[e][t];
    }

    return SUCCESS;
}

int Hyb4_CopyAssetEQ(HYB4_TREE_DATA          *tn,
                     HYB3_TREE_DATA       *tt,
                     ASSET_EQ             *asset)
{
    int t;
    int NbTP = tt->NbTP;

    for (t=0; t<=NbTP; ++t)
    {
        asset->Fwd[t]      = tt->FwdEq[t];
        asset->SpotVol[t]  = tt->SpotEqVol[t];
        asset->Center[t]   = tt->EqMidNode[t];
    }

    return SUCCESS;
}

int Hyb4_CopyAssetFX(HYB4_TREE_DATA          *tn,
                     HYB3_TREE_DATA       *tt,
                     ASSET_FX             *asset)
{
    int t;
    int NbTP = tt->NbTP;

    for (t=0; t<=NbTP; ++t)
    {
        asset->Fwd[t]      = tt->FwdFx[t];
        asset->SpotVol[t]  = tt->SpotFxVol[t];
        asset->Center[t]   = tt->FxMidNode[t];
    }

    return SUCCESS;
}
