/*****************************************************************************/
/*        Drifts and transition probabilities.                               */
/*****************************************************************************/
/*        HYB4_LATTICE.C                                                     */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hyb4_lib.h"

#ifndef INIT_GRID
#define INIT_GRID(Q,M,V,i) (((fabs(Q)) > (QCUTOFF)) ? ((M)*(exp((Q)*(V)*(i)))) : ((M)*(V)*(i)))
#endif

#ifndef UPDT_GRID
#define UPDT_GRID(Q,Grid,J) (((fabs(Q)) > (QCUTOFF)) ? ((Grid)*(J)) : ((Grid)+(J)))
#endif

#ifndef CALC_PROB
#define CALC_PROB(_p, _JumpCoeff, _c)        \
                                             \
    _p->u = 0.5 * (_JumpCoeff + _c + _c * _c); \
    _p->d = _p->u - _c;                      \
    _p->m = 1 - _p->u - _p->d;

#endif

int GetAssetParamIR(double    *MRFactor,
                    double    *y1,
                    double    *y2,
                    HYB4_TREE_DATA *tree_data,
                    int       *t,
                    int       *dim)
{
    *MRFactor = - tree_data->Length[*t] * ((ASSET_IR *)(tree_data->asset[*dim]))->MR;
    return SUCCESS;
}

int GetAssetParamEQ(double    *Fwd,
                    double    *y1,
                    double    *y2,
                    HYB4_TREE_DATA *tree_data,
                    int       *t,
                    int       *dim)
{
    return SUCCESS;
}

int GetAssetParamFX(double    *y0,
                    double    *y1,
                    double    *y2,
                    HYB4_TREE_DATA *tree_data,
                    int       *t,
                    int       *dim)
{
    return SUCCESS;
}

int InitializeAssetValueIR(double    *X,
                           double    *PJ,
                           double    *Discount0,
                           double    *Discount1,
                           HYB4_DEV_DATA  *dev_data,
                           HYB4_TREE_DATA *tree_data,
                           int       *t,
                           int       *dim,
                           int       *indexMin,
                           int       *indexMax)
{
    ASSET_IR *asset = (ASSET_IR *)(tree_data->asset[*dim]);

    double    QLeft         =  asset->QLeft;
    double    QRight        =  asset->QRight;
    double    VolBbq        =  asset->VolBbq;
    double    FwdRateA      =  asset->FwdRateA[*t];
    double    MLeft         =  asset->MLeft[*t];
    double    MRight        =  asset->MRight[*t];
    double    SLeft         =  asset->SLeft[*t];
    double    SRight        =  asset->SRight[*t];
    double    ZFRatio1      =  asset->ZFRatio1[*t];
    double    ZCenter       =  asset->ZCenter[*t];

    double    RateJumpLeft;
    double    RateJumpRight;
    double    Zidx;
    int       Mid;

    int a;
    
    double    GridIr;

    RateJumpLeft = RateJumpRight = *PJ * VolBbq * FwdRateA;

    if (fabs(QLeft)  > QCUTOFF)
    {
        RateJumpLeft  = exp(*PJ * QLeft  * VolBbq);
    }

    if (fabs(QRight) > QCUTOFF)
    {
        RateJumpRight = exp(*PJ * QRight * VolBbq);
    }

    Mid = (int) ceil(- ZCenter / *PJ) - 1;
    Mid = MIN ( MAX (Mid, *indexMin - 1), *indexMax);

    Zidx = ZCenter + *X + *PJ * *indexMin;

    GridIr = INIT_GRID(QLeft, MLeft, VolBbq, Zidx);

    for (a=*indexMin; a<=Mid; ++a)
    {
        Discount0[a] = 1 / (SLeft + GridIr);
        Discount1[a] = ZFRatio1 * Discount0[a];
  
        GridIr = UPDT_GRID(QLeft, GridIr, RateJumpLeft);
    }

    Zidx = ZCenter + *X + *PJ * (Mid + 1);
    
    GridIr = INIT_GRID(QRight, MRight, VolBbq, Zidx);

    for (a=Mid+1; a<=*indexMax; ++a)
    {
        Discount0[a] = 1 / (SRight + GridIr);
        Discount1[a] = ZFRatio1 * Discount0[a];
       
        GridIr = UPDT_GRID(QRight, GridIr, RateJumpRight);
    }

    return SUCCESS;
}

int InitializeAssetValueEQ(double    *X,
                           double    *PJ,
                           double    *Eq,
                           double    *EqNone,
                           HYB4_DEV_DATA  *dev_data,
                           HYB4_TREE_DATA *tree_data,
                           int       *t,
                           int       *dim,
                           int       *indexMin,
                           int       *indexMax)
{
    ASSET_EQ *asset    = (ASSET_EQ *)(tree_data->asset[*dim]);

    double    Fwd      = asset->Fwd[*t];
    double    Factor   = exp(*PJ);
    double    Center   = asset->Center[*t];
    double    EqValue  = Fwd * exp(Center + *X + *PJ * *indexMin);

    int a;

    for (a=*indexMin; a<=*indexMax; ++a)
    {
        Eq[a]    = EqValue;
        EqValue *= Factor;
    }

    return SUCCESS;
}

int InitializeAssetValueFX(double    *X,
                           double    *PJ,
                           double    *Fx,
                           double    *FxNone,
                           HYB4_DEV_DATA  *dev_data,
                           HYB4_TREE_DATA *tree_data,
                           int       *t,
                           int       *dim,
                           int       *indexMin,
                           int       *indexMax)
{
    ASSET_FX *asset    = (ASSET_FX *)(tree_data->asset[*dim]);

    double    Fwd      = asset->Fwd[*t];
    double    Factor   = exp(*PJ);
    double    Center   = asset->Center[*t];
    double    FxValue  = Fwd * exp(Center + *X + *PJ * *indexMin);

    int a;

    for (a=*indexMin; a<=*indexMax; ++a)
    {
        Fx[a]    = FxValue;
        FxValue *= Factor;
    }

    return SUCCESS;
}

int UpdateDriftNONE(double *a0,
                    double *a1,
                    double *a2,
                    double *b0,
                    double *b1,
                    double *b2,
                    double *DF0,
                    double *DF1,
                    double *Drift,
                    int    *i,
                    double *A,
                    double *M,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int *t)
{
    return SUCCESS;
}

int UpdateDriftIR(double *a0,
                  double *a1,
                  double *a2,
                  double *MRFactor,
                  double *b1,
                  double *b2,
                  double *DF0,
                  double *DF1,
                  double *Drift,
                  int    *i,
                  double *A,
                  double *M,
                  HYB4_DEV_DATA  *dev_data,
                  HYB4_TREE_DATA *tree_data,
                  int *t)
{
    *Drift += (*MRFactor) * (*i) * (*M);
    return SUCCESS;
}

int UpdateDriftEQF(double *a0,
                   double *a1,
                   double *a2,
                   double *b0,
                   double *b1,
                   double *b2,
                   double *DF0,
                   double *DF1,
                   double *Drift,
                   int    *i,
                   double *A,
                   double *M,
                   HYB4_DEV_DATA  *dev_data,
                   HYB4_TREE_DATA *tree_data,
                   int *t)
{
    return SUCCESS;
}

int UpdateDriftEQD(double *a0,
                    double *a1,
                    double *a2,
                    double *b0,
                    double *b1,
                    double *b2,
                    double *DF0,
                    double *DF1,
                    double *Drift,
                    int    *i,
                    double *A,
                    double *M,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int *t)
{
    return SUCCESS;
}

int UpdateDriftFX(double *a0,
                    double *a1,
                    double *a2,
                    double *b0,
                    double *b1,
                    double *b2,
                    double *DF0,
                    double *DF1,
                    double *Drift,
                    int    *i,
                    double *A,
                    double *M,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int *t)
{
    return SUCCESS;
}

int UpdateDriftIRD_FX(double *a0,
                    double *a1,
                    double *a2,
                    double *b0,
                    double *b1,
                    double *b2,
                    double *DF0,
                    double *DF1,
                    double *Drift,
                    int    *i,
                    double *A,
                    double *M,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int *t)
{
    *Drift -= log(*DF1) * (*A);

    return SUCCESS;
}

int UpdateDriftIRF_FX(double *a0,
                    double *a1,
                    double *a2,
                    double *b0,
                    double *b1,
                    double *b2,
                    double *DF0,
                    double *DF1,
                    double *Drift,
                    int    *i,
                    double *A,
                    double *M,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int *t)
{
    *Drift += log(*DF1) * (*A);

    return SUCCESS;
}

int UpdateDriftIR_EQ(double *a0,
                     double *a1,
                     double *a2,
                     double *b0,
                     double *b1,
                     double *b2,
                     double *DF0,
                     double *DF1,
                     double *Drift,
                     int    *i,
                     double *A,
                     double *M,
                     HYB4_DEV_DATA  *dev_data,
                     HYB4_TREE_DATA *tree_data,
                     int *t)
{
    *Drift -= log(*DF1) * (*A);

    return SUCCESS;
}

int FinalizeDriftIR(double    *AV,
                    double    *Drift,
                    int       *Shift,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int       *t)
{
    return SUCCESS;
}

int FinalizeDriftEQF(double    *AV,
                    double    *Drift,
                    int       *Shift,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int       *t)
{
    return SUCCESS;
}

int FinalizeDriftEQD(double    *AV,
                    double    *Drift,
                    int       *Shift,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int       *t)
{
    return SUCCESS;
}

int FinalizeDriftFX(double    *AV,
                    double    *Drift,
                    int       *Shift,
                    HYB4_DEV_DATA  *dev_data,
                    HYB4_TREE_DATA *tree_data,
                    int       *t)
{
    return SUCCESS;
}

/* configuration of the lattice code, depending on which asset type */
/* appears in each dimension; called once before pricing begins     */

int Hyb4_ConfigureLatticeCode(HYB4_TREE_DATA *tree_data, LATTICE_PROG *LP)
{
    int           *AssetType = tree_data->AssetType;
    int            NbAssetOn = tree_data->NbAssetOn;

    int            n;
    int            z;

    for (n=0; n<NbAssetOn; ++n)
    {
        switch (AssetType[n])
        {
            case IRF:

                LP->LC[n].GetAssetParam        = &GetAssetParamIR;
                LP->LC[n].InitializeAssetValue = &InitializeAssetValueIR;
                LP->LC[n].FinalizeDrift        = &FinalizeDriftIR;

                for (z=n; z<4; ++z)
                {
                    switch (AssetType[z])
                    {
                        case IRF: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR;     break;
                        case IRD: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR;     break;
                        case FX:  LP->LC[n].UpdateDrift[z] = &UpdateDriftIRF_FX; break;
                        case EQF: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR_EQ;  break;
                        case EQD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE;   break;
                        default:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE;   break;
                    }
                }

                break;

            case IRD:

                LP->LC[n].GetAssetParam        = &GetAssetParamIR;
                LP->LC[n].InitializeAssetValue = &InitializeAssetValueIR;
                LP->LC[n].FinalizeDrift        = &FinalizeDriftIR;

                for (z=n; z<4; ++z)
                {
                    switch (AssetType[z])
                    {
                        case IRF: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR;     break;
                        case IRD: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR;     break;
                        case FX:  LP->LC[n].UpdateDrift[z] = &UpdateDriftIRD_FX; break;
                        case EQF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE;   break;
                        case EQD: LP->LC[n].UpdateDrift[z] = &UpdateDriftIR_EQ;  break;
                        default:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE;   break;
                    }
                }

                break;

            case FX:

                LP->LC[n].GetAssetParam        = &GetAssetParamFX;
                LP->LC[n].InitializeAssetValue = &InitializeAssetValueFX;
                LP->LC[n].FinalizeDrift        = &FinalizeDriftFX;

                for (z=n; z<4; ++z)
                {
                    switch (AssetType[z])
                    {
                        case IRF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case IRD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case FX:  LP->LC[n].UpdateDrift[z] = &UpdateDriftFX;   break;
                        case EQF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case EQD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        default:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                    }
                }

                break;

            case EQF:

                LP->LC[n].GetAssetParam        = &GetAssetParamEQ;
                LP->LC[n].InitializeAssetValue = &InitializeAssetValueEQ;
                LP->LC[n].FinalizeDrift        = &FinalizeDriftEQF;

                for (z=n; z<4; ++z)
                {
                    switch (AssetType[z])
                    {
                        case IRF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case IRD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case FX:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case EQF: LP->LC[n].UpdateDrift[z] = &UpdateDriftEQF;  break;
                        case EQD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        default:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                    }
                }

                break;

            case EQD:

                LP->LC[n].GetAssetParam        = &GetAssetParamEQ;
                LP->LC[n].InitializeAssetValue = &InitializeAssetValueEQ;
                LP->LC[n].FinalizeDrift        = &FinalizeDriftEQD;

                for (z=n; z<4; ++z)
                {
                    switch (AssetType[z])
                    {
                        case IRF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case IRD: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case FX:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case EQF: LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                        case EQD: LP->LC[n].UpdateDrift[z] = &UpdateDriftEQD;  break;
                        default:  LP->LC[n].UpdateDrift[z] = &UpdateDriftNONE; break;
                    }
                }

                break;
        }
    }

    return SUCCESS;
}

int Hyb4_Lattice(HYB4_DEV_DATA    *dev_data,
                 HYB4_TREE_DATA   *tree_data,
                 LATTICE_PROG  *LP,
                 int            t,
                 int            T)
{ 
    GET_ASSET_PARAM         GetAssetParam_0        = LP->LC[0].GetAssetParam;
    GET_ASSET_PARAM         GetAssetParam_1        = LP->LC[1].GetAssetParam;
    GET_ASSET_PARAM         GetAssetParam_2        = LP->LC[2].GetAssetParam;
    GET_ASSET_PARAM         GetAssetParam_3        = LP->LC[3].GetAssetParam;

    INITIALIZE_ASSET_VALUE  InitializeAssetValue_0 = LP->LC[0].InitializeAssetValue;
    INITIALIZE_ASSET_VALUE  InitializeAssetValue_1 = LP->LC[1].InitializeAssetValue;
    INITIALIZE_ASSET_VALUE  InitializeAssetValue_2 = LP->LC[2].InitializeAssetValue;
    INITIALIZE_ASSET_VALUE  InitializeAssetValue_3 = LP->LC[3].InitializeAssetValue;

    UPDATE_DRIFT            UpdateDrift_0_0        = LP->LC[0].UpdateDrift[0];
    UPDATE_DRIFT            UpdateDrift_0_1        = LP->LC[0].UpdateDrift[1];
    UPDATE_DRIFT            UpdateDrift_0_2        = LP->LC[0].UpdateDrift[2];
    UPDATE_DRIFT            UpdateDrift_0_3        = LP->LC[0].UpdateDrift[3];

    UPDATE_DRIFT            UpdateDrift_1_1        = LP->LC[1].UpdateDrift[1];
    UPDATE_DRIFT            UpdateDrift_1_2        = LP->LC[1].UpdateDrift[2];
    UPDATE_DRIFT            UpdateDrift_1_3        = LP->LC[1].UpdateDrift[3];

    UPDATE_DRIFT            UpdateDrift_2_2        = LP->LC[2].UpdateDrift[2];
    UPDATE_DRIFT            UpdateDrift_2_3        = LP->LC[2].UpdateDrift[3];

    UPDATE_DRIFT            UpdateDrift_3_3        = LP->LC[3].UpdateDrift[3];

    FINALIZE_DRIFT          FinalizeDrift_0        = LP->LC[0].FinalizeDrift;
    FINALIZE_DRIFT          FinalizeDrift_1        = LP->LC[1].FinalizeDrift;
    FINALIZE_DRIFT          FinalizeDrift_2        = LP->LC[2].FinalizeDrift;
    FINALIZE_DRIFT          FinalizeDrift_3        = LP->LC[3].FinalizeDrift;

    int NbAssetOn = tree_data->NbAssetOn;

    int Asset_0_On = (NbAssetOn > 0);
    int Asset_1_On = (NbAssetOn > 1);
    int Asset_2_On = (NbAssetOn > 2);
    int Asset_3_On = (NbAssetOn > 3);

    int Asset_1_Off = !Asset_1_On;
    int Asset_2_Off = !Asset_2_On;
    int Asset_3_Off = !Asset_3_On;

    static int    I0 = 0;
    static int    I1 = 1;
    static int    I2 = 2;
    static int    I3 = 3;

    double JumpCoeff = tree_data->JumpCoeff[t];

    int iIndexMin = 0;
    int iIndexMax = 0;

    int jIndexMin = 0;
    int jIndexMax = 0;

    int kIndexMin = 0;
    int kIndexMax = 0;

    int LIndexMin = 0;
    int LIndexMax = 0;

    int NodeOffset0 = 0;
    
    int *NodeOffset1 = NULL;
    int  NodeOffset1i = 0;
    
    int **NodeOffset2 = NULL;
    int  *NodeOffset2i = NULL;
    int   NodeOffset2ij = 0;
    
    int ***NodeOffset3 = NULL;
    int  **NodeOffset3i = NULL;
    int   *NodeOffset3ij = NULL;
    int    NodeOffset3ijk = 0;
    
    int iMin = 0;
    int iMax = 0;

    int *jMin = NULL;
    int *jMax = NULL;
    int  jMini = 0;
    int  jMaxi = 0;

    int **kMin = NULL;
    int **kMax = NULL;
    int  *kMini = NULL;
    int  *kMaxi = NULL;
    int   kMinij = 0;
    int   kMaxij = 0;

    int ***LMin = NULL;
    int ***LMax = NULL;
    int  **LMini = NULL;
    int  **LMaxi = NULL;
    int   *LMinij = NULL;
    int   *LMaxij = NULL;
    int    LMinijk = 0;
    int    LMaxijk = 0;

    int iOutMin = 0;
    int iOutMax = 0;

    int *jOutMin = NULL;
    int *jOutMax = NULL;

    int **kOutMin = NULL;
    int **kOutMax = NULL;

    int ***LOutMin = NULL;
    int ***LOutMax = NULL;

    int iShiftDn;
    int iShiftMd;
    int iShiftUp;

    int jShiftDn;
    int jShiftMd;
    int jShiftUp;

    int kShiftDn;
    int kShiftMd;
    int kShiftUp;

    double J00 = 0;
    double J10 = 0, J11 = 0;
    double J20 = 0, J21 = 0, J22 = 0;
    double J30 = 0, J31 = 0, J32 = 0, J33 = 0;

    double PJ00 = 0;
    double PJ10 = 0, PJ11 = 0;
    double PJ20 = 0, PJ21 = 0, PJ22 = 0;
    double PJ30 = 0, PJ31 = 0, PJ32 = 0, PJ33 = 0;

    double A0 = 0;
    double A1 = 0;
    double A2 = 0;
    double A3 = 0;

    double B00 = 0;
    double B10 = 0, B11 = 0;
    double B20 = 0, B21 = 0, B22 = 0;
    double B30 = 0, B31 = 0, B32 = 0, B33 = 0;

    double M00 = 0;
    double M10 = 0, M11 = 0;
    double M20 = 0, M21 = 0, M22 = 0;
    double M30 = 0, M31 = 0, M32 = 0, M33 = 0;

    double D00 = 0;
    double D10 = 0, D11 = 0;
    double D20 = 0, D21 = 0, D22 = 0;
    double D30 = 0, D31 = 0, D32 = 0, D33 = 0;

    double Drift0 = 0;
    double DriftF = 0;
    double Drift1 = 0;
    double Drift2 = 0;
    double Drift3 = 0;

    double Drift0i = 0;
    double DriftFi = 0;
    double Drift1i = 0;
    double Drift2i = 0;
    double Drift3i = 0;

    double Drift1ij = 0;
    double Drift2ij = 0;
    double Drift3ij = 0;

    double Drift2ijk = 0;
    double Drift3ijk = 0;

    double Drift3ijkL = 0;

    double X0=0;
    double X1=0;
    double X2=0;
    double X3=0;

    double X1i = 0;
    double X2i = 0;
    double X3i = 0;

    double X2ij = 0;
    double X3ij = 0;
    
    double X3ijk = 0;

    TPROB_0  *p0 = NULL;
    TPROB_0  *p0i = NULL;

    TPROB_0  *pF = NULL;
    TPROB_0  *pFi = NULL;

    TPROB_0  *p1i = NULL;
    TPROB_0  *p1ij = NULL;

    TPROB_0  *p2ij = NULL;
    TPROB_0  *p2ijk = NULL;

    TPROB_0  *p3ijk = NULL;
    TPROB_0  *p3ijkL = NULL;

    TPROB_1  *r1i = NULL;
    TPROB_1  *r1ij = NULL;

    TPROB_2  *r2ij = NULL;
    TPROB_2  *r2ijk = NULL;

    int *Shift0 = NULL;
    int *Shift0i = NULL;

    int *ShiftF = NULL;
    int *ShiftFi = NULL;

    int *Shift1i = NULL;
    int *Shift1ij = NULL;

    int *Shift2ij = NULL;
    int *Shift2ijk = NULL;

    int *Shift3ijk = NULL;
    int *Shift3ijkL = NULL;

    double    *AV00; double   *AV00i;
    double    *AV01; double   *AV01i;

    double   *AV10i; double  *AV10ij;
    double   *AV11i; double  *AV11ij;

    double  *AV20ij; double *AV20ijk;
    double  *AV21ij; double *AV21ijk;

    double *AV30ijk; double *AV30ijkL;
    double *AV31ijk; double *AV31ijkL;

    double y00=0; double y01=0; double y02=0;
    double y10=0; double y11=0; double y12=0;
    double y20=0; double y21=0; double y22=0;
    double y30=0; double y31=0; double y32=0;

    int assetType;

    int i, j, k, L;

    int xT = tree_data->xT;

    double elapsed;
    clock_t start, end;

    start = clock();

    /* space-invariant operations */

    if (Asset_0_On)
    {
        PJ00 = tree_data->J00[t-1];

        J00 = tree_data->J00[t];

        A0 = 1 / J00;

        M00 = PJ00 / J00;

        B00 = 1;

        D00 = (PJ00 - J00) / J00;

        Drift0 = tree_data->Drift0[t];
        DriftF = tree_data->DriftF[t];

        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        if (t < T)
        {
            iOutMin = tree_data->iOutMin[t+1];
            iOutMax = tree_data->iOutMax[t+1];
        }

        assetType = tree_data->AssetType[0];

        (*GetAssetParam_0)(&y00, &y01, &y02, tree_data, &t, &I0);

        if (Asset_1_On)
        {
            PJ10 = tree_data->J10[t-1];
            PJ11 = tree_data->J11[t-1];

            J10 = tree_data->J10[t];
            J11 = tree_data->J11[t];

            A1 = 1 / J11;
    
            M10 = PJ10 / J11;
            M11 = PJ11 / J11;

            B10 = J10 / J11;
            B11 = 1;

            D10 = (PJ10 - J10) / J11;
            D11 = (PJ11 - J11) / J11;

            Drift1 = tree_data->Drift1[t];
    
            jMin = tree_data->jMin[t];
            jMax = tree_data->jMax[t];

            if (t < T)
            {
                jOutMin = tree_data->jOutMin[t+1];
                jOutMax = tree_data->jOutMax[t+1];
            }

            assetType = tree_data->AssetType[1];

            (*GetAssetParam_1)(&y10, &y11, &y12, tree_data, &t, &I1);

            if (Asset_2_On)
            {
                PJ20 = tree_data->J20[t-1];
                PJ21 = tree_data->J21[t-1];
                PJ22 = tree_data->J22[t-1];

                J20 = tree_data->J20[t];
                J21 = tree_data->J21[t];
                J22 = tree_data->J22[t];

                A2 = 1 / J22;
        
                M20 = PJ20 / J22;
                M21 = PJ21 / J22;
                M22 = PJ22 / J22;

                B20 = J20 / J22;
                B21 = J21 / J22;
                B22 = 1;

                D20 = (PJ20 - J20) / J22;
                D21 = (PJ21 - J21) / J22;
                D22 = (PJ22 - J22) / J22;

                Drift2 = tree_data->Drift2[t];
        
                kMin = tree_data->kMin[t];
                kMax = tree_data->kMax[t];

                if (t < T)
                {
                    kOutMin = tree_data->kOutMin[t+1];
                    kOutMax = tree_data->kOutMax[t+1];
                }

                assetType = tree_data->AssetType[2];

                (*GetAssetParam_2)(&y20, &y21, &y22, tree_data, &t, &I2);
                
                if (Asset_3_On)
                {
                    PJ30 = tree_data->J30[t-1];
                    PJ31 = tree_data->J31[t-1];
                    PJ32 = tree_data->J32[t-1];
                    PJ33 = tree_data->J33[t-1];

                    J30 = tree_data->J30[t];
                    J31 = tree_data->J31[t];
                    J32 = tree_data->J32[t];
                    J33 = tree_data->J33[t];

                    A3 = 1 / J33;

                    M30 = PJ30 / J33;
                    M31 = PJ31 / J33;
                    M32 = PJ32 / J33;
                    M33 = PJ33 / J33;

                    B30 = J30 / J33;
                    B31 = J31 / J33;
                    B32 = J32 / J33;
                    B33 = 1;

                    D30 = (PJ30 - J30) / J33;
                    D31 = (PJ31 - J31) / J33;
                    D32 = (PJ32 - J32) / J33;
                    D33 = (PJ33 - J33) / J33;

                    Drift3 = tree_data->Drift3[t];

                    LMin = tree_data->LMin[t];
                    LMax = tree_data->LMax[t];

                    if (t < xT)
                    {
                        LOutMin = tree_data->LOutMin[t+1];
                        LOutMax = tree_data->LOutMax[t+1];
                    }

                    assetType = tree_data->AssetType[3];

                    (*GetAssetParam_3)(&y30, &y31, &y32, tree_data, &t, &I3);
                }
            }
        }
    }

    /* loop through the grid */

    if (Asset_0_On)
    {
        NodeOffset0 = tree_data->NodeOffset0[t];
        NodeOffset1 = tree_data->NodeOffset1[t];
        NodeOffset2 = tree_data->NodeOffset2[t];
        NodeOffset3 = tree_data->NodeOffset3[t];

        AV00 = dev_data->AssetValue0[0] + NodeOffset0;
        AV01 = dev_data->AssetValue0[1] + NodeOffset0;

        (*InitializeAssetValue_0)(&X0, &PJ00,
                                  AV00, AV01,
                                  dev_data,
                                  tree_data,
                                  &t,
                                  &I0,
                                  &iMin,
                                  &iMax);

        if ((Asset_1_Off) &&
            (t == T))
        {
            return SUCCESS;
        }

        if (t < T)
        {
            iIndexMin = iOutMin + 1;
            iIndexMax = iOutMax - 1;
        }

        Shift0 = dev_data->Shift0 + NodeOffset0;
        ShiftF = dev_data->ShiftF + NodeOffset0;

        p0 = dev_data->p0 + NodeOffset0;
        pF = dev_data->pF + NodeOffset0;

        for (i=iMin; i<=iMax; ++i)
        {
            AV00i = AV00 + i;
            AV01i = AV01 + i;

            Shift0i = Shift0 + i;
            ShiftFi = ShiftF + i;

            p0i = &(p0[i]);
            pFi = &(pF[i]);

            Drift0i = Drift0 + D00 * i;
            DriftFi = DriftF + D00 * i;

            (*UpdateDrift_0_0)(&y00, &y01, &y02,
                               &y00, &y01, &y02,
                               AV00i, AV01i,
                               &Drift0i, &i, &A0, &M00,
                               dev_data,
                               tree_data,
                               &t);

            (*UpdateDrift_0_0)(&y00, &y01, &y02,
                               &y00, &y01, &y02,
                               AV00i, AV01i,
                               &DriftFi, &i, &A0, &M00,
                               dev_data,
                               tree_data,
                               &t);

            if (Asset_1_On)
            {
                X1i = X1 + PJ10 * i;

                Drift1i = Drift1 + D10 * i;

                (*UpdateDrift_0_1)(&y00, &y01, &y02,
                                   &y10, &y11, &y12,
                                   AV00i, AV01i,
                                   &Drift1i, &i, &A1, &M10,
                                   dev_data,
                                   tree_data,
                                   &t);

                Drift1i -= Drift0i * B10;

                if (Asset_2_On)
                {
                    X2i = X2 + PJ20 * i;

                    Drift2i = Drift2 + D20 * i;

                    (*UpdateDrift_0_2)(&y00, &y01, &y02,
                                       &y20, &y21, &y22,
                                       AV00i, AV01i,
                                       &Drift2i, &i, &A2, &M20,
                                       dev_data,
                                       tree_data,
                                       &t);

                    Drift2i -= Drift0i * B20;

                    if (Asset_3_On)
                    {
                        X3i = X3 + PJ30 * i;

                        if (t < xT)
                        {
                            Drift3i = Drift3 + D30 * i;

                            (*UpdateDrift_0_3)(&y00, &y01, &y02,
                                               &y30, &y31, &y32,
                                               AV00i, AV01i,
                                               &Drift3i, &i, &A3, &M30,
                                               dev_data,
                                               tree_data,
                                               &t);

                            Drift3i -= Drift0i * B30;
                        }
                    }
                }
            }

            *Shift0i = NEAR_INT(Drift0i);
            *ShiftFi = NEAR_INT(DriftFi);

            if (t < T) /* VEZI */
            {
                *Shift0i = MIN(MAX(iIndexMin - i, *Shift0i), iIndexMax - i);
                *ShiftFi = MIN(MAX(iIndexMin - i, *ShiftFi), iIndexMax - i);
            }

            (*FinalizeDrift_0)(AV00i,
                               &Drift0i,
                                Shift0i,
                               dev_data,
                               tree_data,
                               &t);

            (*FinalizeDrift_0)(AV00i,
                               &DriftFi,
                                ShiftFi,
                               dev_data,
                               tree_data,
                               &t);

            Drift0i -= *Shift0i;
            DriftFi -= *ShiftFi;

            CALC_PROB(p0i, JumpCoeff, Drift0i);
            CALC_PROB(pFi, JumpCoeff, DriftFi);

            if (Asset_1_On)
            {
                jMini = jMin[i];
                jMaxi = jMax[i];

                NodeOffset1i = NodeOffset1[i];

                if (Asset_2_On)
                {
                    kMini = kMin[i];
                    kMaxi = kMax[i];

                    NodeOffset2i = NodeOffset2[i];

                    if (Asset_3_On && (t <= xT))
                    {
                        LMini = LMin[i];
                        LMaxi = LMax[i];

                        NodeOffset3i = NodeOffset3[i];
                    }
                }

                AV10i = dev_data->AssetValue1[0] + NodeOffset1i;
                AV11i = dev_data->AssetValue1[1] + NodeOffset1i;

                (*InitializeAssetValue_1)(&X1i, &PJ11,
                                          AV10i, AV11i,
                                          dev_data,
                                          tree_data,
                                          &t,
                                          &I1,
                                          &jMini,
                                          &jMaxi);

                if ((Asset_2_Off) &&
                    (t >= T))
                {
                    continue;
                }

                iShiftUp = (iShiftMd = (iShiftDn = i + *Shift0i - 1) + 1) + 1;

                if (t < T)
                {
                    jIndexMin =                jOutMin[iShiftDn] ;
                    jIndexMin = MAX(jIndexMin, jOutMin[iShiftMd]);
                    jIndexMin = MAX(jIndexMin, jOutMin[iShiftUp]) + 1;

                    jIndexMax =                jOutMax[iShiftDn] ;
                    jIndexMax = MIN(jIndexMax, jOutMax[iShiftMd]);
                    jIndexMax = MIN(jIndexMax, jOutMax[iShiftUp]) - 1;

                    if (jIndexMin > jIndexMax)
                    {
                        DR_Error("Problem building the tree: IndexMin > IndexMax.");
                    }
                }

                Shift1i = dev_data->Shift1 + NodeOffset1i;

                p1i = dev_data->p1 + NodeOffset1i;

                if (Asset_2_On)
                {
                    r1i = dev_data->r1 + NodeOffset1i;
                }

                for (j=jMini; j<=jMaxi; ++j)
                {
                    AV10ij = AV10i + j;
                    AV11ij = AV11i + j;

                    Shift1ij = Shift1i + j;

                    p1ij = &(p1i[j]);

                    if (Asset_2_On)
                    {
                        r1ij = &(r1i[j]);
                    }

                    Drift1ij = Drift1i + D11 * j;

                    (*UpdateDrift_1_1)(&y10, &y11, &y12,
                                       &y10, &y11, &y12,
                                       AV10ij, AV11ij,
                                       &Drift1ij, &j, &A1, &M11,
                                       dev_data,
                                       tree_data,
                                       &t);

                    if (Asset_2_On)
                    {
                        X2ij = X2i + PJ21 * j;

                        Drift2ij = Drift2i + D21 * j;

                        (*UpdateDrift_1_2)(&y10, &y11, &y12,
                                           &y20, &y21, &y22,
                                           AV10ij, AV11ij,
                                           &Drift2ij, &j, &A2, &M21,
                                           dev_data,
                                           tree_data,
                                           &t);

                        Drift2ij -= Drift1ij * B21;

                        if (Asset_3_On)
                        {
                            X3ij = X3i + PJ31 * j;

                            if (t < xT)
                            {
                                Drift3ij = Drift3i + D31 * j;

                                (*UpdateDrift_1_3)(&y10, &y11, &y12,
                                                   &y30, &y31, &y32,
                                                   AV10ij, AV11ij,
                                                   &Drift3ij, &j, &A3, &M31,
                                                   dev_data,
                                                   tree_data,
                                                   &t);

                                Drift3ij -= Drift1ij * B31;
                            }
                        }
                    }

                    *Shift1ij = NEAR_INT(Drift1ij);

                    if (t < T) /* VEZI */
                    {
                        *Shift1ij = MIN(MAX(jIndexMin - j, *Shift1ij), jIndexMax - j);
                    }
                        
                    (*FinalizeDrift_1)(AV10ij,
                                       &Drift1ij,
                                        Shift1ij,
                                       dev_data,
                                       tree_data,
                                       &t);

                    Drift1ij -= *Shift1ij;

                    CALC_PROB(p1ij, JumpCoeff, Drift1ij);

                    if (Asset_2_On)
                    {
                        kMinij = kMini[j];
                        kMaxij = kMaxi[j];

                        NodeOffset2ij = NodeOffset2i[j];

                        if (Asset_3_On && (t <= xT))
                        {
                            LMinij = LMini[j];
                            LMaxij = LMaxi[j];
                            
                            NodeOffset3ij = NodeOffset3i[j];
                        }

                        AV20ij = dev_data->AssetValue2[0] + NodeOffset2ij;
                        AV21ij = dev_data->AssetValue2[1] + NodeOffset2ij;

                        (*InitializeAssetValue_2)(&X2ij, &PJ22,
                                                  AV20ij, AV21ij,
                                                  dev_data,
                                                  tree_data,
                                                  &t,
                                                  &I2,
                                                  &kMinij,
                                                  &kMaxij);

                        if ((Asset_3_Off) &&
                            (t >= T))
                        {
                            continue;
                        }

                        r1ij->dd = p0i->d * p1ij->d;
                        r1ij->dm = p0i->d * p1ij->m;
                        r1ij->du = p0i->d * p1ij->u;
                        r1ij->md = p0i->m * p1ij->d;
                        r1ij->mm = p0i->m * p1ij->m;
                        r1ij->mu = p0i->m * p1ij->u;
                        r1ij->ud = p0i->u * p1ij->d;
                        r1ij->um = p0i->u * p1ij->m;
                        r1ij->uu = p0i->u * p1ij->u;

                        jShiftUp = (jShiftMd = (jShiftDn = j + *Shift1ij - 1) + 1) + 1;

                        if (t < T)
                        {
                            kIndexMin =                kOutMin[iShiftDn][jShiftDn] ;
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftDn][jShiftMd]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftDn][jShiftUp]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftMd][jShiftDn]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftMd][jShiftMd]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftMd][jShiftUp]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftUp][jShiftDn]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftUp][jShiftMd]);
                            kIndexMin = MAX(kIndexMin, kOutMin[iShiftUp][jShiftUp]) + 1;

                            kIndexMax =                kOutMax[iShiftDn][jShiftDn] ;
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftDn][jShiftMd]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftDn][jShiftUp]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftMd][jShiftDn]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftMd][jShiftMd]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftMd][jShiftUp]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftUp][jShiftDn]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftUp][jShiftMd]);
                            kIndexMax = MIN(kIndexMax, kOutMax[iShiftUp][jShiftUp]) - 1;
                        
                            if (kIndexMin > kIndexMax)
                            {
                                DR_Error ("Problem in building the tree (IndexMin > IndexMax)!");
                            }
                        }

                        Shift2ij = dev_data->Shift2 + NodeOffset2ij;

                        p2ij = dev_data->p2 + NodeOffset2ij;

                        if (Asset_3_On && (t < xT))
                        {
                            r2ij = dev_data->r2 + NodeOffset2ij;
                        }

                        for (k=kMinij; k<=kMaxij; ++k)
                        {
                            AV20ijk = AV20ij + k;
                            AV21ijk = AV21ij + k;

                            Shift2ijk = Shift2ij + k;

                            p2ijk = &(p2ij[k]);

                            if (Asset_3_On && (t < xT))
                            {
                                r2ijk = &(r2ij[k]);
                            }

                            Drift2ijk = Drift2ij + D22 * k;

                            (*UpdateDrift_2_2)(&y20, &y21, &y22,
                                               &y20, &y21, &y22,
                                               AV20ijk, AV21ijk,
                                               &Drift2ijk, &k, &A2, &M22,
                                               dev_data,
                                               tree_data,
                                               &t);
 
                            if (Asset_3_On)
                            {
                                X3ijk = X3ij + PJ32 * k;

                                if (t < xT)
                                {
                                    Drift3ijk = Drift3ij + D32 * k;

                                    (*UpdateDrift_2_3)(&y20, &y21, &y22,
                                                       &y30, &y31, &y32,
                                                       AV20ijk, AV21ijk,
                                                       &Drift3ijk, &k, &A3, &M32,
                                                       dev_data,
                                                       tree_data,
                                                       &t);

                                    Drift3ijk -= Drift2ijk * B32;
                                }
                            }

                            *Shift2ijk = NEAR_INT(Drift2ijk);
 
                            if (t < T) /* VEZI */
                            {
                                *Shift2ijk = MIN(MAX(kIndexMin - k, *Shift2ijk), kIndexMax - k);
                            }

                            (*FinalizeDrift_2)(AV20ijk,
                                               &Drift2ijk,
                                                Shift2ijk,
                                               dev_data,
                                               tree_data,
                                               &t);

                            Drift2ijk -= *Shift2ijk;

                            CALC_PROB(p2ijk, JumpCoeff, Drift2ijk);

                            if ((Asset_3_On) && (t <= xT))
                            {
                                LMinijk = LMinij[k];
                                LMaxijk = LMaxij[k];

                                NodeOffset3ijk = NodeOffset3ij[k];

                                AV30ijk = dev_data->AssetValue3[0] + NodeOffset3ijk;
                                AV31ijk = dev_data->AssetValue3[1] + NodeOffset3ijk;

                                (*InitializeAssetValue_3)(&X3ijk, &PJ33,
                                                          AV30ijk, AV31ijk,
                                                          dev_data,
                                                          tree_data,
                                                          &t,
                                                          &I3,
                                                          &LMinijk,
                                                          &LMaxijk);

                                if (t >= xT)
                                {
                                    continue;
                                }

                                r2ijk->ddd = r1ij->dd * p2ijk->d;
                                r2ijk->ddm = r1ij->dd * p2ijk->m;
                                r2ijk->ddu = r1ij->dd * p2ijk->u;
                                r2ijk->dmd = r1ij->dm * p2ijk->d;
                                r2ijk->dmm = r1ij->dm * p2ijk->m;
                                r2ijk->dmu = r1ij->dm * p2ijk->u;
                                r2ijk->dud = r1ij->du * p2ijk->d;
                                r2ijk->dum = r1ij->du * p2ijk->m;
                                r2ijk->duu = r1ij->du * p2ijk->u;
                                r2ijk->mdd = r1ij->md * p2ijk->d;
                                r2ijk->mdm = r1ij->md * p2ijk->m;
                                r2ijk->mdu = r1ij->md * p2ijk->u;
                                r2ijk->mmd = r1ij->mm * p2ijk->d;
                                r2ijk->mmm = r1ij->mm * p2ijk->m;
                                r2ijk->mmu = r1ij->mm * p2ijk->u;
                                r2ijk->mud = r1ij->mu * p2ijk->d;
                                r2ijk->mum = r1ij->mu * p2ijk->m;
                                r2ijk->muu = r1ij->mu * p2ijk->u;
                                r2ijk->udd = r1ij->ud * p2ijk->d;
                                r2ijk->udm = r1ij->ud * p2ijk->m;
                                r2ijk->udu = r1ij->ud * p2ijk->u;
                                r2ijk->umd = r1ij->um * p2ijk->d;
                                r2ijk->umm = r1ij->um * p2ijk->m;
                                r2ijk->umu = r1ij->um * p2ijk->u;
                                r2ijk->uud = r1ij->uu * p2ijk->d;
                                r2ijk->uum = r1ij->uu * p2ijk->m;
                                r2ijk->uuu = r1ij->uu * p2ijk->u;

                                kShiftUp = (kShiftMd = (kShiftDn = k + *Shift2ijk - 1) + 1) + 1;
    
                                LIndexMin =                LOutMin[iShiftDn][jShiftDn][kShiftDn] ;
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftDn][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftDn][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftMd][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftMd][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftMd][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftUp][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftUp][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftDn][jShiftUp][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftDn][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftDn][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftDn][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftMd][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftMd][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftMd][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftUp][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftUp][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftMd][jShiftUp][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftDn][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftDn][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftDn][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftMd][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftMd][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftMd][kShiftUp]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftUp][kShiftDn]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftUp][kShiftMd]);
                                LIndexMin = MAX(LIndexMin, LOutMin[iShiftUp][jShiftUp][kShiftUp]) + 1;

                                LIndexMax =                LOutMax[iShiftDn][jShiftDn][kShiftDn] ;
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftDn][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftDn][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftMd][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftMd][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftMd][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftUp][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftUp][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftDn][jShiftUp][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftDn][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftDn][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftDn][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftMd][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftMd][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftMd][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftUp][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftUp][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftMd][jShiftUp][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftDn][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftDn][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftDn][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftMd][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftMd][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftMd][kShiftUp]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftUp][kShiftDn]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftUp][kShiftMd]);
                                LIndexMax = MIN(LIndexMax, LOutMax[iShiftUp][jShiftUp][kShiftUp]) - 1;

                                Shift3ijk = dev_data->Shift3 + NodeOffset3ijk;

                                p3ijk = dev_data->p3 + NodeOffset3ijk;

                                for (L=LMinijk; L<=LMaxijk; ++L)
                                {
                                    AV30ijkL = AV30ijk + L;
                                    AV31ijkL = AV31ijk + L;
    
                                    Shift3ijkL = Shift3ijk + L;

                                    p3ijkL = &(p3ijk[L]);

                                    Drift3ijkL = Drift3ijk + D33 * L;

                                    (*UpdateDrift_3_3)(&y30, &y31, &y32,
                                                       &y30, &y31, &y32,
                                                       AV30ijkL, AV31ijkL,
                                                       &Drift3ijkL, &L, &A3, &M33,
                                                       dev_data,
                                                       tree_data,
                                                       &t);

                                    *Shift3ijkL = NEAR_INT(Drift3ijkL);

                                    if (t < T) /* VEZI */
                                    {
                                        *Shift3ijkL = MIN(MAX(LIndexMin - L, *Shift3ijkL), LIndexMax - L);
                                    }

                                    (*FinalizeDrift_3)(AV30ijkL,
                                                       &Drift3ijkL,
                                                        Shift3ijkL,
                                                       dev_data,
                                                       tree_data,
                                                       &t);

                                    Drift3ijkL -= *Shift3ijkL;

                                    CALC_PROB(p3ijkL, JumpCoeff, Drift3ijkL);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    end = clock();
    elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
    tree_data->LatticeTime += elapsed;

    return SUCCESS;
}
