/****************************************************************************/
/*      HYB4_MODEL.h                                                        */
/****************************************************************************/

#include "esl.h"

#ifndef _HYB4_MODEL_H
#define _HYB4_MODEL_H

#define HYB4_BARRIER_TOL 0.0

enum
{
    DN, MD, UP
};

enum
{
    IRF, IRD, FX, EQF, EQD
};

typedef struct
{
    double  MR;
    double  QLeft;
    double  QRight;
    double  VolBbq;
    double *FwdRateA;
    double *MLeft;
    double *MRight;
    double *SLeft;
    double *SRight;
    double *ZFRatio1;
    double *ZFRatio2;
    double *ZCenter;
    double *SpotVol;

    int Fx;

} ASSET_IR;

typedef struct
{
    double *Fwd;
    double *SpotVol;
    double *Center;

    int IrDom;

} ASSET_EQ;

typedef struct
{
    double *Fwd;
    double *SpotVol;
    double *Center;

    int IrDom;
    int IrFor;

} ASSET_FX;

typedef struct
{
    /* slices */

    double        *AssetValue0[2];
    double        *AssetValue1[2];
    double        *AssetValue2[2];
    double        *AssetValue3[2];

    double        *Aux0;
    double        *Aux1;
    double        *Aux2;
    double        *Aux3;

    int           *Shift0;
    int           *ShiftF;
    int           *Shift1;
    int           *Shift2;
    int           *Shift3;

    TPROB_0       *p0;
    TPROB_0       *pF;
    TPROB_0       *p1;
    TPROB_0       *p2;
    TPROB_0       *p3;

    TPROB_1       *r1;
    TPROB_2       *r2;

    double        *t1;
    double        *t2;
    double        *t3;

} HYB4_DEV_DATA;

typedef struct
{
    int    *NodeOffset0;
    int   **NodeOffset1;
    int  ***NodeOffset2;
    int ****NodeOffset3;

    int     NbAsset;
    int     NbAssetOn;
    int     NbTP;

    long   *TPDate;

    int     AssetType[4];
    void   *asset[4];

    int     Width[4];
    int     HalfWidth[4];

    int     xWidth[4];
    int     xHalfWidth[4];

    int     MaxIndex[4];

    int    *iMin;
    int    *iMax;
    int   **jMin;
    int   **jMax;
    int  ***kMin;
    int  ***kMax;
    int ****LMin;
    int ****LMax;

    int    *iOutMin;
    int    *iOutMax;

    int   **jOutMin;
    int   **jOutMax;

    int  ***kOutMin;
    int  ***kOutMax;

    int ****LOutMin;
    int ****LOutMax;

    double *Drift0;
    double *DriftF;
    double *Drift1;
    double *Drift2;
    double *Drift3;

    double *Length;
    double *LengthJ;
    double *JumpCoeff;

    double *Time;

    /* spatial jumps */

    double *J00;

    double *J10;
    double *J11;
    
    double *J20;
    double *J21;
    double *J22;
    
    double *J30;
    double *J31;
    double *J32;
    double *J33;

    /* correlations */

    double *R10;
    
    double *R20;
    double *R21;
    
    double *R30;
    double *R31;
    double *R32;

    double LatticeTime;

    /* moments */

    double *M0;
    double *M1;
    double *M2;
    double *M3;

    double *V00;
    double *V10;
    double *V11;
    double *V20;
    double *V21;
    double *V22;
    double *V30;
    double *V31;
    double *V32;
    double *V33;

    /* cutoff date and index for memory thrift */

    long   xDate;
    int    xT;

    /* multithreading */

    int    NbThread;

    /* empirical state probabilities */

    double *StateProb;

    /* empirical limits */

    double *LimLoMin;
    double *LimLoMax;
    double *LimHiMin;
    double *LimHiMax;

} HYB4_TREE_DATA;

typedef int (*GET_ASSET_PARAM)         (double *, double *, double *, HYB4_TREE_DATA *, int *, int *);
typedef int (*INITIALIZE_ASSET_VALUE)  (double *, double *, double *, double *, HYB4_DEV_DATA *, HYB4_TREE_DATA *, int *, int *, int *, int *);
typedef int (*UPDATE_DRIFT)            (double *, double *, double *, double *, double *, double *, double *, double *, double *, int *, double *, double *, HYB4_DEV_DATA *, HYB4_TREE_DATA *, int *);
typedef int (*FINALIZE_DRIFT)          (double *, double *, int *, HYB4_DEV_DATA *, HYB4_TREE_DATA *, int *);

typedef struct
{
    GET_ASSET_PARAM        GetAssetParam;
    INITIALIZE_ASSET_VALUE InitializeAssetValue;
    UPDATE_DRIFT           UpdateDrift[4];
    FINALIZE_DRIFT         FinalizeDrift;

} LATTICE_CODE;

typedef struct
{
    LATTICE_CODE LC[4];

} LATTICE_PROG;

#define DISC_4D_CUPS 5 /* DEV in 4-D using the currency protected foreign IR */

#endif /* HYB4_MODEL.H */
