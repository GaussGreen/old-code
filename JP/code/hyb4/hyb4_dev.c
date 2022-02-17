/*****************************************************************************/
/*        Discounted expected value.                                         */
/*****************************************************************************/
/*        HYB4_DEV.C                                                         */
/*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hyb4_lib.h"

int Hyb4_Dev(TSLICE          Price,         /* (I/O) Slice to be discounted   */
             int             t,             /* (I) Current time period        */
             int             T,             /* (I) Total number of periods    */
             int             DCurve,        /* (I) Disc curve to use (0,1,2)  */
             int             DMode,         /* (I) Nb of dims of dev required */
             HYB4_DEV_DATA  *dev_data,      /* (I) Probas, shifts, limits, etc*/
             HYB4_TREE_DATA *tree_data)     /* (I) Structure of tree data     */
{
    int ThisNodeOffset0;

    int *ThisNodeOffset1;
    int  ThisNodeOffset1i;
    
    int **ThisNodeOffset2;
    int  *ThisNodeOffset2i;
    int   ThisNodeOffset2ij;
    
    int ***ThisNodeOffset3;
    int  **ThisNodeOffset3i;
    int   *ThisNodeOffset3ij;
    int    ThisNodeOffset3ijk;
    
    int    NextNodeOffset0;
    int   *NextNodeOffset1;
    int  **NextNodeOffset2;
    int ***NextNodeOffset3;
    
    int iMin = 0;
    int iMax = 0;

    int iMinNext = 0;
    int iMaxNext = 0;

    int *jMin = NULL;
    int *jMax = NULL;
    int  jMini = 0;
    int  jMaxi = 0;

    int *jMinNext = NULL;
    int *jMaxNext = NULL;
    int  jMinNexti = 0;
    int  jMaxNexti = 0;

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

    int jOutMini = 0;
    int jOutMaxi = 0;

    int **kOutMin = NULL;
    int **kOutMax = NULL;

    int ***LOutMin = NULL;
    int ***LOutMax = NULL;

    double val;

    TPROB_0  *p0;
    TPROB_0   p0i;

    TPROB_0  *pF;
    TPROB_0   pFi;

    TPROB_0  *p1i;
    TPROB_0   p1ij;

    TPROB_0  *p2ij;
    TPROB_0   p2ijk;

    TPROB_0  *p3ijk;
    TPROB_0   p3ijkL;

    TPROB_1  *r1i;
    TPROB_1   r1ij;

    TPROB_2  *r2ij;
    TPROB_2   r2ijk;

    double *t1;
    double *t2;
    double *t3;

    int *Shift0;
    int *ShiftF;
    int *Shift1i;
    int *Shift2ij;
    int *Shift3ijk;

    double *DF;
    double *DFi;

    double DFiVal;
    double DFij;

    double *Aux0;
    double *Aux1i;
    double *Aux2ij;
    double *Aux3ijk;

    double a0;
    double a1;
    double a2;

    double a00;
    double a01;
    double a02;
    double a10;
    double a11;
    double a12;
    double a20;
    double a21;
    double a22;

    double a000;
    double a001;
    double a002;
    double a010;
    double a011;
    double a012;
    double a020;
    double a021;
    double a022;
    double a100;
    double a101;
    double a102;
    double a110;
    double a111;
    double a112;
    double a120;
    double a121;
    double a122;
    double a200;
    double a201;
    double a202;
    double a210;
    double a211;
    double a212;
    double a220;
    double a221;
    double a222;

    double *x;

    double *x0;
    double *x1;
    double *x2;

    double *x00;
    double *x01;
    double *x02;
    double *x10;
    double *x11;
    double *x12;
    double *x20;
    double *x21;
    double *x22;

    double *x000;
    double *x001;
    double *x002;
    double *x010;
    double *x011;
    double *x012;
    double *x020;
    double *x021;
    double *x022;

    double *x100;
    double *x101;
    double *x102;
    double *x110;
    double *x111;
    double *x112;
    double *x120;
    double *x121;
    double *x122;

    double *x200;
    double *x201;
    double *x202;
    double *x210;
    double *x211;
    double *x212;
    double *x220;
    double *x221;
    double *x222;

    int jLowerNext;
    int jUpperNext;

    int kLowerNext;
    int kUpperNext;

    int LLowerNext;
    int LUpperNext;

    int i, j, k, L;

    int i0, i1, i2;
    int j0, j1, j2;
    int k0, k1, k2;
    int L0, L1, L2;

    int Asset_0_On = (DMode > 0);
    int Asset_1_On = (DMode > 1);
    int Asset_2_On = (DMode > 2);
    int Asset_3_On = (DMode > 3);

    if (t >= T)
    {
        return(SUCCESS);
    }

    if (Asset_0_On)
    {
        iMin = tree_data->iMin[t];
        iMax = tree_data->iMax[t];

        iMinNext = tree_data->iMin[t+1];
        iMaxNext = tree_data->iMax[t+1];

        iOutMin = tree_data->iOutMin[t+1];
        iOutMax = tree_data->iOutMax[t+1];

        if (Asset_1_On)
        {
            jMin = tree_data->jMin[t];
            jMax = tree_data->jMax[t];

            jMinNext = tree_data->jMin[t+1];
            jMaxNext = tree_data->jMax[t+1];

            jOutMin = tree_data->jOutMin[t+1];
            jOutMax = tree_data->jOutMax[t+1];

            if (Asset_2_On)
            {
                kMin = tree_data->kMin[t];
                kMax = tree_data->kMax[t];

                kOutMin = tree_data->kOutMin[t+1];
                kOutMax = tree_data->kOutMax[t+1];

                if (Asset_3_On)
                {
                    LMin = tree_data->LMin[t];
                    LMax = tree_data->LMax[t];

                    LOutMin = tree_data->LOutMin[t+1];
                    LOutMax = tree_data->LOutMax[t+1];
                }
            }
        }
    }

    if (DMode == DISC_1D_NOCUPS)
    {
        ThisNodeOffset0 = tree_data->NodeOffset0[t];
        NextNodeOffset0 = tree_data->NodeOffset0[t+1];

        pF     = dev_data->pF                  + ThisNodeOffset0;
        ShiftF = dev_data->ShiftF              + ThisNodeOffset0;
        DF     = dev_data->AssetValue0[DCurve] + ThisNodeOffset0;
        Aux0   = dev_data->Aux0                + ThisNodeOffset0;

        x      = Price                         + NextNodeOffset0;

        val = x[iMinNext];

        for (i=iOutMin; i<iMinNext; ++i)
        {
            x[i] = val;
        }

        val = x[iMaxNext];

        for (i=iMaxNext+1; i<=iOutMax; ++i)
        {
            x[i] = val;
        }

        for (i=iMin; i<=iMax; ++i)
        {
            pFi = pF[i];

            i2 = (i1 = (i0 = i + ShiftF[i] - 1) + 1) + 1;

            x0 = x + i0;
            x1 = x + i1;
            x2 = x + i2;

            Aux0[i] = pFi.d * (*x0) +
                      pFi.m * (*x1) +
                      pFi.u * (*x2);

            Aux0[i] *= DF[i];
        }

        if (Hyb4_CopySlice(Price,
                           dev_data->Aux0,
                           1,
                           t,
                           tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }
    else if (DMode == DISC_2D_NOCUPS)
    {
        ThisNodeOffset0 = tree_data->NodeOffset0[t];
        NextNodeOffset0 = tree_data->NodeOffset0[t+1];

        ThisNodeOffset1 = tree_data->NodeOffset1[t];
        NextNodeOffset1 = tree_data->NodeOffset1[t+1];

        x = Price + NextNodeOffset1[iMinNext];

        val = x[jMaxNext[iMinNext]];

        for (i=iOutMin; i<iMinNext; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            x = Price + NextNodeOffset1[i];

            for (j=jOutMini; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        for (i=iMinNext; i<=iMaxNext; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            jMinNexti = jMinNext[i];
            jMaxNexti = jMaxNext[i];

            x = Price + NextNodeOffset1[i];

            val = x[jMinNexti];

            for (j=jOutMini; j<jMinNexti; ++j)
            {
                x[j] = val;
            }

            val = x[jMaxNexti];

            for (j=jMaxNexti+1; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        x = Price + NextNodeOffset1[iMaxNext];

        val = x[jMaxNext[iMaxNext]];

        for (i=iMaxNext+1; i<=iOutMax; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            x = Price + NextNodeOffset1[i];

            for (j=jOutMini; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        p0     = dev_data->p0                  + ThisNodeOffset0;
        Shift0 = dev_data->Shift0              + ThisNodeOffset0;
        DF     = dev_data->AssetValue0[DCurve] + ThisNodeOffset0;

        for (i=iMin; i<=iMax; ++i)
        {
            ThisNodeOffset1i = ThisNodeOffset1[i];

            Shift1i = dev_data->Shift1 + ThisNodeOffset1i;
            Aux1i   = dev_data->Aux1   + ThisNodeOffset1i;
            p1i     = dev_data->p1     + ThisNodeOffset1i;

            p0i     = p0[i];
            DFiVal  = DF[i];

            a0 = p0i.d;
            a1 = p0i.m;
            a2 = p0i.u;

            i2 = (i1 = (i0 = i + Shift0[i] - 1) + 1) + 1;

            x0 = Price + NextNodeOffset1[i0];
            x1 = Price + NextNodeOffset1[i1];
            x2 = Price + NextNodeOffset1[i2];

            jLowerNext =                 jOutMin[i0] ;
            jLowerNext = MAX(jLowerNext, jOutMin[i1]);
            jLowerNext = MAX(jLowerNext, jOutMin[i2]);

            jUpperNext =                 jOutMax[i0] ;
            jUpperNext = MIN(jUpperNext, jOutMax[i1]);
            jUpperNext = MIN(jUpperNext, jOutMax[i2]);

            t1 = dev_data->t1 - jLowerNext;

            for (j=jLowerNext; j<=jUpperNext; ++j)
            {
                t1[j] = a0 * x0[j] +
                        a1 * x1[j] +
                        a2 * x2[j];
            }

            jMini = jMin[i];
            jMaxi = jMax[i];

            for (j=jMini; j<=jMaxi; ++j)
            {
                p1ij = p1i[j];

                j2 = (j1 = (j0 = j + Shift1i[j] - 1) + 1) + 1;

                Aux1i[j] = p1ij.d * t1[j0] +
                           p1ij.m * t1[j1] +
                           p1ij.u * t1[j2];

                Aux1i[j] *= DFiVal;
            }
        }

        if (Hyb4_CopySlice(Price,
                           dev_data->Aux1,
                           2,
                           t,
                           tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }
    else if (DMode == DISC_2D_CUPS)
    {
        ThisNodeOffset0 = tree_data->NodeOffset0[t];
        NextNodeOffset0 = tree_data->NodeOffset0[t+1];

        ThisNodeOffset1 = tree_data->NodeOffset1[t];
        NextNodeOffset1 = tree_data->NodeOffset1[t+1];

        x = Price + NextNodeOffset1[iMinNext];

        val = x[jMaxNext[iMinNext]];

        for (i=iOutMin; i<iMinNext; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            x = Price + NextNodeOffset1[i];

            for (j=jOutMini; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        for (i=iMinNext; i<=iMaxNext; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            jMinNexti = jMinNext[i];
            jMaxNexti = jMaxNext[i];

            x = Price + NextNodeOffset1[i];

            val = x[jMinNexti];

            for (j=jOutMini; j<jMinNexti; ++j)
            {
                x[j] = val;
            }

            val = x[jMaxNexti];

            for (j=jMaxNexti+1; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        x = Price + NextNodeOffset1[iMaxNext];

        val = x[jMaxNext[iMaxNext]];

        for (i=iMaxNext+1; i<=iOutMax; ++i)
        {
            jOutMini = jOutMin[i];
            jOutMaxi = jOutMax[i];

            x = Price + NextNodeOffset1[i];

            for (j=jOutMini; j<=jOutMaxi; ++j)
            {
                x[j] = val;
            }
        }

        p0     = dev_data->p0     + ThisNodeOffset0;
        Shift0 = dev_data->Shift0 + ThisNodeOffset0;

        for (i=iMin; i<=iMax; ++i)
        {
            ThisNodeOffset1i = ThisNodeOffset1[i];

            DFi     = dev_data->AssetValue1[DCurve] + ThisNodeOffset1i;
            Shift1i = dev_data->Shift1              + ThisNodeOffset1i;
            Aux1i   = dev_data->Aux1                + ThisNodeOffset1i;
            p1i     = dev_data->p1                  + ThisNodeOffset1i;

            p0i = p0[i];

            a0 = p0i.d;
            a1 = p0i.m;
            a2 = p0i.u;

            i2 = (i1 = (i0 = i + Shift0[i] - 1) + 1) + 1;

            x0 = Price + NextNodeOffset1[i0];
            x1 = Price + NextNodeOffset1[i1];
            x2 = Price + NextNodeOffset1[i2];

            jLowerNext =                 jOutMin[i0] ;
            jLowerNext = MAX(jLowerNext, jOutMin[i1]);
            jLowerNext = MAX(jLowerNext, jOutMin[i2]);

            jUpperNext =                 jOutMax[i0] ;
            jUpperNext = MIN(jUpperNext, jOutMax[i1]);
            jUpperNext = MIN(jUpperNext, jOutMax[i2]);

            t1 = dev_data->t1 - jLowerNext;

            for (j=jLowerNext; j<=jUpperNext; ++j)
            {
                t1[j] = a0 * x0[j] +
                        a1 * x1[j] +
                        a2 * x2[j];
            }

            jMini = jMin[i];
            jMaxi = jMax[i];

            for (j=jMini; j<=jMaxi; ++j)
            {
                p1ij = p1i[j];

                j2 = (j1 = (j0 = j + Shift1i[j] - 1) + 1) + 1;

                Aux1i[j] = p1ij.d * t1[j0] +
                           p1ij.m * t1[j1] +
                           p1ij.u * t1[j2];

                Aux1i[j] *= DFi[j];
            }
        }

        if (Hyb4_CopySlice(Price,
                           dev_data->Aux1,
                           2,
                           t,
                           tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }
    else if (DMode == DISC_3D_CUPS)
    {
        ThisNodeOffset0 = tree_data->NodeOffset0[t];
        NextNodeOffset0 = tree_data->NodeOffset0[t+1];

        ThisNodeOffset1 = tree_data->NodeOffset1[t];
        NextNodeOffset1 = tree_data->NodeOffset1[t+1];

        ThisNodeOffset2 = tree_data->NodeOffset2[t];
        NextNodeOffset2 = tree_data->NodeOffset2[t+1];

        /* we don't fill in the borders */

        Shift0 = dev_data->Shift0 + ThisNodeOffset0;

        for (i=iMin; i<=iMax; ++i)
        {
            ThisNodeOffset1i = ThisNodeOffset1[i];
            ThisNodeOffset2i = ThisNodeOffset2[i];

            DFi     = dev_data->AssetValue1[DCurve] + ThisNodeOffset1i;
            Shift1i = dev_data->Shift1              + ThisNodeOffset1i;
            r1i     = dev_data->r1                  + ThisNodeOffset1i;

            i2 = (i1 = (i0 = i + Shift0[i] - 1) + 1) + 1;

            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            for (j=jMini; j<=jMaxi; ++j)
            {
                ThisNodeOffset2ij = ThisNodeOffset2i[j];

                Shift2ij = dev_data->Shift2 + ThisNodeOffset2ij;
                p2ij     = dev_data->p2     + ThisNodeOffset2ij;
                Aux2ij   = dev_data->Aux2   + ThisNodeOffset2ij;

                DFij = DFi[j];

                r1ij = r1i[j];

                a00 = r1ij.dd;
                a01 = r1ij.dm;
                a02 = r1ij.du;

                a10 = r1ij.md;
                a11 = r1ij.mm;
                a12 = r1ij.mu;

                a20 = r1ij.ud;
                a21 = r1ij.um;
                a22 = r1ij.uu;

                j2 = (j1 = (j0 = j + Shift1i[j] - 1) + 1) + 1;

                x00 = Price + NextNodeOffset2[i0][j0];
                x01 = Price + NextNodeOffset2[i0][j1];
                x02 = Price + NextNodeOffset2[i0][j2];
                x10 = Price + NextNodeOffset2[i1][j0];
                x11 = Price + NextNodeOffset2[i1][j1];
                x12 = Price + NextNodeOffset2[i1][j2];
                x20 = Price + NextNodeOffset2[i2][j0];
                x21 = Price + NextNodeOffset2[i2][j1];
                x22 = Price + NextNodeOffset2[i2][j2];

                kLowerNext =                 kOutMin[i0][j0] ;
                kLowerNext = MAX(kLowerNext, kOutMin[i0][j1]);
                kLowerNext = MAX(kLowerNext, kOutMin[i0][j2]);
                kLowerNext = MAX(kLowerNext, kOutMin[i1][j0]);
                kLowerNext = MAX(kLowerNext, kOutMin[i1][j1]);
                kLowerNext = MAX(kLowerNext, kOutMin[i1][j2]);
                kLowerNext = MAX(kLowerNext, kOutMin[i2][j0]);
                kLowerNext = MAX(kLowerNext, kOutMin[i2][j1]);
                kLowerNext = MAX(kLowerNext, kOutMin[i2][j2]);

                kUpperNext =                 kOutMax[i0][j0] ;
                kUpperNext = MIN(kUpperNext, kOutMax[i0][j1]);
                kUpperNext = MIN(kUpperNext, kOutMax[i0][j2]);
                kUpperNext = MIN(kUpperNext, kOutMax[i1][j0]);
                kUpperNext = MIN(kUpperNext, kOutMax[i1][j1]);
                kUpperNext = MIN(kUpperNext, kOutMax[i1][j2]);
                kUpperNext = MIN(kUpperNext, kOutMax[i2][j0]);
                kUpperNext = MIN(kUpperNext, kOutMax[i2][j1]);
                kUpperNext = MIN(kUpperNext, kOutMax[i2][j2]);

                t2 = dev_data->t2 - kLowerNext;

                for (k=kLowerNext; k<=kUpperNext; ++k)
                {
                    t2[k] = a00 * x00[k] +
                            a01 * x01[k] +
                            a02 * x02[k] +
                            a10 * x10[k] +
                            a11 * x11[k] +
                            a12 * x12[k] +
                            a20 * x20[k] +
                            a21 * x21[k] +
                            a22 * x22[k];
                }

                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                for (k=kMinij; k<=kMaxij; ++k)
                {
                    p2ijk = p2ij[k];

                    k2 = (k1 = (k0 = k + Shift2ij[k] - 1) + 1) + 1;

                    Aux2ij[k] = p2ijk.d * t2[k0] +
                                p2ijk.m * t2[k1] +
                                p2ijk.u * t2[k2];

                    Aux2ij[k] *= DFij;
                }
            }
        }

        if (Hyb4_CopySlice(Price,
                           dev_data->Aux2,
                           3,
                           t,
                           tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }
    else if (DMode == DISC_4D_CUPS)
    {
        ThisNodeOffset0 = tree_data->NodeOffset0[t];
        NextNodeOffset0 = tree_data->NodeOffset0[t+1];

        ThisNodeOffset1 = tree_data->NodeOffset1[t];
        NextNodeOffset1 = tree_data->NodeOffset1[t+1];

        ThisNodeOffset2 = tree_data->NodeOffset2[t];
        NextNodeOffset2 = tree_data->NodeOffset2[t+1];

        ThisNodeOffset3 = tree_data->NodeOffset3[t];
        NextNodeOffset3 = tree_data->NodeOffset3[t+1];

        /* we don't fill in the borders */

        Shift0 = dev_data->Shift0 + ThisNodeOffset0;

        for (i=iMin; i<=iMax; ++i)
        {
            ThisNodeOffset1i = ThisNodeOffset1[i];
            ThisNodeOffset2i = ThisNodeOffset2[i];
            ThisNodeOffset3i = ThisNodeOffset3[i];

            DFi     = dev_data->AssetValue1[DCurve] + ThisNodeOffset1i;
            Shift1i = dev_data->Shift1              + ThisNodeOffset1i;

            i2 = (i1 = (i0 = i + Shift0[i] - 1) + 1) + 1;

            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            LMini = LMin[i];
            LMaxi = LMax[i];

            for (j=jMini; j<=jMaxi; ++j)
            {
                ThisNodeOffset2ij = ThisNodeOffset2i[j];
                ThisNodeOffset3ij = ThisNodeOffset3i[j];

                DFij = DFi[j];

                Shift2ij = dev_data->Shift2 + ThisNodeOffset2ij;
                r2ij     = dev_data->r2     + ThisNodeOffset2ij;

                j2 = (j1 = (j0 = j + Shift1i[j] - 1) + 1) + 1;

                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                LMinij = LMini[j];
                LMaxij = LMaxi[j];

                for (k=kMinij; k<=kMaxij; ++k)
                {
                    ThisNodeOffset3ijk = ThisNodeOffset3ij[k];

                    Shift3ijk = dev_data->Shift3 + ThisNodeOffset3ijk;
                    p3ijk     = dev_data->p3     + ThisNodeOffset3ijk;
                    Aux3ijk   = dev_data->Aux3   + ThisNodeOffset3ijk;

                    r2ijk = r2ij[k];

                    a000 = r2ijk.ddd;
                    a001 = r2ijk.ddm;
                    a002 = r2ijk.ddu;
                    a010 = r2ijk.dmd;
                    a011 = r2ijk.dmm;
                    a012 = r2ijk.dmu;
                    a020 = r2ijk.dud;
                    a021 = r2ijk.dum;
                    a022 = r2ijk.duu;

                    a100 = r2ijk.mdd;
                    a101 = r2ijk.mdm;
                    a102 = r2ijk.mdu;
                    a110 = r2ijk.mmd;
                    a111 = r2ijk.mmm;
                    a112 = r2ijk.mmu;
                    a120 = r2ijk.mud;
                    a121 = r2ijk.mum;
                    a122 = r2ijk.muu;

                    a200 = r2ijk.udd;
                    a201 = r2ijk.udm;
                    a202 = r2ijk.udu;
                    a210 = r2ijk.umd;
                    a211 = r2ijk.umm;
                    a212 = r2ijk.umu;
                    a220 = r2ijk.uud;
                    a221 = r2ijk.uum;
                    a222 = r2ijk.uuu;

                    k2 = (k1 = (k0 = k + Shift2ij[k] - 1) + 1) + 1;

                    x000 = Price + NextNodeOffset3[i0][j0][k0];
                    x001 = Price + NextNodeOffset3[i0][j0][k1];
                    x002 = Price + NextNodeOffset3[i0][j0][k2];
                    x010 = Price + NextNodeOffset3[i0][j1][k0];
                    x011 = Price + NextNodeOffset3[i0][j1][k1];
                    x012 = Price + NextNodeOffset3[i0][j1][k2];
                    x020 = Price + NextNodeOffset3[i0][j2][k0];
                    x021 = Price + NextNodeOffset3[i0][j2][k1];
                    x022 = Price + NextNodeOffset3[i0][j2][k2];

                    x100 = Price + NextNodeOffset3[i1][j0][k0];
                    x101 = Price + NextNodeOffset3[i1][j0][k1];
                    x102 = Price + NextNodeOffset3[i1][j0][k2];
                    x110 = Price + NextNodeOffset3[i1][j1][k0];
                    x111 = Price + NextNodeOffset3[i1][j1][k1];
                    x112 = Price + NextNodeOffset3[i1][j1][k2];
                    x120 = Price + NextNodeOffset3[i1][j2][k0];
                    x121 = Price + NextNodeOffset3[i1][j2][k1];
                    x122 = Price + NextNodeOffset3[i1][j2][k2];

                    x200 = Price + NextNodeOffset3[i2][j0][k0];
                    x201 = Price + NextNodeOffset3[i2][j0][k1];
                    x202 = Price + NextNodeOffset3[i2][j0][k2];
                    x210 = Price + NextNodeOffset3[i2][j1][k0];
                    x211 = Price + NextNodeOffset3[i2][j1][k1];
                    x212 = Price + NextNodeOffset3[i2][j1][k2];
                    x220 = Price + NextNodeOffset3[i2][j2][k0];
                    x221 = Price + NextNodeOffset3[i2][j2][k1];
                    x222 = Price + NextNodeOffset3[i2][j2][k2];

                    LLowerNext =                 LOutMin[i0][j0][k0] ;
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j0][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j0][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j1][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j1][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j1][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j2][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j2][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i0][j2][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j0][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j0][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j0][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j1][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j1][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j1][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j2][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j2][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i1][j2][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j0][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j0][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j0][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j1][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j1][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j1][k2]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j2][k0]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j2][k1]);
                    LLowerNext = MAX(LLowerNext, LOutMin[i2][j2][k2]);

                    LUpperNext =                 LOutMax[i0][j0][k0] ;
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j0][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j0][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j1][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j1][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j1][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j2][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j2][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i0][j2][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j0][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j0][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j0][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j1][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j1][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j1][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j2][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j2][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i1][j2][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j0][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j0][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j0][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j1][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j1][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j1][k2]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j2][k0]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j2][k1]);
                    LUpperNext = MIN(LUpperNext, LOutMax[i2][j2][k2]);

                    t3 = dev_data->t3 - LLowerNext;

                    for (L=LLowerNext; L<=LUpperNext; ++L)
                    {
                        t3[L] = a000 * x000[L] +
                                a001 * x001[L] +
                                a002 * x002[L] +
                                a010 * x010[L] +
                                a011 * x011[L] +
                                a012 * x012[L] +
                                a020 * x020[L] +
                                a021 * x021[L] +
                                a022 * x022[L] +
                                a100 * x100[L] +
                                a101 * x101[L] +
                                a102 * x102[L] +
                                a110 * x110[L] +
                                a111 * x111[L] +
                                a112 * x112[L] +
                                a120 * x120[L] +
                                a121 * x121[L] +
                                a122 * x122[L] +
                                a200 * x200[L] +
                                a201 * x201[L] +
                                a202 * x202[L] +
                                a210 * x210[L] +
                                a211 * x211[L] +
                                a212 * x212[L] +
                                a220 * x220[L] +
                                a221 * x221[L] +
                                a222 * x222[L];
                    }

                    LMinijk = LMinij[k];
                    LMaxijk = LMaxij[k];

                    for (L=LMinijk; L<=LMaxijk; ++L)
                    {
                        p3ijkL = p3ijk[L];

                        L2 = (L1 = (L0 = L + Shift3ijk[L] - 1) + 1) + 1;

                        Aux3ijk[L] = p3ijkL.d * t3[L0] +
                                     p3ijkL.m * t3[L1] +
                                     p3ijkL.u * t3[L2];

                        Aux3ijk[L] *= DFij;
                    }
                }
            }
        }

        if (Hyb4_CopySlice(Price,
                           dev_data->Aux3,
                           4,
                           t,
                           tree_data) == FAILURE)
        {
            return FAILURE;
        }
    }

    return SUCCESS;
}
