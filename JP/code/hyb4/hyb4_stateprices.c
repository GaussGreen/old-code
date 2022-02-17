/****************************************************************************/
/*      Calculates state-prices in the tree                                 */
/****************************************************************************/
/*      HYB4_STATEPRICES.C                                                  */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hyb4_lib.h"

int Hyb4_FillGrid_Asset3(HYB4_DEV_DATA  *dev_data,
                         HYB4_TREE_DATA *tree_data,
                         LATTICE_PROG   *LP,
                         int             t)
{
    int    status = FAILURE;

    int ***NodeOffset3 = tree_data->NodeOffset3[t];
    int  **NodeOffset3i;
    int   *NodeOffset3ij;
    int    NodeOffset3ijk;

    int    i, j, k;

    double PJ30 = tree_data->J30[t-1];
    double PJ31 = tree_data->J31[t-1];
    double PJ32 = tree_data->J32[t-1];
    double PJ33 = tree_data->J33[t-1];

    double X3 = 0;

    double X3i;
    double X3ij;
    double X3ijk;

    TSLICE AV30ijk;
    TSLICE AV31ijk;

    int    iMin = tree_data->iMin[t];
    int    iMax = tree_data->iMax[t];

    int   *jMin = tree_data->jMin[t];
    int   *jMax = tree_data->jMax[t];

    int  **kMin = tree_data->kMin[t];
    int  **kMax = tree_data->kMax[t];

    int ***LMin = tree_data->LMin[t];
    int ***LMax = tree_data->LMax[t];

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

    static int I3 = 3;

    INITIALIZE_ASSET_VALUE InitializeAssetValue_3 =
                               LP->LC[3].InitializeAssetValue;

    for (i=iMin; i<=iMax; ++i)
    {
        X3i = X3 + PJ30 * i;

        jMini = jMin[i];
        jMaxi = jMax[i];

        kMini = kMin[i];
        kMaxi = kMax[i];

        LMini = LMin[i];
        LMaxi = LMax[i];

        NodeOffset3i = NodeOffset3[i];

        for (j=jMini; j<=jMaxi; ++j)
        {
            X3ij = X3i + PJ31 * j;

            kMinij = kMini[j];
            kMaxij = kMaxi[j];

            LMinij = LMini[j];
            LMaxij = LMaxi[j];
        
            NodeOffset3ij = NodeOffset3i[j];

            for (k=kMinij; k<=kMaxij; ++k)
            {
                X3ijk = X3ij + PJ32 * k;

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
            }
        }
    }

    status = SUCCESS;

    return status;
}

/*****  Hyb4_UpdateStatePrices4D  ***********************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*/
int Hyb4_UpdateStatePrices4D
    (int              DfDim,                  /* (I) discount dim (-1 for probs) */
     int              t,                      /* (I) current time period         */  
     HYB4_TREE_DATA  *tree_data,              /* (I) Tree data                   */
     HYB4_DEV_DATA   *dev_data,               /* (I/O) dev-data for tree         */
     LATTICE_PROG    *LP,
     TSLICE           StatePr0,               /* (I) initial state price         */
     TSLICE           StatePr1)               /* (O) final State price           */
{
    int      status = FAILURE;

    int DfDim_0 = (DfDim == 0);
    int DfDim_1 = (DfDim == 1);
    int DfDim_2 = (DfDim == 2);
    int DfDim_3 = (DfDim == 3);

    double   Df = 1;

    double   x;

    TPROB_0   *pijk;
    TPROB_0    pijkL;

    TPROB_2   *rij;
    TPROB_2    rijk;

    double   p0x, p1x, p2x;

    double   r000;
    double   r001;
    double   r002;
    double   r010;
    double   r011;
    double   r012;
    double   r020;
    double   r021;
    double   r022;
    double   r100;
    double   r101;
    double   r102;
    double   r110;
    double   r111;
    double   r112;
    double   r120;
    double   r121;
    double   r122;
    double   r200;
    double   r201;
    double   r202;
    double   r210;
    double   r211;
    double   r212;
    double   r220;
    double   r221;
    double   r222;

    double  *StatePr0L = NULL;
    double  *StatePr1L = NULL;

    double  *StatePr000, *StatePr001, *StatePr002;
    double  *StatePr010, *StatePr011, *StatePr012;
    double  *StatePr020, *StatePr021, *StatePr022;

    double  *StatePr100, *StatePr101, *StatePr102;
    double  *StatePr110, *StatePr111, *StatePr112;
    double  *StatePr120, *StatePr121, *StatePr122;

    double  *StatePr200, *StatePr201, *StatePr202;
    double  *StatePr210, *StatePr211, *StatePr212;
    double  *StatePr220, *StatePr221, *StatePr222;

    int      NodeOffset0 = tree_data->NodeOffset0[t];
    int     *NodeOffset1 = tree_data->NodeOffset1[t];
    int    **NodeOffset2 = tree_data->NodeOffset2[t];
    int   ***NodeOffset3 = tree_data->NodeOffset3[t];

    int      NodeOffset1i;
    int     *NodeOffset2i;
    int    **NodeOffset3i; 

    int      NodeOffset2ij;
    int     *NodeOffset3ij;

    int      NodeOffset3ijk;

    int   ***NodeOffset3Next = tree_data->NodeOffset3[t+1];
    int    **NodeOffset3Nexti; 
    int     *NodeOffset3Nextij;
    int      NodeOffset3Nextijk;

    int      iMin = tree_data->iMin[t];
    int      iMax = tree_data->iMax[t];

    int     *jMin = tree_data->jMin[t];
    int     *jMax = tree_data->jMax[t];

    int    **kMin = tree_data->kMin[t];
    int    **kMax = tree_data->kMax[t];

    int   ***LMin = tree_data->LMin[t];
    int   ***LMax = tree_data->LMax[t];

    int      jMini;
    int      jMaxi;
    int     *kMini;
    int     *kMaxi;
    int    **LMini;
    int    **LMaxi;

    int      kMinij;
    int      kMaxij;
    int     *LMinij;
    int     *LMaxij;

    int      LMinijk;
    int      LMaxijk;

    int      iMinNext = tree_data->iMin[t+1];
    int      iMaxNext = tree_data->iMax[t+1];

    int     *jMinNext = tree_data->jMin[t+1];
    int     *jMaxNext = tree_data->jMax[t+1];

    int    **kMinNext = tree_data->kMin[t+1];
    int    **kMaxNext = tree_data->kMax[t+1];

    int   ***LMinNext = tree_data->LMin[t+1];
    int   ***LMaxNext = tree_data->LMax[t+1];

    int      jMinNexti;
    int      jMaxNexti;
    int     *kMinNexti;
    int     *kMaxNexti;
    int    **LMinNexti;
    int    **LMaxNexti;

    int      kMinNextij;
    int      kMaxNextij;
    int     *LMinNextij;
    int     *LMaxNextij;

    int      LMinNextijk;
    int      LMaxNextijk;

    int      i, j, k, L;

    int      i0, i1, i2,
             j0, j1, j2,
             k0, k1, k2,
             L0, L1, L2;
            
    int     *Shift0, *Shift1i, *Shift2ij, *Shift3ijk;

    int      T  = tree_data->NbTP;
    int      xT = tree_data->xT;

    if (DfDim_3)
    {
        DR_Error("DfDim must not be 3.\n");
        goto RETURN;
    }

    if (Hyb4_Lattice(dev_data,
                     tree_data,
                     LP,
                     t,
                     T) == FAILURE)
    {
        goto RETURN;
    }

    if (t == 0)
    {
        (StatePr0 + NodeOffset3[0][0][0])[0] = 1;
    }
    
    if (t != xT)
    {
        for (i = iMinNext; i <= iMaxNext; i++)
        {
            jMinNexti = jMinNext[i];
            jMaxNexti = jMaxNext[i];

            kMinNexti = kMinNext[i];
            kMaxNexti = kMaxNext[i];

            LMinNexti = LMinNext[i];
            LMaxNexti = LMaxNext[i];

            NodeOffset3Nexti = NodeOffset3Next[i];

            for (j = jMinNexti; j <= jMaxNexti; j++)
            {                                   
                kMinNextij = kMinNexti[j];
                kMaxNextij = kMaxNexti[j];

                LMinNextij = LMinNexti[j];
                LMaxNextij = LMaxNexti[j];

                NodeOffset3Nextij = NodeOffset3Nexti[j];

                for (k = kMinNextij; k <= kMaxNextij; k++)
                {
                    LMinNextijk = LMinNextij[k];
                    LMaxNextijk = LMaxNextij[k];

                    NodeOffset3Nextijk = NodeOffset3Nextij[k];

                    StatePr1L = StatePr1 + NodeOffset3Nextijk;

                    for (L = LMinNextijk; L <= LMaxNextijk; L++)
                    {
                        StatePr1L[L] = 0.;
                    }
                }
            }
        }

        Shift0 = dev_data->Shift0 + NodeOffset0;

        for (i = iMin; i <= iMax; i++)
        {
            i2 = (i1 = (i0 = i + Shift0[i] - 1) + 1) + 1;

            jMini = jMin[i];
            jMaxi = jMax[i];

            kMini = kMin[i];
            kMaxi = kMax[i];

            LMini = LMin[i];
            LMaxi = LMax[i];

            NodeOffset1i = NodeOffset1[i];
            NodeOffset2i = NodeOffset2[i];
            NodeOffset3i = NodeOffset3[i];

            Shift1i = dev_data->Shift1 + NodeOffset1i;

            if (DfDim_0)
            {
                Df = (dev_data->AssetValue0[1] + NodeOffset0)[i];
            }

            for (j = jMini; j <= jMaxi; j++)
            {                                   
                j2 = (j1 = (j0 = j + Shift1i[j] - 1) + 1) + 1;

                kMinij = kMini[j];
                kMaxij = kMaxi[j];

                LMinij = LMini[j];
                LMaxij = LMaxi[j];

                NodeOffset2ij = NodeOffset2i[j];
                NodeOffset3ij = NodeOffset3i[j];

                Shift2ij = dev_data->Shift2 + NodeOffset2ij;

                rij      = dev_data->r2     + NodeOffset2ij;

                if (DfDim_1)
                {
                    Df = (dev_data->AssetValue1[1] + NodeOffset1i)[j];
                }

                for (k = kMinij; k <= kMaxij; k++)
                {
                    k2 = (k1 = (k0 = k + Shift2ij[k] - 1) + 1) + 1;

                    LMinijk = LMinij[k];
                    LMaxijk = LMaxij[k];

                    NodeOffset3ijk = NodeOffset3ij[k];

                    rijk = rij[k];

                    r000 = rijk.ddd;
                    r001 = rijk.ddm;
                    r002 = rijk.ddu;
                    r010 = rijk.dmd;
                    r011 = rijk.dmm;
                    r012 = rijk.dmu;
                    r020 = rijk.dud;
                    r021 = rijk.dum;
                    r022 = rijk.duu;

                    r100 = rijk.mdd;
                    r101 = rijk.mdm;
                    r102 = rijk.mdu;
                    r110 = rijk.mmd;
                    r111 = rijk.mmm;
                    r112 = rijk.mmu;
                    r120 = rijk.mud;
                    r121 = rijk.mum;
                    r122 = rijk.muu;

                    r200 = rijk.udd;
                    r201 = rijk.udm;
                    r202 = rijk.udu;
                    r210 = rijk.umd;
                    r211 = rijk.umm;
                    r212 = rijk.umu;
                    r220 = rijk.uud;
                    r221 = rijk.uum;
                    r222 = rijk.uuu;

                    Shift3ijk = dev_data->Shift3 + NodeOffset3ijk;

                    pijk      = dev_data->p3     + NodeOffset3ijk;

                    StatePr0L = StatePr0         + NodeOffset3ijk;

                    StatePr000 = StatePr1 + NodeOffset3Next[i0][j0][k0];
                    StatePr001 = StatePr1 + NodeOffset3Next[i0][j0][k1];
                    StatePr002 = StatePr1 + NodeOffset3Next[i0][j0][k2];
                    StatePr010 = StatePr1 + NodeOffset3Next[i0][j1][k0];
                    StatePr011 = StatePr1 + NodeOffset3Next[i0][j1][k1];
                    StatePr012 = StatePr1 + NodeOffset3Next[i0][j1][k2];
                    StatePr020 = StatePr1 + NodeOffset3Next[i0][j2][k0];
                    StatePr021 = StatePr1 + NodeOffset3Next[i0][j2][k1];
                    StatePr022 = StatePr1 + NodeOffset3Next[i0][j2][k2];
                    StatePr100 = StatePr1 + NodeOffset3Next[i1][j0][k0];
                    StatePr101 = StatePr1 + NodeOffset3Next[i1][j0][k1];
                    StatePr102 = StatePr1 + NodeOffset3Next[i1][j0][k2];
                    StatePr110 = StatePr1 + NodeOffset3Next[i1][j1][k0];
                    StatePr111 = StatePr1 + NodeOffset3Next[i1][j1][k1];
                    StatePr112 = StatePr1 + NodeOffset3Next[i1][j1][k2];
                    StatePr120 = StatePr1 + NodeOffset3Next[i1][j2][k0];
                    StatePr121 = StatePr1 + NodeOffset3Next[i1][j2][k1];
                    StatePr122 = StatePr1 + NodeOffset3Next[i1][j2][k2];
                    StatePr200 = StatePr1 + NodeOffset3Next[i2][j0][k0];
                    StatePr201 = StatePr1 + NodeOffset3Next[i2][j0][k1];
                    StatePr202 = StatePr1 + NodeOffset3Next[i2][j0][k2];
                    StatePr210 = StatePr1 + NodeOffset3Next[i2][j1][k0];
                    StatePr211 = StatePr1 + NodeOffset3Next[i2][j1][k1];
                    StatePr212 = StatePr1 + NodeOffset3Next[i2][j1][k2];
                    StatePr220 = StatePr1 + NodeOffset3Next[i2][j2][k0];
                    StatePr221 = StatePr1 + NodeOffset3Next[i2][j2][k1];
                    StatePr222 = StatePr1 + NodeOffset3Next[i2][j2][k2];
                
                    if (DfDim_2)
                    {
                        Df = (dev_data->AssetValue2[1] + NodeOffset2ij)[k];
                    }

                    for (L = LMinijk; L <= LMaxijk; L++)
                    {
                        L2 = (L1 = (L0 = L + Shift3ijk[L] - 1) + 1) + 1;

                        x = StatePr0L[L];

                        pijkL = pijk[L];

                        p0x = pijkL.d * x * Df; 
                        p1x = pijkL.m * x * Df; 
                        p2x = pijkL.u * x * Df;
                    
                        StatePr000[L0] += r000 * p0x;
                        StatePr001[L0] += r001 * p0x;
                        StatePr002[L0] += r002 * p0x;
                        StatePr010[L0] += r010 * p0x;
                        StatePr011[L0] += r011 * p0x;
                        StatePr012[L0] += r012 * p0x;
                        StatePr020[L0] += r020 * p0x;
                        StatePr021[L0] += r021 * p0x;
                        StatePr022[L0] += r022 * p0x;

                        StatePr100[L0] += r100 * p0x;
                        StatePr101[L0] += r101 * p0x;
                        StatePr102[L0] += r102 * p0x;
                        StatePr110[L0] += r110 * p0x;
                        StatePr111[L0] += r111 * p0x;
                        StatePr112[L0] += r112 * p0x;
                        StatePr120[L0] += r120 * p0x;
                        StatePr121[L0] += r121 * p0x;
                        StatePr122[L0] += r122 * p0x;

                        StatePr200[L0] += r200 * p0x;
                        StatePr201[L0] += r201 * p0x;
                        StatePr202[L0] += r202 * p0x;
                        StatePr210[L0] += r210 * p0x;
                        StatePr211[L0] += r211 * p0x;
                        StatePr212[L0] += r212 * p0x;
                        StatePr220[L0] += r220 * p0x;
                        StatePr221[L0] += r221 * p0x;
                        StatePr222[L0] += r222 * p0x;

                        StatePr000[L1] += r000 * p1x;
                        StatePr001[L1] += r001 * p1x;
                        StatePr002[L1] += r002 * p1x;
                        StatePr010[L1] += r010 * p1x;
                        StatePr011[L1] += r011 * p1x;
                        StatePr012[L1] += r012 * p1x;
                        StatePr020[L1] += r020 * p1x;
                        StatePr021[L1] += r021 * p1x;
                        StatePr022[L1] += r022 * p1x;

                        StatePr100[L1] += r100 * p1x;
                        StatePr101[L1] += r101 * p1x;
                        StatePr102[L1] += r102 * p1x;
                        StatePr110[L1] += r110 * p1x;
                        StatePr111[L1] += r111 * p1x;
                        StatePr112[L1] += r112 * p1x;
                        StatePr120[L1] += r120 * p1x;
                        StatePr121[L1] += r121 * p1x;
                        StatePr122[L1] += r122 * p1x;

                        StatePr200[L1] += r200 * p1x;
                        StatePr201[L1] += r201 * p1x;
                        StatePr202[L1] += r202 * p1x;
                        StatePr210[L1] += r210 * p1x;
                        StatePr211[L1] += r211 * p1x;
                        StatePr212[L1] += r212 * p1x;
                        StatePr220[L1] += r220 * p1x;
                        StatePr221[L1] += r221 * p1x;
                        StatePr222[L1] += r222 * p1x;

                        StatePr000[L2] += r000 * p2x;
                        StatePr001[L2] += r001 * p2x;
                        StatePr002[L2] += r002 * p2x;
                        StatePr010[L2] += r010 * p2x;
                        StatePr011[L2] += r011 * p2x;
                        StatePr012[L2] += r012 * p2x;
                        StatePr020[L2] += r020 * p2x;
                        StatePr021[L2] += r021 * p2x;
                        StatePr022[L2] += r022 * p2x;

                        StatePr100[L2] += r100 * p2x;
                        StatePr101[L2] += r101 * p2x;
                        StatePr102[L2] += r102 * p2x;
                        StatePr110[L2] += r110 * p2x;
                        StatePr111[L2] += r111 * p2x;
                        StatePr112[L2] += r112 * p2x;
                        StatePr120[L2] += r120 * p2x;
                        StatePr121[L2] += r121 * p2x;
                        StatePr122[L2] += r122 * p2x;

                        StatePr200[L2] += r200 * p2x;
                        StatePr201[L2] += r201 * p2x;
                        StatePr202[L2] += r202 * p2x;
                        StatePr210[L2] += r210 * p2x;
                        StatePr211[L2] += r211 * p2x;
                        StatePr212[L2] += r212 * p2x;
                        StatePr220[L2] += r220 * p2x;
                        StatePr221[L2] += r221 * p2x;
                        StatePr222[L2] += r222 * p2x;
                    }
                }
            }
        }
    }

    status = SUCCESS;
    
RETURN:

    return status;
}
