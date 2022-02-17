/****************************************************************************/
/*      Calculates state-prices in the tree                                 */
/****************************************************************************/
/*      STATEPRICES.C                                                       */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"



/*****  Hyb3_BuildStatePrices3D  ***********************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*/
int     Hyb3_UpdateStatePrices3D
    (int              t,                      /* (I) current time period   */  
     MKTVOL_DATA     *mktvol_data,            /* (I) Market vol data       */
     HYB3_TREE_DATA  *tree_data,              /* (I) Tree data             */
     HYB3_DEV_DATA   *dev_data,               /* (I/O) dev-data for tree   */
     TSLICE           StatePr,                /* (I) initial state price   */
     TSLICE           StatePr1)               /* (O) final State price     */
{

    double  *StatePrL  = NULL;      /* pointers state price slices      */
    double  *StatePr1L = NULL;

    double  *StatePr11L, *StatePr12L, *StatePr13L;
    double  *StatePr21L, *StatePr22L, *StatePr23L;
    double  *StatePr31L, *StatePr32L, *StatePr33L;

    double  *Discount2DL;           /* pointer to domestic discount     */
    
    double  *quu,  *qu0,  *qud;     /*  2-D probability slices          */
    double  *q0u,  *q00,  *q0d;
    double  *qdu,  *qd0,  *qdd;

    TPROB_0 *r;

    /* pointer to equity or fx structure */
    double      *ExSpot        = NULL;
    double      *ExMidNode     = NULL;
    double      *ExVol         = NULL;
    double      *FwdEx         = NULL;


    int
       *Top1,     *Bottom1,         /* Limits of the tree (1rst dimension)  */
      **Top2,    **Bottom2,         /* Limits of the tree (2nd dimension)   */
     ***Top3,   ***Bottom3;         /* Limits of the tree (3rd dimension)   */

    int     *Shift1, *Shift2, *Shift3;

    /* local variables for easy referencing */
    double    Quu,   Qu0,  Qud;         /* 2-D probabilities (values)   */
    double    Q0u,   Q00,  Q0d;
    double    Qdu,   Qd0,  Qdd;

    double    Rux,   R0x,  Rdx;            /* 3-D probabilities (values)   */

    int       i0, i1, i2,
              j0, j1, j2,
              k0, k1, k2;
            
    double   Discount_ij,               /* domestic discount            */
             x;

    int      i, j, k;

    int     DCurveD,           /* domestic discount curve */

            T,           /* Total number of period in the hybrids tree*/
            offset,
            status = FAILURE; /* Error status = FAILURE initially          */



    /* Total size of tree timeline */     
    T   = tree_data->NbTP;


    /* Assigment of domestic discount curve */
    DCurveD   = tree_data->CvDisc[1];

       
    Top1    = tree_data->Top1;	        
    Top2    = tree_data->Top2;
    Top3    = tree_data->Top3;	    
    Bottom1 = tree_data->Bottom1;            
    Bottom2 = tree_data->Bottom2;    
    Bottom3 = tree_data->Bottom3;

    if (tree_data->TreeType == TTYPE_FX2IR) /* fx */
    {
        ExSpot        = dev_data->FxSpot;                
        ExMidNode     = tree_data->FxMidNode;
        ExVol         = tree_data->FxVol;
        FwdEx         = tree_data->FwdFx;
    }
    else /* equity */
    {
        /* select which structures we're using: eq or fx */
        ExSpot        = dev_data->EqSpot;
        ExMidNode     = tree_data->EqMidNode;
        ExVol         = tree_data->EqVol;
        FwdEx         = tree_data->FwdEq;
    }
    
    /* in the lattice function, dev_data->FxSpot is assumed to be the 
       spot fx at t+1 (which is naturally the case when going backwards).
       Since we go forward here, we need to set FxSpot before calling lattice */
    if (t < T)
    {
        if (Hyb3_FillGrid_3d(tree_data,
                            ExSpot,       /* equity or Fx points */
                            dev_data->gDash,
                            dev_data->kDashTimesX,
                            dev_data->kVar,
                            FwdEx, 
                            ExMidNode, 
                            ExVol, 
                            t+1,
                            FALSE) == FAILURE)
           {
               goto RETURN;
           }
       
    }
    
    /*  Update tree */
    if (Hyb3_Lattice(dev_data,
        t,
        T,
        mktvol_data,
        tree_data) == FAILURE)
    {
        goto RETURN;
    }  
    

    if (t == 0) /* initialise the first state price */
    {
        StatePrL = StatePr + Hyb3_Node_Offset(3, 0, 0, 0, tree_data);
        StatePrL[0] = 1.0;
    }
    
    /* update State Prices */
    if (t != T)
    {
        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr1L = StatePr1 + Hyb3_Node_Offset (3, i, j, t+1, tree_data);
                
                for (k = Bottom3[t+1][i][j]; k <= Top3[t+1][i][j]; k++)
                {
                    StatePr1L[k] = 0.;
                }
            }  /* for j */
        }  /* for i */
        
        Shift1 = dev_data->Shift1 + Hyb3_Node_Offset(1,0,0,t,tree_data);
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {   
            i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;
            
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            
            Discount2DL = dev_data->Discount_2D[DCurveD] + offset;
            
            quu = dev_data->quu + offset; 
            q0u = dev_data->q0u + offset;
            qdu = dev_data->qdu + offset;
            qu0 = dev_data->qu0 + offset;
            q00 = dev_data->q00 + offset;
            qd0 = dev_data->qd0 + offset;
            qud = dev_data->qud + offset;
            q0d = dev_data->q0d + offset;
            qdd = dev_data->qdd + offset;
            
            Shift2 = dev_data->Shift2 + offset;                
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;
                
                StatePrL = StatePr + Hyb3_Node_Offset (3, i, j, t, tree_data);
                
                offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                
                r = dev_data->r + offset; 
                
                Quu = quu[j]; Qu0 = qu0[j]; Qud = qud[j];
                Q0u = q0u[j]; Q00 = q00[j]; Q0d = q0d[j];
                Qdu = qdu[j]; Qd0 = qd0[j]; Qdd = qdd[j];
                
                Shift3    = dev_data->Shift3 + offset;
                
                StatePr11L = StatePr1 + Hyb3_Node_Offset (3, i0, j0, t+1, tree_data);
                StatePr12L = StatePr1 + Hyb3_Node_Offset (3, i0, j1, t+1, tree_data);
                StatePr13L = StatePr1 + Hyb3_Node_Offset (3, i0, j2, t+1, tree_data);
                StatePr21L = StatePr1 + Hyb3_Node_Offset (3, i1, j0, t+1, tree_data);
                StatePr22L = StatePr1 + Hyb3_Node_Offset (3, i1, j1, t+1, tree_data);
                StatePr23L = StatePr1 + Hyb3_Node_Offset (3, i1, j2, t+1, tree_data);
                StatePr31L = StatePr1 + Hyb3_Node_Offset (3, i2, j0, t+1, tree_data);
                StatePr32L = StatePr1 + Hyb3_Node_Offset (3, i2, j1, t+1, tree_data);
                StatePr33L = StatePr1 + Hyb3_Node_Offset (3, i2, j2, t+1, tree_data);
                
                Discount_ij = Discount2DL[j];
                
                for (k = Bottom3[t][i][j]; k <= Top3[t][i][j]; k++)
                {
                    k0 = (k1 = (k2 = k + Shift3[k] - 1) + 1) + 1;
                    
                    x = StatePrL[k] * Discount_ij;
                    
                    Rux = r[k].u * x; 
                    R0x = r[k].m * x; 
                    Rdx = r[k].d * x;
                    
                    StatePr11L[k0] += Quu * Rux;
                    StatePr11L[k1] += Quu * R0x;
                    StatePr11L[k2] += Quu * Rdx;
                    StatePr12L[k0] += Qu0 * Rux;
                    StatePr12L[k1] += Qu0 * R0x;
                    StatePr12L[k2] += Qu0 * Rdx;
                    StatePr13L[k0] += Qud * Rux;
                    StatePr13L[k1] += Qud * R0x;
                    StatePr13L[k2] += Qud * Rdx;
                    StatePr21L[k0] += Q0u * Rux;
                    StatePr21L[k1] += Q0u * R0x;
                    StatePr21L[k2] += Q0u * Rdx;
                    StatePr22L[k0] += Q00 * Rux;
                    StatePr22L[k1] += Q00 * R0x;
                    StatePr22L[k2] += Q00 * Rdx;
                    StatePr23L[k0] += Q0d * Rux;
                    StatePr23L[k1] += Q0d * R0x;
                    StatePr23L[k2] += Q0d * Rdx;
                    StatePr31L[k0] += Qdu * Rux;
                    StatePr31L[k1] += Qdu * R0x;
                    StatePr31L[k2] += Qdu * Rdx;
                    StatePr32L[k0] += Qd0 * Rux;
                    StatePr32L[k1] += Qd0 * R0x;
                    StatePr32L[k2] += Qd0 * Rdx;
                    StatePr33L[k0] += Qdd * Rux;
                    StatePr33L[k1] += Qdd * R0x;
                    StatePr33L[k2] += Qdd * Rdx;
                    
                }/* for k */
            }/* for j */
        }/* for i */
                
    }/* if t!= T */   
    
    status = SUCCESS;
    
RETURN:
        
    return (status);
    
}  /* Hyb3_UpdateStatePrices3D */
    


/*****  Hyb3_BuildStatePrices2D  ***********************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*       Remark: this function only acts on EQ-data (i.e. TTYPE_EQ1IR )
*/
int     Hyb3_UpdateStatePrices2D
    (int              t,                      /* (I) current time period   */  
     MKTVOL_DATA     *mktvol_data,            /* (I) Market vol data       */
     HYB3_TREE_DATA  *tree_data,              /* (I) Tree data             */
     HYB3_DEV_DATA   *dev_data,               /* (I/O) dev-data for tree   */
     TSLICE           StatePr,                /* (I) initial state price   */
     TSLICE           StatePr1)               /* (O) final State price     */
{

    /* Slices and pointers to slices */

    double  *StatePrL    = NULL;      /* pointers state price slices      */
    double  *StatePr1L   = NULL;

    double  *StatePrM1L;
    double  *StatePrM2L;
    double  *StatePrM3L;

    double  *Discount1DL;           /* pointer to domestic discount     */
    
    TPROB_0 *p;
    TPROB_0 *q;

    int
       *Top1,     *Bottom1,         /* Limits of the tree (1rst dimension)  */
      **Top2,    **Bottom2;         /* Limits of the tree (2nd dimension)   */

    int     *Shift1, *Shift2;

    /* local variables for easy referencing */
    double    Pu ,   P0 ,  Pd ;
    double    Qux,   Q0x,  Qdx;

    int       i0, i1, i2,
              j0, j1, j2;
            
    double   Discount_i,               /* domestic discount            */
             x;

    int      i, j;


    int     DCurveD,           /* domestic discount curve */

            T,           /* Total number of period in the hybrids tree*/
            offset,
            status = FAILURE; /* Error status = FAILURE initially          */



    /* Total size of tree timeline */     
    T   = tree_data->NbTP;


    /* Assigment of domestic discount curve */
    DCurveD   = tree_data->CvDisc[0];

       
    Top1    = tree_data->Top1;	        
    Top2    = tree_data->Top2;
    Bottom1 = tree_data->Bottom1;            
    Bottom2 = tree_data->Bottom2;    


    /* in the lattice function, dev_data->EqSpot is assumed to be the 
    spot eq at t+1 (which is naturally the case when going backwards).
    Since we go forward here, we need to set EqSpot before calling lattice */
    if (t < T)
    {
        if (Hyb3_FillGrid_2d(tree_data,
                             dev_data->EqSpot,
                             dev_data->gDash,
                             dev_data->kDashTimesX,
                             dev_data->kVar,
                             tree_data->FwdEq,
                             tree_data->EqMidNode,
                             tree_data->EqVol,
                             t+1,
                             FALSE ) == FAILURE)
        {
            goto RETURN;
        }
    }
    
    /*  Update tree */
    if (Hyb3_Lattice(dev_data,
                    t,
                    T,
                    mktvol_data,
                    tree_data) == FAILURE)
    {
        goto RETURN;
    }  
    
    if (t == 0) /* initialise the first state price */
    {
        StatePrL = StatePr + Hyb3_Node_Offset(2, 0, 0, 0, tree_data);
        StatePrL[0] = 1.0;
    }
    
    /* update State Prices */
    if (t != T)
    {
        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Hyb3_Node_Offset (2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr1L[j] = 0.;
            }
        }
        
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        
        Shift1      = dev_data->Shift1               + offset;
        Discount1DL = dev_data->Discount_1D[DCurveD] + offset;
        p           = dev_data->p                    + offset;
        
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {   
            i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;
            
            Discount_i = Discount1DL[i];
            
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            
            q = dev_data->q + offset;
            
            Shift2   = dev_data->Shift2 + offset;                
            StatePrL = StatePr          + offset;
            
            StatePrM1L = StatePr1 + Hyb3_Node_Offset (2, i0, 0, t+1, tree_data);
            StatePrM2L = StatePr1 + Hyb3_Node_Offset (2, i1, 0, t+1, tree_data);
            StatePrM3L = StatePr1 + Hyb3_Node_Offset (2, i2, 0, t+1, tree_data);
            
            Pu = p[i].u;
            P0 = p[i].m;
            Pd = p[i].d;
            
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                j0 = (j1 = (j2 = j + Shift2[j] - 1) + 1) + 1;
                
                x = StatePrL[j] * Discount_i;
                
                Qux = q[j].u * x;
                Q0x = q[j].m * x;
                Qdx = q[j].d * x;
                
                StatePrM1L[j0] += Pu * Qux;
                StatePrM1L[j1] += Pu * Q0x;
                StatePrM1L[j2] += Pu * Qdx;
                
                StatePrM2L[j0] += P0 * Qux;
                StatePrM2L[j1] += P0 * Q0x;
                StatePrM2L[j2] += P0 * Qdx;
                
                StatePrM3L[j0] += Pd * Qux;
                StatePrM3L[j1] += Pd * Q0x;
                StatePrM3L[j2] += Pd * Qdx;
                
                
            }/* for j */
        }/* for i */
        
        
    }/* if t!= T */
                     

    status = SUCCESS;

    RETURN:


    return (status);

}  /* Hyb3_UpdateStatePrices3D */
