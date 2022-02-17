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
*       Remark: this function only acts on treetypes:
*        -TTYPE_2IR
*        -TTYPE_EQ1IR
*        -TTYPE_2IR2F1D
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
            
    double   x;

    int      i, j;


    int     DCurve,           /* discount curve */

            T,           /* Total number of period in the hybrids tree*/
            offset,
            status = FAILURE; /* Error status = FAILURE initially          */


    if ( ( tree_data->TreeType != TTYPE_EQ1IR) && ( tree_data->TreeType != TTYPE_2IR) && ( tree_data->TreeType != TTYPE_2IR2F1D))
    {
        DR_Error("UpdateStatePrice2D: only 2D equity and CUPS mode supported in this function.\n");
        goto RETURN;
    }

    /* Total size of tree timeline */     
    T   = tree_data->NbTP;


    /* Assigment of discount curve */
    DCurve  = ((tree_data->TreeType == TTYPE_2IR) ? tree_data->CvDisc[1] : tree_data->CvDisc[0]);

    
    Top1    = tree_data->Top1;	        
    Top2    = tree_data->Top2;
    Bottom1 = tree_data->Bottom1;            
    Bottom2 = tree_data->Bottom2;    


    /* in the lattice function, dev_data->EqSpot is assumed to be the 
    spot eq at t+1 (which is naturally the case when going backwards).
    Since we go forward here, we need to set EqSpot before calling lattice */
    if ((t < T) && (tree_data->TreeType == TTYPE_EQ1IR))
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
        Shift1 = ((tree_data->TreeType == TTYPE_2IR2F1D) ? dev_data->Shift4 : dev_data->Shift1) + offset;
        p      = ((tree_data->TreeType == TTYPE_2IR2F1D) ? dev_data->s      : dev_data->p)    + offset;

        /* we separate between  EQ1IR and CUPS mode as the domestic */
        /* discount factor is in different dimensions               */
        if (tree_data->TreeType == TTYPE_EQ1IR)
        {
            double Discount_i;
            double *Discount1DL;

            Discount1DL = dev_data->Discount_1D[DCurve] + offset;
            
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

        } /* if EQ_TREE */
        else 
        {
            double  *Discount2DL;           /* pointer to domestic discount (2IR), 2nd dim. */

            for (i = Bottom1[t]; i <= Top1[t]; i++)
            {   
                i0 = (i1 = (i2 = i + Shift1[i] - 1) + 1) + 1;

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                Discount2DL = dev_data->Discount_2D[DCurve] + offset;

                q        = ((tree_data->TreeType == TTYPE_2IR) ? dev_data->q      : dev_data->t)      + offset;
                Shift2   = ((tree_data->TreeType == TTYPE_2IR) ? dev_data->Shift2 : dev_data->Shift5) + offset;
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

                    x = StatePrL[j] * Discount2DL[j];

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
        }

    }/* if t!= T */
                     

    status = SUCCESS;

    RETURN:


    return (status);

}  /* Hyb3_UpdateStatePrices2D */



/*****  Hyb3_BuildStatePrices1D  ***********************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*       Remark: this function only acts on 
*        -TTYPE_2IR
*        -TTYPE_FX2IR
*        -TTYPE_1IR
*        -TTYPE_EQ1IR
*/
int     Hyb3_UpdateStatePrices1D
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
    double  *StatePrML   = NULL;
    
    double  *Discount1DL  = NULL;     /* pointer to foreign discount (2IR)*/

    TPROB_0 *s;

    int     *Top1,     *Bottom1;      /* Limits of the tree (1rst dimension)  */

    int     *Shift;

    /* local variables for easy referencing */
    double   Su ,   S0 ,  Sd ;
          
    int      i0, i1, i2;

    double   x;

    int      i;

    int      DCurve,           /* discount curve */

             T,                /* Total number of period in the hybrids tree*/
             offset,
             status = FAILURE; /* Error status = FAILURE initially          */


    if ( ( tree_data->TreeType != TTYPE_2IR) && 
         ( tree_data->TreeType != TTYPE_FX2IR) &&
         ( tree_data->TreeType != TTYPE_1IR) &&
         ( tree_data->TreeType != TTYPE_EQ1IR) )
    {
        DR_Error("UpdateStatePrice1D: only 2IR and 2IR+FX supported in this function.\n");  
        goto RETURN;
    }

    /* Total size of tree timeline */     
    T   = tree_data->NbTP;

    /* Assigment of discount curve */
    DCurve   = tree_data->CvDisc[0];
       
    Top1    = tree_data->Top1;	        
    Bottom1 = tree_data->Bottom1;            
  
    /*  Update tree */
    if (Hyb3_Lattice(dev_data,
                    t,
                    T,
                    mktvol_data,
                    tree_data) == FAILURE)
    {
        goto RETURN;
    }  
    
    offset = tree_data->HalfWidth[0];

    if (t == 0) /* initialise the first state price */
    {
        StatePrL = StatePr + offset;
        StatePrL[0] = 1.0;
    }
    
    /* update State Prices */
    if (t != T)
    {
        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + offset;  
            StatePr1L[i] = 0.;
        }

        Shift = dev_data->Shift4 + offset;
        s     = dev_data->s      + offset;
        
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {   
            i0 = (i1 = (i2 = i + Shift[i] - 1) + 1) + 1;

            Discount1DL = dev_data->Discount_1D[DCurve] + offset;

            StatePrL    = StatePr  + offset;
            StatePrML   = StatePr1 + offset;

            Su = s[i].u;
            S0 = s[i].m;
            Sd = s[i].d;

            x = StatePrL[i] * Discount1DL[i];

            StatePrML[i0] += Su * x;
            StatePrML[i1] += S0 * x;
            StatePrML[i2] += Sd * x;

        }/* for i */
        
    }/* if t!= T */
                     

    status = SUCCESS;

    RETURN:


    return (status);

}  /* Hyb3_UpdateStatePrices1D */


/*****  Hyb3_BuildStatePrices2D_1D  **************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*       updates the 2D+1D mode in the tree
*/
int     Hyb3_UpdateStatePrices2D_1D
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

    double  *Discount3DL;           /* pointer to domestic discount     */

    double  *quu,  *qu0,  *qud;     /*  2-D probability slices          */
    double  *q0u,  *q00,  *q0d;
    double  *qdu,  *qd0,  *qdd;

    TPROB_0   *r;


    int
        *Top1,     *Bottom1,         /* Limits of the tree (1rst dimension)  */
        **Top2,    **Bottom2,         /* Limits of the tree (2nd dimension)   */
        ***Top3,   ***Bottom3;         /* Limits of the tree (3rd dimension)   */

    int     *Shift1, *Shift2, *Shift3;

    /* local variables for easy referencing */
    double    Quu,   Qu0,  Qud;         /* 2-D probabilities (values)   */
    double    Q0u,   Q00,  Q0d;
    double    Qdu,   Qd0,  Qdd;

    double    Rux,   R0x,  Rdx;        /* 3-D probabilities (values)   */

    int       i0, i1, i2,
        j0, j1, j2,
        k0, k1, k2;

    double   x;

    int      i, j, k;

    int     DCurveD,           /* domestic discount curve */

        T,           /* Total number of period in the hybrids tree*/
        offset,
        status = FAILURE; /* Error status = FAILURE initially          */


    if (tree_data->TreeType != TTYPE_2IR2F1D )
    {
        DR_Error("Hyb3_UpdateStatePrices2D_1D: only implemented for Hyb2+1!");
        goto RETURN;
    }

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

                Discount3DL = dev_data->Discount_3D[DCurveD] + offset;

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


                for (k = Bottom3[t][i][j]; k <= Top3[t][i][j]; k++)
                {
                    k0 = (k1 = (k2 = k + Shift3[k] - 1) + 1) + 1;

                    x = StatePrL[k] * Discount3DL[k];

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

}  /* Hyb3_UpdateStatePrices2D_1D */



/*****  Hyb3_BuildStatePrices  **************************************/
/*
*       updates the state prices from t to t+1 going forward in the tree.
*       basically calls the correct updating mode
*/
int     Hyb3_UpdateStatePrices
        (int             t,                      /* (I) current time period   */  
        MKTVOL_DATA     *mktvol_data,            /* (I) Market vol data       */
        HYB3_TREE_DATA  *tree_data,              /* (I) Tree data             */
        HYB3_DEV_DATA   *dev_data,               /* (I/O) dev-data for tree   */
        TSLICE           StatePr,                /* (I) initial state price   */
        TSLICE           StatePr1)               /* (O) final State price     */
{
    int status = FAILURE;


    switch (tree_data->TreeType)
    {
    case TTYPE_2IR:
    case TTYPE_EQ1IR:
    case TTYPE_1IR2F:
        Hyb3_UpdateStatePrices2D(t, mktvol_data, tree_data, dev_data, StatePr, StatePr1 );
        break;
    case TTYPE_2IR2F1D:
        Hyb3_UpdateStatePrices2D_1D(t, mktvol_data, tree_data, dev_data, StatePr, StatePr1 );
        break;
    case TTYPE_FX2IR:
    case TTYPE_EQD2IR:
    case TTYPE_EQF2IR:
    case TTYPE_EQC2IR:
        Hyb3_UpdateStatePrices3D(t, mktvol_data, tree_data, dev_data, StatePr, StatePr1 );
        break;
    case TTYPE_EQDFX2IR:
    case TTYPE_EQFFX2IR:
        /* not implemented yet */
        DR_Error("Hyb3_UpdateStatePrice: currently not implemented for 4D mode!");
        goto RETURN;
        break;
    case TTYPE_1IR:
        Hyb3_UpdateStatePrices1D(t, mktvol_data, tree_data, dev_data, StatePr, StatePr1 );
        break;
    default:
        DR_Error("Hyb3_UpdateStatePrice: internal error: mode not supported!");
        goto RETURN;
    }

    status = SUCCESS;

RETURN:

    return (status);

}
