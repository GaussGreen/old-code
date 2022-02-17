/****************************************************************************/
/*      Position the nodes in the lattice.                                  */
/****************************************************************************/
/*      LATTICE.c                                                           */
/****************************************************************************/

/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"




/**                        MACROS FOR IR Q MODEL                             
   The following macros deal with the definition of the Q-mapping function   
   in the normal, Q=0, case.  The  INIT macro deals with initialisation at   
   top of a loop and the  UPDT  macro deals with the update of the grid at   
   the bottom of the loop.                                                 */

#ifndef INIT_GRID
#define INIT_GRID(Q,M,V,i) (((fabs(Q)) > (QCUTOFF)) ? ((M)*(exp((Q)*(V)*(i)))) : ((M)*(V)*(i)))
#endif

#ifndef UPDT_GRID
#define UPDT_GRID(Q,Grid,J) (((fabs(Q)) > (QCUTOFF)) ? ((Grid)*(J)) : ((Grid)+(J)))
#endif



/*******************************************************/
/*
 ******************************************************/
int Hyb3_SolveForShift(
                  double  *R,           /**<(O) = (rup - rdown)                             */
                  long     *nextn,      /**<(O) node shift in 3rd D                         */
                  long     prevn,       /**<(I) prev guess for node shift                   */
                  double  c1,           /**<(I) multiplier of rup in the first mom linear eq*/
                  double  c2,           /**<(I) multiplier of r0  in the first mom linear eq*/
                  double  c3,           /**<(I) multiplier of rd in the first mom linear eq */
                  double  JumpCoeff,    /**<(I)                                             */
                  double  TargetFwd,    /**<(I)  one period forward FX                      */
                  int     *NoSolution
                  )
{

    int status = FAILURE;
    double  rloc;
    
    double  a,b,c;
    double  discriminant;
    double  dummy;
    
    
    a     = (c1 + c3 - 2 * c2);
    b     = (c1 - c3);
    c     = 2 * (c2 - TargetFwd) + JumpCoeff * a;
    
    
    if (fabs(a) < TINY)
    {
        if(fabs(b) < TINY)
        {               
            *R = 0.0;
            *nextn = prevn;
            *NoSolution = TRUE;
        }
        else
        {   
            rloc    = - c/b;

            dummy   = prevn + rloc;
            if((dummy > DR_SHRT_MAX) || (dummy < DR_SHRT_MIN))
            {
                *R      = 0.0;
                *nextn  = prevn;
                *NoSolution = TRUE;
            }
            
            else
            {
                *R      = dummy;
                *nextn  = NEAR_INT(dummy);
                *NoSolution = FALSE;
            }

        }
    }
    else
    {
        discriminant =  b * b - 4 * a * c;
        
        if(discriminant < TINY) /* no solution */
        {
            *R     = 0.0;
            *nextn = prevn; 
            *NoSolution = TRUE;
        }
        else
        {
            double root1,root2;
            double  sqrdiscr;
            double  a_x_2;

            a_x_2 = 2 * a;

            sqrdiscr = sqrt(discriminant);
            root1 = (-b + sqrdiscr)/(a_x_2);
            root2 = (-b - sqrdiscr)/(a_x_2);

            rloc = (fabs(root1) < fabs(root2)) ? root1 : root2;
            dummy = prevn + rloc;

            if((dummy > DR_SHRT_MAX) || (dummy < DR_SHRT_MIN))
            {
                *R          = 0.0;
                *nextn      = prevn;
                *NoSolution = TRUE;
            }

            
            else
            {
                *R    = dummy;
                *nextn = NEAR_INT(dummy);
                *NoSolution = FALSE;
            }
        }
    }

    status = SUCCESS;
    return(status);

}


/*****  Hyb3_Lattice   ***********************************************************/
/**
*       Position the nodes in the lattice: calculate the one period discount
*       factor and the  probabilities at each node in the tree, according to
*       the tree type chosen by  the caller  (i.e. one of six types:  2ccy's
*       with CUPS, FX+2IR, EqDom+2IR, EqFor+2IR, EqComp+2IR and Eq+1IR).
*
*       This routine may work in two or three dimensions depending on chosen
*       tree type -  this implies that  some local variables may be declared
*       and not used.
*
*****************************************************************************/

int     Hyb3_Lattice(
                HYB3_DEV_DATA  *dev_data,     /**< (O) Shifts, probs, assets   */
                int             t,            /**< (I) Current time period     */
                int             T,            /**< (I) Total nb of per in tree */
                MKTVOL_DATA     *mktvol_data, /**< (I) Market vol data         */
                HYB3_TREE_DATA  *tree_data)   /**< (I) Structure of tree data  */
{
    int    Top1,      Bottom1;           /* Tree limits (IR foreign)       */
    int   *Top2,     *Bottom2;           /* Tree limits (IR domestic)      */
    int  **Top3,    **Bottom3;           /* Tree limits (third asset)      */
    int    OutTop1,   OutBottom1;        /* Tree limits (IR foreign)       */
    int   *OutTop2,  *OutBottom2;        /* Tree limits (IR domestic)      */
    int  **OutTop3, **OutBottom3;        /* Tree limits (third asset)      */


    int   *Shift1L;
    int   *Shift2L;
    int   *Shift3L;
    int   *Shift4L;
    int   *Shift5L;

    int   Mid[2];                     /* Mid point, dist goes left->right  */
          
    int    i, j, k;                   /* Node indices                      */
    int    l, lMin, lMax;             /* Node branching shifts             */
    int    m, mMin, mMax;             /* Node branching shifts             */
    int    mMinNoCups, mMaxNoCups;    /* Shift lim. in no cups 2IR2F1D mode*/
    int    n, nMin, nMax;             /* Node branching shifts             */


    int    offset;

    int    NbIRatesOn;               /* 1 for TTYPE_EQ1IR, 2 for all others */

    int    status = FAILURE;         /* Status = FAILURE initially          */


    /* Grid building variables (for efficiency) */
    double    GridIrF;                /* Grid pts foreign interest rate      */
    double    GridIrD;                /* Grid pts for domestic interest rate */

    /* Approximation for zero coupon diffusion, 2nd and 3rd curves */
    double    ZDRatio[2]={0,0};       /* Ratio of Zeros for domestic         */
    double    ZFRatio[2];             /* Ratio of Zeros for foreign          */
                                      /* [1]=index/COF [2]=index/risk        */
        
    double    Length;                 /* Current time step in years          */
    double    LengthJ;                /* Time step used for J-size calc      */
    double    Beta1, Beta2, Beta3;    /* Mean reversion terms                */
    double    Jump[6];                /* Jump sizes at current period        */
    double    PrevJump[6];            /* Jump sizes at prev period           */


    /* In the following variables, 0=foreign, 1=domestic                 */
    double    QLeft[2];               /* Q values for left  portion of dist  */
    double    QRight[2];              /* Q values for right portion of dist  */
    double    FwdShift[2];

    double    VolBbq[2];              /* Sigma used in bone mapping          */
    double    QMid[2];                /* Average q parameter                 */

    double    MLeft[2];               /* M and S are coefficients in the     */
    double    MRight[2];              /* mapping function                    */
    double    SLeft[2];
    double    SRight[2];
        
    double    FwdRateA[2];            /* Forward rate adjusted for Q shift   */
    double    FwdRate[2];             /* Forward rate                        */

    double    ZCenter[2];             /* IR grid center as calibrated        */
    double    RateJumpLeft[2];        /* Jump for IR grid, left portion      */
    double    RateJumpRight[2];       /* Jump for IR grid, right portion     */
    double    Zidx[2];                /* Auxiliary in IR grid construction   */
    double    JumpCoeff;                    
    double    d[6], du;

    register double   Pi;              /* -> First dimension                 */
    register double   Si;              /* -> First dimension (no CUPS)       */
    register double   Qi, Qij;         /* -> Second dimension                */
    register double   Ui, Uij;         /* -> Second dimension (no CUPS)      */
    register double   Ri, Rij, Rijk;   /* -> Third dimension                 */
    register double   pu,   p0,   pd;  /*       Local variables used         */
    register double   su,   s0,   sd;  /*            to store the            */
    register double   qu,   q0,   qd;  /*           probabilities            */
    register double   tu,   t0,   td;
    register double   ru,   r0,   rd;
   
    register double   Quu,Qu0,Qud;
    register double   Q0u,Q00,Q0d;
    register double   Qdu,Qd0,Qdd;

    

    double    *quuL,  *qu0L,  *qudL;   /*          1-D Pointers              */
    double    *q0uL,  *q00L,  *q0dL;   /*              used                  */
    double    *qduL,  *qd0L,  *qddL;   /*               for                  */
                                       /*            efficiency              */
    TPROB_0   *pL;                     
    TPROB_0   *sL;
    TPROB_0   *qL;
    TPROB_0   *tL;
    TPROB_0   *rL;



    double    *Discount_1D[3];           /*  Local variables using for smart   */
    double    *Discount_2D[3];           /*       addressing of slices         */
    double    *Discount_3D[3];
    
    int      CvDiffF = -999;
    int      CvIdx1F = -999;
    int      CvIdx2F = -999;
    int      CvDiscF = -999;
    int      CvDiffD = -999;
    int      CvIdx1D = -999;
    int      CvIdx2D = -999;
    int      CvDiscD = -999;


   
    /* Establish number of interest rates present in current mode */
    if ((tree_data->TreeType == TTYPE_EQ1IR)  ||
        (tree_data->TreeType == TTYPE_1IR)    ||
        (tree_data->TreeType == TTYPE_1IR2F))
    {
        NbIRatesOn = 1;
    }
    else
    {
        NbIRatesOn = 2;
    }

    /*  APPROXIMATION ON THE DIFFUSION OF THE COF/INDEX/RISK CURVES  */
    /*  We must store deterministic ratios of zero coupon bonds      */

    /* First dimension is always an interest rate */     
    CvDiffF = tree_data->CvDiff[0];    /* Internal assigments of */
    CvIdx1F = tree_data->CvIdx1[0];    /*  zero curves  */
    CvIdx2F = tree_data->CvIdx2[0];
    CvDiscF = tree_data->CvDisc[0];

    ZFRatio[0] = (1. + tree_data->FwdRate[0][CvDiffF][t]) 
               / (1. + tree_data->FwdRate[0][CvIdx1F][t]);
    ZFRatio[1] = (1. + tree_data->FwdRate[0][CvDiffF][t]) 
               / (1. + tree_data->FwdRate[0][CvIdx2F][t]);

    /* Unless the tree is running EQ + IR, the 2nd dim will be an IR */
    if (NbIRatesOn == 2)
    {

        CvDiffD = tree_data->CvDiff[1];  
        CvIdx1D = tree_data->CvIdx1[1];  
        CvIdx2D = tree_data->CvIdx2[1];  
        CvDiscD = tree_data->CvDisc[1];

        ZDRatio[0] = (1. + tree_data->FwdRate[1][CvDiffD][t]) 
                   / (1. + tree_data->FwdRate[1][CvIdx1D][t]);
        ZDRatio[1] = (1. + tree_data->FwdRate[1][CvDiffD][t]) 
                   / (1. + tree_data->FwdRate[1][CvIdx2D][t]);  
    }     
    

    /* Now process smile parameters for each IR present */
    for (i=0; i<NbIRatesOn; i++) 
    {
        QLeft[i]    = mktvol_data[i].QLeft;
        QRight[i]   = mktvol_data[i].QRight;
        FwdShift[i] = mktvol_data[i].FwdShift;

        QMid[i]   = (QLeft[i] + QRight[i])/2;

        VolBbq[i] = (1. + FwdShift[i]) / (1. + QMid[i] * FwdShift[i]);

        FwdRateA[i]  = tree_data->FwdRate[i][tree_data->CvDiff[i]][t] 
                        / (1. + FwdShift[i]);

        MLeft[i] = MRight[i] = FwdRateA[i];
        SLeft[i] = SRight[i] = 1 + FwdRateA[i];
        if (fabs(QLeft[i]) > QCUTOFF)
        {
            MLeft[i] /= QLeft[i];
            SLeft[i] -= FwdRateA[i] / QLeft[i];      
        }
        if (fabs(QRight[i]) > QCUTOFF)
        {
            MRight[i] /= QRight[i];
            SRight[i] -= FwdRateA[i] / QRight[i];      
        }
    }



    /*                       MAIN CODE                             */ 
    /* Now deal with one of four cases: one or two interest rates, */
    /* one interest rate and one equity or three factor where the  */
    /* 3rd asset may be equity or foreign exchange.                */



    /***  I - ONE INTEREST RATE   ***/
    if (tree_data->TreeType == TTYPE_1IR)
    {


        if (t == T)                      /* We  do not  need */
        {                                /* anything  at the */
            status = SUCCESS;            /* back of the tree.*/
            return (status);

        }

        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;   
            
    
        /* Current period coefficients [t] */
        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];
     
        Jump[0] = tree_data->Aweight[0][t] * du;   
        d[0] = (PrevJump[0] - Jump[0]) / Jump[0]; 
        Beta1  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];
          
        Top1    = tree_data->Top1[t];
        Bottom1 = tree_data->Bottom1[t];
   
        OutTop1    = tree_data->OutTop1[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
 
    

        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID             */
        /* First index of FwdRate is ccy (i.e. domestic or foreign) */
        /* and the second is zc ( i.e. COF, Index, Risk)            */
       

        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        Discount_1D[0] = dev_data->Discount_1D[0] + offset;
        Discount_1D[1] = dev_data->Discount_1D[1] + offset;
        Discount_1D[2] = dev_data->Discount_1D[2] + offset;


        /* Single IR: one level loop, innermost j'size is Jo */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[0] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[0]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[0]);
        
        
        /* Calculate mid point */
        Zidx[0]  = ZCenter[0]; 
        Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[0]) - 1;
        Mid[0]   = MINMAX (Mid[0], Bottom1 - 1, Top1);
        Zidx[0] += (PrevJump[0]) * Bottom1;
        

        /* SINGLE INTEREST RATE */

        /* Left for single IR */
        GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);

        for (i = Bottom1; i <= Mid[0]; i++)                   
        {                   
            /* Deal with foreign IR */
            Discount_1D[CvDiffF][i] = 1. / (SLeft[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] =ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] =ZFRatio[1] * Discount_1D[CvDiffF][i];   
      
            GridIrF = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);
        
        }  /* for i - Left portion, foreign IR */
    

        /* Right for foreign */
        Zidx[0] = ZCenter[0] + PrevJump[0] * (Mid[0] + 1);
        GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);


        for (i = Mid[0]+1; i <= Top1; i++)                   
        {                   
         
            /* Deal with foreign IR */
            Discount_1D[CvDiffF][i] = 1. / (SRight[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] =ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] =ZFRatio[1] * Discount_1D[CvDiffF][i];    
           
            GridIrF = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);
        
        }  /* for i */

        
            
        /*  (3) CALCULATION OF PROBABILITIES */


        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }


        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

        sL = dev_data->s + offset;
        Shift4L = dev_data->Shift4 + offset;

        for (i = Bottom1; i <= Top1; i++)    /* Single dimension  */
        {

            /* Drift of the foreign currency without CUPS */
            Si = (d[0] - Beta1) * i;

            l = NEAR_INT(Si);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift4L[i] = l;  
            Si -= l;

            sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
            sL[i].d = sd = su - Si;
            sL[i].m = s0 = 1. - su - sd;
 
        }  /* for i */
    }



    /***  II - TWO INTEREST RATES (CUPS)  ***/
    else if (tree_data->TreeType == TTYPE_2IR)
    {

        double  Drift1;   /* Auxiliary for determ. part of drift   */


        if (t == T)                      /* We  do not  need */
        {                                /* anything  at the */
            status = SUCCESS;            /* back of the tree.*/
            return (status);

        }

        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    

        /* No need to apply limits to jump sizes according to correlation*/
        /* level because the I/O manager will force correlations < 0.95  */
    
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;   
        PrevJump[1] = tree_data->Aweight[1][t-1] * du;   
        PrevJump[2] = tree_data->Aweight[2][t-1] * du;       
    
        /* Current period coefficients [t] */
        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        FwdRate[1] = tree_data->FwdRate[1][CvDiffD][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];
        ZCenter[1] = tree_data->IrZCenter[1][t];
    
        Jump[0] = tree_data->Aweight[0][t] * du;   
        Jump[1] = tree_data->Aweight[1][t] * du;   
        Jump[2] = tree_data->Aweight[2][t] * du;                

        d[0] = (PrevJump[0] - Jump[0]) / Jump[0];
        d[1] = (PrevJump[1] - Jump[1]) / Jump[2];
        d[2] = (PrevJump[2] - Jump[2]) / Jump[2];
    
        Beta1  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];
        Beta2  = mktvol_data[1].Beta[0] * Length * PrevJump[1] / Jump[2];
        Beta3  = mktvol_data[1].Beta[0] * Length * PrevJump[2] / Jump[2];
    
        Drift1 = tree_data->DriftCUPS[0][t] / Jump[0]; /* Foreign, CUPS    */
        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
    
        OutTop1    = tree_data->OutTop1[t+1];
        OutTop2    = tree_data->OutTop2[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
        OutBottom2 = tree_data->OutBottom2[t+1];
    

        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID             */
        /* First index of FwdRate is ccy (i.e. domestic or foreign) */
        /* and the second is zc ( i.e. COF, Index, Risk)            */
       

        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        Discount_1D[0] = dev_data->Discount_1D[0] + offset;
        Discount_1D[1] = dev_data->Discount_1D[1] + offset;
        Discount_1D[2] = dev_data->Discount_1D[2] + offset;


        /* Foreign IR: one level loop, innermost j'size is Jo */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[0] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[0]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[0]);
        
        /* Domestic IR: two level loop, innermost j'size is jump J2 */
        RateJumpLeft[1] = RateJumpRight[1] = PrevJump[2] * VolBbq[1] * FwdRateA[1];
        if (fabs(QLeft[1]) >QCUTOFF) RateJumpLeft[1] =exp(QLeft[1]*VolBbq[1]*PrevJump[2]);
        if (fabs(QRight[1])>QCUTOFF) RateJumpRight[1]=exp(QRight[1]*VolBbq[1]*PrevJump[2]);

        
        /* Calculate mid point */
        Zidx[0]  = ZCenter[0]; 
        Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[0]) - 1;
        Mid[0]   = MINMAX (Mid[0], Bottom1 - 1, Top1);
        Zidx[0] += (PrevJump[0]) * Bottom1;
        

        /* FOREIGN CURRENCY */

        /* Left for foreign */
        GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);

        for (i = Bottom1; i <= Mid[0]; i++)                   
        {                   
            /* Deal with foreign IR */
            Discount_1D[CvDiffF][i] = 1. / (SLeft[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] =ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] =ZFRatio[1] * Discount_1D[CvDiffF][i];   
      
            GridIrF = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);
        
        }  /* for i - Left portion, foreign IR */
    

        /* Right for foreign */
        Zidx[0] = ZCenter[0] + PrevJump[0] * (Mid[0] + 1);
        GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);


        for (i = Mid[0]+1; i <= Top1; i++)                   
        {                   
         
            /* Deal with foreign IR */
            Discount_1D[CvDiffF][i] = 1. / (SRight[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] =ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] =ZFRatio[1] * Discount_1D[CvDiffF][i];    
           
            GridIrF = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);
        
        }  /* for i */



        /* DOMESTIC CURRENCY */

        for (i=Bottom1; i<=Top1; i++)
        {

            /* Prepare to address domestic IR 2-D slices */                
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);   
                                                           
            Discount_2D[0] = dev_data->Discount_2D[0] + offset;
            Discount_2D[1] = dev_data->Discount_2D[1] + offset;
            Discount_2D[2] = dev_data->Discount_2D[2] + offset;


            /* Calculate Mid point */
            Zidx[1]  = ZCenter[1] + PrevJump[1] * i;
            Mid[1]   = (int) ceil(-Zidx[1] / PrevJump[2]) - 1;
            Mid[1]   = MINMAX (Mid[1], Bottom2[i] - 1, Top2[i]);
            Zidx[1] += (PrevJump[2]) * Bottom2[i];


            /* Left for domestic */
            GridIrD = INIT_GRID(QLeft[1], MLeft[1], VolBbq[1], Zidx[1]);

            for (j = Bottom2[i]; j <= Mid[1]; j++)            
            {  
                Discount_2D[CvDiffD][j] = 1. / (SLeft[1] + GridIrD);
                Discount_2D[CvIdx1D][j] = ZDRatio[0] * Discount_2D[CvDiffD][j];
                Discount_2D[CvIdx2D][j] = ZDRatio[1] * Discount_2D[CvDiffD][j];	
    
                GridIrD = UPDT_GRID(QLeft[1], GridIrD, RateJumpLeft[1]);         
                    
            }  /* for j */


            /* Right for domestic */
            Zidx[1] = ZCenter[1] + PrevJump[1] * i 
                    + PrevJump[2] * (Mid[1] + 1);
            GridIrD = INIT_GRID(QRight[1], MRight[1], VolBbq[1], Zidx[1]);

            for (j = Mid[1]+1; j <= Top2[i]; j++)            
            {
                Discount_2D[CvDiffD][j] = 1. / (SRight[1] + GridIrD);
                Discount_2D[CvIdx1D][j] = ZDRatio[0] * Discount_2D[CvDiffD][j];
                Discount_2D[CvIdx2D][j] = ZDRatio[1] * Discount_2D[CvDiffD][j];	
    
                GridIrD = UPDT_GRID(QRight[1], GridIrD, RateJumpRight[1]);
                    
            }  /* For j */
        
        } /* For i */

            


        /*  (3) CALCULATION OF PROBABILITIES */


        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }


        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

        sL = dev_data->s + offset;
        pL = dev_data->p + offset;

        Shift1L = dev_data->Shift1 + offset;
        Shift4L = dev_data->Shift4 + offset;


        for (i = Bottom1; i <= Top1; i++)    /* First dim: foreign rate */
        {

            /* Drift of first dim including mean reversion */
            Pi = Drift1 + (d[0] - Beta1) * i; 
          
            /* Part of 2nd dim drift which depend on 1st dim */  
            Qi = (d[1]-Beta2)*i - Pi*Jump[1]/ Jump[2];

            /* And drift of the foreign currency without CUPS */
            Si = (d[0] - Beta1) * i;


            l = NEAR_INT(Si);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift4L[i] = l;  
            Si -= l;

            sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
            sL[i].d = sd = su - Si;
            sL[i].m = s0 = 1. - su - sd;

            l = NEAR_INT(Pi);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift1L[i] = l;  
            Pi -= l;


            pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
            pL[i].d = pd = pu - Pi;
            pL[i].m = p0 = 1. - pu - pd;

            mMax =            OutTop2[i+l-1];
            mMax = MIN (mMax, OutTop2[i+l  ]);
            mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

            mMin =            OutBottom2[i+l-1];
            mMin = MAX (mMin, OutBottom2[i+l  ]);
            mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Hyb3_Lattice: problem in building the tree (mMin > mMax)!");
                goto RETURN;
            }

            /* Initialise 1-D pointers to beginning of 2-D blocks */
            offset = Hyb3_Node_Offset (2, i, 0, t, tree_data);

            qL = dev_data->q + offset;

            Shift2L = dev_data->Shift2 + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++) 
            {
                /* Second dim: domestic rate */
                Qij = Qi + (d[2] - Beta3) * j;
                
                m = NEAR_INT(Qij);
                m = MINMAX (mMin - j, m, mMax - j);
                Shift2L[j] = m; 
                Qij -= m;

                qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                qL[j].d = qd = qu - Qij;
                qL[j].m = q0 = 1. - qu - qd;

            }  /* for j */

        }  /* for i */

    }


    /***  III - ONE INTEREST RATE PLUS EQUITY  ***/

    else if (tree_data->TreeType == TTYPE_EQ1IR) 
    {
        double  *ExMidNode     = tree_data->EqMidNode;
        double  *TmpEqSlice    = NULL;
        double  *EqL           = NULL;

        double  *gDashL        = NULL;
        double  *kDashTimesxL  = NULL;
        double  *kVarL         = NULL;

        double   EqVar, Drift4, Drift3, Drift6, LnRatiosi, NormalisedEq;
        int      PrevIdx, SameIdx, SameParams;

        int      NoChangeInK = FALSE;
        int      ChangeInK   = FALSE;

        double   DeflatedFwdRatio = 0;
        double   CurrFwd = 0;

        int      count;
        int      j0,j1,j2;
        double   c1,c2,c3;
        
        long     nextm = 0;
        long     currm,temp;
        double   QLoc = 0;

        double   QGuess = 0;
        int      mGuess = 0;
        int      NoSolution = FALSE;

        int      iterate1D = 0;
        int      iterate2D = 0;

        /* variable for efficient access to Inner tree dims*/
        int           InnerTop1,      InnerBottom1;
        int          *InnerTop2,     *InnerBottom2;

        int         InnerTop2i,InnerBottom2i;

        int      Idx = tree_data->SmileIndex[t-1];

        double   a1  = tree_data->A1[Idx];
        double   a2  = tree_data->A2[Idx];
        double   a3  = tree_data->A3[Idx];
        
        int      IsLognormal = FALSE;
 
        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    
        /* No need to apply limits to jump sizes according to correlation */
        /* level because the I/O manager will force correlations < 0.95   */        
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;   
        PrevJump[1] = tree_data->Aweight[1][t-1] * du;   
        PrevJump[2] = tree_data->Aweight[2][t-1] * du;   
    
        /* Current period coefficients [t] */
        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];
    
        Jump[0] = tree_data->Aweight[0][t] * du;   
        Jump[1] = tree_data->Aweight[1][t] * du;   
        Jump[2] = tree_data->Aweight[2][t] * du;                

        d[0] = (PrevJump[0] - Jump[0]) / Jump[0];
        d[1] = (PrevJump[1] - Jump[1]) / Jump[2];
        d[2] = (PrevJump[2] - Jump[2]) / Jump[2];
    
        Beta1  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];

        /* Vars that are used not Node dependent and used whether K */
        /* has changed or not.                                      */

        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
		 
        OutTop1    = tree_data->OutTop1[t+1];
        OutTop2    = tree_data->OutTop2[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
        OutBottom2 = tree_data->OutBottom2[t+1];

        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID */

        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1, 0, 0, t, tree_data);

        Discount_1D[0] = dev_data->Discount_1D[0] + offset;
        Discount_1D[1] = dev_data->Discount_1D[1] + offset;
        Discount_1D[2] = dev_data->Discount_1D[2] + offset;

        /* Single IR: one factor, one loop */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[0] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[0]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[0]);

        Zidx[0]  = ZCenter[0]; 
        Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[0]) - 1;
        Mid[0]   = MINMAX (Mid[0], Bottom1 - 1, Top1);
        Zidx[0] += (PrevJump[0]) * Bottom1;


        /* INTEREST RATE DIMENSION */
        GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);
        

        /* Left for single interest rate */
        for (i = Bottom1; i <= Mid[0]; i++)                   
        {	          
            Discount_1D[CvDiffF][i] = 1. / (SLeft[0] + GridIrF);
            Discount_1D[CvIdx1F][i] = ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] = ZFRatio[1] * Discount_1D[CvDiffF][i];        
    
            GridIrF = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);
             
        }  /* for i */


        /* Right for single interest rate */
        Zidx[0] = ZCenter[0] + PrevJump[0] * (Mid[0] + 1);
        GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);

        for (i = Mid[0]+1; i <= Top1; i++)                   
        {	               
            Discount_1D[CvDiffF][i] = 1. / (SRight[0] + GridIrF);
            Discount_1D[CvIdx1F][i] = ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] = ZFRatio[1] * Discount_1D[CvDiffF][i];	
    
            GridIrF  = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);
    
        }  /* for i */

        /* SECOND DIMENSION: EQUITY */

        /* Swap Next and current Eq slices */
        TmpEqSlice           = dev_data->NextEqSpot;
        dev_data->NextEqSpot = dev_data->EqSpot;
        dev_data->EqSpot     = TmpEqSlice;

        if (Hyb3_FillGrid_2d(tree_data,
                             dev_data->EqSpot,
                             dev_data->gDash, 
                             dev_data->kDashTimesX, 
                             dev_data->kVar,   
                             tree_data->FwdEq,
                             tree_data->EqMidNode,
                             tree_data->EqVol,
                             t,
                             TRUE) == FAILURE)
        {
            goto RETURN;
        }

        if (Hyb3_FillGridEqFwd_2d(tree_data,dev_data,t) == FAILURE)
        {
            goto RETURN;
        }

        if (t == T)                      /* We  don't  need  the */
        {                                /* probabilities at the */
            status = SUCCESS;            /* back of the tree.    */
            goto RETURN;
        }

        /*  (3) CALCULATION OF PROBABILITIES */
        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;
        if (lMin > lMax)
        {
            DR_Error ("Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }

        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        Discount_1D[CvDiscF] = dev_data->Discount_1D[CvDiscF] + offset;

        pL = dev_data->p + offset;
        sL = dev_data->s + offset;

        Shift1L = dev_data->Shift1 + offset;
        Shift4L = dev_data->Shift4 + offset;

        InnerTop1    = tree_data->InnerTreeTop1[t];
        InnerTop2    = tree_data->InnerTreeTop2[t];
        InnerBottom1 = tree_data->InnerTreeBottom1[t];
        InnerBottom2 = tree_data->InnerTreeBottom2[t];
        
        Drift4 = (ExMidNode[t] - ExMidNode[t+1])/Jump[2];
        
        EqVar  = 0.5 * tree_data->SpotEqVol[t] * tree_data->SpotEqVol[t] * Length / Jump[2];

        DeflatedFwdRatio = tree_data->ZeroCoupon[0][CvDiscF][t+1] /
                           tree_data->ZeroCoupon[0][CvDiscF][t];

        Drift3 = log(DeflatedFwdRatio) / Jump[2];

        DeflatedFwdRatio *= tree_data->FwdEq[t+1] /
                            tree_data->FwdEq[t];

        Idx     = tree_data->SmileIndex[t];
        PrevIdx = tree_data->SmileIndex[t-1];
            
        a1 = tree_data->A1[Idx];
        a2 = tree_data->A2[Idx];
        a3 = tree_data->A3[Idx];
            
        IsLognormal = ((a2 < A2CUTOFF) && (fabs(a1) < TINY));

        SameIdx     = (Idx == PrevIdx);
        SameParams  = ((fabs(a1 - tree_data->A1[PrevIdx]) < TINY) &&
                       (fabs(a2 - tree_data->A2[PrevIdx]) < TINY) &&
                       (fabs(a3 - tree_data->A3[PrevIdx]) < TINY));

        NoChangeInK = (SameIdx || SameParams);
        ChangeInK   = !NoChangeInK;

        /* to improve running time, treat log-normal EQ seperately */

        if(IsLognormal && NoChangeInK)
        {
            for (i = Bottom1; i <= Top1; i++) /* First dimension: foreign rate */
            {
                /* Drift of first dim including mean reversion */
                Pi = (d[0] - Beta1) * i; 
       
                LnRatiosi = Drift3 - log(Discount_1D[CvDiscF][i]) / Jump[2];
      
                /* Part of 2nd dim drift which depend on 1st dim */  
                Qi = Drift4 + d[1]*i - Pi*Jump[1]/Jump[2];

                l = NEAR_INT(Pi);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift1L[i] = l;  
                Shift4L[i] = l;
                Pi -= l;

                pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
                pL[i].d = pd = pu - Pi;
                pL[i].m = p0 = 1. - pu - pd;

                sL[i].u = pu;
                sL[i].d = pd;
                sL[i].m = p0;

                mMax =            OutTop2[i+l-1];
                mMax = MIN (mMax, OutTop2[i+l  ]);
                mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

                mMin =            OutBottom2[i+l-1];
                mMin = MAX (mMin, OutBottom2[i+l  ]);
                mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

                if (mMin > mMax)
                {
                    DR_Error ("Lattice: problem in building the tree (mMin > mMax)!");
                    goto RETURN;
                }

                /* Initialise 1-D pointers to beginning of 2-D blocks */

                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                qL = dev_data->q + offset;
                
                Shift2L  = dev_data->Shift2 + offset;
                
                for (j = Bottom2[i]; j <= Top2[i]; j++) /* Second dim: domestic rate */
                {
                    Drift6 = LnRatiosi - EqVar;

                    Qij = Qi + Drift6 + d[2] * j;
            
                    m = NEAR_INT(Qij);
                    m = MINMAX (mMin - j, m, mMax - j); 

                    Shift2L[j] = m; 
                    Qij -= m;

                    qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                    qL[j].d = qd = qu - Qij;
                    qL[j].m = q0 = 1. - qu - qd;

                }  /* for j */
            }  /* for i */
        }/* Lognormal equity */
        else
        {   
            double JRatio1;
            double MixJRatio1,MixJRatio2;
            double Kchange = - ExMidNode[t]/Jump[2];
            double Kchangei = 0;
            double Kchangeij = 0;

            double ZeroRatioi;

            double *NextEq0 = NULL;
            double *NextEq1 = NULL;
            double *NextEq2 = NULL;

            JRatio1    = - Jump[1]/Jump[2];

            MixJRatio1 = - PrevJump[1]/Jump[2];
            MixJRatio2 = - PrevJump[2]/Jump[2];

            for (i = Bottom1; i <= Top1; i++)    /* First dimension: foreign rate */
            {
                if (ChangeInK)
                {
                    Kchangei = Kchange + i*MixJRatio1;
                }

                InnerTop2i = InnerBottom2i = 0;

                iterate1D = 0;

                if((i <= InnerTop1) && (i>= InnerBottom1))
                {
                    iterate1D     = 1;
                    InnerTop2i    = InnerTop2[i];
                    InnerBottom2i = InnerBottom2[i];
                }

               /* Drift of first dim including mean reversion */
                Pi = (d[0] - Beta1) * i; 
       
                /* Part of 2nd dim drift which depend on 1st dim */  
                Qi = Drift4 + d[1]*i + Pi*JRatio1;

                l = NEAR_INT(Pi);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift1L[i] = l;  
                Shift4L[i] = l;
                Pi -= l;

                pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
                pL[i].d = pd = pu - Pi;
                pL[i].m = p0 = 1. - pu - pd;

                sL[i].u = pu;
                sL[i].d = pd;
                sL[i].m = p0;
       
                mMax =            OutTop2[i+l-1];
                mMax = MIN (mMax, OutTop2[i+l  ]);
                mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

                mMin =            OutBottom2[i+l-1];
                mMin = MAX (mMin, OutBottom2[i+l  ]);
                mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

                if (mMin > mMax)
                {
                    DR_Error ("Lattice: problem in building the tree (mMin > mMax)!");
                    goto RETURN;
                }

                /* Initialise 1-D pointers to beginning of 2-D blocks */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
    
                EqL = dev_data->EqSpot + offset;

                qL = dev_data->q + offset;

                Shift2L  = dev_data->Shift2 + offset;

                /* local variable pointing to the drift functions */
                gDashL       = dev_data->gDash       + offset;
                kDashTimesxL = dev_data->kDashTimesX + offset;
                kVarL        = dev_data->kVar        + offset;

                NextEq0 = dev_data->NextEqSpot + Hyb3_Node_Offset(2,i+l+1,0,t+1,tree_data);
                NextEq1 = dev_data->NextEqSpot + Hyb3_Node_Offset(2,i+l,  0,t+1,tree_data);
                NextEq2 = dev_data->NextEqSpot + Hyb3_Node_Offset(2,i+l-1,0,t+1,tree_data);

                LnRatiosi  = Drift3 - log(Discount_1D[CvDiscF][i]) / Jump[2];
                ZeroRatioi = 1 / Discount_1D[CvDiscF][i];
      
                for (j = Bottom2[i]; j <= Top2[i]; j++)
                {
                    iterate2D = 0;

                    if(iterate1D)
                    {
                        if((j <= InnerTop2i) && (j >= InnerBottom2i))
                        {
                            iterate2D = 1;
                        }
                    }

                    NormalisedEq = EqL[j] / tree_data->FwdEq[t];

                    Drift6 = LnRatiosi*kDashTimesxL[j] - gDashL[j] * EqVar;

                    if (ChangeInK)
                    {
                        Kchangeij = Kchangei + j*MixJRatio2 + kVarL[j]/Jump[2];
                        Drift6   += Kchangeij;
                    }

                    CurrFwd = EqL[j] * DeflatedFwdRatio * ZeroRatioi;
                    Qij     = Qi + Drift6 + d[2] * j;
            
                    m = NEAR_INT(Qij);

                    /********************************************/
                    /* initialise variables used to iteratively */
                    /* improve the fit to FwdEq                 */
                    /********************************************/
                    
                    mGuess = currm = nextm = m;
                    QGuess = QLoc  = Qij;

                    /*************************************************/
                    /* Try to improve QGuess using only              */
                    /* exisiting grid points.                        */
                    /* if at any iteration:                          */
                    /*    1- no solution found or                    */
                    /*    2- the solution takes us outside limits    */
                    /* THEN revert to QGuess                         */
                    /*************************************************/

                    for(count = 1 ; count <= (int) (MAX_FX_ITER * iterate2D) ; count++)
                    {

                        j1 = j  + currm;
                        j0 = j1 + 1;
                        j2 = j1 - 1;
                    
                        if((j1 > mMax) || (j1 < mMin))
                        {  
                            nextm = currm - 1L;
                            break;
                        }

                        /* Coeffs for first moment equation */

                        c1 = pu * NextEq0[j0] + p0 * NextEq1[j0] + pd * NextEq2[j0];
                        c2 = pu * NextEq0[j1] + p0 * NextEq1[j1] + pd * NextEq2[j1];
                        c3 = pu * NextEq0[j2] + p0 * NextEq1[j2] + pd * NextEq2[j2];

                        if(Hyb3_SolveForShift(&QLoc,&nextm,currm,c1,c2,c3,JumpCoeff,CurrFwd,&NoSolution) == FAILURE)
                        {
                            DR_Error("pb solving shift\n");
                            goto RETURN;
                        }
                        
                    
                        if(NoSolution)
                        {
                            nextm = currm - 1L;
                            break;
                        }
                        else
                        {
                    
                            if (nextm == currm)
                            {
                                break;
                            }
                            else
                            {   
                                temp  = currm;
                                currm = nextm;
                                nextm = temp;
                            }
                        }
                    } /* for count */

                    /******************************************************/
                    /* IF we exit the loop                                */
                    /* 1- without having found a solution OR              */
                    /* 2- the solution takes us outside the initial nodes */
                    /* THEN revert to QGuess                              */
                    /******************************************************/

                    if(nextm != currm)
                    {
                        QLoc  = QGuess;
                        nextm = mGuess;
                    }
                    else
                    {
                        if(nextm != mGuess)
                        {
                            nextm = mGuess;
                            QLoc  = QGuess;
                        }  
                    }


                    nextm      = MINMAX (mMin - j, nextm, mMax - j); 
                    Qij        = QLoc - nextm;
                    Shift2L[j] = nextm; 

                    qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                    qL[j].d = qd = qu - Qij;
                    qL[j].m = q0 = 1. - qu - qd;

                }  /* for j */
            }  /* for i */
        }
    }

    /***  IV - ONE INTEREST RATES TWO FACTOR (CUPS)  ***/
    else if (tree_data->TreeType == TTYPE_1IR2F)
    {

       
        if (t == T)                      /* We  do not  need */
        {                                /* anything  at the */
            status = SUCCESS;            /* back of the tree.*/
            return (status);

        }

        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    

        /* No need to apply limits to jump sizes according to correlation*/
        /* level because the I/O manager will force correlations < 0.95  */
    
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;   
        PrevJump[1] = tree_data->Aweight[1][t-1] * du;   
        PrevJump[2] = tree_data->Aweight[2][t-1] * du;       
    
        /* Current period coefficients [t] */
        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];

    
        Jump[0] = tree_data->Aweight[0][t] * du;   
        Jump[1] = tree_data->Aweight[1][t] * du;   
        Jump[2] = tree_data->Aweight[2][t] * du;                

        d[0] = (PrevJump[0] - Jump[0]) / Jump[0];
        d[1] = (PrevJump[1] - Jump[1]) / Jump[2];
        d[2] = (PrevJump[2] - Jump[2]) / Jump[2];
    
        Beta1  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];
        Beta2  = mktvol_data[0].Beta[1] * Length * PrevJump[1] / Jump[2];
        Beta3  = mktvol_data[0].Beta[1] * Length * PrevJump[2] / Jump[2];
    
        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
    
        OutTop1    = tree_data->OutTop1[t+1];
        OutTop2    = tree_data->OutTop2[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
        OutBottom2 = tree_data->OutBottom2[t+1];
    

        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID             */
        /* First index of FwdRate is ccy (i.e. domestic or foreign) */
        /* and the second is zc ( i.e. COF, Index, Risk)            */
       

        /* Grid jumps in both part of distribution */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[2] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[2]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[2]);
        
        for (i=Bottom1; i<=Top1; i++)
        {

            /* Prepare to address domestic IR 2-D slices */                
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);   
                                                           
            Discount_2D[0] = dev_data->Discount_2D[0] + offset;
            Discount_2D[1] = dev_data->Discount_2D[1] + offset;
            Discount_2D[2] = dev_data->Discount_2D[2] + offset;


            /* Calculate Mid point */
            Zidx[0]  = ZCenter[0] +(PrevJump[0]+ PrevJump[1]) * i;
            Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[2]) - 1;
            Mid[0]   = MINMAX (Mid[0], Bottom2[i] - 1, Top2[i]);
            Zidx[0] += (PrevJump[2]) * Bottom2[i];


            /* Left part of distribution */
            GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);

            for (j = Bottom2[i]; j <= Mid[0]; j++)            
            {  
                Discount_2D[CvDiffF][j] = 1. / (SLeft[0] + GridIrF);
                Discount_2D[CvIdx1F][j] = ZFRatio[0] * Discount_2D[CvDiffF][j];
                Discount_2D[CvIdx2F][j] = ZFRatio[1] * Discount_2D[CvDiffF][j];	
    
                GridIrF = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);         
                    
            }  /* for j */


            /* Right part of distribution */
            Zidx[0] = ZCenter[0] + (PrevJump[0]+PrevJump[1]) * i 
                    + PrevJump[2] * (Mid[0] + 1);
            GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);

            for (j = Mid[0]+1; j <= Top2[i]; j++)            
            {
                Discount_2D[CvDiffF][j] = 1. / (SRight[0] + GridIrF);
                Discount_2D[CvIdx1F][j] = ZFRatio[0] * Discount_2D[CvDiffF][j];
                Discount_2D[CvIdx2F][j] = ZFRatio[1] * Discount_2D[CvDiffF][j];	
    
                GridIrF = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);
                    
            }  /* For j */
        
        } /* For i */

            


        /*  (3) CALCULATION OF PROBABILITIES */


        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }


        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

        sL = dev_data->s + offset;

        Shift4L = dev_data->Shift4 + offset;


        for (i = Bottom1; i <= Top1; i++)    /* First dim: first factor of IR */
        {

            /* Drift of first dim including mean reversion */
            Si =  (d[0] - Beta1) * i; 
          
            /* Part of 2nd dim drift which depend on 1st dim */  
            Ui = (d[1]-Beta2)*i - Si*Jump[1]/ Jump[2];


            l = NEAR_INT(Si);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift4L[i] = l;  
            Si -= l;


            sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
            sL[i].d = sd = su - Si;
            sL[i].m = s0 = 1. - su - sd;

            mMax =            OutTop2[i+l-1];
            mMax = MIN (mMax, OutTop2[i+l  ]);
            mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

            mMin =            OutBottom2[i+l-1];
            mMin = MAX (mMin, OutBottom2[i+l  ]);
            mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Hyb3_Lattice: problem in building the tree (mMin > mMax)!");
                goto RETURN;
            }

            /* Initialise 1-D pointers to beginning of 2-D blocks */
            offset = Hyb3_Node_Offset (2, i, 0, t, tree_data);

            tL = dev_data->t + offset;

            Shift5L = dev_data->Shift5 + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++) 
            {
                /* Second dim: second factor of IR */
                Uij = Ui + (d[2] - Beta3) * j;
                
                m = NEAR_INT(Uij);
                m = MINMAX (mMin - j, m, mMax - j);
                Shift5L[j] = m; 
                Uij -= m;

                tL[j].u = tu = .5 * (JumpCoeff + Uij + Uij*Uij);
                tL[j].d = td = tu - Uij;
                tL[j].m = t0 = 1. - tu - td;

            }  /* for j */

        }  /* for i */

    }
    
    /***  V - TWO IR:TWO FACTOR FOREIGN IR, ONE FACTOR DOM IR ***/
    else if (tree_data->TreeType == TTYPE_2IR2F1D)
    {

        double  Drift1,
                Drift2;        /* Auxiliary for determ. part of drift   */
        double  Beta11,        /* 3D Beta coefficient */
                Beta21,
                Beta31,
                Beta22,
                Beta32,
                Beta33;


        if (t == T)                      /* We  do not  need */
        {                                /* anything  at the */
            status = SUCCESS;            /* back of the tree.*/
            return (status);

        }

        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    

        /* No need to apply limits to jump sizes according to correlation*/
        /* level because the I/O manager will force correlations < 0.95  */
    
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;   
        PrevJump[1] = tree_data->Aweight[1][t-1] * du;   
        PrevJump[2] = tree_data->Aweight[2][t-1] * du;   
        PrevJump[3] = tree_data->Aweight[3][t-1] * du;   
        PrevJump[4] = tree_data->Aweight[4][t-1] * du;   
        PrevJump[5] = tree_data->Aweight[5][t-1] * du;   

    
        /* Current period coefficients [t] */
        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        FwdRate[1] = tree_data->FwdRate[1][CvDiffD][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];
        ZCenter[1] = tree_data->IrZCenter[1][t];
    
        Jump[0] = tree_data->Aweight[0][t] * du;   
        Jump[1] = tree_data->Aweight[1][t] * du;   
        Jump[2] = tree_data->Aweight[2][t] * du;   
        Jump[3] = tree_data->Aweight[3][t] * du;   
        Jump[4] = tree_data->Aweight[4][t] * du;   
        Jump[5] = tree_data->Aweight[5][t] * du;   


        d[0] = (PrevJump[0] - Jump[0]) / Jump[0];
        d[1] = (PrevJump[1] - Jump[1]) / Jump[2];
        d[2] = (PrevJump[2] - Jump[2]) / Jump[2];
        d[3] = (PrevJump[3] - Jump[3]) / Jump[5];
        d[4] = (PrevJump[4] - Jump[4]) / Jump[5];
        d[5] = (PrevJump[5] - Jump[5]) / Jump[5];


        Beta11  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];
        Beta21  = mktvol_data[0].Beta[1] * Length * PrevJump[1] / Jump[2];
        Beta22  = mktvol_data[0].Beta[1] * Length * PrevJump[2] / Jump[2];
        Beta31  = mktvol_data[1].Beta[0] * Length * PrevJump[3] / Jump[5];
        Beta32  = mktvol_data[1].Beta[0] * Length * PrevJump[4] / Jump[5];
        Beta33  = mktvol_data[1].Beta[0] * Length * PrevJump[5] / Jump[5];

    
        Drift1 = tree_data->DriftCUPS[0][t] / Jump[0];      /* Foreign first factor CUPS drift*/
        Drift2 = tree_data->DriftCUPS[1][t] / Jump[2];    /* Foreign second factor CUPS drift*/
        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
    
        OutTop1    = tree_data->OutTop1[t+1];
        OutTop2    = tree_data->OutTop2[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
        OutBottom2 = tree_data->OutBottom2[t+1];
    
        Top3    = tree_data->Top3[t];
        Bottom3 = tree_data->Bottom3[t];
        OutTop3    = tree_data->OutTop3[t+1];
        OutBottom3 = tree_data->OutBottom3[t+1];
        
        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID             */
        /* First index of FwdRate is ccy (i.e. domestic or foreign) */
        /* and the second is zc ( i.e. COF, Index, Risk)            */
       




        /* Foreign IR: two level loop, innermost j'size is J2 */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[2] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[2]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[2]);
        
        /* Domestic IR: three level loop, innermost j'size is jump J5 */
        RateJumpLeft[1] = RateJumpRight[1] = PrevJump[5] * VolBbq[1] * FwdRateA[1];
        if (fabs(QLeft[1]) >QCUTOFF) RateJumpLeft[1] =exp(QLeft[1]*VolBbq[1]*PrevJump[5]);
        if (fabs(QRight[1])>QCUTOFF) RateJumpRight[1]=exp(QRight[1]*VolBbq[1]*PrevJump[5]);
        
        for (i =Bottom1; i<=Top1; i++)
        {

            /* Prepare to address slices */
            offset = Hyb3_Node_Offset(2, i, 0, t, tree_data);

            Discount_2D[0] = dev_data->Discount_2D[0] + offset;
            Discount_2D[1] = dev_data->Discount_2D[1] + offset;
            Discount_2D[2] = dev_data->Discount_2D[2] + offset;

            /* Calculate mid point */
            Zidx[0]  = ZCenter[0] + (PrevJump[0] + PrevJump[1]) * i; 
            Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[2]) - 1;
            Mid[0]   = MINMAX (Mid[0], Bottom2[i] - 1, Top2[i]);
            Zidx[0] += (PrevJump[2]) * Bottom2[i];
            /* FOREIGN CURRENCY */

            /* Left for foreign */
            GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);

            for (j = Bottom2[i]; j <= Mid[0]; j++)                   
            {                   
                /* Deal with foreign IR */
                Discount_2D[CvDiffF][j] = 1. / (SLeft[0] + GridIrF); /* Foreign */
                Discount_2D[CvIdx1F][j] =ZFRatio[0] * Discount_2D[CvDiffF][j];
                Discount_2D[CvIdx2F][j] =ZFRatio[1] * Discount_2D[CvDiffF][j];   
      
                GridIrF = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);
        
            }  /* for j - Left portion, foreign IR */
    

            /* Right for foreign */
            Zidx[0] = ZCenter[0] + (PrevJump[0]+PrevJump[1]) * i 
                    + PrevJump[2] * (Mid[0] + 1);
            GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);


            for (j = Mid[0]+1; j <= Top2[i]; j++)                   
            {                   
         
                /* Deal with foreign IR */
                Discount_2D[CvDiffF][j] = 1. / (SRight[0] + GridIrF); /* Foreign */
                Discount_2D[CvIdx1F][j] =ZFRatio[0] * Discount_2D[CvDiffF][j];
                Discount_2D[CvIdx2F][j] =ZFRatio[1] * Discount_2D[CvDiffF][j];    
           
                GridIrF = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);
        
            }  /* for j */
        } /* for i */



        /* DOMESTIC CURRENCY */

        for (i=Bottom1; i<=Top1; i++)
        {
            for (j=Bottom2[i]; j <= Top2[i];j++)
            {

                /* Prepare to address domestic IR 2-D slices */                
                offset = Hyb3_Node_Offset(3, i, j, t, tree_data);   
                                                           
                Discount_3D[0] = dev_data->Discount_3D[0] + offset;
                Discount_3D[1] = dev_data->Discount_3D[1] + offset;
                Discount_3D[2] = dev_data->Discount_3D[2] + offset;


                /* Calculate Mid point */
                Zidx[1]  = ZCenter[1] + PrevJump[3] * i + PrevJump[4] * j;
                Mid[1]   = (int) ceil(-Zidx[1] / PrevJump[5]) - 1;
                Mid[1]   = MINMAX (Mid[1], Bottom3[i][j] - 1, Top3[i][j]);
                Zidx[1] += (PrevJump[5]) * Bottom3[i][j];


                /* Left for domestic */
                GridIrD = INIT_GRID(QLeft[1], MLeft[1], VolBbq[1], Zidx[1]);

                for (k = Bottom3[i][j]; k <= Mid[1]; k++)            
                {  
                    Discount_3D[CvDiffD][k] = 1. / (SLeft[1] + GridIrD);
                    Discount_3D[CvIdx1D][k] = ZDRatio[0] * Discount_3D[CvDiffD][k];
                    Discount_3D[CvIdx2D][k] = ZDRatio[1] * Discount_3D[CvDiffD][k];
    
                    GridIrD = UPDT_GRID(QLeft[1], GridIrD, RateJumpLeft[1]);         
                    
                }  /* for k */


                /* Right for domestic */
                Zidx[1] = ZCenter[1] + PrevJump[3] * i + PrevJump[4] * j 
                        + PrevJump[5] * (Mid[1] + 1);
                GridIrD = INIT_GRID(QRight[1], MRight[1], VolBbq[1], Zidx[1]);

                for (k = Mid[1]+1; k <= Top3[i][j]; k++)            
                {
                    Discount_3D[CvDiffD][k] = 1. / (SRight[1] + GridIrD);
                    Discount_3D[CvIdx1D][k] = ZDRatio[0] * Discount_3D[CvDiffD][k];
                    Discount_3D[CvIdx2D][k] = ZDRatio[1] * Discount_3D[CvDiffD][k];	
    
                    GridIrD = UPDT_GRID(QRight[1], GridIrD, RateJumpRight[1]);
                    
                }  /* For k */
        
            } /* For j */
        } /* For i */

            


        /*  (3) CALCULATION OF PROBABILITIES */


        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }


        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);

        sL = dev_data->s + offset;
        pL = dev_data->p + offset;

        Shift1L = dev_data->Shift1 + offset;
        Shift4L = dev_data->Shift4 + offset;


        for (i = Bottom1; i <= Top1; i++)    /* First dim: foreign rate */
        {

            /* Drift of first dim including mean reversion */
            Pi = Drift1 + (d[0] - Beta11) * i; 
          
            /* Part of 2nd and 3rd dim drift which depend on 1st dim */  
            Qi = Drift2 + (d[1] - Beta21) * i - Pi * Jump[1] / Jump[2];
            Ri = (d[3] - Beta31) * i - Pi * Jump[3] / Jump[5];

            /* And drift of two factors of the foreign currency without CUPS */
            Si = (d[0] - Beta11) * i;
            Ui = (d[1] - Beta21) * i - Si * Jump[1] / Jump[2];

            l = NEAR_INT(Si);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift4L[i] = l;  
            Si -= l;

            sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
            sL[i].d = sd = su - Si;
            sL[i].m = s0 = 1. - su - sd;

            mMaxNoCups =            OutTop2[i+l-1];
            mMaxNoCups = MIN (mMaxNoCups, OutTop2[i+l  ]);
            mMaxNoCups = MIN (mMaxNoCups, OutTop2[i+l+1]) - 1;

            mMinNoCups =            OutBottom2[i+l-1];
            mMinNoCups = MAX (mMinNoCups, OutBottom2[i+l  ]);
            mMinNoCups = MAX (mMinNoCups, OutBottom2[i+l+1]) + 1;
            if (mMinNoCups > mMaxNoCups)
            {
                DR_Error ("Hyb3_Lattice: problem in building the tree (mMinNoCups > mMaxNoCups)!");
                goto RETURN;
            }



            l = NEAR_INT(Pi);          /* Node shift */
            l = MINMAX (lMin - i , l, lMax - i);
            Shift1L[i] = l;  
            Pi -= l;


            pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
            pL[i].d = pd = pu - Pi;
            pL[i].m = p0 = 1. - pu - pd;

            mMax =            OutTop2[i+l-1];
            mMax = MIN (mMax, OutTop2[i+l  ]);
            mMax = MIN (mMax, OutTop2[i+l+1]) - 1;

            mMin =            OutBottom2[i+l-1];
            mMin = MAX (mMin, OutBottom2[i+l  ]);
            mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;



            if (mMin > mMax)
            {
                DR_Error ("Hyb3_Lattice: problem in building the tree (mMin > mMax)!");
                goto RETURN;
            }

            /* Initialise 1-D pointers to beginning of 2-D blocks */
            offset = Hyb3_Node_Offset (2, i, 0, t, tree_data);

            qL = dev_data->q + offset;
            tL = dev_data->t + offset;
            quuL = dev_data->quu + offset; 
            qu0L = dev_data->qu0 + offset; 
            qudL = dev_data->qud + offset;
            q0uL = dev_data->q0u + offset; 
            q00L = dev_data->q00 + offset; 
            q0dL = dev_data->q0d + offset;
            qduL = dev_data->qdu + offset; 
            qd0L = dev_data->qd0 + offset; 
            qddL = dev_data->qdd + offset;

            Shift2L = dev_data->Shift2 + offset;
            Shift5L = dev_data->Shift5 + offset;

            for (j = Bottom2[i]; j <= Top2[i]; j++) 
            {
                /* Second dim: 2F of foreign rate */
                Qij = Qi + (d[2] - Beta22) * j;
                /* Part of 3rd dim drift which depend on 2nd dim */
                Rij = Ri + (d[4] - Beta32) * j -Qij *Jump[4] / Jump[5];
                /* And drift of the second factor of foreign currency without CUPS */
                Uij = Ui + (d[2] - Beta22) * j;

                /*Prob of 2nd factor Foreign IR without CUPS */

                m = NEAR_INT(Uij);
                m = MINMAX (mMinNoCups - j, m, mMaxNoCups - j);
                Shift5L[j] = m; 
                Uij -= m;
                
                tL[j].u = tu = .5 * (JumpCoeff + Uij + Uij*Uij);
                tL[j].d = td = tu - Uij;
                tL[j].m = t0 = 1. - tu - td;

                /*Prob of 2nd factor Foreign IR with CUPS */
                m = NEAR_INT(Qij);
                m = MINMAX (mMin - j, m, mMax - j);
                Shift2L[j] = m; 
                Qij -= m;

                qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                qL[j].d = qd = qu - Qij;
                qL[j].m = q0 = 1. - qu - qd;
                quuL[j] = pu * qu; qu0L[j] = pu * q0; qudL[j] = pu * qd;
                q0uL[j] = p0 * qu; q00L[j] = p0 * q0; q0dL[j] = p0 * qd;
                qduL[j] = pd * qu; qd0L[j] = pd * q0; qddL[j] = pd * qd;
                    




                nMax =            OutTop3[i+l-1][j+m-1];
                nMax = MIN (nMax, OutTop3[i+l-1][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l-1][j+m+1]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m-1]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l  ][j+m+1]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m-1]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m  ]);
                nMax = MIN (nMax, OutTop3[i+l+1][j+m+1]) - 1;
                    
                nMin =            OutBottom3[i+l-1][j+m-1];
                nMin = MAX (nMin, OutBottom3[i+l-1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l-1][j+m+1]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m-1]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l  ][j+m+1]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m-1]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m  ]);
                nMin = MAX (nMin, OutBottom3[i+l+1][j+m+1]) + 1;
                    
                    
                    if (nMin > nMax)
                    {
                        DR_Error ("Hyb3_Lattice: problem in building the tree (nMin > nMax)!");
                        goto RETURN;
                    }
                    
                    
                    /* Initialise 1-D pointers to beginning of 3-D blocks */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    rL = dev_data->r + offset;
                    
                    Shift3L = dev_data->Shift3   + offset;
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Rijk         = Rij + (d[5] - Beta33) * k  ;

                        n            = NEAR_INT(Rijk);
                        n            = MINMAX (nMin - k, n, nMax - k);
                        Shift3L[k]   = n;
                        Rijk        -= n;

                        rL[k].u = ru = .5 * (JumpCoeff + Rijk + Rijk * Rijk);
                        rL[k].d = rd = ru - Rijk;
                        rL[k].m = r0 = 1. - ru - rd;

                    } /* for k */

     

            }  /* for j */

        }  /* for i */
    }

    /***  VI - THREE FACTOR TREE WITH EQUITY OR FX  ***/
    else 
        
    {
        /* "Ex" stands for "Equity or Fx" */

        int         EquityOff = (tree_data->TreeType == TTYPE_FX2IR);
        int         EquityOn  = !EquityOff;
        int         EquityDom = EquityOn && !(tree_data->TreeType == TTYPE_EQF2IR);
        int         EquityFor = EquityOn && !EquityDom;
        int         EquityOffOrEquityDom = (EquityOff || EquityDom);
        int         ChecksOn  = tree_data->CalcCheckSlices;

        int         Idx;
        int         PrevIdx;

        double      *TmpSlice      = NULL;
        double      *ExSpot        = NULL;
        double      *ExMidNode     = NULL;
        double      *SpotExVol     = NULL;
        double      *ExVol         = NULL;

        double      Drift1;   /* Auxiliary to store deterministic drift */
        double      Drift3;   /* to assets 1 (CUPS) and asset 3 (EQ/FX) */
        double      Drift4;
        double      Drift6;

        double      JRatio1,JRatio3,JRatio4;

        double      ExpJump3;
        double      ExpJump4;
        double      ExpJump5;

        double      a1;
        double      a2;
        double      a3;

        double      *NextExSpot = NULL;
        double      *FwdEx      = NULL;

        double      *NextEx00 = NULL;
        double      *NextEx01 = NULL;
        double      *NextEx02 = NULL;
        double      *NextEx10 = NULL;
        double      *NextEx11 = NULL;
        double      *NextEx12 = NULL;
        double      *NextEx20 = NULL;
        double      *NextEx21 = NULL;
        double      *NextEx22 = NULL;

        double      *ExL      = NULL;
        
        double      CurrFwd;
        double      ExVar;
        double      CupsDriftRaw;
        double      CupsDrift;
        double      DeflatedFwdRatio;
                          
        int         k0,k1,k2;
        
        double      ZeroRatioi,ZeroRatioij;
        double      LnRatiosi,LnRatiosij;
        double      *gDashL;
        double      *kDashTimesxL;
        double      *kVarL;
        
        int         SameParams, SameIdx;
        int         NoChangeInK = FALSE;
        int         ChangeInK   = FALSE;
        int         IsLognormal = FALSE;

        int         count;
        double      c1,c2,c3;
        
        long        nextn = 0;
        long        currn,temp;
        double      RLoc = 0;

        double      RGuess = 0;
        int         nGuess = 0;
        int         NoSolution = FALSE;

        int         iterate1D = 0;
        int         iterate2D = 0;
        int         iterate3D = 0;
        
        /* variable for efficient access to Inner tree dims*/
        int           InnerTop1,      InnerBottom1;
        int          *InnerTop2,     *InnerBottom2;
        int         **InnerTop3,    **InnerBottom3;

        int         InnerTop2i,InnerBottom2i;
        int         InnerTop3ij,InnerBottom3ij;

        /* (1) COEFFICIENTS */

        /* Previous period coefficients [t-1] */
        LengthJ = tree_data->LengthJ[t-1];
        du      = sqrt (JUMPCOEFF * LengthJ);
    
        /* No need to apply limits to jump sizes according to correlation */
        /* level because the I/O manager will force correlations < 0.95   */
    
        PrevJump[0] = tree_data->Aweight[0][t-1] * du;
        PrevJump[1] = tree_data->Aweight[1][t-1] * du;
        PrevJump[2] = tree_data->Aweight[2][t-1] * du;
        PrevJump[3] = tree_data->Aweight[3][t-1] * du;
        PrevJump[4] = tree_data->Aweight[4][t-1] * du;
        PrevJump[5] = tree_data->Aweight[5][t-1] * du;

        ExpJump3 = exp(PrevJump[3]);
        ExpJump4 = exp(PrevJump[4]);
        ExpJump5 = exp(PrevJump[5]);
    
    
        /* Current period coefficients [t] */

        LengthJ = tree_data->LengthJ[t];
        du      = sqrt (JUMPCOEFF * LengthJ);
        Length  = tree_data->Length[t];
        JumpCoeff = Length/(LengthJ * JUMPCOEFF);
        FwdRate[0] = tree_data->FwdRate[0][CvDiffF][t];
        FwdRate[1] = tree_data->FwdRate[1][CvDiffD][t];
        ZCenter[0] = tree_data->IrZCenter[0][t];
        ZCenter[1] = tree_data->IrZCenter[1][t];


        Jump[0] = tree_data->Aweight[0][t] * du;   
        Jump[1] = tree_data->Aweight[1][t] * du;   
        Jump[2] = tree_data->Aweight[2][t] * du;   
        Jump[3] = tree_data->Aweight[3][t] * du;   
        Jump[4] = tree_data->Aweight[4][t] * du;   
        Jump[5] = tree_data->Aweight[5][t] * du; 

    
        JRatio1 = - Jump[1]/Jump[2];
        JRatio3 = - Jump[3]/Jump[5];
        JRatio4 = - Jump[4]/Jump[5];

        /* Drifts and coefficients used to calculate probas. */
        /* (The IR drifts (1, 2 and 4) have been calibrated  */
        /* beforehand. The 3rd asset is calculated here and  */
        /* comprises a term  for tree centre shift plus the  */
        /* term due to the lognormality of the asset.)       */

        d[0] = (PrevJump[0] - Jump[0]) / Jump[0];
        d[1] = (PrevJump[1] - Jump[1]) / Jump[2];
        d[2] = (PrevJump[2] - Jump[2]) / Jump[2];
        d[3] = (PrevJump[3] - Jump[3]) / Jump[5];
        d[4] = (PrevJump[4] - Jump[4]) / Jump[5];
        d[5] = (PrevJump[5] - Jump[5]) / Jump[5];
    
        Beta1  = mktvol_data[0].Beta[0] * Length * PrevJump[0] / Jump[0];
        Beta2  = mktvol_data[1].Beta[0] * Length * PrevJump[1] / Jump[2];
        Beta3  = mktvol_data[1].Beta[0] * Length * PrevJump[2] / Jump[2];
    
        Drift1 = tree_data->DriftCUPS[0][t] / Jump[0]; /* Foreign, CUPS    */
        
        Top1    = tree_data->Top1[t];
        Top2    = tree_data->Top2[t];
        Top3    = tree_data->Top3[t];
        Bottom1 = tree_data->Bottom1[t];
        Bottom2 = tree_data->Bottom2[t];
        Bottom3 = tree_data->Bottom3[t];
    
        OutTop1    = tree_data->OutTop1[t+1];
        OutTop2    = tree_data->OutTop2[t+1];
        OutTop3    = tree_data->OutTop3[t+1];
        OutBottom1 = tree_data->OutBottom1[t+1];
        OutBottom2 = tree_data->OutBottom2[t+1];
        OutBottom3 = tree_data->OutBottom3[t+1];

        /* (2) CALCULATION OF DISCOUNT FACTORS AND GRID             */
        /* First index of FwdRate is ccy (i.e. domestic or foreign) */
        /* and the second is zc (i.e. COF, Index, Risk)             */
    

        /* Prepare to address slices */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        Discount_1D[0] = dev_data->Discount_1D[0] + offset;
        Discount_1D[1] = dev_data->Discount_1D[1] + offset;
        Discount_1D[2] = dev_data->Discount_1D[2] + offset;



        /* Foreign IR: one level loop, innermost j'size is Jo */
        RateJumpLeft[0] = RateJumpRight[0] = PrevJump[0] * VolBbq[0] * FwdRateA[0];
        if (fabs(QLeft[0]) >QCUTOFF) RateJumpLeft[0]=exp(QLeft[0]*VolBbq[0]*PrevJump[0]);
        if (fabs(QRight[0])>QCUTOFF) RateJumpRight[0]=exp(QRight[0]*VolBbq[0]*PrevJump[0]);
        
        /* Domestic IR: two level loop, innermost j'size is jump J2 */
        RateJumpLeft[1] = RateJumpRight[1] = PrevJump[2] * VolBbq[1] * FwdRateA[1];
        if (fabs(QLeft[1]) >QCUTOFF) RateJumpLeft[1] =exp(QLeft[1]*VolBbq[1]*PrevJump[2]);
        if (fabs(QRight[1])>QCUTOFF) RateJumpRight[1]=exp(QRight[1]*VolBbq[1]*PrevJump[2]);

        /* Calculate Mid point (distribution switch) */
        Zidx[0]  = ZCenter[0]; 
        Mid[0]   = (int) ceil(-Zidx[0] / PrevJump[0]) - 1;
        Mid[0]   = MINMAX (Mid[0], Bottom1 - 1, Top1);
        Zidx[0] += (PrevJump[0]) * Bottom1;



        /* ASSET 1:  FOREIGN INTEREST RATE */

        /* Left for foreign */
        GridIrF = INIT_GRID(QLeft[0], MLeft[0], VolBbq[0], Zidx[0]);
        

        for (i = Bottom1; i <= Mid[0]; i++)
        {    
            
            /* Deal with foreign IR first */
            Discount_1D[CvDiffF][i] = 1. / (SLeft[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] = ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] = ZFRatio[1] * Discount_1D[CvDiffF][i];

            GridIrF  = UPDT_GRID(QLeft[0], GridIrF, RateJumpLeft[0]);
        
        }  /* for i */


        /* Right for Foreign */
        Zidx[0] = ZCenter[0] + PrevJump[0] * (Mid[0] + 1);
        GridIrF = INIT_GRID(QRight[0], MRight[0], VolBbq[0], Zidx[0]);
        
        for (i = Mid[0]+1; i <= Top1; i++)   
        {
            
            /* Deal with foreign IR first */
            Discount_1D[CvDiffF][i] = 1. / (SRight[0] + GridIrF); /* Foreign */
            Discount_1D[CvIdx1F][i] = ZFRatio[0] * Discount_1D[CvDiffF][i];
            Discount_1D[CvIdx2F][i] = ZFRatio[1] * Discount_1D[CvDiffF][i];

            GridIrF  = UPDT_GRID(QRight[0], GridIrF, RateJumpRight[0]);

        }  /* for i */


        /* ASSET 2:  DOMESTIC INTEREST RATE */

        for (i = Bottom1; i <= Top1; i++)   
        {
            /* Prepare to address slices */
            offset = Hyb3_Node_Offset(2,i,0,t,tree_data);
            Discount_2D[0] = dev_data->Discount_2D[0] + offset;
            Discount_2D[1] = dev_data->Discount_2D[1] + offset;
            Discount_2D[2] = dev_data->Discount_2D[2] + offset;

            /* Calculate Mid point */
            Zidx[1]  = ZCenter[1] + PrevJump[1] * i;
            Mid[1]   = (int) ceil(-Zidx[1] / PrevJump[2]) - 1;
            Mid[1]   = MINMAX (Mid[1], Bottom2[i] - 1, Top2[i]);
            Zidx[1] += (PrevJump[2]) * Bottom2[i];

            /* Left for domestic */
            GridIrD = INIT_GRID(QLeft[1], MLeft[1], VolBbq[1], Zidx[1]); 

            for (j = Bottom2[i]; j <= Mid[1]; j++)
            {
                /* Domestic grid */
                Discount_2D[CvDiffD][j] = 1. / (SLeft[1] + GridIrD); /* Domestic */
                Discount_2D[CvIdx1D][j] = ZDRatio[0]*Discount_2D[CvDiffD][j];
                Discount_2D[CvIdx2D][j] = ZDRatio[1]*Discount_2D[CvDiffD][j];
                
                GridIrD  = UPDT_GRID(QLeft[1], GridIrD, RateJumpLeft[1]);
                    
            }  /* for j */

            /* Right for domestic */
            Zidx[1] = ZCenter[1] + PrevJump[1] * i 
                    + PrevJump[2] * (Mid[1] + 1);
            GridIrD = INIT_GRID(QRight[1], MRight[1], VolBbq[1], Zidx[1]); 

            for (j = Mid[1]+1; j <= Top2[i]; j++) 
            {

                /* Domestic grid */
                Discount_2D[CvDiffD][j] = 1. / (SRight[1] + GridIrD); /* Domestic */
                Discount_2D[CvIdx1D][j] = ZDRatio[0]*Discount_2D[CvDiffD][j];
                Discount_2D[CvIdx2D][j] = ZDRatio[1]*Discount_2D[CvDiffD][j];
                  
                GridIrD  = UPDT_GRID(QRight[1],GridIrD,RateJumpRight[1]);
                  
            }  /* for j */

        } /* For i*/

        /* ASSET 3: EQUITY or FX ("EX") */

        /* switch between equity and fx */

        if (EquityOn)
        {
            /* Swap Next and current Eq slices */

            TmpSlice             = dev_data->NextEqSpot;
            dev_data->NextEqSpot = dev_data->EqSpot;
            dev_data->EqSpot     = TmpSlice;

            /* select which structures we're using: eq or fx */

            ExSpot        = dev_data->EqSpot;
            NextExSpot    = dev_data->NextEqSpot;

            ExMidNode     = tree_data->EqMidNode;
            SpotExVol     = tree_data->SpotEqVol;
            ExVol         = tree_data->EqVol;
            FwdEx         = tree_data->FwdEq;

            /* compute some deterministic values */

            if (EquityDom)
            {
                DeflatedFwdRatio = tree_data->ZeroCoupon[1][CvDiscD][t+1] /
                                   tree_data->ZeroCoupon[1][CvDiscD][t];
                CupsDriftRaw = 0;
            }
            else /* Equity is foreign  */
            {
                DeflatedFwdRatio = tree_data->ZeroCoupon[0][CvDiscF][t+1] /
                                   tree_data->ZeroCoupon[0][CvDiscF][t];

                CupsDriftRaw = - tree_data->Rho[5][t] *
                                 tree_data->SpotEqVol[t] *
                                 tree_data->SpotFxVol[t] *
                                 Length;
            }

            Drift3 = log(DeflatedFwdRatio) / Jump[5];

            DeflatedFwdRatio *= tree_data->FwdEq[t+1] /
                                tree_data->FwdEq[t];
        }
        else /* fx */
        {
            /* Swap Next and current Fx slices */

            TmpSlice             = dev_data->NextFxSpot;
            dev_data->NextFxSpot = dev_data->FxSpot;
            dev_data->FxSpot     = TmpSlice;

            /* select which structures we're using: eq or fx */

            ExSpot        = dev_data->FxSpot;                
            NextExSpot    = dev_data->NextFxSpot;

            ExMidNode     = tree_data->FxMidNode;
            SpotExVol     = tree_data->SpotFxVol;
            ExVol         = tree_data->FxVol;
            FwdEx         = tree_data->FwdFx;

            /* compute some deterministic values */

            DeflatedFwdRatio = 1;
            Drift3           = log(tree_data->FwdFx[t] /
                                   tree_data->FwdFx[t+1]) / Jump[5];
            CupsDriftRaw = 0;
        }

        CupsDrift = CupsDriftRaw / Jump[5];


        if (Hyb3_FillGrid_3d(tree_data,
                             ExSpot,
                             dev_data->gDash,
                             dev_data->kDashTimesX,
                             dev_data->kVar,
                             FwdEx,
                             ExMidNode,
                             ExVol,
                             t,
                             TRUE) == FAILURE)
        {
            goto RETURN;
        }


        if (EquityOn)
        {
            if (Hyb3_FillGridEqFwd_3d(tree_data,dev_data,t) == FAILURE)
            {
                goto RETURN;
            }
        }


        if (t == T)
        {   
            if(ChecksOn)
            {
                if(Hyb3_CopySlice(dev_data->FwdFX,
                             dev_data->FxSpot,
                             3,
                             T,
                            tree_data) == FAILURE) goto RETURN;
            
                      
                if(Hyb3_SetSlice(dev_data->DomZero,
                            2,
                            1.0,
                            T,
                            tree_data) == FAILURE) goto RETURN;
            }

            status = SUCCESS;            
            return (status);
        }

        /*  (3) CALCULATION OF PROBABILITIES */

        /* Store limits of outer ellipse */
        lMax = OutTop1    - 1;
        lMin = OutBottom1 + 1;
        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Lattice: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }

        /* Prepare to address 1-D slice */
        offset = Hyb3_Node_Offset(1,0,0,t,tree_data);
        Discount_1D[CvDiscF] = dev_data->Discount_1D[CvDiscF] + offset;

        pL = dev_data->p + offset;
        sL = dev_data->s + offset;

        Shift1L = dev_data->Shift1 + offset;
        Shift4L = dev_data->Shift4 + offset;
            
        InnerTop1    = tree_data->InnerTreeTop1[t];
        InnerTop2    = tree_data->InnerTreeTop2[t];
        InnerTop3    = tree_data->InnerTreeTop3[t];
        InnerBottom1 = tree_data->InnerTreeBottom1[t];
        InnerBottom2 = tree_data->InnerTreeBottom2[t];
        InnerBottom3 = tree_data->InnerTreeBottom3[t];
        
        Drift4 = (ExMidNode[t] - ExMidNode[t+1]) / Jump[5];

        ExVar  = 0.5 * SpotExVol[t] * SpotExVol[t] * Length / Jump[5];

        /* Determine whether smile function has changed */
        /* between t-1 and t                            */
        /* Either within same input smile period or  */
        /* two consecutive periods define same smile */

        Idx     = tree_data->SmileIndex[t];
        PrevIdx = tree_data->SmileIndex[t-1];
        
        a1 = tree_data->A1[Idx];
        a2 = tree_data->A2[Idx];
        a3 = tree_data->A3[Idx];
        
        IsLognormal = ((a2 < A2CUTOFF) && (fabs(a1) < TINY));

        SameIdx    = (Idx == PrevIdx);
        SameParams = ((fabs(a1 - tree_data->A1[PrevIdx]) < TINY) &&
                      (fabs(a2 - tree_data->A2[PrevIdx]) < TINY) &&
                      (fabs(a3 - tree_data->A3[PrevIdx]) < TINY));

        NoChangeInK = (SameIdx || SameParams);
        ChangeInK   = !NoChangeInK;
      
        /* to improve running time, treat log-normal EX seperately */

        if(IsLognormal && NoChangeInK)
        {
            for (i = Bottom1; i <= Top1; i++)    /* First dimension: foreign rate */
            {
                /* Drift of first dim including mean reversion */
                Pi = Drift1 + (d[0] - Beta1) * i; 
                
                /* Part of 2nd and 3rd dim drifts which depend on 1st dim */  
                Qi = (d[1]-Beta2)*i + Pi*JRatio1;
                Ri = Drift4 + d[3]*i + Pi*JRatio3;
                
                /* And drift of the foreign currency without CUPS */
                Si = (d[0] - Beta1) * i;
                
                l = NEAR_INT(Si);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift4L[i] = l;  
                Si -= l;
                
                sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
                sL[i].d = sd = su - Si;
                sL[i].m = s0 = 1. - su - sd;
                
                l = NEAR_INT(Pi);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift1L[i] = l;  
                Pi -= l;
                
                pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
                pL[i].d = pd = pu - Pi;
                pL[i].m = p0 = 1. - pu - pd;
                
                mMax =            OutTop2[i+l-1];
                mMax = MIN (mMax, OutTop2[i+l  ]);
                mMax = MIN (mMax, OutTop2[i+l+1]) - 1;
                
                mMin =            OutBottom2[i+l-1];
                mMin = MAX (mMin, OutBottom2[i+l  ]);
                mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;
                
                if (mMin > mMax)
                {
                    DR_Error ("Hyb3_Lattice: problem in building the tree (mMin > mMax)!");
                    goto RETURN;
                }
                
                /* Initialise 1-D pointers to beginning of 2-D blocks */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                qL = dev_data->q + offset;

                quuL = dev_data->quu + offset; 
                qu0L = dev_data->qu0 + offset; 
                qudL = dev_data->qud + offset;
                q0uL = dev_data->q0u + offset; 
                q00L = dev_data->q00 + offset; 
                q0dL = dev_data->q0d + offset;
                qduL = dev_data->qdu + offset; 
                qd0L = dev_data->qd0 + offset; 
                qddL = dev_data->qdd + offset;
                
                Shift2L = dev_data->Shift2 + offset;
                Discount_2D[CvDiscD] = dev_data->Discount_2D[CvDiscD] + offset;

                /* Prepare interest rate differential for drift */

                LnRatiosi = Drift3;
                
                if (EquityOff)
                {
                    LnRatiosi += log(Discount_1D[CvDiscF][i]) / Jump[5];
                }
                else
                {
                    if (EquityFor)
                    {
                        LnRatiosi -= log(Discount_1D[CvDiscF][i]) / Jump[5];
                    }
                }

                for (j = Bottom2[i]; j <= Top2[i]; j++) 
                {
                    Qij = Qi + (d[2] - Beta3) * j;
                    Rij = Ri + d[4]*j + Qij*JRatio4;
                    
                    m = NEAR_INT(Qij);
                    m = MINMAX (mMin - j, m, mMax - j);
                    Shift2L[j] = m; 
                    Qij -= m;
                    
                    qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                    qL[j].d = qd = qu - Qij;
                    qL[j].m = q0 = 1. - qu - qd;
                    
                    quuL[j] = pu * qu; qu0L[j] = pu * q0; qudL[j] = pu * qd;
                    q0uL[j] = p0 * qu; q00L[j] = p0 * q0; q0dL[j] = p0 * qd;
                    qduL[j] = pd * qu; qd0L[j] = pd * q0; qddL[j] = pd * qd;
                    
                    nMax =            OutTop3[i+l-1][j+m-1];
                    nMax = MIN (nMax, OutTop3[i+l-1][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l-1][j+m+1]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m-1]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m+1]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m-1]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m+1]) - 1;
                    
                    nMin =            OutBottom3[i+l-1][j+m-1];
                    nMin = MAX (nMin, OutBottom3[i+l-1][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l-1][j+m+1]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m-1]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m+1]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m-1]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m+1]) + 1;
                    
                    
                    if (nMin > nMax)
                    {
                        DR_Error ("Hyb3_Lattice: problem in building the tree (nMin > nMax)!");
                        goto RETURN;
                    }
                    
                    
                    /* Initialise 1-D pointers to beginning of 3-D blocks */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    rL = dev_data->r + offset;
                    
                    Shift3L = dev_data->Shift3   + offset;
                    ExL     = ExSpot             + offset;

                    /* complete interest rate differential for drift calculation */

                    LnRatiosij = LnRatiosi;

                    if (EquityOffOrEquityDom)
                    {
                        LnRatiosij -= log(Discount_2D[CvDiscD][j]) / Jump[5];
                    }

                    /* Third dim: EX */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        Drift6       = LnRatiosij - ExVar + CupsDrift;
                        Rijk         = Rij + d[5] * k  + Drift6;

                        n            = NEAR_INT(Rijk);
                        n            = MINMAX (nMin - k, n, nMax - k);
                        Shift3L[k]   = n;
                        Rijk        -= n;

                        rL[k].u = ru = .5 * (JumpCoeff + Rijk + Rijk * Rijk);
                        rL[k].d = rd = ru - Rijk;
                        rL[k].m = r0 = 1. - ru - rd;

                    } /* for k */
                }  /* for j */
            }  /* for i */
        }/* if log-normal EX */
        else
        {
            double MixJRatio3 = - PrevJump[3]/Jump[5];
            double MixJRatio4 = - PrevJump[4]/Jump[5];
            double MixJRatio5 = - PrevJump[5]/Jump[5];

            double Kchange = - ExMidNode[t] / Jump[5];
            double Kchangei=0,Kchangeij=0,Kchangeijk;

            for (i = Bottom1; i <= Top1; i++)    /* First dimension: foreign rate */
            {
                if (ChangeInK)
                {
                    Kchangei = Kchange + i*MixJRatio3;
                }

                /* to avoid warnings*/
                InnerTop2i = InnerBottom2i = 0;

                iterate1D = 0;
                if((i <= InnerTop1) && (i>= InnerBottom1))
                {
                    iterate1D = 1;
                    InnerTop2i     = InnerTop2[i];
                    InnerBottom2i = InnerBottom2[i];
                }

                /* Drift of first dim including mean reversion */
                Pi = Drift1 + (d[0] - Beta1) * i; 
                
                /* Part of 2nd and 3rd dim drifts which depend on 1st dim */  
                Qi = (d[1]-Beta2)*i + Pi*JRatio1;
                Ri = Drift4 + d[3]*i + Pi*JRatio3;
                
                /* And drift of the foreign currency without CUPS */
                Si = (d[0] - Beta1) * i;
                
                l = NEAR_INT(Si);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift4L[i] = l;  
                Si -= l;
                
                sL[i].u = su = .5 * (JumpCoeff + Si + Si*Si);
                sL[i].d = sd = su - Si;
                sL[i].m = s0 = 1. - su - sd;
                
                l = NEAR_INT(Pi);          /* Node shift */
                l = MINMAX (lMin - i , l, lMax - i);
                Shift1L[i] = l;  
                Pi -= l;
                
                pL[i].u = pu = .5 * (JumpCoeff + Pi + Pi*Pi);
                pL[i].d = pd = pu - Pi;
                pL[i].m = p0 = 1. - pu - pd;
                
                mMax =            OutTop2[i+l-1];
                mMax = MIN (mMax, OutTop2[i+l  ]);
                mMax = MIN (mMax, OutTop2[i+l+1]) - 1;
                
                mMin =            OutBottom2[i+l-1];
                mMin = MAX (mMin, OutBottom2[i+l  ]);
                mMin = MAX (mMin, OutBottom2[i+l+1]) + 1;
                
                if (mMin > mMax)
                {
                    DR_Error ("Hyb3_Lattice: problem in building the tree (mMin > mMax)!");
                    goto RETURN;
                }
                
                /* Initialise 1-D pointers to beginning of 2-D blocks */
                offset = Hyb3_Node_Offset(2,i,0,t,tree_data);

                qL = dev_data->q + offset;

                quuL = dev_data->quu + offset; 
                qu0L = dev_data->qu0 + offset; 
                qudL = dev_data->qud + offset;
                q0uL = dev_data->q0u + offset; 
                q00L = dev_data->q00 + offset; 
                q0dL = dev_data->q0d + offset;
                qduL = dev_data->qdu + offset; 
                qd0L = dev_data->qd0 + offset; 
                qddL = dev_data->qdd + offset;
                
                Shift2L = dev_data->Shift2 + offset;
                Discount_2D[CvDiscD] = dev_data->Discount_2D[CvDiscD] + offset;
                
                /* Prepare Ratio of zeros for CurrFwdEx[t, t+dt] */
                /* Prepare interest rate differential for drift  */

                LnRatiosi = Drift3;
                
                if (EquityOff)
                {
                    ZeroRatioi = Discount_1D[CvDiscF][i];
                    LnRatiosi += log(Discount_1D[CvDiscF][i]) / Jump[5];
                }
                else
                {
                    if (EquityDom)
                    {
                        ZeroRatioi = 1;
                    }
                    else
                    {
                        ZeroRatioi = 1 / Discount_1D[CvDiscF][i];
                        LnRatiosi -= log(Discount_1D[CvDiscF][i]) / Jump[5];
                    }
                }
                
                for (j = Bottom2[i]; j <= Top2[i]; j++) 
                {
                    if (ChangeInK)
                    {
                        Kchangeij = Kchangei + j*MixJRatio4;
                    }

                    /* to avoid warnings*/
                    InnerTop3ij     = InnerBottom3ij = 0;

                    iterate2D = 0;
                    if(iterate1D)
                    {
                        if((j <= InnerTop2i) && (j >= InnerBottom2i))
                        {
                            iterate2D = 1;
                            InnerTop3ij     = InnerTop3[i][j];
                            InnerBottom3ij  = InnerBottom3[i][j];
                        }
                    }

                    Qij = Qi + (d[2] - Beta3) * j;
                    Rij = Ri + d[4]*j + Qij*JRatio4;
                    
                    m = NEAR_INT(Qij);
                    m = MINMAX (mMin - j, m, mMax - j);
                    Shift2L[j] = m; 
                    Qij -= m;
                    
                    qL[j].u = qu = .5 * (JumpCoeff + Qij + Qij*Qij);
                    qL[j].d = qd = qu - Qij;
                    qL[j].m = q0 = 1. - qu - qd;
                    
                    quuL[j] = Quu = pu * qu; qu0L[j] = Qu0 = pu * q0; qudL[j] = Qud = pu * qd;
                    q0uL[j] = Q0u = p0 * qu; q00L[j] = Q00 = p0 * q0; q0dL[j] = Q0d = p0 * qd;
                    qduL[j] = Qdu = pd * qu; qd0L[j] = Qd0 = pd * q0; qddL[j] = Qdd = pd * qd;
                    
                    nMax =            OutTop3[i+l-1][j+m-1];
                    nMax = MIN (nMax, OutTop3[i+l-1][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l-1][j+m+1]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m-1]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l  ][j+m+1]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m-1]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m  ]);
                    nMax = MIN (nMax, OutTop3[i+l+1][j+m+1]) - 1;
                    
                    nMin =            OutBottom3[i+l-1][j+m-1];
                    nMin = MAX (nMin, OutBottom3[i+l-1][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l-1][j+m+1]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m-1]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l  ][j+m+1]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m-1]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m  ]);
                    nMin = MAX (nMin, OutBottom3[i+l+1][j+m+1]) + 1;
                    
                    if (nMin > nMax)
                    {
                        DR_Error ("Hyb3_Lattice: problem in building the tree (nMin > nMax)!");
                        goto RETURN;
                    }
                    
                    /* Initialise 1-D pointers to beginning of 3-D blocks */
                    offset = Hyb3_Node_Offset(3,i,j,t,tree_data);
                    rL = dev_data->r + offset;
                    
                    Shift3L     = dev_data->Shift3      + offset;
                    ExL         = ExSpot                + offset;
                    gDashL      = dev_data->gDash       + offset;
                    kDashTimesxL= dev_data->kDashTimesX + offset;
                    kVarL       = dev_data->kVar        + offset;

                    
                    /* if outside inner tree on IR dims, no need to access NextEx slices */
                    if(iterate1D && iterate2D)
                    {
                        NextEx00 = NextExSpot + Hyb3_Node_Offset(3,i+l+1,j+m+1,t+1,tree_data);
                        NextEx01 = NextExSpot + Hyb3_Node_Offset(3,i+l+1,j+m  ,t+1,tree_data);
                        NextEx02 = NextExSpot + Hyb3_Node_Offset(3,i+l+1,j+m-1,t+1,tree_data);
                        NextEx10 = NextExSpot + Hyb3_Node_Offset(3,i+l  ,j+m+1,t+1,tree_data);
                        NextEx11 = NextExSpot + Hyb3_Node_Offset(3,i+l  ,j+m  ,t+1,tree_data);
                        NextEx12 = NextExSpot + Hyb3_Node_Offset(3,i+l  ,j+m-1,t+1,tree_data);
                        NextEx20 = NextExSpot + Hyb3_Node_Offset(3,i+l-1,j+m+1,t+1,tree_data);
                        NextEx21 = NextExSpot + Hyb3_Node_Offset(3,i+l-1,j+m  ,t+1,tree_data);
                        NextEx22 = NextExSpot + Hyb3_Node_Offset(3,i+l-1,j+m-1,t+1,tree_data);
                    }
                    
                    /* complete ratio of zeros for FwdEx[t,t+dt]                 */
                    /* complete interest rate differential for drift calculation */

                    ZeroRatioij = ZeroRatioi;
                    LnRatiosij  = LnRatiosi;

                    if (EquityOffOrEquityDom)
                    {
                        ZeroRatioij /= Discount_2D[CvDiscD][j];
                        LnRatiosij  -= log(Discount_2D[CvDiscD][j]) / Jump[5];
                    }
                    
                    /* Third dim: EX */
                    for (k = Bottom3[i][j]; k <= Top3[i][j]; k++)
                    {
                        iterate3D = 0;
                        if(iterate1D && iterate2D)
                        {
                            if((k <= InnerTop3ij) && (k >= InnerBottom3ij))
                            {
                                iterate3D = 1;
                            }
                        }
                     
                        Drift6 = kDashTimesxL[k] * LnRatiosij - gDashL[k] * ExVar + CupsDrift;


                        if (ChangeInK)
                        {
                            Kchangeijk = Kchangeij + k*MixJRatio5 + kVarL[k]/Jump[5];
                            Drift6    += Kchangeijk;
                        }

                        CurrFwd = ExL[k] * DeflatedFwdRatio * ZeroRatioij
                                         * exp(CupsDriftRaw / kDashTimesxL[k]);

                        Rijk    = Rij + d[5] * k + Drift6;

                        n       = NEAR_INT(Rijk);


                        /********************************************/
                        /* initialise variables used to iteratively */
                        /* improve the fit to FwdFX                 */
                        /********************************************/
                        
                        nGuess = currn = nextn = n;
                        RGuess = RLoc  = Rijk;

                        /*************************************************/
                        /* Try to improve RGuess using only              */
                        /* existing grid points.                         */
                        /* if at any iteration:                          */
                        /*    1- no solution found or                    */
                        /*    2- the solution takes us outside limits    */
                        /* THEN revert to RGuess                         */
                        /*************************************************/

                        if (EquityOffOrEquityDom)
                        {
                            for(count = 1 ; count <= (int) (MAX_FX_ITER * iterate3D) ; count++)
                            {

                                k1 = k  + currn;
                                k0 = k1 + 1;
                                k2 = k1 - 1;
                        
                                if((k1 > nMax) || (k1 < nMin))
                                {  
                                    nextn = currn - 1L;
                                    break;
                                }

                                /* Coeffs for first moment equation */

                                c1 = Quu * NextEx00[k0] + Qu0 * NextEx01[k0] + Qud * NextEx02[k0] +
                                     Q0u * NextEx10[k0] + Q00 * NextEx11[k0] + Q0d * NextEx12[k0] +
                                     Qdu * NextEx20[k0] + Qd0 * NextEx21[k0] + Qdd * NextEx22[k0];
                        
                        
                                c2 = Quu * NextEx00[k1] + Qu0 * NextEx01[k1] + Qud * NextEx02[k1] +
                                     Q0u * NextEx10[k1] + Q00 * NextEx11[k1] + Q0d * NextEx12[k1] +
                                     Qdu * NextEx20[k1] + Qd0 * NextEx21[k1] + Qdd * NextEx22[k1];
                        
                        
                                c3 = Quu * NextEx00[k2] + Qu0 * NextEx01[k2] + Qud * NextEx02[k2] +
                                     Q0u * NextEx10[k2] + Q00 * NextEx11[k2] + Q0d * NextEx12[k2] +
                                     Qdu * NextEx20[k2] + Qd0 * NextEx21[k2] + Qdd * NextEx22[k2];
                        
                                if(Hyb3_SolveForShift(&RLoc,&nextn,currn,c1,c2,c3,JumpCoeff,CurrFwd,&NoSolution) == FAILURE)
                                {
                                    DR_Error("pb solving shift\n");
                                    goto RETURN;
                                }
                            
                        
                                if(NoSolution)
                                {
                                    nextn = currn - 1L;
                                    break;
                                }
                                else
                                {
                        
                                    if (nextn == currn)
                                    {
                                        break;
                                    }
                                    else
                                    {   
                                        temp  = currn;
                                        currn = nextn;
                                        nextn = temp;
                                    }
                                }
                            } /* for count */

                            /******************************************************/
                            /* IF we exit the loop                                */
                            /* 1- without having found a solution OR              */
                            /* 2- the solution takes us outside the initial nodes */
                            /* THEN revert to RGuess                              */
                            /******************************************************/

                            if(nextn != currn)
                            {
                                RLoc  = RGuess;
                                nextn = nGuess;
                            }
                            else
                            {
                                if(nextn != nGuess)
                                {
                                    nextn = nGuess;
                                    RLoc  = RGuess;
                                }  
                            }
                        }

                        nextn      = MINMAX (nMin - k, nextn, nMax - k);
                        Rijk       = RLoc - nextn;
                        Shift3L[k] = nextn;

                        rL[k].u = ru = .5 * (JumpCoeff + Rijk + Rijk * Rijk);
                        rL[k].d = rd = ru - Rijk;
                        rL[k].m = r0 = 1. - ru - rd;

                    } /* for k */
                }  /* for j */
            }  /* for i */            
        }/* EX NOT log-normal */

        if(ChecksOn)
        {
            if(Hyb3_Dev(dev_data->FwdFX,
                   t,
                   T,
                   CvDiscD,
                   DISC_3D_CUPS,
                   dev_data,
                   tree_data) == FAILURE) goto RETURN;


            if(Hyb3_Dev(dev_data->DomZero,
                   t,
                   T,
                   CvDiscD,
                   DISC_2D_CUPS,
                   dev_data,
                   tree_data) == FAILURE) goto RETURN;


            if(t==0)
            {
                int     Offset3D;
                int     Offset2D;
                char    ErrorMsg[MAXBUFF];
                double  TreeForZero,DetermFZero;
                double  TreeDomZero,DetermDZero;

                Offset3D = Hyb3_Node_Offset(3,0,0,t,tree_data);
                Offset2D = Hyb3_Node_Offset(2,0,0,t,tree_data);

                TreeForZero = (dev_data->FwdFX+Offset3D)[0]/((dev_data->FxSpot + Offset3D)[0]);
                DetermFZero = tree_data->ZeroCoupon[0][tree_data->CvDisc[0]][T];

                TreeDomZero = (dev_data->DomZero + Offset2D)[0];
                DetermDZero = tree_data->ZeroCoupon[1][tree_data->CvDisc[1]][T];
         
                
                if(fabs(TreeDomZero - DetermDZero) > DOMZERO_TOL)
                {
                    sprintf(ErrorMsg,"DomDisc calibration failed: Tree = %lf ; Determ = %lf\n",
                            TreeDomZero,DetermDZero);

                    DR_Error(ErrorMsg);
                }

                if(fabs(TreeForZero - DetermFZero) > FWDFX_TOL)
                {					
                    sprintf(ErrorMsg,"ForDisc calibration failed: Tree = %lf ; Determ = %lf\n",
                            TreeForZero,DetermFZero);

                    DR_Error(ErrorMsg);
                }
            }/* if t==0 */
        }/* if Price checking is ON */
    } /* End if (CUPS only) else if (EQ+IR) else 3-D mode */
    
    




    status = SUCCESS;

RETURN:

    return (status);

}  /* Hyb3_Lattice */


