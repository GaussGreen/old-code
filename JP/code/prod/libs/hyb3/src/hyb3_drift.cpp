/****************************************************************************/
/*      Building of the interest rate tree: calculation of forward rates,   */
/*      interest rate drift at each time step.                      */
/****************************************************************************/
/*      DRIFT.c                                                             */
/****************************************************************************/


/*
$Header$
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

#define  ADJSMOOTH_FAC  0.5     /* to smooth the adjusted CUPS drift         */
#define  CUPSADJ_CUTOFF 1.0     /* cutoff percentage for adjusted CUPS drift */
#define  MINYEARSTEP    1.0     /* minimum step between points to smooth     */


/******** Hyb3_SmoothCurve *************************************************/
/**
        given a curve crv, the function returns a smoothed approximation
        to it (in terms of cubic spline)  
        (1) assumes smootheCurve to be allocated
        (2) assumes curve to be allocated on 0..nPts (included)
 ***************************************************************************/

int Hyb3_SmoothCurve(double*     smoothCrv,  /**<(O) Smoothed Curve             */
                     double*     crvT,       /**<(I) input curve times          */
                     double*     crvPt,      /**<(I) input curve to be smoothed */
                     long        nPts,       /**<(I) number of initial crv pts  */
                     int*        index,      /**<(I) index for smile parameters ("critical points") */
                     long        step)       /**<(I) number of steps between smoothing pts */
{
    int status = FAILURE;
    int  i, n, nSPL =0;
    double inder, finder ,a, b, h;
    double * AcurveT = NULL;
    double * AcurveVal = NULL;
    double * AcurveSPL = NULL;

    AcurveT   = (double *)DR_Array(DOUBLE, 0, nPts);
    AcurveVal = (double *)DR_Array(DOUBLE, 0, nPts);
    AcurveSPL = (double *)DR_Array(DOUBLE, 0, nPts);
    if ((AcurveT   == NULL)||
        (AcurveVal == NULL)||
        (AcurveSPL == NULL)) 
    {

        DR_Error("Hyb3_SmoothCurve:: Memory allocation failure");
        goto RETURN;
    }

    for (i = 0 ; i <= nPts ; ++i )
    {
        if ( ( (i % step) == 0) || (index[i] != index[i-1] ) )
        {
            if ( (i>0) && (crvT[i] - AcurveT[nSPL-1]) < MINYEARSTEP ) continue; /* min distance */
            AcurveT[nSPL]   = crvT[i];
            AcurveVal[nSPL] = crvPt[i];
            nSPL++;
        }        
    }
    /* add the final point to the spline curve if necessary */
    if ( fabs( AcurveT[nSPL-1] - crvT[nPts]) > TINY )
    {
        AcurveT[nSPL]   = crvT[nPts];
        AcurveVal[nSPL] = crvPt[nPts];
        nSPL++;
    }

    /*Calculation of initial and final derivative and initialising spline */
    inder  = (crvPt[1] - crvPt[0])/(crvT[1]-crvT[0]);
    finder  = (crvPt[nPts] - crvPt[nPts-1])/(crvT[nPts]-crvT[nPts-1]);
    
    if(SplineInterp1dInit(AcurveT, AcurveVal, nSPL, inder, finder, AcurveSPL )!=SUCCESS)
    {
        DR_Error("Hyb3_SmoothCurve:: Pb in Spline Initialisation in smoothing function");
        goto RETURN;
    }
    
    n = 0;
    for (i= 0 ; i <=nPts; ++i)
    {
        if ( (crvT[i]-AcurveT[n+1]) > TINY) n++;

        h = AcurveT[n+1] - AcurveT[n];
        a = (AcurveT[n+1] - crvT[i])/h;
        b = 1.0 - a;

        smoothCrv[i] = a*AcurveVal[n]+b*AcurveVal[n+1]+((a*a*a-a)*AcurveSPL[n]+
                       (b*b*b-b)*AcurveSPL[n+1])*(h*h)/6.0;
    }
    
    /* finally, interpolate the scaling function smoothly */
    status = SUCCESS;
RETURN:

    Free_DR_Array(AcurveT,   INT, 0, nPts);
    Free_DR_Array(AcurveVal, INT, 0, nPts);
    Free_DR_Array(AcurveSPL, INT, 0, nPts);

    return status;
    
}






/******** Hyb3_SimpleFwdCurve****************************************************/
/**
  
       Given NbTP points, it outputs (NbTP - 1) "Forward rates * TimePeriod",
       where FwdRate[i] spans (TPDate[i],TPDate[i+1]).
       
  
  
  
 ***************************************************************************/

int Hyb3_SimpleFwdCurve(
            double  *FwdRate,       /**<(O) Interpolated Fwd rates.           */
            T_CURVE const* crv,
            long    *TPDate,        /**<(I) date of time point                */
            int     NbTP)           /**<(I)total number of time points        */
{
    double z1;      /* discount factor at start of period*/
    double z2;      /* discount factor at end of period  */
    int     i;

    if (NbTP < 1)
    {
        DR_Error("invalid input Nb of time points "
                    "in Hyb3_SimpleFwdCurve: should be at least 1\n");
        return FAILURE;
    }

    if (FwdRate == NULL || TPDate == NULL)
    {   
        DR_Error("Invalid pointer inputs to Hyb3_SimpleFwdCurve\n");
        return FAILURE;
    }

    z1 = GetZeroPrice(TPDate[0], crv);
    for (i = 0; i < NbTP - 1; i++)
    {
        z2 = GetZeroPrice(TPDate[i+1], crv);
        if (z1 < TINY || z2 < TINY) 
            return FAILURE;

        FwdRate[i] = z1/z2 - 1.0;
        z1 = z2;
    }
    return SUCCESS;    
}





/*****  Hyb3_Forward_Curve  *******************************************************/
/**
 *   Interpolate  zero and  calculate  forward  at each  node point i.e.
 *   every Length[i] interval.When the value date is not the same as the 
 *   zero curve value date  we need to interpolate  or  even extrapolate 
 *   backward.
 */
int    	Hyb3_Forward_Curve (
        double*         ZeroCoupon, /* (O) Zero coupon price       */
        double*         ZeroRate,   /* (O) Zero coupon rate        */
        double*         FwdRate,    /* (O) Forward rate            */
        T_CURVE const* crv,   /* (I) Zero curve              */
        long const*     TPDate,     /* (I) Date of each time point */
        int             NbTP)       /* (I) Nb of time points       */
{
    int     i;
    double  inv_mat;            /* 1/Maturity of zero coupon bond   */
    int     status = FAILURE;   /* Error status = FAILURE initially */

    double  ZPrice_i,           /* ZeroPrice from ZValueDate to TPDate[i]   */
            ZPrice_i1;          /* ZeroPrice from ZValueDate to TPDate[i+1] */

    if (NbTP <= 0L) goto RETURN;

    /* compute the 1-period forwards */
    ZPrice_i = GetZeroPrice(TPDate[0], crv);
    if (ZPrice_i < TINY) goto RETURN;

    for (i = 0; i <= NbTP; i++)
    {
        ZPrice_i1 = GetZeroPrice(TPDate[i+1], crv);
        if (ZPrice_i1 < TINY) goto RETURN;

        FwdRate[i] = ZPrice_i/ZPrice_i1 - 1.0;

        ZPrice_i = ZPrice_i1;
    }

    /* compute ZeroCoupon price and ZeroRate that are based from TPDate[0] */
    ZeroCoupon[0] = 1.0;
    for (i = 0; i <= NbTP; i++)
    {
        ZeroCoupon[i+1] = ZeroCoupon[i] / (1.0 + FwdRate[i]);

        inv_mat  = 365.0 / Daysact (TPDate[0], TPDate[i+1]);
        ZeroRate[i+1] = pow(ZeroCoupon[i+1], -inv_mat) - 1.0;
    }
    ZeroRate[0] = ZeroRate[1];

    status = SUCCESS;

RETURN:

    if (status == FAILURE)
        DR_Error("Hyb3_Forward_Curve: failed.");
    return (status);

}




/*****  Hyb3_Find_Drift_1D   *******************************************************/
/**
 *       Calculate the interest rate drift in terms of positioning the centre
 *       of the tree.
 *
 *       We take Ron Levin's method and solve for the drift directly,  rather
 *       than calibrate the probabilities.
 */

int   Hyb3_Find_Drift_1D(double    *ZCenter,   /**< (O) Middle node                 */
                    double    *ZeroCoupon,/**< (I) Zero prices (diffused curve)*/
                    double    *FwdRate,   /**< (I) Forward rate                */
                    double     QLeft,
                    double     QRight,
                    double     FwdShift,
                    double     Beta,      /**< (I) Mean reversion coefficient  */
                    double    *Hyb3_SpotVol,   /**< (I) Spot vol of IR at each node */
                    int       *Top,       /**< (I) Upper limit of the tree     */
                    int       *Bottom,    /**< (I) Lower limit of the tree     */
                    int       *OutTop,    /**< (I) Top of outer ellipse        */
                    int       *OutBottom, /**< (I) Bottom of outer ellipse     */
                    int        NbTP,      /**< (I) Total number of time points */
                    double    *Length,    /**< (I) Length of t steps (ACT/365) */
                    double    *LengthJ,   /**< (I) Length for jump size calc   */
                    HYB3_TREE_DATA *tree_data)
{

 

    

    TSLICE   StatePr;          /* State Prices at t                    */
    TSLICE   Discount;         /* 1 period discount factor             */
    TSLICE   StatePr1;         /* State Prices at t+1                  */

    double  *StatePrL;         /* Local slice pointers for convenience */
    double  *DiscountL;        /*     in addressing the slices         */
    double  *StatePr1L;
    double  *EDevPriceL;       /* St prices kept for express DEV'ing   */



    double  Jump;              /* Size of the jump (in log space)           */
    double  PreviousJump;      /* Jump size at previous period              */
    double  RateJumpLeft;
    double  RateJumpRight;
    double  JumpCoeff, du;     /* Jump coefficients                         */
    double  Grid;              /* Current grid point                        */

    double  BetaLocal;         /* = Beta * Length[t]                        */
    double  Pi;                /* Drift at the current node                 */
    double  d;                 /* Drift due to the change in jump size      */
    double  Zt=0;              /* Center of the tree in X-space             */
    double  DelZt;             /* Correction to Zt for current time point   */
    double  Zidx;
    double  P[3];              /* 2rd degree polynomial coefficients        */
    double  p;                 /* Total probability                         */
    double  x;                 /*                                           */
    
    double  FwdRateA;          /* Fwd rate adjusted                         */
    double  MLeft, SLeft;      /* Multiple coeff and shift for grid point   */
    double  MRight, SRight;

    double  VolBbq;            /* Sigma used in bone mapping                */
    double  QMid;              /* Average q parameter                       */
    double  QSh;               /* q parameter for shift adjustment          */

    double  Q, U;              /* Q and U weight coeff for Z-polynomial     */ 
    double  D0;                /* Consecutive powers of discount factor     */ 
    
    double  aL00;               /* Coeff of Z-polynomial across all states   */
    double  aL10, aL11 ;
    double  aL20, aL21, aL22 ;
    double  aR00;               /* Coeff of Z-polynomial across all states   */
    double  aR10, aR11 ;
    double  aR20, aR21, aR22 ;
    
    double  pu, pd, p0;        /* Probabilities                             */

    int     EDevOn = FALSE;    /* To indicate whether EDev tool is used     */
    int     EDevIdx = 0;       /* Counter of dates for express DEV          */

    int     i;                 /* Node index                                */
    int     l, lMin, lMax;     /* Node branching shift                      */
    int     Mid;
    int     t;                 /* Time step index                           */
    int     offset;            /* For efficient slice addressing            */
    int     status = FAILURE;  /* Error status = FAILURE initially          */



    StatePr  = Hyb3_Alloc_Slice (tree_data, 1);
    Discount = Hyb3_Alloc_Slice (tree_data, 1);
    StatePr1 = Hyb3_Alloc_Slice (tree_data, 1);

    if (   (StatePr  == NULL)
        || (Discount == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Hyb3_Drift_1D: could not allocate memory!");
        goto RETURN;        
    }


    /* Check if Express DEV tool is being used */
    if (tree_data->NbEDevDates > 0)
    {

        /* Make sure we are running 1-D mode */
        if (tree_data->TreeType != TTYPE_1IR)
        {
            DR_Error("Hyb3_Drift_1D: Express DEV tool only available for 1-D CET "
                     "runs.\n");
            goto RETURN;
        }

        EDevOn = TRUE;
        EDevIdx = 0;

        /* Allocate array of pointers to state price slices */
        tree_data->EDevStPrice = (TSLICE *)DR_Array(DOUBLE_PTR, 
                                                    0, 
                                                    tree_data->NbEDevDates-1);
        if (tree_data->EDevStPrice == NULL)
        {
            goto RETURN;
        }


        /* Initialise pointers to NULL for safe freeing */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            tree_data->EDevStPrice[i] = NULL;
        }

        /* Finally allocate the slices themselves in preparation */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            /* Check on TreeType has ensured that the current run */
            /* is in 1-D mode, so 1-factor allocation.            */
            tree_data->EDevStPrice[i] = Hyb3_Alloc_Slice(tree_data, 1);
            if (tree_data->EDevStPrice[i] == NULL)
            {
                goto RETURN;
            }
        }

    }


    QMid = (QLeft + QRight) / 2; 

    StatePrL = StatePr + Hyb3_Node_Offset(1, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;


    for (t = 0; t <= NbTP; t++)
    {   

        /* Prepare for addressing */
        StatePr1L = StatePr1 + Hyb3_Node_Offset(1, 0, 0, t+1, tree_data);
        StatePrL  = StatePr  + Hyb3_Node_Offset(1, 0, 0, t,   tree_data);
        DiscountL = Discount + Hyb3_Node_Offset(1, 0, 0, t,   tree_data);


        /*  Precompute jumps and probabilities coefficients */
        du           = sqrt (JUMPCOEFF * LengthJ[t-1]) ;
        PreviousJump = Hyb3_SpotVol[t-1] * du ;
                                                                                
        du   = sqrt (JUMPCOEFF * LengthJ[t]);
        Jump = Hyb3_SpotVol[t] * du ;
        d    = PreviousJump / Jump - 1.;
        
        BetaLocal = Beta * Length[t] * PreviousJump / Jump;

        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        FwdRateA  = FwdRate[t] / (1. + FwdShift);

        VolBbq = (1. + FwdShift) / (1. + QMid * FwdShift);


        
        MLeft        = MRight        = FwdRateA;
        SLeft        = SRight        = 1 + FwdRateA;
        RateJumpLeft = RateJumpRight = PreviousJump * VolBbq * FwdRateA;
        if (fabs(QLeft) > QCUTOFF)
        {
            MLeft       /= QLeft;
            SLeft       -= FwdRateA / QLeft;      
            RateJumpLeft = exp (QLeft * VolBbq * PreviousJump);
        }
        if (fabs(QRight) > QCUTOFF)
        {
            MRight       /= QRight;
            SRight       -= FwdRateA / QRight;      
            RateJumpRight = exp (QRight * VolBbq * PreviousJump);
        }


        /* Set initial Zt */
        if (t == 0)
        {
            QSh = (FwdShift > 0) ? QRight : QLeft;
            if (fabs(QSh) > QCUTOFF)
            {
                Zt = log(1. + QSh * FwdShift) / (QSh * VolBbq);
            }
            else
            {
                Zt = FwdShift / VolBbq;
            }
        }


        /*   Set up the grid points  */
        /* LEFT part of distribution */
        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid   = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);
        Zidx += (PreviousJump) * Bottom[t];

        if (fabs(QLeft) > QCUTOFF)
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid); 
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */
        Zidx = Zt  + (PreviousJump) * (Mid + 1);

        if (fabs(QRight) > QCUTOFF)
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid); 
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;
            for (i = Mid +1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);     
                Grid += RateJumpRight;
            }
        }

        

        /*  Calculate polynomial coefficients and drift  */
        for (i = 0; i < 3; i++)                     
        {
            P[i] = 0.; 
        }


        Q = QLeft * VolBbq;
        U = (FwdRateA * (1. - QLeft) - QLeft) * VolBbq; 

        aL00 =  1.;
        aL10 = -Q;          aL11 = -U;
        aL20 =  0.5*Q*Q;    aL21 =  1.5*U*Q;    aL22 =   U*U;

        Q = QRight * VolBbq;
        U = (FwdRateA * (1. - QRight) - QRight) * VolBbq; 
    
        aR00 =  1.;
        aR10 = -Q;          aR11 = -U;
        aR20 =  0.5*Q*Q;    aR21 =  1.5*U*Q;    aR22 =   U*U;

        Zidx = Zt; 
        Mid  = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid  = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);

        /* From BOTTOM to MID */
        for (i = Bottom[t]; i <= Mid; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aL10 + aL11 * x) * D0;

            P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  
        }

        /* From MID to TOP */
        for (i = Mid + 1; i <= Top[t]; i++)
        {
            x   = DiscountL[i];
            D0  = StatePrL[i] * x;
            P[0] += D0;

            P[1] += (aR10 + aR11 * x) * D0;

            P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  
        }


        P[0] -= ZeroCoupon[t+1];


        /* Newton-Raphson to find the root of the polynomial. */
        /* We use the previous value as a first guess.        */

        DelZt = NR_Poly (0., P, 2);

        if (fabs(DelZt * VolBbq) > 1.)
        {
            DR_Error ("Problem in calculation of drift (Hyb3_Drift_1D)!");
            goto RETURN;
        }

        Zt += DelZt;
        ZCenter[t] = Zt;

        /* OLD CODE */
        /* MidNode[t] = FwdRate[t] * exp(Zt);*/


        /*   Recalculate the grid using the middle node  */
        /* LEFT part of distribution */
        Zidx  = Zt; 
        Mid   = (int) ceil(-Zidx / PreviousJump) - 1;
        Mid   = MIN ( MAX (Mid, Bottom[t] - 1), Top[t]);
        Zidx += (PreviousJump) * Bottom[t];

        if (fabs(QLeft) > QCUTOFF)
        {
            Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);        
                Grid *= RateJumpLeft;
            }           
        }
        else
        {
            Grid = MLeft * VolBbq * Zidx;
            for (i = Bottom[t]; i <= Mid; i++)
            {
                DiscountL[i] = 1. / (SLeft + Grid);
                Grid += RateJumpLeft;
            }
        }

        /* Right part of distribution */
        Zidx = Zt + (PreviousJump) * (Mid + 1);

        if (fabs(QRight) > QCUTOFF)
        {
            Grid = MRight * exp (QRight * VolBbq * Zidx);
            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                Grid *= RateJumpRight;
            }           
        }
        else
        {
            Grid = MRight * VolBbq * Zidx;
            for (i = Mid + 1; i <= Top[t]; i++)
            {
                DiscountL[i] = 1. / (SRight + Grid);
                Grid += RateJumpRight;
            }
        }
       


        /*  Compute state prices for next period  */

        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {
            StatePr1L[i] = 0.;                        
        }
                

        /* Branching has to be inside outer ellipse at next period */
        lMax = OutTop[t+1] - 1;                                             
        lMin = OutBottom[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Drift_1D: problem in building the tree (lMin > lMax)!");
            goto RETURN;
        }


        for (i = Bottom[t]; i <= Top[t]; i++)
        {                        
            Pi = (d - BetaLocal) * i;
                                                                                
            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                                                                                
            Pi -= l;
                        
                        
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;

            x = StatePrL[i] * DiscountL[i];

            StatePr1L[i+1+l] += pu * x;
            StatePr1L[i  +l] += p0 * x;
            StatePr1L[i-1+l] += pd * x;
        }


        /*   Check that state prices add up to 1 */

        p = 0.;

        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {
            p += StatePr1L[i];
        }

        if (fabs (p - ZeroCoupon[t+1]) > 0.0001)
        {
            DR_Error ("Hyb3_Drift_1D: Drift calibration has failed.\n"
                      "          Please increase sigma or ppy.\n"
                      "          (State prices don't add up to 1) \n");
   
            goto RETURN;
        }
    



        /* Express DEV tool: store state prices if necessary   */
        if ((EDevOn) && (EDevIdx<tree_data->NbEDevDates))
        {
            if (tree_data->TPDate[t+1] == tree_data->EDevDate[EDevIdx])
            {

                /* Store corresponding state prices */
                offset = Hyb3_Node_Offset(1, 0, 0, t+1, tree_data);
                EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;

                for (i = Bottom[t+1]; i <= Top[t+1]; i++)
                {
                    EDevPriceL[i] = StatePr1L[i];
                }

                /* Increment the EDev date index */
                EDevIdx++;
            }
        }




        /*  Permute values for the next time point  */

        StatePrL = StatePr + Hyb3_Node_Offset(1, 0, 0, t+1, tree_data);

        for (i = Bottom[t+1]; i <= Top[t+1]; i++)
        {                                                               
            StatePrL[i] = StatePr1L[i];
        }
    }  /* for t */                                         

                                               
               
    /* Final check on number of EDev dates processed */
    if (EDevOn)
    {
        if (EDevIdx != tree_data->NbEDevDates)
        {
            DR_Error("Express DEV:  Error in processing  the DEV dates.\n"
                     "Please ensure that dates are ordered, that there\n"
                     "is no repetition,that value date is not included\n"
                     "and that the express DEV dates are also critical\n" 
                     "dates .\n");
            goto RETURN;
        }
    }




    status = SUCCESS;


  RETURN:

    
    Hyb3_Free_Slice (StatePr,  tree_data, 1);
    Hyb3_Free_Slice (Discount, tree_data, 1);
    Hyb3_Free_Slice (StatePr1, tree_data, 1);


    return (status);


}  /* Hyb3_Find_Drift_1D */


/*****  Hyb3_Find_Drift_2D   ***************************************************/
/*
 *       Calculate the interest rate drift for a two factor model.
 *       (from Fix3_Drift_2D_Classic)
 */
int     Hyb3_Find_Drift_2D (double    *ZCenter,   /**< (O) Middle node                 */
                    double    *ZeroCoupon,/**< (I) Zero prices (diffused curve)*/
                    double    *FwdRate,   /**< (I) Forward rate                */
                    double     QLeft,
                    double     QRight,
                    double     FwdShift,
                    double     Beta[3],      /**< (I) Mean reversion coefficient  */
                    double    **Aweight,   /**< (I) Spot vol of IR at each node */
                    int       *Top1,       /**< (I) Upper limit of the tree     */
                    int       *Bottom1,    /**< (I) Lower limit of the tree     */
                    int       *OutTop1,    /**< (I) Top of outer ellipse        */
                    int       *OutBottom1, /**< (I) Bottom of outer ellipse     */
                    int       **Top2,       /**< (I) Upper limit of the tree     */
                    int       **Bottom2,    /**< (I) Lower limit of the tree     */
                    int       **OutTop2,    /**< (I) Top of outer ellipse        */
                    int       **OutBottom2, /**< (I) Bottom of outer ellipse     */
                    int        NbTP,      /**< (I) Total number of time points */
                    double    *Length,    /**< (I) Length of t steps (ACT/365) */
                    double    *LengthJ,   /**< (I) Length for jump size calc   */
                    HYB3_TREE_DATA *tree_data)
{

    /* Local values of tree data structure */
    /* Use diffused zero and fwd curves    */


    /* Other variables */
    TSLICE  StatePr;            /* State Prices at t                         */
    TSLICE  Discount;           /* 1 period discount factor                  */
    TSLICE  StatePr1;           /* State Prices at t+1                       */

                                                                              
    double  *StatePrL;           /* Local slice pointers                      */
    double  *DiscountL;                                                       
    double  *StatePr1L;                                                       
    double  *StatePr2L;                                                       
    double  *StatePr3L;  
    double  *EDevPriceL;         /* St prices kept for express DEV'ing        */
   
     
    
    double  *TempPointer;                                                 
                                                                              
    double  Jump[3];             /* Jump sizes at current time point          */
    double  PreviousJump[3];     /* Jump sizes at previous time point         */
    double  RateJumpLeft2;       /* Rate space jump size for last dimension   */
    double  RateJumpRight2;      /* Rate space jump size for last dimension   */
    double  JumpCoeff, d[3], du; /* Jump coefficients                         */
    double  Grid;                /* Grid points                               */
                                                                              
    double  Beta1, Beta2, Beta3; /* Modified mean reversion coefficients      */
    double  Pi, Qi, Qij;         /* Total drift                               */
    double  p;                   /* Total probability                         */
    double  Zt=0;                /* Center of the tree in X-space             */
    double  DelZt;               /* Correction to Zt for current time point   */
    double  Zidx;                /* Zt index i,j,k adjusted                   */
    double  P[3];                /* 2nd degree polynomial coefficients        */
    double  x;                   /*                                           */
    double  FwdRateA;            /* Fwd rate adjusted                         */
    double  MLeft, SLeft;        /* Multiple coeff and shift for grid point   */
    double  MRight, SRight;
    double  VolBbq;              /* Sigma used in bone mapping                */
    double  QMid;                /* Average q parameter                       */
    double  QSh;                 /* q parameter for shift adjustment          */
    double  Q, U;                /* Q and U weight coeff for Z-polynomial     */ 
    double  D0;                  /* Consequtive power of discount factor      */ 
    double  aL00;                /* Coeff of Z-polynomial across all states   */
    double  aL10, aL11 ;
    double  aL20, aL21, aL22 ;
    double  aR00;               
    double  aR10, aR11 ;
    double  aR20, aR21, aR22 ;


    double  pu, pd, p0;          /* Probabilities in first dimension          */
    double  qu, qd, q0;          /* Probabilities in second dimension         */

    int     EDevOn = FALSE;      /* To indicate whether EDev tool is used     */
    int     EDevIdx = 0;         /* Counter of dates for express DEV          */
                                                                              
    int     offset;
    int     Mid;                 /* Mid of distribution index                 */
    int     i, j, j1, j2, j3;    /* Node indices                              */
    int     l, lMin, lMax;       /* Node branching shifts                     */ 
    int     m, mMin, mMax;                                                    
    int     t;                   /* Time step index                           */
    int     status = FAILURE;    /* Error status = FAILURE initially          */


    StatePr  = Hyb3_Alloc_Slice (tree_data,2); /*TWO DIMENSIONAL SLICE?*/
    Discount = Hyb3_Alloc_Slice (tree_data,2);
    StatePr1 = Hyb3_Alloc_Slice (tree_data,2);

    if (   (StatePr  == NULL)
        || (Discount == NULL)
        || (StatePr1 == NULL))
    {
        DR_Error ("Hyb3_Drift_2D: could not allocate memory!");
        goto FREE_MEM_AND_RETURN;        
    }


    /* Check if Express DEV tool is being used */
    if (tree_data->NbEDevDates > 0)
    {

        /* Make sure we are running 1-IR mode with 2 factors*/
        if (tree_data->TreeType != TTYPE_1IR2F)
        {
            DR_Error("Hyb3_Drift_2D: Express DEV tool only available for 1-IR CET "
                     "runs.\n");
            goto FREE_MEM_AND_RETURN;
        }
        EDevOn = TRUE;
        EDevIdx = 0;

        /* Allocate array of pointers to state price slices */
        tree_data->EDevStPrice = (TSLICE *)DR_Array(DOUBLE_PTR, 
                                                     0, 
                                                     tree_data->NbEDevDates-1);
        if (tree_data->EDevStPrice == NULL)
        {
            goto FREE_MEM_AND_RETURN;
        }


        /* Initialise pointers to NULL for safe freeing */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            tree_data->EDevStPrice[i] = NULL;
        }

        /* Finally allocate the slices themselves in preparation */
        for (i=0; i<tree_data->NbEDevDates; i++)
        {
            /* Check on TreeType has ensured that the current run */
            /* is in 2-D mode, so 2-factor allocation.            */
            tree_data->EDevStPrice[i] = Hyb3_Alloc_Slice(tree_data, 2);
            if (tree_data->EDevStPrice[i] == NULL)
            {
                goto FREE_MEM_AND_RETURN;
            }
        }

    }

    QMid = (QLeft + QRight) / 2;

    StatePrL = StatePr + Hyb3_Node_Offset(2, 0, 0, 0, tree_data);
    StatePrL[0] = 1.;

    for (t = 0; t < NbTP; t++)
    {                                   
        du = sqrt (JUMPCOEFF * LengthJ[t-1]);
        
        PreviousJump[0] = Aweight[0][t-1] * du;
        PreviousJump[1] = Aweight[1][t-1] * du;
        PreviousJump[2] = Aweight[2][t-1] * du;
        
        du = sqrt (JUMPCOEFF * LengthJ[t]);
            
        Jump[0] = Aweight[0][t] * du;
        Jump[1] = Aweight[1][t] * du;
        Jump[2] = Aweight[2][t] * du;
        
        d[0] = (PreviousJump[0] - Jump[0]) / Jump[0];
        d[1] = (PreviousJump[1] - Jump[1]) / Jump[2];
        d[2] = (PreviousJump[2] - Jump[2]) / Jump[2];
            
        Beta1  = Beta[0] * Length[t] * PreviousJump[0] / Jump[0];
        Beta2  = Beta[1] * Length[t] * PreviousJump[1] / Jump[2];
        Beta3  = Beta[1] * Length[t] * PreviousJump[2] / Jump[2];
        
        JumpCoeff = Length[t]/(LengthJ[t]*JUMPCOEFF);

        FwdRateA  = FwdRate[t] / (1. + FwdShift);
        VolBbq = (1. + FwdShift) / (1. + QMid * FwdShift);

        MLeft         = MRight         = FwdRateA;
        SLeft         = SRight         = 1 + FwdRateA;
        RateJumpLeft2 = RateJumpRight2 = PreviousJump[2] * VolBbq * FwdRateA;

        if (fabs(QLeft) > QCUTOFF)
        {
            MLeft        /= QLeft;
            SLeft        -= FwdRateA / QLeft;      
            RateJumpLeft2 = exp (QLeft * VolBbq * PreviousJump[2]);
        }
        if (fabs(QRight) > QCUTOFF)
        {
            MRight        /= QRight;
            SRight        -= FwdRateA / QRight;      
            RateJumpRight2 = exp (QRight * VolBbq * PreviousJump[2]);
        }

        if (t == 0)
        {
            QSh = (FwdShift > 0) ? QRight : QLeft;
            if (fabs(QSh) > QCUTOFF)
            {
                Zt = log(1. + QSh * FwdShift) / (QSh * VolBbq);
            }
            else
            {
                Zt = FwdShift / VolBbq;
            }
        }

        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt 
                  + (PreviousJump[0] + PreviousJump[1]) * i;
            Mid   = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (PreviousJump[2]) * Bottom2[t][i];

            DiscountL = Discount + Hyb3_Node_Offset (2, i, 0, t, tree_data);
    
            if (fabs(QLeft) > QCUTOFF)
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (PreviousJump[0] + PreviousJump[1]) * i
                 + (PreviousJump[2]) * (Mid + 1);

            if (fabs(QRight) > QCUTOFF)
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */


                
        for (i = 0; i < 3; i++)                     
        {      
            P[i] = 0.; 
        }

        Q = QLeft * VolBbq;
        U = (FwdRateA * (1. - QLeft) - QLeft) * VolBbq; 

        aL00 =  1.;
        aL10 = -Q;          aL11 = -U;
        aL20 =  0.5*Q*Q;    aL21 =  1.5*U*Q;    aL22 =   U*U;

        Q = QRight * VolBbq;
        U = (FwdRateA * (1. - QRight) - QRight) * VolBbq; 

        aR00 =  1.;
        aR10 = -Q;          aR11 = -U;
        aR20 =  0.5*Q*Q;    aR21 =  1.5*U*Q;    aR22 =   U*U;

        
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {    
            StatePrL  = StatePr  + Hyb3_Node_Offset (2, i, 0, t, tree_data);
            DiscountL = Discount + Hyb3_Node_Offset (2, i, 0, t, tree_data);
                                                
            Zidx = Zt 
                 + (PreviousJump[0] + PreviousJump[1]) * i;
            Mid  = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid  = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);

            for (j = Bottom2[t][i]; j <= Mid ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aL10 + aL11 * x) * D0;
            
                P[2] += (aL20 + (aL21 + aL22 * x) * x) * D0;                  

            }

            for (j = Mid + 1; j <= Top2[t][i] ; j++)
            {
                
                x   = DiscountL[j];
                D0  = StatePrL[j] * x;
                P[0] += D0;

                P[1] += (aR10 + aR11 * x) * D0;
            
                P[2] += (aR20 + (aR21 + aR22 * x) * x) * D0;                  

            }
        }  /* for i */
        

        P[0] -= ZeroCoupon[t+1];

        /* Newton-Raphson to find the root of the polynomial. */
        /* We use the previous value as a first guess.        */

        DelZt = NR_Poly (0., P, 2);


        /* Check power series expansion was valid */ 
        if (fabs(DelZt * VolBbq) > 1.)
        {
            DR_Error ("Hyb3_Drift_2D: change in drift is too large");
            goto FREE_MEM_AND_RETURN;
        }
        


        Zt += DelZt;
        ZCenter[t] = Zt;


        /*   Recalculate the grid using the middle node  */
        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {
            /* LEFT part of distribution */

            Zidx  = Zt 
                  + (PreviousJump[0] + PreviousJump[1]) * i;
            Mid   = (int) ceil(-Zidx / PreviousJump[2]) - 1;
            Mid   = MIN ( MAX (Mid, Bottom2[t][i] - 1), Top2[t][i]);
            Zidx += (PreviousJump[2]) * Bottom2[t][i];

            DiscountL = Discount + Hyb3_Node_Offset (2, i, 0, t, tree_data);
    
            if (fabs(QLeft) > QCUTOFF)
            {
                Grid = MLeft * exp (QLeft * VolBbq * Zidx);
            
                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid *= RateJumpLeft2;
                }
            }
            else
            {
                Grid = MLeft * VolBbq * Zidx;

                for (j = Bottom2[t][i]; j <= Mid; j++)
                {
                    DiscountL[j] = 1. / (SLeft + Grid);
                
                    Grid += RateJumpLeft2;
                }
            }

            /* Right part of distribution */

            Zidx = Zt 
                 + (PreviousJump[0] + PreviousJump[1]) * i
                 + (PreviousJump[2]) * (Mid + 1);

            if (fabs(QRight) > QCUTOFF)
            {
                Grid = MRight * exp (QRight * VolBbq * Zidx);
            
                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);
                
                    Grid *= RateJumpRight2;
                }
            }
            else
            {
                Grid = MRight * VolBbq * Zidx;

                for (j = Mid + 1; j <= Top2[t][i]; j++)
                {
                    DiscountL[j] = 1. / (SRight + Grid);

                    Grid += RateJumpRight2;
                }
            }

        }  /* for i */

        /*
         *  Compute state prices for next period 
         */

        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {                                   
                StatePr1L[j] = 0.;
            }
        }  /* for i */
                
                        
        lMax = OutTop1[t+1] - 1;
        lMin = OutBottom1[t+1] + 1;

        if (lMin > lMax)
        {
            DR_Error ("Hyb3_Drift_2D: problem in building the tree (lMin > lMax)!");
            goto FREE_MEM_AND_RETURN;
        }


        for (i = Bottom1[t]; i <= Top1[t]; i++)
        {                        
            Pi = (d[0] - Beta1) * i;                
            Qi = (d[1] - Beta2) * i - Pi * Jump[1] / Jump[2];

            l = NEAR_INT (Pi);
            l = MIN (MAX (lMin - i , l), lMax - i);
                
            Pi -= l;
                                
            pu = .5 * (JumpCoeff + Pi + Pi * Pi);
            pd = pu - Pi;
            p0 = 1. - pu - pd;
                                
                                
            mMax =            OutTop2[t+1][i+l-1];
            mMax = MIN (mMax, OutTop2[t+1][i+l  ]);
            mMax = MIN (mMax, OutTop2[t+1][i+l+1]) - 1;

            mMin =            OutBottom2[t+1][i+l-1];
            mMin = MAX (mMin, OutBottom2[t+1][i+l  ]);
            mMin = MAX (mMin, OutBottom2[t+1][i+l+1]) + 1;

            if (mMin > mMax)
            {
                DR_Error ("Hyb3_Drift_2D: problem in building the tree (mMin > mMax)!");
                goto FREE_MEM_AND_RETURN;
            }


            StatePrL  = StatePr  + Hyb3_Node_Offset (2, i, 0, t, tree_data);
            DiscountL = Discount + Hyb3_Node_Offset (2, i, 0, t, tree_data);
                                        
            StatePr1L = StatePr1 + Hyb3_Node_Offset (2, i+1+l, 0, t+1, tree_data);
            StatePr2L = StatePr1 + Hyb3_Node_Offset (2, i  +l, 0, t+1, tree_data);
            StatePr3L = StatePr1 + Hyb3_Node_Offset (2, i-1+l, 0, t+1, tree_data);

            
            for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
            {
                Qij = Qi + (d[2] - Beta3) * j;

                m = NEAR_INT (Qij);
                m = MIN (MAX (mMin - j, m), mMax - j);
                
                Qij -= m;
                                
                qu = .5 * (JumpCoeff + Qij + Qij * Qij);
                qd = qu - Qij;
                q0 = 1. - qu - qd;

                        
                j1 = j + 1 + m;                 
                j2 = j     + m;
                j3 = j - 1 + m;

                x = StatePrL[j] * DiscountL[j];

                StatePr1L[j1] += pu * qu * x;
                StatePr1L[j2] += pu * q0 * x;
                StatePr1L[j3] += pu * qd * x;
                StatePr2L[j1] += p0 * qu * x;
                StatePr2L[j2] += p0 * q0 * x;
                StatePr2L[j3] += p0 * qd * x;
                StatePr3L[j1] += pd * qu * x;
                StatePr3L[j2] += pd * q0 * x;
                StatePr3L[j3] += pd * qd * x;
            }
        }  /* for i */

                           
        p = 0.;

        /*
         *   Check that state prices add up to 1
         */

        for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
        {
            StatePr1L = StatePr1 + Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
            
            for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
            {
                p += StatePr1L[j];
            }
        }  /* for i */

        if (fabs (p - ZeroCoupon[t+1]) > 0.0001)
        {
            DR_Error ("Hyb3_Drift_2D: state prices don't add up to 1: "
                        "increase sigma or Ppy!");
            goto FREE_MEM_AND_RETURN;
        }


        /* Express DEV tool: store state prices if necessary   */
        if ((EDevOn) && (EDevIdx<tree_data->NbEDevDates))
        {
            if (tree_data->TPDate[t+1] == tree_data->EDevDate[EDevIdx])
            {

                /* Store corresponding state prices */
                for (i = Bottom1[t+1]; i <= Top1[t+1]; i++)
                {
                    offset =  Hyb3_Node_Offset(2, i, 0, t+1, tree_data);
                    StatePr1L = StatePr1 + offset;
                    EDevPriceL = tree_data->EDevStPrice[EDevIdx] + offset;
            
                    for (j = Bottom2[t+1][i]; j <= Top2[t+1][i]; j++)
                    {
                        EDevPriceL[j] = StatePr1L[j];
                    }
                }  /* for i */

                /* Increment the EDev date index */
                EDevIdx++;
            }
        }

    
        
        /* Permute variables in preparation for next time step */          
        TempPointer = StatePr;
        StatePr     = StatePr1;
        StatePr1    = TempPointer;


    }  /* for t */ 
    
    
    /* Final check on number of EDev dates processed */
    if (EDevOn)
    {
        if (EDevIdx != tree_data->NbEDevDates)
        {
            DR_Error("Express DEV:  Error in processing  the DEV dates.\n"
                     "Please ensure that dates are ordered, that there\n"
                     "is no repetition,that value date is not included\n"
                     "and that the express DEV dates are also critical\n" 
                     "dates .\n");
            goto FREE_MEM_AND_RETURN;
        }
    }
    
                                        

                                                
    status = SUCCESS;

    FREE_MEM_AND_RETURN:

    Hyb3_Free_Slice (StatePr,  tree_data,2);
    Hyb3_Free_Slice (Discount, tree_data,2);
    Hyb3_Free_Slice (StatePr1, tree_data,2);
        
    return (status);

}  /* Hyb3_Find_Drift_2D*/




/*****  Hyb3_CUP_Drift  ********************************************************/
/**
 *  Adjust IR drift for currency protection, for several factor IR
 *
 */
int    Hyb3_CUP_Drift(
                 double    *Drift[],    /* (O) Cups drift for each fact.     */
                 double    **Aweight,   /* (I) Aweight at each time step     */
                 int       NbFactor,    /* (I) Number of factor              */
                 double    *SpotFxVol,  /* (I) Inst fx vol at each time step */
                 double    *Rho[],      /* (I) Corr. of fx with each factor  */
                 double    *Length,     /* (I) Length of each time step      */
                 int        NbTP)       /* (I) Total number of time points   */
 {

    int     t,i,j;
    int     status = FAILURE;    /** Error status = FAILURE initially */
    

    
    
    for (t = 0; t <= NbTP; t++)
    {
        for (i=0; i< NbFactor; i++)
        {
            Drift[i][t] = SQUARE(Aweight[i*(i+1)/2][t]);
            for (j=1; j<=i;j++)
            {
                Drift[i][t]+= SQUARE(Aweight[i*(i+1)/2+j][t]);
            }
            Drift[i][t] = sqrt(Drift[i][t]); /*Now, Drift[i][t] contains the 
                                              sptvol of fact. i at time step t */

            Drift[i][t] *= - Rho[i][t] * SpotFxVol[t] * Length[t];
        } /* For i */
    } /* for t */

    status = SUCCESS;

    return (status);

}  /* Hyb3_CUP_Drift */

/*****  FX_CUPS_FXSmile ***************************************************/
/*
 *  calculate the FX smile adjusted CUPS drift 
 *  currently, the function relies heavily on the 
 *  time consuming FXMomentMatching
 */
int     Hyb3_CUPSAdjust_FXSmile(HYB3_TREE_DATA  *tree_data,     /* (O) Tree data             */
                                MKTVOL_DATA     *mktvol_data,   /* (I) Market vol data       */
                                FX_DATA         *fx_data)       /* (I) FX Data for flat smile test */
{
    int status = FAILURE;
    int isLN = TRUE, i;
    double *scale1st = NULL, *smoothScale = NULL;
    double *newDrift       = NULL;
    double *adjZeroCoupons = NULL;
    double *adjust         = NULL;
    double *timeLine       = NULL ; 

    
    /* Local var. for 1D drift cups */
    double *DriftCUPSL = tree_data->DriftCUPS[0];

    /* FILE * fp = NULL;
    char fName[200]; 
    static int itCnt = 0; */


    /* allocate vector for a first moment matching vector */
    scale1st    = (double *)DR_Array(DOUBLE, -1, tree_data->NbTP+1);
    smoothScale = (double *)DR_Array(DOUBLE, -1, tree_data->NbTP+1);

    newDrift   = (double *)DR_Array(DOUBLE, 0, tree_data->NbTP);
    adjZeroCoupons = (double *)DR_Array(DOUBLE, 0, tree_data->NbTP + 1);
    adjust = (double *)DR_Array(DOUBLE, 0, tree_data->NbTP);
    timeLine = (double *)DR_Array(DOUBLE, 0, tree_data->NbTP);
    if ((scale1st       == NULL)||
        (smoothScale    == NULL)||
        (newDrift       == NULL)||
        (adjZeroCoupons == NULL)||
        (adjust         == NULL)||
        (timeLine       == NULL)) 
    {

        DR_Error("Hyb3_CUPSAdjust_FXSmile: Memory allocation failure");
        goto RETURN;
    }



    /* test if the market contains smile information                  */
    /* if not there is no need for any adjustment (LN Cups is correct) */
    for ( i = 0 ; i < fx_data->NbSmilePt; ++i)
    {
        if ( (fabs(fx_data->a1[i]) > TINY ) || (fabs(fx_data->a2[i]) > TINY ) 
            || (fabs(fx_data->a3[i]) > TINY ) ) 
        {
            isLN = FALSE; 
        }
    }
    if (isLN) 
    {
        status = SUCCESS;
        goto RETURN;
    }


    /* calculate the adjustment factor first */
    if ( Hyb3_FX_MomentMatching( mktvol_data, tree_data, scale1st ) == FAILURE)
    {
        goto RETURN;
    }


    /* smooth the scale function to avoid big swings in CUPS drift */
    for (i = 0 ; i <= tree_data->NbTP; ++i)
    {
        if (i == 0)
        {
            timeLine[i] = 0.0;
        }
        else
        {
            timeLine[i] = timeLine[i-1] + tree_data->Length[i-1];
        }
    }
    if (Hyb3_SmoothCurve(smoothScale, timeLine,scale1st, 
                tree_data->NbTP, tree_data->SmileIndex, 50) != SUCCESS )
    {
        DR_Error("Failed to smooth the scaling curve!");
        goto RETURN;
    }

    
    /* find the adjustment for the CUPS drift */
    for (i = 0 ; i <= tree_data->NbTP+1 ; ++i )
    {
        /* only applicable for foreign IR under FX */
        adjZeroCoupons[i] = tree_data->ZeroCoupon[0][tree_data->CvDiff[0]][i]   
                    *smoothScale[ MIN(i,tree_data->NbTP) ];
    }

    /* invoke the drift-calculation for the modified forwards */
    /* Get one factor IR grid drift */
    if (Hyb3_Find_Drift_1D(newDrift,
                            adjZeroCoupons,
                            tree_data->FwdRate[0][tree_data->CvDiff[0]],
                            mktvol_data[0].QLeft,
                            mktvol_data[0].QRight,
                            mktvol_data[0].FwdShift,
                            mktvol_data[0].Beta[0],
                            tree_data->SpotVol[0],
                            tree_data->Top1,
                            tree_data->Bottom1,
                            tree_data->OutTop1,   
                            tree_data->OutBottom1,
                            tree_data->NbTP,
                            tree_data->Length,
                            tree_data->LengthJ,
                            tree_data) == FAILURE)
    {
        goto RETURN;
    }
    

    /* debug information (to be deleted later ) */
  /*  sprintf(fName,"CUPSDrift_%d.dat",itCnt++);
    fp = fopen(fName, "w"); 
    printf("#no.  Length   newDrift  oldDrift  scale1st smoothed  adjZeroCoupon newDrift IRCenter\n"); */
    for (i = 0 ; i <= tree_data->NbTP ; ++i )
    {
        double dt = tree_data->Length[i];
        
        newDrift[i] -= tree_data->IrZCenter[0][i];               
        adjust[i] = newDrift[i] - ( (i > 0) ? (newDrift[i-1]*exp( -mktvol_data[0].Beta[0]*dt ) ) : 0.0 );
    }


    /* smooth the drift if necessary (stability criterion, see CUPS drift document for explanation) */
    /* increasing number of standard deviation in the tree helps to avoid this smoothing */
    for (i = 0 ; i < tree_data->NbTP ; ++i ) 
    {
        double driftLN = DriftCUPSL[i];

        /* cap and floor maximal variations on the drift (to avoid instabilities) */
        /* this value is a bit arbitrary and needs more testing */
        if (fabs(adjust[i]) > ADJSMOOTH_FAC*fabs(driftLN))
        {
            /* average over consecutive adjustment for stability reason */
            adjust[i+1] = (adjust[i] = 0.5*(adjust[i+1]+adjust[i]) );


            /* cut the adjustment to be no more than CUTOFF of the original drift -> stability */
            adjust[i]   = MAX( -CUPSADJ_CUTOFF*fabs(driftLN),
                                MIN(CUPSADJ_CUTOFF*fabs(driftLN), adjust[i] ));

            adjust[i+1] = MAX( -CUPSADJ_CUTOFF*fabs(DriftCUPSL[i+1]),
                                MIN(CUPSADJ_CUTOFF*fabs(DriftCUPSL[i+1]), adjust[i+1] ));

            /* the next lines are kept for testing purpose (to be deleted later) */
            
            /* fprintf(fp,"%d  %6.5f  %12.10f %12.10f  %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n"
                 ,i,tree_data->Length[i],DriftCUPSL[i]+adjust[i],DriftCUPSL[i],
                  scale1st[ MIN(i,tree_data->NbTP) ],smoothScale[ MIN(i,tree_data->NbTP) ],adjZeroCoupons[i],newDrift[i],tree_data->IrZCenter[0][i],
                 DriftCUPSL[i]+0.5*(adjust[i+1]+adjust[i]));
 
             fprintf(fp,"%d  %6.5f  %12.10f %12.10f  %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n"
                 ,i+1,tree_data->Length[i+1],DriftCUPSL[i+1]+adjust[i+1],DriftCUPSL[i+1],
                 scale1st[ MIN(i+1,tree_data->NbTP) ],smoothScale[ MIN(i+1,tree_data->NbTP) ],adjZeroCoupons[i+1],newDrift[i+1],tree_data->IrZCenter[0][i+1],
                 DriftCUPSL[i+1]+0.5*(adjust[i+1]+adjust[i]) ); */
                

            DriftCUPSL[i]   += adjust[i]; 
            DriftCUPSL[i+1] += adjust[i+1]; 

            ++i;
        }
        else
        {
            
       /*       fprintf(fp, "%d  %6.5f  %12.10f %12.10f  %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n"
                 ,i,tree_data->Length[i],DriftCUPSL[i]+adjust[i],DriftCUPSL[i],
                 scale1st[ MIN(i,tree_data->NbTP) ],smoothScale[ MIN(i,tree_data->NbTP) ],adjZeroCoupons[i],newDrift[i],tree_data->IrZCenter[0][i],
                 DriftCUPSL[i]+adjust[i]); */

            DriftCUPSL[i] += adjust[i]; 
        }

    }

    status = SUCCESS;
RETURN:
    /* if (fp != NULL) fclose(fp);  */

    Free_DR_Array(scale1st, DOUBLE, -1, tree_data->NbTP+1);
    Free_DR_Array(smoothScale, DOUBLE, -1, tree_data->NbTP+1);
    Free_DR_Array(newDrift, DOUBLE, 0, tree_data->NbTP);
    Free_DR_Array(adjZeroCoupons, DOUBLE, 0, tree_data->NbTP+1);
    Free_DR_Array(adjust, DOUBLE, 0, tree_data->NbTP);
    Free_DR_Array(timeLine, DOUBLE, 0, tree_data->NbTP);
    
    return (status);
}


/*****  FX_MomentMatching  ***************************************************/
/*
 *  Adjust first and second moment of the FX-spot to fit "externally"  
 *  specified FX fwd and ATM volatility
 *  the code is an adopted version of the FXseries
 */
int     Hyb3_FX_MomentMatching(MKTVOL_DATA     *mktvol_data,   /* (I) Market vol data       */
                          HYB3_TREE_DATA       *tree_data ,    /* (I) Tree data             */
                          double               *scale1st)      /* (O) moment matching scaling */
{

    HYB3_DEV_DATA     dev_data;

    /* Slices and pointers to slices */
    TSLICE       StatePr   = NULL;  /* Variable containing state prices at t   */  
    TSLICE       StatePr1 = NULL;   /* Variable containing state prices at t+1 */  

    double  *StatePrL  = NULL;      /* pointers state price slices      */
    double  *TempPointer = NULL;

    double  *SpotFXL;               /* pointer to spot fx               */
    
    int
       *Top1,     *Bottom1,         /* Limits of the tree (1rst dimension)  */
      **Top2,    **Bottom2,         /* Limits of the tree (2nd dimension)   */
     ***Top3,   ***Bottom3;         /* Limits of the tree (3rd dimension)   */


    /* local variables for easy referencing */
    int      i, j, k,
             IdxBuffer;

    int     DCurveD,           /* domestic discount curve */
            ExerciseFlag,      /* TRUE if current date is an exercise date */

            T,           /* Total number of period in the hybrids tree*/
            t,          /* Current time period                       */
            status = FAILURE; /* Error status = FAILURE initially          */


    /* FX Moment matching variables */
    double DetermFwdFX;        /* local (temporary) variable for FX-FWD     */
    double DetermDZero;        /* domestic zero coupon bond					*/
    double DetermFZero;        /* domestic zero coupon bond					*/

    double TreeFwdFX;          /* local (temporary) variable for tree FX-FWD    */
    double TreeDomZero;        /* local (temporary) variable for tree zero bond */


    /* initialise dev structure */
    Hyb3_Dev_Init(&dev_data);

    /* Total size of tree timeline */     
    T   = tree_data->NbTP;


    /* Assigment of domestic discount curve */
    DCurveD   = tree_data->CvDisc[1];

    IdxBuffer         = -1;    

    StatePr  = Hyb3_Alloc_Slice (tree_data,3);
    StatePr1 = Hyb3_Alloc_Slice (tree_data,3);
    if ((StatePr  == NULL) || (StatePr1 == NULL))
    {
        DR_Error ("could not allocate memory!");
        goto RETURN;
    }
       
    if (Hyb3_Dev_Alloc(&dev_data, tree_data) == FAILURE)
    {
        goto RETURN;
    }

    StatePrL = StatePr + Hyb3_Node_Offset(3, 0, 0, 0, tree_data);
    StatePrL[0] = 1.0;
       
    Top1    = tree_data->Top1;	        
    Top2    = tree_data->Top2;
    Top3    = tree_data->Top3;	    
    Bottom1 = tree_data->Bottom1;            
    Bottom2 = tree_data->Bottom2;    
    Bottom3 = tree_data->Bottom3;



    /* to avoid checks in lattice.s, 
       we use our own check here for forward going tree */
    tree_data->CalcCheckSlices = FALSE;
    scale1st[0] = 1.0; /* initialise the first scaling point */

    for (t = 0; t <= T; t++)
    {     
        /*needs checking what flag is general for this purpose..*/
        ExerciseFlag = TRUE; /* tree_data->TPType[0][t]; */
    
        /* updates state prices (and fills the FX grid) */
        Hyb3_UpdateStatePrices(t,mktvol_data,tree_data,&dev_data,StatePr,StatePr1);
               

        if (( t > 0) && ExerciseFlag) /* only future time points need to be rescaled */
        {
            /* scale the average of the FX (first moment) */
            DetermFwdFX = tree_data->FwdFx[t];
            DetermFZero = tree_data->ZeroCoupon[0][tree_data->CvDisc[0]][t];
            DetermDZero = tree_data->ZeroCoupon[1][tree_data->CvDisc[1]][t];
            
            
            TreeFwdFX   = 0.0;
            TreeDomZero = 0.0;
            for (i = Bottom1[t]; i <= Top1[t]; i++)
            {
                for (j = Bottom2[t][i]; j <= Top2[t][i]; j++)
                {
                    StatePrL  = StatePr + Hyb3_Node_Offset (3, i, j, t, tree_data);
                    SpotFXL   = dev_data.FxSpot + Hyb3_Node_Offset (3, i, j, t, tree_data);
                    
                    for (k = Bottom3[t][i][j]; k <= Top3[t][i][j]; k++)
                    {
                        TreeFwdFX   += StatePrL[k] *  SpotFXL[k];
                        TreeDomZero += StatePrL[k];
                    }
                }  /* for j */
            }  /* for i */
            /* do the actual scaling */
            if ( IS_EQUAL( TreeFwdFX, 0.0) )
            {
                DR_Error ("TreeFwdFX is zero and would lead to division by zero (internal error)!");
                goto RETURN;
            }
            scale1st[t] = DetermFwdFX * DetermDZero/TreeFwdFX;
            
        } /* if (t>0) */

        TempPointer = StatePr;
        StatePr     = StatePr1;
        StatePr1    = TempPointer;

    }  /* End of loop forwards in tree */

    status = SUCCESS;

    RETURN:

    Hyb3_Free_Slice (StatePr,  tree_data, 3);
    Hyb3_Free_Slice (StatePr1, tree_data, 3);

    Hyb3_Dev_Free    (&dev_data, tree_data);

                                    
    return (status);

} /* FX_Momentmatching */
