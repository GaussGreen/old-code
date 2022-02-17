/****************************************************************************/
/*      Calculation of forward stock price, forward equity volatility.       */
/****************************************************************************/
/*      STOCK.c                                                             */
/****************************************************************************/


/*
$Header$
*/



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cupslib.h"

/*****  Hyb3_Get_ZeroBorrow  *****************************************************/
/* Borrow Rate is:
 *  - flat before first benchmark,
 *  - linear between benchmarks,
 *  - extrapolated after the last benchmark
 * Input Rates are assumed Act/365, always
 */
int Hyb3_Get_ZeroBorrow(
        double *ZeroRate,       /* (O) Zero Rate                            */
        double *ZeroCost,       /* (O) Borrowing cost (factor)              */
        long    ValueDate,      /* (I) Value date                           */
        long    EndDate,        /* (I) Date for Zero Borrow                 */
        int     NbBorrow,       /* (I) Number of points in curve            */
        long   *BorrowDate,     /* (I) Borrowing curve dates                */
        double *Borrow)         /* (I) Borrowing cost zero rate             */
{
    double
        Days;                   /* Maturity of the zero in days             */
    int
        k,
        status = FAILURE;

    if (NbBorrow <= 0)
    {
        DR_Error("Borrowing Cost of Equity Missing.\n");
        goto RETURN;
    }
    if ((NbBorrow == 1) ||
        (EndDate <= BorrowDate[0]))
    {
        *ZeroRate = Borrow[0];
    }
    else
    {
        k = NbBorrow - 2;
        while ((EndDate < BorrowDate[k]) && (k > 0))   /* k is the index of the zero maturing just before the maturity. 	*/
            k--;

        dlinterp(EndDate,
                 ZeroRate,
                 BorrowDate[k],
                 BorrowDate[k+1],
                 Borrow[k],
                 Borrow[k+1]);
    }

    Days = Daysact (ValueDate, EndDate);
    *ZeroCost = exp(- *ZeroRate * Days / 365.);

    status = SUCCESS;
RETURN:
    return (status);
}

/*****  Hyb3_Forward_Stock  ******************************************************/
/*
*       Interpolate forward stock price at each node point i.e. every
*       Length[i] interval.There are 4 different ways to do the inter-
*       polation to cope with commodities and equity.
*       In addition, this function positions the centre of the eq dim
*       on the equity forward.  If running on 3F mode  (i.e. no  Q(t)
*       smile available) then the  centre is later corrected to be on
*       the centre of the distribution.
*
*       See also function Hyb3_Eq_Shift_Centre().
*
*/
int	Hyb3_Forward_Stock(
        double   *FwdEq,          /* (O) Fwd stock price at eachtime point   */
        long     *NodeSettleDate, /* (O) Stock settlement date for time point*/
        double   *NodeSettleTime, /* (O) Time from time point to stock settl */
        long     Today,           /* (I) Date of equity spot quote           */
        double   Spot,            /* (I) Spot stock price                    */
        int      NbDiv,           /* (I) Number of dividend / forward prices */
        double   *Div,            /* (I) Dividends / Forward prices          */
        long     *DivDate,        /* (I) Dividends / forwards dates          */
        char     *DivType,        /* (I) Dividend / forward type             */
        int      NbSettle,        /* (I) Nb of settlement dates <0 if rolling*/
        long     *LastTrading,    /* (I) Last trading date of settl period   */
        long     *SettleDate,     /* (I) Corresponding settlement date       */
        char     SettleType,      /* (I) Settlement type:'F'ixed or 'R'olling*/
        int      NbBorrow,        /* (I) Number of points in borrowing curve */
        long     *BorrowDate,     /* (I) Borrowing curve dates               */
        double   *Borrow,         /* (I) Borrowing cost zero rate            */

        T_CURVE const* zcrv,      /* (I) Zero curve                        */

        long     *TPDate,         /* (I) Date of each node                   */
        double   *Length,         /* (I) Length of each time step            */
        double   *ZeroCoupon,     /* (I) Zero coupon price at time point     */
        int      NbTP)            /* (I) Total number of time points         */
{
        int
                NbTPAux;
        long
                *TPDateAux = NULL;
        double
                *FwdEqAux = NULL;
        double
                x, y, z, b, dt;
        int
                i1,        /* Node index of the previous volatility point */
                i2,        /* Node index of the current volatility point  */
                i3,        /* Node index of the next volatility point     */
                i, k, v,
                status = FAILURE; /* Error status = FAILURE initially     */
        char
                ErrorMsg[MAXBUFF]; /* Error message to be sent to DR_Error*/

        /* 	Calculate settlement date and time for each node. */
        if (SettleType == 'R')     /* Rolling settlement */
        {
            for (i = 0; i <= NbTP; i++)
            {
                /* Take into account week-end but not holidays */
                NodeSettleDate[i] = Nxtwkday(TPDate[i],(long) NbSettle);
                NodeSettleTime[i] = Daysact(TPDate[i],NodeSettleDate[i]) /365.;

            }  /* for i */
        }
        else  /* Fixed settlement */
        {
            i1 = i2 = 0;

            for (k = 0; k < NbSettle; k++)
            {
                /* Look for node index i2 falling after curr last trading dt */
                while ((TPDate[i2] <= LastTrading[k]) && (i2 <= NbTP))	
                        i2++;

                for (i = i1; i < i2; i++)
                {
                    NodeSettleDate[i] = SettleDate[k];
                    NodeSettleTime[i] = Daysact(TPDate[i],NodeSettleDate[i])/365.;

                }  /* for i */

                i1 = i2;

            }  /* for k */

            /* If there are some nodes left after the last settlement */
            /* date just assume it settles immediately.               */
            for (i = i1; i <= NbTP; i++)
            {
                NodeSettleDate[i] = TPDate[i];
                NodeSettleTime[i] = 0.;
            }
        }

        /* construct TPDateAux = {Today} U {TPDate} U
                                 {DivDate > Today, <= TPDate[NbTP]} */

        TPDateAux = (long   *) DR_Array(LONG,   0, 1 + NbTP+1 + NbDiv - 1);
        FwdEqAux  = (double *) DR_Array(DOUBLE, 0, 1 + NbTP+1 + NbDiv - 1);

        if ((TPDateAux == NULL) ||
            (FwdEqAux  == NULL))
        {
            DR_Error("Could not allocate memory for aux timeline (hyb3_stock.c).\n");
            goto RETURN;
        }

        if (TPDate[0] < Today)
        {
            DR_Error("Today must not be after ValueDate.\n");
            goto RETURN;
        }

        TPDateAux[0] = Today;

        k = 0;

        while ((k < NbDiv) &&
               (DivDate[k] <= Today))
        {
            ++k;
        }

        i = 0; /* TP    */
        v = 1; /* TPAux */

        while (i <= NbTP)
        {
            while ((k < NbDiv) &&
                   (DivDate[k] <= TPDate[i]))
            {
                if (v == 1 + NbTP+1 + NbDiv)
                {
                    DR_Error("Error: not enough memory for TP + Div.\n");
                    goto RETURN;
                }

                TPDateAux[v] = DivDate[k];
                ++k;
                ++v;
            }

            if (TPDateAux[v-1] != TPDate[i])
            {
                if (v == 1 + NbTP+1 + NbDiv)
                {
                    DR_Error("Error: not enough memory for TP + Div.\n");
                    goto RETURN;
                }

                TPDateAux[v] = TPDate[i];
                ++v;
            }

            ++i;
        }

        NbTPAux = v;

        /* Verify that we really have a sorted union.          */
        /* It's still quicker than appending and then sorting. */

        for (v=1; v<NbTPAux; ++v)
        {
            if (TPDateAux[v] <= TPDateAux[v-1])
            {
                DR_Error("Coding error: aux timeline not sorted "
                         "(Hyb3_Forward_Stock).\n");
                goto RETURN;
            }
        }

        v = 0;
        
        for (i=0; i<=NbTP; ++i)
        {
            while ((v < NbTPAux) &&
                   (TPDateAux[v] != TPDate[i]))
            {
                ++v;
            }

            if (v == NbTPAux)
            {
                DR_Error("Coding error: aux timeline does not contain "
                         "original timeline. (Hyb3_Forward_Stock).\n");
                goto RETURN;
            }
        }

        v = 0;
        
        for (k=0; k<NbDiv; ++k)
        {
            if (DivDate[k] <= Today)
            {
                continue;
            }

            if (DivDate[k] > TPDate[NbTP])
            {
                break;
            }

            while ((v < NbTPAux) &&
                   (TPDateAux[v] != DivDate[k]))
            {
                ++v;
            }

            if (v == NbTPAux)
            {
                DR_Error("Coding error: aux timeline does not contain "
                         "all divs > Today. (Hyb3_Forward_Stock).\n");
                goto RETURN;
            }
        }

        /* Calculate the (aux) deflated forward curve. */

        z = GetZeroPrice(NodeSettleDate[0], zcrv);

        if (z < TINY)
        {
            goto RETURN;
        }

        FwdEqAux[0] = Spot * z;

        z = GetZeroPrice(TPDate[0], zcrv);

        if (z < TINY)
        {
            goto RETURN;
        }

        i1 = i2 = 0;

        for (k = 0; k < NbDiv; k++)
        {
            /* If current type is continuous div then the following are too */
            if (DivType[k] == 'C')
                break;

            /* Look for the node index i2 corresponding to the current div */
            /* It exists as dividend dates are critical dates.             */
            while ((TPDateAux[i2] < DivDate[k]) && (i2 < NbTPAux-1))  
                i2++;

            if (i2 == i1)
                continue;

            /* Forward stock is flat between two dividend dates */
            for (i = i1 + 1; i <= i2; i++)
            {
                FwdEqAux[i] = FwdEqAux[i1];
            }

            /* Need this test for i2 = NbTP and TPDate[i2]<DivDate[k] */
            if (TPDateAux[i2] == DivDate[k])
            {
                switch (DivType[k])
                {
                    case 'D' : /* FIXED $ DIV */
                    {
                        if (Hyb3_Get_ZeroBorrow(&x,
                                           &b,
                                           zcrv->ValueDate,
                                           DivDate[k],
                                           NbBorrow,
                                           BorrowDate,
                                           Borrow) == FAILURE)  { goto RETURN; }

                        if (GetZeroPriceRate(&x,
                                      &y,
                                      DivDate[k],
                                      zcrv) == FAILURE)  { goto RETURN; }

                        y /= z;

                        /* Subtract PV of $ div from previous ex-div price */
                        FwdEqAux[i2] -= Div[k] * y / b;

                        break;
                    }
                    case 'Y' : /* FIXED DIV YIELD */
                    {
                        /* 'Y' means a proportional drop of the Forward,
                         * irrespective of Settlement, irrelevant to that
                         * degree of precision */
                        FwdEqAux[i2] *= (1.0 - Div[k]);

                        break;
                    }
                    default:
                    {
                        DR_Error("Internal Error, wrong dividend type.\n");
                        goto RETURN;
                    }
                }  /* switch */
            }

            if (FwdEqAux[i2] < 0.0)
            {
                sprintf(ErrorMsg, "Stock Forward goes negative on the %ld.\n",
                        DivDate[k]);
                DR_Error (ErrorMsg);
                goto RETURN;
            }
        
            i1 = i2;

        }  /* for k */

         /* There are some continuous dividends left */
        if (k < NbDiv)
        {
            while ((TPDateAux[i2] < DivDate[k]) && (i2 < NbTPAux-1))
                    i2++;

            /* Flat forward from last discrete dividend date */
            /* until current continuous dividend date.	     */
            for (i = i1 + 1; i <= i2; i++)
            {
                FwdEqAux[i] = FwdEqAux[i1];
            }  /* for i */
        }

        i3 = 0;

        for (; k < NbDiv; k++)
        {
            /* Only point left: we fill nodes until last */
            if (k == NbDiv - 1)
            {
                i3 = NbTPAux-1;
            }
            else /* Some continuous dividends left afterwards */
            {
                while((TPDateAux[i3] < DivDate[k+1]) && (i3 < NbTPAux-1))
                    i3++;
            }


            for (i = i2 + 1; i <= i3; i++)
            {
                dt = (double)Daysact(TPDateAux[i-1], TPDateAux[i]) / 365;

                /* Continuous dividends are not input in % */
                FwdEqAux[i] = FwdEqAux[i-1] * exp (- Div[k] * dt);
            }

            i2 = i3;

        }  /* while */

        /* Flat forward stock after the last discrete div or forward price */
        /* Only ever get there if there are no continuous divs             */
        for (i = i2 + 1; i < NbTPAux; i++)  
        {
            FwdEqAux[i] = FwdEqAux[i2];
        }

        /* Take the subsample for the tree timeline */
        /* and forward the value.                   */

        v=0;

        for (i=0; i<=NbTP; ++i)
        {
            while ((v < NbTPAux) &&
                   (TPDateAux[v] != TPDate[i]))
            {
                ++v;
            }

            if (v == NbTPAux)
            {
                DR_Error("Coding error in div schedule setup.\n");
                goto RETURN;
            }

            FwdEq[i] = FwdEqAux[v];

            /* Forward-Value to settlement date */
            if (GetZeroPriceRate(&x,             /* Get zero corresponding to spot */
                         &y,             /* settlement date                */
                         NodeSettleDate[i],
                         zcrv) != SUCCESS)
            {
                goto RETURN;
            }
            /* Pay borrow up to the settlement date */
            if (Hyb3_Get_ZeroBorrow(&x,
                               &b,
                               zcrv->ValueDate,
                               NodeSettleDate[i],
                               NbBorrow,
                               BorrowDate,
                               Borrow) == FAILURE)  { goto RETURN; }
            FwdEq[i] *= b / y;

        }  /* for i */

        status = SUCCESS;

        RETURN:

        Free_DR_Array(TPDateAux, LONG,   0, 1 + NbTP+1 + NbDiv - 1);
        Free_DR_Array(FwdEqAux,  DOUBLE, 0, 1 + NbTP+1 + NbDiv - 1);

        return (status);

}  /* Hyb3_Forward_Stock */



/*****	Hyb3_Composite_Eq  *******************************************************/
/*
*	Adjust forward,  volatility,  correlations and tree center of composite 
*   equity. The adjustements are based on the transformation of the foreign
*   equity in a domestic asset by multiplying by the fx rate.
*/
int	Hyb3_Composite_Eq ( double	*FwdEq,          /* (I/O) Forward stock at each node */
                        double	*EqMidNode,      /* (I/O) Center of the equity tree  */
                        double	*SpotEqVol,      /* (I/O) Instantaneous equity vol   */
                        double **Rho,            /* (I/O) Instantaneous correlations */
                        long    *NodeSettleDate, /* (I) Stock settlement dates       */
                        T_CURVE const* fcrv,     /* Foreign zero curve               */

                        T_CURVE const* dcrv,     /* Domestic zero curve              */

                        double	*FwdFx,          /* (I) Forward fx                   */
                        double	*SpotFxVol,      /* (I) Instantaneous fx volatility  */
                        double  *Length,         /* (I) Length of each time step     */
                        int     NbTP)            /* (I) Total number of time steps   */
{



    double  FxFwdtoSettle;    /* Fx forward maturing on settlement date */
    double  VolD;             /* Foreign (i.e. composite) equity vol    */
    double  VarF;             /* Total variance of domestic equity      */
    double  VarD;             /* Total variance of foreign equity       */
    double  x, y;
    int     i;
    int     status = FAILURE; /* Error status = FAILURE initially  */
        

    VarF = VarD = 0;

    for (i = 0; i <= NbTP; i++)
    {
        /* Calculate fx forward maturing on settlement date. It will */
        /* be used to adjust the fwd stock (we cannot use FwdFx[i]   */
        /* as it matures on current date, not settlement date).      */

        FxFwdtoSettle = FwdFx[0]; /* i.e. spot fx */

        if (GetZeroPriceRate(&x,
                     &y,
                     NodeSettleDate[i],
                     fcrv) != SUCCESS)
        {
            goto RETURN;
        
        }

        FxFwdtoSettle *= y;

        if (GetZeroPriceRate(&x,
                     &y,
                     NodeSettleDate[i],
                     dcrv) != SUCCESS)
        {
            goto RETURN;
        
        }

        FxFwdtoSettle /= y;

        FwdEq[i] *= FxFwdtoSettle;


        /* Multiply the tree center by the fx fwd and change the lognormal */
        /* correction to take into account the new volatility.             */

        EqMidNode[i] *= FwdFx[i];
        EqMidNode[i] *= exp (-.5 * (VarD - VarF));


        /* Calculate composite volatility and new correlations */

        VolD = sqrt (SQUARE(SpotEqVol[i]) + SQUARE(SpotFxVol[i])
                             + 2. * Rho[5][i] * SpotEqVol[i] * SpotFxVol[i]);

        Rho[3][i] = (Rho[3][i] * SpotEqVol[i] + Rho[1][i] * SpotFxVol[i]) / VolD;
        Rho[4][i] = (Rho[4][i] * SpotEqVol[i] + Rho[2][i] * SpotFxVol[i]) / VolD;
        Rho[5][i] = (Rho[5][i] * SpotEqVol[i] +             SpotFxVol[i]) / VolD;


        VarF += SQUARE(SpotEqVol[i]) * Length[i];
        VarD += SQUARE(VolD) * Length[i];


        SpotEqVol[i] = VolD;

    }  /* for i */	
                                   

    SpotEqVol[-1] = SpotEqVol[0];

    Rho[3][-1] += Rho[2][0];
    Rho[4][-1] += Rho[1][0];
    Rho[5][-1] += Rho[2][0];



    status = SUCCESS;

    RETURN:

    return (status);

}  /* Hyb3_Composite_Eq */



/*****  Hyb3_CUPSEqMidNode  *******************************************************/
/*
 *	Adjust equity mid-node for currency protection.
 *
 */
int Hyb3_CUPSEqMidNode(double   *EqMidNode, /*(I/O) Mid node to be adjusted       */
                  double   *SpotEqVol, /*(I) Inst eq vol at each time step   */
                  double   *SpotFxVol, /*(I) Inst fx vol at each time step 	 */
                  double   *Rho,       /*(I) Correlation of fx with equity   */      
                  double   *Length,    /*(I) Length of each time step        */
                  int       NbTP)      /*(I) Total number of nodes           */
{

    
    double 
        integral = 0.0;
    int
        i,
        status = FAILURE;    /* Error status = FAILURE initially */


    for (i = 0; i <= NbTP; i++)
    {
        
        EqMidNode[i] += integral; 
        integral -= Rho[i] * SpotFxVol[i] * SpotEqVol[i] * Length[i];
 
    }  /* for i */


    status = SUCCESS;

    return (status);

}  /* End of Hyb3_CUPSEqMidNode() */



/*****    Hyb3_MultiFac_Spot_EqVol  ********************************************************/
/*
 *    Interpolate the correlation  curve at each node  in the tree. 
 *    Calculate  equity spot volatility at each node in the tree by
 *    bootstrapping input benchmark option volatility curve, taking
 *    into account interest rate volatility.
 */
 
int   Hyb3_MultiFac_Spot_EqVol 
    (double  *SpotEqVol,   /* (O) Spot equity volatility at each time step   */
     double  *FwdEq,       /* (I) Forward equity price at each time step     */
     /* Model IR related quantities */
     int      NbFac,       /* (I) Number of factors for the considered rate  */
     double  *Alpha,       /* (I) Value of alpha,either domestic or foreign  */
     double  *Beta,        /* (I) Mean reversion coefficient                 */
     /* IR rates characteristics */
     double  *SpotVol,     /* (I) Spot vol of short rate at each time step   */
     double  *FwdRate,     /* (I) Forward rate at each node in the tree      */
     double  *RhoIrEq,
     /*General tree quantities */
     long    *TPDate,      /* (I) Date of each node                          */
     long     NbTP,
   
     int      NbImpVol,    /* (I) Number of points in equity volatility curve*/
     long    *VolDate,     /* (I) Benchmark equity option expiration dates   */
     double  *EqVol,       /* (I) Benchmark equity option volatilities       */
     /*EQS*/   
     /*Additional equity smile parameters*/
     int      EqCutOffFlag,   /* (I) true if cut off allowed       */
     int      EqCutOffLast,   /* (I) true if cut off at last level */
     double   EqCutOffLevel,  /* (I) user-given or from cups.h     */
     double  *A_Eq,           /* (I) EQ smile param, [0,NbTP-2]    */
     double  *B_Eq,           /* (I) EQ smile param, [0,NbTP-2]    */
     double  *C_Eq            /* (I) EQ smile param, [0,NbTP-2]    */
     )
{

    double  Poly[3],
            EqVolPrev,
            epsilon, 
            tau_0, t_low, delta_t, slope, /* times used for interpolation of implied vols   */ 
            temp1, temp2,              
            Sigma = 0.0,                  /* y Spot Vol in current bucket                  */
            discriminant;                 /* discriminant of sec. order polynomial          */
   
    int    
            i = 0, 
            k,
            NbNoEqInt,
            NbWithEqInt,
            status = FAILURE;  /* Error status = FAILURE initially          */
    char
            ErrorMsg[MAXBUFF]; /* Error message to be sent to DR_Error      */
    double  *DomA;
    double  *NoEqInt;
    double  *WithEqInt;
    double  *DeltaTime;
    double  *Time;
    double  *a;
    double  *b;
    double  *c;
    double  *volA;
    double  *loc_vol_x;
    double  *loc_vol_xx;
    double  *SpotEqVolSoFar;
    double  *EqVolExt;
    double  *Perturb;
    double  *EqVolUnSmile;

    int      kVol,
             TPIntBefore,           /* integration time index of prev.Vol point */
             TPIntSoFar,
             TPFirstImpVol,         /* time index for first implied vol         */
             TPLastImpVol,          /* time index for last implied vol          */
             NextTP,                /* temp variable used for TPDate offset     */ 
             VolTooLow,             /* flag for spot vol criterion              */
             l,n,
             CurrentNbInt;          /* used to index integrals not involving y rate */
          
    long IsSmileE = TRUE;           /* Flag for equity smile                         */

/*Additional perliminary variables to match the multifac code */

   /* Quick checks first */
/**************************/
    if (NbTP < 2)
    {
        DR_Error ("Invalid number of time points in MultiFac_Spot_EqVol: "
                    "should be at least 2\n");
        return (status);
    }

    if (NbImpVol < 1)
    {
        DR_Error ("Invalid number of Vol points in MultiFac_Spot_EqVol2: "
                    "should be at least 1");
        return (status);
    }

    if( SpotVol          == NULL ||
        Alpha            == NULL ||
        FwdRate          == NULL ||
        VolDate          == NULL ||
        EqVol            == NULL || 
        Beta             == NULL ||
        TPDate           == NULL)
    {   
        DR_Error("Invalid pointer inputs to MultiFac_Spot_EqVol\n");
        return (status);
    }
    
    NbNoEqInt = NbFac  + 
                (NbFac - 1) * (NbFac) / 2;

    NbWithEqInt = (NbFac);


    NoEqInt           = (double*) DR_Array (DOUBLE,0,NbNoEqInt - 1);
    WithEqInt         = (double*) DR_Array (DOUBLE,0,NbWithEqInt - 1);
    DomA              = (double*) DR_Array (DOUBLE,0,NbFac - 1);
    DeltaTime         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    Time              = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    a                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    b                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    c                 = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA              = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_x         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_xx        = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    SpotEqVolSoFar    = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    EqVolExt          = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    EqVolUnSmile      = (double*) DR_Array (DOUBLE,0,NbTP - 1); 
    Perturb           = (double*) DR_Array (DOUBLE,0,NbTP - 1); 

    if( NoEqInt         == NULL ||
        WithEqInt       == NULL ||
        DomA            == NULL ||
        DeltaTime       == NULL ||
        Time            == NULL ||
        a               == NULL ||
        b               == NULL ||
        c               == NULL ||
        volA            == NULL ||
        loc_vol_x       == NULL ||
        loc_vol_xx      == NULL ||
        SpotEqVolSoFar  == NULL ||
        EqVolExt        == NULL ||
        EqVolUnSmile    == NULL ||
        Perturb         == NULL)
    {
        DR_Error("Memory allocation failed in MultiFac_Spot_EqVol2\n");
        status = SUCCESS;
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1;i++)
    {   
    
        if ( TPDate[i+1] <= TPDate[i])
        {
            sprintf(ErrorMsg,"TPDate[%d] and TPDate[%d] in input time line"
                      "are not in a strictly ascending order\n",i,i+1);
            DR_Error (ErrorMsg);
            return (status);
        }
    
        /* time step between time points, and corresponding times in years */
        DeltaTime[i]  =    Daysact(TPDate[i],TPDate[i+1]) / 365.;
        Time[i]       =    Daysact(TPDate[0],TPDate[i]) / 365.; 
        a[i]          =    A_Eq[i];
        b[i]          =    B_Eq[i];
        c[i]          =    C_Eq[i];

        /* atm slope wrt log-fwd-money */  
        loc_vol_x[i]    =   c[i] * a[i];  
        /* atm curvature wrt log-fwd-money */  
        loc_vol_xx[i]   =   c[i] * c[i] * b[i];

    }

    /* we have already checked that NbTp >= 2 */
    Time[NbTP-1] = Daysact(TPDate[0],TPDate[NbTP-1]) / 365.; 

    /* find the first time point on or after the FIRST imp vol date */
    
    TPFirstImpVol = GetDLOffset(NbTP,
                                TPDate,
                                VolDate[0],
                                CbkHIGHER);
    
    if (TPFirstImpVol < 0) 
    {
        sprintf(ErrorMsg, "Vol Integration date %ld is beyond the"
            "input time line (MultiFac_Spot_EqVol2)\n",VolDate[0]);
        DR_Error(ErrorMsg);
        goto RETURN;
    }

    TPLastImpVol = GetDLOffset(NbTP,
                               TPDate,
                               VolDate[NbImpVol-1],
                               CbkHIGHER);

    if (TPLastImpVol < 0)
    {
        sprintf(ErrorMsg, "Vol Integration date %ld is beyond the"
            "input time line (MultiFac_Spot_EqVol)\n",VolDate[NbImpVol-1]);
        DR_Error(ErrorMsg);
        goto RETURN;
    }
/*End additional variables */


    /*EQS*/
    /******************************************/
    /* remove 1st order EQ smile contribution */
    /******************************************/
    
    /* first estimate zero expiry limit of implied vol = EqVolUnSmile[0]              */ 
    /* the estimate is obtained from Eq. 8 of "Long-dated FX smile model"             */
    /* by assuming EQ spot vol constant and approximating RHS                         */

    /* time to first option expiry */
    tau_0 =  Daysact(TPDate[0],VolDate[0]) / 365.;
    
    /* approximation assumes const EQ spot vol until first expiry */
    epsilon = 2. * tau_0 * EqVol[0] * EqVol[0] * 
              ( loc_vol_xx[0] - 0.5 * loc_vol_x[0] * loc_vol_x[0] ) / 3.0;

    if (epsilon < -0.5) /* sqr root of negative number */
    {
        sprintf(ErrorMsg, "MultiFac_Spot_EqVol: can't calculate initial EQ spot vol\n");
        DR_Error(ErrorMsg);
        goto RETURN;        
    }  
    else if (fabs(epsilon) < TINY) 
    {
        EqVolUnSmile[0] = EqVol[0];
    }
    else 
    {
        EqVolUnSmile[0] = EqVol[0] * sqrt( (- 1. + sqrt(1 + 2. * epsilon) ) / epsilon);
    }

    /* extend the implied vols to full time line */    
    t_low     = 0.0;               /* initialise as today                               */
    delta_t   = tau_0;             /* first step size equals first expiry               */
    EqVolPrev = EqVolUnSmile[0];   /* guess for implied vol extrapolated to zero expiry */
    slope     = EqVol[0] * EqVol[0] - EqVolUnSmile[0] * EqVolUnSmile[0];    

    for (i = 0, k = 0 ; i <= NbTP - 1; i++)
    {   
        if (k <= NbImpVol - 1)
        {
            if (TPDate[i] >=  VolDate[k])       
            {
                if (k < NbImpVol - 1)
                {
                    EqVolPrev =  EqVol[k];
                    slope     =  EqVol[k+1] * EqVol[k+1] - EqVol[k] * EqVol[k];
                    t_low    +=  delta_t;
                    /* time step for next update */
                    delta_t   =  Daysact(VolDate[k],VolDate[k+1]) / 365.;
                    k++;
                }
                else 
                {
                    /*  flat after last implied vol */
                    slope     = 0.0;
                    EqVolPrev = EqVol[NbImpVol-1];
                    k++; /* i.e. k = NbImpVols */
                }    
            }
            /* interpolate linearly on squares of vols - in the short end, in particular */ 
            /* this is more accurate than interpolation on the vols                      */
            EqVolExt[i] = sqrt( EqVolPrev * EqVolPrev + slope * (Time[i] - t_low) / delta_t );
        }
        else /* flat after last implied vol */
        {         
            EqVolExt[i] = EqVol[NbImpVol - 1];
        }
    }

    /* initialise dummy for integration */
    temp1   = EqVolUnSmile[0] * EqVolUnSmile[0] * DeltaTime[0];
    volA[1] = EqVolUnSmile[0];

    /* only volA is used outside this loop */
    for (i = 1; i < TPLastImpVol; i++)
    {   
        /* if smile then use the bootstrapping routine described in "Long-dated FX" */
        if (IsSmileE)
        {    
            if (Hyb3_CalcPerturb(Perturb,                 /* (O) [0,i]...[i-1,i]                       */ 
                                 i,                       /* (I) number of time points                 */
                                 EqVolUnSmile,            /* (I) implied eq vol                        */
                                 loc_vol_x,               /* (I) slope of local vol wrt log-money      */
                                 loc_vol_xx,              /* (I) curvature of local vol wrt log-money  */
                                 Time,                    /* (I) time grid, Time[0] = 0                */
                                 DeltaTime)               /* (I) intervals in grid                     */
                                 == FAILURE) goto RETURN;

            /* EqVol[i] etc is the implied vol at Time[i] - not Time[i+1] */
            temp2   =    ( Time[i+1] * EqVolExt[i+1] * EqVolExt[i+1] - 
                           Time[i] * EqVolExt[i] * EqVolExt[i] ) / 
                         ( 1.0 + Perturb[0] );
        }

        /* otherwise, no perturbation in the implied volatility surface due to the equity smile */
        else
        {
            /* EqVol[i] etc is the implied vol at Time[i] - not Time[i+1] */
            temp2   =    ( Time[i+1] * EqVolExt[i+1] * EqVolExt[i+1] - 
                           Time[i] * EqVolExt[i] * EqVolExt[i] );
        }

        if (temp2 < 0.0)
        {
            DR_Error("MultiFac_Spot_EqVol: negative eq vol in calibration\n");
            goto RETURN;
        }
        else
        {
            EqVolUnSmile[i] = sqrt(temp2/DeltaTime[i]);
        }
        temp1   +=  EqVolUnSmile[i] * EqVolUnSmile[i] * DeltaTime[i];

        /* Time[i] > 0 because i > 0 and TPDate is strictly increasing  */
        /* volA[i] refers to Time[i]                                    */
        volA[i+1]  =  sqrt(temp1/Time[i+1]); 
    }
/************* 1st order EQ smile contribution removed ******************/
/****************/
/*EQS*/

        /*     Bootstrap equity volatility curve. */
        /*     Change to multifac format          */

    TPIntBefore = TPIntSoFar = 0;
    Sigma = 0;

    for (kVol = 0; kVol < NbImpVol; kVol++)
    {

/*New reading of TPIntSoFar, TPIntBefore */
/*******************************************************************/
        if ( VolDate[kVol] < TPDate[0])
        {
            DR_Error("VolIntegrationdate cannot be before TPDate[0].\n");
            goto RETURN;
        }

        /* check that VolIntegrationDate is on input time line         */

        /* - GetDLOffset searches for the new TPIntSoFar               */
        /*   in the range {TPDate[TPIntSoFar],TPDate[NbTP-1]}          */
        /* - since function returns offset relative to TPIntSoFar,     */
        /*   we add it to old value of TPIntSoFar                      */
        NextTP = GetDLOffset(NbTP - TPIntSoFar,
                             TPDate + TPIntSoFar,
                             VolDate[kVol],
                             CbkEXACT);

        TPIntSoFar += NextTP;

        if (NextTP < 0)
        {
            sprintf(ErrorMsg, "Vol Integration date %ld is not on the"
                    "input time line (or dates are not in ascending order)"
                    "(MultiFac_Spot_EqVol)\n",VolDate[kVol]);
            DR_Error(ErrorMsg);
            goto RETURN;
        }
        
       /* check that Voldate is on input time line                    */

       /* the following should not normally happen as fx   */ 
       /* dates are input using New_Eq_Input_W,            */
       /* which checks for strictly ascending order and    */
       /* eliminates dates before EqBaseDate=TPDate[0]     */

        if (TPIntSoFar == TPIntBefore)
        {
            continue;
        }

/*****************************************************************/

        /*EQS**********************/
        /*Initalize all quantiites */
            

        /* first,some initialisations*/
        /*DomA, domestic or foreing factor */
        for (i = 0; i < NbFac;i++)
        {
            DomA[i] = 0.;
        } /* for i*/

        for (i = 0; i < NbNoEqInt; i++)
        {
            NoEqInt[i] = 0.;
        } /* for i*/
        
        for (i = 0; i < NbWithEqInt; i++)
        {
            WithEqInt[i] = 0.;
        } /* for i*/

        for (i = 0; i < 3; i++)
        {
            Poly[i] = 0.;
        } /* for i*/
        /* prepare A-factors    */

        /*New TPIntSoFar-TPIntBefore ************************/
        for (i = TPIntSoFar - 1; i >= TPIntBefore; i--)
        {
            CurrentNbInt = -1;
            
            /* domestic A's*/
            for (l = 0; l < NbFac; l++)
            {
                DomA[l] *= exp(-Beta[l] * DeltaTime[i]);                 /*B*/
                DomA[l] += FwdRate[i] * Hyb3_ExpDecay(Beta[l],DeltaTime[i]);/*B*/

            }/* for l */
            
            /* variance of domestic factors*/
            for (l = 0;l < NbFac;l++)
            {   
                CurrentNbInt++;
                NoEqInt[CurrentNbInt] += Alpha[l] * Alpha[l] * 
                                         SpotVol[i] * SpotVol[i]
                                         * DomA[l] * DomA[l] * DeltaTime[i];
            }
                 
            /*covariance between Eq and domestic*/
            for (l = 0; l < NbFac; l++)
            {
                WithEqInt[l] += 2.0 * SpotVol[i] * Alpha[l] * 
                                DomA[l] * RhoIrEq[l] * DeltaTime[i];

            }/* for l*/
        }   /* for i*/


        for (i = TPIntBefore - 1; i >= 0; i--)
        {               
            CurrentNbInt = -1;

          /* domestic A's*/
            for (l = 0; l < NbFac; l++)
            {
                DomA[l] *= exp(-Beta[l] * DeltaTime[i]);
                DomA[l] += FwdRate[i] * Hyb3_ExpDecay(Beta[l],DeltaTime[i]);

            }/* for l */
                                    
            /* variance of each domestic factor*/
            for (l = 0; l < NbFac; l++)
            {   
                CurrentNbInt++;
                NoEqInt[CurrentNbInt] += Alpha[l] * Alpha[l] * 
                                         SpotVol[i] * SpotVol[i]
                                         * DomA[l] * DomA[l] * DeltaTime[i];
            }
                        
            /* constant coefficients due to covar fx/dom*/
            
            for (l = 0; l < NbFac; l++)
            {
                Poly[0] += 2.0 * SpotVol[i] * Alpha[l] * DomA[l] *
                                                 SpotEqVolSoFar[i] 
                                               * RhoIrEq[l] * DeltaTime[i];
            }

            /* constant coeffs due to fx variance*/
            Poly[0] += SpotEqVolSoFar[i] * SpotEqVolSoFar[i] * DeltaTime[i];
        }/* for i */


/*********************************New version of the second part of the bootstrapping *****/
        for (l = 0; l < NbNoEqInt; l++)
        {
            Poly[0] += NoEqInt[l];
        }

        /* 1) we have assumed TPDate[0] = EqBaseDate                         */
        /* 2) volA is defined on the "finer" time line  {TPDate}             */
        /*    so we use the index TPIntSoFar, for which                      */  
        /*    TPdate[TPIntSoFar] = VolDate[kVol]                             */
        
        Poly[0] -= volA[TPIntSoFar] * volA[TPIntSoFar] * 
                            Daysact(TPDate[0],TPDate[TPIntSoFar]) / 365.;

        for (l = 0; l < NbWithEqInt; l++)
        {
            Poly[1] += WithEqInt[l];
        }
        
        /* notice that TPDate[TPIntBefore] < VolIntegrationDate[kVol]        */
        /* i.e (Poly[2] non zero)                                            */
        /* because of earlier check on TPIntBefore==TPIntSoFar               */

        Poly[2] = (double) Daysact(TPDate[TPIntBefore],VolDate[kVol])
                           / 365.;

        discriminant = Poly[1] * Poly[1] - 4.* Poly[2] * Poly[0];

        VolTooLow = (discriminant < TINY);
  
        if (discriminant >= TINY)
        {       
			/* using the larger root ensures a continuous spot fx vol        */ 
            Sigma  = (-Poly[1] + sqrt(discriminant)) / 2. / Poly[2];
                            
            /* Evaluate whether spotvol satisfies appropriate criterion */
            if((EqCutOffLast) || (!EqCutOffFlag))
            {
                VolTooLow  = (EqVol[kVol] > EqCutOffLevel * Sigma);
            }
            else
            {
                VolTooLow = (Sigma < EqCutOffLevel);
            }
        }
        /* now implement cutoff actions as the case may be */

        if (VolTooLow)
        {
            if (EqCutOffFlag)
            {
                if (EqCutOffLast)
                {
                    if(kVol == 0)
                    {
                        DR_Error("Unable to bootstrap at least one FX\n"
                                "vol point and therefore unable to\n"
                                "cut off at last spot vol level.\n");

                        goto RETURN;
                    }
                    else
                    {   
                        for (n = kVol; n < NbImpVol; n++)
                        {
                            SpotEqVol[n] = SpotEqVol[kVol-1];
                        }
                    }
                }/* end of EqCutOffLast=true */

                else
                {                           
                    for(n = kVol; n < NbImpVol; n++)
                    {                            
                        SpotEqVol[n] = EqCutOffLevel;
                    }
                }

                break;

            }/* end of EqCutOffFlag */

            else
            {
                sprintf(ErrorMsg,"Problem in calculating eq vol at %ld:less"
                    " than %4.2lf%% \n(level determined by ratio to fwd vol)"
                    ".\n(MultiFac_Eq_Spot_Vol)\n",
                    VolDate[kVol],100.*EqVol[kVol]/EqCutOffLevel);
                DR_Error(ErrorMsg);
                goto RETURN;
            }
        } /* end of VolTooLow */

        SpotEqVol[kVol] = Sigma;

        for (i = TPIntBefore; i < TPIntSoFar; i++)
        {
            SpotEqVolSoFar[i] = Sigma;

        }/* for i*/

        TPIntBefore  = TPIntSoFar;

    
    }/*for kVol*/



    status = SUCCESS;

RETURN:

    Free_DR_Array (NoEqInt,          DOUBLE,0,NbNoEqInt - 1);
    Free_DR_Array (WithEqInt,        DOUBLE,0,NbWithEqInt - 1);
    Free_DR_Array (DomA,             DOUBLE,0,NbFac - 1);
    Free_DR_Array (DeltaTime,        DOUBLE,0,NbTP - 1);
    Free_DR_Array (Time,             DOUBLE,0,NbTP - 1);
    Free_DR_Array (a,                DOUBLE,0,NbTP - 1);     
    Free_DR_Array (b,                DOUBLE,0,NbTP - 1);    
    Free_DR_Array (c,                DOUBLE,0,NbTP - 1);   
    Free_DR_Array (volA,             DOUBLE,0,NbTP - 1);
    Free_DR_Array (loc_vol_x,        DOUBLE,0,NbTP - 1);   
    Free_DR_Array (loc_vol_xx,       DOUBLE,0,NbTP - 1);  
    Free_DR_Array (SpotEqVolSoFar,   DOUBLE,0,NbTP - 1);
    Free_DR_Array (EqVolExt,         DOUBLE,0,NbTP - 1);
    Free_DR_Array (EqVolUnSmile,     DOUBLE,0,NbTP - 1);   
    Free_DR_Array (Perturb,          DOUBLE,0,NbTP - 1);  

    return (status);

}  /* Hyb3_Eq_Spot_Vol */



/*************** Hyb3_MultiFac_EqVol2 *********************************************/
/* Outputs as many implied eq vols as there are Vol Expiry dates (i.e NbVol).   */
/* The correlation factors are with respect to the exponential factors          */ 
/*                                                                              */
/* NOTES:                                                                       */
/* 1-Assumes that TPDate[0]=fx_data.Today.                                      */
/* 2-Assumes and checks that each VolDate and EqRstDate is on the time          */
/* line.                                                                        */
/* 3-VolDate and EqMatDate must be entered in a strictly ascending              */
/* order.                                                                       */
/*                                                                              */
/*****************************************************************************/

int Hyb3_MultiFac_EqVol2(
					double  *ImpVol,           /*(O) Implied FX Vol for FX Expiries       */
					double  *SpotVolDom,       /*(I) Instant.dom int.rate vol             */
					double  *AlphaDom,         /*(I) factor weights                       */
					double  *BetaDom,          /*(I) domestic mean reversions             */
					double  *DomFwdRate,       /*(I) dom rates used in Christian's A facs */
                    int      DomNbFac,         /*(I) Nb factors for Dom yield curve       */
					double  *RhoEqDomFac,      /*(I) Corr. between fx and Dom factors     */
					double  *SpotEqVol,        /*(I) Inst fx vol                          */
					double  *A_Eq,
					double  *B_Eq,
					double  *C_Eq,
					long    *EqMatDate,        /*(I) maturity date                        */
					long    *EqRstDate,        /*(I) fx reset date                        */
					int      NbVol,            /*(I) Nb of output ImpVols                 */
					long    *TPDate,           /*(I) Dates on time line for vol integrat. */
                    int      NbTP)             /*(I) Nb of time points                    */
{

    int           i, k,l, status=FAILURE;

    int
                  NbNoFxInt,           /* Nb of integrals not containing spot fx vol  */
                  NbWithFxInt,         /* Nb of integrals due to covar fx/Dom,fx/For  */
		          IndexRst,IndexMat,   /* FX reset and maturity indices on timeline   */
                  CurrentNbInt;        /* used to index integrlas not involving Fx    */

    double
                  *NoFxInt        = NULL,    /* Integrals not containing spot fxvol         */
                  *WithFxInt      = NULL,    /* Integrals due to covar Fx/Dom and Fx/For    */
                  *DomA           = NULL,
                  *DeltaTime      = NULL,
				  *Time           = NULL,
	          	  *a              = NULL,    /* short for A_Fx                        */
         		  *b              = NULL,    /* short for B_Fx                        */
            	  *c              = NULL,    /* short for C_Fx                        */
				  *volA           = NULL,    /* lowest order approximation            */
				  *volA2          = NULL,    /* square of lowest order approx         */
        		  *volA_x         = NULL,    /* slope wrt log moneyness               */
        		  *volA_xx        = NULL,    /* curvature wrt log moneyness           */
				  *loc_vol_x      = NULL,    /* local vol slope                       */
				  *loc_vol_xx     = NULL;    /* local vol curvature                   */ 
	
    double        
		          varFX,
                  volB_sum = 0.;
            
    char
                  ErrorMsg[MAXBUFF];
    

    /* Quick checks first */
    if (DomNbFac != 1 && DomNbFac != 2 && DomNbFac != 3)
     
    {
        DR_Error("Invalid Number of domestic yield curve factors "
                 "in MultiFac_FxVol: should be 1,2"
                 " or 3\n");
        return (status);
    }


    if( NbTP < 2)
    {
        DR_Error ("Invalid number of time points in MultiFac_FxVol: "
                  "should be at least 2\n");
        return (status);
    }

    if( NbVol < 1)
    {
        DR_Error ("Invalid number of Vol points in MultiFac_fxVol: "
                  "should be at least 1");
        return (status);
    }

    if (ImpVol             == NULL ||
        SpotVolDom         == NULL ||
        AlphaDom           == NULL ||
        BetaDom            == NULL ||
        DomFwdRate         == NULL ||
        RhoEqDomFac        == NULL ||
        SpotEqVol          == NULL ||
        EqRstDate          == NULL ||
        TPDate             == NULL)

    {
        DR_Error ("invalid pointer inputs to MultiFx_FxVol\n");
        return (status);
    }



    NbNoFxInt = (DomNbFac ) + 
                (DomNbFac - 1) * (DomNbFac ) / 2;

    NbWithFxInt = (DomNbFac );


    NoFxInt         = (double*) DR_Array (DOUBLE,0,NbNoFxInt - 1);
    WithFxInt       = (double*) DR_Array (DOUBLE,0,NbWithFxInt -1);
    DomA            = (double*) DR_Array (DOUBLE,0,DomNbFac - 1);
    DeltaTime       = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    Time            = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    a               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    b               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    c               = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA            = (double*) DR_Array (DOUBLE,0,NbTP - 1);
	volA2           = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA_x          = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    volA_xx         = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_x       = (double*) DR_Array (DOUBLE,0,NbTP - 1);
    loc_vol_xx      = (double*) DR_Array (DOUBLE,0,NbTP - 1);

    if( NoFxInt      == NULL ||
        WithFxInt    == NULL ||
        DomA         == NULL ||
    
        DeltaTime    == NULL ||
		Time         == NULL ||
		a            == NULL ||
		b            == NULL ||
		c            == NULL ||
		volA         == NULL ||
		volA2        == NULL ||
		volA_x       == NULL ||
		volA_xx      == NULL ||
		loc_vol_x    == NULL ||
		loc_vol_xx   == NULL)
    {
        goto RETURN;
    }

    for (i = 0; i < NbTP - 1; i++)
    {   

        if (TPDate[i+1] <= TPDate[i])
        {

            sprintf(ErrorMsg,"TPDate[%d] and TPDate[%d] in input time line"
                             "are not in a strictly"
                             " ascending order\n",i,i+1);
            DR_Error (ErrorMsg);
            return (status);

        }

	    /* time step between time points, and corresponding times in years */
        DeltaTime[i] = Daysact(TPDate[i],TPDate[i+1]) / 365.;
  	    Time[i] = Daysact(TPDate[0],TPDate[i]) / 365.; 

		/* change eq smile arrays to short hand */
		a[i]    =    A_Eq[i];
		b[i]    =    B_Eq[i];
		c[i]    =    C_Eq[i];
		
		/* atm slope wrt log-fwd-money */  
		loc_vol_x[i]    =   c[i] * a[i];  
		/* atm curvature wrt log-fwd-money */  
		loc_vol_xx[i]   =   c[i] * c[i] * b[i];
			         
    }

    Time[NbTP-1] =  Daysact(TPDate[0],TPDate[NbTP-1]) / 365.; 

	IndexMat = 0;  /* TPDate[IndexMat] = FxMatDate[k] */

    for (k = 0; k < NbVol; k++)
    {
       	
        /* check that FxRstDate[i] <= FxMatDate[i] */
        if (EqRstDate[k] > EqMatDate[k])
        {

            DR_Error ("EQ reset date is after maturity date\n");
                return(status);

        }

        /* find maturity date on time line */
        IndexMat = GetDLOffset(NbTP,
                               TPDate,
                               EqMatDate[k],
                               CbkEXACT);		

        if (IndexMat < 0)
        {

            sprintf(ErrorMsg, "Option expiry date #%ld is not"
                              " on the input time line\n",EqMatDate[k]);
            DR_Error(ErrorMsg);
            goto RETURN;

        }
        
		/* reset date = base date */
        if (EqRstDate[k] == TPDate[0])
        {

            ImpVol[k] = 0.0;
            continue; /* next k */

        }

        /* find FxRstDate[k] on time line*/        
        IndexRst = GetDLOffset(NbTP,
                               TPDate,
                               EqRstDate[k],
                               CbkEXACT);	

        if (IndexRst < 0)
        {

            sprintf(ErrorMsg, "Option expiry date #%ld is not"
                              " on the input time line\n",
                              EqRstDate[k]);
            DR_Error(ErrorMsg);
            goto RETURN;

        }

    	/*   now TPDate[IndexRst] == FxRstDate[k] */
        
        if (IndexRst > IndexMat )    
        {

            DR_Error("IndexRst is strictly greater than"
                     "IndexMat\n");
            goto RETURN;

        }

        /* Some initialisations */
        varFX = 0.;

        for(l = 0; l < DomNbFac; l++)
        {

           DomA[l]      =   0.;

        }/* for i*/
        

        for(l = 0; l < NbNoFxInt; l++)
        {

            NoFxInt[l]   =   0.;

        }

        for(l = 0; l < NbWithFxInt; l++)
        {
        
            WithFxInt[l] =   0.;

        }
        /******************************************/
        /*      initialisations are complete      */
        /******************************************/  
        
        for(i = IndexMat - 1; i >= IndexRst; i--)
        {
                        
            /* domestic A's*/

            for (l = 0; l < DomNbFac; l++)
            {

				DomA[l]   *=   exp(-BetaDom[l] * DeltaTime[i]);
                DomA[l]   +=   DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            

        }/* for i */

        for (i = IndexRst - 1; i >= 0;i --)
        {

            CurrentNbInt     = 0;

           /* domestic A's*/
            for (l = 0; l < DomNbFac; l++)
            {

                DomA[l]   *=   exp(-BetaDom[l] * DeltaTime[i+1]);
                DomA[l]   +=   DomFwdRate[i] * Hyb3_ExpDecay(BetaDom[l],DeltaTime[i]);

            }/* for l */
            
      
            /* variance of domestic factors */
            for (l = 0; l < DomNbFac; l++)
            {   

                CurrentNbInt++;
                NoFxInt[CurrentNbInt-1]    +=    AlphaDom[l] * AlphaDom[l] * 
                                                 SpotVolDom[i] * SpotVolDom[i]
                                                 * DomA[l] * DomA[l] * DeltaTime[i];

            }
            
            /* covariance between domestic factors*/
            
            

         
            
          
            
            /* covariance between Fx and domestic */
            for (l = 0; l < DomNbFac; l++)
            {
                WithFxInt[l]             +=   2. * SpotEqVol[i] * SpotVolDom[i]
                                              * AlphaDom[l] * DomA[l] * 
                                              RhoEqDomFac[l] * DeltaTime[i];

            }/* for l*/

        

			/* integral of "composite" vol square over (t_i,T) */
            varFX      +=  SpotEqVol[i] * SpotEqVol[i] * DeltaTime[i];  
			volA2[i]    =  varFX; 
            
     		for(l = 0; l < NbNoFxInt; l++)
            {
				
				volA2[i]    +=  NoFxInt[l];
				
            }

			for(l = 0; l < NbWithFxInt; l++)
            {
				
				volA2[i]    +=  WithFxInt[l];
				
            }

        } /* i */
		/* finished IR contribution */ 
		
		/* Note that volA[i] = \Sigma_A(t_i,T), where T = Time[IndexRst] */
		for (i = IndexRst - 1; i >= 0; i--)
        {
			volA[i]      = sqrt( volA2[i] / (Time[IndexRst] - Time[i]) );
		} 

	    /* set up array of curvatures */
        if (Hyb3_CalcVolA_xx(volA_xx,      /* (O)                                      */
 			                 IndexRst,     /* (I) number of time points                */
                             SpotEqVol,    /* (I) spot fx vol                          */
                             loc_vol_x,    /* (I) slope of local vol wrt log-money     */
				             loc_vol_xx,   /* (I) curvature of local vol wrt log-money */
                             Time,         /* (I) time grid, evaluation time = 0       */
				             DeltaTime)    /* (I) intervals in grid                    */
                             == FAILURE) goto RETURN;

		volB_sum = 0.;
		for (i = IndexRst - 1; i >= 0; i--)
        {
			/* integrate result */
			/* using volA[i-1] because volA[i] corresponds to Time[i+1] */

			/* the RHS equals T - t integrated over (t_i,t_{i+1}) */
			volB_sum +=  SpotEqVol[i] * SpotEqVol[i] * volA[i] * volA_xx[i] 
				                      * 0.5 * ( - Time[i+1] * Time[i+1] + Time[i] * Time[i] 
							          + 2.0 * Time[IndexRst] * (Time[i+1] - Time[i]) );
        }   /* for i */
        
		/* Since the implied vol approximation is unsafe for large pertubations, an error   */
		/* is returned in those cases - the user should then lower implied vol curvature.   */                                                                  
		/* Having the cutoff at 100% for vol^2 terms is somewhat arbitary                   */ 
        if( fabs(volB_sum) > 3.0 * volA2[0])
        {
            sprintf(ErrorMsg, "Volatility or its curvature is too large, bootstrapping" 
                              "fails at expiry #%ld\n",EqMatDate[k]);

            DR_Error(ErrorMsg);
            goto RETURN;
        }

		ImpVol[k] = sqrt( (volA2[0] + volB_sum) / (Time[IndexRst] - Time[0]) ); 
        
        if( ImpVol[k] < TINY)
        {
            sprintf(ErrorMsg, "Problem in calculating Implied EQ Vol at"
                              "Expiry #%ld\n",EqMatDate[k]);

            DR_Error(ErrorMsg);
            goto RETURN;
        }

    }/* for k*/

    status = SUCCESS;

RETURN:

    Free_DR_Array (NoFxInt,     DOUBLE,0,NbNoFxInt - 1);
    Free_DR_Array (WithFxInt,   DOUBLE,0,NbWithFxInt - 1);
    Free_DR_Array (DomA,        DOUBLE,0,DomNbFac - 1);
    Free_DR_Array (DeltaTime,   DOUBLE,0,NbTP - 1);
    Free_DR_Array (Time,        DOUBLE,0,NbTP - 1);
	Free_DR_Array (a,           DOUBLE,0,NbTP - 1);     
	Free_DR_Array (b,           DOUBLE,0,NbTP - 1);    
	Free_DR_Array (c,           DOUBLE,0,NbTP - 1);   
	Free_DR_Array (volA,        DOUBLE,0,NbTP - 1);
	Free_DR_Array (volA2,       DOUBLE,0,NbTP - 1);       	 
	Free_DR_Array (volA_x,      DOUBLE,0,NbTP - 1);     
	Free_DR_Array (volA_xx,     DOUBLE,0,NbTP - 1);         
	Free_DR_Array (loc_vol_x,   DOUBLE,0,NbTP - 1);   
	Free_DR_Array (loc_vol_xx,  DOUBLE,0,NbTP - 1);  
	
    return (status);

}/* Hyb3_Multifac_FxVol2*/




/*****  Hyb3_Eq_Vol  **********************************************************************
*
*  Computes the composite Eq Vols at every time point of the tree timeline.
*  
*  Effectively, it calls the MultiFac_EqVol2 which 
*  a) uses the analytical approximation based on Henrik's paper
*  b) is the exact same function as in srm3
*
*
*******************************************************************************************/

int	Hyb3_Eq_Vol(double   *EqVolCurve,        /* (O) Eq volatility at each time step */
		   double   *SpotEqVol,         /* (I) Inst eq vol at each time step   */
           double   *RhoCurve,       /* (I) Correlations at each time step  */
           double   *SpotVol,        /* (I) Inst vol of short rate at node  */
           double   *FwdRateDom,        /* (I) Pass Forward Domestic           */
           double   DomAlpha[3],        /* (I) dom factor weights              */ 
           long     *VolIntegrationDate,/* (I) expiry of option                */
           double   *A1C,               /* (I) FX smile param, [0,NbTP-2]      */
           double   *A2C,               /* (I) FX smile param, [0,NbTP-2]      */
           double   *A3C,               /* (I) FX smile param, [0,NbTP-2]      */
           double   BetaDom[3],         /* (I) Mean reversion coefficient (dom)*/
           long     *TPDate,			/* (I) date of each pt.                */
           int      NbTP)               /* (I) Total nb of nodes or time steps */
{

	int
        k,
        status = FAILURE;        /* Error status = FAILURE initially */

	double *SDomSpotVol = NULL;

	SDomSpotVol    =   (double *)DR_Array(DOUBLE, -1, NbTP+1);

    if (SDomSpotVol == NULL)
    {
        DR_Error("Unable to allocate memory.");
        goto RETURN;
    }

	/*Change spot vol format to match srm3 definitions (division by alpha factors)*/
	for (k = 0 ; k < NbTP; k++)
    {
        SDomSpotVol[k] = SpotVol[k] / DomAlpha[0];
    }


	if(Hyb3_MultiFac_EqVol2( EqVolCurve,
						SDomSpotVol,            
						DomAlpha,        
						BetaDom,         
						FwdRateDom,
						1,                                    
						RhoCurve,                 
						SpotEqVol,      
						A1C,
						A2C,
						A3C,
						VolIntegrationDate,       
						VolIntegrationDate,        
						(NbTP+1),                    
						TPDate,          
						(NbTP+1)) != SUCCESS){
	
		goto RETURN;
	}

    /*Add -1 point*/
	EqVolCurve[-1] = EqVolCurve[0];


	status = SUCCESS;

	RETURN:
	if(SDomSpotVol!=NULL){
		Free_DR_Array(SDomSpotVol, DOUBLE, -1, NbTP+1);
    }

    if (status == FAILURE)
    {
        DR_Error("Eq_Vol: failed");
    }
	return (status);
}




/***** Hyb3_Get_TreeEqSpotVol *****************************************************
* CAREFUL: The Input IR vols are assumed to be in an "Aweight" form 
* i.e "alpha*spotvol". 
* 
* 1-Assumes that TPDate ranges from TPDate[-1] to TPDate[NbTP - 1]
*   The last meaningful TreeSpotVol is TreeSpotVol[NbTP - 2], 
*   The first meaningful TreeSpotVol is TreeSpotVol[0].
*   TreeSpotVol is then extended by:
*   TreeSpotVol[NbTP - 1] ,which is set to TreeSpotVol[NbTP - 2],and
*   TreeSpotVol[-1], which is set to TreeSpotVol[0];
*
* 2-The array of InpCompVol is offset by one, i.e, the first
*   effective expiry is InpVompVol[1], the last one is InpCompVol[NbCompVol]
*
* 3-The array of InpSpotVol is NOT offset by one, i.e, 
*   the first spot vol is InpSpotVol[0], and the last one is 
*   InpSpotVol[NbInpSpotVol - 1]  
*
*****************************************************************************/
int Hyb3_Get_TreeEqSpotVol(double    *TreeSpotVol, /*(O)Bootst'd spotvol on tree    */
                      double    *FwdEq,       /*(I)Forward stock for equity    */
                      long      *TPDate,      /*(I)tree time pts               */
                      long      NbTP,         /*(I)Last TP is TPDate[NbTP-1]   */
                      double    *InpCompVol,  /*(I)input composite vols        */
                      long      *InpCompVDate,/*(I)option expiry dates         */
                      int       NbCompVol,    /*(I)Number of composite vols    */
                      double    *InpSpotVol,  /*(I)Input Spot Eq vols          */
                      long      *InpSpotVDate,/*(I)Input SpotVol Dates         */
                      int       NbInpSpotVol, /*(I)Nb of Input spot vols       */
                      double    *DomSpotVol,  /*(I)dom ir (tree) Aweight       */
                      double    *DomAlpha,    /*(I)dom factor weights          */ 
                      double    *DomFwdRate,  /*(I)dom (tree) fwd rates        */
                      double    *RhoEqDom,    /*(I)                            */
                      int       EqCutOffFlag, /*(I)                            */
                      int       EqCutOffLast, /*(I)                            */
                      double    EqCutOffLevel,/*(I)                            */
                      double    *DomBeta,     /*(I)                            */
                      int       DomNbFac,     /*(I)                            */
                      long      Today,        /*(I)                            */
                      int       EqBootStrapMode,
					  double    *A_Eq,        /* (I) FX smile param, [0,NbTP-2]*/
					  double    *B_Eq,        /* (I) FX smile param, [0,NbTP-2]*/
					  double    *C_Eq)        /* (I)                           */
{


    int status = FAILURE;
    long    PseudoCompVDate[MAXNBDATE+1];/*Stores interpolated expiry dates  */
    double  PseudoCompVol[MAXNBDATE+1];  /*stores interpolated composite vols*/
    long    *RawSpotVDate   = NULL;      /*spot vol dates before extension   */
    double  *RawSpotVol     = NULL;      /*spotvols bef. extension on TPDate */
    int     NbRawSpotVol    = 0;     /*total Nb of spot vols bef. extension  */
    int     NbPseudoCompVol = 0;     /*Nb of valid interpolated compositevols*/


    int     i,k;
    long    TPSofar, TPBefore;
    double  Vol;
    double  *LDomVol = NULL;        /* dom Aweight/alpha*/  
    double  *LForVol = NULL;        /* for Aweight/alpha*/

    /* basic checks */
    if (TreeSpotVol == NULL ||
        TPDate      == NULL ||
        DomSpotVol  == NULL ||
        DomAlpha    == NULL ||
        DomBeta     == NULL ||
        DomFwdRate  == NULL ||
        RhoEqDom    == NULL )
    {
        DR_Error("invalid input pointers to Get_TreeEqSpotVol\n");
        goto RETURN;
    }
    
    if (NbTP < 2) goto RETURN;

    /* Do the nil calibration case and exit*/
    if (EqBootStrapMode == EQ_CONSTANT_SPOT_VOL)
    {
        for (k = -1 ; k < NbTP ; k++)
        {
            TreeSpotVol[k] = EqCutOffLevel;
        }

        return (SUCCESS);
    }

    if ((NbInpSpotVol < 0) || (NbCompVol < 0)) goto RETURN;
    if ((NbCompVol > 0) && ((InpCompVol == NULL) || (InpCompVDate == NULL)))
    {        
        goto RETURN;
    }


    if ((NbInpSpotVol > 0) && ((InpSpotVol == NULL) || (InpSpotVDate == NULL)))
    {
        goto RETURN;
    }
    
    /* Check that composite vol dates have been read in using Fx_Input_W*/
    /* which offsets vols and dates by 1 and sets VolDate[0] to Today   */
    if (NbCompVol > 0)
    { 
        if(InpCompVDate[0] != Today)
        {
            DR_Error("First Composite Vol Date (%ld) is not Today (%ld)\n",
                     InpCompVDate[0],
                     Today);
            goto RETURN;
        }
        if(!IS_EQUAL(InpCompVol[0],InpCompVol[1])) goto RETURN;
    }
    
    
    if (TPDate[0] != Today)
    {
        DR_Error("TPDate[0] is not Today \n");
        goto RETURN;
    }

    /* Prepare Ir vols */
    if (DomAlpha[0] < TINY) goto RETURN;

    LDomVol = (double*) DR_Array(DOUBLE,0,NbTP - 2);

    if (LDomVol == NULL) goto RETURN;

    for (k = 0 ; k < NbTP - 1 ; k++)
    {
        LDomVol[k] = DomSpotVol[k] / DomAlpha[0];

    }


    /* Find Pseudo composite vols and expiries */
    /*************************************************/
    /* Input composite vols are offset by 1, whereas */
    /* Pseudocomposite vols are offset by 0          */
    /*************************************************/

    TPSofar = TPBefore = 0;
    NbPseudoCompVol    = 0;

    for (k = 1; k <= NbCompVol ; k++) /*input composite vols are offset by 1*/
    {
        /* Look for the tree time  point  TPSoFar  falling  */ 
        /* immediately before the current volatility input  */         
        TPSofar = NbTP - 1;
        while ((TPDate[TPSofar] > InpCompVDate[k]) && (TPSofar > 0))  
                    TPSofar--;        	                   
        
        /* Ignore this point if it is falling onto the previous one */        
        if (TPSofar == TPBefore)	                               
            continue;

        /* Interpolate fx volatility at current node as the input */
        /* volatility points may not be falling on a node.        */    
        dlinterp (TPDate[TPSofar],     &Vol,                     
                  InpCompVDate[k-1] ,  InpCompVDate[k],	           
                  InpCompVol[k-1]   ,  InpCompVol[k]);
        
        PseudoCompVDate[NbPseudoCompVol] = TPDate[TPSofar];
        PseudoCompVol[NbPseudoCompVol]   = Vol;

        NbPseudoCompVol++;
        TPBefore = TPSofar;
    }

    /* quick checks*/    
    if ((NbPseudoCompVol > NbCompVol) || (NbPseudoCompVol < 0)) goto RETURN;    
    if ((NbInpSpotVol > 0) && (NbPseudoCompVol > 0))
    {
        if(PseudoCompVDate[NbPseudoCompVol-1] >= InpSpotVDate[0]) goto RETURN;
    }
    
    if (NbPseudoCompVol > 0)
    {
        if(PseudoCompVDate[0] <= Today) goto RETURN;
    }

    if (NbPseudoCompVol > 1)
    {
        for (i = 0 ; i < NbPseudoCompVol - 1; i++)
        {
            if (PseudoCompVDate[i] >= PseudoCompVDate[i+1]) goto RETURN;
        }
    }

    NbRawSpotVol = NbPseudoCompVol + NbInpSpotVol;
    if (NbRawSpotVol < 1) 
    {
        DR_Error("No vols to bootstrap (no spot vols and "
                "no composite vols left after interpolation!)\n");
        goto RETURN;
    }
    
    RawSpotVDate    = (long*)   DR_Array(LONG,  0, NbRawSpotVol - 1);
    RawSpotVol      = (double*) DR_Array(DOUBLE,0, NbRawSpotVol - 1);

    if ((RawSpotVDate == NULL) || (RawSpotVol == NULL)) goto RETURN;


    /* fill in vol dates and spot vols */
    for (k = 0 ; k < NbPseudoCompVol; k++)
    {
        RawSpotVDate[k] = PseudoCompVDate[k];
    }
  
    if (NbPseudoCompVol > 0)
    {   
    

        if (Hyb3_MultiFac_Spot_EqVol(RawSpotVol,
                                FwdEq,
                                DomNbFac,
                                DomAlpha,
                                DomBeta,
                                LDomVol,
                                DomFwdRate,
                                RhoEqDom,
                                TPDate,
                                NbTP,
                                NbPseudoCompVol,
                                PseudoCompVDate,
                                PseudoCompVol,
                                EqCutOffFlag,
                                EqCutOffLast,
                                EqCutOffLevel,
		    		 			A_Eq,
			 					B_Eq,
								C_Eq) == FAILURE) goto RETURN;  

        
       

    }/* if there are composite vols to bootstrap*/
    
    /* Append inputspot vol dates and vols*/
    for (k = NbPseudoCompVol ; k < NbRawSpotVol ; k++)
    {
        RawSpotVDate[k] = InpSpotVDate[k - NbPseudoCompVol];
        RawSpotVol[k]   = InpSpotVol[k - NbPseudoCompVol];
    }
    
    /*EQS*/
    /*Before, nb of points for the extend function was not NbTP*/
    if(ExtendSpotVol(TreeSpotVol,
                     NbRawSpotVol,
                     RawSpotVDate,
                     RawSpotVol,
                     (NbTP + 1),
                     TPDate) == FAILURE) goto RETURN;


    /* extend by one */
    TreeSpotVol[-1] = TreeSpotVol[0];
    status = SUCCESS;

RETURN:
    
    if (NbRawSpotVol >= 1)
    {
        Free_DR_Array(RawSpotVDate,LONG,0,NbRawSpotVol - 1);
        Free_DR_Array(RawSpotVol, DOUBLE,0,NbRawSpotVol - 1);
    }
    if (NbTP >= 2)
    {
        Free_DR_Array(LDomVol, DOUBLE, 0, NbTP - 2);
        Free_DR_Array(LForVol, DOUBLE, 0, NbTP - 2);
    }
    
    return(status);
}



/*****  not used with VHT smile ***************************************************/
/*****  Hyb3_Eq_Shift_Centre  *****************************************************/
/*
 *      Given the initial positioning of the centre of the equity dimension
 *      on the equity forward, this routine shifts the centre so that it is
 *      placed on the centre of the lognormal distribution,  rather than on 
 *      the equity forward.
 *
 *      This  "shifting" operation is kept as a separate module, because in 
 *      the Q(t) smile formulation for equity, we actually want to centre
 *      the equiy dimension on the forward.
 *
 */
int	Hyb3_Eq_Shift_Centre(
        double   *EqMidNode,      /* (I/O) Center of the equity tree         */
        double   *SpotEqVol,      /* (I) Spot equity volatility at time pt   */
        double   *Length,         /* (I) Length of each time step            */
        int       NbTP)           /* (I) Total number of time points         */
{

        double    Vol2;      
        int       i;           
        int       status = FAILURE;   /*  Error status = FAILURE initially   */


        /*	Shift center of equity tree */
        Vol2 = 0.;
        for (i = 0; i <= NbTP; i++)
        {       
                
            EqMidNode[i] *= exp (-.5 * Vol2);
            Vol2 += SpotEqVol[i] * SpotEqVol[i] * Length[i];
 
        } 

        status = SUCCESS;
        return (status);

}  /* Hyb3_Eq_Shift_Centre */


/*****  not used with VHT smile **************************************************/
/*****	Hyb3_Eq_Spot_Vol  ********************************************************/
/*
 *	Interpolate the correlation  curve at each node  in the tree. 
 *  Calculate  equity spot volatility at each node in the tree by
 *	bootstrapping input benchmark option volatility curve, taking
 *	into account interest rate volatility.
 */
int   Hyb3_Eq_Spot_Vol 
    (double  *SpotEqVol,   /* (O) Spot equity volatility at each time step   */
     double  *EqSmiQ,      /* (I) Single Q as function of time               */
     double  *FwdEq,       /* (I) Forward price at each time step            */
     double  *RhoCurve,    /* (I) Correlation at each time step              */
     double  *SpotVol,     /* (I) Spot vol of short rate at each time step   */
     double  *FwdRate,     /* (I) Forward rate at each node in the tree      */
     double  *Length,      /* (I) Length of each time step                   */
     long    *TPDate,      /* (I) Date of each node                          */
     long     ValueDate,   /* (I) Value date                                 */
     int      NbVol,       /* (I) Number of points in equity volatility curve*/
     long    *VolDate,     /* (I) Benchmark equity option expiration dates   */
     double  *EqVol,       /* (I) Benchmark equity option volatilities       */
     double   AnchorStrike,/* (I) Strike at which to map to Q measure        */
     double   Beta,        /* (I) Mean reversion coefficient                 */
     int      NbTP)        /* (I) Total number of nodes (or time steps)      */
{

    double
            Sigma,             /* Equity spot volatility in current bucket  */
            Vol,               /* Current equity volatility                 */
            T,                 /* Time to current calculation in years      */
            atmPr,             /* Black-Scholes option price                */
            strike,            /* strike at which to match the Q vols       */
            I[4],              /* 8 integrals needed to calculate spot vol  */
            B,                 /* A coefficient in Christian's memo         */
            a[3];              /* Coefficients of the 2nd degree polynomial */
    int	
            i1,                /* Node index of the previous vol point      */
            i2,                /* Node index of the current vol point       */
            i = 0, 
            j, 
            k,
            status = FAILURE;  /* Error status = FAILURE initially          */
    char
            ErrorMsg[MAXBUFF]; /* Error message to be sent to DR_Error      */



        /* 	Bootstrap equity volatility curve. */

        i1 = i2 = 0;
        Sigma = 0;

        for (k = 1; k <= NbVol; k++)
        {

            /* Look for the node index i2 falling immediately */
            /* before the current volatility point.           */
            i2 = NbTP;
            while ((TPDate[i2] > VolDate[k]) && (i2 > 0))
                    i2--;
            
            /* Ignore this point if it is falling onto the previous one */
            if (i2 == i1)
                    continue;
            
            dlinterp(TPDate[i2],   &Vol,      /* Interp vol at current node  */
                     VolDate[k-1], VolDate[k],/* as the input vol points may */
                     EqVol[k-1],   EqVol[k]); /* not be falling on a node.   */


            T = Daysact(ValueDate, TPDate[i2])/365.;

            /* Unless running on lognormal case, convert vols to Q meas. */
            if (fabs(EqSmiQ[i2] - 1.0) > TINY)
            {
                /* Lognormal to X-space vol conversion. The single Q */
                /* case is treated  using the general, 2Q  functions */
                /* for convenience.                                  */
                if (AnchorStrike > TINY) /* A strike was entered */
                {
                    strike = AnchorStrike;
                }
                else                     /* Else anchor at spot  */
                {
                    strike = FwdEq[0];
                }

                atmPr = Option_BS2Q(FwdEq[i2],
                                    strike, 
                                    T,
                                    Vol,
                                    'C',1.,1.,0.); /* lognormal call */

                Vol = ImpVol_BS2Q(FwdEq[i2],
                                  strike,
                                  T,
                                  atmPr,
                                  'C',EqSmiQ[i2],EqSmiQ[i2],0., /* Q call */
                                  Vol);

		if ( Vol < 0.0 )
                {
                    DR_Error("Unable to convert stock vol to 2Q measure.\n");
                    goto RETURN;
                }
            }


            B = 0.;
            for (j = 0; j < 4; j++)
                    I[j] = 0.;
            
            for (i = i2 - 1; i >= i1; i--)
            {
                B *= exp (-Beta * Length[i]);
                B += FwdRate[i] * Hyb3_ExpDecay(Beta, Length[i]);
                                    
                I[0] += SpotVol[i] * SpotVol[i] * B * B * Length[i];
                I[3] += 2. * RhoCurve[i] * SpotVol[i] * B * Length[i];
                    
            }  /* for i */

            
            for (; i >= 0; i--)
            {
                B *= exp (-Beta * Length[i]);
                B += FwdRate[i] * Hyb3_ExpDecay(Beta, Length[i]);
                        
                I[0] += SpotVol[i] * SpotVol[i] * B * B * Length[i];

                I[1] += 2.*RhoCurve[i] * SpotEqVol[i] * SpotVol[i]
                        * B * Length[i];
                I[2] += SpotEqVol[i] * SpotEqVol[i] * Length[i];
                    
            }  /* for i */
        
            
            a[2] = (double) Daysact (TPDate[i1], TPDate[i2]) / 365.;
            a[1] = I[3];
            a[0] = - Vol * Vol * Daysact (ValueDate, TPDate[i2]) / 365.;
            
            for (j = 0; j < 3; j++)
                    a[0] += I[j];


            if (a[1]*a[1] - 4.*a[2]*a[0] < 0.)
            {
                sprintf (ErrorMsg, "Problem in bootstrapping equity vol "
                                   "#%d (Hyb3_Eq_Spot_Vol(.))", k);
                DR_Error (ErrorMsg);
                goto RETURN;
            }

            Sigma = (-a[1] + sqrt (a[1]*a[1] - 4.*a[2]*a[0])) / 2. / a[2];
            
            if (Sigma < .02)  /* Arbitrary number: spot volatility > 2% */
            {

                sprintf (ErrorMsg, "Problem in calculating equity vol "
                                   "#%d: less than 2%% (Hyb3_Eq_Spot_Vol(.))", k);
                DR_Error (ErrorMsg);
                goto RETURN;
            }  
        
            /* Spot volatility is constant between benchmark pts   */
            for (i = i1; i < i2; i++)
                    SpotEqVol[i] = Sigma;    
        
            i1 = i2;

        }  /* for k */


        for (; i <= NbTP; i++)       /* If the tree extends beyond the */ 
            SpotEqVol[i] = Sigma;    /* last input vol point we fill   */
                                     /* with the last value found.     */

        /* We need one extra node to build the lattice */
        SpotEqVol[-1] = SpotEqVol[0];


        status = SUCCESS;

        RETURN:

        return (status);

}  /* Hyb3_Eq_Spot_Vol */


/*****  not used with VHT smile **************************************************/
/*****	Hyb3_Eq_Vol_Old***********************************************************/
/*
*	Re-calculate equity volatility at each node in the tree from
*	the spot volatility, the interest rate spot volatilities and
*	correlation curve.
*/
int	Hyb3_Eq_Vol_Old
    (double   *EqVolCurve, /* (O) Equity option volatility at each time step */
     double   *SpotEqVol,  /* (I) Spot equity volatility at each time step   */
     double   *EqSmiQ,     /* (I) Single Q as function of time               */
     double   *FwdEq,      /* (I) Forward price at each time step            */
     double   *RhoCurve,   /* (I) Correlation at each time step              */
     double   *SpotVol,    /* (I) Spot vol of short rate at each time step   */
     double    AnchorStrike,/* (I) Strike at which to map to Q measure       */
     double   *FwdRate,    /* (I) Forward rate at each node in the tree      */
     double   *Length,     /* (I) Length of each time step                   */
     double    Beta,       /* (I) Mean reversion coefficient                 */
     int       NbTP)       /* (I) Total number of nodes (or time steps)      */
{


    double
            I[3],          /* 3 integrals needed to calculate current vol    */
            *Ad = NULL,    /* Intermediate integrals to calculate vol curve  */
            atmPr,         /* Option price used to imply vol                 */
            strike,        /* strike at which to match the Q vols            */
            T;
    int
            i, j, k,
            status = FAILURE;      /* Error status = FAILURE initially     */
    char
            ErrorMsg[MAXBUFF];     /* Error message to be sent to DR_Error */
             

    Ad = DR_Array(DOUBLE, 0, NbTP);
    if (Ad == NULL)
    {
        goto RETURN;
    }
            

    I[2] = 0.;

    for (i = 0; i <= NbTP; i++)
    {
        EqVolCurve[i+1] = 0.;
        I[0] = I[1] = 0.;
        
        for (j = i, T = 0.; j >=0; j--)   /* Double integral */
        {        
            Ad[j] += FwdRate[i] * Hyb3_ExpDecay(Beta, Length[i]) * exp (-Beta * T);
                                                                            
            I[0] += SpotVol[j] * SpotVol[j] * Ad[j] * Ad[j] * Length[j];
            I[1] += 2.*RhoCurve[j]*SpotEqVol[j]*SpotVol[j]*Ad[j]*Length[j];
                                        
            T += Length[j];
                
        }  /* for j */
        
        /* Integral of spot equity volatility: doesn't need the double integral */
        I[2] += SpotEqVol[i] * SpotEqVol[i] * Length[i];   
        
        for (k = 0; k <= 2; k++)     /* The eq vol is the sum of 3 integrals  */
            EqVolCurve[i+1] += I[k]; /* Note the 1 offset for EqVolCurve[i+1] */


        if (EqVolCurve[i+1] < ERROR) /* Shouldn't happen as spot vol is +tive */
        {
            sprintf (ErrorMsg, "Problem in calculating equity vol "
                               "(Q-vol) at node #%d (Hyb3_Eq_Vol(.))", i);
            DR_Error (ErrorMsg);
            goto RETURN;

        }


        EqVolCurve[i+1] /= T;
        EqVolCurve[i+1] = sqrt (EqVolCurve[i+1]);

        if (i==NbTP-1)
            EqVolCurve[0] = EqVolCurve[NbTP];

        /* Unless running on lognormal case, convert vols  */
        if (fabs(EqSmiQ[i] - 1.0) > TINY)
        {
            /* Convert back from Q measure to B&S measure */
            if (AnchorStrike > TINY) /* A strike was entered */
            {
                strike = AnchorStrike;
            }
            else                     /* Else anchor at spot  */
            {
                strike = FwdEq[0];
            }

            atmPr = Option_BS2Q(FwdEq[i],
                                strike,
                                T,
                                EqVolCurve[i+1],
                                'C',EqSmiQ[i],EqSmiQ[i],0.);

            EqVolCurve[i+1] = ImpVol_BS2Q(FwdEq[i],
                                strike,
                                T,
                                atmPr,
                                'C',1.,1.,0.,
                                EqVolCurve[i+1]);

	    if ( EqVolCurve[i+1] < 0.0 )
            {
                DR_Error("Unable to convert Q stock vol to lognormal measure.\n");
                goto RETURN;
            }
        }

    }  /* for i */


    status = SUCCESS;

    RETURN:

    Free_DR_Array(Ad, DOUBLE, 0, NbTP);
        
    return (status);

}  /* Hyb3_Eq_Vol_Old */



