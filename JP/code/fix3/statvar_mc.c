/****************************************************************************/
/*      Bootstrapping of swaption volatility into spot volatility.          */
/****************************************************************************/
/*      STATVAR.c                                                           */
/****************************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fix123head.h"


/* Fix3_gasdev from NR */
static double Fix3_gasdev(long *idum);


/*****  Fix3_IndexLimits_Classic  *****************************************/
/*
*      Find the maximum and minimum yields at NBSTDEV move               
*/
int     Fix3_IndexLimits_Classic (
          double*       LoYield,              /* (O) Min yield at NBSTDEV move    */
          double*       HiYield,              /* (O) Max yield at NBSTDEV move    */
          long          NbStDev,              /* (I) Number of std for state var  */
          long          RateReset,            /* (I) Time to expiration           */
          long          SwapSt,               /* (I) Underlying swap start        */
          int           IndexMat,             /* (I) Index maturity               */
          char          Freq,                 /* (I) Day count convention         */    
          char          DCC,                  /* (I) Frequency of underlying rate */
          MKTVOL_DATA   *mktvol_data,         /* (I) Volatility data              */
          T_CURVE const* crv)                 /* (I) zero curve                   */
{

    double ParYld, Annuity;
    double Time, Vol;
    long   MatDate;
    long   ValueDate = mktvol_data->BaseDate;

    static	char	routine[] = "Fix3_IndexLimits";
    int status = FAILURE;

    if (SwapSt < ValueDate)
    {
        DR_Error ("%s: Swap start date %ld < value date %ld",
                  routine,
                  SwapSt, ValueDate);
        goto RETURN;
    }

	if (ParYield(&ParYld,
               &Annuity,
               crv,
               SwapSt,
               IndexMat,
               DCC,
               Freq) == FAILURE)
    {
        goto RETURN;
    }
            
    MatDate = Nxtmth(SwapSt,
                     IndexMat,
                     1L);

    if (Fix3_IndexVol (&Vol, 
                       SwapSt,
                       MatDate,
                       Freq, 
                       DCC,
                       mktvol_data,
                       crv) == FAILURE)
    {
        goto RETURN;
    }
    /* Maximal and minimal yields at NBSTDEV move */
    Time = Daysact(ValueDate, RateReset) / 365.;
    Time = MAX(Time, 1/365.);  /* to treat 1% vol cases */
    Vol *= NbStDev * sqrt(Time);

    *LoYield = ParYld * exp(-Vol);
    *HiYield = ParYld * exp(+Vol);

    status = SUCCESS;

    RETURN:

    return (status);

}/*Fix3_Index_Limits_Classic */

/*****  Fix3_DoubleArrayFloorIdx  ************************************************/
/*
*      Find the floor index "i" in a double array y[n] given search value x
*      such that y[i] <= x.
*      if x < y[0], i = -1
*/
int    Fix3_DoubleArrayFloorIdx(double const* y, /* (I) Double ascending array */
                                int           n, /* (I) Size of array          */
                                double        x) /* (I) Search value           */
{
        register int jl, ju, jm;
        int    j;    /* return result */

        jl = (-1);
        ju = n;
        while (ju-jl > 1) {
                jm=(ju+jl) >> 1;
                if (x >= y[jm])
                        jl=jm;
                else
                        ju=jm;
        }
        j=jl;

        return(j);
}




/*****  Fix3_DoubleVectSort  ************************************************/
/*
*      Sort input double array in ascending order.
*/
int     Fix3_DoubleVectSort(double *x, int n)
{
static	char	routine[] = "Fix3_DoubleVectSort";

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;
#define M 7
#define NSTACK 50

	unsigned	long i,ir=n,j,k,l=1;
	int		jstack=0,
			istack[NSTACK+3];
	double		*arr, a, temp;

	/* NRC convention */
	arr = x - 1;

	for (;;) {
		if (ir-l < M) {
			for (j=l+1;j<=ir;j++) {
				a=arr[j];
				for (i=j-1;i>=1;i--) {
					if (arr[i] <= a) break;
					arr[i+1]=arr[i];
				}
				arr[i+1]=a;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(arr[k],arr[l+1])
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[l]) {
				SWAP(arr[l+1],arr[l])
			}
			i=l+1;
			j=ir;
			a=arr[l];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j]);
			}
			arr[l]=arr[j];
			arr[j]=a;
			jstack += 2;
			if (jstack > NSTACK) {
				DR_Error("%s: NSTACK too small.\n",
					routine);
				return(FAILURE);
			}
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}

	return(SUCCESS);

#undef M
#undef NSTACK
#undef SWAP

}    /* Fix3_DoubleVectSort */




/*****  Fix3_DoubleQuadraticInterp  **********************************************/
/*
*      Special quardratic interpolation.  
*      Input X(i) must be in ascending order
*      and duplicated elements in both X(i) and Y(i) are allowed,
*/
int    Fix3_DoubleQuadraticInterp(
    double const* xa,     /* (I) array of X(i) (i=0,..,m-1) */
    double const* ya,     /* (I) array of Y(i) (i=0,..,m-1) */
    int           m,      /* (I) arrays size */
    double        x,      /* (I) point to intepolate */
    double*       y)      /* (O) interpolated value */
{
    static  char    routine[] = "Fix3_DoubleQuadraticInterp";
    
    int    i;
    double d0, d1, d2;
    char   interpType;

    /* Input check
     * if X(i)=X(j), then Y(i)=Y(j)    
     * Duplication of X(i) is allowed 
     */
    for (i=1; i<m; i++)
    {
        if ( IS_EQUAL(xa[i], xa[i-1]) &&
            !IS_EQUAL(ya[i], ya[i-1]))
        {
            DR_Error("%s: X[%d] = X[%d], but Y[%d] != Y[%d]",
                     routine,
                     i-1, i, i-1, i);
            return FAILURE;
        }
    }

    
    /* Flat on both ends */
    if (x >= xa[m-1])
    {
        *y = ya[m-1];
        return SUCCESS;
    }
    
    if (x <= xa[0])
    {
        *y = ya[0];
        return SUCCESS;
    }

    /* Interpolate in between */

    /* 1. Two points only */
    if (m == 2)
    {
        linterp (x, 
                 y, 
                 xa[0], xa[1], 
                 ya[0], ya[1]);

        return SUCCESS;
    }
    
    /* Search the corresponding index */
    i = Fix3_DoubleArrayFloorIdx(xa, m, x);

    /* 2. Between two flat points */
    if (IS_EQUAL(xa[i], xa[i+1])) 
    {
        *y = ya[i];
        return SUCCESS;
    }

    /* 3. Linear if both left and right buckets are flat,
     *    otherwise quadratic interpolation involving current period and
     *    either the one on the right or on the left
     */
    if (i == m-2)        /* Right-most bucket */
    {
        if (IS_EQUAL(xa[i-1], xa[i]))  /* left bucket is flat */
        {
            interpType = 'L';
        }
        else    /* quadratic involving i-1, i, i+1 */
        {
            interpType = 'Q';
            i = i-1;
        }
    }
    else if (i == 0)    /* Left-most bucket */
    {
        if (IS_EQUAL(xa[i+1], xa[i+2]))  /* right bucket is flat */
        {
            interpType = 'L';
        }
        else    /* quadratic involving i, i+1, i+2 */
        {
            interpType = 'Q';
        }
    }
    else                /* In the middle  */
    {
        if (!IS_EQUAL(xa[i+1], xa[i+2]))  /* right bucket is NOT flat */
        {
            interpType = 'Q';
        }
        else if (!IS_EQUAL(xa[i-1], xa[i]))  /* left bucket is NOT flat */
        {
            interpType = 'Q';
            i = i-1;
        }
        else        /* Both right and left buckets are flat */
        {
            interpType = 'L';
        }
    }


    /* Interpolate based on interp type */

    if (interpType == 'L')     /* Linear */
    {
        linterp (x, 
                 y, 
                 xa[i], xa[i+1], 
                 ya[i], ya[i+1]);

        return SUCCESS;
    }
    else                         /* Quadratic */
    {
        d0 = 1. / ((xa[i]-xa[i+1])
                  *(xa[i]-xa[i+2]));
        d1 = 1. / ((xa[i+1]-xa[i])
                  *(xa[i+1]-xa[i+2]));
        d2 = 1. / ((xa[i+2]-xa[i])
                  *(xa[i+2]-xa[i+1]));
 
        sqinterp(xa[i],
                 xa[i+1],
                 xa[i+2],
                 ya[i],  
                 ya[i+1],
                 ya[i+2],
                 d0,
                 d1,
                 d2,
                 x, 
                 y);

        return SUCCESS;
    }

}    /* Fix3_DoubleQuadraticInterp */




/*****  Fix3_Initialize_IR_SIM ****************************************************/
/*
*       Initialize IR_SIM structure
*/
int     Fix3_Initialize_IR_SIM (
          IR_SIM  *ir_sim,               /* (I/O) IR simulation data         */
          long    BaseDate,              /* (I) Simulation base date         */
          int     NbPaths)               /* (I) Nb of paths                  */
{
    int    p;

    ir_sim->NbPaths  = NbPaths;

    ir_sim->BaseDate = BaseDate;
    ir_sim->CurrDate = BaseDate;

    if (NbPaths > MAXNBSTATES)
    {
        DR_Error("Fix3_Initialize_IR_SIM: NbPaths (%d) > %d the maximum of "
                 "state variables allowed",
                 NbPaths,
                 MAXNBSTATES);
        return (FAILURE);
    }

    for (p=0; p<MAXNBSTATES; p++)
        ir_sim->X[p] = 0.0;

    /* Dummy 3M rate, not used.
     * This will be updated by actual simulated swap rates 
     * along the way
     */
    ir_sim->SwapSt   = BaseDate;
    ir_sim->SwapMat  = Nxtmth(BaseDate, 3, 1L);
    ir_sim->Freq     = 'Q';
    ir_sim->DCC      = '0';

    /* Seed for Fix3_gasdev */
    ir_sim->Seed     = 0L;

    return (SUCCESS);
}


    
    
/*****  Fix3_SwapYield_MC *********************************************************/
/*
*       Swap yield distributions along the path
*/
int     Fix3_SwapYield_MC (
          double  SwapYield[MAXNBSTATES],/* (O) Swap yield slice             */ 
          double  *MaxYield,             /* (O) Max of yield with NbStDev    */ 
          double  *MinYield,             /* (O) Min of yield with NbStDev    */ 
          int     NbStDev,               /* (I) NbStDev for Max/Min yield    */
          long    SwapSt,                /* (I) Rate1 start                  */
          long    SwapMat,               /* (I) Rate1 maturity               */
          char    Freq,                  /* (I) Rate1 frequency              */
          char    DCC,                   /* (I) Rate1 Day count convention   */
          IR_SIM  *ir_sim,               /* (I) IR simulation data           */
          MKTVOL_DATA *mktvol_data,      /* (I) Volatility data              */
          T_CURVE const* crv)            /* (I) zero curve                   */
{
static  char    routine[] = "Fix3_SwapYield_MC";

    int     status = FAILURE;  /* Error status = FAILURE initially      */

    double  Corr,
            Var1,
            Var2;              /* variances and correlation btw two rates  */

    double  Q,
            QMid,
            ParYield,
            Annuity,
            FwdParYield,
            St,
            Center,
            Sigma;

    int     p;

    /* temporary values for volatility data param*/ 
    double  QLeft;                 /* (I) Left Q mapping coefficient   */
    double  QRight;                /* (I) Right Q mapping coefficient  */
    double  FwdShift;              /* (I) Fwd shift mapping coefficient*/
   


    /* set temporary parameters to mktvol_data parameters */
    QLeft       = mktvol_data->QLeft;
    QRight      = mktvol_data->QRight;
    FwdShift    = mktvol_data->FwdShift;
    
    /* Check simulation time line */
    if (SwapSt < ir_sim->CurrDate)
    {
        DR_Error ("%s: Swap state date %8ld < current simulation date %8ld",
                  routine,
                  SwapSt, ir_sim->CurrDate);
        goto RETURN;
    }  
    
    
    /* Compute the variance and correlation between two rates */
    if (Fix3_IndexCovar (
              &Corr,        
              &Var1,
              &Var2,
              SwapSt,      
              SwapMat,    
              Freq,      
              DCC,      
              ir_sim->SwapSt,  
              ir_sim->SwapMat,
              ir_sim->Freq,  
              ir_sim->DCC,  
              mktvol_data,
              crv) != SUCCESS)
    {
        goto RETURN;
    }


    /* Compute forward par yield */
    if (ParYieldFromDates(
                  &ParYield,
                  &Annuity,
                  SwapSt,
                  SwapMat,
                  DCC,
                  Freq,
                  'F',
                  crv) != SUCCESS)
    {
        goto RETURN;
    }



    /* Compute the drift and vol of swap */
    QMid   = (QLeft + QRight) / 2.;
    FwdParYield = ParYield / (1. + QMid * FwdShift);
    St = sqrt(Var1) * (1 + FwdShift) / (1 + QMid * FwdShift);


    /* Drift */    
    if (ConvexityC_BS2Q(St,
                        QLeft,
                        QRight,
                        FwdShift,
                        &Center) == FAILURE)
    {
        goto RETURN;
    }

    /* Max/Min yield of given NbStDev */
    Sigma = sqrt(Var1) * NbStDev;

    if (IS_Q(QLeft))
    {
        *MinYield = FwdParYield * (1. - 1. / QLeft
                  + exp(-QLeft * Sigma + QLeft * Center) / QLeft);
    }
    else
    {
        *MinYield = FwdParYield * (1. - Sigma);
    }

    if (IS_Q(QRight))
    {
        *MaxYield = FwdParYield * (1. - 1. / QRight
                  + exp(QRight * Sigma + QRight * Center) / QRight);
    }
    else
    {
        *MaxYield = FwdParYield * (1. + Sigma);
    }


    /* Compute the conditional yield distribution given the
     * known swap rate distribution and correlation of OU grid
     * with the current swap rate.
     *
     * This is a simplified mkt3 model.
     */
    for (p=0; p<ir_sim->NbPaths; p++)
    {
        ir_sim->X[p] = Corr * ir_sim->X[p]
                     + sqrt(1. - Corr * Corr) * Fix3_gasdev(&(ir_sim->Seed));
        
        Sigma = sqrt(Var1) * ir_sim->X[p];

        if (Sigma >= -Center)
            Q = QRight;
        else
            Q = QLeft;

        if (IS_Q(Q))
        {
            SwapYield[p] = FwdParYield * (1. - 1. / Q
                         + exp(Q * Sigma + Q * Center) / Q);
        }
        else
        {
            SwapYield[p] = FwdParYield * (1. + Sigma);
        }
    }


    /* Update the current swap info */
    ir_sim->CurrDate  = SwapSt;
    ir_sim->SwapSt    = SwapSt;
    ir_sim->SwapMat   = SwapMat;
    ir_sim->DCC       = DCC;
    ir_sim->Freq      = Freq;

    
    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
        DR_Error("%s: failed", routine);

    return (status);

}  /* Fix3_SwapYield_MC */
    





/*****  Fix3_IndexCovar *********************************************************/
/*
*       Compute variances and correlation between two rate indices.
*/
int     Fix3_IndexCovar (
          double  *Corr,                 /* (O) Correlation btw Rate1 & Rate2*/
          double  *Var1,                 /* (O) Variance of Rate1            */
          double  *Var2,                 /* (O) Variance of Rate2            */
          long    SwapSt1,               /* (I) Rate1 start                  */
          long    SwapMat1,              /* (I) Rate1 maturity               */
          char    Freq1,                 /* (I) Rate1 frequency              */
          char    DCC1,                  /* (I) Rate1 Day count convention   */
          long    SwapSt2,               /* (I) Rate2 start                  */
          long    SwapMat2,              /* (I) Rate2 maturity               */
          char    Freq2,                 /* (I) Rate2 frequency              */
          char    DCC2,                  /* (I) Rate2 Day count convention   */
          MKTVOL_DATA *mktvol_data,      /* (I) Volatility data              */
          T_CURVE const* crv)            /* (I) zero curve                   */
{

static  char    routine[] = "Fix3_IndexCovar";
    int     status = FAILURE;  /* Error status = FAILURE initially      */

    

    double  deltaT;              /* Time between SwapSt1 and SwapSt2    */

    double  Cov;

    double  L1[3][3],
            L2[3][3],
            LC[3][3];          /* Integrals of factor spot vol         */

    double  B1[3], 
            B2[3];             /* B in Christian memo                   */

    int     i, j;

    /* temporary variables for volatility parameters */
    int     NbFactor;
    double  *Beta;

    /* set temporary variables to mktvol_data parameters*/
    NbFactor = mktvol_data->NbFactor;
    Beta     = mktvol_data->Beta;

    

   
    /* Find B coefficients of each index */
    if (Fix3_BFactor (B1,
                 SwapSt1,
                 SwapMat1,
                 DCC1,
                 Freq1,
                 mktvol_data,
                 crv) != SUCCESS || 
        Fix3_BFactor (B2,
                 SwapSt2,
                 SwapMat2,
                 DCC2,
                 Freq2,
                 mktvol_data,
                 crv) != SUCCESS)
    {
        goto RETURN;
    }


    /* Vol integral up to SwapSt1 */
    Fix3_LCoef ( L1,
            SwapSt1,
            mktvol_data); 

    /* Vol integral up to SwapSt2 */
    Fix3_LCoef ( L2,
            SwapSt2,
            mktvol_data); 


    /* Covariance up to MIN(SwapSt1, SwapSt2) */
    if (SwapSt1 <= SwapSt2)
    {
        deltaT = (double) Daysact(SwapSt1, SwapSt2) / 365.;
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++) 
        {
            LC[i][j] = L1[i][j];
        }    
    }
    else
    {
        deltaT = (double) Daysact(SwapSt2, SwapSt1) / 365.;
        for(i = 0; i < 3; i++)
        for(j = 0; j < 3; j++) 
        {
            LC[i][j] = L2[i][j];
        }    
    }
    

    *Var1 = 0.0;
    *Var2 = 0.0;
    Cov   = 0.0;
    /* Compute variaces and covariance */
    for (i=0; i<NbFactor; i++)
    for (j=0; j<NbFactor; j++)
    {
        *Var1 += B1[i] * B1[j] * L1[i][j];
        *Var2 += B2[i] * B2[j] * L2[i][j];
        Cov   += B1[i] * B2[j] * LC[i][j] * exp (-Beta[j] * deltaT);

    }


    if (*Var1 < SQUARE(TINY) ||
        *Var2 < SQUARE(TINY) ||
        Cov   < SQUARE(TINY))
    {
        *Corr = 0.0;
    } 
    else
    {
        *Corr = Cov / sqrt((*Var1) * (*Var2));

        if (fabs(*Corr) > 1.0 + TINY)
        {
            DR_Error("%s: correlation (%lf) > 1", 
                     routine,
                     *Corr);
            goto RETURN;
        }    
    }
        

    status = SUCCESS;

    RETURN:

    if (status == FAILURE)
        DR_Error("%s: failed", routine);

    return (status);

}  /* Fix3_IndexCovar */





/*****  Fix3_LCoef    *********************************************************/
/* integrand: exp(-(beta[i]+beta[j])*(T-t))*sigma^2(t)
*
*       Integrate the spot vol.
*/
int     Fix3_LCoef (
          double  L[3][3],               /* (O) L coef                       */
          long    SwapSt,                /* (I) Underlying swap start        */
          MKTVOL_DATA *mktvol_data)      /* (I) Volatility data              */
          
{
    double  VolT[MAXNBDATE];   /* Expiries in years                     */
    double  t;                 /* Time between two consecutive expiries */

    long    EndDate;           /* End of current integration bucket     */
    int     i, j;

    /* temporary values for volatility parameters */
    int     CalibFlag;             /* (I) Index calibration flag       */
    int     NbVol;                 /* (I) Number of spot vol points    */
    long    VolBaseDate;           /* (I) Volatility curve base date   */
    long    *VolDate;              /* (I) Spot vol dates               */
    double  Aweight[6][MAXNBDATE]; /* (I) Aweight curve                */
    int     NbFactor;              /* (I) Number of factors            */
    double  *Beta;                 /* (I) Mean reversions              */
    /* Initialization */
    for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++) 
    {
        L[i][j] = 0.;
    }    

    
    /* set temporary values to mktvol_data parameters */
    CalibFlag   = mktvol_data->CalibFlag;
    NbVol       = mktvol_data->NbVol;
    VolBaseDate = mktvol_data->BaseDate;
    VolDate     = mktvol_data->VolDate;
    NbFactor    = mktvol_data->NbFactor;
    Beta        = mktvol_data->Beta;

    if (VolBaseDate >= SwapSt)
        return SUCCESS;

    /* set temporary Aweight to mktvol_data Aweight */
    if (CalibFlag == TRUE)
    {
        for (i = 0; i < NbVol; i++)
        {
            Aweight[0][i] = mktvol_data->Aweight[0][i]; 
            if (NbFactor > 1)
            {
                Aweight[1][i] = mktvol_data->Aweight[1][i];
                Aweight[2][i] = mktvol_data->Aweight[2][i];
            }
            if (NbFactor > 2)
            {
                Aweight[3][i] = mktvol_data->Aweight[3][i];
                Aweight[4][i] = mktvol_data->Aweight[4][i];
                Aweight[5][i] = mktvol_data->Aweight[5][i];
            }

        }
    }
    else
    {
        Aweight[0][0] = mktvol_data->Aweight[0][0]; 
        if (NbFactor > 1)
        {
            Aweight[1][0] = mktvol_data->Aweight[1][0];
            Aweight[2][0] = mktvol_data->Aweight[2][0];
        }
        if (NbFactor > 2)
        {
            Aweight[3][0] = mktvol_data->Aweight[3][0];
            Aweight[4][0] = mktvol_data->Aweight[4][0];
            Aweight[5][0] = mktvol_data->Aweight[5][0];
        }
    }

    /* Const spot vol */
    if (CalibFlag == FALSE)
    {
        t  = Daysact (VolBaseDate, SwapSt) / 365.;

        L[0][0] = Aweight[0][0] * Aweight[0][0] * Fix3_ExpDecay(2. * Beta[0], t) * t;
    
        if (NbFactor > 1) 
        {
            L[1][1]  = Aweight[1][0] * Aweight[1][0] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[1][1] += Aweight[2][0] * Aweight[2][0] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1]  = Aweight[0][0] * Aweight[1][0] * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
            L[1][0]  = Aweight[1][0] * Aweight[0][0] * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;

        }
        
        if (NbFactor > 2)
        {
            L[2][2]  = Aweight[3][0] * Aweight[3][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[4][0] * Aweight[4][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[5][0] * Aweight[5][0] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2]  = Aweight[0][0] * Aweight[3][0] * Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[2][0]  = Aweight[3][0] * Aweight[0][0] * Fix3_ExpDecay ((Beta[2]+Beta[0]), t) * t;
            L[1][2]  = Aweight[1][0] * Aweight[3][0] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[1][2] += Aweight[2][0] * Aweight[4][0] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[2][1]  = Aweight[3][0] * Aweight[1][0] * Fix3_ExpDecay ((Beta[2]+Beta[1]), t) * t;
            L[2][1] += Aweight[4][0] * Aweight[2][0] * Fix3_ExpDecay ((Beta[2]+Beta[1]), t) * t;

        }
    
        return (SUCCESS);
    }


    /* Time-dependent vol term structure */
    /* Vol integral up to SwapSt */
    for (i = 0; i < NbVol; i++)
    {
        /* End of integration bucket: if not last bucket, then stop at end */
        /* of bucket or expiry, whichever is smaller; if last bucket, stop */
        /* at expiry (as spot volatility is constant after last bucket).   */

        EndDate = (i < NbVol-1) ? MIN(SwapSt, VolDate[i]) : SwapSt;

        /* Time to expiry in years */
        VolT[i] = Daysact (VolBaseDate, EndDate) / 365.;

        t = ((i == 0) ? VolT[i] : (VolT[i] - VolT[i-1]));

        /* 
         *  Update integrals
         */
        L[0][0] *= exp (-2. * Beta[0] * t);
        L[0][0] += Aweight[0][i] * Aweight[0][i] * Fix3_ExpDecay (2. * Beta[0], t) * t;

        if (NbFactor > 1)
        {
            L[1][1] *= exp (-2. * Beta[1] * t);
            L[1][1] += Aweight[1][i] * Aweight[1][i] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[1][1] += Aweight[2][i] * Aweight[2][i] * Fix3_ExpDecay (2. * Beta[1], t) * t;
            L[0][1] *= exp (-(Beta[0] + Beta[1]) * t);
            L[0][1] += Aweight[0][i] * Aweight[1][i] * Fix3_ExpDecay ((Beta[0]+Beta[1]), t) * t;
            L[1][0] *= exp (-(Beta[1] + Beta[0]) * t);
            L[1][0] += Aweight[1][i] * Aweight[0][i] * Fix3_ExpDecay ((Beta[1]+Beta[0]), t) * t;
        }

        if (NbFactor > 2)
        {
            L[2][2] *= exp (-2. * Beta[2] * t);
            L[2][2] += Aweight[3][i] * Aweight[3][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[4][i] * Aweight[4][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[2][2] += Aweight[5][i] * Aweight[5][i] * Fix3_ExpDecay (2. * Beta[2], t) * t;
            L[0][2] *= exp (-(Beta[0] + Beta[2]) * t);
            L[0][2] += Aweight[0][i] * Aweight[3][i] * Fix3_ExpDecay ((Beta[0]+Beta[2]), t) * t;
            L[2][0] *= exp (-(Beta[2] + Beta[0]) * t);
            L[2][0] += Aweight[3][i] * Aweight[0][i] * Fix3_ExpDecay ((Beta[2]+Beta[0]), t) * t;
            L[1][2] *= exp (-(Beta[1] + Beta[2]) * t);
            L[1][2] += Aweight[1][i] * Aweight[3][i] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[1][2] += Aweight[2][i] * Aweight[4][i] * Fix3_ExpDecay ((Beta[1]+Beta[2]), t) * t;
            L[2][1] *= exp (-(Beta[2] + Beta[1]) * t);
            L[2][1] += Aweight[3][i] * Aweight[1][i] * Fix3_ExpDecay ((Beta[2]+Beta[1]), t) * t;
            L[2][1] += Aweight[4][i] * Aweight[2][i] * Fix3_ExpDecay ((Beta[2]+Beta[1]), t) * t;
        }

        /* End of integration loop */
        if (EndDate == SwapSt)
        {
            break;
        }
    }  /* for i */

    return (SUCCESS);
}




/*  Constants for Rand2  */
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/*************************************************************/
/**    UNIFORM RANDOM GENERATORS FROM NUMERICAL RECIPES     **/
/*************************************************************/
double Fix3_ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

/* Constants for Ran2 */
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX




/*************************************************************/
/**               NORMAL RANDOM GENERATORS                  **/
/*************************************************************/
/*********************/
/* BOX MULLER METHOD */ 
/*********************/
static
double Fix3_gasdev(long *idum)
{
	static int iset=0;
	static double gset;
	double fac,rsq,v1,v2;

	if  (iset == 0) {
		do {
			v1=2.0*Fix3_ran2(idum)-1.0;
			v2=2.0*Fix3_ran2(idum)-1.0;
			rsq=v1*v1+v2*v2;
		} while (rsq >= 1.0 || fabs(rsq) < TINY);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return (v2*fac);
	} else {
		iset=0;
		return (gset);
	}
}


