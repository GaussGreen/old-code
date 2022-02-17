/*********************************************************************************
 * Transition matrix for defaults
 *
 ********************************************************************************/

#include "transition.h"
#include "lapackinterface.h"

#define  INVSQRT2PI             0.398942280401433   /* 1/sqrt(2*pi) */
//#define  PI                     3.141592653589793

static int      one = 1;

/* computed from InitLFac */
static double   LFac1 = 0.13742649;
static double   LFac2 = 0.058854782;

/*********************************************************************************
 * Leverage Factor
 *
 ********************************************************************************/
static void InitLFac()
{
    double dbeta = 0.01;
    double dm = 0.1, m;
    int    nBeta = 99;
    int    i;
    double p = 1e-4, K;
    double b1, b2, beta;
    double sumFsinPi, sumFsin2Pi, sum;
    
    sumFsinPi = 0.0;
    sumFsin2Pi = 0.0;

    K = NormCumInv(p);
    
    for(i = 0;i < nBeta;i++)
    {
        beta = dbeta * i;
        b1   = sqrt(beta);
        b2   = sqrt(1.0 - beta);
        sum = 0.0;
        m = -5.0;
        while(m < 5.0)
        {
            sum += exp(-m*m/2) * pow(1 - NormCum((K - b1*m)/b2), 100);
            m += dm;
        }

        sum *= (INVSQRT2PI * dm);
        sum = 1.0 - sum;
        
        sum /= (100*p);
        
        sum -= (1 - beta + 0.01*beta);
        
        sumFsinPi += sum * sin(PI * beta);
        sumFsin2Pi += sum * sin(2.0 * PI * beta);        
    }

    //LFac[0] = 2.0 * sumFsinPi * dbeta;
    //LFac[1] = 2.0 * sumFsin2Pi * dbeta;
}

static double CalcAlpha(
    double          beta,                 /* (I) correlation for CreditMetrics  */
    int             numTotal)             /* (I) total num of tranchelets       */        
{
    double f;
    
    f = (1.0 - beta + beta/numTotal) + LFac1 * sin(PI * beta)
        + LFac2 * sin(2.0 * PI * beta);

    return numTotal * (1.0 - f )/(numTotal - 1.0);
}


/*********************************************************************************
 * Calculate residual expected loss
 *
 ********************************************************************************/
int ResidualLoss(
    double          *loss,                /* (O) residual expected loss         */
    TDate           date,                 /* (I) a specific date to observe     */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "ResidualLoss";

    int       status  = FAILURE;
    double    P1, P2;
    int       numName = effCurves->numName;
    int       numTotal = effCurves->numTotal;
    
    if(GtoDiscountDate(date,
                       effCurves->pCurves[numName - 2],
                       GTO_FLAT_FORWARDS,
                       &P1) != SUCCESS)
    {
        goto RETURN;
    }
       
    if(GtoDiscountDate(date,
                       effCurves->pCurves[numName - 1],
                       GTO_FLAT_FORWARDS,
                       &P2) != SUCCESS)
    {
        goto RETURN;
    }

    *loss = (1.0-P2)/(1.0-P1) * (numTotal - numName +1) + (numName - 1);
    if(*loss > numTotal) *loss = numName-1;
    
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * default transition matrix from startDate to endDate with catastrophe rate.
 *
 ********************************************************************************/
static int CalcLambda(double *lambda,
                      double *scale,
                      double *p1,
                      double *p2,
                      int    numName)
{
    static char      routine[] = "CalcLambda";
    
    int       status  = FAILURE;
    double    sum = 1.0;

    sum -= p1[numName-1];
    
    *lambda = (p2[numName-1] - p1[numName-1])/sum;

    if(*lambda < 0.0) *lambda = 0.0;
    *scale = 1.0;

    status = SUCCESS;
    return status;
}
/*********************************************************************************
 * default transition matrix from startDate to endDate with catastrophe rate.
 *
 ********************************************************************************/
static int  TransMatCat(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "TransMatCat";

    int       status  = FAILURE;
    int       i, k, ntmp;
    Mat       tmpMat;
    TDate     date1, date2;
    TDateList *tranDates = NULL;
    double    dt, lambda[MAX_NB],dtmp, presum, x;
    double    **P = NULL, **p = NULL, *rhs = NULL, *a = NULL, *p1 = NULL;
    TDate     baseDate = effCurves->today;
    int       numName2, numName;
    double    cRate;

    numName  = effCurves->numName;
    numName2 = numName * numName;
    
    /* initialize */
    MatSetSize(mat, numName);
    MatSetSize(&tmpMat, numName);
 
    for(i = 0;i < mat->size2;i++)
    {
        mat->ptr[i]        = 0.0;
        tmpMat.ptr[i]     = 0.0;
    }
    for(i = 0;i < numName;i++)
    {
        mat->ptr[i * mat->size + i] = 1.0;
    }

    if(startDate >= endDate)
    {
        return SUCCESS;
    }

    /* calculate */
    tranDates  = GtoNewDateList(
        startDate,
        endDate,
        &(effCurves->step),
        FALSE);

    if(tranDates == NULL)
    {
        CRXError("%s: Construction of date list failed.\n",routine);
        goto RETURN;
    }
    
    p = (double **)malloc(tranDates->fNumItems * sizeof(double *));
    if(p == NULL)
    {
        CRXError("%s: Memory allocation failed.\n",routine);
        goto RETURN;
    }
    p[0] = NULL;
    p[0] = (double *)malloc(tranDates->fNumItems * numName * sizeof(double));

    P = (double **)malloc(tranDates->fNumItems * sizeof(double *));
    if(P == NULL)
    {
        CRXError("%s: Memory allocation failed.\n",routine);
        goto RETURN;        
    }

    P[0] = NULL;
    P[0] = (double *)malloc(tranDates->fNumItems * numName * sizeof(double));

    rhs = (double *)malloc(numName * sizeof(double));
    a = (double *)malloc(numName * numName * sizeof(double));

    if(!p[0] || !P[0] || !rhs || !a)
    {
        CRXError("%s: Memory allocation failed.\n",routine);
        goto RETURN;
    }

    for(i = 0;i < tranDates->fNumItems;i++)
    {
        P[i] = &P[0][i * numName];
        p[i] = &p[0][i * numName];        
    }
        
    for(i = 0;i < tranDates->fNumItems;i++)
    {
        date1  = tranDates->fArray[i];
        for(k = 0;k < numName;k++)
        {
            if(GtoDiscountDate(date1,
                               effCurves->pCurves[k],
                               GTO_FLAT_FORWARDS,
                               &P[i][k]) != SUCCESS)
            {
                goto RETURN;
            }

            if(k == numName -1)
            {
                P[i][k] = 1.0;
            }
            
            p[i][k] = (k > 0)?P[i][k] - P[i][k-1]:P[i][k];
        }
    }
    
    for(i = 0;i < tranDates->fNumItems - 1;i++)
    {
        date1  = tranDates->fArray[i];
        date2  = tranDates->fArray[i+1];

        if(GtoDayCountFraction(baseDate,date1, GTO_ACT_365F,&dt) != SUCCESS)
        {
            goto RETURN;
        }

        if(GtoDiscountDateForward(
               date1,
               date2,
               effCurves->cCurve,
               GTO_FLAT_FORWARDS,
               &cRate) != SUCCESS)
        {
            goto RETURN;
        }

        cRate = 1.0/cRate - 1.0;
        
        //CalcLambda(&dtmp,&dtmp1, p[i],p[i+1],numName);

        for(k = 0;k < numName;k++)
        {
            rhs[k] = p[i+1][k] - p[i][k];
            lambda[k] = cRate;

            if(lambda[k] > 1.0) lambda[k] = 1.0;
            if(lambda[k] < SMALL) lambda[k] = 0.0;
        }
        p1  = p[i];

        presum = 0.0;
        PtrSet(tmpMat.ptr,tmpMat.size2,0.0);
        
        for(k = 0;k < numName - 1;k++)
        {
            if(p1[k] < SMALL)
            {
                x = 1.0;                  /* is it reasonable? */
                dtmp = x - lambda[k];                
            } else {
                x = - (rhs[k] - presum)/p1[k];

                /* Warning, stability, good enough? */
                if(x > 1.0) x = 1.0;
                if(x < -0.0) x = 0.0;
                
                dtmp = x - lambda[k];

                if(x < lambda[k]) dtmp = x;

                presum = dtmp * p1[k];
            }

            ntmp = k * numName;

            
            tmpMat.ptr[ntmp + k] = 1.0 - x;            
            if(k < numName -1)
            {
                tmpMat.ptr[ntmp + k + 1] += dtmp;
            }

            tmpMat.ptr[ntmp + numName -1] += x - dtmp;
        }

        tmpMat.ptr[numName2 -1] = 1.0;
        if(MatMultUT(mat, mat, &tmpMat) != SUCCESS)
        {
            goto RETURN;
        }
    }
    
    status = SUCCESS;
RETURN:
    if(tranDates != NULL) GtoFreeDateList(tranDates);

    if(p != NULL)
    {
        if(p[0] != NULL) free(p[0]);        
        free(p);
    }
    if(P != NULL)
    {
        if(P[0] != NULL) free(P[0]);        
        free(P);
    }
    
    if(rhs != NULL) free(rhs);
    if(a != NULL) free(a);

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}


/*********************************************************************************
 * aux routine for exp bootstrapping
 *
 ********************************************************************************/
static int FuncRes(
    double     *res,                      /* (O) return value                   */
    double     *eval,                     /* (I) array of eigenvals             */
    double      x,                        /* (I) trial eigenval                 */
    int         idx,                      /* (I) idx of trial eigenval          */
    double     *p1,                       /* (I) vector 1                       */
    double     *p2,                       /* (I) vector 2                       */
    double     *lambda,                   /* (I) array of lambdas               */
    int         numName)                  /* (I) num of tranchelets             */
{
    static char routine[] = "FuncRes";
    
    int       status  = FAILURE;
    int       i;
    int       n = idx + 1;
    Mat       B;

    MatSetSize(&B,n);

    *res = 0.0;
    
    PtrSet(B.ptr,B.size2,0.0);
    for(i = 0;i < n-1;i++)
    {
        B.ptr[i * n + i] = eval[i];
        B.ptr[i * n + i + 1] = (-eval[i])* ( 1 - lambda[i]);
    }
    B.ptr[B.size2 - 1] = x;
    
    if(MatExp(
           B.ptr,
           6,
           1.0,
           B.ptr,
           n) != SUCCESS)
    {
        goto RETURN;
    }
    
    for(i = 0;i < B.size;i++)
    {
        *res += p1[i]*B.ptr[i*B.size + idx];
    }

    *res -= p2[idx];

    status = SUCCESS;
RETURN:
    
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * default transition matrix kernel from startDate to endDate
 *
 ********************************************************************************/
int  TransMatExpKernel(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves)           /* (I) effective curves               */    
{
    static char  routine[] = "TransMatExpKernel";
    
    int       status  = FAILURE;
    int       i, istart, idx, j;
    double    *p1 = NULL, *p2 = NULL, dtmp1, dtmp2;
    double    *evec = NULL, *ptr = NULL, dT, *eval = NULL;
    double    x1, x2, x, res1, res2, res;
    double    *lhs = NULL;
    Mat       B;
    double    err;
    int       maxIter = 10;
    double    tol = 1e-9;
    double    lambda[MAX_NB];
    int       numName = effCurves->numName;
    
    /* initialize */
    MatSetSize(&B, numName);
    PtrSet(B.ptr,B.size2,0.0);
    
    if(startDate >= endDate)
    {
        return SUCCESS;
    }

    p1 = (double *)malloc(numName * sizeof(double));
    p2 = (double *)malloc(numName * sizeof(double));
    lhs = (double *)malloc(numName * sizeof(double));
    ptr = (double *) malloc(numName * sizeof(double));
    eval = (double *) malloc(numName * sizeof(double));
    evec = (double *) malloc(numName * numName * sizeof(double));

    if(!p1 || !p2 || !lhs || !ptr || !eval || !evec)
    {
        CRXError("%s: Memory allocation failed.\n", routine);
        goto RETURN;
    }
    
    /* generate p's */
    dtmp1 = 0.0;
    for(i = 0;i < numName;i++)
    {
        if(GtoDiscountDate(startDate,
                           effCurves->pCurves[i],
                           GTO_FLAT_FORWARDS,
                           &dtmp2) != SUCCESS)
        {
            goto RETURN;
        }

        p1[i] = dtmp2 - dtmp1;
        dtmp1 = dtmp2;
    }
    dtmp1 = 0.0;
    for(i = 0;i < numName;i++)
    {
        if(GtoDiscountDate(endDate,
                           effCurves->pCurves[i],
                           GTO_FLAT_FORWARDS,
                           &dtmp2) != SUCCESS)
        {
            goto RETURN;
        }

        p2[i] = dtmp2 - dtmp1;
        dtmp1 = dtmp2;
    }

    if(GtoDayCountFraction(startDate, endDate, GTO_ACT_365F,&dT) != SUCCESS)
    {
        goto RETURN;
    }

    
    //CalcLambda(&dtmp,&dtmp1, p1,p2,numName);
    for(i = 0;i < numName;i++)
    {
        lambda[i] = 0.0;
    }

    /* incremently solve */
    eval[0] = log(p2[0]/p1[0]);
    PtrSet(&evec[0],numName,0.0);
    evec[0] = 1.0;

    B.ptr[0] = eval[0];
    B.ptr[1] = (-eval[0]) * (1-lambda[0] );
    B.ptr[numName-1] += (-eval[0])* lambda[0];
    for(idx = 1;idx < numName-1;idx++)
    {
        istart = idx * numName;
        
        x1 = 0.0;
        if(FuncRes(&res1, eval, x1, idx, p1, p2, lambda, numName) != SUCCESS)
        {
            CRXError("%s: FuncRes failed.\n",routine);
            goto RETURN;
        }

        x2 = -1.0;
        if(FuncRes(&res2, eval, x2, idx, p1, p2, lambda, numName) != SUCCESS)
        {
            CRXError("%s: FuncRes failed.\n",routine);
            goto RETURN;
        }
            

        err = 2 * tol;
		x = (x1+x2)/2.0;
        for(j = 0;j < maxIter;j++)
        {
            
            x = x1 - res1 * (x2 - x1)/(res2 - res1);
            if(FuncRes(&res, eval,x, idx,p1, p2, lambda, numName) != SUCCESS)
            {
                CRXError("%s: FuncRes failed.\n",FAILURE);
                goto RETURN;
            }
                
            if(fabs(res) < tol) break;

            x1 = x2;
            res1 = res2;
            x2 = x;
            res2 = res;
        }

        /* update matrix */
        eval[idx] = x;

        B.ptr[istart+idx] = x;
        B.ptr[istart+idx+1] = (-x)*(1 - lambda[idx]);
        
        B.ptr[istart+numName-1] += (-x)*lambda[idx];

    }

    mat->size = B.size;
    mat->size2 = B.size2;

    if(GtoDayCountFraction(startDate, endDate, GTO_ACT_365F,&dT) != SUCCESS)
    {
        goto RETURN;
    }
    
    PtrScale(B.ptr,B.size2,1.0/dT);
    PtrCopy(mat->ptr,B.ptr,B.size2);
    
    status = SUCCESS;
RETURN:
    if(p1 != NULL) free(p1);
    if(p2 != NULL) free(p2);
    if(lhs != NULL) free(lhs);
    if(eval != NULL) free(eval);    
    if(evec != NULL) free(evec);
    if(ptr != NULL) free(ptr);

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}


/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/
int  TransMatExp(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "TransMatExp";
        
    int       status  = FAILURE;
    double    dT;
    Mat       B;
    int       numName = effCurves->numName;
    
    MatSetSize(&B, numName);

    if(GtoDayCountFraction(startDate, endDate,GTO_ACT_365F, &dT) != SUCCESS)
    {
        goto RETURN;
    }

    if(TransMatExpKernel(
           &B,
           startDate,
           endDate,
           effCurves) != SUCCESS)
    {
        goto RETURN;
    }

    PtrScale(B.ptr,B.size2,dT);

    mat->size = B.size;
    mat->size2 = B.size2;
    
    if(MatExp(
           mat->ptr,
           6,
           1.0,
           B.ptr,
           numName) != SUCCESS)
    {
        goto RETURN;
    }
    
    status = SUCCESS;
RETURN:

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/    
static int BKQuickMat(
    Mat             *mat,                 /* (O) mat                            */
    double          *slope,               /* (I) conditional slope              */
    double          kFac,                 /* (I) rotation                       */
    double          decay,                /* (I) drift                          */
    double          alpha,                /* (I) alpha                          */    
    double          *eLoss,               /* (I) array of expected loss         */
    double          *pTime,               /* (I) array of time                  */
    int             numStep,              /* (I) num of timesteps               */
    int             numName,              /* (I) num of tranchelets             */
    int             numTotal,             /* (I) total num of tranchelets       */
    double          *sum1,                /* (I) sum of exponential decay       */
    double          *sum2,                /* (I) sum of exponential decay       */
    double          *tail)                /* (I) tail component                 */
{
    static char      routine[] = "BKQuickMat";
    int       status  = FAILURE;
    int       i, j, k;
    double    drift, dL, c, l, dt, sum, ftd;
    Mat       tmpMat;
    
    MatSetSize(&tmpMat, numName);
    PtrSet(mat->ptr,mat->size2,0.0);


    for(i = 0;i < numName;i++)
    {
        mat->ptr[i*numName + i] = 1.0;
    }
    
    for(i = 0;i < numStep - 1;i++)
    {
        dL = eLoss[i+1] - eLoss[i];
        dt = pTime[i+1] - pTime[i];

        PtrSet(tmpMat.ptr,tmpMat.size2,0.0);
        
        for(j = 0;j < numName;j++)
        {
            l = (j + 1.0)/numTotal;
            
            //drift = (numTotal-j)* (dL + slope[j]*l * dt) /(1.0 + decay * pTime[i]);
            drift = (numTotal-j)* (dL + alpha*l * dt) /(1.0 + alpha * pTime[i]);

            c = slope[j] /(numTotal-j) * drift;
            ftd = c + (1.0 - slope[j]) * drift;

            tmpMat.ptr[j*numName + j] = 1.0 - ftd;

            sum = -ftd;
            for(k = j+1;k < numName;k++)
            {
                tmpMat.ptr[j*numName + k] = ftd * tail[k-j-1];
                sum += tmpMat.ptr[j*numName + k];
            }

            tmpMat.ptr[j*numName+numName-1] += -sum;
            
        }
        tmpMat.ptr[numName*numName-1] = 1.0;
        
        if(MatMultUT(mat, mat, &tmpMat) != SUCCESS)
        {
            goto RETURN;
        }
    }

    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

static int BKFuncRes(
    double          *res,                 /* (O) result                         */
    double          *lossDist,            /* (O) lost dist at end date          */
    double          g,                    /* (I) trial gamma                    */
    int             index1,               /* (I) gamma index                    */
    int             index2,               /* (I) gamma index                    */    
    double          *p1,                  /* (I) loss dist pdf at t1            */
    double          *payoff,              /* (I) payoff vector                  */
    double          loss,                 /* (I) target loss                    */    
    Mat             *mat,                 /* (O) mat                            */
    double          *slope,               /* (I) conditional slope              */
    double          kFac,                 /* (I) rotation                       */
    double          decay,
    double          alpha,                /* (I) alpha                          */
    double          *eLoss,               /* (I) array of expected loss         */
    double          *pTime,               /* (I) array of dt                    */    
    int             numStep,              /* (I) num of timesteps               */
    int             numName,              /* (I) num of tranchelets             */
    int             numTotal,             /* (I) total num of tranchelets       */
    double          *sum1,                /* (I) sum of exponential decay       */
    double          *sum2,                /* (I) sum of exponential decay       */
    double          *tail)                /* (I) tail component                 */    
{
    static char      routine[] = "BKFuncRes";
    int       status  = FAILURE;
    int       i, j;
    double    sum;
    double    drift = 0.0;

    for(i = index1;i < index2;i++)
    {
        slope[i] = g;
    }

    
    if(BKQuickMat(mat,
                  slope,
                  kFac,
                  decay,
                  alpha,
                  eLoss,
                  pTime,
                  numStep,
                  numName,
                  numTotal,
                  sum1,
                  sum2,
                  tail) != SUCCESS)
    {
        goto RETURN;
    }

    *res = 0.0;
    sum = 0.0;
    for(i = 0;i <= index2;i++)
    {
        for(j = 0;j <= index2;j++)
        {
            *res += p1[i] * mat->ptr[i*numName+j] * payoff[j];
        }

    }

    *res /= loss;
    *res -= 1.0;

    for(i = 0;i < numName;i++)
    {
        lossDist[i] = mat->ptr[i];
    }
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/
int  TransMat(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "TransMat";
    int       status  = FAILURE;
    
    if((effCurves->step).prd <= 0)
    {
        CRXError("%s: Timestep <= 0.\n", routine);
        goto RETURN;
    }
    
    status = TransMatCat(
        mat,
        startDate,
        endDate,
        effCurves);
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}
    
/*********************************************************************************
 * calculate Mat * p,M is from startDate to endDate, p is at the endDate
 * it gets the average vector p at the startDate
 *
 ********************************************************************************/
int  MatDiscount(
    double          *price,               /* (I/O) price vector                 */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves,           /* (I) effective curves               */
    TCurve          *irCurve)             /* (I) ir curve                       */    
{
    static char      routine[] = "MatAverage";
   
    int       status  = FAILURE;
    Mat       mat;
    double    priceL[MAX_NB];
    double    discFac;
    int       numName = effCurves->numName;
    
    if(TransMat(&mat,
                startDate,
                endDate,
                effCurves) != SUCCESS)
    {
        goto RETURN;
    }

    if(MatMultVecUT(priceL,&mat,price) != SUCCESS)
    {
        goto RETURN;
    }
    
    PtrCopy(price,priceL,numName);

    if(GtoDiscountDateForward(startDate,
                              endDate,
                              irCurve,
                              GTO_FLAT_FORWARDS,
                              &discFac) != SUCCESS)
    {
        goto RETURN;
    }

    PtrScale(price,numName, discFac);
    
    status = SUCCESS;
RETURN:

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Calibrate DTM to a given loss distritubions on a specific tenor
 * 
 *
 ********************************************************************************/
static int  DTMCalibSingle(
    DTMParm         *dtmParm,             /* (O) DTM param structure            */
    TDate           today,                /* (I) today                          */
    int             index,                /* (I) target mat index               */
    TargetLoss      *target,              /* (I) target tranche structure       */
    TCurve          *portCurve,           /* (I) average spread                 */
    double          rho,                  /* (I) CM correlation                 */
    double          decay,                /* (I) tail decay                     */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "DTMCalibSingle";
    int       status  = FAILURE;
    
    int       i, idx,j, index1,index2, k;
    TDate     date1;
    double    *eLoss = NULL, *pTime = NULL;
    double    survProb, oldSurvProb;
    double    g1, g2, g0, res1, res2, res0;
    double    pSlope[MAX_NB], lossDist0[MAX_NB];
    TDateList *tranDates = NULL;
    int       maxIter;
    double    gStep;
    double    tol = 0.001;
    int       numTStep;
    double    alpha;
    Mat       mat;
    double    tail[MAX_NB],sum1[MAX_NB],sum2[MAX_NB], ek, dtmp;
    double    payoff[MAX_NB][MAX_NB];
    TDate     startDate;
    int       numName = effCurves->numName;
    int       numTotal = effCurves->numTotal;
    
    MatSetSize(&mat, numName);
    PtrSet(mat.ptr,mat.size2,0.0);

    
    maxIter = 20;
    alpha = CalcAlpha(rho, numTotal);

    /* loss dist at start date */
    if(index == 0)                        
    {
        PtrSet(lossDist0,numName, 0.0);
        lossDist0[0] = 1.0;
        startDate = today;
    } else {
        for(i = 0;i < numName;i++)
        {
            lossDist0[i] = target->lossDist[index-1][i];
        }
        startDate = target->maturityDate[index-1];
    }
    
    /* calculate */
    tranDates  = GtoNewDateList(
        startDate,
        target->maturityDate[index],
        &(effCurves->step),
        FALSE);
    
    if(tranDates == NULL)
    {
        CRXError("%s: Construction of date list failed.\n",routine);
        goto RETURN;
    }

    numTStep = tranDates->fNumItems;
    
    eLoss = (double *)malloc(numTStep * sizeof(double));
    pTime = (double *)malloc(numTStep * sizeof(double));    

    if(!eLoss || !pTime)
    {
        CRXError("%s: Memory allocation failed.\n",routine);
        goto RETURN;
    }

    for(i = 0;i < numName;i++)
    {
        tail[i] = 0.0;//debug! pow(rho, i * 
    }


    
    oldSurvProb = 1.0;
    for(i = 0;i < numTStep;i++)
    {
        date1  = (i > 0)?tranDates->fArray[i]:startDate;
        
        if(GtoDayCountFraction(today,date1, GTO_ACT_365F,&pTime[i]) != SUCCESS)
        {
            goto RETURN;
        }
        
        if(GtoDiscountDate(date1,
                           portCurve,
                           GTO_FLAT_FORWARDS,
                           &survProb) != SUCCESS)
        {
            goto RETURN;
        }

        if(i > 0)
        {
            eLoss[i] = eLoss[i-1] + (oldSurvProb - survProb);
        } else {
            eLoss[i] = (oldSurvProb - survProb);
        }
        oldSurvProb = survProb;
    }
    
    /* init */
    PtrSet(pSlope,numName,0.4);

    
    gStep = 0.1;

    /* payoff */
    for(i = 0;i < target->numStrike;i++)
    {
        if(SetTranchePayoff(
               payoff[i],
               0.0,
               target->strikes[i],
               1.0,
               effCurves) != SUCCESS)
        {
            goto RETURN;
        }

        for(j = 0;j < numName;j++)
        {
            payoff[i][j] = 1.0 - payoff[i][j];
        }
    }

    /* calibrate */
    for(j = 0;j < 1;j++)
    {
        decay = 0.0;
        for(i = 0;i < numName;i++)
        {
            decay += pSlope[i];
        }
        decay /= numName;

        for(i = 0;i < target->numStrike;i++)
        {
            idx = 0;
            
            index1 = (i > 0)?target->index[i-1]:0;
            index2 = target->index[i];
            if(index2 > numName) index2 = numName - 1;

            if(index2 < index1) continue; /* not calibrated */
            g1 = (index1 > 0)?pSlope[index1-1]:0.0;
            
            if(BKFuncRes(&res1,
                         &(target->lossDist[index][0]),
                         g1,
                         index1,
                         index2,
                         lossDist0,
                         payoff[i],
                         target->loss[i][index],
                         &mat,
                         pSlope,
                         decay,//debug! kFac
                         decay,
                         alpha,
                         eLoss,
                         pTime,
                         numTStep,
                         numName,
                         numTotal,
                         sum1,
                         sum2,
                         tail) != SUCCESS)
            {
                goto RETURN;
            }
            
            if(fabs(res1) < tol)              /* already converged */
            {
                for(k = index1;k < index2;k++)
                {
                    pSlope[k] = g1;
                }
                continue;
            }
            
            
            g2 = g1 + gStep;
            if(BKFuncRes(&res2,
                         &(target->lossDist[index][0]),
                         g2,
                         index1,
                         index2,
                         lossDist0,
                         payoff[i],
                         target->loss[i][index],
                         &mat,
                         pSlope,
                         decay, //debug! kFac,
                         decay,
                         alpha,                     
                         eLoss,
                         pTime,
                         numTStep,
                         numName,
                         numTotal,
                         sum1,
                         sum2,
                         tail) != SUCCESS)
            {
                goto RETURN;
            }
            
            while(idx < maxIter)
            {
                g0 = g1 - res1 * (g2 - g1)/(res2 - res1);
                
                if(g0 < -1.0) g0 = -1.0;
                if(g0 > 1.0) g0 = 1.0;
                
                if(BKFuncRes(&res0,
                             &(target->lossDist[index][0]),
                             g0,
                             index1,
                             index2,
                             lossDist0,
                             payoff[i],
                             target->loss[i][index],
                             &mat,
                             pSlope,
                             decay, //debug! kFac,
                             decay,
                             alpha,
                             eLoss,
                             pTime,
                             numTStep,
                             numName,
                             numTotal,
                             sum1,
                             sum2,
                             tail) != SUCCESS)
                {
                    goto RETURN;
                }
                
                if(fabs(res0) <= tol)
                {
                    break;
                }
                
                idx++;
                
                g1 = g2;
                res1 = res2;
                g2 = g0;
                res2 = res0;
            }

            for(k = index1;k < index2;k++)
            {
                pSlope[k] = g0;
            }
        }
    }

    /* fill out dtm parm */
    dtmParm->numName    = numName;

    dtmParm->step       = effCurves->step;
    dtmParm->numTotal   = numTotal;
    dtmParm->portCurve  = portCurve;
    dtmParm->baseDate   = today;
    dtmParm->alpha      = alpha;
    dtmParm->kFac       = decay; //debug! kFac;
    dtmParm->decay      = decay;
    dtmParm->numT       = 1;

    for(i = 0;i < target->numDate;i++)
    {
        dtmParm->facDate[i] = target->maturityDate[i];
    }
    
    dtmParm->numDate        = target->numDate;
    
    for(i = 0;i < numName;i++)
    {
        dtmParm->fac[i]    = pSlope[i];
        dtmParm->sum1[i]   = sum1[i];
        dtmParm->sum2[i]   = sum2[i];
        dtmParm->tail[i]   = tail[i];
    }
        
    status = SUCCESS;
RETURN:
    if(tranDates != NULL) GtoFreeDateList(tranDates);

    if(eLoss != NULL)free(eLoss);
    if(pTime != NULL)free(pTime);

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Calibrate DTM to a given loss distritubions on a series of tenors
 * 
 *
 ********************************************************************************/
int  DTMCalib(
    DTMParms        *dtmParm,             /* (O) DTM param structure            */
    TDate           today,                /* (I) today                          */
    TargetLoss      *target,              /* (I) target tranche structure       */
    TCurve          *portCurve,           /* (I) average spread                 */
    double          rho,                  /* (I) CM correlation                 */
    double          decay,                /* (I) tail decay                     */
    EffCurves       *effCurves)           /* (I) effective curves               */    
{
    static char      routine[] = "DTMCalib";
    int       status  = FAILURE;
    int       i;

    for(i = 0;i < target->numDate;i++)
    {
        if(DTMCalibSingle(
               &dtmParm->parms[i],
               today,
               i,
               target,
               portCurve,
               rho,
               decay,
               effCurves) != SUCCESS)
        {
            goto RETURN;
        }
    }
    
    status = SUCCESS;    
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Calc DTM given dtmParm
 * 
 *
 ********************************************************************************/
int  DTMCalc(
    Mat             *mat,                 /* (O) transition matrix              */
    DTMParms        *dtmParms,            /* (I) DTM param structure            */
    TDate           startDate,            /* (I) today                          */
    TDate           endDate)              /* (I) target tranche mat date        */
{
    static char      routine[] = "DTMCalc";
    int       status  = FAILURE;

    int       i, numTStep, index1, index2;
    TDate     date1;
    double    *eLoss = NULL, survProb, oldSurvProb, *pTime = NULL;    
    TDateList *tranDates = NULL;
    DTMParm   *dtmParm = NULL;
    Mat       tmpMat;
    int       numName, numDate;
    TDateInterval step;
    TDate     baseDate;
    TCurve    *portCurve = NULL;
    TDate     *facDate   = NULL;
    
    numName  = dtmParms->parms[0].numName;
    numDate  = dtmParms->parms[0].numDate;
    facDate  = dtmParms->parms[0].facDate;
    step     = dtmParms->parms[0].step;
    baseDate = dtmParms->parms[0].baseDate;
    portCurve = dtmParms->parms[0].portCurve;

    
    MatSetSize(&tmpMat, numName);
    MatSetSize(mat, numName);
    PtrSet(mat->ptr,mat->size2,0.0);
    for(i = 0;i < numName;i++)
    {
        mat->ptr[i*numName + i] = 1.0;
    }
    
    /* calculate */
    tranDates  = GtoNewDateList(
        startDate,
        endDate,
        &step,
        FALSE);
    
    if(tranDates == NULL)
    {
        CRXError("%s: Construction of date list failed.\n",routine);
        goto RETURN;
    }

    numTStep = tranDates->fNumItems;
    
    eLoss = (double *)malloc(numTStep * sizeof(double));
    pTime = (double *)malloc(numTStep * sizeof(double));    

    if(!eLoss || !pTime)
    {
        CRXError("%s: Memory allocation failed.\n",routine);
        goto RETURN;
    }

    oldSurvProb = 1.0;
    for(i = 0;i < numTStep;i++)
    {
        date1  = tranDates->fArray[i];
        
        if(GtoDayCountFraction(baseDate,date1, GTO_ACT_365F,&pTime[i]) != SUCCESS)
        {
            goto RETURN;
        }
        
        if(GtoDiscountDate(date1,
                           portCurve,
                           GTO_FLAT_FORWARDS,
                           &survProb) != SUCCESS)
        {
            goto RETURN;
        }

        if(i == 0)
        {
            eLoss[i] = 1.0 - survProb;
        } else {
            eLoss[i] = eLoss[i-1] + (oldSurvProb - survProb);
        }
        
        oldSurvProb = survProb;
    }

    for(i = 0;i < numDate;i++)
    {
        if(startDate <= facDate[i]) break;
    }


    if(i == numDate)
    {
        dtmParm = &dtmParms->parms[numDate - 1];
        if(BKQuickMat(
               mat,
               dtmParm->fac,
               dtmParm->kFac,
               dtmParm->decay,
               dtmParm->alpha,
               eLoss,
               pTime,
               numTStep,
               dtmParm->numName,
               dtmParm->numTotal,
               dtmParm->sum1,
               dtmParm->sum2,
               dtmParm->tail) != SUCCESS)
        {
            goto RETURN;
        }
        
    } else {
        
        for(index1 = 0;index1 < numTStep;index1++)
        {
            for(index2 = index1; index2 < numTStep;index2++)
            {
                if(tranDates->fArray[index2] <= facDate[i]){
                    continue;
                }
                break;
            }

            if(index2 >= numTStep) index2 = numTStep - 1;


            dtmParm = &dtmParms->parms[i];            
            if(BKQuickMat(
                   &tmpMat,
                   dtmParm->fac,
                   dtmParm->kFac,
                   dtmParm->decay,
                   dtmParm->alpha,
                   &eLoss[index1],
                   &pTime[index1],
                   index2 - index1 + 1,
                   dtmParm->numName,
                   dtmParm->numTotal,
                   dtmParm->sum1,
                   dtmParm->sum2,
                   dtmParm->tail) != SUCCESS)
            {
                goto RETURN;
            }
            
            if(MatMultUT(mat, mat, &tmpMat) != SUCCESS)
            {
                goto RETURN;
            }
            
            
            index1 = index2;
            i++;
            if(i >= numDate) i = numDate-1;
            if(index2 == numTStep-1) break;
        }
    }

    
    status = SUCCESS;
RETURN:
    if(tranDates != NULL) GtoFreeDateList(tranDates);
    if(eLoss != NULL)free(eLoss);
    if(pTime != NULL)free(pTime);

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }
    
    return SUCCESS;
}
