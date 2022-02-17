/*********************************************************************************
 * Base Correlation Smile Calibration using M-Mapping
 *
 ********************************************************************************/

#include "transition.h"
#include "lintrp.h"

#define SQRTPI  1.772453850905516027     /* sqrt(pi) */

/*********************************************************************************
 * Test Func
 *
 ********************************************************************************/
int     CMTest(
    double          *res,                 /* (O) res                            */
    int             n,                    /* (I) n                              */
    int             N,                    /* (I) N                              */    
    double          rho,                  /* (I) rho                            */
    double          prob,                 /* (I) prob                           */
    char            flag)                 /* (I) flag                           */
{
    static char      routine[] = "CMTest";
    int              status  = FAILURE;
    int              i;
    double           y1 = -5, y2 = 5, dy = 0.001, y;
    double           sum, K, sum0;
    double           sr, sr2;
    double           pz, ppz, facy, facn;

    K  = NormCumInv(prob);
    sr = sqrt(rho);
    sr2 = sqrt(1-rho);
    
    y = y1;
    sum = 0.0;
    sum0 = 0.0;
    while(y <= y2)
    {
        pz  = NormCum((K - sr * y)/sr2);
        if(pz >= 1 - SMALL)
        {
            y += dy;
            continue;
        }

        ppz = pow(1-pz,N);
        facy = exp(-y*y/2) * dy;
        facn = 1;

        for(i = 0;i <=n;i++)
        {
            sum += facy * facn * ppz * i;//* (n-i);
            sum0 += facy * facn * ppz;
            ppz = ppz * pz/(1-pz);
            facn = facn * (N-i)/(i+1);
        }
        y += dy;
    }


    //sum /= sum0;
    sum /= (sqrt(2)*SQRTPI);
    sum0 /= (sqrt(2)*SQRTPI);

    if(flag == 'C')
    {
        *res = sum0;
    } else {
        *res = sum;
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
 * Function to compute loss cdf and pdf
 *
 ********************************************************************************/
static double lossCDF(
    double          K1,                   /* (I) K                              */
    double          K2,                   /* (I) K                              */
    BCSmile         *bc);                 /* (I/O) x                            */            


/*********************************************************************************
 * Credit Metrics Large Portfolio (LPHA) approximation PDF
 *
 ********************************************************************************/
static int  LPHAPDF(
    double          *pdf,                 /* (O) cdf                            */
    double          x,                    /* (I) x                              */
    double          beta,                 /* (I) base correlation               */
    double          p)                    /* (I) default probability            */
{
    static char      routine[] = "LPHACDFInv";
    int              status  = FAILURE;
    double           n1, n2;
    
    if(x < 0.0 || x > 1.0)
    {
        CRXError("%s: x should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(p < 0.0 || p > 1.0)
    {
        CRXError("%s: p should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(beta < 0.0 || beta > 1.0)
    {
        CRXError("%s: beta should be > 0 and < 1.", routine);        
        goto RETURN;
    }

    if(beta < SMALL)
    {
        beta = SMALL;
    }         

    n1 = NormCumInv (x);
    n2 = NormCumInv (p);

    *pdf = sqrt((1-beta)/beta) *
        exp(0.5*n1*n1 - 0.5/beta*pow((n2 - sqrt(1-beta)*n1),2.0));
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Credit Metrics Large Portfolio (LPHA) approximation CDF
 *
 ********************************************************************************/
static int  LPHACDF(
    double          *cdf,                 /* (O) cdf                            */
    double          x,                    /* (I) x                              */
    double          beta,                 /* (I) base correlation               */
    double          p)                    /* (I) default probability            */
{
    static char      routine[] = "LPHACDFInv";
    int              status  = FAILURE;
    double           n1, n2;
    

    if(x < 0.0 || x > 1.0)
    {
        CRXError("%s: x should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(p < 0.0 || p > 1.0)
    {
        CRXError("%s: p should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(beta < 0.0 || beta > 1.0)
    {
        CRXError("%s: beta should be > 0 and < 1.", routine);        
        goto RETURN;
    }

    if(beta < SMALL)
    {
        beta = SMALL;
    }         

    n1 = NormCumInv (x);
    n2 = NormCumInv (p);
    
    *cdf = NormCum((sqrt(1-beta) * n1 - n2)/sqrt(beta));
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

static int LPHAELoss(
    double          *loss,                /* (O) expected loss                  */
    double          K1,                   /* (I) K                              */
    double          K2,                   /* (I) K                              */
    double          beta,                 /* (I) base correlation               */
    double          p)                    /* (I) default probability            */
{
    static char      routine[] = "LPHAlossExpected";
    int              status  = FAILURE;
    
    double    sum = 0.0, sum2 = 0.0, dtmp;
    double    step = 0.05;
    double    K;
    double    k1, k2, np, sbeta;

    if(K1 < SMALL) K1 = SMALL;
    if(K2 < SMALL) K2 = SMALL;
    if(K1 > 1 - SMALL) K1 = 1-SMALL;
    if(K2 > 1 - SMALL) K2 = 1-SMALL;
    
    k1 = NormCumInv(K1);
    k2 = NormCumInv(K2);    

    np = NormCumInv(p);
    sbeta = sqrt(1-beta);
    
    K = NormCumInv(SMALL);
    while(K < NormCumInv(1-SMALL))
    {
        
        dtmp = (K2-NormCum(K)) * exp(- 0.5/beta * pow(np - sbeta * K,2.0));
        
        if(K < k2 && K >= k1)
        {
            sum += dtmp;
        }

        sum2 += dtmp;
        K   += step;
    }

    //sum2 = sum2 * sqrt(1/beta-1)/sqrt(2)/SQRTPI * step;

    *loss = sum * sqrt(1/beta-1)/sqrt(2)/SQRTPI * step;
    //*loss = sum/sum2 * p;
        

    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Credit Metrics Large Portfolio (LPHA) approximation CDF Inv
 *
 ********************************************************************************/
static int  LPHACDFInv(
    double          *x,                   /* (O) x                              */
    double          cdf,                  /* (I) cdf                            */
    double          beta,                 /* (I) base correlation               */
    double          p)                    /* (I) default probability            */
{
    static char      routine[] = "LPHACDFInv";
    int              status  = FAILURE;
    double           n2;
    
    if(p < 0.0 || p > 1.0)
    {
        CRXError("%s: p should be > 0 and < 1.", routine);
        goto RETURN;
    }
    
    if(cdf < 0.0 || cdf > 1.0)
    {
        CRXError("%s: cdf should be > 0 and < 1.", routine);
        goto RETURN;
    }
    
    if(beta < 0.0 || beta > 1.0)
    {
        CRXError("%s: beta should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(beta > 1 - SMALL)
    {
        beta = SMALL;
    }         

    n2 = NormCumInv (p);
    
    *x = (NormCumInv(cdf) * sqrt(beta) + n2)/sqrt(1-beta);
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * auxi func.
 *
 ********************************************************************************/
static double getSum(
    int             num,                  /* (I) index                          */
    double          *pX,                  /* (I) strike array 1                 */
    double          *b1,                  /* (I) number array 1                 */
    double          *b2,                  /* (I) number array 2                 */
    double          *b,                   /* (I) number array 2                 */
    int             numGrid)              /* (I) num of interpolation grid      */
{
    double       sum;
    int          i;
    int          mask = 0;
    double       x2, x, xy, y, dtmp, k1, k2;

    mask |= 0x0001;


    for(i = 0;i < numGrid;i++)
    {
        if((mask & num) == 1)
        {
            b[i] = b1[i];
        } else {
            b[i] = b2[i];
        }

        num >>= 1;
    }

    x2 = 0.0;
    x  = 0.0;
    xy = 0.0;
    y  = 0.0;
    for(i = 0;i < numGrid;i++)
    {
        dtmp = pX[i];

        x  += dtmp;
        x2 += dtmp * dtmp;
        y  += b[i];
        xy += dtmp * b[i];
    }

    k1 = (numGrid * xy - x * y)/(x2 *numGrid - x *x);
    k2 = (x2 * y - x * xy)/(x2 *numGrid - x *x);
    
    sum = 0;    
    for(i = 0;i < numGrid;i++)
    {
        dtmp = pX[i];
        sum +=  pow(b[i] - k1 * dtmp - k2,2.0);;
    }
    
    return sum;
}

/*********************************************************************************
 * Credit Metrics Large Portfolio (LPHA) approximation CDF
 *
 ********************************************************************************/
static int  LPHACDF2Corr(
    double          *beta,                /* (I) base correlation               */
    double          cdf,                  /* (I) cdf                            */
    double          x,                    /* (I) x                              */
    double          p,                    /* (I) default probability            */
    BCSmile         *bc,                  /* (I) x                              */
    char            flag)                 /* (I) flag                           */
{
    static char      routine[] = "LPHACDF2Corr";
    int              status  = FAILURE;
    double           np, nc, nx,  a, b, c, delta;
    double           b1, b2, pb1[MAX_NB], pb2[MAX_NB], pb3[MAX_NB], pd[MAX_NB];
    double           t1, t2, t, tmpX, tmpCDF, pX[MAX_NB];
    int              i, numGrid = 8, metric, maxi, maxmetric, idx, idx0;
    double           sum, maxsum;
    
    if(cdf < 0.0 || cdf > 1.0)
    {
        CRXError("%s: cdf should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(x < 0.0 || x > 1.0)
    {
        CRXError("%s: x should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(p < 0.0 || p > 1.0)
    {
        CRXError("%s: p should be > 0 and < 1.", routine);
        goto RETURN;
    }

    t = bc->beta[0];
    for(i = bc->numStrike-1;i >= 0;i--)
    {
        if(bc->strike[i] <= x)
        {
            t = bc->beta[i];
            break;
        }
    }
    t = sqrt(t);
    

    np = NormCumInv(p);

    i = 0;
    idx = 0;
    idx0 = 0;
    while(i < numGrid)
    {
        tmpX = (i+1.0)/(numGrid+1);

        if(fabs(tmpX-x) < SMALL)
        {
            idx0 = idx;
        } else {
            
            if(idx == 0 && tmpX > x)
            {
                pX[idx] = x;
                idx0 = idx;
                idx++;
                continue;
            } else if(tmpX < x && tmpX + 1.0/(numGrid+1) > x)
            {
                pX[idx] = tmpX;
                idx++;
                
                pX[idx] = x;
                idx0 = idx;            
                idx++;
                i++;

                continue;
            }
        }
        
        pX[idx] = tmpX;
        i++;
        idx++;

    }

    numGrid += (idx - i);
    
    for(i = 0;i < numGrid;i++)
    {
        tmpX = pX[i];
        tmpCDF = lossCDF(0,tmpX,bc);
        if(tmpCDF > 1.0) tmpCDF = 1.0-SMALL;
        nc = NormCumInv(tmpCDF);        
        nx = NormCumInv(tmpX);           
        a = nx * nx + nc * nc;
        b = 2 * np * nc;
        c = np * np - nx * nx;
        
        delta = b * b - 4 * a * c;
        pd[i] = delta;
        
        if(delta < 0.0)                       /* minimize */
        {
            if(fabs(-b/2/a) + fabs(1+b/2/a) <= 1.0)
            {
                *beta = - b/2/a;            
            } else {
                if(-b/2/a < 0.0)
                {
                    *beta = 0.0;
                } else {
                    *beta = 1.0;
                }
            }

            pb1[i] = *beta;
            pb2[i] = *beta;
        } else {
            if(a < SMALL)
            {
                if(fabs(np) < SMALL)
                {
                    b1 = t;
                    b2 = t;
                } else {
                    b1 = 1.0;
                    b2 = 1.0;
                }
            } else {
                b1 = (-b + sqrt(delta))/2/a;
                b2 = (-b - sqrt(delta))/2/a;
            }
            
            t1 = fabs(b1) + fabs(1-b1);
            t2 = fabs(b2) + fabs(1-b2);
            
            
            if(t1 <= 1.0)
            {
                if(t2 <= 1.0)
                {
                    if(fabs(b1 - t)  < fabs(b2 - t))
                    {
                        *beta = b1;
                    } else {
                        *beta = b2;
                    }

                    pb1[i] = b1;
                    pb2[i] = b2;
                } else {
                    *beta = b1;

                    pb1[i] = b1;
                    pb2[i] = b1;
                }
                
                
            } else {
                if(t2 <= 1.0)
                {
                    *beta = b2;

                    pb1[i] = b2;
                    pb2[i] = b2;
                    
                } else {
                    *beta = t;//to be added
                }
            }
        }
    }

    metric = 1;
    maxi = -1;
    maxsum = -1.0;
    maxmetric = -1;
    for(i = 0;i < pow(2,numGrid);i++)
    {
        sum = getSum(metric,pX,pb1,pb2,pb3,numGrid);
        
        if(sum < maxsum || maxsum < 0.0)
        {
            maxsum = sum;
            maxi   = i;
            maxmetric = metric;
        }
        
        metric += 1;
    }

    sum = getSum(maxmetric,pX, pb1,pb2,pb3, numGrid);

    if(flag == 'U')
    {
        *beta = (pb1[idx0] < pb2[idx0])? pb2[idx0]:pb1[idx0];
    } else if(flag == 'L')
    {
        *beta = (pb1[idx0] < pb2[idx0])? pb1[idx0]:pb2[idx0];
    } else if(flag == 'D')
    {
        *beta = pd[idx0];
    } else {
        *beta = pb3[idx0];
    }
    
    *beta *= (*beta);
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Credit Metrics Large Portfolio (LPHA) approximation CDF
 *
 ********************************************************************************/
static int  LPHAELoss2Corr(
    double          *beta,                /* (I) base correlation               */
    double          loss,                 /* (I) loss                           */
    double          x,                    /* (I) x                              */
    double          p,                    /* (I) default probability            */
    BCSmile         *bc)                  /* (I) x                              */
{
    static char      routine[] = "LPHAELoss2Corr";
    int              status  = FAILURE;
    double           b1, b2, b;
    double           t1, t2, t;
    double           err, errTOL = 1e-8, idx, maxIter = 10;
    
    if(loss < 0.0 || loss > 1-bc->prob)
    {
        CRXError("%s: put should be > 0 and < %f.", routine, 1-bc->prob);
        goto RETURN;
    }

    if(x < 0.0 || x > 1.0)
    {
        CRXError("%s: x should be > 0 and < 1.", routine);
        goto RETURN;
    }

    if(p < 0.0 || p > 1.0)
    {
        CRXError("%s: p should be > 0 and < 1.", routine);
        goto RETURN;
    }

    b1 = 0.4;
    if(LPHAELoss(
           &t1,
           0,
           x,
           b1,
           p) != SUCCESS)
    {
        goto RETURN;
    }
    
    b2 = b1 * 1.1;
    if(LPHAELoss(
           &t2,
           0,
           x,
           b2,
           p) != SUCCESS)
    {
        goto RETURN;
    }

    b = b2;
    err = fabs(t2 - loss);
    idx = 0;
    while(err > p * errTOL && idx < maxIter)
    {
        if(fabs(t2 - t1) < SMALL) break;
        
        b = b1 - (b2 - b1)/(t2 - t1) * (t1 - loss);

        if(b < 0.0) b = b2/2;
        if(b > 1.0) b = b2 + (1-b2)/2;
        
        if(LPHAELoss(
               &t,
               0,
               x,
               b,
               p) != SUCCESS)
        {
            goto RETURN;
        }

        err = fabs(t - loss);
        if(err < p * errTOL) break;

        b1 = b2;
        t1 = t2;
        b2 = b;
        t2 = t;
        
        idx++;
    }
    
    *beta = b;
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}
    
static double lossPDF(
    double          K,                    /* (I) K                              */
    BCSmile         *bc)                  /* (I/O) x                            */        
{
    double    pdf;
    int       i = 0;
    double    rho = bc->rho;
    double    np = NormCumInv(bc->p);
    double    ns;
    
    pdf = 1.0/sqrt(2)/SQRTPI*sqrt((1-rho)/rho)*exp(-0.5*pow(sqrt(1-rho)*K-np,2.0)/rho);
    
    while(i < bc->numStrike)
    {
        if(i == 0)
        {
            pdf *= exp(bc->A[i]*K);
        } else {
            ns = bc->normStrike[i-1];
        
            if(ns < K)
            {
                pdf *= exp(bc->A[i]*(K - ns));
            } else {
                break;
            }
        }
        
        i++;
    }

    pdf *= bc->C;
    
    return pdf;
}

/*********************************************************************************
 * Function to compute loss cdf and pdf
 *
 ********************************************************************************/
static double lossCDFOld(
    double          K1,                   /* (I) K                              */
    double          K2,                   /* (I) K                              */
    BCSmile         *bc)                  /* (I/O) x                            */            
{

    double    sum = 0.0, sum2 = 0.0, dtmp;
    double    step = 0.05;
    double    K,cK;
    double          k1, k2;

    if(K1 < SMALL) K1 = SMALL;
    if(K2 < SMALL) K2 = SMALL;
    if(K1 > 1 - SMALL) K1 = 1-SMALL;
    if(K2 > 1 - SMALL) K2 = 1-SMALL;
    
    k1 = NormCumInv(K1);
    k2 = NormCumInv(K2);    

    K = NormCumInv(SMALL);
    while(K < NormCumInv(1-SMALL))
    {
        cK   = NormCum(K);
        dtmp = lossPDF(K,bc);
        
        if(K < k2 && K >= k1)
        {
            sum += dtmp;
        }

        sum2 += dtmp;
        K   += step;
    }
    
    sum = sum * step;///sum2 * bc->prob;
    
    return sum;
}

static double lossCDF(
    double          K1,                   /* (I) K                              */
    double          K2,                   /* (I) K                              */
    BCSmile         *bc)                  /* (I/O) x                            */            
{
    double          cdf;
    int             i;
    double          step = 0.01;
    double          k1, k2;
    double          a, b, c, u1, u2, y1, y2, oy2, m;
    double          rho = bc->rho;
    double          np = NormCumInv(bc->p);
    double          small = NormCumInv(SMALL);
    double          large = NormCumInv(1-SMALL);
    double          fac1, fac2;

    if(rho < SMALL) rho = SMALL;
    
    if(K1 < SMALL) K1 = SMALL;
    if(K2 < SMALL) K2 = SMALL;

    if(K1 > 1 - SMALL) K1 = 1-SMALL;
    if(K2 > 1 - SMALL) K2 = 1-SMALL;
    
    k1 = NormCumInv(K1);
    k2 = NormCumInv(K2);    
    
    cdf = 0.0;

    y1 = 0.0;
    y2 = 0.0; oy2 = 0.0;
    a = 0.0;  b   = 0.0;
    fac1 = sqrt(1-rho);
    fac2 = sqrt((1-rho)/rho);
    for(i = 0;i < bc->numStrike;i++)
    {
        a   += bc->A[i];
        if(i == 0)
        {
            b = 0;
        } else {
            b   -= bc->A[i] * bc->normStrike[i-1];
        }
        
        if(i == 0)
        {
            y1 = NormCumInv(SMALL);
        } else {
            y1 = oy2;
        }
        
        if(i < bc->numStrike-1)
        {
            y2 = bc->normStrike[i];
        } else {
            y2 = NormCumInv(1.0-SMALL);
        }

        if(y2 > k2) y2 = k2;

        if(y2 > k1)
        {
            if(y1 < k1) y1 = k1;
            
            m = rho/(1-rho)*a + np/fac1;
            c = -(1-rho)/2/rho*(np*np/(1-rho)-m*m);
            
            u1 = fac2 * (y1 - m);
            u2 = fac2 * (y2 - m);

            if(y1 <= small) u1 = small;
            if(y1 >= large) u1 = large;
            if(y2 <= small) u2 = small;
            if(y2 >= large) u2 = large;
            
            cdf += (NormCum(u2) - NormCum(u1))*exp(b+c);
            
        }

        oy2 = y2;
        if(y2 == k2) break;
    }

    cdf *= bc->C;    
    return cdf;
}

static double lossExpected(
    double          K1,                   /* (I) K                              */
    double          K2,                   /* (I) K                              */
    BCSmile         *bc)                  /* (I/O) x                            */    
{
    double    sum = 0.0, sum2 = 0.0, dtmp;
    double    step = 0.05;
    double    K,cK;
    double          k1, k2;

    if(K1 < SMALL) K1 = SMALL;
    if(K2 < SMALL) K2 = SMALL;
    if(K1 > 1 - SMALL) K1 = 1-SMALL;
    if(K2 > 1 - SMALL) K2 = 1-SMALL;
    
    k1 = NormCumInv(K1);
    k2 = NormCumInv(K2);    

    K = NormCumInv(SMALL);
    while(K < NormCumInv(1-SMALL))
    {
        cK   = NormCum(K);
        dtmp = lossPDF(K,bc) * (K2-cK);
        
        if(K < k2 && K >= k1)
        {
            sum += dtmp;
        }

        sum2 += dtmp;
        K   += step;
    }
    
    sum = sum * step;///sum2 * bc->prob;
    
    return sum;
}

/*********************************************************************************
 * Base Corr Smile Calibration Util
 *
 ********************************************************************************/
static int    calcTarget(
    double          *y,                   /* (O) y                              */
    int             idx,                  /* (I) strike idx                     */
    double          target,               /* (I) target                         */
    BCSmile         *bc)                  /* (I/O) x                            */
{    
    static char      routine[] = "calcTarget";
    int              status  = FAILURE;

    *y        = lossExpected(0,bc->K2,bc);
    *y        -= target;
    
    status = SUCCESS;
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * Base Corr Smile Calibration
 *
 ********************************************************************************/
//#define MCHANGE(x) (exp(-x))
//#define MUNCHANGE(x) (-log(x))

#define MCHANGE(x) (x)
#define MUNCHANGE(x) (x)

int      BCSmileCalib(
    BCSmile          *bc)                 /* (I/O) x                            */
{
    static char      routine[] = "BCSmileCalib";
    int              status = FAILURE;
    int              i, j, k, idx;
    double           target;
    int              maxPStep = 20, maxCStep = 10, maxIter = 10;
    double           err, errTOL = 1e-6;
    double           M1, y1, M2, y2, M, y;
    double           p[MAX_NB], eLoss[MAX_NB], C[MAX_NB], E[MAX_NB], oldSum[MAX_NB];
    int              BiSec, i1, i2;
    double           norm;
    double           dtmp;

    /* convert to target losses */
    for(i = 0;i < bc->numStrike;i++)
    {

        
        if(LPHAELoss(&bc->cdf[i],
                     0.0,
                     bc->strike[i],
                     bc->beta[i],
                     bc->prob) != SUCCESS)
        {
            goto RETURN;
        }

        bc->normStrike[i] = NormCumInv(bc->strike[i]);
        
        if(i > 0)
        {
            if(bc->cdf[i] < bc->cdf[i-1])
            {
                CRXError("%s: CDF not increasing.", routine);
                goto RETURN;
            }
        }
    }

    norm = 1.0;
    
    /* more than one input */
    bc->A[0] = 1.0;
    bc->A[1] = 0.0;
    bc->C    = 1.0;
    bc->rho  = bc->beta[0];

    for(k = 0;k < bc->numStrike;k++) /* loop on strikes */
    {
        bc->A[k] = 0.0;
    }
    
    BiSec = 0;
    for(i = 0;i < maxCStep;i++)            /* loop on A[0] */
    {
        C[i] = 1.0;
        E[i] = 1.0;
        
        if(i == 0) bc->p = bc->prob;

        /* init  */
        for(j = 0;j < maxPStep;j++)
        {
            p[j] = bc->prob;
            eLoss[j] = 0.0;
            oldSum[j] = 0.0;
        }
        
        for(j = 0;j < maxPStep;j++)        /* loop on p    */
        {
            
            for(k = 0;k < bc->numStrike;k++) /* loop on strikes */
            {

                target = bc->cdf[k];
                
                bc->K2 = bc->strike[k];
                bc->K1 = (k > 0)? bc->strike[k-1]:0.0;
                
                bc->A[k] = 0.0;
                if(k > 0)bc->A[k] = bc->A[k-1];
                M1 = bc->A[k];
                if(calcTarget(
                       &y1,
                       k,
                       target,
                       bc) != SUCCESS)
                {
                    goto RETURN;
                }

                bc->A[k] = 1.0;
                M2 = bc->A[k];
                if(calcTarget(
                       &y2,
                       k,
                       target,
                       bc) != SUCCESS)
                {
                    goto RETURN;
                }
                
                err = fabs(y2);
                idx = 0;
                M = M2;
                while(err > errTOL && idx < maxIter)
                {

                    if(fabs(y2 - y1) < SMALL) break;
                    
                    M = M1 - (M2 - M1)/(y2 - y1)*y1;
                    if(M < -10.0 || M > 10.0) M = M2/1.1;
                    
                    bc->A[k] = M;
                    if(calcTarget(
                           &y,
                           k,
                           target,
                           bc) != SUCCESS)
                    {
                        goto RETURN;
                    }

                    err = fabs(y);
                    if(err < errTOL) break;

                    M1 = M2;
                    y1 = y2;
                    M2 = M;
                    y2 = y;
                    
                    idx++;
                }
                
                bc->A[k]       = M;
                
            }

            eLoss[j] = exp(-(lossExpected(0,1.0,bc)/(1-bc->prob)-1.0));
            
            /* update p to match expected loss*/
            if(fabs(eLoss[j] - 1.0) < 1e-2) break;

            p[j] = bc->rho;

            if(j == 0)
            {
                bc->rho *= 1.01;
            } else {
                if(fabs(eLoss[j]/eLoss[j-1] - 1.0) < SMALL) break;
                
                bc->rho = p[j-1] - (p[j] - p[j-1])/(eLoss[j]-eLoss[j-1]) *
                    (eLoss[j-1] - 1.0);

            }
            if(bc->rho > 1.0 || bc->rho < 0.0) bc->rho = p[j] * 0.5;
                            
        }

        /* update A[0] to match normality */
        C[i] = bc->C;
        E[i] = lossCDF(0,1.0,bc);
        
        if(fabs(E[i] - 1.0) < errTOL || i == maxCStep - 1) break;

        if(i == 0)
        {
            bc->C *= 1.01;
        } else {
            if(fabs(E[i]/E[i-1]-1.0) < SMALL) break;

            if(BiSec == 0 && (E[i] - 1.0) * (E[i-1] - 1.0) < 0.0)
            {
                BiSec = 1;
                i1 = i-1;
                i2 = i-2;
            }

            if(BiSec == 0 || BiSec == 1)
            {
                bc->C = C[i-1]-(C[i] - C[i-1])/(E[i] - E[i-1]) * (E[i-1] - 1.0);
            } else {
                if((E[i] - 1.0) * (E[i1] - 1.0) < 0.0)
                {
                    bc->C = 0.5 * (C[i1] + C[i]);
                    i2 = i;
                } else {
                    bc->C = 0.5 * (C[i2] + C[i]);
                    i1 = i;
                }
            }
        }


        if(bc->C < 0.1) bc->C = C[i] * 0.99;
        if(bc->C > 5.0) bc->C = C[i] * 1.01;

        dtmp = lossCDF(0,1.0,bc);
    }
    
    bc->calib = 1;
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}


/*********************************************************************************
 * Base Corr Smile Interoplation
 *
 ********************************************************************************/
int      BCSmileInterp(
    double           *res,                /* (O) interpolated base corr         */
    double           strike,              /* (I) strike                         */
    BCSmile          *bc,                 /* (I) smile structure                */
    char             flag)                /* (I) 'B': return beta, 'C':  CDF
                                             'P': PDF                           */
{
    static char      routine[] = "BCSmileInterp";
    int              status  = FAILURE;
    double           sum;
    
    if(bc->calib == 0)
    {
        CRXError("%s: Smile not calibrated.", routine);
        goto RETURN;
    }
    
    /* expect loss -> corr */
    if(flag == 'C')
    {
        *res = lossCDF(0,strike,bc);
        if(*res > 1.0) *res /= lossCDF(0,1.0,bc);
    } else if(flag == 'P')
    {
        *res = lossPDF(strike,bc);
    } else {
        sum = lossExpected(0,strike,bc);
        if(sum > (1-bc->prob)) sum/= (lossExpected(0,1.0,bc)/(1-bc->prob));
        if(LPHAELoss2Corr(res,
                        sum,
                        strike,
                        bc->prob,
                        bc) != SUCCESS)
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
