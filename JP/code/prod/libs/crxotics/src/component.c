/*********************************************************************************
 * pricers based on transition matrix
 *
 ********************************************************************************/
#include "component.h"

/*********************************************************************************
 * Set Tranche Payoff Vector
 *
 ********************************************************************************/
int SetTranchePayoff(
    double          *payoff,              /* (O) tranche contingent payoff vec  */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    double          notional,             /* (I) notional of the payoff         */
    EffCurves       *effCurves)           /* (I) effective curves               */
{
    static char      routine[] = "SetTranchePayoff";
    
    int              status  = FAILURE;
    int              i;
    double           dtmp;
    double           recovery = effCurves->recovery;
    int              numName = effCurves->numName;    
    double           size = (1.0 - recovery > upperSize - lowerSize)?upperSize - lowerSize:1.0-recovery;

    
    if(size < SMALL)
    {
        CRXError("%s: Tranche size < 0.0.\n", routine);
        goto RETURN;
    }
    
    size /= notional;

    for(i =0;i < numName;i++)
    {
        payoff[i] = 0.0;
        
        dtmp = effCurves->loss[i];
        if(dtmp > lowerSize)
        {
            payoff[i] = (dtmp - lowerSize)/size;
        }
        if(dtmp > upperSize)
        {
            payoff[i] = notional;
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
 * protection leg
 * PAY_AS_YOU_GO
 ********************************************************************************/
int MatrixProtection(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           startDate,            /* (I) protection start date          */
    TDate           endDate,              /* (I) protection end date            */
    double          *prot,                /* (I) protection payoff vector       */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve)             /* (I) ir curve                       */
{
    static char      routine[] = "MatrixProtection";
    
    int       status  = FAILURE;
    int       i,j;
    Mat       mat;
    double    pPrice[MAX_NB], tmpPayoff[MAX_NB];
    TDate     date1,date2;
    double    discFac, discFac2;
    TDateList *tranDates = NULL;
    int       numName = effCurves->numName;
    
    if(startDate < today) startDate = today;
    
    /* calculate */
    tranDates  = GtoNewDateList(
        startDate,
        endDate,
        &(effCurves->step),
        FALSE);
    
    if(tranDates == NULL)
    {
        CRXError("%s: Date list construction failed.\n", routine);
        goto RETURN;
    }

    PtrSet(pPrice,numName,0.0);
    PtrSet(tmpPayoff,numName,0.0);
    
    for(i = tranDates->fNumItems-1;i > 0;i--)
    {
        date1  = tranDates->fArray[i-1];
        date2  = tranDates->fArray[i];

        if(date1 < startDate)
        {
            date1 = startDate;
        }

        
        if(TransMat(&mat,
                    date1,
                    date2,
                    effCurves) != SUCCESS)
        {
            goto RETURN;
        }

        if(GtoDiscountDateForward(date1,
                                  date2,
                                  irCurve,
                                  GTO_FLAT_FORWARDS,
                                  &discFac) != SUCCESS)
        {
            goto RETURN;
        }
        if(GtoDiscountDateForward(date1,
                                  date2,
                                  irCurve,
                                  GTO_FLAT_FORWARDS,
                                  &discFac2) != SUCCESS)
        {
            goto RETURN;
        }
        discFac2 = discFac;

        for(j = 0;j < numName;j++)
        {
            tmpPayoff[j] = discFac * pPrice[j] + discFac2 * prot[j];
        }
        
        if(MatMultVecUT(pPrice,&mat,tmpPayoff) != SUCCESS)
        {
            goto RETURN;
        }

        for(j = 0;j < numName;j++)
        {
            pPrice[j] = pPrice[j] - discFac2 * prot[j];
        }
        

        if(date1 == startDate)break;
    }

    if(GtoDiscountDateForward(startDate,
                              valueDate,
                              irCurve,
                              GTO_FLAT_FORWARDS,
                              &discFac) != SUCCESS)
    {
        goto RETURN;
    }
    
    for(i = 0;i < numName;i++)
    {
        //price[i] = (pPrice[i] + prot[i])/discFac; /* include loss at today */
        price[i] = pPrice[i]/discFac; /* include loss at today */        
    }
    
    
    status = SUCCESS;
RETURN:
    if (tranDates != NULL) GtoFreeDateList(tranDates);

    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}

/*********************************************************************************
 * annuity pricer based on transition matrix
 *
 ********************************************************************************/
int  MatrixAnnuity(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) issueDate                      */
    double          cpnRate,              /* (I) coupon rate                    */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    double          notional,             /* (I) notional                       */
    int             payPrincipal,         /* (I) =1, pay principal at maturity  */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve)             /* (I) ir curve                       */
{
    static char         routine[] = "MatrixAnnuity";
    
    int                 status  = FAILURE;
    int                 i, j;
    TDateList           *tranDates = NULL;
    PayOff              *pPayoff = NULL;
    double              dtmp, dt;
    double              size;
    double              recovery = effCurves->recovery;
    int                 numName = effCurves->numName;

    size = (1.0 - recovery > upperSize - lowerSize)?upperSize - lowerSize:1.0-recovery;    
    if(issueDate < today) issueDate = today;
    
    /* calculate */
    tranDates  = GtoNewDateList(
        issueDate,
        maturityDate,
        &cpnInterval,
        FALSE);
    
    if(tranDates == NULL)
    {
        CRXError("%s: Date list construction failed.\n", routine);
        goto RETURN;
    }

    pPayoff = (PayOff *) calloc(1, sizeof(PayOff));
    if(pPayoff == NULL)
    {
        CRXError("%s: Memory allocation for PayOff failed.\n", routine);
        goto RETURN;
    }
    
    pPayoff->numDate = tranDates->fNumItems - 1;
    
    for(i = 1;i < tranDates->fNumItems;i++)
    {
        pPayoff->payDate[i-1] = tranDates->fArray[i];

        if(GtoDayCountFraction(tranDates->fArray[i-1],
                               tranDates->fArray[i],
                               DCC, &dt) != SUCCESS)
        {
            goto RETURN;
        }


        for(j = 0;j < numName;j++)
        {
            dtmp = effCurves->loss[j];

            dtmp = (upperSize - dtmp)/size;

            if(dtmp > 1.0) dtmp = 1.0;
            if(dtmp < 0.0) dtmp = 0.0;
            pPayoff->payoff[i-1][j] = cpnRate * dtmp * dt * notional;

        }
    }
    pPayoff->cutoffDate = pPayoff->payDate[pPayoff->numDate - 1];

    if(MatrixPrice(price,
                   valueDate,
                   today,
                   pPayoff,
                   effCurves,
                   irCurve) != SUCCESS)
    {
        goto RETURN;
    }
                   

    status = SUCCESS;
RETURN:
    if (tranDates != NULL) GtoFreeDateList(tranDates);
    if (pPayoff != NULL) free(pPayoff);
    
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * fwd par spread distribution based on transition matrix
 *
 ********************************************************************************/
int  FwdSpdDistribution(
    MatState        *res,                 /* (O) mat state                      */
    TDate           fwdDate,              /* (I) fwd date                       */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) issueDate                      */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve)             /* (I) ir curve                       */
{
    static char      routine[] = "FwdSpdDistribution";

    int       status  = FAILURE;
    Mat       mat;
    int       numName = effCurves->numName;

    if(MatrixParSpread(
           res,
           fwdDate,
           today,
           issueDate,
           maturityDate,
           cpnInterval,
           DCC,
           lowerSize,
           upperSize,
           effCurves,
           irCurve) != SUCCESS)
    {
        goto RETURN;
    }

    
    if(TransMat(&mat,
                today,
                fwdDate,
                effCurves) != SUCCESS)
    {
        goto RETURN;
    }

    PtrCopy(res->prob,mat.ptr,numName);
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }

    return status;
}

/*********************************************************************************
 * fwd par spread based on transition matrix
 *
 ********************************************************************************/
int  MatrixParSpread(
    MatState        *res,                 /* (O) mat state                      */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) maturityDate                   */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve)             /* (I) ir curve                       */
{
    static char      routine[] = "MatrixParSpread";
        
    int       status  = FAILURE;
    int       i;
    double    annuity[MAX_NB], prot[MAX_NB], spd[MAX_NB];
    double    payoff[MAX_NB];
    Mat       mat;
    int       numName = effCurves->numName;
    
    if(today >= maturityDate){
        CRXError("%s: Valuedate >= maturity date.\n",routine);
        goto RETURN;
    }
    if(issueDate < valueDate)
    {
        issueDate = valueDate;
    }

    if(MatrixAnnuity(
           annuity,
           valueDate,
           valueDate,
           issueDate,
           maturityDate,
           1.0,
           cpnInterval,
           DCC,
           lowerSize,
           upperSize,
           1.0,
           0,
           effCurves,
           irCurve) != SUCCESS)
    {
        goto RETURN;
    }

    if(SetTranchePayoff(
           payoff,
           lowerSize,
           upperSize,
           1.0,
           effCurves) != SUCCESS)
    {
        goto RETURN;
    }

    if(MatrixProtection(
           prot,
           valueDate,
           valueDate,
           issueDate,
           maturityDate,
           payoff,
           effCurves,
           irCurve) != SUCCESS)
    {
        goto RETURN;
    }

    for(i = 0;i < numName;i++)
    {
        if(annuity[i] < SMALL)
        {
            spd[i] = 0.0;
            //spd[i] = prot[i] / annuity[i];            
        } else {
            spd[i] = prot[i] / annuity[i];
        }
    }

    if(TransMat(&mat,
                today,
                valueDate,
                effCurves) != SUCCESS)
    {
        goto RETURN;
    }

    PtrCopy(res->prob,mat.ptr,numName);
    PtrCopy(res->annuity,annuity,numName);
    PtrCopy(res->prot,prot,numName);
    PtrCopy(res->spd,spd,numName);
    
    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}

/*********************************************************************************
 * pricer based on transition matrix
 *
 ********************************************************************************/
int  MatrixPrice(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    PayOff          *pPayoff,             /* (I) payoff matrix                  */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve)             /* (I) ir curve                       */
{
    static char      routine[] = "MatrixPrice";
        
    int       status  = FAILURE;
    int       i,j;
    Mat       mat;
    double    pPrice[MAX_NB], tmpPayoff[MAX_NB];
    TDate     startDate, endDate;
    double    discFac;
    int       numName = effCurves->numName;
    
    if(pPayoff->numDate < 1){
        CRXError("%s: Num. of payment dates < 1.\n",routine);
        goto RETURN;
    }
    
    if(pPayoff->payDate[pPayoff->numDate-1] < today)
    {
        *price = 0.0;
        status = SUCCESS;
        goto RETURN;
    }    
    
    /* create date mesh */
    for(j = 0;j < numName;j++)
    {
        pPrice[j] = 0.0;
    }

    for(i = pPayoff->numDate -1;i >= 0;i--)
    {
        /* init */
        endDate = pPayoff->payDate[i];
        
        if(i == 0) startDate = today;
        if(i > 0) startDate = (pPayoff->payDate[i-1] >= today)? pPayoff->payDate[i-1]:today;

        if(startDate > endDate)
        {
            CRXError("%s: Payoff date list not in ascending order!\n", routine);
            goto RETURN;
        }
        
        /* payoff */
        for(j = 0;j < numName;j++)
        {
            tmpPayoff[j] = pPayoff->payoff[i][j] + pPrice[j];
        }

        if(pPayoff->cutoffDate < endDate)
        {
            if(TransMat(&mat,
                        startDate,
                        pPayoff->cutoffDate,
                        effCurves) != SUCCESS)
            {
                goto RETURN;
            }
            
            if(MatMultVecUT(pPrice,&mat,tmpPayoff) != SUCCESS)
            {
                goto RETURN;
            }
        } else {            
            if(TransMat(&mat,
                        startDate,
                        endDate,
                        effCurves) != SUCCESS)
            {
                goto RETURN;
            }
            
            if(MatMultVecUT(pPrice,&mat,tmpPayoff) != SUCCESS)
            {
                goto RETURN;
            }
        }
        
        if(GtoDiscountDateForward(startDate,
                                  endDate,
                                  irCurve,
                                  GTO_FLAT_FORWARDS,
                                  &discFac) != SUCCESS)
        {
            goto RETURN;
        }

        for(j = 0;j < numName;j++)
        {
            pPrice[j] *= discFac;
        }

        if(startDate == today)break;
    }

    if(GtoDiscountDateForward(today,
                              valueDate,
                              irCurve,
                              GTO_FLAT_FORWARDS,
                              &discFac) != SUCCESS)
    {
        goto RETURN;
    }
        
    for(i = 0;i < numName;i++)
    {
        price[i] = pPrice[i]/discFac;
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
 * extend payoff to contain exercise dates
 *
 ********************************************************************************/
int ExtendPayoff(
    PayOff     *pPayoffNew,               /* (O) extended payoff                */
    int        *exerFlag,                 /* (O) exer flag, 0: not exer, 1: yes */
    PayOff     *pPayoffOld,               /* (I) original payoff                */
    TDate      *pDate,                    /* (I) exercise dates                 */
    int         numDate,                  /* (I) num. of exer. dates            */
    int         numName)                  /* (I) num of names                   */
{
    static char routine[] = "ExtendPayoff";
    
    int       status  = FAILURE;
    int       i, j, idx;

    if(pPayoffOld->numDate < 0)
    {
        CRXError("%s: Payoff num dates < 0.\n", routine);
        goto RETURN;
    }

    pPayoffNew->cutoffDate = pPayoffOld->cutoffDate;
    
    idx = 0;
    j = 0;
    
    if(pPayoffOld->numDate == 0)
    {
        while(idx < numDate)
        {
            pPayoffNew->payDate[j] = pDate[idx];
            PtrSet(&pPayoffNew->payoff[j][0],numName,0.0);
            exerFlag[j] = 1;
            j++;
            idx++;
        }

    } else {
        for(i = 0;i < pPayoffOld->numDate;i++)
        {
            while(pDate[idx] < pPayoffOld->payDate[i] && idx < numDate)
            {
                pPayoffNew->payDate[j] = pDate[idx];
                PtrSet(&pPayoffNew->payoff[j][0],numName,0.0);
                exerFlag[j] = 1;
                j++;
                idx++;
            }
            
            if(pDate[idx] == pPayoffOld->payDate[i] && idx < numDate)
            {
                pPayoffNew->payDate[j] = pDate[idx];
                PtrCopy(&pPayoffNew->payoff[j][0],&pPayoffOld->payoff[i][0],numName);
                exerFlag[j] = 1;
                j++;
                idx++;
            } else {
                pPayoffNew->payDate[j] = pPayoffOld->payDate[i];
                PtrCopy(&pPayoffNew->payoff[j][0],&pPayoffOld->payoff[i][0],numName);
                exerFlag[j] = 0;
                j++;
            }
        }
    }
    
    pPayoffNew->numDate = j;

    status = SUCCESS;
RETURN:
    if(status != SUCCESS)
    {
        CRXError("%s: Failed.\n", routine);
    }
    return status;
}

/*********************************************************************************
 * Insert payoff at a certain pos
 *
 ********************************************************************************/
int PayoffInsert(
    PayOff    *pPayoff,                   /* (I/O) payoff                       */
    int       pos,                        /* (I)   pos to insert into           */
    TDate     date,                       /* (I)   payoff date                  */
    int       numName)                    /* (I)   num of tranchelets           */
{
    int i,j;
    
    for(i = pPayoff->numDate;i >= pos+1;i--)
    {
        pPayoff->payDate[i] = pPayoff->payDate[i-1];
        for(j = 0;j < numName;j++)
        {
            pPayoff->payoff[i][j] = pPayoff->payoff[i-1][j];
        }
    }
    pPayoff->payDate[pos]  = date;
    
    for(j = 0;j < numName;j++)
    {
        pPayoff->payoff[pos][j] = 0.0;
    }

    pPayoff->numDate++;
    
    return SUCCESS;
}

