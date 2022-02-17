/***************************************************************
 * Module:      cranalyst
 * File:        indexcmcds.c
 * Description: composite fixing adjustment for index CMCDS
 ***************************************************************/

#include <common/include/drutils.h>
#include <ctype.h>
#include <indexcmcds.h>
#include "imsl.h"
#include "alimsl.h"
#include "interp.h"
#include "crxutil.h"
#include "crcrv.h"
#include "alib/dtlist.h"
#include "alib/mgpimsl.h"


/** CDS parameter **
 */
typedef struct
{
    TDate         today;            /**(I) Trade date for new curve         */
    TDate         valDate;          /**(I) Spot-date for new curve          */
    TDateInterval cdsFreq;          /**(I) CDS fee payment frequency        */
    TDayCountConv cdsDCC;           /**(I) CDS fee daycount                 */
    KStubLocation dateStub;         /**(I) CDS fee stub type                */
    KAccrualConv  payAccFlag;       /**(I) Pay accrued on default or not    */
    KProtPayConv  protPayType;      /**(I) Pay on default, or at maturity   */
    TDateInterval delay;            /**(I) Payment delay                    */
    double        recovery;         /**(I) Recovery rate                    */
    TDateInterval frequency;        /**(I) Protection leg integration freq. */
    TCurve        *discZC;          /**(I) Interest rate zero curve         */
    TCurve        *sprdZC;          /**(I) Credit zero curve                */
    int           status;           /**(O) status                           */
} TCDSParams;



/** Par spread calculation parameters. Used by the solver
 */
typedef struct
{
    int           nbFwds;
    int           nPoints;
    TDate         today;            /**(I) Trade date for new curve         */
    TDate         valDate;          /**(I) Spot-date for new curve          */
    TDate         *cdsStDate;       /**(I) Start date for CDS               */
    TDate         *cdsMatDate;      /**(I) CDS maturity                     */
    TDateInterval cdsFreq;          /**(I) CDS fee payment frequency        */
    TDayCountConv cdsDCC;           /**(I) CDS fee daycount                 */
    KStubLocation dateStub;         /**(I) CDS fee stub type                */
    KAccrualConv  payAccFlag;       /**(I) Pay accrued on default or not    */
    KProtPayConv  protPayType;      /**(I) Pay on default, or at maturity   */
    TDateInterval delay;            /**(I) Payment delay                    */
    double        recovery;         /**(I) Recovery rate                    */
    TDateInterval frequency;        /**(I) Protection leg integration freq. */
    TCurve        *discZC;          /**(I) Interest rate zero curve         */
    TCurve        *sprdZC;          /**(I) Credit zero curve                */
    TDate         *spotDates;       /**(I) Target dates                     */
    double        *targetFwdRates;  /**(I) Target solution                  */
    double        *actFwdRates;     /**(O) Solver's solution                */
    TCurve        *tmpZC;           /**(O) temp credit ZC                   */
    int           status;           /**(O) status                           */
} TCDSZCParams;

/**  global variable for solver **/
TCDSZCParams p;

/** imsl passin functions  **/
double dummy(int n, double *x)
{
    return 0;
}

void grad(int n, 
          double *spotRates, 
          double *g);

double matchFwdRates(int n, 
                     double *spotRates);


/** Helper function to get fwd and annuity **/
int getFwdAndAnnuity(TCDSParams   *q,             /** (I) parameters */
                     TDate         start,         /** (I) start date */
                     TDate         end,           /** (I) end date   */
                     double       *fwd,           /** (O) forward    */
                     double       *annuity);      /** (O) annuity    */

/** Helper function to generate zero dates **/
int generateSpotDates(int nbFwds,
                      int *nPoints,
                      TDate* fwdSttDates,
                      TDate* fwdEndDates,
                      TDate** spotDates);

/**********************************************************************************
 ** Composite fixing adjustment for cmcds index.
 ** Return SUCCESS/FAILURE.
 */
int CrxFwdsToCrv(TDate         today,            /**<(I) Trade date                           */
                 TDate         valDate,          /**<(I) Spot-date                            */
                 const TCurve  *discZC,          /**<(I) Interest rate zero curve             */
                 int           nbFwds,           /**<(I) Number of fwds                       */
                 TDate         *fwdSttDates,     /**<(I) Start dates of fwds to be adjusted   */
                 TDate         *fwdEndDates,     /**<(I) Fwd end dates                        */
                 double        *targetFwdRates,  /**<(I) Target fwd rates                     */
                 double        recovery,         /**<(I) Recovery rate                        */
                 TDateInterval cdsFreq,          /**<(I) CDS fee payment frequency            */
                 TDayCountConv cdsDCC,           /**<(I) CDS fee daycount                     */
                 KStubLocation dateStub,         /**<(I) CDS fee stub type                    */
                 KAccrualConv  payAccFlag,       /**<(I) Pay accrued on default or not        */
                 KProtPayConv  protPayType,      /**<(I) Pay on default, or at maturity       */
                 TDateInterval frequency,        /**<(I) Protection leg integration freq.     */
                 TDateInterval delay,            /**<(I) Payment delay                        */
                 int           nPoints,          /**<(I) number of spot dates                 */
                 TDate         *spotDates,       /**<(I) spot dates                           */
                 double        *inGuessRates,    /**<(I) guesss rates, size of nPoints        */
                 double        *actFwdRates,     /**<(O) actual fwd rates                     */
                 TCurve        **newSprdZC)      /**<(O) New Clean spread zero curve          */
{
    static char routine[] = "CrxFwdsToCrv";
    int         status    = FAILURE;

    double *solveRates = NULL;
    double fValue;
    double *guessRates = NULL;
    int i;


    /** Do some basic consistency check */
    if(nbFwds < 1) {
        DR_Error("%s: number of fwds must be greater than zero.", routine);
        goto done;
    }

    if(nPoints < 1) 
    {
        DR_Error("%s: number of spot dates must be greater than zero.", routine);
        goto done;
    }

    if(!discZC) 
    {
        DR_Error("%s: interest rate curve is NULL.", routine);
        goto done;
    }

    if(!inGuessRates)
    {
        guessRates =  (double *) MALLOC (nPoints * sizeof(double));
        for(i = 0; i < nPoints; i++)
        {
            guessRates[i] = 0.01;
        }
    } else {
        guessRates = inGuessRates;
    }

    if(fwdSttDates[0] < today)
    {
        DR_Error("%s: Input fwd start date is before today", routine);
        goto done;
    }

    for(i = 0; i < nbFwds; i++) 
    {
        if(fwdSttDates[i] >= fwdEndDates[i]) 
        {
            DR_Error("%s: fwd start date must be before end date.", routine);
            goto done;
        }

        if(targetFwdRates[i] <= 0)
        {
            DR_Error("%s: target fwd rate [%d] less than 0: %f ", routine, i, targetFwdRates[i]);
            goto done;
        }

        if(i)
        {
            if(fwdSttDates[i] < fwdSttDates[i-1])
            {
                DR_Error("%s: fwd start date must be increasing.", routine);
                goto done;
            }

            if(fwdEndDates[i] < fwdEndDates[i-1])
            {
                DR_Error("%s: fwd end date must not be decreasing.", routine);
                goto done;
            }
        }
    }

    /** fill in members in p **/
    p.nbFwds = nbFwds;
    p.nPoints = nPoints;
    p.today = today;
    p.valDate = valDate;
    p.cdsStDate = fwdSttDates;
    p.cdsMatDate = fwdEndDates;
    p.cdsFreq = cdsFreq;
    p.cdsDCC = cdsDCC;
    p.dateStub = dateStub;
    p.payAccFlag = payAccFlag;
    p.protPayType = protPayType;   
    p.delay = delay;         
    p.recovery = recovery;      
    p.frequency = frequency;     
    p.spotDates = spotDates;
    p.discZC = (TCurve*)discZC;
    p.sprdZC = NULL;
    p.targetFwdRates = targetFwdRates;
    p.actFwdRates = actFwdRates;

    /** new curve **/
    p.tmpZC = GtoNewTCurve(today, nPoints, 1, GTO_ACT_365F);


    /** Solve for new spot rates **/
    GtoImslCmathInit();
    solveRates = imsl_d_min_uncon_multivar(matchFwdRates,
                                           nPoints,
                                           IMSL_XGUESS, guessRates,
                                           IMSL_FVALUE, &fValue,
                                           IMSL_GRAD, grad,
                                           0);

    if(solveRates) {
        DR_Error("%s: TotalError=%13.7e %", routine, fValue);
        DR_Error(" i, startDate,   endDate, targetFwdRates,  actFwdRates[i],   diff[i]");

        for(i=0; i < nbFwds; i++) {
            DR_Error("%2d %9d %10d %16.8e %16.8e %16.8e", i, fwdSttDates[i], fwdEndDates[i], 
                     p.targetFwdRates[i], actFwdRates[i], p.actFwdRates[i]-p.targetFwdRates[i]);
        }

        DR_Error(" i, date,     spotRates");
        for(i=0; i < nPoints; i++) {
            DR_Error("%2d %10d %16.8e", i, p.spotDates[i], solveRates[i]);
        }
    }

    /** check results */
    for(i = 0; i < nPoints; i++)
    {
        if(solveRates != NULL && solveRates[i] <= 0 )
        {
            DR_Error("%s: Negagive spot rate[%d] %f%%", routine, i, solveRates[i]*100);
            p.status = FAILURE;
        }
    }

    /*    if(solveRates == NULL || p.status == FAILURE) */
    if(!solveRates)
    {
        DR_Error("%s: Can not solve for spot rates", routine);
        goto done;
    }

    (*newSprdZC) = GtoCopyCurve(p.tmpZC);
    if(!(*newSprdZC))
    {
        DR_Error("%s: New curve creation failed", routine);
        goto done;
    }

    for(i=0; i<nPoints; i++)
    {
        (*newSprdZC)->fArray[i].fRate = solveRates[i];
    }
    

    status = SUCCESS;
 done:
    GtoFreeTCurve(p.tmpZC);
    FREE_ARRAY(solveRates);
    if(!inGuessRates) 
    {
        FREE_ARRAY(guessRates);
    }

    return status;
}


/**********************************************************************************
 ** Composite fixing adjustment for cmcds index.
 ** Return SUCCESS/FAILURE.
 */
int CrxCompositeCMCDSCrv(TDate         today,            /**<(I) Trade date                           */
                         TDate         valDate,          /**<(I) Spot-date                            */
                         TDateInterval cmcdsTenor,       /**<(I) CMCDS tenor, e.g. 3Y                 */
                         const TCurve  *discZC,          /**<(I) Interest rate zero curve             */
                         const TCurve  *sprdZC,          /**<(I) Index clean curve                    */
                         int           fwdAdj,           /**<(I) Y/N flag for fwd adjustment          */
                         int           matAdj,           /**<(I) Y/N flag for maturity adjustment     */
                         int           nbFwds,           /**<(I) Number of fwds                       */
                         TDate         *fwdSttDates,     /**<(I) Start dates of fwds to be adjusted   */
                         TDate         *fwdEndDates,     /**<(I) Fwd end dates                        */
                         double        recovery,         /**<(I) Recovery rate                        */
                         TDateInterval cdsFreq,          /**<(I) CDS fee payment frequency            */
                         TDayCountConv cdsDCC,           /**<(I) CDS fee daycount                     */
                         KStubLocation dateStub,         /**<(I) CDS fee stub type                    */
                         KAccrualConv  payAccFlag,       /**<(I) Pay accrued on default or not        */
                         KProtPayConv  protPayType,      /**<(I) Pay on default, or at maturity       */
                         TDateInterval frequency,        /**<(I) Protection leg integration freq.     */
                         TDateInterval delay,            /**<(I) Payment delay                        */
                         TCurve        **newSprdZC)      /**<(O) New Clean spread zero curve          */
{
    static char routine[] = "CrxCompositeCMCDSCrv";
    int         status    = FAILURE;
    
    int         nPoints;
    double      *actFwdRates = NULL;
    TDate       *fwdEndNoRollDates = NULL;
    TDate       *newStartDates = NULL;
    TDate       *newEndDates = NULL;
    TDate       *spotDates = NULL;
    double      *guessRates = NULL;
    TDateInterval adjInterval;
    double        fwd1, fwd2, ann1, ann2;
    double       *adjusts = NULL;
    double       *adjusted = NULL;
    double       mats;
    int          i;
    double       xtemp;
    double       tenor;
    TDate        dateTemp;
    KFeeLeg_D    *fl = NULL;
    TCurve       *cc = NULL;
    TCurve       *ccStub = NULL;
    TCDSParams   q;

    /** Do some basic consistency check */
    if(nbFwds < 1) {
        DR_Error("%s: number of fwds must be greater than zero.", routine);
        goto done;
    }

    if(!discZC) {
        DR_Error("%s: interest rate curve is NULL.", routine);
        goto done;
    }
    
    if(!sprdZC) {
        DR_Error("%s: credit rate curve is NULL.", routine);
        goto done;
    }
    
    if(fwdSttDates[0] < today){
        DR_Error("%s: Input fwd start date is before today", routine);
        goto done;
    }

    for(i = 0; i < nbFwds; i++) 
    {
        if(fwdSttDates[i] >= fwdEndDates[i]) 
        {
            DR_Error("%s: fwd start date must be before end date.", routine);
            goto done;
        }

        if(i)
        {
            if(fwdSttDates[i] < fwdSttDates[i-1])
            {
                DR_Error("%s: fwd start date must be increasing.", routine);
                goto done;
            }

            if(fwdEndDates[i] < fwdEndDates[i-1])
            {
                DR_Error("%s: fwd end date must not be decreasing.", routine);
                goto done;
            }
        }
    }

    /** get adjustment tenor **/
    GtoDateIntervalToYears(&cmcdsTenor, &tenor);
    GtoDateIntervalToYears(&cdsFreq, &xtemp);
    GtoYearsToDateInterval(tenor-xtemp, &adjInterval);

    /** generate spot dates **/
    if(generateSpotDates(nbFwds,
                         &nPoints,
                         fwdSttDates,
                         fwdEndDates,
                         &spotDates) == FAILURE)
    {
        DR_Error("%s: failed to generate spot dates.", routine);
        goto done;

    }

    /** fwd dates **/
    fwdEndNoRollDates = (TDate *) MALLOC (nbFwds * sizeof(TDate));
    newStartDates = (TDate *) MALLOC (nbFwds * sizeof(TDate));
    newEndDates = (TDate *) MALLOC (nbFwds * sizeof(TDate));
    guessRates =  (double *) MALLOC (nPoints * sizeof(double));

    fwdEndNoRollDates[0] = fwdEndDates[0];
    
    for(i = 0; i < nbFwds; i++) 
    {
        newStartDates[i] = fwdSttDates[i];
        newEndDates[i] = fwdEndDates[i];
        if(i) {
            fwdEndNoRollDates[i] = newEndDates[i-1];
        }
    }

    for(i = 0; i < nPoints; i++)
    {
        GtoInterpRate(spotDates[i], (TCurve*)sprdZC, INTERP_METHOD, &(guessRates[i]));
    }


    /** fill in members in q */
    q.today = today;
    q.valDate = valDate;
    q.cdsFreq = cdsFreq;
    q.cdsDCC = cdsDCC;
    q.dateStub = dateStub;
    q.payAccFlag = payAccFlag;
    q.protPayType = protPayType;   
    q.delay = delay;         
    q.recovery = recovery;      
    q.frequency = frequency;     
    q.discZC = (TCurve*)discZC;
    q.sprdZC = (TCurve*)sprdZC;


    /** maturity and forward adjustments
     */
    adjusts = (double*) MALLOC(nbFwds * sizeof(double));
    adjusted = (double*) MALLOC(nbFwds * sizeof(double));

	for(i = 0; i < nbFwds; i++)
	{
        getFwdAndAnnuity(&q, newStartDates[i], newEndDates[i], &(adjusted[i]), &ann2);
	}

    getFwdAndAnnuity(&q, today, newStartDates[0], &fwd1, &ann1);
    GtoDateFromDateAndOffset(newStartDates[0], &adjInterval, 1, &dateTemp);
    getFwdAndAnnuity(&q, today, dateTemp, &fwd2, &ann2);

    if(fwdAdj) {
        adjusts[0] = fwd1 * ann1 / (ann2 - ann1);
    } else {
        adjusts[0] = 0;
    }
    adjusted[0] += adjusts[0];
    mats = 0;

    for(i = 1; i < nbFwds; ++i)
    {
        adjusts[i] = 0;

        /** t -> T */
		if(matAdj)
		{
            getFwdAndAnnuity(&q, newStartDates[i], fwdEndNoRollDates[i], &fwd1, &ann1);
		    getFwdAndAnnuity(&q, newStartDates[i], newEndDates[i], &fwd2, &ann2);
			mats += fwd1 - fwd2;  /** maturity adj */
            adjusts[i] += mats;
		}

        /** 0 -> t */
        if(fwdAdj)
        {
            getFwdAndAnnuity(&q, today, newStartDates[i], &fwd1, &ann1);
            GtoDateFromDateAndOffset(newStartDates[i], &adjInterval, 1, &dateTemp);
            getFwdAndAnnuity(&q, today, dateTemp, &fwd2, &ann2);
            adjusts[i] += (fwd1*ann1)/(ann2-ann1);
        }

        /** adjusted fwds */
        adjusted[i] += adjusts[i];  
    }
    

    actFwdRates = (double *) MALLOC (nbFwds * sizeof(double));
    if(!actFwdRates) {
        DR_Error("%s: unable to allocate memory for actFwdRates", routine);
        goto done;
    }


    if(CrxFwdsToCrv(q.today,
                    q.valDate,
                    q.discZC,
                    nbFwds,
                    newStartDates,
                    newEndDates,
                    adjusted,
                    q.recovery,
                    q.cdsFreq,
                    q.cdsDCC,
                    q.dateStub,
                    q.payAccFlag,
                    q.protPayType,
                    q.frequency,
                    q.delay,
                    nPoints,
                    spotDates,
                    guessRates,
                    actFwdRates,
                    newSprdZC) == FAILURE)
    {
        DR_Error("%s: unable to generate credit clean curve to match fwds", routine);
    }


    status = SUCCESS;

 done:
    FREE_ARRAY(spotDates);
    FREE_ARRAY(fwdEndNoRollDates);
    FREE_ARRAY(guessRates);
    FREE_ARRAY(adjusts);
    FREE_ARRAY(adjusted);
    FREE_ARRAY(actFwdRates);
    FREE_ARRAY(newStartDates);
    FREE_ARRAY(newEndDates);

    return status;
}


/** Helper function, return fwd spread and annuity
 */
int getFwdAndAnnuity(TCDSParams   *q,             /** (I) parameters */
                     TDate         start,         /** (I) start date */
                     TDate         end,           /** (I) end date   */
                     double       *fwd,           /** (O) forward    */
                     double       *annuity)       /** (O) annuity    */
{
    return CrxFwdParCDSSpread(fwd,
                              annuity,
                              q->today,
                              q->valDate,
                              start,
                              end,
                              q->cdsFreq,
                              q->cdsDCC,
                              q->dateStub,
                              q->payAccFlag,
                              q->protPayType,
                              q->delay,
                              q->recovery,
                              q->frequency,
                              q->discZC,
                              q->sprdZC);
}
                     
/** imsl required functions  **/
double matchFwdRates(int n, 
                     double *spotRates)
{
    int i;
    double errSum = 0;
    double err;
    double annuity;
    double *parSprd;

    for(i = 0; i < n; ++i)
    { 
        if(spotRates[i] <= 0 )
        {
            spotRates[i] = 0.00001;
        }
        p.tmpZC->fArray[i].fDate = p.spotDates[i];
        p.tmpZC->fArray[i].fRate = spotRates[i];
    }

    for(i = 0, parSprd = p.actFwdRates; i < p.nbFwds; ++i, ++parSprd)
    {
        CrxFwdParCDSSpread(parSprd,
                           &annuity,
                           p.today,
                           p.valDate,
                           p.cdsStDate[i],
                           p.cdsMatDate[i],
                           p.cdsFreq,
                           p.cdsDCC,
                           p.dateStub,
                           p.payAccFlag,
                           p.protPayType,
                           p.delay,
                           p.recovery,
                           p.frequency,
                           p.discZC,
                           p.tmpZC);
        err = (*parSprd) - p.targetFwdRates[i];
        errSum += err * err;
    }

    //errSum = sqrt(errSum);
    return errSum;
}

void grad(int n, 
          double *spotRates, 
          double *g)
{
    int i;
    int j;
    double parSprd, annuity;
    double err, sumErr;
    static const double SHIFT = 0.0001;
    TCurve *tmpZC;

    tmpZC = GtoNewTCurve(p.today, n, 1, GTO_ACT_365F);

    for(i = 0; i < n; ++i)
    { 
        tmpZC->fArray[i].fDate = p.spotDates[i];
        tmpZC->fArray[i].fRate = spotRates[i];
    }

    for(j = 0; j < n; j++)
    {
        g[j] = 0;
        sumErr = 0;

        tmpZC->fArray[j].fRate += SHIFT;
        
        for(i = 0; i < p.nbFwds; ++i)
        {
            err = 2*(p.actFwdRates[i] - p.targetFwdRates[i]);
            sumErr += err * err;

            CrxFwdParCDSSpread(&parSprd,
                               &annuity,
                               p.today,
                               p.valDate,
                               p.cdsStDate[i],
                               p.cdsMatDate[i],
                               p.cdsFreq,
                               p.cdsDCC,
                               p.dateStub,
                               p.payAccFlag,
                               p.protPayType,
                               p.delay,
                               p.recovery,
                               p.frequency,
                               p.discZC,
                               tmpZC);
            g[j] += err * (parSprd - p.actFwdRates[i]) / SHIFT;
        }

        //g[j] /= 0.5 * sqrt(sumErr + 1.0e-60);

        tmpZC->fArray[j].fRate -= SHIFT;
    }

    GtoFreeTCurve(tmpZC);
}


/** Helper function: generate spot dates **/
int generateSpotDates(int nbFwds,
                      int *nPoints,
                      TDate* fwdSttDates,
                      TDate* fwdEndDates,
                      TDate** spotDates)
{
    static char routine[] = "generatSpotDates";
    TDateList *tmp1 = NULL;
    TDateList *tmp2 = NULL;
    TDateList *tmp = NULL;
    int i;
    int status = FAILURE;
    
    FREE_ARRAY(*spotDates);
    *spotDates = NULL;


    //tmp2 has all the end dates
    tmp2 = GtoNewDateListFromDates(fwdEndDates, nbFwds);
    if(!tmp2) {
        DR_Error("%s: Cannot create datelist form fwdEndDates", routine);
        goto done0;
    }

    //tmp1 has all the start dates
    tmp1 = GtoNewDateListFromDates(fwdSttDates, nbFwds);
    if(!tmp1) {
        DR_Error("%s: Cannot create datelist form fwdSttDates", routine);
        goto done0;
    }
    
    //Merge
    tmp = GtoMergeDateListDuplicates(tmp1, tmp2);
    if(!tmp) {
        DR_Error("%s: Cannot create merge datelist", routine);
        goto done0;
    }
    GtoFreeDateList(tmp1);
    tmp1 = NULL;
    GtoFreeDateList(tmp2);
    tmp2 = NULL;


    //tmp1 has distince dates
    tmp1 = GtoNewDateListDistinct(tmp);
    if(!tmp1) {
        DR_Error("%s: Cannot create distinct datelist", routine);
        goto done0;
    }
    GtoFreeDateList(tmp);
    tmp = NULL;
    

    //fill in spotDates with tmp
    *nPoints = tmp1->fNumItems;
    *spotDates = (TDate *) MALLOC (tmp1->fNumItems * sizeof(TDate));
    for(i = 0; i < tmp1->fNumItems; i++) {
        (*spotDates)[i] = tmp1->fArray[i];
    }

    GtoFreeDateList(tmp1);
    tmp1 = NULL;

    status = SUCCESS;
 done0:
    GtoFreeDateList(tmp);
    GtoFreeDateList(tmp1);
    GtoFreeDateList(tmp2);

    return status;
}
