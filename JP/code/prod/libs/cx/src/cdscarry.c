 /*
*************************************************************************
**FILE NAME: cdscarry.c
**
** CDS Carry functions
**
**
*************************************************************************
*/

#include "cdscarry.h"

#include <math.h>

#include "cds.h"
#include "cdsbootstrap.h"

#include <alib/rtbrent.h>

#include <cxutils/include/cxmacros.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/dateutils.h>
#include <cxutils/include/zerocurve.h>

typedef struct _BREAKEVEN_PARAMS
{
  TDate horizonDate;
  TCurve* discCurve;
  TDate startDate;
  TDate date1;
  TDate date2;
  double notional1;
  double notional2;
  double spread1;
  double fee1;
  double fee2;
  CxTRecoveryCurve* recoveryCurve;
  CxTCdsConventions*   cdsConventions;
  CxTCalendar* calendar;
  double target;
}BREAKEVEN_PARAMS;


static double FwZeroPrice
(TCurve        *discCurve,
 TDate         startDate,
 TDate         endDate);

static double FwSurvPorb
(CxTCreditCurve        *spreadCurve,
 TDate         startDate,
 TDate         endDate);

static int BreakevenSolver(
    double spread,
    BREAKEVEN_PARAMS *params,
    double *priceDiff);

CxTCdsCarryOutputs* CxCdsCarry
(TDate           today,
 TDate           valueDate,
 TDate           startDate,
 TDate           endDate,
 TDate            horizonDate,
 double          couponRate,
 TDateInterval  *dateInterval,
 CxTStubType      stubType,
 CxTDayCountConv  paymentDcc,
 CxTBadDayConv    badDayConv,
 CxTCalendar     *calendar,
 TCurve         *discCurve,
 TBoolean        protectStart)
{
  static char routine[] = "CxCdsCarry";
  int status = FAILURE;
  CxTFeeLeg *fl = NULL;


  /*We discount the flows back to value date*/
  double valueDatePV;
  int i;

  /*Outputs*/
  double couponsPV = 0.0;
  double horizonAccrual = 0.0;
  double carry;
  CxTCdsCarryOutputs *carryOutputs = NULL;
  
  valueDatePV =  FwZeroPrice(discCurve,today,valueDate);

  REQUIRE (horizonDate > startDate);
  /*Other requirements will be made by the called functions*/

    fl = CxCdsFeeLegMake (startDate, endDate, CX_ACCRUAL_PAY_NONE,
                          dateInterval, stubType, 1.0,
                          1.0, paymentDcc, badDayConv,
                          calendar, protectStart);
    if (fl == NULL)
        goto done; /* failure */

    for (i=0 ; i < fl->nbDates ; ++i)
      {
	double thisCouponPV;
	double accTime;
	double discount;


	if (fl->accEndDates[i] <= horizonDate)
	  {
	accTime = CxDayCountFraction (fl->accStartDates[i],
				      fl->accEndDates[i],
				      paymentDcc);
	discount = FwZeroPrice(discCurve,
			       today,
			       fl->accEndDates[i]);
	thisCouponPV = accTime * couponRate * discount / valueDatePV;
	couponsPV += thisCouponPV;
	  }


	if (fl->accEndDates[i] > horizonDate)
	  {
	accTime = CxDayCountFraction (fl->accStartDates[i],
				      horizonDate,
				      paymentDcc);
	discount = FwZeroPrice(discCurve,
			       today,
			       fl->accEndDates[i]);

	horizonAccrual = accTime * couponRate * discount / valueDatePV;
	i=fl->nbDates;
	  }
      }
    carry = couponsPV + horizonAccrual;
    carryOutputs = CxCdsCarryOutputsMake(couponsPV,
					 horizonAccrual,
					 carry);
    
    status = SUCCESS;

 done:

  if (status != SUCCESS)
    {
       GtoErrMsgFailure (routine);
       CxCdsCarryOutputsFree(carryOutputs);
       carryOutputs = NULL;
    }

    CxFeeLegFree (fl);
    return carryOutputs;

}
/* */
int CxCdsSlideBreakeven(
			 TDate horizonDate,
			 TCurve *discCurve, /* as of Horizon Date */
			 TDate startDate,
			 TDate maturityDate1,
			 TDate maturityDate2,
			 double signedNotional1,
			 double signedNotional2,
			 double spread1, /* the spread you suppose maturityDate1 spread is at */
			 double fee1, /* The fee you are paying */
			 double fee2,
			 CxTRecoveryCurve *recoveryCurve,
			 CxTCdsConventions  *cdsConventions, /* conventions for the CDS pricing and bootstrapping */
			 CxTCalendar *calendar,
			 double target,
			 double *spread2) /* We will solve for the slide to be equal to target */
{
  static char routine[] = "CdsSlideBreakeven";
  
  int status = FAILURE;
  
  BREAKEVEN_PARAMS params;

    /* Parameters for root finder */
    double boundLo       = 0.0;
    double boundHi       = 1.0;
    int    numIterations = 100;
    double guess         = spread1;
    double initialXStep  = 0.0001; /*1 bp*/
    double initialFDeriv = 0.0;
    double xacc          = 2; /* will not put restrictions on x accuracy as spread between 0 and 1 */
    double facc          = MIN (ABS(signedNotional1),ABS(signedNotional2)) * 1e-8;

    REQUIRE(maturityDate1 > horizonDate);
    REQUIRE(maturityDate2 > horizonDate);
    
  params.horizonDate = horizonDate;
  params.discCurve = discCurve;
  params.startDate =  startDate;
  params.date1 = maturityDate1;
  params.date2 = maturityDate2;
  params.notional1 = signedNotional1;
  params.notional2 = signedNotional2;
  params.spread1 = spread1;
  params.fee1 = fee1;
  params.fee2 = fee2;
  params.recoveryCurve = recoveryCurve;
  params.cdsConventions = cdsConventions;
  params.calendar =  calendar;
  params.target = target;

  status = GtoRootFindBrent ((TObjectFunc)BreakevenSolver,
			     &params,
			     boundLo,
			     boundHi,
			     numIterations,
			     guess,
			     initialXStep,
			     initialFDeriv,
			     xacc,
			     facc,
			     spread2);

 done:
    if (status != SUCCESS)
      GtoErrMsgFailure (routine);
  
  return status;
    
}



			 
static double FwZeroPrice
(TCurve        *discCurve,
 TDate         startDate,
 TDate         endDate)
{
  double df = 1.0 ;

  if (discCurve != NULL)
    
     df = CxForwardZeroPrice(discCurve,
			       startDate,
			       endDate);

  return df;
}

static double FwSurvPorb
(CxTCreditCurve        *spreadCurve,
 TDate         startDate,
 TDate         endDate)
{
  double df = 1.0 ;

  if (spreadCurve != NULL)
    
    df = CxForwardZeroPrice(spreadCurve->tc,
			    startDate,
			    endDate);
  
  return df;
}

static int BreakevenSolver(
    double spread,
    BREAKEVEN_PARAMS *params,
    double *priceDiff)
{
  int status = FAILURE;
  TDate *endDates = NULL;
  double *spreads = NULL;
  CxTCreditCurve *spreadCurve = NULL;
  double price1 = 0.0;
  double price2 = 0.0;
  
  endDates = NEW_ARRAY(TDate,2);
  endDates[0]=params->date1;
  endDates[1]=params->date2;
  
  spreads = NEW_ARRAY(double,2);
  spreads[0]=params->spread1;
  spreads[1]=spread;
  
  spreadCurve =  CxCdsBootstrap (
				 params->horizonDate,
				 params->discCurve,
				 params->startDate,
				 params->horizonDate,
				 2,
				 endDates,
				 spreads,
				 NULL,
				 NULL,
				 params->recoveryCurve,
				 params->cdsConventions->payAccOnDefault,
				 params->cdsConventions->couponInterval,
				 params->cdsConventions->paymentDCC,
				 params->cdsConventions->stubType,
				 params->cdsConventions->curveType,
				 params->cdsConventions->timestep,
				 params->cdsConventions->smoothInterval,
				 params->cdsConventions->protectStart,
				 params->cdsConventions->isPriceClean,
				 params->cdsConventions->delay,
				 params->cdsConventions->badDayConv,
				 params->calendar);
  
  if (spreadCurve == NULL)
    goto done; /* Failure */
  
  CxCdsPrice (
	      params->horizonDate, 
	      params->horizonDate,
	      params->startDate,
	      params->date1,
	      params->cdsConventions->delay,
	      params->fee1,
	      params->cdsConventions->payAccOnDefault,
	      params->cdsConventions->couponInterval,
	      params->cdsConventions->stubType,
	      params->cdsConventions->paymentDCC,
	      params->cdsConventions->badDayConv,
	      params->calendar,
	      params->discCurve,
	      spreadCurve,
	      params->recoveryCurve,
	      params->cdsConventions->protectStart,
	      params->cdsConventions->isPriceClean,
	      &price1);

  CxCdsPrice (
	      params->horizonDate, 
	      params->horizonDate,
	      params->startDate,
	      params->date2,
	      params->cdsConventions->delay,
	      params->fee2,
	      params->cdsConventions->payAccOnDefault,
	      params->cdsConventions->couponInterval,
	      params->cdsConventions->stubType,
	      params->cdsConventions->paymentDCC,
	      params->cdsConventions->badDayConv,
	      params->calendar,
	      params->discCurve,
	      spreadCurve,
	      params->recoveryCurve,
	      params->cdsConventions->protectStart,
	      params->cdsConventions->isPriceClean,
	      &price2);


   *priceDiff = params->notional1 * price1 + params->notional2 * price2 - params->target;

   status = SUCCESS;
 
 done:

   FREE (endDates);
   FREE(spreads);
   CxCreditCurveFree(spreadCurve);
   

   return status;
 
}
