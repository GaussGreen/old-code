/*
*************************************************************************
**FILE NAME: bond.c
**
** Bond testing functions
** Quick and not too clean implementation
** This is still under test
*************************************************************************
*/

#include "bond.h"

#include <math.h>

#include "cds.h"
#include "cdsbootstrap.h"
#include "cdscarry.h" /* only for discount functions*/
#include "contingentleg.h"
#include "recovery.h"
#include "feeleg.h"

#include <alib/bondflow.h>
#include <alib/bondcnst.h>


#include <cxutils/include/cxmacros.h>
#include <cxutils/include/datelist.h>
#include <cxutils/include/dateutils.h>
#include <cxutils/include/zerocurve.h>

int  CxBondTest(TBond* bond,
		double* output)
{
  *output =2;
  return SUCCESS;
}


static double FwZeroPrice
(TCurve        *discCurve,
 TDate         startDate,
 TDate         endDate);

static double FwSurvProb
(CxTCreditCurve        *spreadCurve,
 TDate         startDate,
 TDate         endDate);


int CxCashFlowListPV(TCashFlowList * cfl,
		     TDate valueDate,
		     TCurve* discCurve,
		     CxTCreditCurve *spreadCurve,
		     double         *pv)
{
    static char routine[] = "CxCashFlowListPV";
    int         status    = FAILURE;
    int i;

    REQUIRE(cfl != NULL);
    REQUIRE(valueDate != 0);
    
    *pv = 0;
    
    for (i = 0; i < cfl->fNumItems ; i++)
      {
	double thisCouponPV;
	double survival;
	double discount;
	survival = FwSurvProb (spreadCurve,
			       valueDate,
			       cfl->fArray[i].fDate);
        discount = FwZeroPrice (discCurve,
				valueDate,
				cfl->fArray[i].fDate);
	thisCouponPV = cfl->fArray[i].fAmount * survival * discount;
	*pv +=thisCouponPV ;

      }
    
    status = SUCCESS;

 done:
    
 if(status ==FAILURE)
    GtoErrMsgFailure(routine);

    return status;

}

int CxBondFeeLegPV(TBond* bond,
		   TDate valueDate,
		   TDate startDate,
		   TCurve* discCurve,
		   CxTCreditCurve *spreadCurve,
		   CxTRecoveryCurve *recoveryCurve,
		   TBoolean isClean,
		   double         *pv)
{

  static char routine[] = "CxBondFeeLegPV";
  int         status    = FAILURE;
  TCashFlowList* cfl = NULL;
  TCashFlowList* cflWithoutPrincipal = NULL;
  TBondCFLParams* cflParams = NULL;
  TBondTradeData * bondTradeData = NULL;
  double bondAccrualOnDefault = 0.0;

  REQUIRE(valueDate != 0);
  REQUIRE(startDate !=0);
  REQUIRE(bond !=NULL);
  REQUIRE (bond->maturityDate > startDate);
  REQUIRE(startDate >= valueDate);

  cflParams = GtoStringToBondCFLParams("I,0");
  bondTradeData = GtoNewTBondTradeData(startDate);
  if (cflParams == NULL || bondTradeData == NULL) goto done;
  
  /* Bond cash flow object including the accrued interest
     (if isClean = TRUE and principal*/
  cfl = GtoBondCashFlows(bond,
			 bondTradeData,
			 isClean);
  if (cfl == NULL) goto done;
  ASSERT(cfl->fArray[0].fDate >= startDate);

  /* Bond Cash flow object without accrued interest nor principal*/
  cflWithoutPrincipal = GtoBondCashFlowInterval2(bond,
						 bondTradeData,
						 NULL,
						 cflParams);
  if(cflWithoutPrincipal == NULL) goto done;
  
  if (CxCashFlowListPV(cfl,
		       valueDate,
		       discCurve,
		       spreadCurve,
		       pv) == FAILURE) goto done;
    

  if(recoveryCurve !=NULL)
    {
      if (cflWithoutPrincipal ==NULL) goto done;
      if(CxBondFeeLegAccrualsOnDefaultPV(cflWithoutPrincipal,
					 valueDate,
					 discCurve,
					 spreadCurve,
					 recoveryCurve,
					 &bondAccrualOnDefault)==FAILURE) goto done;
    }
  *pv += bondAccrualOnDefault;
  
  status = SUCCESS;
  
 done:
  
  GtoFreeTBondTradeData(bondTradeData);
  GtoFreeCFL(cfl);
  GtoFreeCFL(cflWithoutPrincipal);
  GtoFreeBondCFLParams(cflParams);
  if (status != SUCCESS)
    {
      GtoErrMsgFailure (routine);
    }
    

  return status;
    
}


int CxBondFeeLegAccrualsOnDefaultPV
(TCashFlowList* cfl,
 TDate valueDate,
 TCurve* discCurve,
 CxTCreditCurve *spreadCurve,
 CxTRecoveryCurve* recoveryCurve,
 double *bondAccrualOnDefault)
{
  static char routine[] = "CxBondFeeLegAccrualsOnDefault";
  int         status    = FAILURE;
  int i;
    
  TDate accStartDate = valueDate;

  REQUIRE(cfl !=NULL);
  REQUIRE(valueDate !=0);
  REQUIRE(discCurve !=NULL);
  REQUIRE(spreadCurve !=NULL);
  REQUIRE(recoveryCurve !=NULL);
  *bondAccrualOnDefault = 0.0;
  
  if(cfl->fNumItems == 0)
    {
      status = SUCCESS;
      goto done;
    }
  
 
 for (i=0; i<cfl->fNumItems;i++)
   {
     double amount = cfl->fArray[i].fAmount;
     TDate accEndDate=cfl->fArray[i].fDate;
     if(accEndDate > valueDate)
       {
     double averageRecovery;
     double accrual;
     
     if(CxRecoveryCurveAverage(recoveryCurve,
			       accStartDate,
			       accEndDate,
			       &averageRecovery)==FAILURE)
       goto done;
     
    
     
  if (AccrualOnDefaultPVWithTimeLine (valueDate,
                                      accStartDate,
				      accEndDate,
				      0, /* delay */
				      amount,
				      discCurve,
				      spreadCurve,
				      NULL,
				      &accrual) != SUCCESS)
    goto done; /* failure */
  *bondAccrualOnDefault += accrual * averageRecovery;
  accStartDate=accEndDate;  
       }

   }
 
 
  
 status = SUCCESS;
 done:

 if(status ==FAILURE)
    GtoErrMsgFailure(routine);

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

static double FwSurvProb
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

/**************************************************
*** Contingent Leg 
*************************************************/
static CxTRecoveryCurve* lossCurveMake
(CxTRecoveryCurve* recoveryCurve);


/*
 * Makes a contingent leg for a vanilla bond
 *
 * Protection starts either at the beginning of startDate (protectStart=True)
 * or at the end of startDate.
 *
 * Protection ends at the end of endDate.
 *
 * Notional is the amount of notional protected.
 * Delay is the delay in payment after default measured in days.
 */
CxTContingentLeg* CxBondContingentLegMake
(TDate     startDate,
 TBond*     bond,  
 double    notional, 
 long      delay,    
 TBoolean  protectStart)
{
    static char routine[] = "CxBondContingentLegMake";
    int         status    = FAILURE;
    
    CxTContingentLeg *cl = NULL;
    
    REQUIRE(startDate !=0);
    REQUIRE(bond != NULL);
    REQUIRE (bond->maturityDate > startDate);
    REQUIRE (delay >= 0);

    cl = CxContingentLegMakeEmpty(1);
    if (cl == NULL)
        goto done; /* failure */

    /* cl->startDate is defined as giving protection from end of startDate.
       So if we want to protect on the start date, we need to move this
       date forward by one. */
    if (protectStart)
    {
        cl->startDate = startDate-1;
    }
    else
    {
        cl->startDate = startDate;
    }
    cl->dates[0]     = bond->maturityDate;
    cl->notionals[0] = notional;
    cl->payType      = CX_PROT_PAY_DEF;
    cl->payDelay     = delay;

    status = SUCCESS;

 done:

    if (status != SUCCESS)
    {
        CxContingentLegFree(cl);
        cl = NULL;
        GtoErrMsgFailure(routine);
    }

    return cl;
}


/*
 * computes the PV for a contingent leg for a vanilla CDS
 */
int CxBondContingentLegPV
(TDate             today,
 TDate             valueDate,
 TDate             startDate,
 TBond             * bond,
 double            notional,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart,
 double           *pv)
{
    static char routine[] = "CxBondContingentLegPV";
    int         status    = FAILURE;

    CxTContingentLeg *cl = NULL;
    CxTRecoveryCurve* lossCurve = NULL;

    REQUIRE(valueDate !=0);
    REQUIRE(startDate !=0);
    REQUIRE(recoveryCurve !=NULL);
    REQUIRE(discCurve !=NULL);
    REQUIRE(spreadCurve !=NULL);
    REQUIRE(bond !=NULL);
    
    lossCurve = lossCurveMake(recoveryCurve);
    
    cl = CxBondContingentLegMake (startDate, bond, notional, delay,
                                 protectStart);
    if (cl == NULL)
        goto done; /* failure */

    if (CxContingentLegPV (cl, today, valueDate, discCurve, spreadCurve,
                           lossCurve, pv) != SUCCESS)
        goto done; /* failure */

    status = SUCCESS;

 done:

    if (status != SUCCESS)
        GtoErrMsgFailure (routine);

    CxRecoveryCurveFree(lossCurve);
    CxContingentLegFree (cl);

    return status;
}

/*Takes a recovery curve and will given back a loss recovery by returning
  (1 - Recovery). This allows us to use the CDS analytics for calculating
  contingent leg integrals*/
static CxTRecoveryCurve* lossCurveMake
(CxTRecoveryCurve* recoveryCurve)
{
  static char routine[]="lossCurveMake";
  int status = FAILURE;
  
  int i;

  CxTRecoveryCurve* lossCurve = NULL;

  REQUIRE(recoveryCurve !=NULL);
  
  lossCurve = CxRecoveryCurveCopy(recoveryCurve);
  
  if (lossCurve ==NULL) goto done;
  
  for (i=0;i<recoveryCurve->numItems;i++)
    {
      lossCurve->recoveryRates[i] = 1.0 - recoveryCurve->recoveryRates[i];
    }

  status = SUCCESS;
  
 done:
  if (status == FAILURE)
    {
      CxRecoveryCurveFree(lossCurve);
      lossCurve = NULL;
      GtoErrMsgFailure (routine);
    }
  
  return lossCurve;

}


/**************************************
** Bond Price
**************************************/
int CxBondPrice
(TDate             today,
 TDate             valueDate,
 TDate             startDate,
 TBond*            bond,
 double            notional,
 long              delay,
 TCurve           *discCurve,
 CxTCreditCurve   *spreadCurve,
 CxTRecoveryCurve *recoveryCurve,
 TBoolean          protectStart,
 TBoolean              isClean,
 double           *pv)
{
 static char routine[] = "CxBondPrice";
    int         status    = FAILURE;

    double feeLegPV;
    double contingentLegPV;

    if (CxBondContingentLegPV( today,
			       valueDate,
			       startDate,
			       bond,
			       notional,
			       delay,
			      discCurve,
			      spreadCurve,
			      recoveryCurve,
			       protectStart,
			       &contingentLegPV)== FAILURE) goto done;

    if(  CxBondFeeLegPV(bond,
			valueDate,
			startDate,
			discCurve,
			spreadCurve,
			recoveryCurve,
			isClean,
			&feeLegPV) == FAILURE) goto done;

    *pv = contingentLegPV + feeLegPV * notional;

    status = SUCCESS;
    
 done:
    if (status != SUCCESS)
      GtoErrMsgFailure (routine);
    
    return status;
    
}
