// capitalwrap.cpp : Defines the entry point for the console application.
//

/*------------------------------------------------------------------
C FILE:         capital_wrap.c

CREATED BY:     Neil Yang - July 2000

PURPOSE:        DR Wrapper to price credit default swap
---------------------------------------------------------------------- */
//#include "kwrapcds.h"

//#include "stdafx.h"
#include "stdio.h"
#include "iostream.h"
#include "string.h"
#include "stdlib.h"
#include "capital.h"
#include "cerror.h"

/*----------------------------------------------------------------------------
FUNCTION:       main
                                                                                                                                                                                                                                                                       
CREATED BY:     Neil Yang - Feb 2000

DESCRIPTION:    run valuation of credit default swap
----------------------------------------------------------------------------*/
int maintemp();
int main()
{
  int i;
  for(i=0; i<2 ;i++)
  {
    printf("#Run = %d \n",i);
    maintemp();
  }

  return -1;
}


int maintemp()
{
    int status = -1;
    
       
	double spotPrice[] = {1,7.3125};
	long   valueDateDr = 20000317;
	long   valueDate = 145957;
//	DrDateToTDate(valueDateDr, &valueDate);
	long   maturityDateDr = 20021019;
	long   maturityDate = valueDate+184;
//	DrDateToTDate(maturityDateDr,&maturityDate);
	long   divDates[] = {1,valueDate+3};
	double divRates[] = {1,0.0};
	long   repoDates[] = {2, 145988, 153262};
	double repoRates[] = {2, 0.005, 0.005};
	long   swapDates[] = {1, maturityDate};
	double swapRates[] = {1, 0.07};
	long   volDates[] = {2, maturityDate+356,maturityDate+780};
	double volRates[] = {2, 0.29,0.30};
	double volShift[] = {1, 0};

	long instType[] = {1, BOND};
	double  notional[] = {1,-100};
	double recoveryRate[] = {1, 0.0};
	long cfStartDates[] = {1, valueDate};
	long cfEndDates[] = {1, maturityDate};
	long cfDates[] = {1, maturityDate};
	double cfCoupon[] = {1,-10};
	double cfAmort[] = {1, -100};
	long claimDates[] = {2,valueDate, maturityDate};
	double claimAmounts[] = {2,1,1};

	long exerStartDates[] = {1, 145919};
	long exerEndDates[] = {1, 146479};
	double exerStartStrikes[] = {1, 1.0};
	double exerEndStrikes[] = {1,1.0};
	long   optionDir[] = {1,-1};
	long opotionType[] = {1, HY_CAPITAL_CALL};
	double optionBarrier[] = {1,100.0};
	long exerType[] = {1,1};

	double lim1[] = {1, 0.75};
	double lim2[] = {1,1};
	double vollim[] = {1, 0};
	double x[] = {1,1};
	double lim[] = {1,0.75};
	double dps[] = {1, 9.07};
	double beta[] = {1, 0.25};


	
	long assetProcessType[] = {1,1};
	long ppy[] = {1,100};
	long ValueDates[] = {2, 145957,145957};
	char outChar[10];
	double outputs[40];


	
   
	HYMCapitalWrapper(spotPrice,
					  divDates,
					  divRates,
					  dps,
					  repoDates,
					  repoRates,
					  swapDates,
					  swapRates,
					  volDates,
					  volRates,
					  volShift,
					  instType,
					  notional,
					  recoveryRate,
					  cfStartDates,
					  cfEndDates,
					  cfDates,
					  cfCoupon,
					  cfAmort,
					  claimDates,
					  claimAmounts,
					  exerStartDates,
					  exerEndDates,
					  exerStartStrikes,
					  exerEndStrikes,
					  optionDir,
					  opotionType,
					  optionBarrier,
					  exerType,
					  lim1,
					  lim2,
					  vollim,
					  x,
					  lim,
					  beta,
					  assetProcessType,
					  ppy,
					  ValueDates,
					  outChar,
					  outputs);

	printf("price = %f\n, delta = %f\n,gamma = %f\n,vega = %f\n,accr=%f\n,annuity = %f\n",
		outputs[1],outputs[3],outputs[4],outputs[5],outputs[7],outputs[13]);
    
        status = 0;


  

    return(status);
}


						
