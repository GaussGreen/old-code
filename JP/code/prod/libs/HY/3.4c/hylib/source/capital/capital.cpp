/*-----------------------------------------------------------------------
HEADER FILE:    capital.c

CREATED BY:     7/14/2000 Neil Yang

PURPOSE:   Hy model wrapper for 
$Header$
---------------------------------------------------------------------- 
Proprietary software, whether developed for Morgan by in-house
staff or external contractors, must not be copied or removed from Morgan
premises without the approval of a Senior Vice President and Audit.

This proprietary software has been developed strictly for Morgan's own
internal use.  Any use or misuse, intentional or otherwise which contradicts
or places this policy in jeopardy is strictly forbidden.  Also, any actions or 
inactions resulting in a breach of this goal will be dealt with harshly.

Do Not Distribute this to anyone outside the Analytics Library
Group without authorization from its management. 

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */
#include "cerror.h"
#include "kexception.h"
#include "kratecurve.h"
#include "ddmap.h"
#include "genfunc.h"
#include "instruments.h"
#include "defswap.h"
#include "macros.h"
#include "capital.h"
#include "strikes.h"   //StrikeClass, OptionContext
#include "lintrp.h"		//GtoLinInterpLongPoint1, HY3.4v

#define MAX_LENGTH 200

int HYMCapitalWrapperCheckInputs
(double*	spotPrice,	/*	1	(I)	equity spot price  */
 double*	divRefSpot,	/*		(I) ref spot price to convert from div amount to div yield	*/
 long*		divDates,	/*	2	(I) projected dividend dates*/
 double*	dividends,	/*	3	(I) projected div yields if divRefSpot==0, otherwise div amount */
 double*	dps,	    /*	4	(I) */
 long*		repoDates,	/*	5	(I) repo curve dates */
 double*	repoRates,	/*	6	(I) repo curve rates*/
 long*		swapDates,	/*	7	(I) swap curve rates*/
 double*	swapRates,	/*	8	(I) swap curve dates*/
 double*	volRefSpot,	/*		(I) ref spot price to convert from stock vol to asset vol	*/
 long*		volDates,	/*	9	(I) asset vol curve dates*/
 double*	volRates,	/*	10	(I) asset vol rates if volRefSpot==0, otherwise stock vols*/
 double*    volShift,   /*  11  (I) vol shift */
 /*
  *  instrument description 
  */
  long*		instType,			/*	12	(I) instrument type */
  double*   notional,			/*  13  (I)  */
  double*	recoveryRate,		/*	14	(I) */
  long*     instCFAccStartDates, /* 15 */ 
  long*     instCFAccEndDates,   // 16
  long*		instCFDates,		/*	17	(I) */
  double*	instCFCoupons,		/*	18	(I) */
  double*	instCFAmorts,		/*	19	(I) */
  long*		claimCFDates,		/*	20	(I) issueDate, claimPardate*/
  double*	claimCFAmounts,		/*	21	(I) issuePrice */
  /*
   *  instrument optionality
   */
  long*		exerStartDates,		/*	22	(I) */
  long*		exerEndDates,		/*	23	(I) */
  double*	exerStartStrikes,	/*	24	(I) */
  double*	exerEndStrikes,		/*	25	(I) */
  long*		optionDirections,	/*	26	(I) */
  long*     optionType,         /*  27  (I) HY_CAPITAL_CALL,HY_CAPITAL_PUT */
  double*   option_barrier,     /*  28  (I) */
  long*	    exerTypes,			/*	29	(I) american=1, european =0 */
  /*
   * Model parameter
   */
   double*	lim1,			/*	30	(I) HY Model Parameter */
   double*	lim2,			/*	31	(I) HY Model Parameter */
   double*	vollim,			/*	32	(I) HY Model Parameter */
   double*	x,				/*	33	(I) HY Model Parameter */
   double*	lim,			/*	34	(I) HY Model Parameter */ 
   double*	beta,			/*	35	(I) */
  /*
   * other values !
   */
   long*	assetProcessType,	/*	36	(I) */
   long*	ppy,				/*	37	(I) */
   long*	valueDates);		/*	38	(I) */
 

HY_EXPORT char*  HY_version(){return (char*)"HY V3.4d";}


HY_EXPORT int HYMCapitalWrapper
(double*	spotPrice,	/*	1	(I)	equity spot price  */
 double*	divRefSpot,	/*		(I) ref spot price to convert from div amount to div yield	*/
 long*		divDates,	/*	2	(I) projected dividend dates*/
 double*	dividends,	/*	3	(I) projected div yields if divRefSpot==0, otherwise div amount*/
 double*	dps,	    /*	4	(I) */
 long*		repoDates,	/*	5	(I) repo curve dates */
 double*	repoRates,	/*	6	(I) repo curve rates*/
 long*		swapDates,	/*	7	(I) swap curve rates*/
 double*	swapRates,	/*	8	(I) swap curve dates*/
 double*	volRefSpot,	/*		(I) ref spot price to convert from stock vol to asset vol	*/
 long*		volDates,	/*	9	(I) asset vol curve dates*/
 double*	volRates,	/*	10	(I) asset vol rates if volRefSpot==0, otherwise stock vols*/
 double*    volShift,   /*  11  (I) vol shift */
 /*
  *  instrument description 
  */
  long*		instType,			/*	12	(I) instrument type */
  double*   notional,			/*  13  (I)  */
  double*	recoveryRate,		/*	14	(I) */
  long*     instCFAccStartDates, /* 15 */ 
  long*     instCFAccEndDates,   // 16
  long*		instCFDates,		/*	17	(I) */
  double*	instCFCoupons,		/*	18	(I) */
  double*	instCFAmorts,		/*	19	(I) */
  long*		claimCFDates,		/*	20	(I) issueDate, claimPardate*/
  double*	claimCFAmounts,		/*	21	(I) issuePrice */
  /*
   *  instrument optionality
   */
  long*		exerStartDates,		/*	22	(I) */
  long*		exerEndDates,		/*	23	(I) */
  double*	exerStartStrikes,	/*	24	(I) */
  double*	exerEndStrikes,		/*	25	(I) */
  long*		optionDirections,	/*	26	(I) */
  long*     optionType,         /*  27  (I) HY_CAPITAL_CALL,HY_CAPITAL_PUT */
  double*   option_barrier,     /*  28  (I) */
  long*	    exerTypes,			/*	29	(I) american=1, european =0 */
  /*
   * Model parameter
   */
   double*	lim1,			/*	30	(I) HY Model Parameter */
   double*	lim2,			/*	31	(I) HY Model Parameter */
   double*	vollim,			/*	32	(I) HY Model Parameter */
   double*	x,				/*	33	(I) HY Model Parameter */
   double*	lim,			/*	34	(I) HY Model Parameter */ 
   double*	beta,			/*	35	(I) */
  /*
   * other values !
   */
   long*	assetProcessType,	/*	36	(I) */
   long*	ppy,				/*	37	(I) */
   long*	valueDates,			/*	38	(I) */
   char*	outputStrings,		/*	39  (O) Some text output */
   double*	outputNumbers)		/*	40	(O) 40ish doubles from model*/
{
    static char routine[] = "HYMCapitalWrapper";
	int status = FAILURE;
//	KRateCurve   irCurve(valueDates[2],&swapDates[1],&swapRates[1],swapDates[0],1,1)
	int i;
	KDate  *kdividentDates = 0;
	double *dividendYields = 0;
	double *assetVols = 0;
	LogNormalFunction *assetProcessFunc = 0;
	KDate  *kPayDates = 0;
	KDate  *kAccStartDates = 0;
	KDate  *kAccEndDates = 0;
	double *couponAmounts = 0;
	double *amortAmounts = 0;
	Instrument *deal = 0;
	OptionContext  *oc = 0;
	OptionContext  *ocConvert = 0;
	StrikesClass   *sc = 0;
	StrikesClass   *scSoft = 0;
	StrikesClass   *scConvert = 0;
	bool isCVOption = false;				//HY3.4v
	long  longShort, convertLongShort, exerType, exerTypeConvert;

	KRateCurve  *strike = 0;
	TDate   strikeDates[MAX_LENGTH];
	double  strikes[MAX_LENGTH];
	TDate   softStrikeDates[MAX_LENGTH];
	double  softStrikes[MAX_LENGTH];
	TDate   convertDates[MAX_LENGTH];
	double  convertStrikes[MAX_LENGTH];
	long    numStrikes=0;
	long    numSoftStrikes=0;
	long    numConvertStrikes = 0;

#ifdef DEBUG_K
	GtoErrMsgOn();
	GtoErrMsg((char*)"%s:\n",(char*) "Begining of the Kapital Wrapper");
#endif

	// check inputs
	if(HYMCapitalWrapperCheckInputs(spotPrice,
									divRefSpot,
									divDates,
									dividends,	
									dps,	    
									repoDates,	
									repoRates,
									swapDates,	
									swapRates,	
									volRefSpot,
									volDates,	
									volRates,
									volShift,   
									instType,	
									notional,     
									recoveryRate,	
									instCFAccStartDates, 
									instCFAccEndDates,   
									instCFDates,	
									instCFCoupons,
									instCFAmorts,	
									claimCFDates,	
									claimCFAmounts,	
									exerStartDates,	
									exerEndDates,	
									exerStartStrikes,	
									exerEndStrikes,		
									optionDirections,
									optionType,        
									option_barrier,     
									exerTypes,			
									lim1,			
									lim2,			
									vollim,			
									x,	
									lim,	
									beta,		
									assetProcessType,
									ppy,		
									valueDates) == FAILURE)
	{
		goto error;
	}

#ifdef DEBUG_K
	GtoErrMsg((char*)"%s:\n",(char*) "Inputs check OK");
#endif

	// set a local scope so that the code can compile
	try{
	
	kdividentDates = new KDate[divDates[0]];
	dividendYields = new double[divDates[0]];
	if(divRefSpot[1] == 0.0)
	{
		for(i = 0; i<divDates[0]; i++)
		{
			kdividentDates[i] = divDates[i+1];
			dividendYields[i] = dividends[i+1];
		}
	}
	else
	{
		for(i = 0; i<divDates[0]; i++)
		{
			kdividentDates[i] = divDates[i+1];
			dividendYields[i] = dividends[i+1]/divRefSpot[1];
		}
	}

	DDMap  divident(divDates[0], kdividentDates, dividendYields);
	KRateCurve  krepoCurve(valueDates[1],&repoDates[1],&repoRates[1],repoDates[0]);
	KRateCurve  kzeroCurve(valueDates[1],&swapDates[1],&swapRates[1],swapDates[0]);

	AssetToStockMapping1 ats(x[1],lim[1],dps[1]);

	assetVols = new double[volDates[0]];
	if(volRefSpot[1] == 0.0)
	{
		for(i = 0; i<volDates[0]; i++)
			assetVols[i]= volRates[i+1];
	}
	else
	{
		for(i = 0; i<volDates[0]; i++)
			assetVols[i]= ats.stockVolToAssetVol(volRefSpot[1],volRates[i+1]);
	}
//	KRateCurve  kvolCurve1(valueDates[1],&volDates[1],assetVols,volDates[0]);		//HY3.4v
	

  if(assetProcessType[1] ==1)
	{
		assetProcessFunc = new LogNormalFunction();
	}

//	AssetToStockMapping1 ats(x[1],lim[1],dps[1]);

	int     numDates = instCFDates[0];
	kPayDates = new KDate[numDates];
	kAccStartDates = new KDate[numDates];
	kAccEndDates = new KDate[numDates];
	couponAmounts = new double[numDates];
	amortAmounts = new double[numDates];

#ifdef DEBUG_K
	GtoErrMsg((char*)"%s:\n",(char*) "Market construction OK");
#endif

	for(i=0; i<numDates;i++)
	{
		kPayDates[i] = instCFDates[i+1];
		kAccStartDates[i] = instCFAccStartDates[i+1];
		kAccEndDates[i] = instCFAccEndDates[i+1];
		if(instType[1] == CDS || instType[1] == CDSOPT)
		{
			couponAmounts[i] = -instCFCoupons[i+1]/notional[1];
			amortAmounts[i] = instCFAmorts[i+1]/notional[1];
			
		}
		else
		{
			couponAmounts[i] = instCFCoupons[i+1]/notional[1];
			amortAmounts[i] = instCFAmorts[i+1]/notional[1];
		}

/*		if(couponAmounts[i] < 0 || amortAmounts[i] <0)
		{
			GtoErrMsg("%s: cashflow i = %d has the wrong sign.\n",routine,i);
			goto error;
		}   */
	}

	TDate maturityDate = instCFDates[instCFDates[0]];
	TDate tempDate;

	KRateCurve  strike1(valueDates[1], &maturityDate,&lim1[1],1);
	KRateCurve  strike2(valueDates[1], &maturityDate,&lim2[1],1);

	double couponAccPercent;
	if(instType[1] == CDS || instType[1] == CDSOPT)
	{
		couponAccPercent = 1.0;
	}
	else
	{
		couponAccPercent = recoveryRate[1];
	}
	CashFlowList  cfList(numDates,kPayDates,kAccStartDates,kAccEndDates,
	couponAmounts,amortAmounts,couponAccPercent,claimCFDates[1],claimCFDates[2],claimCFAmounts[1]);

	
//	return 0.0 for expired trades
	if(instType[1] == EQUITYOPT)
	{
		tempDate = exerEndDates[(int) exerEndDates[0]];
	}
	else
	{
		tempDate = maturityDate;
	}

	if(tempDate <= valueDates[1])				//HY3.4c
	{
		outputNumbers[0] = 40;

		outputNumbers[1] = 0.0;		//dirty price
		outputNumbers[2] = 0.0;		//clean price
		outputNumbers[3] = 0.0;		// delta
		outputNumbers[4] = 0.0;		//gamma
		outputNumbers[5] = 0.0;		//vega

		outputNumbers[7] = 0.0;		//accrued 
		outputNumbers[9] = 0.0;		//stock vol 
		outputNumbers[13] = 0.0;	// annuity
		
		outputNumbers[18] = 0.0;					//initial asset value
		outputNumbers[19] = dps[1];					//dps
		outputNumbers[20] = recoveryRate[1];
		outputNumbers[15] = 0.0;					//asset vol

		status = SUCCESS;
		goto error;
	}

	double period;

	switch(instType[1])
	{
	case  CDS:
		deal = new DefaultSwap(&strike1,&strike2,recoveryRate[1],cfList);
		period = (maturityDate-valueDates[1])/365.0;
		tempDate = maturityDate;
		break;

	case  CDSOPT:	//this block has been totally rewritten in HY3.3.2v
		{
			long		numStrikesOrig = exerEndDates[0] - 1;	

			TDate		exerStartDatesOrig[MAX_LENGTH];			
			TDate		exerEndDatesOrig[MAX_LENGTH];			
			double		exerStartStrikesOrig[MAX_LENGTH];		
			double		exerEndStrikesOrig[MAX_LENGTH];			
			long		optionDirectionsOrig[MAX_LENGTH];		
			long		optionTypeOrig[MAX_LENGTH];				
			double		option_barrierOrig[MAX_LENGTH];			
			long	    exerTypesOrig[MAX_LENGTH];				

			//undo add/shift of strike in Kapital
			for(i=0; i<numStrikesOrig; i++)							
			{
				exerStartDatesOrig[i] = exerStartDates[i+2];		
				exerEndDatesOrig[i] = exerEndDates[i+2];			

				exerStartStrikesOrig[i] = exerStartStrikes[i+1];	
				exerEndStrikesOrig[i] = exerEndStrikes[i+1];		
				optionDirectionsOrig[i] = optionDirections[i+1];	
				optionTypeOrig[i] = optionType[i+1];				
				option_barrierOrig[i] = option_barrier[i+1];		
				exerTypesOrig[i] = exerTypes[i+1];					
			}
		
			//sort out different strikes
			for(i=0; i<numStrikesOrig; i++)
			{
					strikeDates[numStrikes] = exerStartDatesOrig[i];		
					strikes[numStrikes] = exerStartStrikesOrig[i];			
					numStrikes++;											

					if (exerEndDatesOrig[i] == maturityDate)
					{
						strikeDates[numStrikes] = exerEndDatesOrig[i];			
					}
					else
					{
						strikeDates[numStrikes] = exerEndDatesOrig[i] - 1;			
					}
					strikes[numStrikes] = exerEndStrikesOrig[i];			
					numStrikes++;											

					longShort = optionDirectionsOrig[i];
					exerType = exerTypesOrig[i];
			}
					
			sc = new StrikesClass(strikeDates,strikes,numStrikes,"L");
			oc = new OptionContext(longShort,exerType? true:false,sc);

//			sc = new StrikesClass(&exerEndDates[1],&exerEndStrikes[1],exerEndDates[0],"S");
//			oc = new OptionContext(optionDirections[1],exerTypes[1]? true:false,sc);

			deal = new DefaultSwap(&strike1,&strike2,recoveryRate[1],cfList,oc);
			period = (maturityDate-valueDates[1])/365.0;
			tempDate = maturityDate;
			break;
		}

	case  BOND:

		deal = new Bond(&strike1,&strike2,recoveryRate[1],cfList);
		period = (maturityDate-valueDates[1])/365.0;
		tempDate = maturityDate;
		break;

	case  BONDOPT:		//this block has been totally rewritten in HY3.3.2v
		{
			long		numStrikesOrig = exerEndDates[0] - 1;	

			TDate		exerStartDatesOrig[MAX_LENGTH];			
			TDate		exerEndDatesOrig[MAX_LENGTH];			
			double		exerStartStrikesOrig[MAX_LENGTH];		
			double		exerEndStrikesOrig[MAX_LENGTH];			
			long		optionDirectionsOrig[MAX_LENGTH];		
			long		optionTypeOrig[MAX_LENGTH];				
			double		option_barrierOrig[MAX_LENGTH];			
			long	    exerTypesOrig[MAX_LENGTH];				

			//undo add/shift of strike in Kapital
			for(i=0; i<numStrikesOrig; i++)							
			{
				exerStartDatesOrig[i] = exerStartDates[i+2];		
				exerEndDatesOrig[i] = exerEndDates[i+2];			

				exerStartStrikesOrig[i] = exerStartStrikes[i+1];	
				exerEndStrikesOrig[i] = exerEndStrikes[i+1];		
				optionDirectionsOrig[i] = optionDirections[i+1];	
				optionTypeOrig[i] = optionType[i+1];				
				option_barrierOrig[i] = option_barrier[i+1];		
				exerTypesOrig[i] = exerTypes[i+1];					
			}
				
			//sort out different strikes
			for(i=0; i<numStrikesOrig; i++)
			{
				switch(optionTypeOrig[i])
				{
				case HY_CAPITAL_CLEAN:
					strikeDates[numStrikes] = exerStartDatesOrig[i];		
					strikes[numStrikes] = exerStartStrikesOrig[i];			
					numStrikes++;											

					if (exerEndDatesOrig[i] == maturityDate)
					{
						strikeDates[numStrikes] = exerEndDatesOrig[i];			
					}
					else
					{
						strikeDates[numStrikes] = exerEndDatesOrig[i] - 1;			
					}
					strikes[numStrikes] = exerEndStrikesOrig[i];			
					numStrikes++;											

					longShort = optionDirectionsOrig[i];
					exerType = exerTypesOrig[i];
					break;
				case HY_CAPITAL_SOFTCALL:
					softStrikeDates[numSoftStrikes] = exerStartDatesOrig[i];	
					softStrikes[numSoftStrikes] = exerStartStrikesOrig[i];		
					numSoftStrikes++;											

					if (exerEndDatesOrig[i] == maturityDate)
					{
						softStrikeDates[numSoftStrikes] = exerEndDatesOrig[i];			
					}
					else
					{
						softStrikeDates[numSoftStrikes] = exerEndDatesOrig[i] - 1;			
					}
					softStrikes[numSoftStrikes] = exerEndStrikesOrig[i];		
					numSoftStrikes++;										

					break;
				case HY_CAPITAL_SHARES:
					convertDates[numConvertStrikes] = exerStartDatesOrig[i];			
					convertStrikes[numConvertStrikes] = exerStartStrikesOrig[i];		
					numConvertStrikes++;												

					if (exerEndDatesOrig[i] == maturityDate)
					{
						convertDates[numConvertStrikes] = exerEndDatesOrig[i];			
					}
					else
					{
						convertDates[numConvertStrikes] = exerEndDatesOrig[i] - 1;			
					}
					convertStrikes[numConvertStrikes] = exerEndStrikesOrig[i];		
					numConvertStrikes++;												

					convertLongShort = optionDirectionsOrig[i];
					exerTypeConvert = exerTypesOrig[i];
					break;
				default:
					GtoErrMsg("Wrong strike type %d.\n",optionTypeOrig[i]);
					goto error;
				}
			}

			if(numStrikes != 0)
			{
				sc = new StrikesClass(strikeDates,strikes,numStrikes,"L");
			}
			if(numSoftStrikes != 0)
			{
				scSoft = new StrikesClass(softStrikeDates,softStrikes,numSoftStrikes,"L");
			}
			if(numConvertStrikes != 0)
			{
				scConvert = new StrikesClass(convertDates,convertStrikes,numConvertStrikes,"L");
			}
			if(numStrikes != 0)
			{
				oc = new OptionContext(longShort,exerType? true:false,sc,scSoft);
			}
			if(numConvertStrikes != 0)
			{
				ocConvert = new OptionContext(convertLongShort,exerTypeConvert? true:false,scConvert);
				isCVOption = true;			//HY3.4v
			}
			deal = new Bond(&strike1,&strike2,recoveryRate[1],cfList,0,oc, ocConvert);
			period = (maturityDate-valueDates[1])/365.0;
			tempDate = maturityDate;
			break;
		}
	case  EQUITYOPT:
		strike = new KRateCurve(valueDates[1], &exerEndDates[1],&exerEndStrikes[1],exerEndDates[0]);
		if(optionType[1] == HY_CAPITAL_CALL)
		{
			deal = new CallOptionT(strike, exerTypes[1]? true:false);
		}
		else
		{
			deal = new PutOptionT(strike, exerTypes[1]? true:false);
		}
		tempDate = exerEndDates[(int) exerEndDates[0]];
		period = (tempDate-valueDates[1])/365.0;

	    beta[1] = 0.0;		//HY3.3.1v
		break;
	}
	
	long volSize = volDates[0];			//HY3.4v
	double volRate;						//HY3.4v
	if(tempDate <= volDates[1])
	{
		volRate = assetVols[0];				//HY3.4v
	}
	else if(tempDate >= volDates[volSize])			//HY3.4v
	{
		volRate = assetVols[volSize - 1];				//HY3.4v
	}
	else
	{
		int Gto = GtoLinInterpLongPoint1(&volDates[1], sizeof(long), volDates[0],assetVols, sizeof(double),
											tempDate, NULL, &volRate);		//HY3.4v
	}
	
	//	double volRate = volRateTemp;												//HY3.4v
//	double volRate = kvolCurve1.get_rate(tempDate);

#ifdef DEBUG_K
	
	GtoErrMsg((char*)"status= %lf \n",volRate);
#endif
	
	KRateCurve  kvolCurve(valueDates[1],&tempDate,&volRate,1);

// HY3.4v
	TDate volDate_opt = tempDate;
	double assetVol_opt = volRate - volShift[1];
	if(assetVol_opt <= 0.0)
	{
		GtoErrMsg((char*)"%s: assetVol_opt %d is zero or negative.\n",routine,assetVol_opt);
		goto error;
	}

	KRateCurve  kvolCurve_opt(valueDates[1],&volDate_opt,&assetVol_opt,1);

// end HY3.4v

#ifdef DEBUG_K
	
	GtoErrMsg((char*)"%s:\n",(char*) "Cashflow construction OK");
#endif
	

	long ppy2 = (long)(MAX(ppy[1],ppy[1]*sqrt(period))/period);
//	long ppy2 = ppy[1];


	KValarray<double> p = GeneralPricer(spotPrice[1],
					    &divident,
					    &krepoCurve,
					    &kzeroCurve,
					    &kvolCurve,
					    &kvolCurve_opt,		//HY3.4v
					    deal,
					    &ats,
					    assetProcessFunc,
					    ppy2,
					    beta[1],
					    valueDates[2],
						isCVOption);
	outputNumbers[0] = 40;
	
	outputNumbers[1] = p[0]*notional[1];		//dirty price
	outputNumbers[2] = (p[0]-p[6])*notional[1]; //clean price
	outputNumbers[3] = p[1]*notional[1];		// delta
	outputNumbers[4] = p[2]*notional[1];		//gamma
	outputNumbers[5] = p[3]*notional[1];		//vega

	outputNumbers[7] = p[6]*notional[1];		//accrued 
	outputNumbers[9] = p[5];					//stock vol 
	outputNumbers[13] = p[4];					// annuity
	
	outputNumbers[18] = p[7];					//initial asset value
	outputNumbers[19] = dps[1];					//dps
	outputNumbers[20] = recoveryRate[1];
	outputNumbers[15] = kvolCurve.get_rate(maturityDate);   //asset vol

	
	}

	catch(KException &e)
	{
		GtoErrMsg((char*)"%s:",(char*) e.what());
	}

	status = SUCCESS;

#ifdef DEBUG_K
	
	GtoErrMsg((char*)"%s:\n",(char*) "Pricing OK");
#endif

error:

	delete[] kdividentDates;
	delete[] dividendYields;
	delete[] assetVols;
	delete assetProcessFunc;
	delete[] kPayDates;
	delete[] kAccStartDates;
	delete[] kAccEndDates;
	delete[] couponAmounts;
	delete[] amortAmounts;
	delete deal;
	delete oc;
	delete ocConvert;
	delete strike;

#ifdef DEBUG_K
	
	GtoErrMsg((char*)"%s status= %d \n",(char*) "End of the Kapital Wrapper",status);
#endif

	return status;
}



int HYMCapitalWrapperCheckInputs
(double*	spotPrice,	/*	1	(I)	equity spot price  */
 double*	divRefSpot,	/*		(I) ref spot price to convert from div amount to div yield	*/
 long*		divDates,	/*	2	(I) projected dividend dates*/
 double*	dividends,	/*	3	(I) projected dividend yields or amount*/
 double*	dps,	    /*	4	(I) */
 long*		repoDates,	/*	5	(I) repo curve dates */
 double*	repoRates,	/*	6	(I) repo curve rates*/
 long*		swapDates,	/*	7	(I) swap curve rates*/
 double*	swapRates,	/*	8	(I) swap curve dates*/
 double*	volRefSpot,	/*		(I) ref spot price to convert from stock vol to asset vol	*/
 long*		volDates,	/*	9	(I) asset vol curve dates*/
 double*	volRates,	/*	10	(I) asset vol curve rates*/
 double*    volShift,   /*  11  (I) vol shift */
 /*
  *  instrument description 
  */
  long*		instType,		/*	12	(I) instrument type */
  double*   notional,       /*  13  (I)  */
  double*	recoveryRate,	/*	14	(I) */
  long*     instCFAccStartDates, /* 15 */ 
  long*     instCFAccEndDates,   // 16
  long*		instCFDates,	    /*	17	(I) */
  double*	instCFCoupons,	    /*	18	(I) */
  double*	instCFAmorts,	    /*	19	(I) */
  long*		claimCFDates,	    /*	20	(I) issueDate, claimPardate*/
  double*	claimCFAmounts,	    /*	21	(I) issuePrice */
  /*
   *  instrument optionality
   */
  long*		exerStartDates,		/*	22	(I) */
  long*		exerEndDates,		/*	23	(I) */
  double*	exerStartStrikes,	/*	24	(I) */
  double*	exerEndStrikes,		/*	25	(I) */
  long*		optionDirections,	/*	26	(I) */
  long*     optionType,         /*  27  (I) HY_CAPITAL_CALL,HY_CAPITAL_PUT */
  double*   option_barrier,     /*  28  (I) */
  long*	    exerTypes,			/*	29	(I) american=1, european =0 */
  /*
   * Model parameter
   */
   double*	lim1,			    /*	30	(I) HY Model Parameter */
   double*	lim2,			    /*	31	(I) HY Model Parameter */
   double*	vollim,			    /*	32	(I) HY Model Parameter */
   double*	x,		            /*	33	(I) HY Model Parameter */
   double*	lim,	            /*	34	(I) HY Model Parameter */ 
   double*	beta,			    /*	35	(I) */
  /*
   * other values !
   */
   long*	assetProcessType,	/*	36	(I) */
   long*	ppy,			    /*	37	(I) */
   long*	valueDates)		    /*	38	(I) */
{
	static char routine[] = "HYMCapitalWrapperCheckInputs";
	int status = FAILURE;
	int i,numDates;
	
	if(spotPrice[0] != 1)
	{
		GtoErrMsg((char*) "%s: spotPrice is not a scaler.\n",routine);
		goto error;
	}

	if(divRefSpot[0] != 1)
	{
		GtoErrMsg((char*) "%s: divRefSpot is not a scaler.\n",routine);
		goto error;
	}

	if(divDates[0] != dividends[0])
	{
		GtoErrMsg((char*)"%s: number of divDates %d != number of dividends %d.\n",routine,
			divDates[0],dividends[0]);
		goto error;
	}
	
	if(dps[0] != 1)
	{
		GtoErrMsg((char*)"%s: dps %d is not a scaler.\n",routine,dps[0]);
		goto error;
	}
	
	if(repoDates[0] != repoRates[0])
	{
		GtoErrMsg((char*)"%s: number of repoDates %d != number of repoRates %d.\n",routine,
			repoDates[0],repoRates[0]);
		goto error;
	}

	if(swapDates[0] != swapRates[0])
	{
		GtoErrMsg((char*)"%s: number of swapDates %d != number of swapRates %d.\n",routine,
			swapDates[0],swapRates[0]);
		goto error;
	}

	
	if(volDates[0] != volRates[0])
	{
		GtoErrMsg((char*)"%s: number of volDates %d != number of volRates %d.\n",routine,
			volDates[0],volRates[0]);
		goto error;
	}

	
	if(volRefSpot[0] != 1)
	{
		GtoErrMsg((char*)"%s: volRefSpot %d is not a scaler.\n",routine);
		goto error;
	}

	if(volShift[0] != 1)
	{
		GtoErrMsg((char*)"%s: volShift %d is not a scaler.\n",routine,volShift[0]);
		goto error;
	}


	if(instType[0] != 1)
	{
		GtoErrMsg((char*)"%s: instType %d is not a scaler.\n",routine,instType[0]);
		goto error;
	}

	
	if(notional[0] != 1)
	{
		GtoErrMsg((char*)"%s: notional %d is not a scaler.\n",routine,notional[0]);
		goto error;
	}
	
	if(recoveryRate[0] != 1)
	{
		GtoErrMsg((char*)"%s: recoveryRate %d is not a scaler.\n",routine,recoveryRate[0]);
		goto error;
	}

	
	if(instCFAccStartDates[0] != instCFAccEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of instCFAccStartDates %d != number of instCFAccEndDates %d.\n",routine,
			instCFAccStartDates[0],instCFAccEndDates[0]);
		goto error;
	}

	if(instCFDates[0] != instCFAccEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of instCFDates %d != number of instCFAccEndDates %d.\n",routine,
			instCFDates[0],instCFAccEndDates[0]);
		goto error;
	}

	
	if(instCFCoupons[0] != instCFAccEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of instCFCoupons %d != number of instCFAccEndDates %d.\n",routine,
			instCFCoupons[0],instCFAccEndDates[0]);
		goto error;
	}
	
	if(instCFAmorts[0] != instCFAccEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of instCFAmorts %d != number of instCFAccEndDates %d.\n",routine,
			instCFAmorts[0],instCFAccEndDates[0]);
		goto error;
	}

	
	if(claimCFDates[0] != claimCFAmounts[0])
	{
		GtoErrMsg((char*)"%s: number of claimCFDates %d != number of claimCFAmounts %d.\n",routine,
			claimCFDates[0],claimCFAmounts[0]);
		goto error;
	}
	

	if(exerStartDates[0] != exerEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of exerStartDates %d != number of exerEndDates %d.\n",routine,
			exerStartDates[0],exerEndDates[0]);
		goto error;
	}
	
	if(exerStartStrikes[0] != exerEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of exerStartStrikes %d != number of exerEndDates %d.\n",routine,
			exerStartStrikes[0],exerEndDates[0]);
		goto error;
	}
	
	if(exerEndStrikes[0] != exerEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of exerEndStrikes %d != number of exerEndDates %d.\n",routine,
			exerEndStrikes[0],exerEndDates[0]);
		goto error;
	}

	if(optionDirections[0] != exerEndDates[0])
	{
		GtoErrMsg((char*)"%s: number of optionDirections %d != number of exerEndDates %d.\n",routine,
			optionDirections[0],exerEndDates[0]);
		goto error;
	}

	if(optionType[0] != optionDirections[0])
	{
		GtoErrMsg((char*)"%s: number of optionDirections %d != number of optionType %d.\n",routine,
			optionDirections[0],optionType[0]);
		goto error;
	}

	if(option_barrier[0] != optionType[0])
	{
		GtoErrMsg((char*)"%s: number of option_barrier %d != number of optionType %d.\n",routine,
			option_barrier[0],optionType[0]);
		goto error;
	}

	if(exerTypes[0] != option_barrier[0])
	{
		GtoErrMsg((char*)"%s: number of option_barrier %d != number of exerTypes %d.\n",routine,
			option_barrier[0],exerTypes[0]);
		goto error;
	}

	if(lim1[0] != 1)
	{
		GtoErrMsg((char*)"%s: lim1 %d is not a scaler.\n",routine,lim1[0]);
		goto error;
	}

	if(lim2[0] != 1)
	{
		GtoErrMsg((char*)"%s: lim2 %d is not a scaler.\n",routine,lim2[0]);
		goto error;
	}
	
	if(vollim[0] != 1)
	{
		GtoErrMsg((char*)"%s: vollim %d is not a scaler.\n",routine,vollim[0]);
		goto error;
	}
	
	if(x[0] != 1)
	{
		GtoErrMsg((char*)"%s: x %d is not a scaler.\n",routine,x[0]);
		goto error;
	}

	if(lim[0] != 1)
	{
		GtoErrMsg((char*)"%s: lim %d is not a scaler.\n",routine,lim[0]);
		goto error;
	}

	if(beta[0] != 1)
	{
		GtoErrMsg((char*)"%s: beta %d is not a scaler.\n",routine,beta[0]);
		goto error;
	}
	
	if(assetProcessType[0] != 1)
	{
		GtoErrMsg((char*)"%s: assetProcessType %d is not a scaler.\n",routine,assetProcessType[0]);
		goto error;
	}

	
	if(ppy[0] != 1)
	{
		GtoErrMsg((char*)"%s: ppy %d is not a scaler.\n",routine,ppy[0]);
		goto error;
	}

	if(valueDates[0] != 2)
	{
		GtoErrMsg((char*)"%s: valueDates %d is not 2 elements.\n",routine,valueDates[0]);
		goto error;
	}

	numDates = instCFCoupons[0];
	for(i=0; i<numDates;i++)
	{
		double tempCoupon;
		double tempAmort;

		if(instType[1] == CDS || instType[1] == CDSOPT)
		{
			tempCoupon = -instCFCoupons[i+1]/notional[1];
			tempAmort = instCFAmorts[i+1]/notional[1];
			
		}
		else
		{
			tempCoupon = instCFCoupons[i+1]/notional[1];
			tempAmort = instCFAmorts[i+1]/notional[1];
		}

		if(tempCoupon < 0 || tempAmort <0)
		{
			GtoErrMsg((char*)"%s: cashflow i = %d has the wrong sign.\n",routine,i);
			goto error;
		}
	}

	status = SUCCESS;

error:

	return status;

}
