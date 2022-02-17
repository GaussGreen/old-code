/*-----------------------------------------------------------------------
  SOURCE FILE:    capital.h
  
  CREATED BY:     Neil Yang  July 2000
  
  PURPOSE:        header file to test interface to HY model
  
   ---------------------------------------------------------------------- */
#ifndef HY_CAPITAL_H
#define HY_CAPITAL_H

#ifdef __cplusplus
extern "C" {
#endif


/* Option Types */
#define LONG			1
#define SHORT			-1

#define HY_CAPITAL_CALL      1
#define HY_CAPITAL_PUT       2
#define HY_CAPITAL_CLEAN     1
#define HY_CAPITAL_SHARES    2
#define HY_CAPITAL_SOFTCALL  3

/* Instrument Types */
#define CDS				1
#define BOND			2
#define EQUITYOPT		3
#define BONDOPT			4
#define CDSOPT			5
#define	EQUITY			6
#define CASH			7

#ifdef _WIN32
   #define  HY_EXPORT    _declspec( dllexport )
#else
   #define  HY_EXPORT 
#endif

HY_EXPORT char*  HY_version();

HY_EXPORT int HYMCapitalWrapper
(double*	spotPrice,	/*	1	(I)	equity spot price  */
 double*	divRefSpot,	/*		(I) ref spot price to convert from div amount to div yield	*/
 long*		divDates,	/*	2	(I) projected dividend dates*/
 double*	dividends,	/*	3	(I) projected div yields if divRefSpot==0, otherwise div amount  */
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
  long*		instType,		/*	12	(I) instrument type */
  double*   notional,       /*  13  (I)  */
  double*	recoveryRate,	/*	14	(I) */
  long*     instCFAccStartDates, /* 15 */ 
  long*     instCFAccEndDates,   // 16
  long*		instCFDates,	/*	17	(I) */
  double*	instCFCoupons,	/*	18	(I) */
  double*	instCFAmorts,	/*	19	(I) */
  long*		claimCFDates,	/*	20	(I) issueDate, claimPardate*/
  double*	claimCFAmounts,	/*	21	(I) issuePrice */
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
   double*	x,		/*	33	(I) HY Model Parameter */
   double*	lim,	/*	34	(I) HY Model Parameter */ 
   double*	beta,			/*	35	(I) */
  /*
   * other values !
   */
   long*	assetProcessType,	/*	36	(I) */
   long*	ppy,			/*	37	(I) */
   long*	valueDates,		/*	38	(I) */
   char*	outputStrings,	/*	39  (O) Some text output */
   double*	outputNumbers); /*	40	(O) 40ish doubles from model*/
	/*
	 *  Output description
	 *
	 *
	 *
	 *
	 *
	 *
	 *
	 *
	 */


#ifdef __cplusplus
}
#endif

#endif








