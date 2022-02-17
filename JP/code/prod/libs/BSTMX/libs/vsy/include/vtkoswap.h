/***************************************************************
   HEADER FILE:		vtkoswapl.h 

   CREATED BY:      	David Liu,	02/23/2000

   PURPOSE:		Basis knock-out swap based on basis tree.
 ***************************************************************/
#ifndef BS_KOSWAP_H
#define BS_KOSWAP_H

#include "bastypes.h"


#ifdef  __cplusplus
extern "C" {
#endif


GTO_EXPORT(int)
DVTreeKnockOutBasisSwapL(
					//
					// --- INSTRUMENT
					//

					/* --- Receive Leg 	       */
    TDate  *resetDatesRecL,		/* (I) Reset dates             */
    TDate  *resetEffDatesRecL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesRecL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesRecL,		/* (I) Accrual end dates       */
    TDate  *payDatesRecL,		/* (I) Payment dates           */
	double *spreadsRecL,        /* (I) Spreads                 */
    double *notionalsRecL,		/* (I) Notionals               */
    
					/* --- Pay Leg 		       */
    TDate  *resetDatesPayL,		/* (I) Reset dates             */
    TDate  *resetEffDatesPayL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesPayL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesPayL,		/* (I) Accrual end dates       */
    TDate  *payDatesPayL,		/* (I) Payment dates           */
	double *spreadsPayL,        /* (I) Spreads                 */
    double *notionalsPayL,		/* (I) Notionals               */

					/* --- Receive & pay legs      */
    char   *payDayCountStrsL,		/* (I) Dcc for payment legs  (X2) */  
    char   *stubConvsL,			/* (I) Stub conventions      (X2) */  

					/* --- Index of receive, pay, 
					       & knock-out	       */
    double *indexWeightsL,		/* (I) Index weights (0=fixed) */
    double *indexRateSpreadsL,		/* (I) Index rate spreads      */
    char   *indexFreqsL,		/* (I) Index frequency         */
    char   *indexMaturitiesL,		/* (I) Maturity of index       */
    char   *indexDayCountStrsL,		/* (I) Dcc string for index    */  
    char   *indexCurvesL,		/* (I) Index zero curve  (X3)  */  
    char   *discIndexCurvesL,		/* (I) discounting for each leg (X3)*/

					/* --- Knock-out 	       */
    char   *kIOTypesL,			/* (I) ko type, ko window, smooth */
    TDate  *kIODatesL,			/* (I) Knock-Out dates */
    TDate  *kIOEffDatesL,		/* (I) Knock-Out Effective Dates */
    TDate  *kIOSettleDatesL,		/* (I) Knock-Out settlement Dates */
    double *kIOLoBarrierRatesL,		/* (I) KO Low barrier definition */
    double *kIOHiBarrierRatesL,		/* (I) KO High barrier definition*/
    double *kIORebateAmountsL,		/* (I) KO Rebate amounts  */

					//
					// ENVIRONMENT + MODEL
					//

					/* --- Zero curvs 	       */
    TDate  *baseDatesL,			/* (I) Today and 4 value dates(X5)*/
    char   *zeroFrequenciesL,		/* (I) 'A'nual or 'S'emi (X4) */
    char   *zeroDCCsL,			/* (I) ACT/360, ACT/365F (X4)*/
                                
    TDate  *zeroDiffDatesL,		/* (I) Diffuse zero curve dates */  
    double *zeroDiffRatesL,		/* (I) Diffuse zero curve rates */  
                                                                  
    TDate  *zero1DatesL,		/* (I) Index zero curve dates 1 */  
    double *zero1RatesL,		/* (I) Index zero curve rates 1 */  
                                                                  
    TDate  *zero2DatesL,		/* (I) Index zero curve dates 2 */  
    double *zero2RatesL,		/* (I) Index zero curve rates 2 */  
                                                                    
    TDate  *zeroBDatesL,		/* (I) Basis curve dates        */  
    double *zeroBRatesL,		/* (I) Basis curve rates        */  


    double *zeroBInfoL,			/* (I) Extra param for basis crv (X4)
					 * [1] Ref basis Libor index
					 * [2] Ref basis Disc  index
					 * [3] Basis type (Spread, Percentage)
					 * [4] Delay shift (intvl) 
					 * [5] Basis DCC 0, 3, 5
					 * [6] Libor DCC 0, 3, 5 */


					/* --- Volatility curvs        */
    TDate  *irVolDatesL,		/* (I) IR Volatility dates */
    TDate  *irVolMatsL,			/* (I) IR Volatility underlying mats */
    int    *irVolFreqL,			/* (I) IR Volatility frequencies */
    double *irVolsL,			/* (I) IR Vols (base or spot)    */

    TDate  *basisVolDatesL,		/* (I) Basis volatility dates */
    double *basisVolsL,			/* (I) Basis spot vols        */
    char   *volTypeL,			/* (I) Vol type ('N'orm, 'L'ognormal)*/

					/* --- Model Parameters		*/
    double *irMrParamL,			/* (I) IR mr parameters		*/
    double *irSmileParamL,		/* (I) IR smile parameters	*/
    double *bsMrParamL,			/* (I) Basis mr parameters	*/
    double *bsSmileParamL,		/* (I) Basis smile parameters	*/
    double *corrIrBsL,			/* (I) IR/Basis corelation	*/

					/* --- Reset bank		*/
    char   *resetBankIndsL,		/* (I) Rate index 'R', 'P', and 'K' */
    TDate  *resetBankDatesL,		/* (I) Reset dates	       */
    double *resetBankRatesL,		/* (I) Reset Rates	       */

    int	   *debugLevelL,		/* (I) 0=no debug, >0 debug  */

    double *outputsL); 			/* (O) outputs
					   [1] KO Price
					   [2] Swap Price
					   [3] Fugit
					   [4] Knock probability */



/** Wrapper function for initializing error logging into global buffer.
 */
GTO_EXPORT(void)	DVErrorLogInit (void);
 
 
/**  Wrapper function for returning the size of error buffer
 */
GTO_EXPORT(long)	DVErrorBufferSize(void);
 
 
/** Wrapper function for retrieving error messages from global buffer.
 *
 *  Returns FAILURE if outBufLen is less than the length
 *  of the accumulated error message string
 */
GTO_EXPORT(int)		DVRetrieveErrorMessages(
	long    outBufLen,      /* (I) length of the input buffer. */
	char    *outBuffer);	/* (I/O) output buffer. */



#ifdef  __cplusplus
};
#endif

#endif



