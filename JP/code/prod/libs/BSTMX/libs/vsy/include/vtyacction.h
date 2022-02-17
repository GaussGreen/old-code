/***************************************************************
   HEADER FILE:		vtyacction.h 

   CREATED BY:      	David Liu,	06/15/2000

   PURPOSE:		Swaption with yacc formula support.
 ***************************************************************/
#ifndef BS_YACCTION_H
#define BS_YACCTION_H

#include "bastypes.h"


#ifdef  __cplusplus
extern "C" {
#endif


GTO_EXPORT(int)
DVTreeYacctionL(
					//
					// --- INSTRUMENT

					/* --- Receive Leg 	       */
    TDate  *resetDatesRecL,		/* (I) Reset dates             */
    TDate  *resetEffDatesRecL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesRecL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesRecL,		/* (I) Accrual end dates       */
    TDate  *payDatesRecL,		/* (I) Payment dates           */
    double *notionalsRecL,		/* (I) Notionals               */
    char   *formulaRecL,		/* (I) Payment formula	       */
    
					/* --- Pay Leg 		       */
    TDate  *resetDatesPayL,		/* (I) Reset dates             */
    TDate  *resetEffDatesPayL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesPayL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesPayL,		/* (I) Accrual end dates       */
    TDate  *payDatesPayL,		/* (I) Payment dates           */
    double *notionalsPayL,		/* (I) Notionals               */
    char   *formulaPayL,		/* (I) Payment formula	       */

					/* --- Receive & pay conventions  */
    char   *payDayCountStrsL,		/* (I) Dcc for payment legs  (X2) */  
    char   *stubConvsL,			/* (I) Stub conventions      (X2) */  

					/* --- Receive, pay indices    */
    double *indexWeightsL,		/* (I) Index weights (0=fixed) */
    double *indexRateSpreadsL,		/* (I) Index rate spreads      */
    char   *indexFreqsL,		/* (I) Index frequency         */
    char   *indexMaturitiesL,		/* (I) Maturity of index       */
    char   *indexDayCountStrsL,		/* (I) Dcc string for index    */  
    char   *indexCurvesL,		/* (I) Index zero curve  (X2)  */  
    char   *discIndexCurvesL,		/* (I) Discounting for each leg (X2)*/

					/* --- Principle cash flows 	   */
    TDate  *cashDatesL,			/* (I) Cash flow Dates 		   */
    double *cashAmountsL,		/* (I) Cash amounts  		   */
    char   *discCashCurvesL,		/* (I) Discount curves		   */

					/* --- Option schedule 	           */
    char   *optTypeL,			/* (I) Call, Put, T-Call, T-Put    */
    TDate  *optNotifDatesL,		/* (I) Option notification dates   */
    TDate  *optSettleDatesL,		/* (I) Option settlement Dates     */
    TDate  *optStrikePayDatesL,		/* (I) Option strike payment Dates */
    double *optStrikesL,		/* (I) Option strikes  		   */
    double *optNotionalsL,		/* (I) Notionals for option strikes*/
    char   *discOptCurveL,		/* (I) Discounting for the option  */

					//
					// ENVIRONMENT + MODEL
					//

					/* --- Zero curvs 	       */
    TDate  *baseDatesL,			/* (I) Today and 4 value dates(X5)*/
    char   *zeroFrequenciesL,		/* (I) 'A'nual or 'S'emi (X4) */
    char   *zeroDCCsL,			/* (I) ACT/365F (X4)*/
                                
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
    char   *resetBankIndsL,		/* (I) Rate index 'R'ec, 'P'ay  */
    TDate  *resetBankDatesL,		/* (I) Reset dates	        */
    double *resetBankRatesL,		/* (I) Reset Rates	        */

    int	   *debugLevelL,		/* (I) 0=no debug, >0 debug  */

    double *outputsL); 			/* (O) outputs
					   [1] option Price
					   [2] Underlying swap price
					   [3] Fugit
					   [4] Exercise probability */



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

