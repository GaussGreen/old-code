/**************************************************************
//
//	Yacction
//
//	D. Liu 
//
**************************************************************/
#include "kstdinc.h"
#include "kutilios.h"
#include "kmodpar.h"
#include "krate.h"
#include "kstlutil.h"		// Conversion between ptrs and vectors
#include "vtyacction.h"		// Global error log buffer functions
 
                                // Products
#include "vpbundle.h"
#include "vpfleg.h"
#include "vpcashfl.h"
#include "vpoption.h"
                                // VTree Tools
#include "vpipars.h"

#include "vtlprice.h"		// main price routine
 
extern  "C" {
#include "date_sup.h"
 
#include "drlstr.h"             // Strtok
#include "drltime.h"
#include "drlio.h"
#include "drlmem.h"
#include "drlproc.h"
#include "drlts.h"
#include "drlstr.h"
#include "drlgetop.h"           // getopt()
 
#include "dritkwrp.h"           // TDrWrapper
};



extern  int      debugLevel;

static	double	FreqCharToBasis(char freq);
static	String	ZcIdxToName(int idx);
static	String	ZcTypeCharToZcName(const char zcTypeChar);


#define IS_VALID_ARRAY(statementL) 	\
	((statementL) != NULL && (long)(statementL)[0] != 0)

extern "C" {


int YacctionInputs(
					// --- INSTRUMENT

					/* --- Receive Leg 	       */
    TDate  *resetDatesRecL,		/* (I) Reset dates             */
    TDate  *resetEffDatesRecL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesRecL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesRecL,		/* (I) Accrual end dates       */
    TDate  *payDatesRecL,		/* (I) Payment dates           */
    double *notionalsRecL,		/* (I) Notionals               */
    char   *formulasRecL,		/* (I) Payment formula	       */
    
					/* --- Pay Leg 		       */
    TDate  *resetDatesPayL,		/* (I) Reset dates             */
    TDate  *resetEffDatesPayL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesPayL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesPayL,		/* (I) Accrual end dates       */
    TDate  *payDatesPayL,		/* (I) Payment dates           */
    double *notionalsPayL,		/* (I) Notionals               */
    char   *formulasPayL,		/* (I) Payment formula	       */

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
    char   *discCashCurveL,		/* (I) Discount curve		   */

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

    char   *routine);

//--------------------------------------------------------------
//


GTO_EXPORT(int)
DVTreeYacctionL(
					//
					// --- INSTRUMENT
					//

					/* --- Receive Leg 	       */
    TDate  *resetDatesRecL,		/* (I) Reset dates             */
    TDate  *resetEffDatesRecL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesRecL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesRecL,		/* (I) Accrual end dates       */
    TDate  *payDatesRecL,		/* (I) Payment dates           */
    double *notionalsRecL,		/* (I) Notionals               */
    char   *formulasRecL,		/* (I) Payment formula	       */
    
					/* --- Pay Leg 		       */
    TDate  *resetDatesPayL,		/* (I) Reset dates             */
    TDate  *resetEffDatesPayL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesPayL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesPayL,		/* (I) Accrual end dates       */
    TDate  *payDatesPayL,		/* (I) Payment dates           */
    double *notionalsPayL,		/* (I) Notionals               */
    char   *formulasPayL,		/* (I) Payment formula	       */

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
    char   *discCashCurveL,		/* (I) Discount curves		   */

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

    double *outputsL) 			/* (O) outputs
					   [1] Option Price
					   [2] Swap Price
					   [3] Fugit
					   [4] Option exercise probability */
{
static  char		routine[]="DVTreeYacctionL";
	int		status=FAILURE;

	int		i;

			// Basis tree model and market data

	KZCurve		*zcD = NULL,
			*zc1 = NULL,
			*zc2 = NULL,
			*zcB = NULL;

	TDate		today;
	KMarketCurves	market;
	KVolDiag	*irVolDiag = NULL;
	KMrParam	irMrParam;
	KSmileParam	irSmileParam;
	KVolDiag	*bsVolDiag = NULL;
	KMrParam	bsMrParam;
	KSmileParam	bsSmileParam;
	double		corrIrBs = 0e0;

	KResetBank	vpResetBank;
	

			/* --- Two legs of swap */
	SharedPointer<KVPFloatLeg>	vpLegRec;
	SharedPointer<KVPFloatLeg>	vpLegPay;

			/* --- Cash flows */
	SharedPointer<KVPCashFlows>	vpCash;

	SharedPointer<KVPWBundle>	swap;

	SharedPointer<KRate>	floatRateRec,
				floatRatePay;

			/* --- Option  */
	SharedPointer<KVPOption>	optSwap;
	SharedPointer<KVPInstr>		vpRoot;

	KVPOption::Type	optType;

	KRate		floatRateRecNS,
			floatRatePayNS;

			/* --- Basis volatility curve	    */
	bool		isBasis = false;

	double		discZeroShift;
 
	KDayCc		bsDCC, liborDCC;

	KMap(String, double)	results;

    try {

	// Enable logging, redirect error messages to buffer.
	//
	ASSERT_OR_THROW(ARGSIZE(debugLevelL) >= 1)
	debugLevel = debugLevelL[1];
	if (debugLevel > 0)
		DppLogMsgSetFile("run.log");
        DppErrMsgBufferInit();



	//--------------------------------------------------------------
	//
	// Log call arguments
	//
	//--------------------------------------------------------------

	if (debugLevel > 0) {
            if (DrlLilVectLoggingFp(DppLogFp(), "VT_YACCTION",
		/* --- Receive Leg	       */
                DRL_TDATE_L, (void*) resetDatesRecL,      "RESET_DATES_1",
		DRL_TDATE_L, (void*) resetEffDatesRecL,	  "RESET_EFF_DATES_1",
		DRL_TDATE_L, (void*) accStartDatesRecL,	  "ACC_START_DATES_1",
		DRL_TDATE_L, (void*) accEndDatesRecL,     "ACC_END_DATES_1",
		DRL_TDATE_L, (void*) payDatesRecL,	  "PAY_DATES_1",
		DRL_FLOAT_L, (void*) notionalsRecL,	  "NOTIONALS_1",
		DRL_CHAR_BLOCK_L,(void*) formulasRecL, 	  "FORMULAS_1",
    
		/* --- Pay Leg 		       */
                DRL_TDATE_L, (void*) resetDatesPayL,      "RESET_DATES_2",
		DRL_TDATE_L, (void*) resetEffDatesPayL,	  "RESET_EFF_DATES_2",
		DRL_TDATE_L, (void*) accStartDatesPayL,	  "ACC_START_DATES_2",
		DRL_TDATE_L, (void*) accEndDatesPayL,     "ACC_END_DATES_2",
		DRL_TDATE_L, (void*) payDatesPayL,	  "PAY_DATES_2",
		DRL_FLOAT_L, (void*) notionalsPayL,	  "NOTIONALS_2",
		DRL_CHAR_BLOCK_L,(void*) formulasPayL, 	  "FORMULAS_2",

		/* --- Receive & pay legs      */
		DRL_CHAR_BLOCK_L,(void*) payDayCountStrsL, "PAY_DAY_COUNTS",
		DRL_CHAR_BLOCK_L,(void*) stubConvsL,	   "STUB_CONVS",

		/* --- Index of receive, pay       */
		DRL_FLOAT_L, (void*) indexWeightsL,	  "INDEX_WEIGHTS",
		DRL_FLOAT_L, (void*) indexRateSpreadsL,   "INDEX_SPREADS",
		DRL_CHAR_BLOCK_L,(void*) indexFreqsL,	  "INDEX_FREQS",
		DRL_CHAR_BLOCK_L,(void*) indexMaturitiesL,"INDEX_MATURITIES",
		DRL_CHAR_BLOCK_L,(void*) indexDayCountStrsL, "INDEX_DCC",	
		DRL_CHAR_BLOCK_L,(void*) indexCurvesL,	  "INDEX_CURVES",
		DRL_CHAR_BLOCK_L,(void*) discIndexCurvesL,"DISC_INDEX_CURVES",

    		DRL_TDATE_L, (void*)     cashDatesL, 	  "CASH_DATES",
		DRL_FLOAT_L, (void*)  	 cashAmountsL, 	  "CASH_AMOUNTS",
		DRL_CHAR_BLOCK_L,(void*) discCashCurveL,  "DISC_CASH_CURVE",

		/* --- Option 	       */
		DRL_CHAR_BLOCK_L,(void*) optTypeL,	  "OPT_TYPES",
		DRL_TDATE_L, (void*)  	 optNotifDatesL,  "OPT_NOTIF_DATES",
		DRL_TDATE_L, (void*)  	 optSettleDatesL, "OPT_SETTLE_DATES",
    		DRL_TDATE_L, (void*)     optStrikePayDatesL, "OPT_STKPAY_DATES",
		DRL_FLOAT_L, (void*)  	 optStrikesL, 	  "OPT_STRIKES",
		DRL_FLOAT_L, (void*)     optNotionalsL,   "OPT_NOTIONALS",
		DRL_CHAR_BLOCK_L,(void*) discOptCurveL,   "DISC_OPTION_CURVE",

		/*
		 * ENVIRONMENT + MODEL
		 */

		/* --- Zero curvs 	       */
		DRL_TDATE_L, (void*) 	 baseDatesL,	  "ZC_BASE_DATES",
		DRL_CHAR_BLOCK_L,(void*) zeroFrequenciesL,"ZC_FREQS",	
		DRL_CHAR_BLOCK_L,(void*) zeroDCCsL,	  "ZC_DCCS",
                                
		DRL_TDATE_L, (void*)     zeroDiffDatesL,  "ZC_DATES_D",
		DRL_FLOAT_L, (void*)     zeroDiffRatesL,  "ZC_RATES_D",	
                                                                  
		DRL_TDATE_L, (void*)     zero1DatesL, 	  "ZC_DATES_1",
		DRL_FLOAT_L, (void*)     zero1RatesL,	  "ZC_RATES_1",
                                                                  
		DRL_TDATE_L, (void*)     zero2DatesL,	  "ZC_DATES_2",
		DRL_FLOAT_L, (void*)     zero2RatesL,	  "ZC_RATES_2",
                                                                    
		DRL_TDATE_L, (void*)     zeroBDatesL,	  "ZC_DATES_B",
		DRL_FLOAT_L, (void*)     zeroBRatesL,	  "ZC_RATES_B",

		DRL_FLOAT_L, (void*)     zeroBInfoL,	  "ZC_BASIS_INFO",

		/* --- Volatility curvs        */
		DRL_TDATE_L, (void*)     irVolDatesL,	  "IR_VOL_DATES",
		DRL_TDATE_L, (void*)     irVolMatsL,	  "IR_VOL_MATS",
		DRL_LONG_L,  (void*) 	 irVolFreqL,	  "IR_VOL_FREQ",
		DRL_FLOAT_L, (void*)     irVolsL,	  "IR_VOL_RATES",

		DRL_TDATE_L, (void*)     basisVolDatesL,  "BS_VOL_DATES",	
		DRL_FLOAT_L, (void*)     basisVolsL,	  "BS_VOL_RATES",
		DRL_CHAR_BLOCK_L,(void*) volTypeL,	  "BS_VOL_TYPE",

		/* --- Model Parameters		*/
		DRL_FLOAT_L, (void*)     irMrParamL,	  "IR_MR_PARAMS",
		DRL_FLOAT_L, (void*)     irSmileParamL,	  "IR_SMILE_PARAMS",
		DRL_FLOAT_L, (void*)     bsMrParamL,	  "BS_MR_PARAMS",
		DRL_FLOAT_L, (void*)     bsSmileParamL,	  "BS_SMILE_PARAMS",
		DRL_FLOAT_L, (void*)     corrIrBsL,	  "CORR_IR_BS",

		/* --- Reset bank		*/
		DRL_CHAR_BLOCK_L,(void*) resetBankIndsL,  "RESET_BANK_INDEX",
		DRL_TDATE_L, (void*)  	 resetBankDatesL, "RESET_BANK_DATES",
		DRL_FLOAT_L, (void*)  	 resetBankRatesL, "RESET_BANK_RATES",

		/* --- Debugging level 		*/
		DRL_LONG_L,  (void*)     debugLevelL,	  "INT_SCALARS", 	
                0) != SUCCESS)
		throw KFailure("%s: failed to log the inputs.\n", routine);
        }


	//--------------------------------------------------------------
	//
	// Check and log arguments
	//
	//--------------------------------------------------------------

	ASSERT_OR_THROW( ARGSIZE(resetDatesRecL) == ARGSIZE(resetEffDatesRecL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesRecL) == ARGSIZE(accStartDatesRecL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesRecL) == ARGSIZE(accEndDatesRecL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesRecL) == ARGSIZE(payDatesRecL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesRecL) == ARGSIZE(notionalsRecL));
	// We do this because of the 128 size limitation for char blocks
	//
	if (ARGSIZE(resetDatesRecL) < 128) {
	    ASSERT_OR_THROW(ARGSIZE(resetDatesRecL) == ARGSIZE(formulasRecL));
	} else {
	    ASSERT_OR_THROW(ARGSIZE(formulasRecL) == 127);
	}
	
	ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(resetEffDatesPayL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(accStartDatesPayL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(accEndDatesPayL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(payDatesPayL));
	ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(notionalsPayL));
	// We do this because of the 128 size limitation for char blocks
	//
	if (ARGSIZE(resetDatesPayL) < 128) {
	    ASSERT_OR_THROW( ARGSIZE(resetDatesPayL) == ARGSIZE(formulasPayL));
	} else {
	    ASSERT_OR_THROW(ARGSIZE(formulasPayL) == 127);
	}

	ASSERT_OR_THROW( ARGSIZE(payDayCountStrsL)	== 2);
	ASSERT_OR_THROW( ARGSIZE(stubConvsL)		== 2);

	ASSERT_OR_THROW( ARGSIZE(indexWeightsL)		== 2);
	ASSERT_OR_THROW( ARGSIZE(indexRateSpreadsL)	== 2);
	ASSERT_OR_THROW( ARGSIZE(indexFreqsL)		== 2);
	ASSERT_OR_THROW( ARGSIZE(indexMaturitiesL)	== 2);
	ASSERT_OR_THROW( ARGSIZE(indexDayCountStrsL)	== 2);
	ASSERT_OR_THROW( ARGSIZE(indexCurvesL)		== 2);
	ASSERT_OR_THROW( ARGSIZE(discIndexCurvesL)	== 2);

	ASSERT_OR_THROW( ARGSIZE(optNotifDatesL)==ARGSIZE(optSettleDatesL));
	ASSERT_OR_THROW( ARGSIZE(optNotifDatesL)==ARGSIZE(optStrikePayDatesL));
	ASSERT_OR_THROW( ARGSIZE(optNotifDatesL)==ARGSIZE(optStrikesL));
	ASSERT_OR_THROW( ARGSIZE(optNotifDatesL)==ARGSIZE(optNotionalsL));

	ASSERT_OR_THROW( ARGSIZE(baseDatesL)   	   == 5);
	ASSERT_OR_THROW( ARGSIZE(zeroFrequenciesL) == 4);
	ASSERT_OR_THROW( ARGSIZE(zeroDCCsL)        == 4);

	ASSERT_OR_THROW( ARGSIZE(irVolDatesL)      == ARGSIZE(irVolMatsL));
	ASSERT_OR_THROW( ARGSIZE(irVolDatesL)      == ARGSIZE(irVolsL));
	ASSERT_OR_THROW( ARGSIZE(irVolDatesL)      == ARGSIZE(irVolsL));

	ASSERT_OR_THROW( ARGSIZE(irMrParamL)	> 4);
	ASSERT_OR_THROW( ARGSIZE(irSmileParamL)	== 4 );

	ASSERT_OR_THROW( ARGSIZE(outputsL)	>= 4);


	// Input logging
	//
	if (debugLevel > 0){
		YacctionInputs(
					/* --- Receive Leg	       */
    				resetDatesRecL,	
				resetEffDatesRecL,
				accStartDatesRecL,
				accEndDatesRecL,
				payDatesRecL,
				notionalsRecL,
				formulasRecL,
    
					/* --- Pay Leg 		       */
				resetDatesPayL,	
				resetEffDatesPayL,
				accStartDatesPayL,	
				accEndDatesPayL,
				payDatesPayL,
				notionalsPayL,	
				formulasPayL,

					/* --- Receive & pay legs      */
				payDayCountStrsL,
				stubConvsL,

					/* --- Index of receive, pay   */
				indexWeightsL,	
				indexRateSpreadsL,
				indexFreqsL,	
				indexMaturitiesL,
				indexDayCountStrsL,	
				indexCurvesL,	
				discIndexCurvesL,

					/* --- Cash flows */
    				cashDatesL,
    				cashAmountsL,	
    				discCashCurveL,

					/* --- option 	       */
    				optTypeL,
    				optNotifDatesL,
    				optSettleDatesL,
    				optStrikePayDatesL,
    				optStrikesL,
    				optNotionalsL,
    				discOptCurveL,

					//
					// ENVIRONMENT + MODEL
					//

					/* --- Zero curvs 	       */
				baseDatesL,	
				zeroFrequenciesL,	
				zeroDCCsL,
                                
				zeroDiffDatesL,	
				zeroDiffRatesL,	
                                                                  
				zero1DatesL,
				zero1RatesL,
                                                                  
				zero2DatesL,
				zero2RatesL,	
                                                                    
				zeroBDatesL,
				zeroBRatesL,

				zeroBInfoL,

					/* --- Volatility curvs        */
				irVolDatesL,	
				irVolMatsL,
				irVolFreqL,	
				irVolsL,

				basisVolDatesL,	
				basisVolsL,
				volTypeL,

					/* --- Model Parameters		*/
				irMrParamL,
				irSmileParamL,	
				bsMrParamL,
				bsSmileParamL,	
				corrIrBsL,

					/* --- Reset bank		*/
				resetBankIndsL,
				resetBankDatesL,
				resetBankRatesL,

				debugLevelL,	

				routine);
	};


	//--------------------------------------------------------------
	//
	// Construct Market and Model information 
	//
	//--------------------------------------------------------------

	// 
	// Today date
	//
	today = baseDatesL[1];
	market.mToday = today;

	// 
	// Construct zero curves and insert them in market environment
	//
	if(IS_VALID_ARRAY(zeroDiffDatesL)) {
		ASSERT_OR_THROW( ARGSIZE(zeroDiffDatesL) 
						== ARGSIZE(zeroDiffRatesL));
		zcD = new KZCurve(
			baseDatesL[2],
			DppVectorFromArray((int)zeroDiffDatesL[0], 
					   zeroDiffDatesL+1),
			DppVectorFromArray((int)zeroDiffRatesL[0], 
					   zeroDiffRatesL+1),
			FreqCharToBasis(zeroFrequenciesL[WRAP_STR_IDX(1)]),
			KDayCc(&zeroDCCsL[WRAP_STR_IDX(1)]));

		market.Insert(
			KZCurve((TCurve*)*zcD, today),
			KV_DIFF,
			baseDatesL[2],
			ZcTypeCharToZcName('D'));
	}
	else
		throw KFailure("%s: invalid diffuse zero curve.\n",
				routine);

	if(IS_VALID_ARRAY(zero1DatesL)) {
		ASSERT_OR_THROW( ARGSIZE(zero1DatesL) == ARGSIZE(zero1RatesL));
	    	zc1 = new KZCurve(
			baseDatesL[3],
			DppVectorFromArray((int)zero1DatesL[0], zero1DatesL+1),
			DppVectorFromArray((int)zero1RatesL[0], zero1RatesL+1),
			FreqCharToBasis(zeroFrequenciesL[WRAP_STR_IDX(2)]),
			KDayCc(&zeroDCCsL[WRAP_STR_IDX(2)]));

	    	market.Insert(
			KZCurve((TCurve*)*zc1, today),
			KV_IDX1,
			baseDatesL[3],
			ZcTypeCharToZcName('1'));
		market.InsertValueDate(KV_IDX1, baseDatesL[3]);
	}

	if(IS_VALID_ARRAY(zero2DatesL)) {
		ASSERT_OR_THROW( ARGSIZE(zero2DatesL) == ARGSIZE(zero2RatesL));
	     	zc2 = new KZCurve(
			baseDatesL[4],
			DppVectorFromArray((int)zero2DatesL[0], zero2DatesL+1),
			DppVectorFromArray((int)zero2RatesL[0], zero2RatesL+1),
			FreqCharToBasis(zeroFrequenciesL[WRAP_STR_IDX(3)]),
			KDayCc(&zeroDCCsL[WRAP_STR_IDX(3)]));

	    	market.Insert(
			KZCurve((TCurve*)*zc2, today),
			KV_IDX2,
			baseDatesL[4],
			ZcTypeCharToZcName('2'));
		market.InsertValueDate(KV_IDX2, baseDatesL[4]);
	}

	if (IS_VALID_ARRAY(zeroBDatesL) && IS_VALID_ARRAY(zeroBInfoL)) {
		isBasis = true;

		ASSERT_OR_THROW( ARGSIZE(zeroBDatesL) == ARGSIZE(zeroBRatesL));
		ASSERT_OR_THROW( ARGSIZE(zeroBInfoL) >= 4);

		KSpd	bsType;

		zcB = new KZCurve(
			baseDatesL[5],
			DppVectorFromArray((int)zeroBDatesL[0], zeroBDatesL+1),
			DppVectorFromArray((int)zeroBRatesL[0], zeroBRatesL+1),
			FreqCharToBasis(zeroFrequenciesL[WRAP_STR_IDX(4)]),
			KDayCc(&zeroDCCsL[WRAP_STR_IDX(4)]));


		switch((int) (zeroBInfoL[3])) {
		case 0:
			bsType = SUB_SPREAD;
			break;
		case 1:
			bsType = PER_SPREAD;
			break;
		case 2:
			bsType = ADD_SPREAD;
			break;
		default:
			throw KFailure("%s: invalid basis type (%d).\n",
				       routine, ((int) (zeroBInfoL[3])));
		}

		// Specify the basis and libor DCCs
		if (ARGSIZE(zeroBInfoL) >= 6)
		{
			bsDCC    = KDayCc((int)zeroBInfoL[5]);
			liborDCC = KDayCc((int)zeroBInfoL[6]);
		}
		else {
			bsDCC    = KDayCc(GTO_ACT_365);
			liborDCC = KDayCc(GTO_ACT_360);
                }


		// Add to market env
		market.InsertBasis(
			KZCurve((TCurve*)*zcB, today),
			KV_BASIS,
			baseDatesL[5],
			ZcTypeCharToZcName('B'),
			bsType,
			ZcIdxToName((int) zeroBInfoL[1]),
			ZcIdxToName((int) zeroBInfoL[2]),
			zeroBInfoL[4],
			bsDCC,
			liborDCC);


		// If index curve is spread or par spread,
		// insert additional curve name and type to market table
		for (i=1; i<=3; ++i)
		{
		    if (toupper(indexCurvesL[WRAP_STR_IDX(i)]) == 'S')
	    		market.InsertZc(
				KZCurve((TCurve*)*zcB, today),
				KV_SPREAD,
				baseDatesL[5],
				ZcTypeCharToZcName('S'));
		    else if(toupper(indexCurvesL[WRAP_STR_IDX(i)]) == 'P')
	    		market.InsertZc(
				KZCurve((TCurve*)*zcB, today),
				KV_PAR_SPREAD,
				baseDatesL[5],
				ZcTypeCharToZcName('P'));
		}

	}



	// 
	// Volatilities
	//

	/* IR volatility */
	if ((int)irVolDatesL[0] == 1) {
		// Nil calibration $$$
		throw KFailure("%s: nil calibration not yet supported. Sorry.\n",
			routine);
	} else {
		irVolDiag = new KVolDiag(
			DppVectorFromArray((int)irVolDatesL[0], irVolDatesL+1),
			DppVectorFromArray((int)irVolMatsL[0], irVolMatsL+1),
			DppVectorFromArray((int)irVolFreqL[0], irVolFreqL+1),
			DppVectorFromArray((int)irVolsL[0], irVolsL+1));
	}



	/* Basis volatility diag:
	 * interp same at same dates as irVolDiag
	 */
	bsVolDiag = new KVolDiag(*irVolDiag);
	*bsVolDiag = 0e0;	// clear values
	if (isBasis) 
	{
	    bsVolDiag->BasisVolDiag(
		&volTypeL[1],
		DppVectorFromArray((int)basisVolDatesL[0], basisVolDatesL+1),	
		DppVectorFromArray((int)basisVolsL[0], basisVolsL+1),	
		4L, /* VOL FREQ SHOULD BE AN INPUT   */
		*irVolDiag);
	}

	//
	// Model parameters
	//
	irMrParam.ReadLil(irMrParamL);
	irSmileParam.ReadLil(irSmileParamL);

	if (isBasis) 
	{
		ASSERT_OR_THROW( ARGSIZE(bsMrParamL)	> 0);
		ASSERT_OR_THROW( ARGSIZE(bsSmileParamL)	> 0);
		ASSERT_OR_THROW( ARGSIZE(corrIrBsL)	> 0);

		bsMrParam.ReadLil(bsMrParamL);
		bsSmileParam.ReadLil(bsSmileParamL);

		corrIrBs = corrIrBsL[1];
	}



	//--------------------------------------------------------------
	//
	// Construct Product 
	//
	//--------------------------------------------------------------

	/* 
	 * Construct the receive leg
	 */


	/* Construct the KRate  */
	if (IS_ALMOST_ZERO(indexWeightsL[1]))
	{
		floatRateRec = Raw2SharedPointer(
				new KRate(indexRateSpreadsL[1]));	
	}	
	else
	{
		floatRateRec = Raw2SharedPointer(
			new KRate(
			ZcTypeCharToZcName(indexCurvesL[WRAP_STR_IDX(1)]),
			KDateInterval(&indexMaturitiesL[WRAP_STR_IDX(1)]),
			KDateInterval(&indexFreqsL[WRAP_STR_IDX(1)]),
			KDayCc(&indexDayCountStrsL[WRAP_STR_IDX(1)]),
			KDateInterval(0e0),   /* use resetEffDates */
			indexRateSpreadsL[1],
			indexWeightsL[1]));
	}


	vpLegRec = Raw2SharedPointer(
		new KVPFloatLeg(
		"Swap Rec Leg",
		DppVectorFromArray((int)resetDatesRecL[0], resetDatesRecL+1),
		DppVectorFromArray((int)resetEffDatesRecL[0], 
				    resetEffDatesRecL+1),
		DppVectorFromArray((int)accStartDatesRecL[0], 
				    accStartDatesRecL+1),
		DppVectorFromArray((int)accEndDatesRecL[0], accEndDatesRecL+1),
		DppVectorFromArray((int)payDatesRecL[0], payDatesRecL+1),
		DppVectorFromArray((int)notionalsRecL[0], notionalsRecL+1),
		DppStrVectorFromCharBlock((int)resetDatesRecL[0], formulasRecL),
		KDayCc(&payDayCountStrsL[WRAP_STR_IDX(1)]),
		KStubConv(&stubConvsL[WRAP_STR_IDX(1)]),
		floatRateRec,
		ZcTypeCharToZcName(discIndexCurvesL[WRAP_STR_IDX(1)]).c_str()));

	/* 
	 * Construct the pay leg 
	 */

	/* Construct the KRate  */
	if (IS_ALMOST_ZERO(indexWeightsL[2]))
	{
		floatRatePay = Raw2SharedPointer(
			new KRate(indexRateSpreadsL[2]));	
	}	
	else
	{
		floatRatePay = Raw2SharedPointer(
			new KRate(
			ZcTypeCharToZcName(indexCurvesL[WRAP_STR_IDX(2)]),
			KDateInterval(&indexMaturitiesL[WRAP_STR_IDX(2)]),
			KDateInterval(&indexFreqsL[WRAP_STR_IDX(2)]),
			KDayCc(&indexDayCountStrsL[WRAP_STR_IDX(2)]),
			KDateInterval(0e0),   /* use resetEffDates */
			indexRateSpreadsL[2],
			indexWeightsL[2]));
	}

	vpLegPay = Raw2SharedPointer(
		new KVPFloatLeg(
		"Swap Pay Leg",
		DppVectorFromArray((int)resetDatesPayL[0], resetDatesPayL+1),
		DppVectorFromArray((int)resetEffDatesPayL[0], 
				   resetEffDatesPayL+1),
		DppVectorFromArray((int)accStartDatesPayL[0], 
			  	   accStartDatesPayL+1),
		DppVectorFromArray((int)accEndDatesPayL[0], accEndDatesPayL+1),
		DppVectorFromArray((int)payDatesPayL[0], payDatesPayL+1),
		DppVectorFromArray((int)notionalsPayL[0], notionalsPayL+1),
		DppStrVectorFromCharBlock((int)resetDatesPayL[0], formulasPayL),
		KDayCc(&payDayCountStrsL[WRAP_STR_IDX(2)]),
		KStubConv(&stubConvsL[WRAP_STR_IDX(2)]),
		floatRatePay,
		ZcTypeCharToZcName(discIndexCurvesL[WRAP_STR_IDX(2)]).c_str()));


	// 
	// Construct the swap from two legs 
	//
	swap = Raw2SharedPointer(new KVPWBundle("Swap"));

	// Cast to KVPAtom
	SharedPointer<KVPInstr> vpRec;
	SharedPointer<KVPInstr> vpPay;
 
	SharedPointerConvertTo(vpLegRec, vpRec);
	SharedPointerConvertTo(vpLegPay, vpPay);

	swap->AddDepWeight(vpRec, 1e0);
	swap->AddDepWeight(vpPay, -1e0);


	// Add cash flows
	//
	
	if(IS_VALID_ARRAY(cashDatesL))	// ignored if NULL
	{
	    ASSERT_OR_THROW( ARGSIZE(cashDatesL)  == ARGSIZE(cashAmountsL));

	    vpCash = Raw2SharedPointer(
		new KVPCashFlows(
		"Cash Flows",
		DppVectorFromArray((int)cashDatesL[0], cashDatesL+1),
		DppVectorFromArray((int)cashAmountsL[0], cashAmountsL+1),
		ZcTypeCharToZcName(discCashCurveL[WRAP_STR_IDX(1)]).c_str()));

	    SharedPointer<KVPInstr> vpC;
	    SharedPointerConvertTo(vpCash, vpC);

	    swap->AddDepWeight(vpC, 1e0);
	}


	//
	// Build Option instrument
	//

	/* Option type */
	switch (toupper(optTypeL[WRAP_STR_IDX(1)])) {
	case 'C':
		optType = KVPOption::CALL;
		break;
	case 'P':
		optType = KVPOption::PUT;
		break;
	default:
		throw KFailure("%s: invalid option type (%c).\n",
				routine,
				optTypeL[WRAP_STR_IDX(1)]);
	}
 
	optSwap = Raw2SharedPointer(
		new KVPOption(
		"Yacction",
		optType,
		FALSE,
		DppVectorFromArray((int)optNotifDatesL[0], optNotifDatesL+1),
		DppVectorFromArray((int)optSettleDatesL[0], optSettleDatesL+1),
		DppVectorFromArray((int)optStrikePayDatesL[0], 
				   optStrikePayDatesL+1),
		DppVectorFromArray((int)optStrikesL[0], optStrikesL+1),
		DppVectorFromArray((int)optNotionalsL[0], optNotionalsL+1),
		KDateInterval(0e0),	// notif days, used by American option
		ZcTypeCharToZcName(discOptCurveL[WRAP_STR_IDX(1)]).c_str()));



	/* 
	 * Add the dependency of underlying swap 
	 */
	SharedPointer<KVPAtom> vpSwap;
	SharedPointerConvertTo(swap, vpSwap);

	optSwap->AddDep(vpSwap);


	//--------------------------------------------------------------
	//
	// Reset bank
	//
	//--------------------------------------------------------------

	if(IS_VALID_ARRAY(resetBankIndsL))	// ignored if NULL
	{
		ASSERT_OR_THROW(ARGSIZE(resetBankIndsL) 
				== ARGSIZE(resetBankDatesL));
		ASSERT_OR_THROW(ARGSIZE(resetBankIndsL) 
				== ARGSIZE(resetBankRatesL));
 
		/* Insert into reset bank */
		for (i=1; i<=(int)resetBankDatesL[0]; i++)
		{
			switch(toupper(resetBankIndsL[WRAP_STR_IDX(i)])) {
			case 'R':
				floatRateRecNS = *floatRateRec;
				floatRateRecNS.SetSpread(0e0);

				vpResetBank.Insert(floatRateRecNS,
					   resetBankDatesL[i],
					   resetBankRatesL[i]);
				break;
			case 'P':
				floatRatePayNS = *floatRatePay;
				floatRatePayNS.SetSpread(0e0);

				vpResetBank.Insert(floatRatePayNS,
					   resetBankDatesL[i],
					   resetBankRatesL[i]);
				break;
			default:
				throw KFailure("%s: invalid (%d)th index rate "
				       "type (%c).  Choose among (R)eceive and "
				       "(P)ay.\n",
				       routine,
				       i,
				       resetBankIndsL[WRAP_STR_IDX(i)]);
			}
		}
	}



	/* 
	 * Zero shift from value date to today for discount curve.
	 */
	discZeroShift = market.ZeroShift(ZcTypeCharToZcName(
					discOptCurveL[WRAP_STR_IDX(1)]));


	//----------------------------------------------
	// Log preprocessed inputs
	//----------------------------------------------
	if (debugLevel > 0) {
		dppLog << "------------------------------------------" << endl;
		dppLog << "Product:" << endl;

		dppLog << "------------------------------------------" << endl;
                dppLog << "Root Dependency Graph: " << endl;
		optSwap->PutTree(dppLog);

		dppLog << "------------------------------------------" << endl;
                dppLog << "Details: " << endl;
		optSwap->PutRecursive(dppLog);

		dppLog << "------------------------------------------" << endl;
		dppLog << "Tree inputs:" << endl;
		dppLog << "market:\n"       << market << endl;
		dppLog << "irVolDiag:\n"    << *irVolDiag << endl;
		dppLog << "irMrParam:\n"    << irMrParam << endl;
		dppLog << "irSmileParam:\n" << irSmileParam << endl;
		dppLog << "bsVolDiag:\n"    << *bsVolDiag << endl;
		dppLog << "bsMrParam:\n"    << bsMrParam << endl;
		dppLog << "bsSmileParam:\n" << bsSmileParam << endl;
		dppLog << "corrIrBs:\n"     << corrIrBs << endl;
		dppLog << "resetBank:\n"    << vpResetBank << endl;

	}

	//----------------------------------------------
	// Call the main pricing routine
	//----------------------------------------------
	SharedPointerDownCastTo(optSwap, vpRoot);

	KVPRootPrice(
		results, 
 
		vpRoot,
		market, 
		*irVolDiag,
		irMrParam,
		irSmileParam,
		*bsVolDiag, 
		bsMrParam,
		bsSmileParam,
		corrIrBs,
		vpResetBank,
		debugLevel);


	/*
	 * Output results
	 */

	outputsL[1] = results["PV"] / discZeroShift;
	outputsL[2] = results["UNDER"]/discZeroShift;
	outputsL[3] = results["XTFUG"];
	outputsL[4] = results["XTPROB"];


	/* Free memory */

	delete zcD;
	if (zero1DatesL)
		delete zc1;
	if (zero2DatesL)
		delete zc2;
	if (zeroBDatesL)
		delete zcB;

	delete irVolDiag;
	delete bsVolDiag;

	status = SUCCESS;

   }
   catch (KFailure) {


	delete zcD;
	if (zero1DatesL)
		delete zc1;
	if (zero2DatesL)
		delete zc2;
	if (zeroBDatesL)
		delete zcB;


	delete irVolDiag;
	delete bsVolDiag;

	dppErr << format("%s: failed.\n", routine);

	return (status);
   }
   catch (...) {
	
	dppErr << format("%s: failed (uncaught).\n", routine);

	return (status);
   }

	return (status);
 
}

};	// extern "C"


//--------------------------------------------------------------
//
int YacctionInputs(
					// --- INSTRUMENT

					/* --- Receive Leg 	       */
    TDate  *resetDatesRecL,		/* (I) Reset dates             */
    TDate  *resetEffDatesRecL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesRecL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesRecL,		/* (I) Accrual end dates       */
    TDate  *payDatesRecL,		/* (I) Payment dates           */
    double *notionalsRecL,		/* (I) Notionals               */
    char   *formulasRecL,		/* (I) Payment formula	       */
    
					/* --- Pay Leg 		       */
    TDate  *resetDatesPayL,		/* (I) Reset dates             */
    TDate  *resetEffDatesPayL,		/* (I) Reset effective dates   */
    TDate  *accStartDatesPayL,		/* (I) Accrual start dates     */
    TDate  *accEndDatesPayL,		/* (I) Accrual end dates       */
    TDate  *payDatesPayL,		/* (I) Payment dates           */
    double *notionalsPayL,		/* (I) Notionals               */
    char   *formulasPayL,		/* (I) Payment formula	       */

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
    char   *discCashCurveL,		/* (I) Discount curve		   */

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
					 * [4] Delay shift (intvl) */


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

    char   *routine)
{
	int             status    = FAILURE;
	int		iRec, iPay, idx;
	int		numFactorIR, numFactorBS;

	bool		isBasis = false;

	dppLog << format("\n ===================== %s INPUTS: =====================\n", routine);

	dppLog << endl;

	dppLog << format("Today:  %s\n", GtoFormatDate(baseDatesL[1]));

	dppLog << endl;


	dppLog << format("\n ======================== Receive Leg: ========================\n");
	dppLog << format(" Reset       Eff Reset   Accr St     Accr End    Pay Date    Notional     Formula\n");
	for(iRec=0; iRec<=(long)resetDatesRecL[0]-1; iRec++)
	{
		dppLog << format("%10s  ", 
			  GtoFormatDate(resetDatesRecL[iRec+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(resetEffDatesRecL[iRec+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(accStartDatesRecL[iRec+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(accEndDatesRecL[iRec+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(payDatesRecL[iRec+1]));
		dppLog << format(" %lf\n", 
			  notionalsRecL[iRec+1]);
		dppLog << format(" %s\n", 
			  &formulasRecL[WRAP_STR_IDX(iRec+1)]);
	}

	dppLog << endl;

	dppLog << format("\n ======================== Pay Leg: ========================\n");
	dppLog << format(" Reset       Eff Reset   Accr St     Accr End    Pay Date    Notional     Formula\n");
	for(iPay=0; iPay<=(long)resetDatesPayL[0]-1; iPay++)
	{
		dppLog << format("%10s  ", 
			  GtoFormatDate(resetDatesPayL[iPay+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(resetEffDatesPayL[iPay+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(accStartDatesPayL[iPay+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(accEndDatesPayL[iPay+1]));
		dppLog << format("%10s  ", 
			  GtoFormatDate(payDatesPayL[iPay+1]));
		dppLog << format(" %lf\n", 
			  notionalsPayL[iPay+1]);
		dppLog << format(" %s\n", 
			  &formulasPayL[WRAP_STR_IDX(iPay+1)]);
	}

	dppLog << endl;

	dppLog << format("\n ======================== Payment Conventions: ========================\n");
	dppLog << format("                 Receive \t Pay \n");
	dppLog << format("Payment DCC:     %s \t %s\n",
		  &payDayCountStrsL[WRAP_STR_IDX(1)],
		  &payDayCountStrsL[WRAP_STR_IDX(2)]);
	dppLog << format("Stub Convention: %s \t %s\n",
		  &stubConvsL[WRAP_STR_IDX(1)],
		  &stubConvsL[WRAP_STR_IDX(2)]);
	
	dppLog << endl;

	dppLog << format("\n ======================== Rate Indices: ========================\n");
	dppLog << format("                Receive \t Pay \n");
	dppLog << format("Weight:         %lf \t %lf \n",
		  indexWeightsL[1], 
		  indexWeightsL[2]); 
	dppLog << format("Spread:         %lf \t %lf \n",
	      	  indexRateSpreadsL[1], 
		  indexRateSpreadsL[2]); 
	dppLog << format("Freq:           %s \t\t %s \n",
	      	  &indexFreqsL[WRAP_STR_IDX(1)], 
		  &indexFreqsL[WRAP_STR_IDX(2)]); 
	dppLog << format("Maturity:       %s \t\t %s \n",
	      	  &indexMaturitiesL[WRAP_STR_IDX(1)], 
	  	  &indexMaturitiesL[WRAP_STR_IDX(2)]); 
	dppLog << format("Day Count Conv: %s    \t %s \n",
	      	  &indexDayCountStrsL[WRAP_STR_IDX(1)], 
		  &indexDayCountStrsL[WRAP_STR_IDX(2)]); 
	dppLog << format("Index Curve:    %c \t\t %c \n",
	      	  indexCurvesL[WRAP_STR_IDX(1)], 
		  indexCurvesL[WRAP_STR_IDX(2)]); 
	dppLog << format("Disc Curve:     %c \t\t %c \n",
	      	  discIndexCurvesL[WRAP_STR_IDX(1)], 
		  discIndexCurvesL[WRAP_STR_IDX(2)]); 

	dppLog << endl;


	if(IS_VALID_ARRAY(cashDatesL))
	{
	    ASSERT_OR_THROW( ARGSIZE(cashDatesL)  == ARGSIZE(cashAmountsL));

	    dppLog << format("\n ======================== Cash Flow Schedule: ========================\n");
	    dppLog << format("Discount Curve: %s.\n", 
			 &discCashCurveL[WRAP_STR_IDX(1)]);
	    dppLog << endl;
	    dppLog << format("  DATES     CASH \n");
	    for(idx=0; idx<=(long)cashDatesL[0]-1; idx++)
		dppLog << format("%10s  %12.6f\n",
			   GtoFormatDate(cashDatesL[idx+1]),
			   cashAmountsL[idx+1]);
	    dppLog << endl;
	}


	dppLog << format("\n ======================== Option Schedule: ========================\n");
	dppLog << format("Option Type: %s.\n", &optTypeL[WRAP_STR_IDX(1)]);
	dppLog << endl;
	dppLog << format("Discount Curve: %s.\n", 
			 &discOptCurveL[WRAP_STR_IDX(1)]);
	dppLog << endl;
	dppLog << format(" NOTIF       SETTLE      STRIKE_PAY  STRIKE     NOTIONAL\n");
	for(idx=0; idx<=(long)optNotifDatesL[0]-1; idx++) 
		dppLog << format("%10s  %10s  %10s %10.6f   %12.6f\n",
                        GtoFormatDate(optNotifDatesL[idx+1]),
                        GtoFormatDate(optSettleDatesL[idx+1]),
                        GtoFormatDate(optStrikePayDatesL[idx+1]),
                        optStrikesL[idx+1],
                        optNotionalsL[idx+1]);

	dppLog << endl;

	dppLog << format("\n ======================== Zero Curves: ========================\n");
	dppLog << format("Diffuse Curve:\n\n");
	dppLog << format("Value Date: %s\n", GtoFormatDate(baseDatesL[2]));
	dppLog << format("Frequency: %s\n", &zeroFrequenciesL[WRAP_STR_IDX(1)]);
	dppLog << format("Swap DCC:   %s\n", &zeroDCCsL[WRAP_STR_IDX(1)]);

	dppLog << endl;
	dppLog << format("  Dates        Rates\n");
	for(idx=0; idx<= (long)zeroDiffDatesL[0]-1; idx++)
	{
		dppLog << format("%10s   %8.4f%%\n",
			  GtoFormatDate(zeroDiffDatesL[idx+1]),
			  zeroDiffRatesL[idx+1]*1e2);
	}	
                                                                  
	dppLog << endl;

	if(IS_VALID_ARRAY(zero1DatesL))
	{
		dppLog << format("Index Curve 1:\n\n");
		dppLog << format("Value Date: %s\n", GtoFormatDate(baseDatesL[3]));
		dppLog << format("Frequency: %s\n", &zeroFrequenciesL[WRAP_STR_IDX(2)]);
		dppLog << format("Swap DCC:   %s\n", &zeroDCCsL[WRAP_STR_IDX(2)]);

		dppLog << format("\n");
		dppLog << format("  Dates        Rates\n");
		for(idx=0; idx<= (long)zero1DatesL[0]-1; idx++)
		{
			dppLog << format("%10s   %8.4f%%\n",
			  	  GtoFormatDate(zero1DatesL[idx+1]),
			  	  zero1RatesL[idx+1]*1e2);
		}	
	}

	dppLog << endl;

	if(IS_VALID_ARRAY(zero2DatesL))
	{
		dppLog << format("Index Curve 2:\n\n");
		dppLog << format("Value Date: %s\n", GtoFormatDate(baseDatesL[4]));
		dppLog << format("Frequency: %s\n", &zeroFrequenciesL[WRAP_STR_IDX(3)]);
		dppLog << format("Swap DCC:   %s\n", &zeroDCCsL[WRAP_STR_IDX(3)]);

		dppLog << format("\n");
		dppLog << format("  Dates        Rates\n");
		for(idx=0; idx<= (long)zero2DatesL[0]-1; idx++)
		{
			dppLog << format("%10s   %8.4f%%\n",
			  	  GtoFormatDate(zero2DatesL[idx+1]),
			  	  zero2RatesL[idx+1]*1e2);
		}	
	}

	dppLog << endl;

	if(IS_VALID_ARRAY(zeroBDatesL) && IS_VALID_ARRAY(zeroBDatesL))
	{
		isBasis = true;

		dppLog << format("Basis Curve:\n\n");
		dppLog << format("Libor Index Curve: %d\n", (long)zeroBInfoL[1]);
		dppLog << format("Disc Index Curve:  %d\n", (long)zeroBInfoL[2]);
		if (ARGSIZE(zeroBInfoL) >= 6)
		{
		    dppLog << format("Basis DCC :  %s\n", 
				GtoFormatDayCountConv((int)zeroBInfoL[5]));
		    dppLog << format("Libor DCC :  %s\n", 
				GtoFormatDayCountConv((int)zeroBInfoL[6]));
		}

		dppLog << format("Basis Type:        %d\n", (long)zeroBInfoL[3]);
		dppLog << format("Delay Shift:       %lf\n", zeroBInfoL[4]);

		dppLog << format("\n");

		dppLog << format("Value Date:        %s\n", 
			   GtoFormatDate(baseDatesL[5]));
		dppLog << format("Frequency: %s\n", &zeroFrequenciesL[WRAP_STR_IDX(4)]);
		dppLog << format("Swap DCC:          %s\n", 
			   &zeroDCCsL[WRAP_STR_IDX(4)]);

		dppLog << format("\n");
		dppLog << format("  Dates        Rates\n");
		for(idx=0; idx<= (long)zeroBDatesL[0]-1; idx++)
		{
			dppLog << format("%10s   %8.4f%%\n",
			  	  GtoFormatDate(zeroBDatesL[idx+1]),
			  	  zeroBRatesL[idx+1]*1e2);
		}	
	}

	dppLog << endl;

	dppLog << format("\n ======================== Vol Curves: ========================\n");
	
	dppLog << endl;
	dppLog << format("IR Vol Curve:\n");
	
	dppLog << format(" Vol Dates    Maturity    Frequency   Volatility\n");
	for(idx=0; idx<=(long)irVolDatesL[0]-1; idx++)
	{
		dppLog << format("%10s  %10s  %6d     %10.4f%%\n",
                        GtoFormatDate(irVolDatesL[idx+1]),
                        GtoFormatDate(irVolMatsL[idx+1]),
                        irVolFreqL[idx+1],
                        irVolsL[idx+1]*1e2);
	}
	
	dppLog << format("\n");

	if(isBasis)
	{
		dppLog << format("Basis Spread Vol Curve:\n");
		dppLog << format("Vol Type : %s\n", 
				 &volTypeL[WRAP_STR_IDX(1)]);

		dppLog << endl;

		dppLog << format(" Vol Dates    Volatility\n");
		for(idx=0; idx<=(long)basisVolDatesL[0]-1; idx++)
		{
			dppLog << format("%10s %10.4f%%\n",
                        	  GtoFormatDate(basisVolDatesL[idx+1]),
                        	  basisVolsL[idx+1]*1e2);
		}
	}

	dppLog << endl;

	dppLog << format("\n ======================== Model Parameters: ========================\n");

	numFactorIR = (int)irMrParamL[5];
	
	int esz = 5 + 2*numFactorIR + numFactorIR*(numFactorIR-1)/2;
	if (ARGSIZE(irMrParamL) != esz) {
	    throw KFailure("%s: expects array length %d+5 with %d factors "
		"(got %d).\n", 
		routine, esz-5, 
		numFactorIR, ARGSIZE(irMrParamL));
	}

	dppLog << format("\n");
	dppLog << format("IR Model Parameters:\n");
	dppLog << format("PPY:               %d\n", (int)irMrParamL[1]);
	dppLog << format("Num Stdev Cut:     %d\n", (int)irMrParamL[2]);
	dppLog << format("Tree Smoothing:    %d\n", (int)irMrParamL[3]);
	dppLog << format("Backbone Q:        %f\n", irMrParamL[4]);
	dppLog << format("Factor Number:     %d\n", numFactorIR);
	dppLog << format("Mean Reversion:    ");
	for(idx=0; idx<=numFactorIR-1; idx++)
		dppLog << format("%lf\t", irMrParamL[6+idx]);
	dppLog << endl;

	dppLog << format("Factor Weighting:  ");
	for(idx=0; idx<=numFactorIR-1; idx++)
		dppLog << format("%lf\t", irMrParamL[6+numFactorIR+idx]);
	dppLog << endl;
	
	if(numFactorIR > 1)
	{
		dppLog << format("Factor Correlation: ");
		for(idx=0; idx<=numFactorIR*(numFactorIR-1)/2-1; idx++)
			dppLog << format("%lf\t", irMrParamL[6+2*numFactorIR+idx]);
		dppLog << endl;
	}

	dppLog << format("Smile Q1:          %f\n", irSmileParamL[1]);
	dppLog << format("Smile Q2:          %f\n", irSmileParamL[2]);
	dppLog << format("Forward Shift:     %f\n", irSmileParamL[3]);
	dppLog << format("Number of CET iteration: %d\n", (int)irSmileParamL[4]);
	dppLog << endl;


	if(isBasis)
	{
		numFactorBS = (int)bsMrParamL[5];

		esz = 5 + 2*numFactorBS + numFactorBS*(numFactorBS-1)/2;
		if (ARGSIZE(bsMrParamL) != esz) {
	    		throw KFailure("%s: expects at array length %d+5 with "
				       "%d factors (got %d).\n", 
				routine, esz-5, 
				numFactorBS, ARGSIZE(bsMrParamL));
		}
	
		dppLog << format("Basis Model Parameters:\n");
		dppLog << format("Factor Number:     %d\n", numFactorBS);
		dppLog << format("Mean Reversion:    ");
		for(idx=0; idx<=numFactorBS-1; idx++)
			dppLog << format("%lf\t", bsMrParamL[6+idx]);
		dppLog << endl;

		dppLog << format("Factor Weighting:  ");
		for(idx=0; idx<=numFactorBS-1; idx++)
		       dppLog << format("%lf\t", bsMrParamL[6+numFactorBS+idx]);
		dppLog << endl;
	
		if(numFactorBS > 1)
		{
			dppLog << format("Factor Correlation: ");
			for(idx=0; idx<=numFactorBS*(numFactorBS-1)/2-1; idx++)
				dppLog << format("%lf\t", 
					      bsMrParamL[6+2*numFactorBS+idx]);
			dppLog << endl;
		}

		dppLog << format("Smile Q1:          %f\n", bsSmileParamL[1]);
		dppLog << format("Smile Q2:          %f\n", bsSmileParamL[2]);
		dppLog << format("Forward Shift:     %f\n", bsSmileParamL[3]);

		dppLog << endl;

		if(numFactorIR != 0 && numFactorBS != 0)
		{
			dppLog << format("IR/Basis Correlations:\n");
			for(idx=0; idx<=numFactorIR*numFactorBS-1; idx++)
			dppLog << format("%lf\t", corrIrBsL[idx+1]);

			dppLog << endl;
		}
	}

	if(IS_VALID_ARRAY(resetBankIndsL))	// ignored if NULL
	{
		dppLog << format("\n ======================== Reset Bank: ========================\n");
	
		dppLog << format("Index  Reset Dates  Reset Rates\n");
		for (idx=0; idx<=(long)resetBankIndsL[0]-1; idx++)
			dppLog << format("%c      %10s   %8.4f%%\n",
			  	  resetBankIndsL[WRAP_STR_IDX(idx+1)],
			  	  GtoFormatDate(resetBankDatesL[idx+1]),
			  	  resetBankRatesL[idx+1]*1e2);

		dppLog << endl;
	}

	dppLog << endl;

	dppLog << format("DebugLevel:  %d\n", debugLevelL[1]);

	dppLog << format("\n =====================================================\n");
	
	dppLog << endl;

	status = SUCCESS;

	return(status);

}



static
double FreqCharToBasis(char freq)
{
static  char    routine[] = "_ZCFreqencyCharToBasis";
	double	basis;
	
	if (toupper(freq) == 'A')
		basis = (double)1L;
	else if (toupper(freq) == 'S')
		basis = (double)2L;
	else{
		throw KFailure("%s: invalid zero curve basis type (%c). "
                          "Only (A)nnual or (S)emi-annual are allowed.\n",
                          routine,
                          freq);
	}

	return (basis);
}


static
String	ZcIdxToName(int idx)
{
	String	cvName;

	switch (idx) {
	case 1:
		cvName = "Diffuse Curve";
		break;
	case 2:
		cvName = "Index Curve1";
		break;
	case 3:
		cvName = "Index Curve2";
		break;
	case 4:
		cvName = "Basis Curve";
		break;
	}
	
	return cvName;
}

static
String	ZcTypeCharToZcName(const char zcTypeChar)
{
	char	cvTypeName;
	String	cvName;

	cvTypeName = toupper(zcTypeChar);

	switch (cvTypeName) {
	case 'D':
		cvName = "Diffuse Curve";
		break;
	case '1':
		cvName = "Index Curve1";
		break;
	case '2':
		cvName = "Index Curve2";
		break;
	case 'B':
		cvName = "Basis Curve";
		break;
	}

	return cvName;
}


