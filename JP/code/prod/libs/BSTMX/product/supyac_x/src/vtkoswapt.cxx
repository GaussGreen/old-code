/****************************************************************
 * Purpose:	Test Driver
 * Author:	C. Daher
 *
 ****************************************************************/
#include "kstdinc.h"

extern	"C" {
#include "cgeneral.h"
#include "bastypes.h"
#include "argdef.h"
#include "cerrsup.h"
#include "mainlot.inc"		/* Includes main() */

#include "drlvtype.h"		/* DrlLilVectLogging() */
};

#include "vtkoswap.h"
#include "kstdinc.h"


/*--------------------------------------------------------------
 * callRoutine(): Called by driver.
 */
extern	"C" {
void callRoutine(FILE *fp)
{
#define	NOUTPUT	4
	int	i,
		errCode = 0;
	double	outVal[NOUTPUT+1];

	long	errorBufferSize;
	char	*errorBuff = NULL;


	outVal[0] = (double) NOUTPUT;
	for (i=1; i<=NOUTPUT; i++) outVal[i] = 0.;


	fprintf(fp, "Calling Test Driver ...\n");	


	/*
	 * Call the Lotus wrapper.
	 */
	errCode = DVTreeKnockOutBasisSwapL(
			/* --- Receive Leg 	*/
	(TDate  *)	TestArgCacheGet(  0),
	(TDate  *)	TestArgCacheGet(  1),
	(TDate  *)	TestArgCacheGet(  2),
	(TDate  *)	TestArgCacheGet(  3),
	(TDate  *)	TestArgCacheGet(  4),
	(double *)	TestArgCacheGet(  5),
	(double *)	TestArgCacheGet(  6),
			/* --- Pay Leg 		*/
	(TDate  *)	TestArgCacheGet(  7),
	(TDate  *)	TestArgCacheGet(  8),
	(TDate  *)	TestArgCacheGet(  9),
	(TDate  *)	TestArgCacheGet( 10),
	(TDate  *)	TestArgCacheGet( 11),
	(double *)	TestArgCacheGet( 12),
	(double *)	TestArgCacheGet( 13),
			/* --- Receive and Pay Legs  	*/
	(char   *)	TestArgCacheGet( 14),
	(char   *)	TestArgCacheGet( 15),
			/* --- Index receive, pay & ko	*/
	(double *)	TestArgCacheGet( 16),
	(double *)	TestArgCacheGet( 17),
	(char   *)	TestArgCacheGet( 18),
	(char   *)	TestArgCacheGet( 19),
	(char   *)	TestArgCacheGet( 20),
	(char   *)	TestArgCacheGet( 21),
	(char   *)	TestArgCacheGet( 22),
			/* --- 	KIO		*/
	(char   *)	TestArgCacheGet( 23),
	(TDate  *)	TestArgCacheGet( 24),
	(TDate  *)	TestArgCacheGet( 25),
	(TDate  *)	TestArgCacheGet( 26),
	(double *)	TestArgCacheGet( 27),
	(double *)	TestArgCacheGet( 28),
	(double *)	TestArgCacheGet( 29),
			/* --- Zero curve info 	*/
	(TDate  *)	TestArgCacheGet( 30),
	(char   *)	TestArgCacheGet( 31),
	(char   *)	TestArgCacheGet( 32),

	(TDate  *)	TestArgCacheGet( 33),
	(double *)	TestArgCacheGet( 34),
	(TDate  *)	TestArgCacheGet( 35),
	(double *)	TestArgCacheGet( 36),
	(TDate  *)	TestArgCacheGet( 37),
	(double *)	TestArgCacheGet( 38),
	(TDate  *)	TestArgCacheGet( 39),
	(double *)	TestArgCacheGet( 40),

	(double *)	TestArgCacheGet( 41),
			/* ---	Vol Curve info 	*/
	(TDate  *)	TestArgCacheGet( 42),
	(TDate  *)	TestArgCacheGet( 43),
	(int    *)	TestArgCacheGet( 44),
	(double *)	TestArgCacheGet( 45),

	(TDate  *)	TestArgCacheGet( 46),
	(double *)	TestArgCacheGet( 47),
	(char   *)	TestArgCacheGet( 48),
			/* -- Model Info	*/
	(double *)	TestArgCacheGet( 49),
	(double *)	TestArgCacheGet( 50),
	(double *)	TestArgCacheGet( 51),
	(double *)	TestArgCacheGet( 52),
	(double *)	TestArgCacheGet( 53),

			/* --- Reset Bank	*/
	(char   *)	TestArgCacheGet( 54),
	(TDate  *)	TestArgCacheGet( 55),
	(double *)	TestArgCacheGet( 56),
			/* --- Other 		*/
	(int    *)	TestArgCacheGet( 57),
			/* --- Output 		*/
	outVal
	);


	if (errCode == SUCCESS) {
	    fprintf(fp, "\nOUTPUT: SUCCESS\n");
	    DrlLilVectLogging(fp, DRL_FLOAT_L, (void*)outVal, "OUT_PRICE");
	} else {
	    errorBufferSize = DppErrMsgBufferSize();
	    errorBuff = new char[errorBufferSize];
 
	    if(DppErrMsgBufferRetrieve(errorBufferSize,
					errorBuff) != SUCCESS)
	    	fprintf(fp, "Failed to retieve error message from buffer\n");
            else
		fprintf(fp, "%s\n", errorBuff);
	
	    fprintf(fp, "\nOUTPUT: FAILED\n");
	
	    delete [] errorBuff;
	}


	return;
}		


/*--------------------------------------------------------------
 *
 */
ARG_DEF g_ArgInfoArray[]=
{

	{"RESET_DATES_1",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"RESET_EFF_DATES_1",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ACC_START_DATES_1",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ACC_END_DATES_1",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"PAY_DATES_1",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"SPREADS_1",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"NOTIONALS_1",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"RESET_DATES_2",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"RESET_EFF_DATES_2",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ACC_START_DATES_2",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ACC_END_DATES_2",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"PAY_DATES_2",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"SPREADS_2",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"NOTIONALS_2",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"PAY_DAY_COUNTS",	AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"STUB_CONVS",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */

	{"INDEX_WEIGHTS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"INDEX_SPREADS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"INDEX_FREQS",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"INDEX_MATURITIES",	AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"INDEX_DCC",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"INDEX_CURVES",	AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"DISC_INDEX_CURVES",	AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */

	{"KIO_TYPES",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"KIO_DATES",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"KIO_EFF_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"KIO_SETTLE_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"KIO_LO_BARRIERS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"KIO_HI_BARRIERS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"KIO_REBATES",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"ZC_BASE_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ZC_FREQS",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"ZC_DCCS",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */

	{"ZC_DATES_D",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ZC_RATES_D",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"ZC_DATES_1",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ZC_RATES_1",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"ZC_DATES_2",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ZC_RATES_2",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"ZC_DATES_B",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"ZC_RATES_B",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"ZC_BASIS_INFO",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */


	{"IR_VOL_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"IR_VOL_MATS",		AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"IR_VOL_FREQ",		AT_ARRAY(AT_INT)},		/* ( 0) */
	{"IR_VOL_RATES",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"BS_VOL_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"BS_VOL_RATES",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"BS_VOL_TYPE",		AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */

	{"IR_MR_PARAMS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"IR_SMILE_PARAMS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"BS_MR_PARAMS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"BS_SMILE_PARAMS",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */
	{"CORR_IR_BS",		AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"RESET_BANK_INDEX",	AT_ARRAY(AT_CHAR_BLOCK)},	/* ( 0) */
	{"RESET_BANK_DATES",	AT_ARRAY(AT_DATE)},		/* ( 0) */
	{"RESET_BANK_RATES",	AT_ARRAY(AT_DOUBLE)},		/* ( 0) */

	{"INT_SCALARS",		AT_ARRAY(AT_INT)},		/* ( 0) */
};

long g_NumArgs = sizeof(g_ArgInfoArray)/sizeof(ARG_DEF);

}	// extern	"C"
