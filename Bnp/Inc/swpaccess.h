/* ================================================================================
   
   FILENAME: SwpAccess.h                                          
 
   FUNCTION: Header for accessing functions in Sort for Swaps Pricing
   
   ================================================================================ */

#include "swp_h_all.h"

 
#ifndef SWPACCESS_H
#define SWPACCESS_H

char *SwpInitExternalFct(
		DiscFuncType       disc_func,   /* NOT USED FOR THE MOMENT: can be NULL */
		SpreadFuncType     spread_func,
		RateInfoFuncType   rate_info_func,
		VolFuncType        vol_func,
		SABRVolFuncType    sabrvol_func,
		FixingFuncType     fixing_func);

char *SwpInitYC(char *ycName, char *ycType, char *ccyName, long today, long spotdate);


char *SwpSwapUnwind(long value_date, long start, long end, long today,
                     char *cpd, char *basis, double strike, char *firstShortFull,
                     double initial_not, double final_not, char *liborCpdStr,
                     char *LiborBasis, double liborFix, 
                     char *undName, double *pv);

char *SwpSwapCalc(long value_date, long start, long  nfp_or_end, long today,
                   char *cpdStr, char *basisStr, double strike, 
                   double initial_not, double final_not, char *info_message,
                   char *undName, double *answer);

char *SwpSwapPV(long start, long  nfp_or_end,
                 char *cpdStr, char *basisStr, double strike, 
                 char *dnctName,  char *refRateCode, double *pv); 

char *SwpMargin(double pv, long start, long  nfp_or_end,
                 char *cpdStr, char *basisStr,
                 char *ycName, double *margin);

char *SwpFwdRate(long start, long  nfp_or_end,
                 char *cpdStr, char *basisStr,
                 char *ycName, char * refRateCode, double *fwdRate);

#endif