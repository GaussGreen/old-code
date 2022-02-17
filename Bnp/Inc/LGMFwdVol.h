#ifndef LGMFwdVol_h
#define LGMFwdVol_h

/*	Compute prices of options on forward swaps in LGM and AUTOCAL	*/

/*	1.-	LGM	*/
Err LGMFwdSwaptionLGM(long exp, long start, long nfp_or_end, char *cpdStr,
                      char *basisStr, char *recPayStr, double strike,
                      char *refRateCodeStr, char *sUndPtrName, double *price);

/*	2.-	Autocal */
Err LGMFwdSwaptionATC(long exp, long start, long nfp_or_end, char *cpdStr,
                      char *basisStr, char *recPayStr, double strike,
                      char *refRateCodeStr, char *yc_name, LGM_TS *tsPtr,
                      double *price);

#endif