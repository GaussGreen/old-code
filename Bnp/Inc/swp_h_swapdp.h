/* =========================================================================
   FILENAME: swp_h_swapdp.h
   FUNCTION: swp_f_initSwapDP
   PURPOSE:  to initialise a SwapDP from start, end_nfp, compString, basisStr
   ========================================================================= */

#ifndef SWP_H_SWAPDP_H
#define SWP_H_SWAPDP_H

#define     MAXREFRATENAMESIZE     32

typedef enum{
	FWD,
	BKWD,
	LASTSWAPDATEDIR} 
SwapDateDir,SrtSwapDateDir;

typedef struct {
	SrtBasisCode      basis_code;
	Date              start;
	Date              end;
	Date              first_full_fixing;
	int               nfp;
	SrtCompounding    compd;
	SrtSwapDateDir    direction;
	int               spot_lag;
	} SwapDP, SrtSwapDP;


Err swp_f_initSwapDP(	long start, 
						long end_nfp, 
						String compStr, 
                  		String basisStr,
						SwapDP *sdp);

Err swp_f_setSwapDP(	long start, 
						long end_nfp, 
						SrtCompounding comp, 
                  		BasisCode basis,
						SwapDP *sdp);

#endif
