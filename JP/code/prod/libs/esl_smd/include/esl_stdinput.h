#ifndef ESL_STDINPUT_DOT_H
#define ESL_STDINPUT_DOT_H

/** NOTE: This file should be only included through 'esl_stdinput.c'
 */

#include <stdio.h>

#ifdef  __cplusplus
extern "C" {
#endif

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_alloc.h"
#include "esl_date.h"
#include "esl_util.h"

/*****	FindAndSkipComLine  *************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine (FILE         *stream,
                            char const   *Label,
                            char const   *Routine,
                            char const   *FileName);

/*****	Term_Input_W  *******************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int  Term_Input_W (T_CURVE  *t_curve    /** (O) Structure of zero curve data  */
                  ,char const* FileName /** (I) File name including extension */
		  );

/*****  Term_Check_W  *******************************************************/
/**
*       Check validity of DR Wrapper term structure inputs.
*/
int     Term_Check_W (
		T_CURVE  *t_curve  /** (I) Structure of zero curve data */
		);

/*****	BaseVol_Input_W  ****************************************************/
/**
*       Read base volatility input and check validity of input.
*/
int     BaseVol_Input_W (
            MKTVOL_DATA  *mktvol_data   /** (O) Volatility data               */
           ,char   const* FileName      /** (I) File name including extension */
	    );

/*****	SwapVol_Input_W  ****************************************************/
/**
*       Read swaption volatility input and check validity of input.
*/
int     SwapVol_Input_W (int          *NbRows,      /** (O) Volatility matrix */
                         int          *NbCol,
                         long         **Expiry,
                         long         **FwdMat,
                         double       ***VolMatrix,
						 char         **DoMoY,     /** (O) Maturities in days, 
												           months or years   */
                         char     const *FileName  /** (I) File name */
		);

/*****  MktVol_Input_W  ******************************************************/
/**
*       Utility routine converting an index name to a series of option expi-
*       ries, index maturities and volatilities taken from either the base
*       volatility curve or the swaption matrix.
*/
int     MktVol_Input_W (
		MKTVOL_DATA*   mktvol_data, /**< (O) Volatility data      */
                char const*    Index,       /**< (I) Index to calibrate   */
                T_CURVE const* t_curve,     /**< (I) Term structure data  */
                char const*    BaseVolFile, /**< (I) Base vol curve file  */
                char const*    SwapVolFile);/**< (I) Swaption matrix file */

/*****	VolDiag_Input_W  ****************************************************/
/**
 *       Read volatility input
 */
int     VolDiag_Input_W (
            long    *BaseDate    /** (O) Volatility data */
           ,int     *VolUnit
           ,int     *NbVol
           ,long    *VolDate
           ,long    *SwapSt
           ,long    *SwapMat
           ,double  *Vol
           ,char    *VolType
           ,char    *FileName    /** (I) File name including extension */
	   );

/*****  MktVol_Input_Plus_W  ************************************************/
/**
        Utility routine converting an index name to a series of option expi-
        ries, index maturities and volatilities taken from either the base
        volatility curve or the swaption matrix.
 
        Same as MktVol_Input_W except this function allows the reading 
        from ir_voldiag_i.dat
        The reading logic is as follows:
        - treats nil calibration first
        - checks if it's basevol caibration and read from basevol.dat
        - checks availability of d(f)swapvol.dat in the folder
        - if not available, then reads from ir_voldiag_i.dat
        - otherwise, always reads from d(f)swapvol.dat
*/
int MktVol_Input_Plus_W (
	MKTVOL_DATA *mktvol_data     /** (O) Volatility data               */
       ,char       *Index           /** (I) Index to calibrate             */
       ,T_CURVE    *t_curve         /** (I) Term structure data            */
       ,char       *BaseVolFile     /** (I) Base vol curve file            */
       ,char       *SwapVolFile     /** (I) Swaption matrix file           */
       ,char       *SwapVolFileSRM3 /** (I) Swaption matrix file from SRM3 */
		);

/*****	MktVol_Check_W  *****************************************************/
/**
*       Check validity of base volatility inputs.
*/
int     MktVol_Check_W (
		MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
		);

/** 
*   Print probability in output file for wrapper.
*/
void    printProbabilityInFile(
		OPT_OUT_DATA *opt_out_data
		);

#ifdef  __cplusplus
}
#endif


#endif



