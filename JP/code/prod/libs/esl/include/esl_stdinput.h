#ifndef ESL_STDINPUT_DOT_H
#define ESL_STDINPUT_DOT_H

#include <stdio.h>

#include "esl_macros.h"
#include "esl_types.h"
#include "esl_error.h"
#include "esl_alloc.h"
#include "esl_date.h"
#include "esl_util.h"

#ifdef  __cplusplus
extern "C" {
#endif

int     FindAndSkipSectionLine (int, FILE *, char const*, char const*,char const *);


/*****	FindAndSkipComLine  *************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int     FindAndSkipComLine (FILE         *stream,
                            char const   *Label,
                            char const   *Routine,
                            char const   *FileName);

int     FindAndSkipComLine_2 (int, FILE *,char const *,char const *,char const *);

/*****	Term_Input_W  *******************************************************/
/**
*       Read term structure input for DR Wrapper and check validity of input.
*/
int  Term_Input_W (T_CURVE* t_curve    /** (O) Structure of zero curve data  */
                  ,char const* FileName); /** (I) File name including extension */


                    
/*****	Term_Input_New_W  ****************************************************/
/**
*       Read term structure input from ir_curve and check validity of input.
*/
int  Term_Input_New_W (T_CURVE*   crv,            /**< (O) Structure of zero curve data  */
                        char const*      fileName /**< (I) File name including extension */ );

/*****  Term_Check_W  *******************************************************/
/**
*       Check validity of DR Wrapper term structure inputs.
*/
int     Term_Check_W (
		T_CURVE const*    t_curve  /** (I) Structure of zero curve data */
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
		MKTVOL_DATA*                    mktvol_data, /**< (O) Volatility data      */
                char const*             Index,       /**< (I) Index to calibrate   */
                T_CURVE const*    t_curve,     /**< (I) Term structure data  */
                char const*             BaseVolFile, /**< (I) Base vol curve file  */
                char const*             SwapVolFile);/**< (I) Swaption matrix file */

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
           ,char const* FileName    /** (I) File name including extension */
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
	MKTVOL_DATA*             mktvol_data     /** (O) Volatility data               */
    ,char const*             Index           /** (I) Index to calibrate             */
    ,T_CURVE const*    t_curve         /** (I) Term structure data            */
    ,char const*             BaseVolFile     /** (I) Base vol curve file            */
    ,char const*             SwapVolFile     /** (I) Swaption matrix file           */
    ,char const*             SwapVolFileSRM3 /** (I) Swaption matrix file from SRM3 */
		);


/*****  MktVol_Input_New_W  ************************************************/
/**
        Utility routine converting an index name to a series of option expi-
        ries, index maturities and volatilities taken from either the base
        volatility curve or the swaption matrix.
 
        Same as MktVol_Input_W except this function  reads 
        from ir_voldiag_i.dat.  Also, vol point skipping flag is optionally
        in ir_voldiag_i.dat (as opposed to embedded in some calibration index
        specification).
*/
int MktVol_Input_New_W (
	MKTVOL_DATA*            mktvol_data     /** (O) Volatility data                */
    ,T_CURVE *              t_curve         /** (I) Term structure data            */
    ,char const*            VolFileSRM3     /** (I) Swaption matrix file           */
		);

/*****	MktVol_Check_New_W  *****************************************************/
/**
*       Check validity of  volatility inputs.
*/
int     MktVol_Check_New_W (
		MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
		);

/*****VolDiag_Check_W  *****************************************************/
/**
*       Check validity of  volatility inputs.
*/
int     VolDiag_Check_W (long    BaseDate,        /* Volatility data */
                         int     VolUnit,
                         int     NbVol,
                         long    *VolDate,
                         long    *SwapSt,
                         long    *SwapMat,
                         double  *Vol,
                         char    *VolType);    


/*****	MktVol_Check_W  *****************************************************/
/**
*       Check validity of base volatility inputs.
*/
int     MktVol_Check_W (
		MKTVOL_DATA  *mktvol_data /** (O) Volatility data */
		);


int MktInfo(
         MKTVOL_DATA*       mktvol_data,       /** (O) Volatility data         */
         int*               diffusionCurveIdx, /** (0) Index of diffused curve */
         int*               discountCurveIdx,  /** (0) Index of discount curve */
         T_CURVE *          t_curve,           /** (I) Term structure data     */
         char const*        IRInfoFileName);   /** (I) IR info file            */

int SummaryInfo(
         T_CURVE *          t_curve,         /** (I) Term structure data    */
         char const*        FileName);        /** (I) Summary.dat file       */

int ModelInfo(const char* fileName,             /** (I) Name of model info file (e.g. model.dat)           */
              const char* expectedEngineName,   /** (I) Expected engine type to find in model info file    */
              char**      modelChoiceString);   /** (O) Uppercase model choice string.  Callers must free! */

/*****	MktVol_PrintStructure  **********************************************/
/**
*       Print out all fields of the structure to provided stream 
*       to assist debugging
*/
int     MktVol_PrintStructure (
        FILE*        stream,
		MKTVOL_DATA  *mktvol_data
		);

/** 
*   Print probability in output file for wrapper.
*/
void    printProbabilityInFile(
		OPT_OUT_DATA *opt_out_data
		);

void EslPrintToStringZeroCurve2(T_CURVE const* crv, char** curveString);

void    EslPrintZeroCurve(T_CURVE const* crv,  
                          FILE*          file);
void    EslPrintZeroCurveAndDiscFactors(T_CURVE const* crv,  
                                        int includeDiscFactors,     /* (I) TRUE or FALSE */
                                        FILE*          file);
/* MAW 2.0 compliant version */
void EslPrintZeroCurve2(T_CURVE const* crv, FILE* file);
#ifdef  __cplusplus
}
#endif






#ifdef  __cplusplus
// Market environment input routines - stream versions

#include <iosfwd>
#include <string>

/// Read zero curve input from the stream. Function expects the stream
/// to contain three curves each optionally enclosed in
/// <zero></zero> tags.
/// @param  crv     array of three zero curves
/// @param  is      input stream
void EslTermInput(T_CURVE       crv[3],  ///< array of zero curves
                  std::istream& is);     ///< input stream

/// Read zero curve input from files. Function expects to find
/// 'zero.dat', 'disczero.dat' and optionally 'riskzero.dat' files
/// in the current directory.
/// @param  crv     array of three zero curves
void EslTermInput(T_CURVE       crv[3]); ///< array of zero curves - file input

/// Read volatility data from the stream. Function expects the stream
/// to contain basevol and swapvol data optionally enclosed in
/// <basevol></basevol> and <swapvol></swapvol> tags.
/// @param  mktvol  market vol data structure
/// @param  index   name of the calibration index
/// @param  t_curve zero curve data structure
/// @param  is      input stream
void EslMktVolInput(MKTVOL_DATA*      mktvol,  ///< volatility data
                 char const*       index,   ///< index to calibrate
                 T_CURVE const*    t_curve, ///< term structure data
                 std::istream&     is);     ///< input stream

/// Read volatility data from the input files. Function expects to find
/// 'basevol.dat' and 'swapvol.dat' in the current directory.
/// @param  mktvol  market vol data structure
/// @param  index   name of the calibration index
/// @param  t_curve zero curve data structure
void EslMktVolInput(MKTVOL_DATA*      mktvol,  ///< volatility data
                    char const*       Index,   ///< index to calibrate
                    T_CURVE const*    t_curve);///< term structure data

/// Read the <tag> line or the title line from the input stream
/// @param  tag     tag name without enclosing angle brackets
/// @param  is      input stream
/// @return         true if the tag has been found
void EslHeaderInput(std::string const& tag, std::istream& is);

/// Read the </tag> line line from the input stream. Should only be
/// invoked after earlier call to EslHeaderInput returns true as it
/// attempts to read lines until the matching tag is found.
/// @param  tag     tag name without enclosing angle brackets
/// @param  is      input stream
/// @return         true if the tag has been found
void EslFooterInput(std::string const& tag, std::istream& is);

#endif




#endif
