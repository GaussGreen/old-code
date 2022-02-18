/* ==================================================================================

   FILENAME :  srt_h_und_struct.h

   PURPOSE :  the structure to store the underlying information

   ================================================================================== */
#ifndef SRT_H_UND_STRUCT_H
#define SRT_H_UND_STRUCT_H

#include "srt_h_ts.h"
#include "srt_h_types.h"
#include "utconst.h"

/* -----------------------------------------------------------------------------
   THE GENERIC UNDERLYING STORAGE STRUCTURE AND THE UNDERLYING TYPE :
    a SrtUndDesc (stored in the Double linked list)
   ----------------------------------------------------------------------------- */

typedef struct
{
    char              underl_name[SRTBUFSZ];
    char              underl_lbl[SRTBUFSZ]; /* String corresponding to the type */
    SrtUnderlyingType underl_type;
    long              underl_ticker; /* Version number (same as in the list) */
    char*             underl_ccy;

    void* spec_desc;
} SrtUndDesc, SrtUndObj, *SrtUndPtr, *SrtUndObjPtr;

/* -------------------------------------------------------------------------- */

/* ---------------------------- INTEREST RATE --------------------------------- */

typedef struct
{
    char        mdl_lbl[SRTBUFSZ];
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts;
    char        yc_name[SRTBUFSZ];

    /*	Added Toto 29Nov1999	-	possibility to have user specs	*/
    void* spec;

} SrtIrDesc;

typedef struct
{
    int           num_of_factors;
    int           num_of_tenors;
    long*         TenorDates;
    double**      VolMatrix;
    char*         cRefRateCode;
    double*       shifts;
    BGMCalibType  CalibType;
    BGMCorrelType CorrelType;
    double**      histcorrel;
} BGMUndSpec;

typedef struct
{
    int       MaxNumPeriod;
    int       MaturityInPeriod;
    long*     TenorDates;
    double**  TSLibor;
    char*     cRefRateCode;
    double**  Correl;
    long*     IStartLiquidATM;
    long*     IEndLiquidATM;
    int       NumLiquidATM;
    int       NumMatCube;
    int       NumUndCube;
    long*     IMatCube;
    long*     IUndCube;
    double*** ATMSens;
    long*     IStartLiquidSABR;
    long*     IEndLiquidSABR;
    int       NumLiquidSABR;
    double*** AlphaSens;
    double*** RhoSens;
    double    bumpATM;
    double    bumpalpha;
    double    bumprho;
} BGMSABRUndSpec;

/* ---------------------------- EQUITY  --------------------------------- */

typedef struct
{
    char        mdl_lbl[SRTBUFSZ];
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts;
    double      spot;
    char        disc_name[SRTBUFSZ];
    char        dvd_name[SRTBUFSZ];
    char        repo_name[SRTBUFSZ];

    /*	Added Toto 29Nov1999	-	possibility to have user specs	*/
    void* spec;

} SrtEqDesc;

/* ---------------------------- FOREIGN EXCHANGE  ------------------------ */

typedef struct
{
    char        mdl_lbl[SRTBUFSZ];
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts;
    double      spot;
    char        dom_name[SRTBUFSZ];
    char        for_name[SRTBUFSZ];
    char        spread_curve_name[SRTBUFSZ];

    /*	Added Toto 29Nov1999	-	possibility to have user specs	*/
    void* spec;

} SrtFxDesc;

/* --------------------------------- BOND  ------------------------------- */

typedef struct
{
    char        mdl_lbl[SRTBUFSZ];
    SrtMdlDim   mdl_dim;
    SrtMdlType  mdl_type;
    TermStruct* ts;
    /*
            SwapDP          p;
    */
    double coupon;
    double clean_price;
    double yield;
    char   disc_name[SRTBUFSZ];
    char   repo_name[SRTBUFSZ];
    /*
            SrtBndAuxParam  aux_par;
    */

    /*	Added Toto 29Nov1999	-	possibility to have user specs	*/
    void* spec;

} SrtBndDesc;

#endif
