#ifndef HYB3_MAKRET_H
#define HYB3_MAKRET_H

#include <eslhead.h>
#include "cupslib.h"

#ifdef  __cplusplus
extern "C" {
#endif

//-----------------------------------------------------------------------------
//                            Market environment
//-----------------------------------------------------------------------------


typedef struct Hyb3Market_
{
    SWAPVOL_DATA  mSwapVol[2];
    BASEVOL_DATA  mBaseVol[2];

    int mNbZeroCurve[2];
    T_CURVE       mZeroCurve[2][3];

    int                   mModelParamsSupplied[2];  // true/false 
    MODELPARAMETERS_DATA  mModelParameter[2];

    // FX specifics 
    FXVOLATILITY_DATA     mFXVolatility;
    CORRELATION_DATA      mCorrelation;
    FXSMILE_DATA          mFXSmile;



} Hyb3Market;


/*************************************************
 *
 *  Read Raw Market Data
 *
 *  These functions read the "raw" wrapper type
 *  market environment data (<files>.dat) and 
 *  store in market data structures, independently
 *  of deal specifications
 *
 *  Interface functions read different market data
 *  for different wrapper types is....
 *
 *  HYB3Read<Wrapper type>MarketW
 *
 * 
 **************************************************/

/***********************
 *
 *  Interface functions
 *
 ***********************/

/* DRWrapper Type (FX and 2 interest rates */
int HYB3ReadType3MarketW(
    Hyb3Market   *market);


/**********************************
 *
 *  Internal Functions
 *  should only be called from 
 *  interface functions
 *
 **********************************/
/* reads fxvolatility.dat and fxsmile.dat */
int HYB3ReadFXEnvW(
    FXVOLATILITY_DATA  *FXVol,
    FXSMILE_DATA       *FXSmile,
    char const         *FXVolFilename,
    char const         *FXSmileFilename);

int HYB3ReadCorrelation(
    CORRELATION_DATA   *correl,
    char const         *filename);

int HYB3ReadModelParams(
    MODELPARAMETERS_DATA *modelParams,
    int                  NbFactor,        // always 1 factor IR for hyb3
    int                  *paramsSupplied, // (O) true/false
    int                  readSmileParams, // (I) true/false
    char const          *filename);



/********************************************
 *
 *  Populate underlying HYB3 tree structures
 *
 *  Takes model parameters and overwrites
 *  from the deal file and populates
 *  MKTVOL_DATA, FX_DATA, TREE_DATA
 *
 ********************************************/

int HYB3ParamOverwrites(
    MKTVOL_DATA           *mktvol_data,              // (O) Volatility data
    HYB3_TREE_DATA        *tree_data,                // (O) Tree data structure
    int                   NbFactor,                  // (I) Number of factors
    MODELPARAMETERS_DATA* modelParameters_data,
    int                   modelParametersSupplied,
    char                  OverWriteString[5][MAXBUFF]);

int HYB3CorrelOverwrites(
    FX_DATA            *fx_data,                    // (O) Fx data
    char               OverWriteString[3][MAXBUFF], // (I) Overwrite strings
    CORRELATION_DATA   *correlation_data);

int HYB3FXSmileOverwrites(
    FX_DATA            *fx_data,                 // (O) Fx data
    char               OverWriteString[MAXBUFF], // (I) Owrite str
    FXVOLATILITY_DATA  *Volatility_data,         // (I) FXVols file name
    FXSMILE_DATA       *Smile_data,              // (I) FX smile file name
    char               NbFXSmileOWS[MAXBUFF],    // (I)
    char               FXSmileParamOWS[MAXNBDATE][MAXBUFF]);

#ifdef  __cplusplus
}
#endif

#endif
