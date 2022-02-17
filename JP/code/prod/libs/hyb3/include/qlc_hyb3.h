#ifndef QLC_HYB3_H
#define QLC_HYB3_H
// ESL QLib Conversion Utils

#include "qlc_utils.h"
#include <stdio.h>
#include <eslhead.h>
#include "cupsmodl.h"

#ifdef  __cplusplus
extern "C" {
#endif

/***** market *****/

int qlc_printHyb3Market(
    FILE        *f,
    T_CURVE     t_curve[2][3],/* (I) Structure of zero curve data */
    char        *baseVolFilename1,
    char        *swapVolFilename1,
    char        *baseVolFilename2,
    char        *swapVolFilename2,
    char        *FXVolFilename,
    FX_DATA     *fx_data,
    MKTVOL_DATA *mktvol_data,
    int         nbPastFX,
    long        pastFXDates[],
    double      pastFXRates[],
    HYB3_TREE_DATA *tree_data
);

int qlc_printHyb3EQMarket(
	FILE *f,
	T_CURVE *t_curve,/* (I) Structure of zero curve data */
    char   *baseVolFilename,
    char   *swapVolFilename,
    MKTVOL_DATA *mktvol_data,
    HYB3_TREE_DATA *tree_data,
    EQ_DATA *eq_data);

/***** model *****/

int qlc_printHyb3FXModel(
    FILE *f,
    HYB3_TREE_DATA *tree_data,
    const char *termPRNFileName);

int qlc_printHyb3EQModel(
    FILE *f,
    HYB3_TREE_DATA *tree_data,
    const char *termPRNFileName,
    long matDate,
    const char *eqSpotVolOverride);

void qlc_printSimpleEquity(
    FILE *f, 
    const char *tag, 
    EQ_DATA *eq_data, 
    long valueDate);

#ifdef  __cplusplus
}

namespace qlc {

struct FxVol : public BaseCounter {
	const char* name;
	FXVOLATILITY_DATA *fxVol;
	FX_DATA *fx;
	FxVol(const char* name, FXVOLATILITY_DATA *fxVol, FX_DATA *fx) : name(name), fxVol(fxVol), fx(fx) {}
	void print(FILE *f, const std::string &tag);
};


struct FxAsset : public BaseCounter {
	const char *name;
	UntweakableYC *yc1;
	UntweakableYC *yc2;
	FxVol *fvh;
    int nbPastFixings;
    long *pastDate;
    double *pastRate;
	FxAsset(const char *name, UntweakableYC *yc1, UntweakableYC *yc2, FxVol *fvh) 
	: name(name), yc1(yc1), yc2(yc2), fvh(fvh) {}
	void print(FILE *f, const std::string &tag);
};

} /* end of namespace qlc */

#endif

#endif
