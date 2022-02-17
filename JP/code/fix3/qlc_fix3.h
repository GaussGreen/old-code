#ifndef QLC_FIX3_H
#define QLC_FIX3_H
// ESL QLib Conversion Utils

#include <stdio.h>
#include <eslhead.h>
#include "qlc_utils.h"
#include "fix123.h"

#ifdef  __cplusplus
extern "C" {
#endif

void qlc_printFix3Market(
	FILE *f,
	T_CURVE t_curve[3],/* (I) Structure of zero curve data */
	char   *baseVolFilename,
	char   *swapVolFilename,
	FIX3_TREE_DATA *tree_data,
    MKTVOL_DATA *mktvol_data
);

void qlc_printFix3Param(
    FILE *f,
    char *volCalibIndex,
    FIX3_TREE_DATA *tree_data,
    char *termPRNFileName,
    int cetSmoothing,
    int useFix3CB,
    int incValDatePricingEvents);

#ifdef  __cplusplus
}
#endif

#endif
