/*	SrtAccessLon
        Files from London for SrtAccess	*/

#ifndef SRTACCESSLON_H
#define SRTACCESSLON_H

#include "Fx3FBetaDLMUtil.h"
#include "MCEBOptimisation.h"
#include "grf_h_all.h"
#include "srt_h_all.h"

char* SrtGrfnFx3DFxBetaDLMTree(
    char*                      und3dfx,
    int                        numeventdates,
    long*                      eventdates,
    long                       tableauRows,
    long                       tableauCols,
    char***                    tableauStrings,
    int**                      tableauMask,
    long                       auxWidth,
    long*                      auxLen,
    double**                   aux,
    FxBetaDLM_GRFNNumerParams* Params,
    double*                    barrier,
    int*                       bar_col,
    long                       num_stp,
    int*                       num_prod,
    double**                   prod_val);

char* SrtGrfn3DFxBetaDLMMc(
    char*    und3dfx,
    int      numeventdates,
    long*    eventdates,
    long     tableauRows,
    long*    tableauCols,
    char***  tableauStrings,
    int**    tableauMask,
    long     auxWidth,
    long*    auxLen,
    double** aux,

    // Params
    FxBetaDLM_GRFNNumerParams* Params,

    // for Optimisation of exercise boundary
    int        do_optimisation,
    int*       optimise,
    MCEBPARAMS params,
    long*      resRows,

    long      num_paths,
    int*      nb_prod,
    double*** prod_val);

#endif
