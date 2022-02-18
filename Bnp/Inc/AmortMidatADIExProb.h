// prevent multiple inclusions
#pragma once

////////////////////
//	warnings
#pragma warning( \
    disable : 4786)  //"identifier was truncated to '255' characters in the debug information"

// NB: force warnings for unused arguments and local parameters
#pragma warning(1 : 4100 4101)

#include "AmortMidatADIGrid.h"
#include "AmortMidatADIModel.h"
#include "AmortMidatADIUtils.h"

const char* Compute_ExProb_Backward(
    const double*   pdEx_Begin,
    const double*   pdEx_End,
    _Grid*          pGrid,
    const int*      pnSizeX3_Begin,
    const int*      pnSizeX1_Begin,
    const double*** pppdExIndicator_Begin,
    const int*      pnExIndicator_RBegin,
    const int*      pnExIndicator_REnd,
    const int*      pnExIndicator_CBegin,
    const int*      pnExIndicator_CEnd,
    const _Model*   pModel,
    // results
    double* pdExProb);

const char* Compute_ExProb_Forward_(
    const double*   pdT_Begin,
    const double*   pdT_End,
    _Grid*          pGrid,
    const int*      pnSizeX1_Begin,
    const int*      pnSizeX3_Begin,
    const double*   pdVarX1_Begin,
    const double*   pdVarX3_Begin,
    const double*** pppdExIndicator_Begin,
    const int*      pnExIndicator_RBegin,
    const int*      pnExIndicator_REnd,
    const int*      pnExIndicator_CBegin,
    const int*      pnExIndicator_CEnd,
    const _Model*   pModel,
    // results
    double* pdExProb_Begin);

const char* Compute_ExProb(
    const double*   pdEx_Begin,
    const double*   pdEx_End,
    _Grid*          pGrid,
    const int*      pnSizeX3_Begin,
    const int*      pnSizeX1_Begin,
    const double*   pdVarX3_Begin,
    const double*   pdVarX1_Begin,
    const _Model*   pModel,
    const double*** pppdExIndicator_Begin,
    const int*      pnExIndicator_RBegin,
    const int*      pnExIndicator_REnd,
    const int*      pnExIndicator_CBegin,
    const int*      pnExIndicator_CEnd,
    // results
    double* pdExProb,
    int     nUse_Backward);
