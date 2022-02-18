#ifndef _INFLATION_SPREAD_H_
#define _INFLATION_SPREAD_H_

//------------------------------------------------------------------------------------------------------------------------------------
//
// InflationSpread.h
//
//------------------------------------------------------------------------------------------------------------------------------------
// Include Files
#include "utDates.h"
#include "utError.h"
//#include "swp_h_swapdp.h"
//#include "swp_h_DateList.h"

Err srt_InflationSpread(
    char*   csUndName,
    Date    dtExpiry,
    Date*   pdlMaturity,
    int     iNumDates,
    double* dpIntArg,
    double* dpPrice);

Err srt_CalibrateInflationSpread(
    char*   csYCname,
    double  dSpotSpread,
    double  dRevSpeed,
    double  dCorrelation,
    double  dSpreadVol,
    double  dTau,
    Date    dtToday,
    int     iNumDates,
    double* dvFwdSpread,
    double* dvLGMvol,
    Date*   dtvDates,
    double* dvCalibratedLevels);

Err srt_InflationSpreadForward(
    char* csVasicekUnd, char* csLGMUnd, Date dtExpiry, double dCorrelation, double* pdForward);

#endif  // _INFLATION_H_