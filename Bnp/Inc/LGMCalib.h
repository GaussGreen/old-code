#ifndef LGM_CALIB_H
#define LGM_CALIB_H

#include "srt_h_lgmtypes.h"

/* Finds the G(t) values by calibrating to the "1 into k" swaptions */
LGMErr FitAllGFwdSwap(LGM_TS* tsPtr, LGMCalSet* CSPtr, Date* zdatefix, double* zetafix);

#endif