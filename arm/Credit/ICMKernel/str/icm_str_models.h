#ifndef STR_MODEL_H
#define STR_MODEL_H

#include <vector>
// 17783 using namespace std;

#include "ICMKernel\mod\modelmulticurves.h"
#include "ICMKernel\crv\icm_default_str.h"

// --------------------------
// Creation d'un model
// --------------------------

ICM_DefCurvStr* AverageCurveStr(const std::vector<const ICM_DefaultCurve*> &TabCurve );
ICM_DefaultCurve* AverageCurvePWC(ICM_DefaultCurve** TabCurve, int NbCurve);

#endif // STR_MODEL_H