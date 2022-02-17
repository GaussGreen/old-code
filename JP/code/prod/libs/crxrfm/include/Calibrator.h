#include "Magnet/Magnet.h"
#include "General/General.h"

using namespace CM;

namespace CM {

enum CalibrationMode {
    QUASI_NEWTON,
    LM,
    NON_LINEAR
};

extern SymbolMap<CalibrationMode> CalibrationModeSymbolMap;

CM_MANAGED_ENUM_DATA( CalibrationMode, CalibrationModeSymbolMap)

}

Array<double> Calibrate2Fact(
        CalibrationMode     mode,
        int                 n,
        Array<double>&      seeds,
        Array<int>&         masks,
        Array<double>&      targetSpreads, 
        Array<int>&         loss,
        Array<double>&      defaultProb,
        int                 notional,
        Array<double>&      attachment,
        double              noIntervals,
        double              lbound,
        double              ubound,
        double              mbound);

Array<double> Calibrate3FactHomogRR(
       CalibrationMode      mode,
	   int					n,
	   Array<double>		&seeds,
       Array<int>           &masks,
       Array<double>        &targetSpreads,
       double               r2,
	   Array<double>		&Recovery,
	   int					indNotional,
	   double				DefaultProb,
	   Array<double>		&Attachment,
	   double				NoIntervals,
	   double				lbound,							
	   double				ubound,
	   double				mbound);
