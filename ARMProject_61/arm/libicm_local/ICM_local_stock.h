#ifndef ICM_LOCAL_STOCK_H
#define ICM_LOCAL_STOCK_H

#include <ARM\libarm\ARM_result.h>

extern long ICMLOCAL_STOCK   (double	FirmVolIn,
					  double	DividendYieldIn,
					  long TreeModelId,
					  double SpotIn,
					  ARM_result& result,
					  long objId = -1);

extern long ICMLOCAL_STOCKCALLOPTION  (long	StockId,
					  double	StrikePriceIn,
					  double		FirstDateIn,
					  double		LastDateIn,
					  ARM_result& result,
					  long objId = -1);

#endif	// ARM_LOCAL_STOCK_H

// EOF %M%