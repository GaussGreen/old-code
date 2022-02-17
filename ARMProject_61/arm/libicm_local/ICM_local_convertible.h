#ifndef ARM_LOCAL_CONVERTIBLE_H
#define ARM_LOCAL_CONVERTIBLE_H

#include <ARM\libarm\ARM_result.h>


extern long ICMLOCAL_CBOPTION(double		FirstDateIn,
							  double		LastDateIn,
							  int			OptionTypeIn,
							  double		StrikeIn,
							  int			StrikeTypeIn,
							  int			AccruedOnExIn,
							  int			BarrierTypeIn,
							  double		BarrierStrikeIn,
							  int			BarrierStrikeType,
							  ARM_result&	result,
							  long			objId = -1);


extern long ICMLOCAL_CONVERTIBLE(long		DefBondId,
								 long		StockId,
								 double		YieldIn,
								 const VECTOR<long>& CBOptionIds,
								 ARM_result&result,
								 long		objId = -1);

extern long ICMLOCAL_CBCALLABLEASW(long		ConvertId,
								 const VECTOR<long>& CBOptionIds,
								 double		SwapStart,
								 long		irIndexId,
								 double		Spread,
								 double		EarlyCall,
								 double		RecallSpread,
								 bool		hasEuribidClause,
								 bool		CallIfITM,
								 ARM_result&result,
								 long		objId = -1);



#endif