#include "Magnet/Magnet.h"
#include "General/General.h"

using namespace CM;

int GCD(int x, int y );

int LCM(int x,int y );

Array<int> CalcLossUnit(
     int			  NoNames,
	 const Array<int> LossGivenDefault);

Array<int> CalcLossUnit(
     int NoNames, 
     const Array<int>& LossGivenDefault,
     const Array<double>& ExpectedRecovery,
     const Array<double>& Co);

Array<double> CondLossDist(
	 const int				NoNames,
	 const Array<int>		&LossGivenDefault,
	 const Array<double>	&DefaultProb);

Array<double> CondLossDist(
	 const int				NoNames,
	 const Array<int>		&LossGivenDefault,
	 const Array<double>	&DefaultProb,
     int                    unit,
     int                    maxUnits);

Matrix <double> CondLossDistUnits(
     int					TotalLossUnits,
	 Array<double>			&LossUnitsGivenDefault,
	 Array<double>			&DefaultProb);


Matrix<double> CondLossDistMultiper(
	 const int				NoNames,
	 const Array<int>		&LossGivenDefault,
	 const Matrix<double>	&DefaultProb);
