
#include "Magnet/Magnet.h"
#include "General/General.h"

using namespace CM;

void BetaTweakSingle(
	 double		   &Beta,
	 double		   Tweakfact);

void BetaTweak(
	 Array<double> &Beta,
	 double		   Tweakfact);

void BetaTweakV(
	 Array<double> &Beta,
	 double		   Tweakfact,
	 Array<double> &BetaTweak);

Array<double> CorrelationSkewSimple(
	 int		   NoNames,
	 int		   NotionalperName,
	 Array<double> &ExpectedLosses,
	 double		   Beta,
	 Array<double> &Attachment,
	 int		   LossperName,
	 double		   DefaultProb);

Array<double> CorrelationSkew(
	 const int	    Notional,
	 Array<double>	&ExpectedLosses,
	 double			Beta,
	 Array<double> &Attachment,
	 Array<int>		&Loss,
	 Array<double> &DefaultProb);


Array<double> CorrelationSkewGeneral(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Array<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb);

Array<double> CorrelationSkew2(
	 const int	    Notional,
	 Array<double>	&ExpectedLosses,
	 double			Beta,
	 Matrix<double> &Attachment,
	 Array<int>		&Loss,
	 Array<double> &DefaultProb);


Array<double> CorrelationSkewGeneral2(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Matrix<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb);


Array<double> CorrelationSkewGeneralSecant(
	 const int	   Notional,
	 Array<double> &ExpectedLosses,
	 Array<double> &Beta,
	 Array<double> &Attachment,
	 Array<int>	   &Loss,
	 Array<double> &DefaultProb);

