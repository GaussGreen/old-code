

#include "Magnet/Magnet.h"
#include "General/General.h"

using namespace CM;



double Collar(

	 const double Loss,
	 const double TrancheLow,
	 const double TrancheHigh);


Array<double> ExpectedLossTranche(

	 const double	  Not,
	 const int		  LossUnit,
	 Array<double>	  &Attachment,
	 Array<double>	  &LossDist);

Array<double> ExpectedLossTranchelet(

	 const double	  Not,
	 const int		  LossUnit,
	 Matrix<double>	  &Attachment,
	 Array<double>	  &LossDist);



Array<double> ExpectedLossTrancheBuck(

	 double		   Not,
	 double		   LossUnit,
	 int			MaxNoUnits,
	 Array<double> &Attachment, 
	 Matrix<double> &LossDist);



Matrix<double> SpreadTranche(

	 const double  Not,
	 const int     LossUnit,
	 const int	   MaxNoUnits,
	 Array<double> &Attachment,
	 Matrix<double> &LossDist,
	 Array<double>	&Maturity,
	 double			time);