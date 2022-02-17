

#include "Magnet/Magnet.h"
#include "General/General.h"
#include "imsl.h"
#include "imsls.h"




using namespace CM;


double UncondSurvProb(

		 const double   c,
		 const double	mean,
		 const double   stddev,
		 const double   betaH,
		 const double   betaL,
		 const double	threshold);

double UncondSurvProbGen(

		 double   c,
		 double   mean,
		 double   stddev,
		 double   beta,
		 double   DispGen,
		 double   Disp1,
		 double   Disp2,
		 double   theta1,
		 double   theta2);


double UncondSurvProbPolExp(

		double   c,
		double   mean,
		double   stddev,
		double   a0,
		double   a1,
		double   a2,
		double   aexp,
		double   kexp);


double UncondSurvProbDiscrete(

		double			c,
		double		    mean,
		double			stddev,
		Array<double>	&prob,
		Array<double>	&theta);



double UncondSurvProbNfact(

		double			c,
		double			mean,
		double			stddev,
		Array<double>   beta,
		Array<double>	theta);



Matrix<double> OrderBetaRecTh(

		Array<double>   &Beta,
		Array<double>	&theta,
		Array<double>   &Recovery,
		Array<double>	&thetaR);


double UncondSurvProbNBetaNRec(

		double			c,
		double			mean,
		double			stddev,
		Matrix<double>  &Param);



// Calculates loss distribution for portfolio with homogeneous beta


Array<double> LossDistConstBetaHomogeneousSpread(

	 const int      NoNames,
	 double			Beta,
	 double			DefaultProb);



Array<double> LossDistConstBeta(

	 const int      NoNames,
	 double			Beta,
	 Array<int>     &LossGivenDefault,
	 Array<double>  &DefaultProb);



Array<double> LossDist(

	 const int    NoNames,
	 Array<double> &Correlation,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb);




Array<double> LossDistRandomFactorPolExp(

	 const int	   NoNames,
	 double		   a0,
	 double		   a1,
	 double		   a2,
	 double		   aexp,
	 double		   kexp,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb);



Array<double> LossDistRandomNFactorsHomogenSimp(

	 const int	    NoNames,
	 Array<double>  beta,
	 Array<double>  theta,
	 double			DefaultProb,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound);


Array<double> LossDistRandomNFactorsSimp(

	 const int	   NoNames,
	 Array<double> beta,
	 Array<double> theta,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound);


Array<double> LossDistRandomDiscrete(

	 const int	   NoNames,
	 Array<double> probMarket,
	 Array<double> theta,
	 Array<int>    &LossGivenDefault,
	 Array<double> &DefaultProb);


Matrix<double> LossDistRandomFactorGeneralMultper(

	 const int	   NoNames,
	 double		   beta,
	 double		   DispGen,
	 double		   Disp1,
	 double		   Disp2,
	 double		   theta1,
	 double		   theta2,
	 Array<int>    &LossGivenDefault,
	 Matrix<double> &DefaultProb);



Array<double> UncondMomentsNfactor(

	 Array<double>  beta,
	 Array<double>  theta);



Matrix<double> LossDistRandomNFactorsSimpsBuck(

	 const int	   NoNames,
	 double		   fracG1,
	 double		   fracG2,
	 double		   probG1,
	 double		   probG2,
	 double		   probJ,
	 Array<double> beta,
	 Array<double> theta,
	 Array<int>    &Notional,
	 Array<double> &Recovery,
	 Array<double> &DefaultProb,
	 double		    NoIntervals,
	 double			lbound,
	 double			ubound,
	 double			mbound,
	 double			NoUnits);



Array<double> LossDistRandomNFactorsHomogenSimpRR(

	 const int	    NoNames,
	 Array<double>  Beta,
	 Array<double>  Theta,
	 double			ProbJ,
	 double			DefaultProbT,
	 double		    NoIntervals,
	 Array<double>	Recovery,
	 Array<double>	&ThresholdR,
	 double			lbound,
	 double			ubound,
	 double			mbound);


Array<double> LossDistRandomNFactorsSimpRR(
	   int					n,
	   Array<double>		&BetaCo,
	   Array<double>		&Theta,
       Array<double>		&RecoveryCo,
	   Array<double>		&ThresholdR,
       Array<double>        &Beta,
       Array<double>        &ExpectedRecovery,
	   Array<int>   		&Loss,
	   Array<double>		&DefaultProb,
	   double				NoIntervals,
	   double				lbound,							
	   double				ubound,
	   double				mbound,
       int&                 lossunit,
       bool                 linearBeta);


