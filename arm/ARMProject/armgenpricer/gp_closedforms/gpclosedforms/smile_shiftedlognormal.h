#ifndef _GP_CF_SMILE_SHIFTEDLOGNORMAL_H
#define _GP_CF_SMILE_SHIFTEDLOGNORMAL_H

#include "firsttoinc.h"
#include "gpbase/port.h"

#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/vanilla_bs.h"


CC_BEGIN_NAMESPACE(ARM)

class SLNDensity_Approx : public ARM_GP::UnaryFunc<double,double> 
{
public: 
		SLNDensity_Approx(double fwd,double time,double strike1,double target1,double strike2,double target2) :
		itsFwd(fwd),
		itsTime(time),
		itsStrike1(strike1),
		itsTarget1(target1),
		itsStrike2(strike2),
		itsTarget2(target2)
		{
		};

		inline virtual double operator() (double shift) const
		{
			double price2	= BS(itsFwd + shift, itsStrike2 + shift, itsTime, GetVol(shift));
			return price2-itsTarget2;
		};
		inline double GetVol(double shift) const
		{
			bool success(true);
			double sigma = VanillaImpliedVol_BS(itsFwd + shift, itsStrike1 + shift, itsTime, itsTarget1, 1, NULL, &success);
			return sigma;
		};
private:
	double						itsFwd;
	double						itsTime;
	double						itsStrike1;
	double						itsTarget1;
	double						itsStrike2;
	double						itsTarget2;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the Shifted Loognormal Distributions the  
///		 Underlying[0],	//INDEX
///      Underlying[1],	//SIGMA
///		 Underlying[2],	//ALPHA
///		 
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct ShiftedLogNormal_Smile
{
	static double call_option(double f,double K,double tex,double sigma, double alpha);
	static double digital_call_option(double f,double K,double tex,double sigma,  double alpha);
	static double inverse_distribution(double f,double K,double tex, double sigma, double alpha);
	static double gaussian_to_distribution(double f,double K,double tex, double sigma, double alpha);
	static void   volshift_calibration( double f,  double tex, 
										double K1, double price1,
										double K2, double price2,
										/// result
										double& sigma,
										double& alpha);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return Underlying[1];}	
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

