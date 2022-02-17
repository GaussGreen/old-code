
#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpclosedforms/stochasticvol_ln_formula.h"
//#include "gpinfra/gramfunctorStochVolcf.h"



CC_BEGIN_NAMESPACE( ARM )



///////////////////////////////////////////////////////////////////////
///  
///			Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(	double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
	int averagingtype= ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);
	
	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(a);
}

double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(	double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
	int averagingtype=ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);

	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(a);
}



double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(int i,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
	int averagingtype=ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);

	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(i,a);
}

double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(int i,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
	int averagingtype=ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);

	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(i,a);
}


double Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(int i,int j,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
int averagingtype=ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);

	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(i,j,a);
}

double Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(int i,int j,double forward,
							double strike,
							double maturity,
							double drift,
							double volatility,
							double voldrift,
							double volvol,
							double averaging,
							double reset,
							double CallPut,
							int nbhermite
							)
{
int averagingtype=ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST;
	ArgumentList a(forward,strike,maturity,drift,volatility,voldrift,volvol,averaging,reset,averagingtype,CallPut,nbhermite);

	Power_Expression<ARM_CF_StochasticVol_LN_Formula> y;
	return y(i,j,a);
}
///////////////////////////////////////////////////////////////////////
///  
///			End of Export Pricing Functions 
///
///////////////////////////////////////////////////////////////////////



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/