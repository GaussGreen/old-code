
#include "ICMKernel\mod\icm_vasicek_diffusion.h"


double ICM_VasicekDiffusion::FwdZeroCoupon(double t, double T, double Rt)
{
	double flat_volatility =its_Volatility.Elt(0);

	double _m = 1.-exp(-its_MeanReversion*(T-t));
	double Mt = its_LongRate*(T-t)+(Rt-its_LongRate)*(_m)/its_MeanReversion;
	double Vt = -flat_volatility*flat_volatility/(2.*its_MeanReversion*its_MeanReversion*its_MeanReversion);
	Vt *= _m*_m;
	Vt += flat_volatility*flat_volatility/(its_MeanReversion*its_MeanReversion)*(T-t-_m/its_MeanReversion);

	double Esperance = exp(-Mt+0.5*Vt);

	return (Esperance);
}