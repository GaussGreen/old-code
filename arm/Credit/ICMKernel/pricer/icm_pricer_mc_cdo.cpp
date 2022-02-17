#include "ICMKernel\pricer\icm_pricer_mc_cdo.h"

void __stdcall Register(long address)
{
	globalFunc = (FUNCPTR)address;
}


double __stdcall CallFunc()
{
	return globalFunc();
}

double ICM_Pricer_MC_CDO::ComputePrice(qCMPMETH mode)
{
	double value=0.;

	for (int i=0;i<itsnbsimuls;i++)
	{
		value+=PayOff();
	}

	value/=itsnbsimuls;

	return value;

}

double ICM_Pricer_MC_CDO::PayOff()
{
	double value = globalFunc();
	
	return value;
}