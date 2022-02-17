#include "ICMKernel\\random\ICM_RandomInvNormAcklam.h"
#include "ICMKernel\util\icm_macro.h"


ICM_RandomInvNormAcklam::ICM_RandomInvNormAcklam() 
{ 
	Init();
}

void ICM_RandomInvNormAcklam::Init(){

		its_a1 = -3.969683028665376e+01;
		its_a2 =  2.209460984245205e+02;
		its_a3 = -2.759285104469687e+02;
		its_a4 =  1.383577518672690e+02;
		its_a5 = -3.066479806614716e+01;
		its_a6 =  2.506628277459239e+00;
    
		its_b1 = -5.447609879822406e+01;
		its_b2 =  1.615858368580409e+02;
		its_b3 = -1.556989798598866e+02;
		its_b4 =  6.680131188771972e+01;
		its_b5 = -1.328068155288572e+01;
    
		its_c1 = -7.784894002430293e-03;
		its_c2 = -3.223964580411365e-01;
		its_c3 = -2.400758277161838e+00;
		its_c4 = -2.549732539343734e+00;
		its_c5 =  4.374664141464968e+00;
		its_c6 =  2.938163982698783e+00;
    
		its_d1 =  7.784695709041462e-03;
		its_d2 =  3.224671290700398e-01;
		its_d3 =  2.445134137142996e+00;
		its_d4 =  3.754408661907416e+00;
    
	  // Limits of the approximation region.
		its_u_low   = 0.02425;
		its_u_high  = 1.0 - its_u_low;

}
ICM_RandomInvNormAcklam::ICM_RandomInvNormAcklam(const ICM_RandomGenerator& pRandomUnif): ICM_RandomLaws(pRandomUnif)
{
	Init();
}

ICM_RandomInvNormAcklam::ICM_RandomInvNormAcklam(const ICM_RandomInvNormAcklam& ref):ICM_RandomLaws(ref)
{
	Init();
}

// randomGen is not cloned !! 
ARM_Object * ICM_RandomInvNormAcklam::Clone(void)
{
	return new ICM_RandomInvNormAcklam(*this);
}

ICM_RandomInvNormAcklam& ICM_RandomInvNormAcklam::operator=(const ICM_RandomInvNormAcklam& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomInvNormAcklam(); 
		new(this)ICM_RandomInvNormAcklam(ref);
	}
	return *this; 
}


ICM_RandomInvNormAcklam::~ICM_RandomInvNormAcklam() 
{
// nothing
}


