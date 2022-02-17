#include "ICMKernel\\random\icm_randomInvNormMoro.h"
#include "ICMKernel\util\icm_macro.h"


ICM_RandomInvNormMoro::ICM_RandomInvNormMoro() 
{ 
	Init();
}

void ICM_RandomInvNormMoro::Init(){
	its_a1 = 2.50662823884;
	its_a2 = -18.61500062529;
	its_a3 = 41.39119773534;
	its_a4 = -25.44106049637;

	its_b1 = -8.4735109309;
	its_b2 = 23.08336743743;
	its_b3 = -21.06224101826;
	its_b4 = 3.13082909833;

	its_c1 = 0.337475482272615;
	its_c2 = 0.976169019091719;
	its_c3 = 0.160797971491821;
    its_c4 = 2.76438810333863E-02;
	its_c5 = 3.8405729373609E-03;
	its_c6 = 3.951896511919E-04;
    its_c7 = 3.21767881768E-05;
	its_c8 = 2.888167364E-07;
	its_c9 = 3.960315187E-07;
}

ICM_RandomInvNormMoro::ICM_RandomInvNormMoro(const ICM_RandomGenerator& pRandomUnif): ICM_RandomLaws(pRandomUnif) 
{
	Init();
}

ICM_RandomInvNormMoro::ICM_RandomInvNormMoro(const ICM_RandomInvNormMoro& ref):ICM_RandomLaws(ref)
{
	Init();
}

// randomGen is not cloned !! 
ARM_Object * ICM_RandomInvNormMoro::Clone(void)
{
	return new ICM_RandomInvNormMoro(*this);
}

ICM_RandomInvNormMoro& ICM_RandomInvNormMoro::operator=(const ICM_RandomInvNormMoro& ref)
{
	if (this!=&ref)
	{
		this->~ICM_RandomInvNormMoro(); 
		new(this)ICM_RandomInvNormMoro(ref);
	}
	return *this; 
}


ICM_RandomInvNormMoro::~ICM_RandomInvNormMoro() 
{
// nothing
}






