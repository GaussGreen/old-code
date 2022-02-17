
#include "ICMKernel\util\icm_Gen_Student1F.h"
#include "math.h"



ICM_Gen_Student1F::ICM_Gen_Student1F(unsigned int itsSeed0, unsigned int itsNbfactors0, unsigned int itsFreedomDegree0, doubleVector & itsBeta0) : 
			ICM_Gengauss1F(itsSeed0, itsNbfactors0, itsBeta0) 
{
	Init();
	itsFreedomDegree	=	itsFreedomDegree0;

} 



ICM_Gen_Student1F::ICM_Gen_Student1F(ICM_Gen_Student1F& Gen_Student0)
{
	Init();
	
	ICM_Gengauss1F::setParameters(Gen_Student0.getSeed(), Gen_Student0.getNbfactors(), Gen_Student0.getBeta());
	
	SetFreedomDegree(Gen_Student0.GetFreedomDegree());
}



void ICM_Gen_Student1F::setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, unsigned int itsFreedomDegree0, doubleVector & itsBeta0)
{
	ICM_Gengauss1F::setParameters(itsSeed0, itsNbfactors0, itsBeta0); 
	SetFreedomDegree(itsFreedomDegree0);
}


void ICM_Gen_Student1F::CommonFactors()
{
	double	N01;

	ICM_Gengauss1F::CommonFactors();
	
	itsInvSqrtChi2TimesFreedomDegree	=	0.0;
	double	Chi2	=	0.0;

	for (int i=0; i<itsFreedomDegree; i++)
	{
		N01	=	NAG_random_normal(0.,1.);
		Chi2	+=	N01 * N01;
	}

	itsInvSqrtChi2TimesFreedomDegree = sqrt((double)itsFreedomDegree/Chi2);

}

