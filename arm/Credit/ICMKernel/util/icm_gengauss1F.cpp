
#include "ICMKernel/glob/icm_maths.h" 
#include "ICMKernel\util\icm_gengauss1F.h"
#include "math.h"

#include <nags.h>		//	s15abc


ICM_Gengauss1F::ICM_Gengauss1F(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector & itsBeta0) : 
			ICM_Root_Generator(itsSeed0, itsNbfactors0) 
{
	Init() ;
	setBeta(itsBeta0);
	setSeed(itsSeed0);
} 



ICM_Gengauss1F::ICM_Gengauss1F(ICM_Gengauss1F & Gengauss1F0)
{
	Init() ;
	ICM_Root_Generator::setParameters(Gengauss1F0.getSeed(), Gengauss1F0.getNbfactors()) ; 
	setBeta(Gengauss1F0.getBeta()) ;
}

void ICM_Gengauss1F::setSeed(unsigned int seed0)
{
	ICM_Root_Generator::setSeed(seed0);
	NAG_random_init_repeatable(seed0);
}


void ICM_Gengauss1F::setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector & itsBeta0)
{
	ICM_Root_Generator::setParameters(itsSeed0, itsNbfactors0) ; 
	setBeta(itsBeta0) ;
}



void ICM_Gengauss1F::setBeta(doubleVector & itsBeta0) 
{
	itsBeta.resize(itsBeta0.size()) ;
	itsBeta = itsBeta0 ;
	itsCour.resize(itsBeta0.size()) ;
	unsigned int nFactors = getNbfactors() ;
	for(unsigned int j = 0 ; j < nFactors ; j++)
		itsCour[j] = sqrt(1 - itsBeta0[j]*itsBeta0[j]) ;
	
	itsCommonFactors.resize(1);
}


void ICM_Gengauss1F::setCommonFactors(doubleVector & itsCommonFactors0) 
{
	itsCommonFactors.resize(itsCommonFactors0.size()) ;
	itsCommonFactors = itsCommonFactors0 ;
}

void ICM_Gengauss1F::CommonFactors() 
{
	itsCommonFactor = NAG_random_normal(its_CommonFactor_drift, 1.);
	itsCommonFactors[0]	=	itsCommonFactor;
}



doubleVector & ICM_Gengauss1F::getBeta()  
{
	return itsBeta ;
}


double ICM_Gengauss1F::getConditionalDefaultProbability(double level, unsigned int j)
{
	double	TheValue;

	TheValue	=	level - itsBeta[j] * itsCommonFactor;
	TheValue	/=	itsCour[j];

	return	NAG_cumul_normal(TheValue);
}

