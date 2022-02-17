
#include "ICMKernel\util\icm_gengauss_1F_FactorLoading_2Points.h"
#include "ICMKernel\glob\icm_maths.h" 

#include "math.h"

#include <nags.h>		//	s15abc
#include <nagg01.h>		//	g01hac


ICM_Gengauss1F_FactorLoading_2Points::ICM_Gengauss1F_FactorLoading_2Points(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector& itsBeta01, doubleVector& itsBeta02, doubleVector& itsAlpha0) :
			ICM_Gengauss1F(itsSeed0, itsNbfactors0, itsBeta01) 
{
	Init() ;
	setBeta_2(itsBeta02);
	setAlpha(itsAlpha0);
} 



ICM_Gengauss1F_FactorLoading_2Points::ICM_Gengauss1F_FactorLoading_2Points(ICM_Gengauss1F_FactorLoading_2Points & Gengauss1F0)
{
	Init() ;
	ICM_Gengauss1F::setParameters(Gengauss1F0.getSeed(), Gengauss1F0.getNbfactors(), Gengauss1F0.getBeta()) ; 
	
	setBeta_2(Gengauss1F0.getBeta_2());
	setAlpha(Gengauss1F0.getAlpha());
}


void ICM_Gengauss1F_FactorLoading_2Points::setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector& itsBeta01, doubleVector& itsBeta02, doubleVector& itsAlpha)
{
	ICM_Gengauss1F::setParameters(itsSeed0, itsNbfactors0, itsBeta01) ; 
	setBeta_2(itsBeta02) ;
	setAlpha(itsAlpha) ;
}



void ICM_Gengauss1F_FactorLoading_2Points::setBeta_2(doubleVector & itsBeta02) 
{
	itsBeta_2.resize(itsBeta02.size()) ;
	itsBeta_2 = itsBeta02 ;
	itsCour_2.resize(itsBeta02.size()) ;
	unsigned int nFactors = getNbfactors() ;
	for(unsigned int j = 0 ; j < nFactors ; j++)
		itsCour_2[j] = sqrt(1 - itsBeta02[j]*itsBeta02[j]) ;
}


double ICM_Gengauss1F_FactorLoading_2Points::CumulativeDistribution(int	i, double x)
{
	double	alpha;
	double	beta_1;
	double	beta_2;

	doubleVector	TheBeta0	=	getBeta();
	// From Andersen
	// Should I check the index...
	alpha	=	itsAlpha[i];
	beta_1	=	TheBeta0[i];
	beta_2	=	itsBeta_2[i];

	// g01hac(x, y, rho, NAGERR_DEFAULT)
	return	NAG_cumul_normal(x) + NAG_bivariate_normal_dist(alpha, x, beta_1) - NAG_bivariate_normal_dist(alpha, x, beta_2);

}

