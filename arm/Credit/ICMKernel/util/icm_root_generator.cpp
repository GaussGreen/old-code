#include "firsttoinc.h"

#include "ICMKernel\util\icm_root_generator.h"




ICM_Root_Generator::ICM_Root_Generator(unsigned int itsSeed0, unsigned int itsNbfactors0) 
{
	Init() ;
	setParameters(itsSeed0, itsNbfactors0) ; 
} 



ICM_Root_Generator::ICM_Root_Generator(ICM_Root_Generator & Generator0) 
{
	Init() ;
	setParameters(Generator0.getSeed(), Generator0.getNbfactors()) ; 
}




void ICM_Root_Generator::setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0) 
{
	setSeed(itsSeed0) ;
	setNbfactors(itsNbfactors0) ;
}



void ICM_Root_Generator::setSeed(unsigned int itsSeed0) 
{
	itsSeed = itsSeed0 ;
}



unsigned int ICM_Root_Generator::getSeed() 
{
	return itsSeed ;
}



unsigned int ICM_Root_Generator::getNbfactors()  
{
	return itsNbfactors ;
}




