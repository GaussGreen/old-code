#ifndef _ICM_ROOT_GENERATOR_H_
#define _ICM_ROOT_GENERATOR_H_

/*********************************************************************************/
/*! \class  ICM_Generator icm_root_generator.h "icm_root_generator.h"
 *  \author R Marciano
 *	\version 1.0
 *	\date   March 2004
 *	\file   ICM_root_generator.h
 *		\brief Creates un generator MC 
/***********************************************************************************/

#include <vector>


#include "ARMKernel/glob/armglob.h"
#include "ICMKernel\util\icm_utils.h"
#include "ICMKernel\util\icm_macro.h"

typedef std::vector<double> doubleVector ;

class ICM_Root_Generator : public ARM_Object
{

private :
	
	unsigned int itsSeed ;

	unsigned int itsNbfactors ; 
	

protected :

	void Init()
	{
		itsSeed = 0 ;
		itsNbfactors = 0 ;
	}

	void setNbfactors(unsigned int nFactors0)
	{
		itsNbfactors = nFactors0 ;
	}


public : 

	ICM_Root_Generator() { Init() ;}
	
	ICM_Root_Generator(unsigned int itsSeed0, unsigned int itsNbfactors0) ;

	ICM_Root_Generator(ICM_Root_Generator & Generator0) ;

	void setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0) ;

	virtual void setSeed(unsigned int seed0) ;

	unsigned int getSeed() ;

	unsigned int getNbfactors() ;

	virtual double generateRandom(unsigned int j)
	{
		return 0. ;
	}

	virtual double getConditionalDefaultProbability(double level, unsigned int j)
	{
		return 0. ;
	}

	virtual void Reset() {};
	
	virtual void CommonFactors(){}
	

	virtual doubleVector& getBeta() {ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
	virtual doubleVector& getCour() {ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }
	virtual doubleVector& getCommonFactors() {ICMTHROW(ERR_INVALID_OPERATION,"Not Implemented"); }

	virtual double	GetCommonFactor() {return 0.0;}
	virtual void setCommonFactor_Drift(double data){}
	virtual double	GetIdiosyncraticValue() {return 0.0;}

	virtual double CumulativeDistribution(int i, double x) {return 0.0;}

	~ICM_Root_Generator() 
	{}


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_Root_Generator * Generator = (ICM_Root_Generator *) src ;

		itsSeed      = Generator->getSeed() ;
		itsNbfactors = Generator->getNbfactors() ;
	}


	// -------------
	//	Copy Method 
	// -------------
	
	void Copy(const ARM_Object* src)
	{
		ARM_Object::Copy(src) ;
 
		BitwiseCopy(src) ;
	}

	// --------------
	//	Clone Method
	// --------------

	ARM_Object * Clone(void)
	{
		ICM_Root_Generator * theClone = new ICM_Root_Generator() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

};


#endif 