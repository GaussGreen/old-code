#ifndef _ICM_GENGAUSS1F_H_
#define _ICM_GENGAUSS1F_H_

/*********************************************************************************/
/*! \class  ICM_gengauss1F ICM_gengauss1F.h "ICM_gengauss1F.h"
 *  \author R Marciano
 *	\version 1.0
 *	\date   March 2004
 *	\file   ICM_gengauss1F.h
 *		\brief Creates a MC generator with one factor gaussian conditional independance
/*****************************************************************************************/

#include "ICMKernel/util/icm_root_generator.h"
#include "ICMKernel/glob/icm_maths.h"



class ICM_Gengauss1F : public ICM_Root_Generator
{

private :
	
	double	itsCommonFactor;
	double	its_CommonFactor_drift;

	vector<double> itsBeta;
	vector<double> itsCour;
	
	vector<double>	itsCommonFactors;

	double	itsIdiocsyncraticValue;

protected: 

	void Init()
	{
		itsCommonFactor = 0.;
		its_CommonFactor_drift	=	0.0;
		itsIdiocsyncraticValue	=	0.0;

		itsBeta.clear();
		itsCour.clear();

		itsCommonFactors.clear();
	}


public : 

	ICM_Gengauss1F() { Init() ;}
	
	ICM_Gengauss1F(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector & itsBeta0) ;

	ICM_Gengauss1F(ICM_Gengauss1F & Gengauss1F0) ;

	void setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector & itsBeta0) ;

	void setBeta(doubleVector & itsBeta0) ;
	void setCommonFactorsBeta(doubleVector & itsCommonFactors0);

	virtual void setSeed(unsigned int seed0);
	void setCommonFactor_Drift(double data) {its_CommonFactor_drift = data;}

	virtual void CommonFactors() ;

	virtual void Reset() {NAG_random_init_repeatable(getSeed());};

	doubleVector& getBeta() ;
	doubleVector& getCour() {return itsCour;}
	doubleVector& getCommonFactors() {return itsCommonFactors;}

	double	GetCommonFactor() {return itsCommonFactor;}
	double	GetCommonFactor_Drift() {return its_CommonFactor_drift;}
	double	GetIdiosyncraticValue() {return itsIdiocsyncraticValue;}

	void	setCommonFactors(doubleVector & itsCommonFactors0);

	double generateRandom(unsigned int j) 
	{
		itsIdiocsyncraticValue	=	NAG_random_normal(0.,1.);
		
		return itsBeta[j] * itsCommonFactor + itsCour[j] * itsIdiocsyncraticValue;
	}

	double getConditionalDefaultProbability(double level, unsigned int j);

	~ICM_Gengauss1F() 
	{
	itsBeta.clear();
	itsCour.clear();
	itsCommonFactors.clear();
	}


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_Gengauss1F * Gengauss1F = (ICM_Gengauss1F *) src ;

		setSeed(Gengauss1F->getSeed()) ;
		setNbfactors(Gengauss1F->getNbfactors()) ;
		setBeta(Gengauss1F->getBeta()) ;
		setCommonFactor_Drift(Gengauss1F->GetCommonFactor_Drift());
		setCommonFactors(Gengauss1F->getCommonFactors()) ;
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
		ICM_Gengauss1F * theClone = new ICM_Gengauss1F() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

};


#endif 