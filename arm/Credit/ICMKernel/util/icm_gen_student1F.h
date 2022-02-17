#ifndef _ICM_Gen_Student1F_H_
#define _ICM_Gen_Student1F_H_

/*********************************************************************************/
/*! \class  ICM_Gen_Student1F ICM_Gen_Student1F.h "ICM_Gen_Student1F.h"
 *  \author L Jacquel
 *	\version 1.0
 *	\date   May 2005
 *	\file   ICM_Gen_Student1F.h
 *		\brief Creates a MC generator with one factor student conditional independance
/*****************************************************************************************/

#include "ICMKernel/util/ICM_Gengauss1F.h"


class ICM_Gen_Student1F : public ICM_Gengauss1F
{

public:
	// SOME PARAMETERS
	unsigned int	GetFreedomDegree() {return itsFreedomDegree;}
	void	SetFreedomDegree(unsigned int data) {itsFreedomDegree = data;}

private :
	
	unsigned int		itsFreedomDegree;
	double				itsInvSqrtChi2TimesFreedomDegree;

protected: 

	void Init()
	{
		itsFreedomDegree	=	4;
		itsInvSqrtChi2TimesFreedomDegree	=	0.0;
	}


public : 

	ICM_Gen_Student1F() {Init();}
	
	ICM_Gen_Student1F(unsigned int itsSeed0, unsigned int itsNbfactors0, unsigned int itsFreedomDegree0, doubleVector & itsBeta0) ;

	ICM_Gen_Student1F(ICM_Gen_Student1F & Gen_Student1F0) ;

	void setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, unsigned int itsFreedomDegree0, doubleVector & itsBeta0) ;

	virtual void CommonFactors();

	virtual void Reset() {NAG_random_init_repeatable(getSeed());};


	double	GetChi2Factor() {return itsInvSqrtChi2TimesFreedomDegree;}

	double generateRandom(unsigned int j) 
	{
		doubleVector	TheBeta, TheCour;
		double	TheCommonFactor;

		TheBeta	=	getBeta();
		TheCour	=	getCour();
		TheCommonFactor	=	GetCommonFactor();

		return (TheBeta[j] * TheCommonFactor + TheCour[j] * NAG_random_normal(0.,1.)) * itsInvSqrtChi2TimesFreedomDegree;
	}

	~ICM_Gen_Student1F() 
	{
	}


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object * src)
	{
	    ICM_Gen_Student1F * Gen_Student1F = (ICM_Gen_Student1F *) src ;

		itsFreedomDegree	=	Gen_Student1F->itsFreedomDegree;
		itsInvSqrtChi2TimesFreedomDegree	=	Gen_Student1F->itsInvSqrtChi2TimesFreedomDegree;
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
		ICM_Gen_Student1F * theClone = new ICM_Gen_Student1F() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

};


#endif 