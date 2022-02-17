#ifndef _ICM_Gengauss1F_FactorLoading_2Points_FACTORLOADING_2POINTS_H_
#define _ICM_Gengauss1F_FactorLoading_2Points_FACTORLOADING_2POINTS_H_

/*********************************************************************************/
/*! \class  ICM_Gengauss1F_FactorLoading_2Points ICM_Gengauss1F_FactorLoading_2Points.h "ICM_Gengauss1F_FactorLoading_2Points.h"
 *  \author L. JACQUEL
 *	\version 1.0
 *	\date   May 2005
 *	\file   ICM_Gengauss1F_FactorLoading_2Points.h
 *		\brief Creates a MC generator with a 2 Points Factor Loadings Model
/*****************************************************************************************/

#include "ICMKernel/util/ICM_Gengauss1F.h"


class ICM_Gengauss1F_FactorLoading_2Points : public ICM_Gengauss1F
{

private :
	
	vector<double> itsAlpha;
	vector<double> itsBeta_2;
	vector<double> itsCour_2;

protected: 

	void Init()
	{
		itsAlpha.clear();
		itsBeta_2.clear();
		itsCour_2.clear();
	}


public : 

	ICM_Gengauss1F_FactorLoading_2Points() {Init();}
	
	ICM_Gengauss1F_FactorLoading_2Points(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector& itsBeta01, doubleVector& itsBeta02, doubleVector& itsAlpha0);

	ICM_Gengauss1F_FactorLoading_2Points(ICM_Gengauss1F_FactorLoading_2Points& Gengauss1F0);

	void setParameters(unsigned int itsSeed0, unsigned int itsNbfactors0, doubleVector& itsBeta01, doubleVector& itsBeta02, doubleVector& itsAlpha);

	void setBeta_2(doubleVector& itsBeta02);
	void setAlpha(doubleVector& Alpha) {itsAlpha	=	Alpha;}

	doubleVector& getBeta_2() {return itsBeta_2;}
	doubleVector& getCour_2() {return itsCour_2;}
	doubleVector& getAlpha() {return itsAlpha;}

	double generateRandom(unsigned int j)
	{
		double TheCommonFactor	=	GetCommonFactor();
		doubleVector	TheBeta	=	getBeta();
		doubleVector	TheCour	=	getCour();

		if (TheCommonFactor < itsAlpha[j])
			return TheBeta[j] * TheCommonFactor + TheCour[j] * NAG_random_normal(0.,1.);
		else
			return itsBeta_2[j] * TheCommonFactor + itsCour_2[j] * NAG_random_normal(0.,1.);
	}

	double CumulativeDistribution(int i, double x);
	
	~ICM_Gengauss1F_FactorLoading_2Points() 
	{
		itsAlpha.clear();
		itsBeta_2.clear();
		itsCour_2.clear();
	}


	// ----------------------------
	//	Copy of members data
	// ----------------------------
	void BitwiseCopy(const ARM_Object* src)
	{
		ICM_Gengauss1F::BitwiseCopy(src);

		ICM_Gengauss1F_FactorLoading_2Points * Gengauss1F = (ICM_Gengauss1F_FactorLoading_2Points *) src;
		setBeta_2(Gengauss1F->getBeta_2());
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
		ICM_Gengauss1F_FactorLoading_2Points * theClone = new ICM_Gengauss1F_FactorLoading_2Points() ;

		theClone->Copy(this) ;
 
		return(theClone) ;
	}

};


#endif 