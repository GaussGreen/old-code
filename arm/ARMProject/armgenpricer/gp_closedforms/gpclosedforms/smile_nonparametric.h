/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_nonparametric.h
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date fev 2007
 */
 
#ifndef _GP_CF_SMILE_NONPARAMETRIC_H
#define _GP_CF_SMILE_NONPARAMETRIC_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpclosedforms/basic_distributions.h"

#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
/// here is a smile class that uses the following conventions :
/// 
///				For the for the non parametric Distributions the  
///      Underlying.V(0),	//Vector of x
///		 Underlying.V(1),	//Vector of y
///		 Underlying.V(2),	//Vector of y2
///		 Underlying[0],	//Index in the Beginning
///      Underlying[1],	//Index in the End
///		 Underlying[2]	//BEGINFLAG  can be : GAUSSIAN_DISTRIBUTIUON, LOGNORMAL_DISTRIBUTION,
///										GAUSSIAN_VOL,LOGNORMALVOL,POWER_DISTRIBUTION,POWER_VOL, 
/// 	 Underlying[3]	//ENDFLAG  can be : GAUSSIAN_DISTRIBUTIUON, LOGNORMAL_DISTRIBUTION,
///										GAUSSIAN_VOL,LOGNORMALVOL,POWER_DISTRIBUTION,POWER_VOL, 
///      Underlying[4]	//S 
///
///
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct NonParametric_LogSmile
{
	static double call_option(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,
		double S,double T,double strike);

	static double digital_call_option(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike);

	static double inverse_distribution(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
	double index_begin,double index_end,int beginflag,int endflag,double S,double T, double proba);

	static  double gaussian_to_distribution(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double k);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static ArgumentList* HomotheticTransformation2(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
	std::vector<double> Y2P;
};


struct NonParametric_NormalSmile
{
	static double call_option(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,
		double S,double T,double strike);

	static double digital_call_option(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike);

	static double inverse_distribution(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
	double index_begin,double index_end,int beginflag,int endflag,double S,double T, double proba);

	static  double gaussian_to_distribution(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double k);

	/// now come the generic interface to copula:
	static double gaussian_to_distribution(const ArgumentList& Underlying1, double x, double t);
	static double distribution_to_gaussian(const ArgumentList& Underlying1, double x, double t);
	static double quantile(const ArgumentList& Underlying1, double x, double t);
	static double probability_density(const ArgumentList& Underlying, double x, double t);
	static double probability_distribution(const ArgumentList& Underlying, double x, double t);
	static double distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying, double x, double t);
	static double probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, double x, double t);
	static ArgumentList* HomotheticTransformation(const ArgumentList* arg, double positivenumber);
	static ArgumentList* HomotheticTransformation2(const ArgumentList* arg, double positivenumber);
	static int BoundedfromBelow(void) {return 1;}
	static double LowerBoundary(const ArgumentList& Underlying) {return 0.;}
	std::vector<double> Y2P;
};





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

