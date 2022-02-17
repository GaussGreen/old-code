/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file smile_nonparametric.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date	fev 2007
 */
#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"
#include "gpbase/removenagwarning.h"

#include <cmath>
#include <complex>

#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"


#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/nonparametric_spline.h"
#include "gpclosedforms/nonparametric_quantile.h"
#include "gpclosedforms/smile_nonparametric.h"

#include <glob/expt.h>   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000

double NonParametric_LogSmile::call_option(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike)
{
	return BlackSholes_Formula(S,
							 NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike),
							 1.0,
							 strike,
							 T,
							 K_CALL);
}

double NonParametric_LogSmile::digital_call_option(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
						   double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike)
{
	return DigitalBlackSholes_Formula(S,
							 NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike)*sqrt(T),
							 1.0,
							 strike,
							 K_CALL);
}

double NonParametric_LogSmile::inverse_distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
	return NonParametric_LN_Quantile_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,S,T,proba);
}


double NonParametric_LogSmile::gaussian_to_distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
								double index_begin,double index_end,int beginflag,int endflag,double S,double T,double k)
{
	return NonParametric_LN_Quantile_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,S,T,NormalCDF(k));
}

double NonParametric_LogSmile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{

	return NonParametric_LN_Quantile_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
		Underlying[1],Underlying[2],Underlying[3],Underlying[4],t,NormalCDF(x));
}

double NonParametric_LogSmile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
{
	return NormalCDF(NonParametric_LN_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x));
}

double NonParametric_LogSmile::quantile(const ArgumentList& Underlying, double x, double t)
{
	return NonParametric_LN_Quantile_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
		Underlying[1],Underlying[2],Underlying[3],Underlying[4],t,x);
}

double NonParametric_LogSmile::probability_density(const ArgumentList& Underlying, double x, double t)
{
	double d1=NonParametric_LN_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x*0.9995);
	double d2=NonParametric_LN_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x*1.0005);


	return (d2-d1)/(0.001*x);
}

double NonParametric_LogSmile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	return NonParametric_LN_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x);
}

double NonParametric_LogSmile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying,
												 double x, double t)
{
	return 0;
}

double NonParametric_LogSmile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, 
												 double x, double t)
{
	return 0;
}

ArgumentList* NonParametric_LogSmile::HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	return new ArgumentList(Underlying->V(0),Underlying->V(1),Underlying->V(2),
		(*Underlying)[0],(*Underlying)[1],(*Underlying)[2],(*Underlying)[3],positivenumber*(*Underlying)[4]);
}


///////////////////////////////////////////////////////////////////////////////////////






double NonParametric_NormalSmile::call_option(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
		double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike)
{
	return BlackSholes_Formula(S,
							 NonParametric_NormalVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike),
							 1.0,
							 strike,
							 T,
							 K_CALL);
}

double NonParametric_NormalSmile::digital_call_option(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
						   double index_begin,double index_end,int beginflag,int endflag,double S,double T,double strike)
{
	return VanillaDigitalOption_N(S,
        NonParametric_NormalVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike),
							 strike,
							 T,
							 K_CALL);
}

double NonParametric_NormalSmile::inverse_distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
	return NonParametric_N_Quantile_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,S,T,proba);
}


double NonParametric_NormalSmile::gaussian_to_distribution(ARM_GP_Vector* x,ARM_GP_Vector* y,ARM_GP_Vector* y2,
								double index_begin,double index_end,int beginflag,int endflag,double S,double T,double k)
{
	return NonParametric_N_Quantile_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,S,T,NormalCDF(k));
}

double NonParametric_NormalSmile::gaussian_to_distribution(const ArgumentList& Underlying, double x, double t)
{

	return NonParametric_N_Quantile_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
		Underlying[1],Underlying[2],Underlying[3],Underlying[4],t,NormalCDF(x));
}

double NonParametric_NormalSmile::distribution_to_gaussian(const ArgumentList& Underlying, double x, double t)
{
	return NormalCDF(NonParametric_N_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x));
}

double NonParametric_NormalSmile::quantile(const ArgumentList& Underlying, double x, double t)
{
	return NonParametric_N_Quantile_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
		Underlying[1],Underlying[2],Underlying[3],Underlying[4],t,x);
}

double NonParametric_NormalSmile::probability_density(const ArgumentList& Underlying, double x, double t)
{
	double d1=NonParametric_N_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x*0.9995);
	double d2=NonParametric_N_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x*1.0005);


	return (d2-d1)/(0.001*x);
}

double NonParametric_NormalSmile::probability_distribution(const ArgumentList& Underlying, double x, double t)
{
	return NonParametric_N_Distribution_Interpolate(Underlying.V(0),Underlying.V(1),Underlying.V(2),Underlying[0],
							Underlying[1],Underlying[2],Underlying[3],Underlying[4],t, x);
}

double NonParametric_NormalSmile::distribution_to_gaussian_first_derivative(int i,const ArgumentList& Underlying,
												 double x, double t)
{
	return 0;
}

double NonParametric_NormalSmile::probability_distribution_First_Derivative(int i,const ArgumentList& Underlying, 
												 double x, double t)
{
	return 0;
}

ArgumentList* NonParametric_NormalSmile::HomotheticTransformation(const ArgumentList* Underlying, double positivenumber)
{
	int i;
ARM_GP_Vector*  Vp0=new ARM_GP_Vector((Underlying->V(0))->size());
for(i=0;i<(Underlying->V(0))->size()-1;i++) (*Vp0)[i]=(*(Underlying->V(0)))[i];

ARM_GP_Vector*  Vp1=new ARM_GP_Vector((Underlying->V(1))->size());
for(i=0;i<(Underlying->V(1))->size()-1;i++) (*Vp1)[i]=(*(Underlying->V(1)))[i]*positivenumber;

ARM_GP_Vector*  Vp2=new ARM_GP_Vector((Underlying->V(2))->size());
for(i=0;i<(Underlying->V(2))->size()-1;i++) (*Vp2)[i]=(*(Underlying->V(2)))[i]*positivenumber;

	return new ArgumentList(Vp0,Vp1,Vp2,
		(*Underlying)[0],(*Underlying)[1],(*Underlying)[2],(*Underlying)[3],positivenumber*(*Underlying)[4]);
}




CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
