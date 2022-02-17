/* Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file nonparametric_quantile.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date	fev 2007
 */
#include "firsttoinc.h"
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
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/nonparametric_spline.h"
#include "gpclosedforms/nonparametric_quantile.h"
#include "gpclosedforms/spreadoption_nonparametric_formula.h"

#include "expt.h"   // for the exceptions

using namespace std;

CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_EPS 1.0e-13
#define ARM_CF_MAXIT 2000





/// toutes ces fonctions suppose le vecteur des derivées secondes y2 precalculées

double NonParametric_LogVolatility_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike)
{
	double initialstrike=*(x->begin());
	
	if (strike< initialstrike)
	{
		double initialvol=(*(y->begin()));
		double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
		switch (beginflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
			{
				return initialvol*pow(fabs(strike/initialstrike),index_begin);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
			{
				return initialvol;
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
			{
				return initialvol+initialvolDer*initialstrike*log(fabs(initialstrike))/index_begin*
					(pow(fabs(log(fabs(strike))/log(fabs(initialstrike))),index_begin)-1.0);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
			{
				double val= initialvol+initialvolDer*(strike-initialstrike);
				if(val>0) return val; else return initialvol/10000.;
				break;
			}
			
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"NonParametric_LogVolatility_Interpolate::Begin Tail Method  : bad value ");
				break;
			}
		}
	}
	double finalstrike=*(x->end()-1);
	if (strike> finalstrike)
	{
		double finalvol=(*(y->end()-1));
		double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
		switch (endflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
			{
				return finalvol*pow(fabs(strike/finalstrike),index_end);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
			{
				return finalvol;
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
			{
				return finalvol+finalvolDer*finalstrike*log(fabs(finalstrike))/index_end*
					(pow(fabs(log(fabs(strike))/log(fabs(finalstrike))),index_end)-1.0);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
			{
				double val= finalvol+finalvolDer*(strike-finalstrike);
				if(val>0) return val; else return finalvol/10000.;
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_LogVolatility_Interpolate::end Tail Method  : bad value ");
				break;
			}
		}
	}

	return cubicspline_interpolate(*x,*y,*y2,strike);
}

/// the following function is necessary to transmit the slope of the smile et the boundary to the spline algo.
void NonParametric_LogVolatility_TailSlope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	double initialstrike=(*(x->begin()));
	double initialvol=(*(y->begin()));
	double finalstrike=(*(x->end()-1));
	double finalvol=(*(y->end()-1));
	switch (beginflag)
	{
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
		{
			*beginderivative=(index_begin/initialstrike)*initialvol;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
		{
			*beginderivative=0;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
		{
			double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
			*beginderivative=initialvolDer;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
		{
			double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
			*beginderivative=initialvolDer;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH :
		{
			double strike1=*(x->begin()+1);
			double strike2=*(x->begin()+2);
			double initialvolDer=((*(y->begin()+1))-initialvol)/(strike1-initialstrike);
			double initialvolDer2=((*(y->begin()+2))-(*(y->begin()+1)))/(strike2-strike1);
			*beginderivative=initialvolDer-(initialvolDer2-initialvolDer)/(strike2-initialstrike)*(strike1-initialstrike);
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_LogVolatility_Interpolate::Begin Tail Method  : bad value ");
			break;
		}
	}
	switch (endflag)
	{
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
		{
			*endderivative=(index_end/finalstrike)*finalvol;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
		{
			*endderivative=0;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
		{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;

		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
			{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;

			}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH :
			{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;

			}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_LogVolatility_Interpolate::end Tail Method  : bad value ");
			break;
		}
	}
}

/// pour les nouveaux types mixtes qui extrapolent au niveau de la distribution,
/// le traitement est dans cette fonction
double NonParametric_LN_Distribution_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike)
{
	double initialstrike=*(x->begin());
	if (strike<= initialstrike)
	{
		switch (beginflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH :
			{
				double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
					ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
				shift=(shift<initialstrike/10.)?shift:(initialstrike/10.);
				double calculstrike=initialstrike+1.001*shift;
				double vol=NonParametric_LogVolatility_Interpolate(x,y,y2,3,3,beginflag,endflag,calculstrike);  /// on selectionne 3, mais cela n'a pas d'importance
				double shiftedvol=NonParametric_LogVolatility_Interpolate(x,y,y2,3,3,beginflag,endflag,calculstrike+shift);
				double derive=(BlackSholes_Formula(S,shiftedvol,1.0,calculstrike+shift,T,K_CALL)-BlackSholes_Formula(S,vol,1.0,calculstrike,T,K_CALL))/shift;
				double initialDistrib= 1.0+derive;
				return pow(strike/initialstrike,index_begin)*initialDistrib;
				
				break;
			}
			
		default :
			{	
				double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
					ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
				shift=(shift<strike/10.)?shift:(strike/10.);
				double vol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike);
				double shiftedvol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike+shift);
				double derive=(BlackSholes_Formula(S,shiftedvol,1.0,strike+shift,T,K_CALL)-BlackSholes_Formula(S,vol,1.0,strike,T,K_CALL))/shift;
				return  1.0+derive;
				break;
			}
		}
	}
	double finalstrike=*(x->end()-1);
	if (strike>= finalstrike)
	{
		switch (endflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH :
			{
				double finalvol=(*(y->end()-1));
				double finalvol1=(*(y->end()-2));
				double finalstrike1=*(x->end()-2);
				double finalDistrib=1.0+(BlackSholes_Formula(S,finalvol1,1.0,finalstrike1,T,K_CALL)-BlackSholes_Formula(S,finalvol,1.0,finalstrike,T,K_CALL))/(finalstrike1-finalstrike);
				return 1.0-(1.0-finalDistrib)*exp((finalstrike-strike)*index_end/finalstrike);
				break;
			}
			
		default :
			{	
				double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
					ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
				shift=(shift<strike/10.)?shift:(strike/10.);
				double vol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike);
				double shiftedvol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike+shift);
				double derive=(BlackSholes_Formula(S,shiftedvol,1.0,strike+shift,T,K_CALL)-BlackSholes_Formula(S,vol,1.0,strike,T,K_CALL))/shift;
				return  1.0+derive;
				break;
			}
		}
	} 
	else
	{
		switch (beginflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::DISTRIBUTIONSMOOTH :
			{
				double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
					ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
				shift=(shift<strike/10.)?shift:(strike/10.);
				double vol=NonParametric_LogVolatility_Interpolate(x,y,y2,3,3,beginflag,endflag,strike);  /// on selectionne 3, mais cela n'a pas d'importance
				double shiftedvol=NonParametric_LogVolatility_Interpolate(x,y,y2,3,3,beginflag,endflag,strike+shift);
				double derive=(BlackSholes_Formula(S,shiftedvol,1.0,strike+shift,T,K_CALL)-BlackSholes_Formula(S,vol,1.0,strike,T,K_CALL))/shift;
				return  1.0+derive;
				break;
			}
			
		default :
			{	
				double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
					ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
				shift=(shift<strike/10.)?shift:(strike/10.);
				double vol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike-shift);
				double shiftedvol=NonParametric_LogVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike+shift);
				double derive=(BlackSholes_Formula(S,shiftedvol,1.0,strike+shift,T,K_CALL)-BlackSholes_Formula(S,vol,1.0,strike-shift,T,K_CALL))/(2.0*shift);
				return  1.0+derive;
				break;
			}
		}
	}
	
	
	
}


double NonParametric_LN_Quantile_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		std::vector<double>* x_0;
		std::vector<double>* y_0;
		std::vector<double>* y2_0;
		double index_begin_0;
		double index_end_0;
		double S_0;
		double T_0;
		int beginflag_0;
		int endflag_0;
		DistributionToInverse(std::vector<double>* x_1,std::vector<double>* y_1,std::vector<double>* y2_1,
			double index_begin_1,double index_end_1,int beginflag_1,int endflag_1,double S_1,double T_1):
		x_0(x_1),y_0(y_1),y2_0(y2_1),index_begin_0(index_begin_1),index_end_0(index_end_1),
			beginflag_0(beginflag_1),endflag_0(endflag_1),S_0(S_1),T_0(T_1)
		{}
		virtual double operator() (double K0)  const
		{
			return NonParametric_LN_Distribution_Interpolate(x_0,y_0,y2_0,index_begin_0,index_end_0,
				beginflag_0,endflag_0,S_0,T_0,K0);
		}
	};
	
	DistributionToInverse d(x,y,y2,index_begin,index_end,beginflag,endflag,S,T);
	
	
	if(proba<0.0005)			/// cut off down of the proba
	{	
		/// initial guess
		double guess = S/2;
		/// initial Standard dev of the gess
		double stddev=S/4;
		double eps=0.0005;
		double p0=Inverse(d,Inverse::ALWAYSPOSITIVE)(eps,guess,stddev,1e-12);  
		return p0*proba/eps;
	}
	
	else if (proba>0.9995)		/// cut off up of the proba
	{	
		/// initial guess
		double guess = 2*S;
		/// initial Standard dev of the gess
		double stddev=S;
		double eps=0.9995;
		double p0=Inverse(d,Inverse::ALWAYSPOSITIVE)(eps,guess,stddev,1e-12);  
		return p0*proba/eps;
	}
	else
	{
		/// initial guess
		double guess = S*proba*1.5;
		/// initial Standard dev of the gess
		double stddev=guess/2;
		return Inverse(d,Inverse::ALWAYSPOSITIVE)(proba,guess,stddev,1e-12);  // suppose that  Y is always postive
	}
}
	


///////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///          NORMAL VOLATILITY>
///
///
////////////////////////////////////////////////////////////////////////////////////////////////////////////


/// toutes ces fonctions suppose le vecteur des derivées secondes y2 precalculées

double NonParametric_NormalVolatility_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike)
{
	double initialstrike=*(x->begin());
	if (strike< initialstrike)
	{
		double initialvol=(*(y->begin()));
		double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
		switch (beginflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
			{
				return initialvol*pow((1.0+(strike-initialstrike)/index_begin*initialvolDer/initialvol),index_begin);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
			{
				return initialvol;
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
			{
				return initialvol+initialvolDer*initialstrike*log(fabs(initialstrike))/index_begin*
					(pow(fabs(log(fabs(strike))/log(fabs(initialstrike))),index_begin)-1.0);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
			{
				double val= initialvol+initialvolDer*(strike-initialstrike);
				if(val>0) return val; else return initialvol/10000.;
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
					"NonParametric_NormalVolatility_Interpolate::Begin Tail Method  : bad value ");
				break;
			}
		}
	}
	double finalstrike=*(x->end()-1);
	if (strike> finalstrike)
	{
		double finalvol=(*(y->end()-1));
		double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
	
		switch (endflag)
		{
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
			{
				return finalvol*pow((1.0+(strike-finalstrike)/index_end*finalvolDer/finalvol),index_end);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
			{
				return finalvol;
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
			{
				return finalvol+finalvolDer*finalstrike*log(fabs(finalstrike))/index_end*
					(pow(fabs(log(fabs(strike))/log(fabs(finalstrike))),index_end)-1.0);
				break;
			}
		case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
			{
				double val= finalvol+finalvolDer*(strike-finalstrike);
				if(val>0) return val; else return finalvol/10000.;
				break;
			}
		default :
			{	
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_NormalVolatility_Interpolate::end Tail Method  : bad value ");
				break;
			}
		}
	}

	return cubicspline_interpolate(*x,*y,*y2,strike);
}

void NonParametric_NormalVolatility_TailSlope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative )
{
	double initialstrike=(*(x->begin()));
	double initialvol=(*(y->begin()));
	double finalstrike=(*(x->end()-1));
	double finalvol=(*(y->end()-1));
	switch (beginflag)
	{
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
		{
			double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
			*beginderivative=initialvolDer;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
		{
			*beginderivative=0;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
		{
			double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
			*beginderivative=initialvolDer;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
		{
			double initialvolDer=((*(y->begin()+1))-initialvol)/((*(x->begin()+1))-initialstrike);
			*beginderivative=initialvolDer;
			break;
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_NormalVolatility_TailSlope_Compute::Begin Tail Method  : bad value ");
			break;
		}
	}
	switch (endflag)
	{
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYTAILINDEX :
		{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYCONSTANT :
		{
			*endderivative=0;
			break;
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYSLOGTAILINDEX :
		{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;
			
		}
	case ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::VOLATILITYLINEAR :
		{
			double finalvolDer=((*(y->end()-2))-finalvol)/((*(x->end()-2))-finalstrike);
			*endderivative=finalvolDer;
			break;
			
		}
	default :
		{	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"NonParametric_NormalVolatility_TailSlope_Compute::end Tail Method  : bad value ");
			break;
		}
	}
}

double NonParametric_N_Distribution_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike)
{
	double shift=ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::specific_shift(
	ARM_CF_NonParametric_Gaussian_PowerSpreadOption_Formula::K1);
	double vol=NonParametric_NormalVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike);
	double shiftedvol=NonParametric_NormalVolatility_Interpolate(x,y,y2,index_begin,index_end,beginflag,endflag,strike+shift);
	double derive=(VanillaOption_N(S,shiftedvol,strike+shift,T,K_CALL)-VanillaOption_N(S,vol,strike,T,K_CALL))/shift;
	return  1.0+derive;
}

double NonParametric_N_Quantile_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba)
{
class DistributionToInverse : public DoubleToDoubleFunc 
	{
	public: 
		std::vector<double>* x_0;
		std::vector<double>* y_0;
		std::vector<double>* y2_0;
		double index_begin_0;
		double index_end_0;
		double S_0;
		double T_0;
		int beginflag_0;
		int endflag_0;
		DistributionToInverse(std::vector<double>* x_1,std::vector<double>* y_1,std::vector<double>* y2_1,
			double index_begin_1,double index_end_1,int beginflag_1,int endflag_1,double S_1,double T_1):
		x_0(x_1),y_0(y_1),y2_0(y2_1),index_begin_0(index_begin_1),index_end_0(index_end_1),
			beginflag_0(beginflag_1),endflag_0(endflag_1),S_0(S_1),T_0(T_1)
		{}
		virtual double operator() (double K0)  const
		{
			return NonParametric_N_Distribution_Interpolate(x_0,y_0,y2_0,index_begin_0,index_end_0,
				beginflag_0,endflag_0,S_0,T_0,K0);
		}
	};
	
	DistributionToInverse d(x,y,y2,index_begin,index_end,beginflag,endflag,S,T);
	
	int sizex=x->size();
	if(proba<0.0001)			/// cut off down of the proba
	{	
		/// initial guess
		double guess = (*x)[0];
		/// initial Standard dev of the gess
		double stddev=((*x)[sizex-1]-guess)/10.0;
		double eps=0.0001;
		double p0=Inverse(d,Inverse::REAL)(eps,guess,stddev,1e-12);  
		return p0*proba/eps;
	}
	
	else if (proba>0.9999)		/// cut off up of the proba
	{	
		/// initial guess
		double guess = (*x)[sizex-1];
		/// initial Standard dev of the gess
		double stddev=fabs(((*x)[sizex-2]-guess));
		guess+=stddev;
		double eps=0.9999;
		double p0=Inverse(d,Inverse::REAL)(eps,guess,stddev,1e-12);  
		return p0*proba/eps;
	}
	else
	{
		/// initial guess
		double guess = ((*x)[sizex-1]+(*x)[0])/2.0;
		/// initial Standard dev of the gess
		double stddev=((*x)[sizex-1]-(*x)[0])/10.0;
		return Inverse(d,Inverse::REAL)(proba,guess,stddev,1e-12);  // suppose that  Y is always postive
	}
}

CC_END_NAMESPACE()


#undef ARM_CF_EPS
#undef ARM_CF_MAXIT

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
