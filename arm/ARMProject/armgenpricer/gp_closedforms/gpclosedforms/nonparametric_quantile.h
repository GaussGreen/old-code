/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 02/15/2007
 *
 *  basic functions for the closed form framework 
 *
 *	\file nonparametric_quantile.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date fev 2007
 */
 
#ifndef _GP_CF_NONPARAMETRIC_QUANTILE_H
#define _GP_CF_NONPARAMETRIC_QUANTILE_H

#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpvector.h"


#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE(ARM)

void NonParametric_LogVolatility_TailSlope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative );

void NonParametric_NormalVolatility_TailSlope_Compute(std::vector<double>* x,std::vector<double>* y,double index_begin,double index_end,
							int beginflag,int endflag, double* beginderivative,double* endderivative );


double NonParametric_LogVolatility_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike);

double NonParametric_NormalVolatility_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double strike);
							
double NonParametric_LN_Distribution_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike);

double NonParametric_N_Distribution_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag, double S,double T, double strike);

double NonParametric_LN_Quantile_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba);

double NonParametric_N_Quantile_Interpolate(std::vector<double>* x,std::vector<double>* y,std::vector<double>* y2,
							double index_begin,double index_end,
							int beginflag,int endflag,double S,double T, double proba);





CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

