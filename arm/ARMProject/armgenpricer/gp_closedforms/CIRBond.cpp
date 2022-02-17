/*
 * Copyright (c) IXIS CIB February 2007 Paris
 *
/*! \file EqModelParams.cpp
 *
 *  \brief 
 *
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date February 2007
 *
 *	\class ARM_CIRModelParams
 */

#include "gpbase/numericconstant.h"
#include "gpbase/functor.h"

//#include "gpmodels\typedef.h"

#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/CIRBond.h"
#include "gpclosedforms/gamma.h"


#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"

// c++ standard headers
#include <map>
#include  <stdlib.h>

#pragma warning(disable : 4786) 
using std::pair;


CC_BEGIN_NAMESPACE( ARM )



const double DEF_PRECISION	= 1e-5;
const double BOUND_PRECISION= 1e-6;
const double PADE_TOL		= 1e-7;
const double SHIFT			= 1e-6;
const double PADE_PERIOD	= 2*ARM_NumericConstants::ARM_PI;
const double PADE_FREQUENCY = 1.0;

const int	 PER_ANNI		= 365;

////////////////////////////////////////////////////////////////
///	Class  : UnaryFuncLaplaceInverse
///	Routine: operator()
///	Returns: U
///	Action : compute the inverse laplace
/// Algo   : "A flexible Inverse Laplace Transform Algorithm and its Application" J.Ahn, S.Kang, Y. Kwon 
////////////////////////////////////////////////////////////////

template<typename T,typename U>

U UnaryFuncLaplaceInverse<T,U>::operator () ( T x ) const{
	double			S;
	double			R;
	Complex			c(0,ARM_NumericConstants::ARM_PI/itsPeriod);

	ARM_GP_T_Vector<Complex > a(NB_PADE);
	ARM_GP_T_Vector<Complex > d(NB_PADE);
	ARM_GP_T_Vector<Complex > A(NB_PADE+1);	
	ARM_GP_T_Vector<Complex > B(NB_PADE+1);

	ARM_GP_T_Matrix<Complex > e(NB_PADE,NB_PADE);
	ARM_GP_T_Matrix<Complex > q(NB_PADE,NB_PADE);

	Complex	z		= exp(x*c);

	Complex tmp		= itsFrequency;
	Complex post	= itsFunction(tmp);
	Complex prev	= 0.5*post;

	a[0]	= prev ;
	d[0]	= prev ;
	A[0]	= 0.0 ;
	B[0]	= 1.0 ;
	A[1]	= d[0];
	B[1]	= 1.0;

	e(1,0)	= 0.0;
	tmp		= itsFrequency+c;
	post	= itsFunction(tmp);
	q(0,1)	= post/prev;
	prev	= post;

	d[1]	= -q(0,1);
	tmp		= z*d[1];
	A[2]	= A[1] + tmp*A[0];
	B[2]	= B[1] + tmp*B[0];

	e(2,0)	= 0.0;
	tmp		= itsFrequency+2.0*c;
	post	= itsFunction(tmp);
	q(1,1)	= post/prev;
	prev	= post;
	e(0,1)	= e(1,0)+q(1,1)-q(0,1);

	d[2]	= -e(0,1);
	tmp		= z*d[2];
	A[3]	= A[2] + tmp*A[1];
	B[3]	= B[2] + tmp*B[1];

	S = 1.0;
	for( int M=3; M<NB_PADE; M++){
	
		e(M,0)	= 0.0;
		tmp		= itsFrequency +((double) M)*c;
		post	= itsFunction(tmp);
		q(M-1,1)= post/prev;
		prev	= post;

		if( div(M+1,2).rem ==0 ) {
			int k=div(M+1,2).quot;
			for( int j=1; j<k;j++){
				e(2*k-2*j-1,j)	= e(2*k-2*j,j-1)+q(2*k-2*j,j)-q(2*k-2*j-1,j);
				q(2*k-2*j-2,j+1)= q(2*k-2*j-1,j)*e(2*k-2*j-1,j)/e(2*k-2*j-2,j);
			}
			d[M] = -q(0,k);
		}
		else {
			int k=div(M,2).quot;
			for( int j=1; j<k; j++){
				e(2*k-2*j,j)	= e(2*k-2*j+1, j-1) +q(2*k-2*j+1,j) - q(2*k-2*j,j);
				q(2*k-2*j-1,j+1)= q(2*k-2*j,j)*e(2*k-2*j,j)/e(2*k-2*j-1,j);
			}
			e(0,k)	= e(1,k-1)+q(1,k)-q(0,k);
			d[M]	= -e(0,k);

		}
		A[M+1]	= A[M] + z*d[M]*A[M-1];
		B[M+1]	= B[M] + z*d[M]*B[M-1];	

		R = norm(B[M+1]/A[M+1]-B[M]/A[M]);
		if( R< PADE_TOL && R>S ) break;
		S = R;
	}
	
	return exp(itsFrequency*x)*real(A[M-1]/B[M-1])/itsPeriod;
}





////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: ARM_CIRModelParams
///	Returns: void
///	Action : constructor
////////////////////////////////////////////////////////////////

ARM_CIRModelParams::ARM_CIRModelParams( const ARM_ModelParamVector & params ):ARM_ModelParams(params){
	int i;

    const std::vector<double>* schedAt	= GetModelParam( ARM_ModelParamType::Volatility		).ToCurveModelParam().GetCurve()->GetAbscisses();
 	if( schedAt.size() == 0)		ARMTHROW(ERR_INVALID_ARGUMENT, "Volatility Model Params should contaions at least one element" );
	const std::vector<double>* schedBt	= GetModelParam( ARM_ModelParamType::MeanReversion	).ToCurveModelParam().GetCurve()->GetAbscisses();
	if( schedBt.size() == 0)		ARMTHROW(ERR_INVALID_ARGUMENT, "MeanReversion Model Params should contaions at least one element" );
    const std::vector<double>* schedCt	= GetModelParam( ARM_ModelParamType::LongTermVol	).ToCurveModelParam().GetCurve()->GetAbscisses();
 	if( schedCt.size() == 0)		ARMTHROW(ERR_INVALID_ARGUMENT, "LongTermVol Model Params should contaions at least one element" );


	map<int,int> tmpMap;

	for( i=0; i< schedAt.size(); i++)
		tmpMap.insert( pair<int,int> (schedAt[i],0) );

	for( i=0; i< schedBt.size(); i++)
		tmpMap.insert( pair<int,int> (schedBt[i],0) );

	for( i=0; i< schedCt.size(); i++)
		tmpMap.insert( pair<int,int> (schedCt[i],0) );


	std::map<int,int>::iterator it;
	for(it= tmpMap.begin(); it!=tmpMap.end(); it++)
		itsSchedule.push_back( it->first );

	if(itsSchedule[0] > K_NEW_DOUBLE_TOL)
        itsSchedule.insert(itsSchedule.begin(),0.0);

	
	ARM_CurveModelParam	cMeanReversion	( GetModelParam( ARM_ModelParamType::MeanReversion	).ToCurveModelParam()	);
	ARM_CurveModelParam	cLongTermVol	( GetModelParam( ARM_ModelParamType::LongTermVol	).ToCurveModelParam()	);
	ARM_CurveModelParam	cVolatility		( GetModelParam( ARM_ModelParamType::Volatility		).ToCurveModelParam()	);

	itsMapParam.insert(pair<ARM_ModelParamType::ParamNb,ARM_CurveModelParam> ( ARM_ModelParamType::MeanReversion,	cMeanReversion) );
	itsMapParam.insert(pair<ARM_ModelParamType::ParamNb,ARM_CurveModelParam> ( ARM_ModelParamType::LongTermVol,		cLongTermVol) );
	itsMapParam.insert(pair<ARM_ModelParamType::ParamNb,ARM_CurveModelParam> ( ARM_ModelParamType::Volatility,		cVolatility) );

	itsLaplaceParam.isOptimized=false;
}

ARM_CIRModelParams::ARM_CIRModelParams( const ARM_CIRModelParams & rhs ):ARM_ModelParams(rhs){ 

	itsSchedule		= rhs.itsSchedule;
	itsLaplaceMeth	= rhs.itsLaplaceMeth;
	itsLaplaceParam	= rhs.itsLaplaceParam;
	itsMapParam		= rhs.itsMapParam;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptElem
///	Returns: void
///	Action : computeE[ exp(-rho - psi*XT - phi*Int(Xs,s=0..T) / Xt=x) -> Exp( -rho' -psi'*Xt - phi'*Int(Xs,s=0..t))	
//
////////////////////////////////////////////////////////////////

void ARM_CIRModelParams::CptElem(	const double	& time,
									const double	& kappa,
									const double	& theta,
									const double	& sigma,
									Complex			& phi,
									Complex			& psi,
									Complex			& rho){

	Complex lambda	= psi;
	Complex mu		= phi;
	Complex gamma	= sqrt( kappa*kappa+2*sigma*sigma*phi );
	Complex expG	= exp(gamma*time/( ( Complex) PER_ANNI) );
	Complex dem		= lambda*sigma*sigma*(expG-1.0)+gamma-kappa+expG*(gamma+kappa);
	Complex num		= lambda*(gamma+kappa+expG*(gamma-kappa))+ 2.0*mu*(expG-1.0);

	phi	=	2.0*gamma*exp(time/PER_ANNI*(gamma+kappa)/2.0);
	phi	/=	dem;
	phi	=	log(phi);
	phi *=	-2.0/(sigma*sigma);
	phi	*=	kappa*theta;
	
	rho	+=	phi;
	psi =	num/dem;
	phi =	mu;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptGlobal
///	Returns: complex
///	Action : computeE[ exp(-rho - psi*XT - phi*Int(Xs,s=0..T) / X0=x0) 	
//
////////////////////////////////////////////////////////////////

Complex ARM_CIRModelParams::CptGlobal(	const	std::vector<double>	& schedule,	
										const	mapParamType	& mapParam,
										const	Complex			& x0,
										const	Complex			& phi,
										const	Complex			& psi,
										const	Complex			& rho){
	double tEnd;
	double tStart;
	Complex c_phi(phi);
	Complex c_psi(psi);
	Complex c_rho(rho);

	tEnd=schedule[schedule.size()-1];
	for (int  i= 0; i<schedule.size()-1;i++){
		tStart = schedule[schedule.size()-2-i];

		CptElem(	tEnd-tStart,
					mapParam.find(ARM_ModelParamType::MeanReversion)->second.GetValue(tStart),
					mapParam.find(ARM_ModelParamType::LongTermVol  )->second.GetValue(tStart),
					mapParam.find(ARM_ModelParamType::Volatility   )->second.GetValue(tStart),
					c_phi,
					c_psi,
					c_rho);

		tEnd=tStart;
	}
	
	return exp(-c_rho-c_psi*x0);
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptLaplace_Bond
///	Returns: complex
///	Action : compute	E[ exp(- phi*Int(Xs,s=0..T) / X0=x)
////////////////////////////////////////////////////////////////

Complex	ARM_CIRModelParams::CptLaplace_Bond	(	const	std::vector<double>	& schedule,	
												const	mapParamType	& mapParam,
												const	Complex			& x0,
												const	Complex			& x){
	return CptGlobal( schedule, mapParam, x0, x, 0.0, 0.0);
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptExp_Bond
///	Returns: double
///	Action : compute	E[ exp(- phi*Int(Xs,s=0..T) / X0=x0)
////////////////////////////////////////////////////////////////

std::vector<double> ARM_CIRModelParams::GenerateSchedule( const double & time) const{
	std::vector<double> schedule;
	for( int i =0; i< itsSchedule.size(); i++){
		if( itsSchedule[i]< time)
			schedule.push_back(itsSchedule[i]);
		else break;
	}
	schedule.push_back(time);
	
	return schedule;
}

double ARM_CIRModelParams::CptExp_Bond ( const double & time, const double & x0, const double & phi) const {

	return real ( CptLaplace_Bond(	GenerateSchedule(time), itsMapParam, x0, phi) );
}




////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptDensity_Bond
///	Returns: double
///	Action : compute the density of Int(Xs,s=0..T) / X0=x)
////////////////////////////////////////////////////////////////





/****************************************************

	\class		ARM_CIR_InvLaplace
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007

	Use of the simplest method: Pade ( coninued fraction) 
	
*****************************************************/

struct ARM_CIR_InvLaplace:public ARM_GP::UnaryFunc<Complex,Complex>{

	ARM_CIR_InvLaplace(		const std::vector<double>	& schedule, 
							const mapParamType	& mapParam,
							const double		& x0):
								itsSchedule(schedule),
								itsMapParam(mapParam),
								itsX0(x0){}

	virtual ~ARM_CIR_InvLaplace(){}
	virtual Complex	operator()	( Complex  x )	const{
		return ARM_CIRModelParams::CptLaplace_Bond( itsSchedule, itsMapParam, itsX0, x );

	}

private:
	std::vector<double>	itsSchedule;
	mapParamType	itsMapParam;
	double			itsX0;
};



double ARM_CIRModelParams::CptDensity_Bond( const double & time, const double & x0, const double & x) const{

	ARM_CIR_InvLaplace laplace( GenerateSchedule(time), itsMapParam, x0);
	UnaryFuncLaplaceInverse<Complex,Complex> inverse(laplace,PADE_PERIOD, PADE_FREQUENCY);
	return real(inverse(x));
}


double ARM_CIRModelParams::CptDensity_Bond	(	const	std::vector<double>	& schedule,	
												const	mapParamType	& mapParam,
												const	double			& x0,
												const	double			& x){

	ARM_CIR_InvLaplace laplace( schedule, mapParam, x0);
	UnaryFuncLaplaceInverse<Complex,Complex> inverse(laplace,PADE_PERIOD, PADE_FREQUENCY);
	return real(inverse(x));
}
		
////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptDistribSet
///	Returns: DistribSet
///	Action : compute the expectation of the density of Int(Xs,s=0..T) / X0=x)
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////


class InvLaplace : public ARM_GP::UnaryFunc<double,double>{
public: 
	InvLaplace(		const std::vector<double>	&	schedule,	
					const mapParamType	&	mapParam, 
					const double		&	x0,
					const double		&	target):
					itsMapParam(mapParam),	itsSchedule(schedule),	itsX0(x0), itsTarget(target){ };
	
	virtual double operator() (double x) const{ 
	return ARM_CIRModelParams::CptDensity_Bond ( itsSchedule, itsMapParam, itsX0, x) -itsTarget;	}

private:
	std::vector<double>	itsSchedule;
	mapParamType	itsMapParam;
	double			itsX0;
	double			itsTarget;
};


mapParamType Scale( const double & factor, mapParamType map){

	mapParamType::iterator it;
	

	for ( it = map.begin(); it != map.end(); it++){
		if (  it->first== ARM_ModelParamType::LongTermVol )
			( * it->second.GetCurve() )*= factor;
		else if(  it->first== ARM_ModelParamType::Volatility )
			( * it->second.GetCurve() )*= sqrt(factor);
	}
	return map;
}

DistribSet	ARM_CIRModelParams::CptDistribSet(const double & time,	const double & x0, const int& nbDisc, const double & factor ) { 
	itsLaplaceParam.isOptimized		= false;
	std::vector<double>		schedule	= GenerateSchedule(time);
	mapParamType		mapParam	= Scale( factor, itsMapParam);
	const double		epsilon		= SHIFT;

	double tmp, tmpX0;

	tmpX0 = x0 *factor;

	double exp = 1.0-real(ARM_CIRModelParams::CptLaplace_Bond	(schedule, mapParam, tmpX0, epsilon));
	exp	/= epsilon;
	
	double target =  DEF_PRECISION * CptDensity_Bond ( schedule, mapParam, tmpX0, exp);
	double bound  = exp;

	InvLaplace func( schedule, mapParam, tmpX0, target );
	UnaryFuncWithNumDerivative<double> funcWithDeriv(func);
	T_NewtonRaphsonSolverBoundedLenght<UnaryFuncWithNumDerivative<double> > solver(funcWithDeriv,BOUND_PRECISION, BOUND_PRECISION,BOUND_PRECISION);
	solver.SetLength(1.0);
	solver.setInitialGuess(bound);
	bound= solver.Solve();

	double infBound, supBound;

	bool isInf = false;
	if ( bound> exp) {
		supBound = bound;
		infBound = 2*exp-bound>0.0 ? 2*exp-bound:BOUND_PRECISION;
		bound	 = infBound;
		isInf	 = true;

	}
	else {
		infBound = bound;
		supBound = 2*exp-bound;
		bound	 = supBound;

	}
	solver.setInitialGuess(bound);
	bound= solver.Solve();

	if ( isInf )
		infBound= bound>0?bound:0.0;
	else
		supBound = bound;

	GaussLegendre_Coefficients dist( nbDisc, infBound, supBound);

	itsDistribSet.itsWeight		= std::vector<double>(nbDisc);
	itsDistribSet.itsPoint		= std::vector<double>(nbDisc);
	itsDistribSet.itsDensity	= std::vector<double>(nbDisc);

	for (int i = 0; i<  nbDisc; i++){
		tmp = dist.get_point(i);
		itsDistribSet.itsPoint[i]  = tmp/factor;
		itsDistribSet.itsWeight[i] = dist.get_weight(i)/factor;
		itsDistribSet.itsDensity[i]= CptDensity_Bond ( schedule, mapParam, tmpX0, tmp)*factor;
	}

	return itsDistribSet; 
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptWeight_Bond
///	Returns: double
///	Action : compute the expectation of the density of Int(Xs,s=0..T) / X0=x)
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ARM_CIRModelParams::CptWeight_Bond() const{

	int nbDisc = itsDistribSet.itsWeight.size();
	if ( nbDisc<=0 ) 
		ARMTHROW(ERR_INVALID_ARGUMENT, "The Distribution has not be built" );
	double sum = 0.0;
	for ( int i = 0; i<  nbDisc; i++)
		sum	+=itsDistribSet.itsWeight[i]* itsDistribSet.itsDensity[i];
	return sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptExpectation_Bond
///	Returns: double
///	Action : compute the expectation of Int(Xs,s=0..T) / X0=x) by a small shift
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////

double	ARM_CIRModelParams::CptExpectation_Bond	() const{

	double sum = 0.0;
	int nbDisc = itsDistribSet.itsWeight.size();
	if ( nbDisc<=0 ) 
		ARMTHROW(ERR_INVALID_ARGUMENT, "The Distribution has not be built" );

	for ( int i = 0; i<  nbDisc; i++)
		sum		+=itsDistribSet.itsWeight[i]* itsDistribSet.itsDensity[i]*itsDistribSet.itsPoint[i];

	return sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_CIRModelParams
///	Routine: CptVariance_Bond
///	Returns: double
///	Action : compute the variance of Int(Xs,s=0..T) / X0=x)
			
////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
double	ARM_CIRModelParams::CptVariance_Bond () const{

	double	tmp;
	int		nbDisc	= itsDistribSet.itsWeight.size();
	double	sum		= 0.0;
	double	exp		= CptExpectation_Bond();

	if ( nbDisc<=0 ) 
		ARMTHROW(ERR_INVALID_ARGUMENT, "The Distribution has not be built" );


	for ( int i = 0; i<  nbDisc; i++){
		tmp= itsDistribSet.itsPoint[i] - exp ;
		sum		+=itsDistribSet.itsWeight[i]* itsDistribSet.itsDensity[i]*tmp*tmp;
	}
	return sum;
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

