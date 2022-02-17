/*
 * Copyright (c) IXIS CIB February 2007 Paris
 *
/*! \file CIRBOND.h
 *
 *  \brief 
 *
 *	\author  M. Bernardo
 *	\version 1.0
 *	\date February 2007
 *
 *	\Compute the price of Exp( -lambda*XT(x) -mu*int(Xs(x),s=t..T) / Ft) 
 *	and Xt(x) follows cir process:
 */


#ifndef _INGPCLOSEDFORMS_CIRMODELPARAMS_H
#define _INGPCLOSEDFORMS_CIRMODELPARAMS_H

#pragma warning(disable : 4786) 

/// gpbase
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/typedef.h"

#include <math.h>
#include <cmath>
#include <complex>

const int	 NB_PADE		= 100;

CC_BEGIN_NAMESPACE( ARM )
using std::complex;

typedef complex<double> Complex;
/****************************************************

	\class		UnaryFuncLaplaceInverse
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007

  	Use of the simplest method: Pade ( coninued fraction) 

*****************************************************/


template<typename T=Complex,typename U=T>
struct UnaryFuncLaplaceInverse{

	UnaryFuncLaplaceInverse( const ARM_GP::UnaryFunc<T,U>& f,
							 const double & period,
							 const double & frequency):itsFunction(f),itsPeriod(period),itsFrequency(frequency){	}

	virtual ~UnaryFuncLaplaceInverse(){};
	virtual U operator () ( T ) const ;

private:
	const	ARM_GP::UnaryFunc<T,U> & itsFunction;
	double	itsPeriod;
	double	itsFrequency;
};




typedef map<ARM_ModelParamType::ParamNb,ARM_CurveModelParam> mapParamType;

struct ARM_InvLaplace{

	ARM_InvLaplace( const ARM_GP_Vector	& schedule, 
					const mapParamType	& mapParam,
					const double		& x0):itsSchedule(schedule), itsMapParam(mapParam),itsX0(x0){}
	virtual ~ARM_InvLaplace(){	};
	
	virtual double	operator()	( const double & )	const	=0;

protected:

	double			itsX0;
	ARM_GP_Vector	itsSchedule;
	mapParamType	itsMapParam;
};



/****************************************************

	\class		ARM_CIRModelParams
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007

	dXt = kappa*(theta-Xt)dt - sigma*sqrt(Xt) dWt

	theta: LongTermVol
	kappa: MeanReversion
	sigma: Volatility
	
	\aim: compute Exp( -XT(x)*nu		 - mu *Int(Xs(x),s=t..T) / Xt(x)=x) 
				= Exp( -x*Psi(nu,mu,t,T) - kappa*theta*Phi(nu,mu,t,T) )

 
	
	\ Asset		Xt(x)
	\ Bond	ln(Xs(x),s=0..T)

		
*****************************************************/


struct LaplaceParam{

	LaplaceParam( ){ };
	LaplaceParam( const int & n){	itsParam.resize(n); };
	LaplaceParam( const ARM_GP_Vector & param):itsParam(param){ };
	LaplaceParam( const LaplaceParam& rhs):isOptimized(rhs.isOptimized),itsParam(rhs.itsParam){}
	virtual ~LaplaceParam(){};

	bool			isOptimized;
	ARM_GP_Vector	itsParam;
};

struct DistribSet{

	DistribSet(){}
	DistribSet( const double & weight, const double & point): itsWeight(weight), itsPoint(point){};
	virtual ~DistribSet(){};

	ARM_GP_Vector	itsWeight;
	ARM_GP_Vector	itsPoint;
	ARM_GP_Vector	itsDensity;
};


class ARM_CIRModelParams: public ARM_ModelParams{

public:
	ARM_GP_Vector itsSchedule;

	ARM_CIRModelParams( const ARM_ModelParamVector & params = ARM_ModelParamVector() );
	ARM_CIRModelParams(	const ARM_CIRModelParams & rhs );
	virtual ARM_Object*		Clone()	   const {	return new ARM_CIRModelParams(*this);	}

	ASSIGN_OPERATOR		( ARM_CIRModelParams );

	virtual string ExportShortName() const { return "LCIRB";}
	virtual ~ARM_CIRModelParams(){};

	// Calibration
    virtual void	PreProcessing	(ARM_ModelFitter& modelFitter,		int factorNb=0) {};
    virtual void	PostProcessing	(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,	int factorNb=0 ) {};
	
	ARM_GP_Vector	GenerateSchedule		(	const double & time) const;
	
	double			CptExpectation_Bond		() const;
	double			CptVariance_Bond		() const;
	double			CptWeight_Bond			() const;

	double			CptExp_Bond				(	const double & time,	const double & x0, const double & phi)	const;
	double			CptDensity_Bond			(	const double & time,	const double & x0, const double & x)	const;

	DistribSet		CptDistribSet			(	const double & time,	const double & x0, const int & nbDiscr, const double& factor = 1.0) ;

	static Complex	CptLaplace_Bond	(	const	ARM_GP_Vector	& schedule,	
										const	mapParamType	& mapParam,
										const	Complex			& x0,
										const	Complex			& phi);
	
	static double	CptDensity_Bond	(	const	ARM_GP_Vector	& schedule,	
										const	mapParamType	& mapParam,
										const	double			& x0,
										const	double			& x);


private:
	
	string			itsLaplaceMeth;
	LaplaceParam	itsLaplaceParam;
	mapParamType	itsMapParam;
	DistribSet		itsDistribSet;

	static void		CptElem			(	const double			& time,
										const double			& kappa,
										const double			& theta,
										const double			& sigma,
										Complex					& phi,
										Complex					& psi,
										Complex					& rho);

	static Complex	CptGlobal		(	const	ARM_GP_Vector	& schedule,	
										const	mapParamType	& mapParam,
										const	Complex			& x0,
										const	Complex			& phi,
										const	Complex			& psi,
										const	Complex			& rho);

};





CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/


