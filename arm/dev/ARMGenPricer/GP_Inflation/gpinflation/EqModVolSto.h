

/*----------------------------------------------------------------------------*/
 
/*! \file EqModVolSto.h
 * Copyright (c) CDC IXIS CM February 2007 Paris
 *
 *  \HW model on IRS Part... Heston model on each cpi.
 *
 *
 *	\author  Mathieu Bernardo
 *	\version 1.0
 *	\date February 2007
 */

/*----------------------------------------------------------------------------*/


#ifndef _INGPINFLATION_HW1FIRVSINF_H
#define _INGPINFLATION_HW1FIRVSINF_H


/// gpinflation
#include "gpinflation/infcurv.h"

/// gpmodels
#include "gpmodels/HWSV1F.h"
#include "gpmodels/EqModelParams.h"
#include "gpmodels/typedef.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/mktdatamanagerrep.h"

/// gpbase
#include "gpbase/countedptr.h"


/// gpclosedforms
#include "gpclosedforms/CIRBond.h"

/// std
#include <cmath>
#include <complex>
CC_USING_NS(std,complex)

CC_BEGIN_NAMESPACE( ARM )



/****************************************************

	\class		ARM_EQHWSV
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007
	
	\diffusion:

	-> asset:
	dz_t = (m_t - k.z_t)dt + p_t( sqrt(V_t).dW_t + sqrt(v_t).dw_t)
	z_0  = 0

	-> variance:
	dV_t = h_t.( 1  - V_t)dt + q_t sqrt(V_t). ( sqrt(1-r_t²) dW_t + r_t dZ_t)
	V_0  = 1

	-> reconstruction function:
	dm_t = p_t²(V_t+v) - 2k.m_tdt
	m_0	 = 0

	->paramters:
	k	double			MeanReversion
	q_t	vector<double>	VolMeanReversion
	h_t	vector<double>	VolOfVol
	v_t	vector<double>	CompoundVol
	r_t	vector<double>	Correlation
	p_t	vector<double>	Volatility
	
	\NB:
	The short rate is r_t= z_t + f_t
	where f_t is the instantaneous forward rate such that f_t=-d lnB(t,T) /dT

	The correlation is implictly defined under v. 
	This difference is relevant and necessary in the inflation context:


*****************************************************/


typedef ARM_CountedPtr< ARM_EQHWSV_ModelParams		>	ARM_EQHWSV_ModelParamsPtr;
typedef ARM_CountedPtr< ARM_EQHWSV_NumMethods		>	ARM_EQHWSV_NumMethodsPtr;

struct VarianceDistribution{

	VarianceDistribution(	ARM_CIRModelParams		& cirModel, 
							const double			& time, 
							const double			& x0);

	VarianceDistribution(	const ARM_GP_Vector		& weight, 
							const ARM_GP_Vector		& point, 
							const ARM_GP_Vector		& density);

	virtual~VarianceDistribution(){}

	ARM_GP_Vector		itsWeight;
	ARM_GP_Vector		itsPoint;
	ARM_GP_Vector		itsDensity;

};


class ARM_EQHWSV: public ARM_HWSV1F
{        

protected:

	ARM_EQHWSV_ModelParamsPtr	itsModelParams;
	ARM_EQHWSV_NumMethodsPtr	itsNumMethods;

public:
	ARM_EQHWSV( ){ }
	ARM_EQHWSV( ARM_ZeroCurve*, const ARM_EQHWSV_ModelParams*, const ARM_EQHWSV_NumMethods*);
	ARM_EQHWSV( const ARM_EQHWSV & );
	virtual ARM_Object*	Clone	( ) const	{	return new ARM_EQHWSV(*this); }
	virtual ~ARM_EQHWSV(){ }


	// stantard products

	virtual ARM_VectorPtr DiscountFactor( 	const string& curveName, 
											double evalTime, 
											double maturityTime, 
											const ARM_PricingStatesPtr& states) const{ 
		ARM_VectorPtr t;	return t; 
	}


	virtual ARM_VectorPtr VanillaSwaption(	const string& curveName,
											double evalTime,
											double swapResetTime,
											const ARM_GP_Vector& fixNotional,
											const ARM_GP_Vector& floatNotional,
											double floatStartTime,
											double floatEndTime,
											const ARM_GP_Vector& floatResetTimes,
											const ARM_GP_Vector& floatStartTimes,
											const ARM_GP_Vector& floatEndTimes,
											const ARM_GP_Vector& floatIntTerms,
											const ARM_GP_Vector& fixPayTimes,
											const ARM_GP_Vector& fixPayPeriods,
											const ARM_GP_Matrix& strikesPerState,
									        int callPut,
											const ARM_PricingStatesPtr& states,
											bool isConstantNotional = true,
											bool isConstantSpread = true,
											bool isConstantStrike = true) const;


	virtual ARM_VectorPtr VanillaCaplet(	const string& curveName, 
											double evalTime,
											double payTime, 
											double period,
  											double payNotional,
											double fwdResetTime, 
											double fwdStartTime,
 										    double fwdEndTime,
											double fwdPeriod,
											const ARM_GP_Vector& strikesPerState,
									        int capFloor,
											const ARM_PricingStatesPtr& states) const{  ARM_VectorPtr t;	return t; }

	virtual string ExportShortName() const { return "LHWSE";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	virtual bool	ValidateModelParams	( const ARM_ModelParams& params) const;



	void ReduceSwaptionSet	(	const double & evalTime,
								const double & floatEndTime,
								const double & fixPayTime,
								const ARM_GP_Vector& fixNotional,
								const ARM_GP_Vector& floatNotional,
								const ARM_GP_Matrix& strikesPerState,
								bool isConstantNotional,
								bool isConstantSpread,
								bool isConstantStrike ) const;


	static double	CptBasketPrice( const ARM_GP_Vector& polCoef,	const ARM_GP_Vector& expCoef,	const double & x);
	double			CptFrontiere(	const ARM_GP_Vector& polCoef,	const ARM_GP_Vector& expCoef,	double x ) const;

private:

//	VarianceDistribution	CptVarDist		( const double& tStart, const double& tEnd, const double & x0) const;
	ARM_CIRModelParams		GenerateCirModel( const double& time) const ;

	double CptLocalVariance(			const double		& tStart, 
										const double		& tEnd) const;

	void   CptExpPolCoef(				ARM_GP_Vector		& polCoef, 
										ARM_GP_Vector		& expCoef, 
										const ARM_GP_Vector & fixPayPeriods, 
										const ARM_GP_Vector & fixPayTimes, 
										const double		& floatStartTime,
										const double		& evalTime,
										const double		& swapResetTime,
										const double		& strike,
										const double		& variance) const;

	void CptFreezedQuantities(			const double		& evalTime,
										const double		& startTime, 
										const double		& endTime,
										const double		& strike,
										const ARM_GP_Vector	& fixPayTimes,
										const ARM_GP_Vector	& fixPayPeriods,
										double				& newStrike,
										double				& drift) const;

	double CptOswDriftFreezedNonVolSto(	const double		& evalTime,
										const double		& swapResetTime,
										const double		& floatStartTime, 
										const double		& floatEndTime,
										const double		& strike,
										const double		& notional,
										const int			& callPut,
										const ARM_GP_Vector & fixPayTimes,
										const ARM_GP_Vector	& fixPayPeriods )  const;

	double	CptOswCorrelNull(			const	double			& evalTime,
										const	double			& swapResetTime,
										const	double			& floatStartTime, 
										const	double			& floatEndTime,
										const	double			& strike,
										const	double			& notional,
										const	int				& callPut,
										const	ARM_GP_Vector	& fixPayTimes,
										const	ARM_GP_Vector	& fixPayPeriods) const;

};

						





/****************************************************

	\class		ARM_ModifiedAnalyticalRiccati ( ARM_ModifiedRiccati )
	\author		Mathieu Bernardo
	\version	1.0
	\date		February 2007
	
	
 *  AIM: 
 
	The aim of this class (struct) is to define the elementary green function of an affine model.
	The equation differs from ARM_AnalyticalRiccati for the alpha dynamic:
	We solve backard the 2 ode ( cf affine model) by backward propagation of dt amplitude.

 *	SDE:
	
	dX_t=	( e  - kx.y ) dt	+ Ex		 ( sqrt(y).DWy_t +  sqrt(S).DWx_t )
	dY_t=	( 1  - ky.y ) dt	+ Ey.sqrt(y) ( cos(a) .DWy_t +  sin(a).	DWz_t )

	all the brownian motions are independent on to each other

 * 	PDE:

	L	=	(e - kx.y).Dx		+	(1 - ky.y).Dy	
		+	Ex².(y + S)/2.Dxx²	+	Ey².y/2.Dyy²			+  EyEx cos(a).y Dxy

 *  PDE (fourier):

	l	=	-( S.u².Ex²/2 + iu.e) + y ( -u²Ex²/2 + iu kx)	+	(1 - (ky+iu.EyEx.cos(a)).y).Dy	+	Ey².y/2.Dyy²
			
 *  ODE:

	Esp ( Exp( Gamma(t,u) + y_t. Alpha(t,u) / F_{t-dt} ) = Exp( Gamma(t-dt,u) + y_{t-dt}.Alpha(t-dt,u) )
	d Gamma / dt	= -( S.u².Ex²/2 + iu.e ) + Alpha
	d Alpha / dt	=	Ey²/2.Alpha² - (ky+iu.EyEx.cos(a)).Alpha  + ( -u²Ex²/2 + iu kx)

 *  VARIABLES:

	double itsEta		-> e		Limit Value:		Asset		Process
	double itsKappa_x	-> kx		Mean Reversion:		Asset		Process
	double itsEspilon_x -> Ex		Scaling Parameter:	Asset		Process
	double itsVol_x		-> S		Normal Volatility:  Asset		Process

	double itsCorrel	-> cos(a)	correlation			Variance	Process		
	double itsKappa_y	-> ky		Mean Reversion:		Variance	Process
	double itsEspilon_y -> Ey		Volatility:			Variance	Process

	double itsGamma		-> Gamma	Gamma Function
	double itsAlpha		-> Alpha	Alpha Function

	double itsReV		-> -Ex²/2	parameter Re(v) = -Ex²/2.u²
	double itsImV		-> kx		parameter Im(v) = kx.u
	double itsReG		-> -S.Ex²/2	parameter Re(v) = -S.u².Ex²/2
	double itsImG		-> -e		parameter Im(v) = -e.u

 *  NB:
	
	We restrict to the computation on a unique variable complex u.

*****************************************************/
/*struct ARM_ModifiedRiccati : public ARM_AnalyticalRiccati
{        


	ARM_ModifiedRiccati()	{}
	ARM_ModifiedRiccati(	double eta,
							double kappa_x,
							double epsilon_x,
							double vol_x,
							double kappa_y,
							double epsilon_y,
							double correl,
							double dt );

	void SolveSystem( complex<double>& u, complex<double>& alpha, complex<double>& gamma) const;

	double itsReV;
	double itsImV;
	double itsReG;
	double itsImG;
	};

inline void ARM_ModifiedRiccati::SolveSystem( complex<double>& u, complex<double>& alpha, complex<double>& gamma) const{
	
		complex<double> U (-u.imag(),u.real());

		double ReV = itsReV*(u.real()*u.real() - u.imag()*u.imag() ) - itsImV*u.imag();
		double ImV = 2 * itsReV * u.real()*u.imag() + itsImV*u.real() ;
		complex<double> V (ReV,	ImV);
	
		ARM_AnalyticalRiccati::SolveSystem(U,V,alpha,gamma);

		double ReG = itsReG*(u.real()*u.real() - u.imag()*u.imag() ) - itsImG*u.imag();
		double ImG = 2 * itsReG * u.real()*u.imag() + itsImG*u.real() ;
		complex<double> G (ReG,	ImG);

		gamma += itsdt * G;
};

*/


CC_END_NAMESPACE()

#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/



