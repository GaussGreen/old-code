/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HWSV.h
 *
 *  \brief 
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date October 2006
 */


#ifndef _INGPMODELS_HWSV_H
#define _INGPMODELS_HWSV_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

#include "gpinfra/pricingmodelir.h"

#include <vector>
CC_USING_NS(std,vector)

#include <cmath>
#include <complex>
CC_USING_NS(std,complex)
CC_USING_NS(std,sqrt)
CC_USING_NS(std,log)
CC_USING_NS(std,real)
CC_USING_NS(std,imag)
CC_USING_NS(std,exp)

#include "gpnumlib/odefunctions.h"

#include "gpclosedforms/gaussian_integrals.h"

#include <glob/armdef.h>

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
CC_BEGIN_NAMESPACE( ARM )

/// Constants for oscillatory integral computation
const double DEFAULT_PRECISION					= 1.e-5;
const double DEFAULT_PRECISION_SO				= 1.e-5;
const size_t OSCILLATION_NBMAX					= 10000;

/// Default value for closed form formula
const size_t DEFAULT_NBPTS_FIRST_OSCIL_STD		= 64;
const size_t DEFAULT_NBPTS_FIRST_OSCIL_LIMIT_STD= 32;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_STD		= 32;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_LIMIT_STD	= 16;
const double DEFAULT_OSCILLATION_LIMIT_STD		= 250.0;

const size_t DEFAULT_NBPTS_FIRST_OSCIL_STD_SO		= 64;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_STD_SO		= 32;
const size_t DEFAULT_NBPTS_LAST_OSCIL_STD_SO		= 8;
const size_t DEFAULT_NBPTS_FIRST_OSCIL_LIMIT_STD_SO	= 32;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_LIMIT_STD_SO	= 16;
const size_t DEFAULT_NBPTS_LAST_OSCIL_LIMIT_STD_SO	= 8;
const double DEFAULT_OSCILLATION_LIMIT_STD_SO		= 250.0;
const double DEFAULT_OSCILLATION_NEXT_LIMIT_STD_SO	= 750.0;
const double DEFAULT_OSCILLATION_LAST_LIMIT_STD_SO	= 2000.0;

const size_t DEFAULT_NBPTS_FIRST_OSCIL_ENH		= 32;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_ENH		= 16;
const size_t DEFAULT_NBPTS_LAST_OSCIL_ENH		= 8;
const double DEFAULT_OSCILLATION_LIMIT_ENH		= 3.0;
const double DEFAULT_OSCILLATION_NEXT_LIMIT_ENH	= 15.0;
const double DEFAULT_OSCILLATION_LAST_LIMIT_ENH	= 125.0;

const size_t DEFAULT_NBPTS_FIRST_OSCIL_ENH_SO		= 16;
const size_t DEFAULT_NBPTS_NEXT_OSCIL_ENH_SO		= 16;
const size_t DEFAULT_NBPTS_LAST_OSCIL_ENH_SO		= 8;
const double DEFAULT_OSCILLATION_LIMIT_ENH_SO		= 0.3;
const double DEFAULT_OSCILLATION_NEXT_LIMIT_ENH_SO	= 4.0;
const double DEFAULT_OSCILLATION_LAST_LIMIT_ENH_SO	= 15.0;
const double DEFAULT_INVFOURIER_POINT_ENH_SO		= 0.15;

/// Constant for interpolation threshold
const double INTERPOL_EPS = 1.e-6;

/// Default constants for Riccati's Runge-Kutta solver
const double	RK5_ODE_PRECISION		= 1.0e-5;
const double	RK5_TINY				= K_NEW_DOUBLE_TOL;

const double	RK4_NBSTEPS_PER_YEAR	= 2;
const int		RK4_NBSTEPS_MIN			= 3;
const int		RK4_NBSTEPS_MAX			= 1000;
const double	RK4_REF_KAPPA			= 0.04;	 
const double	RK4_FACTOR_KAPPA		= 1000.0;	 
const double	RK4_FACTOR_TIME			= 0.015;	 
const double	RK4_FACTOR_NU			= -1/(0.5*0.5*0.5);

//-----------------------------------------------------------------------------
// \struc ARM_AnalyticalRiccati
// \brief
//  Structure to deal with stepwise constant Riccati equation.
//  For each interval ]Ti,Ti+1] and a given complex numbers u, v and w
//  we will solve the system :
//
//		d(alpha)/dt = A.alpha^2 + (B1+u.B2).alpha + C1.w + C2.v
//		d(gamma)/dt = D.alpha
//
//	Systems are solved backward form TN=T to T0=t with initial conditions :
//
//		alpha(T,u,v)=0
//		gamma(T,u,v)=0
//-----------------------------------------------------------------------------
struct ARM_AnalyticalRiccati
{
	double itsdt;	// Ti+1-Ti

	double itsNu2;	// to take advantage of alpha squared term (=-0.5*Nu^2)

	double itsA;	// alpha squared term
	double itsB1;	// alpha u independent linear term
	double itsB2;	// alpha u dependent linear term
	double itsC1;	// alpha w independent constant term
	double itsC2;	// alpha v independent constant term

	double itsD;	// gamma cross linear term


	ARM_AnalyticalRiccati() : itsC1(0.0) {}

	inline void SolveSystem(complex<double>& u, complex<double>& v,
						complex<double>& alpha, complex<double>& gamma) const
	{
		complex<double> zero(0,0);
		SolveSystem(u,v,zero,alpha, gamma);
	}

	inline void SolveSystem(complex<double>& u, complex<double>& v, complex<double>& w,
							complex<double>& alpha, complex<double>& gamma) const
	{
		register double a2 = 2*itsA;
		complex<double> B(itsB1 + itsB2*u.real(),itsB2*u.imag());
		complex<double> C(itsC2*v.real()+(itsC1 ? itsC1*w.real() : 0.0),
						  itsC2*v.imag()+(itsC1 ? itsC1*w.imag() : 0.0));
		complex<double> unit(1,0);

		complex<double> delta = std::sqrt(B*B-2*a2*C);
		complex<double> alpha1((-B.real()+delta.real())/a2,(-B.imag()+delta.imag())/a2);
		complex<double> alpha2((-B.real()-delta.real())/a2,(-B.imag()-delta.imag())/a2);

		double adt = itsA*itsdt;
		double rhoTerm	= exp((alpha2.real()-alpha1.real())*adt); /// because itsdt = Ti+1-Ti
		double thetaTerm = (alpha2.imag()-alpha1.imag())*adt;
		complex<double> expTerm(rhoTerm*cos(thetaTerm),rhoTerm*sin(thetaTerm));
		complex<double> ratioTerm = (alpha-alpha1)/(alpha-alpha2);
		expTerm *= ratioTerm;

		alpha = (alpha1 - alpha2*expTerm)/(unit-expTerm);

		gamma += itsD*(-alpha1*itsdt - std::log((unit-expTerm)/(unit-ratioTerm))/itsA);
	}
};


//-----------------------------------------------------------------------------
// \class ARM_RiccatiHWSV
// \brief
//  general class to define Riccati system of H&W SV models
//-----------------------------------------------------------------------------
class ARM_RiccatiHWSV : public ARM_ODEFunc
{
protected :
	CC_IS_MUTABLE size_t itsSystemIdx;
	ARM_IntVector itsNbSteps;
	CC_IS_MUTABLE size_t itsNbDerivCalls;
	size_t itsSolverType;
	bool itsIsStdFormula;

public:

	ARM_RiccatiHWSV(size_t solverType=ARM_ODEFunc::RK5Adaptative,bool isStdFormula=false);
	ARM_RiccatiHWSV(const ARM_RiccatiHWSV& rhs);
	virtual ~ARM_RiccatiHWSV() {}

    virtual ARM_Object* Clone() const = 0;

	virtual size_t GetSystemSize() const = 0;
	virtual double GetYf(size_t i) const = 0;
	virtual size_t GetNbPhis() const = 0;

	inline void LocateSystemIdx(double t) const;

	void SetSystemIdx(size_t idx) {itsSystemIdx = idx;}
	int GetNbSteps(size_t i) const {return itsNbSteps[i];}
	ARM_IntVector GetNbSteps() const {return itsNbSteps;}
	void SetNbSteps(const ARM_IntVector& nbSteps) {itsNbSteps=nbSteps;}
	void SetNbSteps(size_t i, int nbSteps) {itsNbSteps[i]=nbSteps;}

	size_t GetNbDerivCalls() const {return itsNbDerivCalls;}
	void ResetNbDerivCalls() {itsNbDerivCalls=0;}

	virtual int ComputeNbSteps(const std::vector<double>& solverParams);

	virtual void derivs(double t, std::vector<double>& yt, std::vector<double>& dyt) const = 0;
};

inline void ARM_RiccatiHWSV::LocateSystemIdx(double t) const
{
/****/
	// New version
	if(itsSystemIdx+1<GetSystemSize() && t > GetYf(0) + K_NEW_DOUBLE_TOL)
		++itsSystemIdx;
	while(itsSystemIdx>0)
	{
		if(t > GetYf(itsSystemIdx-1) + K_NEW_DOUBLE_TOL)
		{
			while(itsSystemIdx+1<GetSystemSize() &&
				  t > GetYf(itsSystemIdx) + K_NEW_DOUBLE_TOL)
				++itsSystemIdx;
			break;
		}
		else
			--itsSystemIdx;
	}
/****/

/****
	// Old Version
	itsSystemIdx += (itsSystemIdx+1 < GetSystemSize() ? 1 : 0);
	while(itsSystemIdx>0)
	{
		if(t > GetYf(itsSystemIdx-1) + K_NEW_DOUBLE_TOL)
			break;
		else
			--itsSystemIdx;
	}
****/
}

class GaussLegendre_Coefficients;

//-----------------------------------------------------------------------------
// \class HWSVNumericals
// \brief
//  class to deal with any numericals in H&W SV models
//-----------------------------------------------------------------------------
class HWSVNumericals
{
public:
    enum RK4ParamsType
    {
        RK4NbStepsPerYear=0,
        RK4NbStepsMin,
		RK4NbStepsMax,
		RK4RefKappa,
		RK4FactorKappa,
		RK4FactorTime,
		RK4FactorNu,
		RK4NbParams
	};

    enum RK5ParamsType
    {
        RK5OdePrecision=0,
        RK5Tiny,
		RK5NbParams
	};


	enum FormulaType
	{
            Heston = 0,
			Lewis,
			Unknown
	};

	enum FormulaParamsType
	{
		LimitStep=0,
		NextLimitStep,
		LastLimitStep,
		FirstNbSteps,
		NextNbSteps,
		LastNbSteps,
		IntegrationPrecision,
		FirstLimitNbSteps,
		NextLimitNbSteps,
		LastLimitNbSteps
	};

	enum FormulaParamsTypeAlias
	{
		InvFourierPoint = FirstLimitNbSteps
	};

	enum IntegrationTraceType
    {
        NbPeriods	=0,
        NbPoints,
		FirstPeriod,
		NextPeriod,
		LastPeriod,
		OscilTraceSize
	};

	static size_t HestonFormulaNbParams;
	static size_t LewisFormulaNbParams;
	static size_t LewisFormulaSONbParams;

	static size_t HestonFctMultiplier;
	static size_t LewisFctMultiplier;

	HWSVNumericals(size_t solverType=ARM_ODEFunc::RK5Adaptative,const std::vector<double>& solverParams=std::vector<double>(0),
		int formulaType=Heston,const std::vector<double>& formulaParams=std::vector<double>(0),
		double maxDecay=0.0,bool isSOFormula=false);

	HWSVNumericals(const HWSVNumericals& rhs);
    HWSVNumericals& operator = (const HWSVNumericals& rhs);
	virtual ~HWSVNumericals() {}

	size_t GetSolverType() const {return itsSolverType;}
	const std::vector<double>& GetSolverParams() const {return itsSolverParams;}
	double GetSolverParam(size_t idx) const {return itsSolverParams[idx];}
	void SetSolverParam(size_t idx, double x) {itsSolverParams[idx]=x;}
	bool IsStdFormula() const {return itsIsStdFormula;}
	void SetIsStdFormula(bool isStdFormula) {itsIsStdFormula=isStdFormula;}
	const std::vector<double>& GetFormulaParams() const {return itsFormulaParams;}
	double GetFormulaParam(size_t idx) const {return itsFormulaParams[idx];}
	double GetMaxDecay() const {return itsMaxDecay;}
	void SetMaxDecay(double maxDecay) {itsMaxDecay=maxDecay;}

	void SetIsSOFormula(bool isSOFormula) {itsIsSOFormula=isSOFormula;}

	static void BuildHestonDefaultParams(std::vector<double>& formulaParams,bool isSOFormula=false);
	static void BuildLewisDefaultParams(std::vector<double>& formulaParams,bool isSOFormula=false);

	void AddOscilTraceNbPoints(int nbPts) {itsOscilTrace[NbPoints] += nbPts;}
	void SetOscilTrace(size_t idx, double value) {itsOscilTrace[idx] = value;}
	void ResetNbStepsTrace() {itsNbStepsTrace.clear();}
	void PushNbStepsTrace(const ARM_IntVector& nbSteps) {itsNbStepsTrace.push_back(nbSteps);}
	void SetScheduleTrace(const std::vector<double>& times) {itsScheduleTrace=times;}


	inline void SolveAnalyticalHestonSystem(const vector< ARM_AnalyticalRiccati >& analyticalDatas,
					 const std::vector<double>& ImU,
					 std::vector<double>& PsiX1,std::vector<double>& PsiY1,
					 std::vector<double>& PsiX0,std::vector<double>& PsiY0,
					 double var=1.0) const;

	inline void SolveAnalyticalLewisSystem(const vector< ARM_AnalyticalRiccati >& analyticalDatas,
					 const std::vector<double>& ImU,
					 std::vector<double>& PsiX,std::vector<double>& PsiY,
					 double var=1.0, double factor=0.5) const;

	inline void IntegrateHestonSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
						 double scalet,double lnK,
						 const std::vector<double>& PsiX1,const std::vector<double>& PsiY1,
						 const std::vector<double>& PsiX0,const std::vector<double>& PsiY0,
						 double& localIntegral_1,double& localIntegral_2) const;

	inline void IntegrateLewisSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
						 double scalet,double lnK,
						 const std::vector<double>& PsiX,const std::vector<double>& PsiY,
						 double& localIntegral,
						 double factor=0.5) const;

	void SolveSystem(ARM_RiccatiHWSV& riccatiSystem,
					 vector<std::vector<double>>& solverVars,
					 std::vector<double>& PsiX1,std::vector<double>& PsiY1,
					 std::vector<double>& PsiX0,std::vector<double>& PsiY0,
					 double var=1.0) const;

	void IntegrateSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
						 double scalet,double lnK,
						 const std::vector<double>& PsiX1,const std::vector<double>& PsiY1,
						 const std::vector<double>& PsiX0,const std::vector<double>& PsiY0,
						 double& localIntegral_1,double& localIntegral_2) const;

	string toString(const string& indent, const string& nextIndent) const;

private:
	size_t			itsSolverType;
	std::vector<double>	itsSolverParams;
	bool			itsIsStdFormula; // Standard = Heston, Enhanced = Lewis
	std::vector<double>	itsFormulaParams;

	double			itsMaxDecay;	// to manage sampling of exp(MRS*t)

	bool			itsIsSOFormula; // Spread option OR caplet,floorlet,swaption

	std::vector<double> itsOscilTrace;
	vector<ARM_IntVector> itsNbStepsTrace;
	std::vector<double> itsScheduleTrace;
};

inline void HWSVNumericals::SolveAnalyticalHestonSystem(const vector< ARM_AnalyticalRiccati >& analyticalDatas,
					 const std::vector<double>& ImU,
					 std::vector<double>& PsiX1,std::vector<double>& PsiY1,
					 std::vector<double>& PsiX0,std::vector<double>& PsiY0,
					 double var) const
{
	complex< double > alpha1,alpha0,gamma1,gamma0,zero(0,0);
	complex< double > u1(-1,0),u0(0,0),v1,v0;
	register double phi,phi2,rhoTerm,thetaTerm;
	complex< double > psi1,psi0;
	for(size_t i=0;i<ImU.size();++i)
	{
		phi = ImU[i];
		u1.imag(-phi);u0.imag(-phi);	// u1=-(1+i.ImU), u0=-i.ImU
		phi2 = phi*phi;
		v1.real(phi2);v1.imag(-phi);	// v=u(1-u) => v1=ImU^2 - i.ImU
		v0.real(phi2);v0.imag(phi);		//			   v0=ImU^2 + i.ImU
		alpha1=zero,gamma1=zero,alpha0=zero,gamma0=zero;
		for(int j=analyticalDatas.size()-1;j>=0;--j)
		{
			analyticalDatas[j].SolveSystem(u1,v1,zero,alpha1,gamma1);
			analyticalDatas[j].SolveSystem(u0,v0,zero,alpha0,gamma0);
		}
		rhoTerm		= exp(gamma1.real()+alpha1.real()*var);
		thetaTerm	= gamma1.imag()+alpha1.imag()*var;
		PsiX1[i]	= cos(thetaTerm)*rhoTerm;
		PsiY1[i]	= sin(thetaTerm)*rhoTerm;

		rhoTerm		= exp(gamma0.real()+alpha0.real()*var);
		thetaTerm	= gamma0.imag()+alpha0.imag()*var;
		PsiX0[i]	= cos(thetaTerm)*rhoTerm;
		PsiY0[i]	= sin(thetaTerm)*rhoTerm;
	}
}

inline void HWSVNumericals::SolveAnalyticalLewisSystem(const vector< ARM_AnalyticalRiccati >& analyticalDatas,
					 const std::vector<double>& ImU,
					 std::vector<double>& PsiX,std::vector<double>& PsiY,
					 double var, double factor) const
{
	complex< double > alpha,gamma,zero(0,0);
	complex< double > u(-factor,0),v(0,0),w(-factor,0);
	register double phi,rhoTerm,thetaTerm;
	complex< double > psi;
	double phi2,f2 = factor*factor;
	for(size_t i=0;i<ImU.size();++i)
	{
		phi = ImU[i];
		phi2 = phi*phi;
		u.imag(phi);			// x=ImU+i.factor => u=i.x=-factor+i.ImU
		if(itsIsSOFormula)
		{
			v.real(phi2-f2);// v=x^2=ImU^2-factor^2+i.2.factor.phi
			v.imag(2*factor*phi);// w=i.x=-factor+i.ImU
			w.imag(phi);
		}
		else
		{
			/// factor=0.5 is always used
			u.real(-0.5);
			v.real(phi2+0.25);	// v=x(x-i)=ImU^2+0.25
			w.real(0.0);
		}

		alpha=zero,gamma=zero;
		for(int j=analyticalDatas.size()-1;j>=0;--j)
		{
			analyticalDatas[j].SolveSystem(u,v,w,alpha,gamma);
		}
		rhoTerm		= exp(gamma.real()+alpha.real()*var);
		thetaTerm	= gamma.imag()+alpha.imag()*var;
		PsiX[i]	= cos(thetaTerm)*rhoTerm;
		PsiY[i]	= sin(thetaTerm)*rhoTerm;
	}
}


inline void HWSVNumericals::IntegrateHestonSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
									 double scalet,double lnK,
									 const std::vector<double>& PsiX1,const std::vector<double>& PsiY1,
									 const std::vector<double>& PsiX0,const std::vector<double>& PsiY0,
									 double& localIntegral_1,double& localIntegral_2) const
{
	double scalew,t,tLnK,cos_,sin_;
	size_t i,nbPts=GL.get_order();
	localIntegral_1=0.0;
	localIntegral_2=0.0;
	for(i=0;i<nbPts;++i)
	{
		t = ImU[i];
		scalew = GL.get_weight(i) * scalet/t;
		tLnK=-t*lnK;
		cos_ = cos(tLnK);
		sin_ = sin(tLnK);

		localIntegral_1 += (PsiY1[i] * cos_ + PsiX1[i] * sin_) * scalew;
		localIntegral_2 += (PsiY0[i] * cos_ + PsiX0[i] * sin_) * scalew;
	}
}


inline void HWSVNumericals::IntegrateLewisSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
									 double scalet,double lnK,
									 const std::vector<double>& PsiX,const std::vector<double>& PsiY,
									 double& localIntegral,
									 double factor) const
{
	double scalew,t,tLnK,tK;
	size_t i,nbPts=GL.get_order();
	localIntegral=0.0;
	if(itsIsSOFormula)
	{
		/// Inversion at u=(t,factor)
		double t2,x,costK,sintK,re,im,f2=factor*factor,_2f=2*factor;
		for(i=0;i<nbPts;++i)
		{
			t = ImU[i];
			t2 = t*t;
			x = 1/(t2+f2);
			scalew = GL.get_weight(i) * scalet * x*x;
			tK=t*lnK;
			costK = cos(tK);
			sintK = sin(tK);
			re = PsiX[i] * costK - PsiY[i] * sintK;
			im = PsiY[i] * costK + PsiX[i] * sintK;
			localIntegral += ((t2-f2)*re + _2f*t*im) * scalew;
		}
	}
	else
	{
		/// Inversion at u=(t,0.5)
		for(i=0;i<nbPts;++i)
		{
			t = ImU[i];
			scalew = GL.get_weight(i) * scalet/(t*t+0.25);
			tLnK=t*lnK;
			localIntegral += (PsiX[i] * cos(tLnK) - PsiY[i] * sin(tLnK)) * scalew;
		}
	}
}


//-----------------------------------------------------------------------------
// \class ARM_HWSV
// \brief
//  Hull & White Stochastic models abstract class
//-----------------------------------------------------------------------------
class ARM_HWSV : public ARM_PricingModelIR
{
private:
	/// For swaption analytical formula correction
	ARM_PricingModelPtr itsHWModel;

	/// To save event dates and keep track of realized variance
	std::vector<double> itsModelSchedule;

protected:
	CC_IS_MUTABLE HWSVNumericals itsNumericals;
	CC_IS_MUTABLE HWSVNumericals itsNumericalsSO;

	CC_IS_MUTABLE vector< ARM_AnalyticalRiccati > itsAnalyticalDatas;

public:
	ARM_HWSV(){}

	ARM_HWSV(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params=NULL,
		int solverType=ARM_ODEFunc::RK5Adaptative,const std::vector<double>& solverParams=std::vector<double>(0),
		int formulaType=HWSVNumericals::Lewis,const std::vector<double>& formulaParams=std::vector<double>(0),
		int formulaTypeSO=HWSVNumericals::Heston,const std::vector<double>& formulaParamsSO=std::vector<double>(0),
		double maxDecay=0.0,double maxDecaySO=0.0);
	ARM_HWSV(const ARM_HWSV& rhs);
	virtual ~ARM_HWSV() {}

	void	SetAnalyticalModel(const ARM_PricingModelPtr& model) {itsHWModel=model;}
		const	ARM_PricingModelPtr& GetAnalyticalModel() const {return itsHWModel;}

	const std::vector<double>& GetModelSchedule() const { return itsModelSchedule; }

	/// Model initialisation
	virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

    // Give local drifts and variances w.r.t. a given schedule
	virtual void NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const;

	void ModelStateLocalVariances(const std::vector<double>& timeSteps,ARM_MatrixVector& localVariances) const {}

     /// General calibration methods
    virtual void Re_InitialiseCalibParams(ARM_ModelFitter& modelFitter){};
    virtual void PreProcessing(ARM_ModelFitter& modelFitter);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter) {};
    virtual void AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter){};
    virtual void AdviseCurrentCalib(ARM_ModelFitter& modelFitter){};
	virtual void ValidateCalibMethod(ARM_CalibMethod& calibMethod);


	/// Question to analyse model
    virtual bool SupportBackwardInduction() const							{return false;}
    virtual bool SupportForwardInduction()  const							{return true;}
	virtual bool SupportAnalyticMarginal()  const							{return false;}
	virtual bool IsMeanRevertingCompatible() const							{return true;}
	virtual bool NeedsToCholeskyDecomposeFactors( ) const					{return false;}
	virtual bool ClosedFormulaSwaptionFlag(	bool isConstantNominal,
											bool isConstantSpread,
											bool isConstantstrike) const	{return true;}




	virtual void PostInit(){};
	virtual std::vector<double>* ComputeModelTimes(const ARM_TimeInfoPtrVector& timeInfos ) {return new std::vector<double>(0);}
// FIXMEFRED: mig.vc8 (24/05/2007 10:46:11):cast
	virtual ARM_VectorPtr ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const {return static_cast<ARM_VectorPtr>(NULL);}


	/// Unimplemented methods
	double VarianceToTime(double var,double minTime,double maxTime) const;
	void TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const;
	void IntegratedLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts ) const;
	void EulerLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const;
	virtual void ModelStateLocalCorrels( const std::vector<double>& timeSteps,ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& map);
	ARM_GP_VectorPtr VolatilitiesAndCorrelationTimesSteps() const;
	virtual void VolatilitiesAndCorrelations( const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& vols,ARM_GP_MatrixPtr& d1Vols,ARM_GP_MatrixPtr& correls ) const;
	double IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const;
	double IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const;
	double VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const;
	ARM_VectorPtr  VanillaSpreadOptionLet(const string& curveName,
														double evalTime,
														int callPut,
														double startTime,
														double endTime,
														double resetTime,
														double payTime,
														double payPeriod,
														double notional,
														double coeffLong,
														double coeffShort,
														const std::vector<double>& strikes,
														double swapLongFloatStartTime,
														double swapLongFloatEndTime,
														const std::vector<double>& swapLongFixPayTimes,
														const std::vector<double>& swapLongFixPayPeriods,
														double swapShortFloatStartTime,
														double swapShortFloatEndTime,
														const std::vector<double>& swapShortFixPayTimes,
														const std::vector<double>& swapShortFixPayPeriods,
														const ARM_PricingStatesPtr& states) const;



	//////////////////////////////////////////////////////////////////////////////////////////////
    /// Standard ARM object support
	//////////////////////////////////////////////////////////////////////////////////////////////
	virtual ARM_Object* Clone() const = 0;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
