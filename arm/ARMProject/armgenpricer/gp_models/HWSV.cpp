/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HWSV.cpp
 *  \brief Markov Stochastic Volatility 1 factor model
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date October 2006
 */

/// this header comes first as it include some preprocessor constants
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestrip.h"


#include "gpmodels/hwsv.h"
#include "gpmodels/modelparamshwsv.h"

#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"

#include "gpnumlib/rungekutta.h"
#include "gpnumlib/argconvdefault.h"

#include "gpcalib/calibmethod.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV
///	Routine: Constructors
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_RiccatiHWSV::ARM_RiccatiHWSV(size_t solverType, bool isStdFormula)
: itsSystemIdx(0), itsSolverType(solverType),
  itsNbDerivCalls(0), itsIsStdFormula(isStdFormula)
{}

ARM_RiccatiHWSV::ARM_RiccatiHWSV(const ARM_RiccatiHWSV& rhs)
: itsSystemIdx(rhs.itsSystemIdx),itsSolverType(rhs.itsSolverType),
  itsNbDerivCalls(rhs.itsNbDerivCalls),itsIsStdFormula(rhs.itsIsStdFormula)
{}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV
///	Routine: ComputeNbSteps
///	Returns: error
///	Action :
int ARM_RiccatiHWSV::ComputeNbSteps(const std::vector<double>& solverParams)
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown ODE solver type" );
}

////////////////////////////////////////////////////
///	Class  : HWSVNumericals
/// Static variables
////////////////////////////////////////////////////

size_t HWSVNumericals::HestonFormulaNbParams	= 10;
size_t HWSVNumericals::LewisFormulaNbParams		= 7;
size_t HWSVNumericals::LewisFormulaSONbParams	= LewisFormulaNbParams+1;
size_t HWSVNumericals::HestonFctMultiplier		= 2;
size_t HWSVNumericals::LewisFctMultiplier		= 1;

////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: BuildHestonDefaultParams (static routine)
///	Returns: 
///	Action : Set default params for Heston like formula
////////////////////////////////////////////////////
void HWSVNumericals::BuildHestonDefaultParams(std::vector<double>& formulaParams,bool isSOFormula)
{
	formulaParams.resize(HestonFormulaNbParams);
	if(isSOFormula)
	{
		formulaParams[LimitStep]			= DEFAULT_OSCILLATION_LIMIT_STD_SO;
		formulaParams[NextLimitStep]		= DEFAULT_OSCILLATION_NEXT_LIMIT_STD_SO;
		formulaParams[LastLimitStep]		= DEFAULT_OSCILLATION_LAST_LIMIT_STD_SO;
		formulaParams[FirstNbSteps]			= DEFAULT_NBPTS_FIRST_OSCIL_STD_SO;
		formulaParams[NextNbSteps]			= DEFAULT_NBPTS_NEXT_OSCIL_STD_SO;
		formulaParams[LastNbSteps]			= DEFAULT_NBPTS_LAST_OSCIL_STD_SO;
		formulaParams[IntegrationPrecision]	= DEFAULT_PRECISION_SO;
		formulaParams[FirstLimitNbSteps]	= DEFAULT_NBPTS_FIRST_OSCIL_LIMIT_STD_SO;
		formulaParams[NextLimitNbSteps]		= DEFAULT_NBPTS_NEXT_OSCIL_LIMIT_STD_SO;
		formulaParams[LastLimitNbSteps]		= DEFAULT_NBPTS_LAST_OSCIL_LIMIT_STD_SO;
	}
	else
	{
		formulaParams[LimitStep]			= DEFAULT_OSCILLATION_LIMIT_STD;
		formulaParams[NextLimitStep]		= DEFAULT_OSCILLATION_LIMIT_STD;
		formulaParams[LastLimitStep]		= DEFAULT_OSCILLATION_LIMIT_STD;
		formulaParams[FirstNbSteps]			= DEFAULT_NBPTS_FIRST_OSCIL_STD;
		formulaParams[NextNbSteps]			= DEFAULT_NBPTS_NEXT_OSCIL_STD;
		formulaParams[LastNbSteps]			= DEFAULT_NBPTS_NEXT_OSCIL_STD;
		formulaParams[IntegrationPrecision]	= DEFAULT_PRECISION;
		formulaParams[FirstLimitNbSteps]	= DEFAULT_NBPTS_FIRST_OSCIL_LIMIT_STD;
		formulaParams[NextLimitNbSteps]		= DEFAULT_NBPTS_NEXT_OSCIL_LIMIT_STD;
		formulaParams[LastLimitNbSteps]		= DEFAULT_NBPTS_NEXT_OSCIL_LIMIT_STD;
	}
}

////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: BuildLewisDefaultParams (static routine)
///	Returns: 
///	Action : Set default params for Lewis like formula
////////////////////////////////////////////////////
void HWSVNumericals::BuildLewisDefaultParams(std::vector<double>& formulaParams,bool isSOFormula)
{
	formulaParams.resize(LewisFormulaNbParams>LewisFormulaSONbParams ? LewisFormulaNbParams : LewisFormulaSONbParams);
	if(isSOFormula)
	{
		formulaParams[LimitStep]			= DEFAULT_OSCILLATION_LIMIT_ENH_SO;
		formulaParams[NextLimitStep]		= DEFAULT_OSCILLATION_NEXT_LIMIT_ENH_SO;
		formulaParams[LastLimitStep]		= DEFAULT_OSCILLATION_LAST_LIMIT_ENH_SO;
		formulaParams[FirstNbSteps]			= DEFAULT_NBPTS_FIRST_OSCIL_ENH_SO;
		formulaParams[NextNbSteps]			= DEFAULT_NBPTS_NEXT_OSCIL_ENH_SO;
		formulaParams[LastNbSteps]			= DEFAULT_NBPTS_LAST_OSCIL_ENH_SO;
		formulaParams[IntegrationPrecision]	= DEFAULT_PRECISION_SO;
		formulaParams[InvFourierPoint]		= DEFAULT_INVFOURIER_POINT_ENH_SO;
	}
	else
	{
		formulaParams[LimitStep]			= DEFAULT_OSCILLATION_LIMIT_ENH;
		formulaParams[NextLimitStep]		= DEFAULT_OSCILLATION_NEXT_LIMIT_ENH;
		formulaParams[LastLimitStep]		= DEFAULT_OSCILLATION_LAST_LIMIT_ENH;
		formulaParams[FirstNbSteps]			= DEFAULT_NBPTS_FIRST_OSCIL_ENH;
		formulaParams[NextNbSteps]			= DEFAULT_NBPTS_NEXT_OSCIL_ENH;
		formulaParams[LastNbSteps]			= DEFAULT_NBPTS_LAST_OSCIL_ENH;
		formulaParams[IntegrationPrecision]	= DEFAULT_PRECISION;
		formulaParams[InvFourierPoint]		= 0.5; // always u = phi + i/2
	}
}


////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: Constructor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
HWSVNumericals::HWSVNumericals(size_t solverType,const std::vector<double>& solverParams,int formulaType,const std::vector<double>& formulaParams,
							   double maxDecay,bool isSOFormula)
:itsSolverType(solverType),itsSolverParams(solverParams),
 itsIsStdFormula(formulaType == Heston),itsFormulaParams(formulaParams),
 itsOscilTrace(std::vector<double>(OscilTraceSize,0.0)),itsIsSOFormula(isSOFormula),
 itsMaxDecay(maxDecay),itsScheduleTrace(0)
{	
	if(itsSolverType==ARM_ODEFunc::RK4Constant)
	{
		if(itsSolverParams.size() < RK4NbParams)
		{
			itsSolverParams.resize(RK4NbParams);
			itsSolverParams[RK4NbStepsPerYear]	= RK4_NBSTEPS_PER_YEAR;
			itsSolverParams[RK4NbStepsMin]		= RK4_NBSTEPS_MIN;
			itsSolverParams[RK4NbStepsMax]		= RK4_NBSTEPS_MAX;
			itsSolverParams[RK4RefKappa]		= RK4_REF_KAPPA;
			itsSolverParams[RK4FactorKappa]		= RK4_FACTOR_KAPPA;
			itsSolverParams[RK4FactorTime]		= RK4_FACTOR_TIME;
			itsSolverParams[RK4FactorNu]		= RK4_FACTOR_NU;
		}
	}
	else if(itsSolverType==ARM_ODEFunc::RK5Adaptative)
	{
		if(itsSolverParams.size() < RK5NbParams)
		{
			itsSolverParams.resize(RK5NbParams);
			itsSolverParams[RK5OdePrecision]	= RK5_ODE_PRECISION;
			itsSolverParams[RK5Tiny]			= RK5_TINY;
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown ODE solver type" );

	if(itsIsStdFormula)
	{
		if(itsFormulaParams.size() < HestonFormulaNbParams)
			BuildHestonDefaultParams(itsFormulaParams,itsIsSOFormula);
	}
	else
	{
		if(itsIsSOFormula)
		{
			if(itsFormulaParams.size() < LewisFormulaSONbParams)
				BuildLewisDefaultParams(itsFormulaParams,itsIsSOFormula);
		}
		else
		{
			if(itsFormulaParams.size() < LewisFormulaNbParams)
				BuildLewisDefaultParams(itsFormulaParams,itsIsSOFormula);
		}

	}
}
////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
HWSVNumericals::HWSVNumericals(const HWSVNumericals& rhs)
{
	itsSolverType		= rhs.itsSolverType;
	itsSolverParams		= rhs.itsSolverParams;
	itsIsStdFormula		= rhs.itsIsStdFormula;
	itsFormulaParams	= rhs.itsFormulaParams;
	itsIsSOFormula		= rhs.itsIsSOFormula;
	itsOscilTrace		= rhs.itsOscilTrace;
	itsMaxDecay			= rhs.itsMaxDecay;
	itsScheduleTrace	= rhs.itsScheduleTrace;
}


////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
HWSVNumericals& HWSVNumericals::operator=(const HWSVNumericals& rhs)
{
	if (&rhs != this)
	{ 
		this->~HWSVNumericals();
		new (this) HWSVNumericals (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : HWSVNumericals
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string HWSVNumericals::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;

	bool isAnalytic = itsMaxDecay!=0.0 && !(itsIsStdFormula && itsIsSOFormula);
	
	if(isAnalytic)
	{
		os << indent << "Analytic Riccati solving\n";
		os << indent << "MaxDecay for exp sampling = " << itsMaxDecay << "\n\n";
	}
	else
	{
		os << indent << "Numerical Riccati solving => Runge-Kutta used\n";
		os << indent << ARM_ArgConvReverse_ODESolverType.GetString(itsSolverType) << " :" << "\n";
		if(itsSolverType == ARM_ODEFunc::RK4Constant)
		{
			os << indent << "\tNbStepsPY   = " << itsSolverParams[RK4NbStepsPerYear] << "\n";
			os << indent << "\tNbStepsMin  = " << static_cast<int>(itsSolverParams[RK4NbStepsMin]) << "\n";
			os << indent << "\tNbStepsMax  = " << static_cast<int>(itsSolverParams[RK4NbStepsMax]) << "\n";
			os << indent << "\tRefKappa    = " << itsSolverParams[RK4RefKappa] << "\n";
			os << indent << "\tFactorKappa = " << itsSolverParams[RK4FactorKappa] << "\n";
			os << indent << "\tFactorTime  = " << itsSolverParams[RK4FactorTime] << "\n";
			os << indent << "\tFactorNu    = " << itsSolverParams[RK4FactorNu] << "\n\n";
		}
		else if(itsSolverType == ARM_ODEFunc::RK5Adaptative)
		{
			os << indent << "\tPrecision = " << itsSolverParams[RK5OdePrecision] << "\n";
			os << indent << "\tTiny Val  = " << itsSolverParams[RK5Tiny] << "\n\n";
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown Riccati solver type" );
	}

	os << indent << (itsIsStdFormula ? "Standard (Heston) " : "Enhanced (Lewis) ") << "closed form formula\n";
	os << indent << "\tLimitStep            = " << itsFormulaParams[LimitStep] << "\n";
	os << indent << "\tNextLimitStep        = " << itsFormulaParams[NextLimitStep] << "\n";
	os << indent << "\tLastLimitStep        = " << itsFormulaParams[LastLimitStep] << "\n";
	os << indent << "\tFirstNbSteps         = " << itsFormulaParams[FirstNbSteps] << "\n";
	os << indent << "\tNextNbSteps          = " << itsFormulaParams[NextNbSteps] << "\n";
	os << indent << "\tLastNbSteps          = " << itsFormulaParams[LastNbSteps] << "\n";
	os << indent << "\tIntegrationPrecision = " << itsFormulaParams[IntegrationPrecision] << "\n";
	if(itsIsStdFormula)
	{
		os << indent << "\tFirstLimitNbSteps    = " << itsFormulaParams[FirstLimitNbSteps] << "\n";
		os << indent << "\tNextLimitNbSteps     = " << itsFormulaParams[NextLimitNbSteps] << "\n";
		os << indent << "\tLastLimitNbSteps     = " << itsFormulaParams[LastLimitNbSteps] << "\n";
	}
	else if(itsIsSOFormula)
	{
		os << indent << "\tFT Inversion point = " << itsFormulaParams[InvFourierPoint] << "\n";
	}

	os << indent << "Oscillatory Integral Trace" << "\n";
	os << indent << "\tNbPeriods    = " << itsOscilTrace[NbPeriods] << "\n";
	os << indent << "\t1st Period   = "	<< itsOscilTrace[FirstPeriod] << "\n";
	os << indent << "\t2nd Period   = "	<< itsOscilTrace[NextPeriod] << "\n";
	os << indent << "\tNext Periods = "	<< itsOscilTrace[LastPeriod] << "\n";
	os << indent << "\tNbPoints     = "	<< itsOscilTrace[NbPoints] << "\n";

	size_t i;
	for(i=0;i<itsNbStepsTrace.size();++i)
		os << indent << itsNbStepsTrace[i].toString(indent,nextIndent);
	os << "\n\n";

	if(itsScheduleTrace.size()>0)
	{
		os << indent << "Model Sampling Trace" << "\n";
		for(i=0;i<itsScheduleTrace.size();++i)
			os << indent << itsScheduleTrace[i] << "\t";
		os << "\n\n";
	}

    return os.str();
}

////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: SolveSystem
///	Returns: 
///	Action : Solve Riccati system and compute characteristic
///			 function Psi(u)=PsiX(u) + i.PsiY(u)
///			 Computation done for u=(1.0,phi) & u=(0,phi)
////////////////////////////////////////////////////
void HWSVNumericals::SolveSystem(ARM_RiccatiHWSV& riccatiSystem,
								vector<std::vector<double>>& solverVars,
								std::vector<double>& PsiX1,std::vector<double>& PsiY1,
								std::vector<double>& PsiX0,std::vector<double>& PsiY0,
								double var) const
{
	size_t systemSize		= riccatiSystem.GetSystemSize();
	size_t phiIdx,nbPhis	= riccatiSystem.GetNbPhis();

	/// Reset trace of number of calls to derivatives function
	riccatiSystem.ResetNbDerivCalls();

	size_t i,nbFcts = nbPhis * 4 *
			(itsIsStdFormula ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);

	if(nbFcts != solverVars[0].size())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			" : internal inconsistency in number of Runge-Kunta to solve simultaneously" );

	std::vector<double> values(nbFcts,0.0);
	std::vector<double> t;

	double yfStart=riccatiSystem.GetYf(systemSize-1),yfEnd;

	int nbSteps;
	if(itsSolverType==ARM_ODEFunc::RK4Constant)
	{
		/// 4th order Runge-Kutta with constant time step number per year
		/// between each model schedule (vol,kappa,nu,rho and lambda are
		/// constant inside each interval)

		/// Compute number of steps for each interval
		int nbStepMax = riccatiSystem.ComputeNbSteps(itsSolverParams);

		t.reserve(nbStepMax+1);

		for(int idx=systemSize-1;idx>=1;--idx)
		{
			/// Set index to compute derivatives on the right interval
			riccatiSystem.SetSystemIdx(idx);
			yfEnd=riccatiSystem.GetYf(idx-1);

			nbSteps=riccatiSystem.GetNbSteps(idx);
			t.resize(nbSteps+1,0.0);

			RkFunction(&values,nbFcts,yfStart,yfEnd,nbSteps,&t,riccatiSystem,solverVars);

			/// Shift to previous interval
			yfStart=yfEnd;
		}
	}
	else if(itsSolverType==ARM_ODEFunc::RK5Adaptative)
	{
		/// 5th order Runge-Kutta with adaptative time steps
		double h1		= -0.3;
		double hmin		= 0.0;
		int kmax		= 1000;
		double dxsav	= 1.0/365.0; // to allow time saving in t[] with a 1day minimal step

		t.resize(kmax,0.0);

		/// One single RK pass from start to end
		int nok = 0;
		int nbad = 0;
		int kount=0;

		/// Initialise index to compute derivatives on the last interval
		/// then derivative function will move index backward
		riccatiSystem.SetSystemIdx(systemSize-1);
		//yfEnd=0.0;
		yfEnd=riccatiSystem.GetYf(0);

		/// initValues[] also contains output values
		odeint(&values,nbFcts,yfStart,yfEnd,itsSolverParams[RK5OdePrecision],
			itsSolverParams[RK5Tiny],h1,hmin,kmax,dxsav,riccatiSystem,solverVars,
			&nok,&nbad,&kount,nbSteps,&t,NULL);

		/// Save nbSteps
		riccatiSystem.SetNbSteps(ARM_IntVector(1,nbSteps));

/****
FILE* f=fopen("c:\\temp\\dumpHWSV.txt","a");
fprintf(f,"RK Schedule : nok=%5d, nbad=%5d\n",nok,nbad);
for(size_t ii=0;ii<kount;++ii)
	fprintf(f,"%9.3lf  ",t[ii]*365);
fprintf(f,"\n");
fclose(f);
****/

	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown Riccati solver type for HWSV model" );


	/// Compute the function Phi(u)=exp(gamma(u)+alpha2(u)*v(t)) for u1=(1.0,phi) & u0=(0.0,phi) 
	/// or function PhiTilde(u)=gammaTilde(u)+alphaTilde2(u)*v(t) for u=(0.0,phi)
	double expRe,Im,ReTilde,ImTilde,cosIm,sinIm;
	size_t fctIncr;
	if(itsIsStdFormula)
	{
		/// Heston like formula : compute fct at u1=(1.0,phi) & u0=(0.0,phi) or both fcts at u=(0.0,phi)
		fctIncr = 4 * HWSVNumericals::HestonFctMultiplier;
		if(var!=1.0)
		{
			for(phiIdx=0,i=0;phiIdx<nbPhis;++phiIdx,i+=fctIncr)
			{
				values[i]	*= var; /// Re(alpha2(u=1+i.phi)) ou Re(alpha2(u=i.phi))
				values[i+1]	*= var; /// Im(alpha2(u=1+i.phi)) ou Im(alpha2(u=i.phi))
				values[i+4]	*= var; /// Re(alpha2(u=i.phi)) ou Re(alphaTilde2(u=i.phi))
				values[i+5]	*= var; /// Im(alpha2(u=i.phi)) ou Im(alphaTilde2(u=i.phi))
			}

		}
		if(itsIsSOFormula)
		{
			for(phiIdx=0,i=0;phiIdx<nbPhis;++phiIdx,i+=fctIncr)
			{
				/// Phi
				expRe=exp(values[i]+values[i+2]);
				Im=values[i+1]+values[i+3];
				cosIm=cos(Im);
				sinIm=sin(Im);

				PsiX1[phiIdx] = expRe*cosIm;
				PsiY1[phiIdx] = expRe*sinIm;


				/// Phi tilde
				ReTilde=values[i+4]+values[i+6];
				ImTilde=values[i+5]+values[i+7];

				PsiX0[phiIdx] = expRe*(cosIm*ReTilde - sinIm*ImTilde);
				PsiY0[phiIdx] = expRe*(cosIm*ImTilde + sinIm*ReTilde);
			}
		}
		else
		{
			for(phiIdx=0,i=0;phiIdx<nbPhis;++phiIdx,i+=fctIncr)
			{
				/// Phi(1+i.phi)
				expRe=exp(values[i]+values[i+2]);
				Im=values[i+1]+values[i+3];
				PsiX1[phiIdx] = expRe*cos(Im);
				PsiY1[phiIdx] = expRe*sin(Im);

				/// Phi(i.phi)
				expRe=exp(values[i+4]+values[i+6]);
				Im=values[i+5]+values[i+7];
				PsiX0[phiIdx] = expRe*cos(Im);
				PsiY0[phiIdx] = expRe*sin(Im);
			}
		}
	}
	else if(itsIsSOFormula)
	{
		/// Lewis like formula : computation at u=(phi,0.5)
		/// To DO
	}
	else
	{
		/// Lewis like formula : computation at u=(phi,0.5)
		fctIncr = 4 * HWSVNumericals::LewisFctMultiplier;
		if(var!=1.0)
		{
			for(phiIdx=0,i=0;phiIdx<nbPhis;++phiIdx,i+=fctIncr)
			{
				values[i]	*= var; /// Re(beta2(u=phi+i/2))
				values[i+1]	*= var; /// Im(beta2(u=phi+i/2))
			}
		}
		for(phiIdx=0,i=0;phiIdx<nbPhis;++phiIdx,i+=fctIncr)
		{
			expRe=exp(values[i]+values[i+2]);
			Im=values[i+1]+values[i+3];
			PsiX1[phiIdx] = expRe*cos(Im);
			PsiY1[phiIdx] = expRe*sin(Im);
		}
	}
}


////////////////////////////////////////////////////
///	Class  : HWSVNumericals
///	Routine: IntegrateSystem
///	Returns: 
///	Action : Integrate Heston or Lewis function of
///			 input Riccati solutions w.r.t. input
///			 Gauss-Legendre sampling
////////////////////////////////////////////////////
void HWSVNumericals::IntegrateSystem(GaussLegendre_Coefficients& GL,const std::vector<double>& ImU,
									 double scalet,double lnK,
									 const std::vector<double>& PsiX1,const std::vector<double>& PsiY1,
									 const std::vector<double>& PsiX0,const std::vector<double>& PsiY0,
									 double& localIntegral_1,double& localIntegral_2) const
{
	if(itsIsStdFormula)
		/// Heston formula
		IntegrateHestonSystem(GL,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,localIntegral_1,localIntegral_2);
	else
		/// Lewis formula
		IntegrateLewisSystem(GL,ImU,scalet,lnK,PsiX1,PsiY1,localIntegral_1);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV::ARM_HWSV(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params,int solverType,const std::vector<double>& solverParams,int formulaType,const std::vector<double>& formulaParams,int formulaTypeSO,const std::vector<double>& formulaParamsSO,
				   double maxDecay,double maxDecaySO)
:ARM_PricingModelIR(zc,params),itsNumericals(solverType,solverParams,formulaType,formulaParams,maxDecay),
 itsNumericalsSO(solverType,solverParams,formulaTypeSO,formulaParamsSO,maxDecaySO),itsHWModel(NULL),
 itsModelSchedule(0)
{}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV::ARM_HWSV(const ARM_HWSV& rhs)
: ARM_PricingModelIR(rhs) 
{
	itsHWModel		= ARM_PricingModelPtr((rhs.itsHWModel != ARM_PricingModelPtr(NULL)) ? (ARM_PricingModel*) rhs.itsHWModel->Clone() : NULL);
	itsNumericals	= rhs.itsNumericals;
	itsNumericalsSO	= rhs.itsNumericalsSO;
	itsModelSchedule= rhs.itsModelSchedule;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HWSV::PreProcessing(ARM_ModelFitter& modelFitter)
{ 
    dynamic_cast<ARM_ModelParamsHWSV*>(GetModelParams())->GenerateSchedule();
}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routine : ValidateCalibMethod
///	Returns : void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_HWSV::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWSV::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : numerical method not set in HWSV models !");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : numeraire not set in the HWSV models!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

		itsModelSchedule = *ptimeSteps;

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: NumMethodStatesGlobalVariances
///	Returns: A vector of global VCV matrixes
///	Action : Compute global VCV matrix of the state 
/// variables between each time step
/// (default implementation => identity matrix)
////////////////////////////////////////////////////
void ARM_HWSV::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double nextStep, step=timeSteps[0];
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : localDrifts.size() != offsetIndex" );
#endif

	variances.resize(nbSteps*(modelNb+1));

	size_t i,nbFactors = FactorCount();

	ARM_GP_Matrix identity(nbFactors,nbFactors,0.0);

	for(i=0;i<nbFactors;++i)
		identity(i,i) = 1.0;

	/// fills the variance
    variances[offsetIndex2+0]=static_cast<ARM_GP_Matrix*>(identity.Clone());

    for(i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
        
		/// [i+1] => variance from 0 -> ti+1
        variances[offsetIndex2+i+1] = static_cast<ARM_GP_Matrix*>(identity.Clone());
        step=nextStep;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: IntegratedLocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_HWSV::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" IntegratedLocalDrifts : unimplemented function for ARM_HWSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_HWSV::EulerLocalDrifts(const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" EulerLocalDrifts : unimplemented function for ARM_HWSV Model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_HWSV::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" TreeStatesToModelStates : unimplemented function for ARM_HWSV Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_HWSV::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" VarianceToTime : unimplemented function for ARM_HWSV Model!");

}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : void
///	Action  : VolatilitiesAndCorrelationTimesSteps for PDE
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_HWSV::VolatilitiesAndCorrelationTimesSteps() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" VolatilitiesAndCorrelationTimesSteps : unimplemented function for ARM_HWSV Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_HWSV::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" VolatilitiesAndCorrelations : unimplemented function for ARM_HWSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: IntegratedBondSquaredVol
///	Returns : double
///	Action  : Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_HWSV::IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" IntegratedBondSquaredVol : unimplemented function for ARM_HWSV Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: IntegratedBondCovariance
///	Returns : double
///	Action  : Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_HWSV::IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" IntegratedBondCovariance : unimplemented function for ARM_HWSV Model!");

}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: VolatilityScalarProduct
///	Returns : double (qui l'eut cru)
///	Action  :  Int_startTime^endTime Gamma(s,bondMaturity) * dW_s 
////////////////////////////////////////////////////
double ARM_HWSV::VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" VolatilityScalarProduct : unimplemented function for ARM_HWSV Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: ModelStateLocalCorrels
///	Returns : void
///	Action  :  
////////////////////////////////////////////////////
void ARM_HWSV::ModelStateLocalCorrels( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& multiassets )
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" ModelStateLocalCorrels : unimplemented function for ARM_HWSV Model!");


}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_HWSV::VanillaSpreadOptionLet(const string& curveName,
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
													const ARM_PricingStatesPtr& states) const
{
    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
		" VanillaSpreadOptionLet : unimplemented function for ARM_HWSV Model!");
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

