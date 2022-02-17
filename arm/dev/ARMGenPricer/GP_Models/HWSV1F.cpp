/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV1F.cpp
 *  \brief Markov Stochastic Volatility 1 factor model
 *
 *	\author  A Triki, J-M Prié
 *	\version 1.0
 *	\date October 2005
 */

/*--------------------------------------------------------------
  --------------------------------------------------------------

			1 factor Hull & White with stochastic vol

dX		= [phi(t) - mrs1.X1]dt + sig(t).sqrt(V(t)).dW
dV		= kappa.[theta(t)-V(t)]dt + nu(t).sqrt(V(t)).dZ

with 

dphi	= [sig^2(t) - 2.mrs].dt

<dW,dZ> = rho(t).dt
--------------------------------------------------------------
--------------------------------------------------------------*/

/// this header comes first as it include some preprocessor constants
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/hwsv1f.h"
#include "gpmodels/hwsv.h"
#include "gpmodels/modelparamshwsv.h"
#include "gpmodels/hw1f.h"


/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/cplx.h"
#include "gpbase/timer.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/irrate.h"


/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/normal.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpnumlib/gaussiananalytics.h"


/// gpnumlib
#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/numfunction.h"
#include "gpnumlib/rungekutta.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/ran2.h"



/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )

const int X_VARIABLE	= 0; // for Xt
const int V_VARIABLE	= 1; // for Vt
const int PHI_VARIABLE	= 2; // for phit

const int DEEP_OTM			= -1;
const int DEEP_ITM			= 1;
const int STD_ATM			= 0;

//#define PRICING_NB_TIMES	1000

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV1F
///	Routine: Constructors & affectation operator
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_RiccatiHWSV1F::ARM_RiccatiHWSV1F(const vector< ARM_SystemData >& systemDatas,size_t solverType, bool isStdFormula)
: ARM_RiccatiHWSV(solverType,isStdFormula),itsSystemDatas(systemDatas)
{
	size_t sysSize = systemDatas.size() > 0 ? systemDatas.size()-1 : 0;
	SetSystemIdx(sysSize);
	SetNbSteps(ARM_IntVector(sysSize,RK4_NBSTEPS_MIN));
}

ARM_RiccatiHWSV1F::ARM_RiccatiHWSV1F(const ARM_RiccatiHWSV1F& rhs)
: ARM_RiccatiHWSV(rhs),itsSystemDatas(rhs.itsSystemDatas)
{}

ARM_RiccatiHWSV1F& ARM_RiccatiHWSV1F::operator=(const ARM_RiccatiHWSV1F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_RiccatiHWSV1F();
		new (this) ARM_RiccatiHWSV1F (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV1F
///	Routine: ComputeNbSteps
///	Returns: 
///	Action : Compute number of steps for each interval
///			 of the model schedule (constant model
///			 parameter inside)
////////////////////////////////////////////////////
int ARM_RiccatiHWSV1F::ComputeNbSteps(const ARM_GP_Vector& solverParams)
{
	size_t nbTimes = itsSystemDatas.size();
	itsNbSteps.resize(nbTimes);
	int i,nbStepMax=solverParams[HWSVNumericals::RK4NbStepsMin];
	double nbStepsPerYear = solverParams[HWSVNumericals::RK4NbStepsPerYear];
	double yfStart = itsSystemDatas[nbTimes-1].itsYf,yfEnd;
	double nuFactor;
	for(i=nbTimes-1;i>=1;--i)
	{
		yfEnd=itsSystemDatas[i-1].itsYf;

		nuFactor = itsSystemDatas[i].itsModA*solverParams[HWSVNumericals::RK4FactorNu];
		if(nuFactor < 1.0)
			nuFactor = 1.0;

		if(itsSystemDatas[i].itsKappaTheta < solverParams[HWSVNumericals::RK4RefKappa])
			itsNbSteps[i] = static_cast<int>(floor(0.5+(yfStart-yfEnd)*nuFactor*(nbStepsPerYear +
			solverParams[HWSVNumericals::RK4FactorKappa] *
			(solverParams[HWSVNumericals::RK4RefKappa]-itsSystemDatas[i].itsKappaTheta))));
		else
			itsNbSteps[i] = static_cast<int>(floor(0.5+(yfStart-yfEnd)*nuFactor*nbStepsPerYear));
		if(itsNbSteps[i]<solverParams[HWSVNumericals::RK4NbStepsMin])
			itsNbSteps[i]=solverParams[HWSVNumericals::RK4NbStepsMin];
		else if(itsNbSteps[i]>solverParams[HWSVNumericals::RK4NbStepsMax])
			itsNbSteps[i]=solverParams[HWSVNumericals::RK4NbStepsMax];

		nbStepMax = (nbStepMax < itsNbSteps[i] ? itsNbSteps[i] : nbStepMax);
		yfStart = yfEnd;
	}
	return nbStepMax;
}


////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV1F
///	Routine: derivs
///	Returns: 
///	Action : Compute derivatives of the Riccati system
////////////////////////////////////////////////////
void ARM_RiccatiHWSV1F::derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const
{
	++itsNbDerivCalls;

	/// For each postive number phi, we compute derivatives of :
	/// for U=(1.0,phi) or (phi,0.5) :
	///		Re(Alpha2)	: [idx+0]
	///		Im(Alpha2)	: [idx+1]
	///		Re(Gamma)	: [idx+2]
	///		Im(Gamma)	: [idx+3]
	///
	/// for U=(0.0,phi) :
	///		Re(Alpha2)	: [idx+4]
	///		Im(Alpha2)	: [idx+5]
	///		Re(Gamma)	: [idx+6]
	///		Im(Gamma)	: [idx+7]


	if(itsSolverType == ARM_ODEFunc::RK5Adaptative)
		/// Search system index (starting with the current one)
		LocateSystemIdx(t);


	size_t fctMultiplier = 4 * (itsIsStdFormula ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);
	size_t i,j,nbPhis = (yt->size()) / fctMultiplier;
	size_t fctIncr = fctMultiplier;

//FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
//fprintf(f,"#26 : t=%8.3lf\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\n",t,
//		(*yt)[fctIncr*26],(*yt)[fctIncr*26+1],(*yt)[fctIncr*26+2],(*yt)[fctIncr*26+3],
//		(*yt)[fctIncr*26+4],(*yt)[fctIncr*26+5],(*yt)[fctIncr*26+6],(*yt)[fctIncr*26+7]);
//fprintf(f,"#27 : t=%8.3lf\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\n\n",t,
//		(*yt)[fctIncr*27],(*yt)[fctIncr*27+1],(*yt)[fctIncr*27+2],(*yt)[fctIncr*27+3],
//		(*yt)[fctIncr*27+4],(*yt)[fctIncr*27+5],(*yt)[fctIncr*27+6],(*yt)[fctIncr*27+7]);
//fprintf(f,"t,Idx (%3d)=\t%8.3lf\t%3d\n",itsSystemDatas.size()-1,t,itsSystemIdx);
//fclose(f);

	double expt			= exp(itsSystemDatas[itsSystemIdx].itsMrs*t);
	double exp2t		= expt*expt;

	double At			= itsSystemDatas[itsSystemIdx].itsModA;

	double kappaTheta	= itsSystemDatas[itsSystemIdx].itsKappaTheta;

	double BtUx0,BtUx1,BtUx05;
	double BtUy,Ctx,Cty,alpha2x0,alpha2y0,alpha2x1,alpha2y1;
	if(itsIsStdFormula)
	{
		/// Standard computation with two ODE system to solve at
		/// point U=(1.0,phi) and U=(0.0,phi)
		BtUx0	= itsSystemDatas[itsSystemIdx].itsModB1 +
				  itsSystemDatas[itsSystemIdx].itsModB2 * expt;
		BtUx1	= BtUx0 + itsSystemDatas[itsSystemIdx].itsModB2x * expt;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y[i] * expt;
			Ctx		= itsSystemDatas[itsSystemIdx].itsModCx[i] * exp2t;
			Cty		= itsSystemDatas[itsSystemIdx].itsModCy[i] * exp2t;


			/// Derivatives for U=(1.0,phi)
			alpha2x1 = (*yt)[j];
			alpha2y1 = (*yt)[j+1];
			(*dyt)[j]	= At * ( alpha2x1*alpha2x1 - alpha2y1*alpha2y1 ) +
						  BtUx1 * alpha2x1 - BtUy * alpha2y1 + Ctx;

			(*dyt)[j+1] = 2*At * alpha2x1*alpha2y1 +
						  BtUx1 * alpha2y1 + BtUy * alpha2x1 - Cty;

			(*dyt)[j+2] = -kappaTheta * alpha2x1;
			(*dyt)[j+3] = -kappaTheta * alpha2y1;


			/// Derivatives for U=(0.0,phi)
			alpha2x0 = (*yt)[j+4];
			alpha2y0 = (*yt)[j+5];
			(*dyt)[j+4]	= At * ( alpha2x0*alpha2x0 - alpha2y0*alpha2y0 ) +
						  BtUx0 * alpha2x0 - BtUy * alpha2y0 + Ctx;

			(*dyt)[j+5] = 2*At* alpha2x0*alpha2y0 +
						  BtUx0 * alpha2y0 + BtUy * alpha2x0 + Cty;

			(*dyt)[j+6] = -kappaTheta * alpha2x0;
			(*dyt)[j+7] = -kappaTheta * alpha2y0;
		}
	}
	else
	{
		/// New computation with on single ODE system to solve at
		/// point U=(phi,0.5)
		BtUx05	= itsSystemDatas[itsSystemIdx].itsModB1 +
				  ( itsSystemDatas[itsSystemIdx].itsModB2 +
				    itsSystemDatas[itsSystemIdx].itsModB2x ) * expt;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y[i] * expt;
			Ctx		= itsSystemDatas[itsSystemIdx].itsModCx[i] * exp2t;


			/// Derivatives for U=(phi,0.5)
			alpha2x1 = (*yt)[j];
			alpha2y1 = (*yt)[j+1];
			(*dyt)[j]	= At * ( alpha2x1*alpha2x1 - alpha2y1*alpha2y1 ) +
						  BtUx05 * alpha2x1 + BtUy * alpha2y1 + Ctx;

			(*dyt)[j+1] = 2*At * alpha2x1*alpha2y1 +
						  BtUx05 * alpha2y1 - BtUy * alpha2x1;

			(*dyt)[j+2] = -kappaTheta * alpha2x1;
			(*dyt)[j+3] = -kappaTheta * alpha2y1;
		} 
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV1F_SO
///	Routine: derivs
///	Returns: 
///	Action : Compute derivatives of the Riccati system
////////////////////////////////////////////////////
void ARM_RiccatiHWSV1F_SO::derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const
{
	++itsNbDerivCalls;

	/// For each postive number phi, we compute derivatives of :
	/// for U=(0,phi) :
	///		Re(Alpha2)		: [idx+0]
	///		Im(Alpha2)		: [idx+1]
	///		Re(Gamma)		: [idx+2]
	///		Im(Gamma)		: [idx+3]
	///
	/// for U=(0,phi) :
	///		Re(AlphaTilde2)	: [idx+4]
	///		Im(AlphaTilde2)	: [idx+5]
	///		Re(GammaTilde)	: [idx+6]
	///		Im(GammaTilde)	: [idx+7]


	/// Search system index (starting with the current one)
	LocateSystemIdx(t);


	size_t fctMultiplier = 4 * (itsIsStdFormula ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);
	size_t i,j,nbPhis = (yt->size()) / fctMultiplier;
	size_t fctIncr = fctMultiplier;

	double expt	= exp(itsSystemDatas[itsSystemIdx].itsMrs*t);
	double expSqt	= expt*expt;

	double At			= itsSystemDatas[itsSystemIdx].itsModA;

	double Bt			= itsSystemDatasSO[itsSystemIdx].itsModB;

	double kappaTheta	= itsSystemDatas[itsSystemIdx].itsKappaTheta;

	double BtUx,BtTildeUx,BtUy,Ctx,Cty,Dtx,Dty,alpha2x,alpha2y,alphaTilde2x,alphaTilde2y;
	if(itsIsStdFormula)
	{
		/// Standard computation with two ODE system to solve at
		/// point U=(0.0,phi)
		BtUx		= itsSystemDatas[itsSystemIdx].itsModB1 +
					  expt * itsSystemDatas[itsSystemIdx].itsModB2;
		BtTildeUx	= expt * itsSystemDatas[itsSystemIdx].itsModB2x;

		Dtx		= itsSystemDatasSO[itsSystemIdx].itsModD * expSqt;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y[i] * expt;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModCx[i] * expSqt;

			Cty		= itsSystemDatas[itsSystemIdx].itsModCy[i] * expSqt;

			Dty		= itsSystemDatasSO[itsSystemIdx].itsModDx[i] * expSqt;


			/// Derivatives for U=(0.0,phi) of psi system
			alpha2x = (*yt)[j];
			alpha2y = (*yt)[j+1];
			(*dyt)[j]	= At * ( alpha2x*alpha2x - alpha2y*alpha2y ) +
						  BtUx * alpha2x - BtUy * alpha2y + Ctx;

			(*dyt)[j+1] = 2*At * alpha2x*alpha2y +
						  BtUx * alpha2y + BtUy * alpha2x + Dty;

			(*dyt)[j+2] = -kappaTheta * alpha2x;
			(*dyt)[j+3] = -kappaTheta * alpha2y;


			/// Derivatives for U=(0.0,phi) of psi tilde system
			alphaTilde2x = (*yt)[j+4];
			alphaTilde2y = (*yt)[j+5];
			(*dyt)[j+4]	= Bt * ( alpha2x*alphaTilde2x - alpha2y*alphaTilde2y )
						  + BtUx * alphaTilde2x - BtUy * alphaTilde2y
						  + BtTildeUx * alpha2x + Dtx;

			(*dyt)[j+5] = Bt * ( alpha2x*alphaTilde2y + alpha2y*alphaTilde2x )
						  + BtUx * alphaTilde2y + BtUy * alphaTilde2x
						  + BtTildeUx * alpha2y - 2*Cty;

			(*dyt)[j+6] = -kappaTheta * alphaTilde2x;
			(*dyt)[j+7] = -kappaTheta * alphaTilde2y;
		}
	}
	else
	{
		/// New computation with on single ODE system to solve at
		/// point U=(phi,0.5)
/**** TO DO !!!!
		BtUx	= itsSystemDatas[itsSystemIdx].itsModB1 +
				  ( itsSystemDatas[itsSystemIdx].itsModB2 +
				    itsSystemDatas[itsSystemIdx].itsModB2x ) * expt;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y[i] * expt;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModCx[i] * expSqt;


			/// Derivatives for U=(phi,0.5)
			alpha2x = (*yt)[j];
			alpha2y = (*yt)[j+1];
			(*dyt)[j]	= At * ( alpha2x*alpha2x - alpha2y*alpha2y ) +
						  BtUx * alpha2x + BtUy * alpha2y + Ctx;

			(*dyt)[j+1] = 2*At * alpha2x*alpha2y +
						  BtUx * alpha2y - BtUy * alpha2x;

			(*dyt)[j+2] = -kappaTheta * alpha2x;
			(*dyt)[j+3] = -kappaTheta * alpha2y;
		} 
****/
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV1F::ARM_HWSV1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params,int solverType,const ARM_GP_Vector& solverParams,
					   int formulaType,const ARM_GP_Vector& formulaParams,
					   int formulaTypeSO,const ARM_GP_Vector& formulaParamsSO,
					   double maxDecay,double maxDecaySO)
:ARM_HWSV(zc,params,solverType,solverParams,formulaType,formulaParams,formulaTypeSO,formulaParamsSO,maxDecay,maxDecaySO),
 itsRealizedVar(0,0),itsVarIdx(0),itsRealizedStates(0)
{	
	ARM_ModelParamVector paramVector(2);
	ARM_CurveModelParam  paramvol = GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam();
	ARM_CurveModelParam  paramMRS = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam();
	paramVector[0] = &paramvol;
	paramVector[1] = &paramMRS;
	ARM_PricingModelPtr hwModel( ARM_PricingModelPtr (static_cast< ARM_PricingModel* >(new ARM_HullWhite1F(  zc ) ) ) );
	hwModel->SetModelParams( ARM_ModelParamsHW1FStd(paramVector) );
	SetAnalyticalModel(hwModel);
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV1F::ARM_HWSV1F(const ARM_HWSV1F& rhs)
: ARM_HWSV(rhs), itsRealizedStates(0) 
{
	itsSystemDatas		= rhs.itsSystemDatas;
    itsFunctionDatas	= rhs.itsFunctionDatas;
	itsAngularShifts1	= rhs.itsAngularShifts1;
	itsPrevValues1		= rhs.itsPrevValues1;
	itsAngularShifts2	= rhs.itsAngularShifts2;
	itsPrevValues2		= rhs.itsPrevValues2;
	itsRealizedVar		= rhs.itsRealizedVar;
	itsVarIdx			= rhs.itsVarIdx;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HWSV1F& ARM_HWSV1F::operator=(const ARM_HWSV1F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_HWSV1F();
		new (this) ARM_HWSV1F (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_HWSV1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F HW Stochastic Volatility Model\n";
    os << indent << "---------------------------------\n\n";

    os << indent << "Numericals for caplet/swaption\n";
    os << indent << "------------------------------\n";
    os << itsNumericals.toString(indent,nextIndent);
    os << indent << "\n\n";

    os << indent << "Numericals for CMS caplet/spread option\n";
    os << indent << "---------------------------------------\n";
    os << itsNumericalsSO.toString(indent,nextIndent);
    os << indent << "\n\n";


    os << ARM_PricingModel::toString(indent);

    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_HWSV1F::Clone() const
{
	return new ARM_HWSV1F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: GetMrs
///	Returns: double
///	Action : fill MRS with 0 limit value
////////////////////////////////////////////////////
double ARM_HWSV1F::GetMrs() const
{
	double mrs = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
    if(fabs(mrs)<=K_NEW_DOUBLE_TOL)
        mrs=(mrs>=0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

	return mrs;
}


///////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWSV1F::FirstPricingStates( size_t bucketSize ) const
{
	size_t nbModelStates=ModelStatesSize();
	size_t nbNumStates=FactorCount();
	ARM_PricingStatesPtr states(new ARM_PricingStates(bucketSize,nbModelStates,0,nbNumStates));
	int nbStates = states->size();

	double LongTermVar0 = 1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		LongTermVar0 = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(0.0);

	/// To save Integ{Event(i)->Event(i+1), sigma(t)^2 * V(t) * exp(2*mrs*t) * dt
	/// for MaxRate keyword
	size_t i,nbEvents = GetModelSchedule().size();
	if(nbEvents==0)
		nbEvents = 1; // one fictive event to be consistent with itsVarIdx=0
	itsVarIdx = 0;
	itsRealizedVar.resize(nbEvents,nbStates);
	itsRealizedStates.resize(nbEvents);

	/// Skip last event because we always need previous event states
	for(i=0;i+1<nbEvents;++i)
		itsRealizedStates[i] = ARM_PricingStatesPtr(new ARM_PricingStates(bucketSize,nbModelStates));

	for(i=0;i<nbStates;++i)
	{
		itsRealizedVar(itsVarIdx,i)=0.0;
		states->SetModelState(i,V_VARIABLE,LongTermVar0); // V(0)=LongTermVar(0)
	}

	return states;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_HWSV1F::ValidateModelParams(const ARM_ModelParams& params) const
{
    if(!params.DoesModelParamExist(ARM_ModelParamType::MeanReversion) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::Volatility) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::Correlation))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "At least 1 Model Param is not of a good type!");
    }
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_HWSV1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    ARM_GP_Vector  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
    case ARM_ModelParamType::Volatility:
	case ARM_ModelParamType::Correlation:
	case ARM_ModelParamType::VolOfVol:
	case ARM_ModelParamType::LongTermVol:
        {
            double date = portfolio->GetAsset(0)->GetResetDates()->Elt(0) - asOfDate;
            tmpdates.push_back(date);
            for(i=1; i<size1; i++) 
            {
                double resetlag = portfolio->GetAsset(i)->GetResetDates()->Elt(0) - asOfDate;
                if(fabs (date - resetlag) > FRMVOL_LAG_THRESHOLD)
                {
                    tmpdates.push_back(resetlag);
                    date = resetlag;
                }
				else
				{
					/// ignore this instrument
					portfolio->SetWeight(0.0,i);
				}
            }
			modelParam->UpdateValues(&tmpdates);
        }
        break;
    case ARM_ModelParamType::MeanReversion:
        {
            double date = portfolio->GetAsset(0)->GetFlowEndDates()->Elt(0) - asOfDate; 
            tmpdates.push_back(date);
            for(i=1; i<size1; i++)
            {
                double startlag = portfolio->GetAsset(i)->GetFlowEndDates()->Elt(0) - asOfDate;
                if(fabs (date - startlag) > FRMVOL_LAG_THRESHOLD)
                {
                    tmpdates.push_back(startlag);
                    date = startlag;
                }
				else
				{
					/// ignore this instrument
					portfolio->SetWeight(0.0,i);
				}
            }
			modelParam->UpdateValues(&tmpdates);
        }
    default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... Model Param Not Supported by MSV1F" );
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of ND matrix
///	Action : Compute local variances of the state variable 
/// between each time step
////////////////////////////////////////////////////
void ARM_HWSV1F::NumMethodStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelRank= GetModelRank();
    double step		= timeSteps[0],
		   nextStep, rho_t;
	size_t offsetIndex	= (nbSteps-1)*modelRank;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelRank+1));
	size_t i,j;

	size_t factorNb = FactorCount();

	ARM_GP_Matrix correlMatrix(factorNb,factorNb,1.0);
	// All the variance is in the numerical method

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// Interpolate Correlation
		rho_t = GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve()->Interpolate(step);
		/// [i] => local variance from ti->ti+1
		for (int l = 0; l < factorNb; ++l)
		{	
			for(j = 0; j < l; ++j)
			{
				correlMatrix(l,j) = correlMatrix(j,l) = rho_t;
			}
		}
		localVariances[offsetIndex+i] = static_cast<ARM_GP_Matrix*>(correlMatrix.Clone());
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_HWSV1F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	bool isLNApprox = true;

#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : time index bigger than max size!");
#endif

	const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
	const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

	double time = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt = (nextTime - time)/K_YEAR_LEN;
	size_t factorsNb= FactorCount();
	size_t statesNb = states->size();
	double X_State,V_State,Phi_State;
	size_t modelNb	= GetModelNb();

	/// Locate nextTime in the model schedule (=events)
	const ARM_GP_Vector& modelSchedule = GetModelSchedule();
	size_t nbEvents = modelSchedule.size();
	for(size_t varIdx=0;varIdx<nbEvents;++varIdx)
	{
		if(modelSchedule[varIdx]>nextTime-K_NEW_DOUBLE_TOL)
			break;
	}
	if(varIdx-itsVarIdx > 1 || varIdx>=nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : pb in finding index in model schedule");

	bool isResetVar = varIdx == itsVarIdx+1;
	if(isResetVar)
	{
		/// Reset realized variance to previous one
		for(size_t i=0;i<statesNb;++i)
			itsRealizedVar(varIdx,i)=itsRealizedVar(itsVarIdx,i);
		itsVarIdx = varIdx;
	}

	/// Constant Rate Mean Reversion
	double Lambda_t		= GetMrs();
	double Lambda_dt	= Lambda_t * dt;

	if(time<K_NEW_DOUBLE_TOL)
			time+=INTERPOL_EPS;

	/// Time dependent Volatility
	double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->Interpolate(time);
	
	/// Possible time depedent Long Term Vol
	double LongTermVol_t = 1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		LongTermVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(time);


	/// Time dependent Vol of Vol
	double VolOfVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(time);	
	
	/// Time dependent Rate Vol Scaling and Rate/Variance correlation
	double rho_t = GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve()->Interpolate(time);	
	double vol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->Interpolate(time);	

	/// Possible Variance Mean Reversion
	double VolMeanReversion_t;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::VolMeanReversion))
		VolMeanReversion_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()->Interpolate(time);
	else
		VolMeanReversion_t  = Lambda_t - vol_t  * rho_t * VolOfVol_t / Lambda_t;

	bool isTerminalZc = (GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc ||
						 GetNumeraire()->GetType() == ARM_Numeraire::TerminalEventZc);

	double coef_var_2 = 1.0;
	if(isTerminalZc)
	{
		double numTime = GetNumeraire()->GetMaturity();
		/// Average beta to avoid time dependent mean reversion
		double avgeBetatTStar = exp(-Lambda_t*(numTime-nextTime)/K_YEAR_LEN)-exp(-Lambda_t*(numTime-time)/K_YEAR_LEN);
		avgeBetatTStar = (1.0 - avgeBetatTStar/(Lambda_t*dt))/Lambda_t;
		coef_var_2 = 1 + VolOfVol_t*rho_t*vol_t*avgeBetatTStar/VolMeanReversion_t;
		VolMeanReversion_t *= coef_var_2;
		if(VolMeanReversion_t < 0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : negative variance mean reversion !");

		if(VolMeanReversion_t < K_NEW_DOUBLE_TOL)
			VolMeanReversion_t = K_NEW_DOUBLE_TOL;
	}
	double VolMeanReversion_dt = VolMeanReversion_t * dt;


	/// Variables of Andersen's diffusion methods
	/// var_1 = exp(-Kappa * dt) 
	/// var_2 = 1-exp(-Kappa * dt) 
	/// var_3 = Epsilon²/Kappa
	/// VolOfVol_2_t = Epsilon²

	double var_1		= exp(-VolMeanReversion_dt);
	double var_2		= 1 - var_1;
	double var_2_2		= LongTermVol_t*var_2/coef_var_2;
	double VolOfVol_2_t	= VolOfVol_t*VolOfVol_t;
	double var_3		= VolOfVol_2_t/VolMeanReversion_t;

	double var_5		= exp(-Lambda_dt);
	double var_5_2		= var_5*var_5;

	double var_4	= var_2*var_3;
	double var_7	= var_1*var_4;
	double var_8	= 0.5*var_2_2*var_4;

	double coefPhi = exp(2*Lambda_t*nextTime/K_YEAR_LEN);

//FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
//fprintf(f,"t=\t%8.3lf\n",time);

	if(isLNApprox)
	{
		double var_7	= 0.5 * var_3 * (1 - var_1*var_1);
		double sdt		= sqrt(dt);

		double KX	= var_5;
		double KP	= var_5_2;
		double KPV	= (1-var_5_2)/(2*Lambda_t);
		double KXV	= sqrt(KPV);

		/// V(t+dt) conditional to V(t) is assumed to be a lognormal
		double dx,dv,mean,var,stdDev,Eta,LastPhi_State;
		for( size_t i=0;i<statesNb; ++i )
		{
			X_State = states->GetModelState(i,modelNb+X_VARIABLE);
			V_State = states->GetModelState(i,modelNb+V_VARIABLE);
			Phi_State = states->GetModelState(i,modelNb+PHI_VARIABLE);

			Eta = Volatility_t * sqrt(V_State);

			if(isResetVar)
			{
				/// Save states at current time for MaxRate keyword computation
				/// varIdx is the index of next event time
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+X_VARIABLE,X_State);
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+V_VARIABLE,V_State);
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+PHI_VARIABLE,Phi_State);
			}

			/// V(t) : lognormal mean & variance are fitted
			dv = states->GetNumMethodState(i,modelNb+V_VARIABLE);
			mean = V_State*var_1 + var_2_2;
			var = log(1 + V_State*var_7/(mean*mean)); // V(s)=V(ti)
			//var	= log(1 + (V_State*var_7 + var_8)/(mean*mean)); // V(s)=Eti(V(s))
			stdDev = sqrt(var);

			V_State = mean*exp(-0.5*var + stdDev*dv);
			if(V_State < 0.0)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Variance process is negative !");
			states->SetModelState(i,modelNb+V_VARIABLE,V_State);

			/// Phi Variable for reconstruction Zc formula : Euler scheme
			LastPhi_State = Phi_State;

			/// H&W X(t) process : Euler Scheme
			dx = states->GetNumMethodState(i,modelNb+X_VARIABLE);
			if(isTerminalZc)
			{
				//KP = exp(-2*Lambda_t*dt)
				//Phi_State = KP*Phi_State + (1-KP)/(2*Lambda_t)*Eta^2;
				//X_State = X_State - Lambda_t * X_State*dt + sdt*dx*Eta ;
				Phi_State = KP*Phi_State + KPV*Eta*Eta;
				X_State = KX*X_State + KXV*Eta*dx;
			}
			else
			{
				//Phi_State = Phi_State + (Eta*Eta - 2*Lambda_t*Phi_State)*dt;
				//X_State = X_State + (Phi- Lambda_t * X_State)*dt + sdt*dx*Eta ;
				Phi_State = Phi_State + (Eta*Eta - 2*Lambda_t*Phi_State)*dt;
				X_State = X_State + (0.5*(Phi_State + LastPhi_State) - Lambda_t * X_State)*dt + sdt*dx*Eta ;
			}
			states->SetModelState(i,modelNb+PHI_VARIABLE,Phi_State);
			states->SetModelState(i,modelNb+X_VARIABLE,X_State);

			/// To save Integ{Event(i)->Event(i+1), sigma(t)^2 * V(t) * exp(2*lambda*t) * dt
			/// for MaxRate keyword computation
			itsRealizedVar(varIdx,i) += KPV*Eta*Eta*coefPhi;

	//		fprintf(f,"#%5d\tx(t),v(t),phi(t)=\t%15.10lf\t%15.10lf\t%15.10lf\n",i,X_State,V_State,Phi_State);
		}
	}
	else
	{
		/// V(t+dt) conditional to V(t) is assumed to be a khi squared
		double blendPhi1= 0.5;
		double blendPhi2= 1-blendPhi1;
		double blendV1	= 0.5;
		double blendV2	= 1-blendV1;

		double var_6	= 1-rho_t*rho_t;
		double var_9	= 1/sqrt(var_6); 

		double var_10	= VolMeanReversion_t/(VolOfVol_t*var_2);

		double XPhiFactor	= (1-var_5)/Lambda_t;
		double phiVFactor	= (1-var_5_2)/(2*Lambda_t)*vol_t*vol_t;
		double XVFactor		= XPhiFactor*vol_t*rho_t;
		double XVOrthFactor	= var_6*phiVFactor;

		/// Phi diffusion coefficients for Phi(ti), V(ti) and V(ti+1)
		double KP	= var_5_2;
		double KPV1	= phiVFactor*blendV1; // V(ti)
		double KPV2	= phiVFactor*blendV2; // V(ti+1)

		/// X diffusion coefficients for X(ti)
		double KX	= var_5;

		/// X diffusion coefficients for Phi(ti) and Phi(ti+1),
		double KXP1		= XPhiFactor*blendPhi1; // Phi(ti)
		double KXP2		= XPhiFactor*blendPhi2; // Phi(ti+1)

		/// X diffusion coefficients for V(ti) and V(ti+1) in V colinear part
		double KXV2		= XVFactor*var_10; // V(ti+1)
		double KXV0		= -KXV2*var_2_2;
		double KXV1		= -KXV2*var_1; // V(ti)

		/// X diffusion coefficients for V(ti) and V(ti+1) in V orthognal part
		double KXVOrth1	= XVOrthFactor*blendV1; // V(ti)
		double KXVOrth2	= XVOrthFactor*blendV2; // V(ti+1)


		double dx,dv,mean,var,psi,a,b,c,u;
		double LastPhi_State,LastV_State;
		for( size_t i=0;i<statesNb; ++i )
		{
			X_State = states->GetModelState(i,modelNb+X_VARIABLE);
			V_State = states->GetModelState(i,modelNb+V_VARIABLE);
			Phi_State = states->GetModelState(i,modelNb+PHI_VARIABLE);

			if(isResetVar)
			{
				/// Save states at current time for MaxRate keyword computation
				/// varIdx is the index of next event time
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+X_VARIABLE,X_State);
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+V_VARIABLE,V_State);
				itsRealizedStates[varIdx-1]->SetModelState(i,modelNb+PHI_VARIABLE,Phi_State);
			}

			/// Correlated draws... then uncorrelate !
			dx = states->GetNumMethodState(i,modelNb+X_VARIABLE);
			dv = states->GetNumMethodState(i,modelNb+V_VARIABLE);
			dv = (dv-rho_t*dx)*var_9;

			/// V variable diffusion : Khi2 approximation
			LastV_State = V_State;
			mean= V_State*var_1 + var_2_2;
			var	= V_State*var_7 + var_8;
			psi = var/(mean*mean);

			if(psi<=1.5)
			{
				c = 2/psi-1;
				b = sqrt(c+sqrt(c*(c+1)));
				a = mean/(1+b*b);
				V_State = b+dv;
				V_State = a*V_State*V_State;
			}
			else
			{
				a = (psi-1)/(psi+1);
				b = mean/(1-a);
				u = NormalCDF(dv); /// restoring the [0,1] uniform draw
				V_State = (u<=a ? 0 : log((1-a)/(1-u))*b);
			}

			if(V_State < 0.0)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Variance process is negative !");
			states->SetModelState(i,modelNb+V_VARIABLE,V_State);

			/// Phi variable diffusion : integration with averaged exp(mrs.t) and
			/// blending of V(ti) & V(ti+1)
			LastPhi_State = Phi_State;
			Phi_State = KP*Phi_State + KPV1*LastV_State + KPV2*V_State;
			states->SetModelState(i,modelNb+PHI_VARIABLE,Phi_State);

			/// To save Integ{Event(i)->Event(i+1), sigma(t)^2 * V(t) * exp(2*mrs*t) * dt
			/// for MaxRate keyword computation
			itsRealizedVar(varIdx,i) += (KPV1*LastV_State + KPV2*V_State)*coefPhi;


			/// X variable diffusion : integration with averaged exp(mrs.t),
			/// reuse of V colinar part and blending of V(ti) & V(ti+1) in
			/// V orthogonal part
			if(isTerminalZc)
				X_State = KX*X_State + KXV0 + KXV1*LastV_State + KXV2*V_State
						  + dx*sqrt(KXVOrth1*LastV_State + KXVOrth2*V_State);
			else
				X_State = KX*X_State + KXP1*LastPhi_State + KXP2*Phi_State
						  + KXV0 + KXV1*LastV_State + KXV2*V_State
						  + dx*sqrt(KXVOrth1*LastV_State + KXVOrth2*V_State);
			states->SetModelState(i,modelNb+X_VARIABLE,X_State);
		}
	//		fprintf(f,"#%5d\tx(t),v(t),phi(t)=\t%15.10lf\t%15.10lf\t%15.10lf\n",i,X_State,V_State,Phi_State);
	}
//fclose(f);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeRiccatiSchedule
///	Returns: 
///	Action : Compute a schedule based on model schedule
///			 that ends at the input time T
////////////////////////////////////////////////////
void ARM_HWSV1F::ComputeRiccatiSchedule(double T, ARM_GP_Vector& schedule, double Tref, bool isSOFormula, double Tstart) const
{
	schedule.empty();
	ARM_GP_Vector modelSched = static_cast<const ARM_ModelParamsHWSV* const>(GetModelParams())->GetSchedule();
	if(modelSched[0] != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : HWSV1F model schedule must start at spot time (t=0)" );


	if(modelSched[modelSched.size()-1]< T - K_NEW_DOUBLE_TOL)
		modelSched.push_back(T);

	/// The schedule starts with Tstart...
	size_t l;
	for(l=0;l<modelSched.size();++l)
	{
		if(modelSched[l]>Tstart+K_NEW_DOUBLE_TOL)
		{
			schedule.push_back(Tstart);
			break;
		}
		else if(Tstart-K_NEW_DOUBLE_TOL <= modelSched[l])
			break;
	}

	///...is made of all model dates and ends with T
	for(;l<modelSched.size() && modelSched[l]<T-K_NEW_DOUBLE_TOL;++l)
		schedule.push_back(modelSched[l]);
	schedule.push_back(T);


	double maxDecay = (isSOFormula ? itsNumericalsSO.GetMaxDecay() : itsNumericals.GetMaxDecay());
	if(maxDecay>0)
	{
		/// Add points to keep |exp(-2*MRS.(Tref-t))-exp(-2*MRS.(Tref-ti))| < MaxDecay for all t in ]ti,ti+1]
		ARM_GP_Vector additionalPoints;
		double factor=-2*GetMrs()/K_YEAR_LEN;
		double lastt=schedule[0],t;
		double expLastt = exp(factor*(Tref-lastt)),expt,h,dh;
		size_t i,j,n;
		for(i=1;i<schedule.size();++i)
		{
			expt=exp(factor*(Tref-schedule[i]));
			h=expt-expLastt;
			n=static_cast<size_t>(ceil(fabs(h)/maxDecay));
			dh = h/n;
			for(h=dh,j=1;j+1<n;++j,h+=dh)
			{
				t = Tref-log(h+expLastt)/factor;
				additionalPoints.push_back(t);
			}
			lastt=schedule[i];
			expLastt=expt;
		}

		/// Merge with additional points
		ARM_GP_Vector tmpSched(schedule.size()+additionalPoints.size());
		CC_NS(std,merge)(schedule.begin(),schedule.end(),additionalPoints.begin(),additionalPoints.end(),tmpSched.begin());
		ARM_GP_Vector::iterator last=CC_NS(std,unique)(tmpSched.begin(),tmpSched.end());
		tmpSched.resize(last-tmpSched.begin());
		schedule = tmpSched;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: InitAnalyticalData
///	Returns: 
///	Action : Initialise model dependent constant set
///			 for Riccati analytical solving
///			 Model parameters are assumed to be stepwise
///			 right constant. Csts saved at [i] are then
///			 related to ]ti, ti+1] interval
///			 expT0 = exp(-mrs*T0/365.0) where
///			 T0 = swap or libor rate start time (lag from asOfDate)
///			 Te = swap or libor rate reset time (lag from asOfDate)
////////////////////////////////////////////////////
void ARM_HWSV1F::InitAnalyticalData(double evalTime, double F0, double Te, double expT0,
									bool isSOFormula, ARM_GP_Vector& schedule, double mu0) const
{
	double Tref = Te;
	if(schedule.size()==0)
		ComputeRiccatiSchedule(Te,schedule,Tref,isSOFormula,evalTime);

    size_t nbTimes = schedule.size();
    itsAnalyticalDatas.resize(nbTimes-1);

	double mrs	= GetMrs();
	double mrsYf = mrs/K_YEAR_LEN;


	double expTref		= exp(mrsYf*Tref);
	double F0expTref	= F0*expTref;
	expT0				*= expTref;
	double cstC2		= 0.5*F0expTref*F0expTref;
	double cstC1		= (isSOFormula ? mu0*expTref*expTref : 0.0);

	const ARM_Curve* nuCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve();
	const ARM_Curve* volCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve();
	const ARM_Curve* rhoCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve();
	const ARM_Curve* kappaCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve();
	ARM_Curve* thetaCurve=NULL;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		thetaCurve = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve();
	double t=schedule[0],nextt,nu,nu2,vol,rho,rhoNuVol,kappa,theta;
	double expt=exp(-mrsYf*(Tref-t)),expNextt,expProxy,expProxyVol,expProxyVol2;
    for(size_t i=0;i<nbTimes-1;++i)
    {
		nextt	= schedule[i+1];
		nu		= nuCurve->Interpolate(nextt);
		nu2		= nu*nu;
		vol		= volCurve->Interpolate(nextt);
		rho		= rhoCurve->Interpolate(nextt);

		rhoNuVol= rho*nu*vol;

		kappa	= kappaCurve->Interpolate(nextt);
		theta	= (thetaCurve!=NULL ? thetaCurve->Interpolate(nextt) : 1.0);

		expNextt	= exp(-mrsYf*(Tref-nextt));
		expProxy	= 0.5*(expt+expNextt);


		itsAnalyticalDatas[i].itsdt	= (nextt-t)/K_YEAR_LEN;
		itsAnalyticalDatas[i].itsNu2= nu2;
		itsAnalyticalDatas[i].itsA	= -0.5*nu2;
		itsAnalyticalDatas[i].itsB1	= kappa + rhoNuVol*(1-expT0*expProxy)/mrs;
		itsAnalyticalDatas[i].itsB2	= rhoNuVol*F0expTref*expProxy;

		expProxyVol		= expProxy*vol;
		expProxyVol2	= expProxyVol*expProxyVol;
		itsAnalyticalDatas[i].itsC2	= cstC2*expProxyVol2;
		if(isSOFormula)
			itsAnalyticalDatas[i].itsC1	= cstC1*expProxyVol2;
		else
			itsAnalyticalDatas[i].itsC1	= 0.0;

		itsAnalyticalDatas[i].itsD	= -kappa*theta;

		t=nextt;
		expt=expNextt;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: InitSystemData
///	Returns: 
///	Action : Initialise model dependent constant set
///			 for Riccati system solving
///			 Model parameters are assumed to be stepwise
///			 right constant. Csts saved at [i] are then
///			 related to ]ti, ti+1] interval
///			 expT0 = exp(-mrs*T0/365.0) where
///			 T0 = swap or libor rate start time (lag from asOfDate)
///			 Te = swap or libor rate reset time (lag from asOfDate)
////////////////////////////////////////////////////
void ARM_HWSV1F::InitSystemData(double evalTime,double F0, double mrs, double Te, double expT0, size_t maxSize,
								bool isStdFormula, bool isSOFormula, double Mu0) const
{
	ARM_GP_Vector schedule;
	ComputeRiccatiSchedule(Te,schedule,Te,isSOFormula,evalTime);

    size_t nbTimes = schedule.size();
    itsSystemDatas.resize(nbTimes);
    //itsSystemDatas.resize(nbTimes-1);

	double cstB	= -expT0/mrs;
	double cstC = 0.5*F0*F0;

	double cstD = - Mu0;
	if(isSOFormula)
		itsSystemDatasSO.resize(nbTimes);
		//itsSystemDatasSO.resize(nbTimes-1);

	const ARM_Curve* nuCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve();
	const ARM_Curve* volCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve();
	const ARM_Curve* rhoCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve();
	const ARM_Curve* kappaCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve();
	ARM_Curve* thetaCurve=NULL;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		thetaCurve = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve();
	double t,teps,nu,vol,rho,rhoNuVol,kappa,theta,vol2,nu2;
    for(size_t i=0;i<nbTimes;++i)
    //for(size_t i=0;i<nbTimes-1;++i)
    {
		t		= schedule[i];
		//t		= schedule[i+1];
		if(t<K_NEW_DOUBLE_TOL)
			teps=1.0e-4;
		else
			teps=t;
		nu		= nuCurve->Interpolate(teps);
		vol		= volCurve->Interpolate(teps);
		rho		= rhoCurve->Interpolate(teps);

		rhoNuVol= rho*nu*vol;

		kappa	= kappaCurve->Interpolate(teps);
		theta	= (thetaCurve!=NULL ? thetaCurve->Interpolate(teps) : 1.0);

		/// To test against Riccati analytical solution
		//kappa = mrs - rhoNuVol/mrs;

		itsSystemDatas[i].itsTime		= t;
		itsSystemDatas[i].itsYf			= t/K_YEAR_LEN;
		itsSystemDatas[i].itsMrs		= mrs;
		itsSystemDatas[i].itsKappaTheta	= kappa*theta;
		nu2 = nu*nu;
		itsSystemDatas[i].itsModA		= -0.5*nu2;
		itsSystemDatas[i].itsModB1		= kappa + rhoNuVol/mrs;
		itsSystemDatas[i].itsModB2		= rhoNuVol*cstB;
		itsSystemDatas[i].itsModB2x		= (isStdFormula ? -rhoNuVol*F0 : -0.5*rhoNuVol*F0);
		vol2 = vol*vol;
		itsSystemDatas[i].itsModC	= cstC*vol2;

		/// Reserve memory space
		itsSystemDatas[i].itsModB2y.reserve(maxSize);
		itsSystemDatas[i].itsModCx.reserve(maxSize);
		itsSystemDatas[i].itsModCy.reserve(maxSize);

		if(isSOFormula)
		{
			itsSystemDatasSO[i].itsModB		= -nu2;
			itsSystemDatasSO[i].itsModD		= cstD*vol2;
			itsSystemDatasSO[i].itsModDx.reserve(maxSize);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: UpdateUDependentSystemData
///	Returns: 
///	Action : Update Riccati constants which depend
///			 on the complex number. Imaginary parts are
///			 input because real parts are already known (0 or 1)
////////////////////////////////////////////////////
void ARM_HWSV1F::UpdateUDependentSystemData(const ARM_GP_Vector& ImU, bool isStdFormula, bool isSOFormula) const
{
    size_t i,nbTimes = itsSystemDatas.size();
	size_t j,nbU = ImU.size();
	double b2x,b2y,c,y,y2,cy,d;
	if(isStdFormula)
	{
		for(i=0;i<nbTimes;++i)
		{
			itsSystemDatas[i].itsModB2y.resize(nbU);
			itsSystemDatas[i].itsModCx.resize(nbU);
			itsSystemDatas[i].itsModCy.resize(nbU);
			if(isSOFormula)
			{
				itsSystemDatasSO[i].itsModDx.resize(nbU);
				d = itsSystemDatasSO[i].itsModD;
			}

			b2x = itsSystemDatas[i].itsModB2x;
			c = itsSystemDatas[i].itsModC;
			for(j=0;j<nbU;++j)
			{
				y = ImU[j];
				itsSystemDatas[i].itsModB2y[j]	= b2x*y;
				cy = c*y;
				itsSystemDatas[i].itsModCx[j]	= cy*y;
				itsSystemDatas[i].itsModCy[j]	= cy;

				if(isSOFormula)
					itsSystemDatasSO[i].itsModDx[j] = d*y;
			}
		}
	}
	else
	{
		for(i=0;i<nbTimes;++i)
		{
			itsSystemDatas[i].itsModB2y.resize(nbU);
			itsSystemDatas[i].itsModCx.resize(nbU);
			if(isSOFormula)
			{
				itsSystemDatasSO[i].itsModDx.resize(nbU);
				d = itsSystemDatasSO[i].itsModD;
			}

			b2y = 2 * itsSystemDatas[i].itsModB2x;
			c = itsSystemDatas[i].itsModC;
			for(j=0;j<nbU;++j)
			{
				y = ImU[j];
				y2 = y*y+0.25;
				itsSystemDatas[i].itsModB2y[j]	= b2y*y;
				itsSystemDatas[i].itsModCx[j]	= c*y2;
				if(isSOFormula)
					itsSystemDatasSO[i].itsModDx[j] = d*y2;
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: InitFunctionData
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HWSV1F::InitFunctionData(double F0,double Te,double expTe) const // exp(mrs*Te)
{
	
	/// Store parameters for each interval ]ti, ti+1]
    /// They are stepwise right constant

	ARM_GP_Vector newSchedule;
	ComputeRiccatiSchedule(Te,newSchedule,Te);

    size_t schedSize = newSchedule.size();
    itsFunctionDatas.resize(schedSize);
    itsPrevValues1.resize(schedSize);
    itsPrevValues2.resize(schedSize);
    itsAngularShifts1.resize(schedSize);
    itsAngularShifts2.resize(schedSize);
    double vol,nu,rho,kappa,theta,nuVol,a;

	double mrs=GetMrs();

	const ARM_Curve* nuCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve();
	const ARM_Curve* volCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve();
	const ARM_Curve* rhoCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam().GetCurve();
	ARM_Curve* thetaCurve=NULL;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		thetaCurve = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve();

	double nextt, expTeF02 = expTe * F0;
	expTeF02 *= expTeF02;
	double t=newSchedule[0],expt=exp(-mrs*(Te-t)/K_YEAR_LEN),expNextt;
    for(size_t i=0;i<schedSize;++i)
    {
		if (i < (schedSize - 1))
		{
			nextt		= newSchedule[i+1];
			expNextt	= exp(-mrs*(Te-nextt)/K_YEAR_LEN);
		}
		else
		{
			nextt		= Te;
			expNextt	= 1.0;
		}

		nu = nuCurve->Interpolate(nextt);
		
		itsFunctionDatas[i].itsVolOfVol_t = nu;

		vol = volCurve->Interpolate(nextt);
		itsFunctionDatas[i].itsVol_t = vol;		

		rho  = rhoCurve->Interpolate(nextt);

		itsFunctionDatas[i].itsRho_t = rho;

        nuVol = vol * nu;
		kappa = mrs - rho * nuVol / mrs;
		itsFunctionDatas[i].itsKappa_t = kappa;

		theta = thetaCurve->Interpolate(nextt);
		itsFunctionDatas[i].itsTheta_t = theta;

		/***
        /// Forward neutral probability existence test !
        if(kappa <= 0.5*nuVol*nuVol)
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                " : inconsistent HW1FSV model parameters (kappa must be superior to 0.5 *(vol*nu)^2)");
        }
		***/


		itsFunctionDatas[i].itsNu_2		= nu * nu;
		a = 0.5 * itsFunctionDatas[i].itsNu_2;
		itsFunctionDatas[i].itsA		= - a;
		itsFunctionDatas[i].itsKappaA	= kappa / a;

		///// Compute B factor
		itsFunctionDatas[i].itsTempB = rho * nu * expTe ;

		///// Compute C factor
		itsFunctionDatas[i].itsTempC = 0.5 * expTeF02;
		

		if (i == (schedSize - 1))
		{
			itsFunctionDatas[i].itsBeta				= 0.;
			itsFunctionDatas[i].itsBetaVolInst		= 0.;
			itsFunctionDatas[i].itsBetaVol			= 0.;
			itsFunctionDatas[i].itsKappaBetaVol		= 0.;
		}
		else
		{
			double beta = (expNextt-expt)/mrs;
			itsFunctionDatas[i].itsBeta			= beta;
			itsFunctionDatas[i].itsBetaVolInst	= vol * beta;

		}
		itsPrevValues1[i]=0.0;
		itsPrevValues2[i]=0.0;
		itsAngularShifts1[i]=0.0;
		itsAngularShifts2[i]=0.0;

		t=nextt;
		expt=expNextt;
    }
	
    for(i=schedSize-1;i> 0;--i)
    {
		itsFunctionDatas[i-1].itsBetaVol = itsFunctionDatas[i].itsBetaVol + itsFunctionDatas[i-1].itsBetaVolInst;
		itsFunctionDatas[i-1].itsKappaBetaVol = itsFunctionDatas[i-1].itsKappa_t * itsFunctionDatas[i-1].itsBetaVolInst;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HWSV1F::UpdateUDependentData(double F0,double Ux,double Uy,double expT0) const // exp (- mrs * T0)
{
	/// Store parameters for each interval ]ti, ti+1]
    /// They are stepwise right constant
    size_t schedSize = itsFunctionDatas.size();
    double deltaX, deltaY,mrs;
	mrs=GetMrs(); // ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion))).GetValueAtPoint(0);   

	double x,y;
	int i;
    for(i=0;i<schedSize;++i)
    {
		///// Compute B
		itsFunctionDatas[i].itsBx =  - itsFunctionDatas[i].itsTempB * ( F0 * Ux + expT0 /mrs) ;   
		itsFunctionDatas[i].itsBy =  - itsFunctionDatas[i].itsTempB *  F0 * Uy;    

		///// Compute C
		itsFunctionDatas[i].itsCx = itsFunctionDatas[i].itsTempC * (Ux - Ux * Ux + Uy * Uy);    
		itsFunctionDatas[i].itsCy = itsFunctionDatas[i].itsTempC * (Uy -2 * Ux * Uy);    
		
		///// Compute Discriminant Delta
		deltaX = itsFunctionDatas[i].itsBx * itsFunctionDatas[i].itsBx - itsFunctionDatas[i].itsBy * itsFunctionDatas[i].itsBy + 2 * itsFunctionDatas[i].itsNu_2 * itsFunctionDatas[i].itsCx;
		deltaY = 2 * (itsFunctionDatas[i].itsBx * itsFunctionDatas[i].itsBy + itsFunctionDatas[i].itsNu_2 * itsFunctionDatas[i].itsCy);

		ARM_Cplx::sqrt (deltaX, deltaY, itsFunctionDatas[i].itsDeltaX, itsFunctionDatas[i].itsDeltaY);

		//// Roots S1 & S2 of the Riccati Equation
		itsFunctionDatas[i].itsRoot_S2_X = (itsFunctionDatas[i].itsBx + itsFunctionDatas[i].itsDeltaX) / itsFunctionDatas[i].itsNu_2;
		itsFunctionDatas[i].itsRoot_S2_Y = (itsFunctionDatas[i].itsBy + itsFunctionDatas[i].itsDeltaY) / itsFunctionDatas[i].itsNu_2;


		itsFunctionDatas[i].itsBetaVolDeltaX = itsFunctionDatas[i].itsBetaVolInst * itsFunctionDatas[i].itsDeltaX;
		itsFunctionDatas[i].itsBetaVolDeltaY = itsFunctionDatas[i].itsBetaVolInst * itsFunctionDatas[i].itsDeltaY;

		//---- compute fact = 0.5 * Nu² / Delta = -A/Delta
		double fact_X, fact_Y;	
		ARM_Cplx::inverse(itsFunctionDatas[i].itsDeltaX, itsFunctionDatas[i].itsDeltaY, fact_X, fact_Y);
		x = -itsFunctionDatas[i].itsA;
		fact_X *= x;
		fact_Y *= x;
		itsFunctionDatas[i].itsFactX = fact_X;
		itsFunctionDatas[i].itsFactY = fact_Y;

		itsFunctionDatas[i].itsFactVolX = fact_X/itsFunctionDatas[i].itsVol_t;
		itsFunctionDatas[i].itsFactVolY = fact_Y/itsFunctionDatas[i].itsVol_t;
    }

	for(i=schedSize-1;i> 0;--i)
    {
		x = - itsFunctionDatas[i].itsBetaVolDeltaX;
		y = - itsFunctionDatas[i].itsBetaVolDeltaY;
		ARM_Cplx::exp(x,y,itsFunctionDatas[i].itsExpBetaVolDeltaX,itsFunctionDatas[i].itsExpBetaVolDeltaY);
		if(itsFunctionDatas[i].itsExpBetaVolDeltaX == 0.0)
		{
			/// Throw an error because further computation will inverse this null quantity !
			/// x=Re(BetaVolDelta) is too large to have an exp(-x) representable by a C++ double
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : exp(-integ{Re(delta).volFwdTe}) too low in HWSV1F model" );
		}

	}
	x = - itsFunctionDatas[0].itsBetaVolDeltaX;
	y = - itsFunctionDatas[0].itsBetaVolDeltaY;
	ARM_Cplx::exp(x,y,itsFunctionDatas[0].itsExpBetaVolDeltaX,itsFunctionDatas[0].itsExpBetaVolDeltaY);	
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HWSV1F::GenerateFunction(double expTe,
								  /// Return
								  double& PsiX,
								  double& PsiY,
								  ARM_GP_Vector& prevValues,
								  ARM_GP_Vector& angularShifts) const
{
    /// evalTime=0 implicitly

	double mrs=GetMrs();

    int i,iFirst=itsFunctionDatas.size() - 1;

    /// Compute constant term at expiry C(Te) such that beta2(Te,Te,u)=0
	double invS2_X, invS2_Y;	
	ARM_Cplx::inverse(itsFunctionDatas[iFirst].itsRoot_S2_X, itsFunctionDatas[iFirst].itsRoot_S2_Y, invS2_X, invS2_Y);
    itsFunctionDatas[iFirst<=1 ? 0 : iFirst-1].itsCstX =
        (itsFunctionDatas[iFirst].itsFactX - invS2_X) / itsFunctionDatas[iFirst].itsVol_t;
	itsFunctionDatas[iFirst<=1 ? 0 : iFirst-1].itsCstY =
        (itsFunctionDatas[iFirst].itsFactY - invS2_Y) / itsFunctionDatas[iFirst].itsVol_t;

    /// Propagate constant term by continuity w.r.t. time schedule
	double re_GRightWithoutExp, im_GRightWithoutExp;
    for(i=iFirst-1;i>0;--i)
    {		
		ARM_Cplx::product(itsFunctionDatas[i].itsCstX, itsFunctionDatas[i].itsCstY, itsFunctionDatas[i].itsExpBetaVolDeltaX, itsFunctionDatas[i].itsExpBetaVolDeltaY, re_GRightWithoutExp, im_GRightWithoutExp);
		re_GRightWithoutExp -= itsFunctionDatas[i].itsFactVolX;
		im_GRightWithoutExp -= itsFunctionDatas[i].itsFactVolY;
		ARM_Cplx::inverse (re_GRightWithoutExp, im_GRightWithoutExp, itsFunctionDatas[i-1].itsCstX, itsFunctionDatas[i-1].itsCstY);
		itsFunctionDatas[i-1].itsCstX += (itsFunctionDatas[i].itsRoot_S2_X * itsFunctionDatas[i].itsVol_t - itsFunctionDatas[i-1].itsRoot_S2_X * itsFunctionDatas[i-1].itsVol_t );
		itsFunctionDatas[i-1].itsCstY += (itsFunctionDatas[i].itsRoot_S2_Y * itsFunctionDatas[i].itsVol_t - itsFunctionDatas[i-1].itsRoot_S2_Y * itsFunctionDatas[i-1].itsVol_t );
		ARM_Cplx::inverse (itsFunctionDatas[i-1].itsCstX, itsFunctionDatas[i-1].itsCstY, itsFunctionDatas[i-1].itsCstX, itsFunctionDatas[i-1].itsCstY);
		itsFunctionDatas[i-1].itsCstX += itsFunctionDatas[i-1].itsFactVolX;
		itsFunctionDatas[i-1].itsCstY += itsFunctionDatas[i-1].itsFactVolY;
	}

	/// Compute G(evalTime)
	double GX, GY, Beta2X, Beta2Y;
	double tempX, tempY; 
	ARM_Cplx::product(itsFunctionDatas[0].itsExpBetaVolDeltaX,itsFunctionDatas[0].itsExpBetaVolDeltaY,itsFunctionDatas[0].itsCstX,itsFunctionDatas[0].itsCstY,tempX,tempY);
	GX =  expTe * (-itsFunctionDatas[0].itsFactVolX +  tempX);
	GY =  expTe * (-itsFunctionDatas[0].itsFactVolY +  tempY);
	
	/// Compute Beta2(evalTime)
	double invGX, invGY;
	ARM_Cplx::inverse(GX, GY, invGX, invGY);
	double temp = itsFunctionDatas[0].itsVol_t / expTe ;
	Beta2X = itsFunctionDatas[0].itsRoot_S2_X * temp + invGX;
	Beta2Y = itsFunctionDatas[0].itsRoot_S2_Y * temp + invGY;

//FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");

	/// Compute Gamma(evalTime). Takecare : the first time step may be > 0
    /// and we integrate from t=0 to t=Te
	double x,y,GammaX =0, GammaY=0, logX, logY, CexpX, CexpY;
	for(i=0;i<iFirst;++i)
	{
		GammaX += itsFunctionDatas[i].itsKappaBetaVol * itsFunctionDatas[i].itsRoot_S2_X;
		GammaY += itsFunctionDatas[i].itsKappaBetaVol * itsFunctionDatas[i].itsRoot_S2_Y;

		ARM_Cplx::product (itsFunctionDatas[i].itsCstX, itsFunctionDatas[i].itsCstY,itsFunctionDatas[i].itsExpBetaVolDeltaX,itsFunctionDatas[i].itsExpBetaVolDeltaY, CexpX, CexpY);
//fprintf(f,"Cexp=\t%15.10lf\t%15.10lf\n",CexpX,CexpY);

		logX = itsFunctionDatas[i].itsExpBetaVolDeltaX;
		logY = itsFunctionDatas[i].itsExpBetaVolDeltaY;
//fprintf(f,"log1=\t%15.10lf\t%15.10lf\n",logX,logY);

		x = -itsFunctionDatas[i].itsFactVolX;
		y = -itsFunctionDatas[i].itsFactVolY;
		ARM_Cplx::product(logX, logY, x, y, logX, logY);
//fprintf(f,"log2=\t%15.10lf\t%15.10lf\n",logX,logY);

		logX += CexpX;
		logY += CexpY;
//fprintf(f,"log3=\t%15.10lf\t%15.10lf\n",logX,logY);

		x = CexpX - itsFunctionDatas[i].itsFactVolX;
		y = CexpY - itsFunctionDatas[i].itsFactVolY;
		ARM_Cplx::quotient (logX, logY, x, y, logX, logY);
//fprintf(f,"log4=\t%15.10lf\t%15.10lf\n",logX,logY);

		ARM_Cplx::log(logX, logY, logX, logY);
		if(ARM_NumericConstants::ARM_PI - ARM_NumericConstants::ARM_PI_BY_4 < prevValues[i] &&
			prevValues[i] < ARM_NumericConstants::ARM_PI &&
			-ARM_NumericConstants::ARM_PI < logY &&
			logY < -ARM_NumericConstants::ARM_PI + ARM_NumericConstants::ARM_PI_BY_4)
		{
			/// Keep record of number of turns
			angularShifts[i] += ARM_NumericConstants::ARM_2_PI;
		}
		else if(-ARM_NumericConstants::ARM_PI < prevValues[i] &&
			prevValues[i] < -ARM_NumericConstants::ARM_PI + ARM_NumericConstants::ARM_PI_BY_4 &&
			ARM_NumericConstants::ARM_PI - ARM_NumericConstants::ARM_PI_BY_4 < logY &&
			logY < ARM_NumericConstants::ARM_PI)
		{
			/// Keep record of number of turns
			angularShifts[i] -= ARM_NumericConstants::ARM_2_PI;
		}

		/// Save new angular part...
		prevValues[i] = logY;

		///... and add shift to maintain continuity
		logY += angularShifts[i];


//fprintf(f,"log5\t%15.10lf\t%15.10lf\n",logX,logY);

		GammaX += itsFunctionDatas[i].itsKappaA * logX;
		GammaY += itsFunctionDatas[i].itsKappaA * logY;
	}
   
	/// Compute charateristic function : psi(t,u)=exp[Gamma(t,u) + Beta2(t,u)] at evalTime
	double logPsiX,logPsiY;
	logPsiX = GammaX + Beta2X;
	logPsiY = GammaY + Beta2Y;
	double tmp = exp(logPsiX);
	PsiX = tmp * cos(logPsiY);
	PsiY = tmp * sin(logPsiY);

/****
fprintf(f,"%15.10lf\t%15.10lf\t%15.10lf\t%15.10lf\n",Beta2X,Beta2Y,GammaX,GammaY);
fclose(f);
****/
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	const ARM_GP_Vector* const timeSteps = GetNumMethod()->GetTimeSteps();

	// Compute the deterministic part
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	double startTime		= (*timeSteps)[timeIdx];
	double endTime			= startTime + dt;

	double zcStart	= ZcCurve->DiscountPrice(startTime/K_YEAR_LEN);
    double zcEnd	= ZcCurve->DiscountPrice(endTime/K_YEAR_LEN);
	double discount = zcEnd/zcStart;

	size_t statesSize		= states->size();
	ARM_GP_Vector* result	= new ARM_GP_Vector(statesSize,0.0);
	size_t modelNb			= GetModelNb();

    /// Simple mapping version  : r(t) = f(0,t) + X(t)
    double Xt;
    dt /= K_YEAR_LEN;
	for( size_t i=0; i<statesSize; ++i )
	{
		Xt = states->GetModelState(i,modelNb+X_VARIABLE);
		(*result)[i] = discount * exp(-Xt*dt);
    }

	return ARM_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
	// Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,zcT) );
    }

    // Volatility computation (ARM_ModelParamsMSV1F class is pure virtual)
    int i,nbStates=states->size();
	if(evalTime >= maturityTime)
		return ARM_VectorPtr(new ARM_GP_Vector(nbStates,1.0));

    double BetatT	= ARM_ModelParamsHW1F::BetatT(GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion),evalTime,maturityTime);
    double BetatT_2 = BetatT*BetatT;
    double zct		= ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);

	if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc  ||
	   GetNumeraire()->GetType() == ARM_Numeraire::TerminalEventZc)
	{
		double numTime = GetNumeraire()->GetMaturity();
		double BetatTStar = BetatT;
		if(numTime != maturityTime)
			BetatTStar = ARM_ModelParamsHW1F::BetatT(GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion),evalTime,numTime);
		BetatT_2 *= 1.0 - 2*BetatTStar/BetatT;

	}
    else if( !(GetNumeraire() == ARM_NumerairePtr(NULL)
			|| GetNumeraire()->GetType() == ARM_Numeraire::Cash) )
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            " : only Cash, Terminal Zc or Terminal Event Zc numeraires are supported by HWSV1F model");
	}



	size_t modelNb = GetModelNb();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
	double temp = 0.;
    for(i=0;i<nbStates;++i)
    {
		(*values)[i]=zcT/zct*exp( -BetatT * states->GetModelState(i,modelNb+X_VARIABLE)
                        - 0.5 * BetatT_2 * states->GetModelState(i,modelNb+PHI_VARIABLE) ); 
	}
    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeCapletF0
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeCapletF0(	double evalTime,
									double fwdStartTime, 
									double fwdEndTime,
									double fwdPeriod,
									const ARM_GP_Vector& strikes,
									const ARM_PricingStatesPtr& states,
									ARM_GP_Vector& newStrikes) const
{	
	size_t stateIdx,nbStates = states->size();

	double mrs = GetMrs();
	double F00 = ( exp(-mrs *(fwdEndTime - evalTime)/K_YEAR_LEN) - exp(-mrs *(fwdStartTime - evalTime)/K_YEAR_LEN) )/mrs;
	ARM_VectorPtr F0(new ARM_GP_Vector(nbStates,F00));

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdStartTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdEndTime,states);

	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		newStrikes[stateIdx]	= (1.0 + fwdPeriod * strikes[stateIdx]) * (*zcEnd)[stateIdx] / (*zcStart)[stateIdx];
		newStrikes[stateIdx]	= 1.0 / newStrikes[stateIdx];
	}

	return F0;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::VanillaCaplet(
		const string& curveName, 
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
		const ARM_PricingStatesPtr& states) const
{


	/// Handle the case of dummy states
	if(states == ARM_PricingStatesPtr(NULL) && evalTime > K_NEW_DOUBLE_TOL)
// FIXMEFRED: mig.vc8 (30/05/2007 16:16:37):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/// If spot evaluation reset the pricing state
	ARM_PricingStatesPtr newStates(states);
	if(evalTime < K_NEW_DOUBLE_TOL)
		newStates = FirstPricingStates(1);
	size_t nbStates = newStates->size();

	if(strikesPerState.size() < 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike is missing in HWSV1F caplet pricing" );

	/// No payment lag adjustement computation
	if(fabs(payTime-fwdEndTime) > 5)
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            " : convexity adjustment not available for HWSV1F model");
	}

	ARM_VectorPtr values;

#ifdef PRICING_NB_TIMES
ARM_Timer timer;
timer.ClockStartTime();

for(size_t tIdx=0;tIdx<PRICING_NB_TIMES;++tIdx)
{
#endif

	ARM_GP_Vector newStrikes(nbStates);
	ARM_VectorPtr F0(ComputeCapletF0(evalTime,fwdStartTime,fwdEndTime,fwdPeriod,strikesPerState,newStates,newStrikes));

	if(nbStates==1 &&  newStrikes[0] <= 0.0)
	{
		if(capFloor==K_FLOOR)
		{
			values = ARM_VectorPtr(new ARM_GP_Vector(1, 0.0));
		}
		else
		{
			ARM_VectorPtr zcStart = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,newStates);
			values = ARM_VectorPtr(new ARM_GP_Vector(1, ((*zcStart)[0]/newStrikes[0])*(newStrikes[0]-1)*payNotional*period/fwdPeriod));
		}
	}
	else
	{
		if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::VolMeanReversion))
		{
			/// Vol MRS is input
			if(itsNumericals.GetMaxDecay()>0)
				/// Exponential terms are sampled then Riccati systems have analytical solutions
				values = ComputeAnalyticalOptionPrice(evalTime,F0,fwdStartTime,fwdResetTime,newStrikes,capFloor,payNotional*period/fwdPeriod,newStates);
			else
				/// Riccatti system are solved using the numerical Runge-Kutta method
				values = ComputeRungeKuttaOptionPrice(evalTime,F0,fwdStartTime,fwdResetTime,newStrikes,capFloor,payNotional*period/fwdPeriod,newStates);
		}
		else if(nbStates==1)
			/// Vol MRS is linked to other model parameters => no time dependency of Riccatti coefficients,
			/// Riccati systems have analytical solutions
			values = ComputeOptionPrice((*F0)[0],fwdStartTime,fwdResetTime,newStrikes[0],capFloor,payNotional*period/fwdPeriod,newStates);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
				" : linked vol MRS version not supported for non asOf evaluation in HWSV1F model");
	}

#ifdef PRICING_NB_TIMES
}
timer.ClockEndTime();
FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
fprintf(f,"Duration = %10.5lf ms\n",timer.GetDuration()*1000.0/PRICING_NB_TIMES);
fclose(f);
#endif

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeSwaptionF0
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeSwaptionF0(	double evalTime,
								double startTime, 
								double endTime,
								const ARM_GP_Vector& fixPayTimes,
								const ARM_GP_Vector& fixPayPeriods,
								int callPut,
								const ARM_GP_Matrix& strikes,
								bool isConstantNotional,
								const ARM_GP_Vector& fixNotional,
								const ARM_GP_Vector& floatNotional,
								size_t refNotionalIdx,
								const ARM_PricingStatesPtr& states,
								ARM_GP_Vector& newStrikes) const
{
	size_t stateIdx,nbStates = states->size();

	double mrs = GetMrs();
	
	double expMrss = exp(-mrs *(startTime - evalTime)/K_YEAR_LEN);
	double expMrse = exp(-mrs *(endTime - evalTime)/K_YEAR_LEN);

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,startTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,endTime,states);

	size_t i,nbPeriods = fixPayTimes.size();

	if(endTime != fixPayTimes[nbPeriods - 1])
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : end time mismatch in bond strip vol computation");

	vector<ARM_VectorPtr> zcFixPay(nbPeriods);
	ARM_GP_Vector expMrs(nbPeriods);
	for(i = 0; i<nbPeriods; ++i)
	{
		zcFixPay[i] = GetDiscountFunctor()->DiscountFactor("",evalTime,fixPayTimes[i],states);
		expMrs[i] = exp(-mrs*(fixPayTimes[i]-evalTime)/K_YEAR_LEN);
	}

	ARM_VectorPtr F0(new ARM_GP_Vector(nbStates,0.0));

	double flow,newStrike,F00;
	if(isConstantNotional)
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx)
		{
			newStrike=0.0,F00=0.0;
			for(i=0;i+1<nbPeriods;++i)
			{
				flow		= strikes(stateIdx,i) * fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F00			+= flow * (expMrs[i] - expMrss);
			}
			flow		= (1.0 + strikes(stateIdx,nbPeriods-1) * fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			newStrike	= (newStrike + flow) / (*zcStart)[stateIdx];
			F00			= (F00 + flow * (expMrse - expMrss)) / (*zcStart)[stateIdx];

			newStrike		= 1.0 / newStrike;
			F00	*= (newStrike/mrs);

			(*F0)[stateIdx]	= F00;
			newStrikes[stateIdx] = newStrike;
		}
	}
	else
	{
		/// Case where first fixed flows have non null notional but floating leg ones are null
		for(stateIdx=0;stateIdx<nbStates;++stateIdx)
		{
			newStrike=0.0,F00=0.0;
			for(i=0;i<refNotionalIdx;++i)
			{
				if(fabs(fixNotional[i]) > K_NEW_DOUBLE_TOL)
				{
					flow		= strikes(stateIdx,i) * fixNotional[i] *
								  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
					newStrike	+= flow;
					F00			+= flow * (expMrs[i] - expMrss);
				}
			}
			for(i=refNotionalIdx;i+1<nbPeriods;++i)
			{
				flow		= floatNotional[i] - floatNotional[i+1] + strikes(stateIdx,i) * fixNotional[i] *
							  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F00			+= flow * (expMrs[i] - expMrss);
			}
			flow		= (floatNotional[nbPeriods-1] + strikes(stateIdx,nbPeriods-1) * fixNotional[nbPeriods-1] *
						  fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			double coefStart = (*zcStart)[stateIdx] * floatNotional[refNotionalIdx];
			newStrike	= (newStrike + flow)/coefStart;
			F00			= (F00 + flow * (expMrse - expMrss))/coefStart;

			newStrike	= 1.0 / newStrike;
			F00	*= (newStrike/mrs);

			(*F0)[stateIdx]	= F00;
			newStrikes[stateIdx] = newStrike;
		}
	}

	return F0;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: UpdateStdHWModel
///	Returns: void
///	Action : Update the deterministic volatility
///			 reference H&W model
////////////////////////////////////////////////////
void ARM_HWSV1F::UpdateStdHWModel(double volFactor) const
{
	ARM_CurveModelParam  paramvol = GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam();
	for(size_t i=0;i<paramvol.size();++i)
		paramvol.SetValueAtPoint(i,paramvol.GetValueAtPoint(i)*volFactor);
	ARM_CurveModelParam  paramMRS = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam();

	ARM_ModelParamVector paramVector(2);
	paramVector[0] = &paramvol;
	paramVector[1] = &paramMRS;
	GetAnalyticalModel()->SetModelParams( ARM_ModelParamsHW1FStd(paramVector) );
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
///			 If the swaption is standard, this method
///          calls ARM_HullWhite::VanillaSwaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::VanillaSwaption(
				const string& curveName,
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
				bool isConstantNotional,
				bool isConstantSpread,
				bool isConstantStrike) const
{
	/// Handle the case of dummy states
	if(states == ARM_PricingStatesPtr(NULL) && evalTime > K_NEW_DOUBLE_TOL)
// FIXMEFRED: mig.vc8 (30/05/2007 16:17:21):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/// If spot evaluation reset the pricing state
	ARM_PricingStatesPtr newStates(states);
	if(evalTime < K_NEW_DOUBLE_TOL)
		newStates = FirstPricingStates(1);
	size_t i,nbStates = newStates->size();

	/// First row is always used
	if(strikesPerState.cols() < 1 || strikesPerState.rows() != nbStates)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike profile is missing in HWSV1F swaption pricing" );

	/// Float & fixed leg ends are equal
	if(fabs(floatEndTime - fixPayTimes[fixPayTimes.size()-1])>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end dates are not matching in HWSV1F swaption pricing" );

	size_t refNotionalIdx=0;
	if(isConstantNotional)
	{
		/// Same notional on float & fixed legs is constant
		if(floatNotional[0] != fixNotional[0] || fabs(floatNotional[0])<K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, " : same non null notional is required on float & fixed legs if bullet in HWSV1F swaption pricing" );
	}
	else
	{
		/// Same notional profile size on float & fixed legs
		if(floatNotional.size() != fixNotional.size())
			ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fix legs are required to be of same frequency in HWSV1F swaption pricing" );

		/// Locate first non null floating leg notional
		for(refNotionalIdx=0;refNotionalIdx<floatNotional.size();++refNotionalIdx)
		{
			if(fabs(floatNotional[refNotionalIdx]) >= K_NEW_DOUBLE_TOL)
				break;
		}
		if(refNotionalIdx >= floatNotional.size()) 
			ARM_THROW( ERR_INVALID_ARGUMENT, " : can't locate a non null notional on floating leg in HWSV1F swaption pricing" );
	}
	double refNotional = floatNotional[refNotionalIdx];

	ARM_VectorPtr values;

#ifdef PRICING_NB_TIMES
ARM_Timer timer;
timer.ClockStartTime();

for(size_t tIdx=0;tIdx<PRICING_NB_TIMES;++tIdx)
{
#endif

	/// Start time = start of the 1st period where floating leg notional is not null
	double startTime = (refNotionalIdx==0 ? floatStartTime : fixPayTimes[refNotionalIdx-1]);

	ARM_GP_Vector newStrikes(nbStates);
	ARM_VectorPtr F0(ComputeSwaptionF0(evalTime,startTime,floatEndTime,fixPayTimes,fixPayPeriods,callPut,
		strikesPerState,isConstantNotional,fixNotional,floatNotional,refNotionalIdx,newStates,newStrikes));

	if(nbStates == 1 && newStrikes[0] <= 0.0)
	{
		if(callPut==K_PUT)
		{
			values = ARM_VectorPtr(new ARM_GP_Vector(1, 0.0));
		}
		else
		{
			ARM_VectorPtr zcStart = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,startTime,newStates);
			values = ARM_VectorPtr(new ARM_GP_Vector(1, ((*zcStart)[0]/newStrikes[0])*(newStrikes[0]-1)*refNotional));
		}
	}
	else
	{
		if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::VolMeanReversion))
		{
			/// Vol MRS is input
			if(itsNumericals.GetMaxDecay()>0)
				/// Exponential terms are sampled then Riccati systems have analytical solutions
				values = ComputeAnalyticalOptionPrice(evalTime,F0,startTime,swapResetTime,newStrikes,callPut,refNotional,newStates);
			else
				/// Riccatti system are solved using the numerical Runge-Kutta method
				values = ComputeRungeKuttaOptionPrice(evalTime,F0,startTime,swapResetTime,newStrikes,callPut,refNotional,newStates);
		}
		else if(nbStates==1)
			/// Vol MRS is linked to other model parameters => no time dependency of Riccatti coefficients,
			/// Riccati systems have analytical solutions
			values = ComputeOptionPrice((*F0)[0],startTime,swapResetTime,newStrikes[0],callPut,refNotional,newStates);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
				" : linked vol MRS version not supported for non asOf evaluation in HWSV1F model");

/***** Formula correction *****/

		/// Refresh linked standard H&W model for deterministic volatility reference price
		/// For deterministic correction, V(t=0) is assumed to be equal to Long Term Vol at t=0
		double LongTermVol0 = 1.0;
		if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		{
			double LongTermVar0 = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
			LongTermVol0 = sqrt(LongTermVar0);
		}
		UpdateStdHWModel(LongTermVol0);

		/// Reference price for correction using standard Hull & White 1F model
		/// Variable notional is internally available
		ARM_VectorPtr refPrices = (dynamic_cast<ARM_HullWhite*> (&*GetAnalyticalModel()))->VanillaSwaption(curveName,
					evalTime,
					swapResetTime,
					fixNotional,
					floatNotional,
					floatStartTime,
					floatEndTime,
					floatResetTimes,
					floatStartTimes,
					floatEndTimes,
					floatIntTerms,
					fixPayTimes,
					fixPayPeriods,
					strikesPerState,
					callPut,
					newStates,
					isConstantNotional,
					isConstantSpread,
					isConstantStrike);

/****/			
		/// Compute standard H&W1F price using same approximation used in stochastic volatility case
		double stdDev=sqrt( ((const ARM_ModelParamsHW1FStd* const)GetAnalyticalModel()->GetModelParams())->StateLocalVariance(evalTime,swapResetTime,0.0) );

		ARM_VectorPtr zcStart = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,startTime,newStates);

		/// Payer swaption => Put on bond, Receiver swaption => Call on bond
		int optType = K_CALL;
		if(callPut==K_CALL)
			optType = K_PUT;
		double F00,approxPrice,stdDev0;
		for(i=0;i<nbStates;++i)
		{
			F00 = (*F0)[i];
			stdDev0 = fabs(F00)*stdDev;
			approxPrice = refNotional/newStrikes[i] * BlackSholes_Formula(1.0,stdDev0,(*zcStart)[i],newStrikes[i],optType);

			/// Correct closed form formula with H&W1F std prices
			(*values)[i] = (*values)[i] + (*refPrices)[i] - approxPrice;
			if((*values)[i] < 0.0)
			{
				if((*values)[i] < -refNotional*0.0001) // 1bp
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, " : price is negative and < -1bp !" );
				}
				else
					(*values)[i]=0.0;
			}
		}
/****/
/****
		/// Compute standard H&W1F price using same approximation used in stochastic volatility case
		double stdDev=fabs((*F0)[0])*sqrt( ((const ARM_ModelParamsHW1FStd* const)GetAnalyticalModel()->GetModelParams())->StateLocalVariance(evalTime,swapResetTime,0.0) );
		ARM_VectorPtr zcStart = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,startTime,newStates);
		ARM_VectorPtr approxPrices(new ARM_GP_Vector(nbStates,0.0));
		if(callPut==K_CALL)
		{
			/// Payer swaption => Put on bond
			for(i=0;i<nbStates;++i)
				(*approxPrices)[i] = refNotional/newStrikes[i] * BlackSholes_Formula(1.0,stdDev,(*zcStart)[i],newStrikes[i],K_PUT);
		}
		else
		{
			/// Receiver swaption => Call on bond
			for(i=0;i<nbStates;++i)
				(*approxPrices)[i] = refNotional/newStrikes[i] * BlackSholes_Formula(1.0,stdDev,(*zcStart)[i],newStrikes[i],K_CALL);
		}

		/// Correct closed form formula with H&W1F std prices
		for(i=0;i<nbStates;++i)
		{
			(*values)[i] = (*values)[i] + (*refPrices)[i] - (*approxPrices)[i];
			if((*values)[i] < 0.0)
			{
				if((*values)[i] < -refNotional*0.0001) // 1bp
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, " : price is negative and < -1bp !" );
				}
				else
					(*values)[i]=0.0;
			}
		}
****/
//***** Formula correction *****
	}


#ifdef PRICING_NB_TIMES
}
timer.ClockEndTime();
FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
fprintf(f,"Duration = %10.5lf ms\n",timer.GetDuration()*1000.0/PRICING_NB_TIMES);
fclose(f);
#endif

	return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeOptionPrice(
											double F0,
											double T0,
											double Te,
											double newStrike,
											int	RecPay,
											double payNotional,
											const ARM_PricingStatesPtr& states) const
{
	
	double PsiX, PsiY, t ;
	
	///// Constant Param
	double mrs=GetMrs();
	double yfT0 = T0/K_YEAR_LEN;
	double yfTe = Te/K_YEAR_LEN;
	double expT0 = exp (- mrs * yfT0);
	double expTe = exp ( mrs * yfTe);
	double lnK = log(newStrike);
    double cos_, sin_;
	double	Integral_1 = 0.0, Integral_2 = 0.0;
    int i;

	double t_min = 1.0e-8;
	InitFunctionData(F0,Te,expTe);

    double abslnK=fabs(lnK);
	double oscilSpeedLimit = ARM_NumericConstants::ARM_2_PI/itsNumericals.GetFormulaParam(HWSVNumericals::LimitStep);
	double oscilSpeed = (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
	bool isSpeedLimit = (oscilSpeed == oscilSpeed);
    double oscilSize = ARM_NumericConstants::ARM_2_PI/oscilSpeed;

	itsNumericals.SetOscilTrace(HWSVNumericals::NbPoints,0);
	itsNumericals.ResetNbStepsTrace();


    /// First oscillation integration
    size_t nbPts = (isSpeedLimit ? itsNumericals.GetFormulaParam(HWSVNumericals::FirstLimitNbSteps)
								: itsNumericals.GetFormulaParam(HWSVNumericals::FirstNbSteps));
    GaussLegendre_Coefficients GL1(nbPts);
    double tmin = 0.0;
    double tmax = oscilSize;
    double scalet=0.5*oscilSize;

    double t0=tmin+scalet;
    double scalew,localIntegral_1=0.0,localIntegral_2=0.0;
//	f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
    for(i=0;i<nbPts;++i)
    {
        t = t0 + GL1.get_point(i) * scalet;
        scalew = GL1.get_weight(i) * scalet/t;
		cos_ = cos(-t*lnK);
		sin_ = sin(-t*lnK);

		UpdateUDependentData(F0,1.0,t,expT0) ;
		GenerateFunction(expTe,PsiX,PsiY,itsPrevValues1,itsAngularShifts1); 
		localIntegral_1 +=  (PsiY * cos_ + PsiX * sin_) * scalew;
//		fprintf(f,"t=\t%15.5lf\tPsi=\t%15.10lf\t%15.10lf\t",t,PsiX,PsiY);

		UpdateUDependentData(F0,0.,t,expT0) ;
		GenerateFunction(expTe,PsiX,PsiY,itsPrevValues2,itsAngularShifts2); 
		localIntegral_2 += (PsiY * cos_ + PsiX * sin_) * scalew;
//		fprintf(f,"%15.10lf\t%15.10lf\n",PsiX,PsiY);
    }
    Integral_1 += localIntegral_1;
    Integral_2 += localIntegral_2;
//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

	itsNumericals.AddOscilTraceNbPoints(nbPts);

    /// Next oscillation integrations
    tmin = tmax;
	nbPts = (isSpeedLimit ? itsNumericals.GetFormulaParam(HWSVNumericals::NextLimitNbSteps)
						 : itsNumericals.GetFormulaParam(HWSVNumericals::NextNbSteps));
    GaussLegendre_Coefficients GL(nbPts);
    size_t iterIdx=1;
	while( iterIdx<OSCILLATION_NBMAX &&
		  (fabs(localIntegral_1)>DEFAULT_PRECISION || fabs(localIntegral_2)>DEFAULT_PRECISION) )
	{
		tmax += oscilSize;
		t0=tmin+scalet;
		localIntegral_1=0.0;
		localIntegral_2=0.0;

		for(i=0;i<nbPts;++i)
		{
			t = t0 + GL.get_point(i) * scalet;
			scalew = GL.get_weight(i) * scalet/t;
			cos_ = cos(-t*lnK);
			sin_ = sin(-t*lnK);

			UpdateUDependentData(F0,1.0,t,expT0) ;
			GenerateFunction(expTe,PsiX,PsiY,itsPrevValues1,itsAngularShifts1); 
			localIntegral_1 += (PsiY * cos_ + PsiX * sin_) * scalew;
//			fprintf(f,"t=\t%15.5lf\tPsi=\t%15.10lf\t%15.10lf\t",t,PsiX,PsiY);

			UpdateUDependentData(F0,0.,t,expT0) ;
			GenerateFunction(expTe,PsiX,PsiY,itsPrevValues2,itsAngularShifts2); 
			localIntegral_2 += (PsiY * cos_ + PsiX * sin_) * scalew;
//			fprintf(f,"%15.10lf\t%15.10lf\n",PsiX,PsiY);
		}

		Integral_1 += localIntegral_1;
		Integral_2 += localIntegral_2;
		tmin = tmax;
		++iterIdx;
//		fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

		itsNumericals.AddOscilTraceNbPoints(nbPts);
	}

//	fclose(f);

	itsNumericals.SetOscilTrace(HWSVNumericals::FirstPeriod,oscilSize);
	itsNumericals.SetOscilTrace(HWSVNumericals::NextPeriod,oscilSize);
	itsNumericals.SetOscilTrace(HWSVNumericals::NbPeriods,iterIdx);

    if(iterIdx >= OSCILLATION_NBMAX)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            " : max iteration reached in option pricing with oscillating integral");
    }


    /// Compute option price
	double pi_1   = 0.5 + Integral_1 / ARM_NumericConstants::ARM_PI;
	double pi_2	  = 0.5 + Integral_2 / ARM_NumericConstants::ARM_PI;
	double expect = pi_1 - newStrike * pi_2 ;
	
	if (RecPay == 1) expect = expect - 1 + newStrike;

	double zcStart	= GetZeroCurve()->DiscountPrice(yfT0);
    ARM_VectorPtr values(new ARM_GP_Vector(1, (zcStart/newStrike)*expect*payNotional));


	/// Reset turns number for further use
	for(i=0;i<itsFunctionDatas.size();++i)
	{
		itsPrevValues1[i]=0.0;
		itsPrevValues2[i]=0.0;
		itsAngularShifts1[i]=0.0;
		itsAngularShifts2[i]=0.0;
	}

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeRungeKuttaOptionPrice(	double evalTime,
														const ARM_VectorPtr& F0,
														double T0,
														double Te,
														const ARM_GP_Vector& newStrikes,
														int	RecPay,
														double payNotional,
														const ARM_PricingStatesPtr& states) const
{
	size_t stateIdx,nbStates = states->size();
	size_t modelIdx	= GetModelNb();
	ARM_VectorPtr values(new ARM_GP_Vector(nbStates,0.0));

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,T0,states);

	double yfT0 = T0/K_YEAR_LEN;
	double mrs=GetMrs();
	double expT0 = exp (- mrs * yfT0);

	/// Set non SO pricing in numerical stuff
	bool isSOFormula = false;
	itsNumericals.SetIsSOFormula(isSOFormula);

	const ARM_GP_Vector& numericals = itsNumericals.GetFormulaParams();

	size_t fctMultiplier = 4 * (itsNumericals.IsStdFormula() ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);

	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	ARM_GP_Vector ImU,PsiX1,PsiY1,PsiX0,PsiY0;
	ImU.reserve(maxNbPts);
	PsiX1.reserve(maxNbPts);
	PsiY1.reserve(maxNbPts);
	PsiX0.reserve(maxNbPts);
	PsiY0.reserve(maxNbPts);

	size_t i,nbPts;
	int nbStepsIncr;
	vector<ARM_GP_Vector> solverVars(RK5_NBTMP_VAR > RK4_NBTMP_VAR ? RK5_NBTMP_VAR : RK4_NBTMP_VAR);
	size_t nbFcts = maxNbPts * fctMultiplier;
	for(i=0;i<solverVars.size();++i)
		solverVars[i].reserve(nbFcts);

	double var0,K,F00,zcStart0,lnK,abslnK;
	double Integral_1,Integral_2;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,localIntegral_1,localIntegral_2=0.0;
	double pi_1,pi_2,expect,pi;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		var0		= states->GetModelState(stateIdx,modelIdx+V_VARIABLE);
		F00			= (*F0)[stateIdx];
		zcStart0	= (*zcStart)[stateIdx];

		K = newStrikes[stateIdx];
		lnK = log(K);
		abslnK=fabs(lnK);

		Integral_1 = 0.0;
		Integral_2 = 0.0;

		InitSystemData(evalTime,F00,mrs,Te,expT0,maxNbPts,itsNumericals.IsStdFormula(),isSOFormula);

		/// 1st oscillation or 1st interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isFirstSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		firstOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// 2nd oscillation or 2nd interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::NextLimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isNextSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		nextOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// Next oscillations or next intervals of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LastLimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isLastSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		lastOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;


		itsNumericals.SetOscilTrace(HWSVNumericals::NbPoints,0);
		itsNumericals.ResetNbStepsTrace();

		/// Initialise a RK solver and its internal variable set (for optimisation purpose)
		ARM_RiccatiHWSV1F riccatiSystem(itsSystemDatas,itsNumericals.GetSolverType(),itsNumericals.IsStdFormula());


	//FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");

	//for(size_t ii=0;ii<itsSystemDatas.size();++ii)
	//	fprintf(f,"(%8.3lf,%8.3lf)\t",itsSystemDatas[ii].itsYf,itsSystemDatas[ii].itsTime);
	//fprintf(f,"\n");
	//fclose(f);

		/// First oscillation integration
		if(itsNumericals.IsStdFormula() && isFirstSpeedLimit)
			nbPts = numericals[HWSVNumericals::FirstLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::FirstNbSteps];
		GaussLegendre_Coefficients GL1(nbPts);
		tmin = 0.0;
		tmax = firstOscilSize;
		scalet=0.5*firstOscilSize;

		ARM_IntVector nbSteps;
		t0=tmin+scalet;
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		PsiX0.resize(nbPts);
		PsiY0.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL1.get_point(i) * scalet;

		if(itsNumericals.GetSolverType()==ARM_ODEFunc::RK4Constant)
		{
			nbStepsIncr=static_cast<int>(floor(itsNumericals.GetSolverParam(HWSVNumericals::RK4FactorTime)*(ImU[nbPts-1]-ImU[0])));
			itsNumericals.SetSolverParam(HWSVNumericals::RK4NbStepsPerYear,itsNumericals.GetSolverParam(HWSVNumericals::RK4NbStepsPerYear) + nbStepsIncr);
		}

		UpdateUDependentSystemData(ImU,itsNumericals.IsStdFormula(),isSOFormula);

		itsNumericals.SolveSystem(riccatiSystem,solverVars,PsiX1,PsiY1,PsiX0,PsiY0,var0);

		itsNumericals.IntegrateSystem(GL1,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,
			localIntegral_1,localIntegral_2);

		Integral_1 += localIntegral_1;
		Integral_2 += localIntegral_2;

	//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);

		itsNumericals.AddOscilTraceNbPoints(nbPts);
		nbSteps=riccatiSystem.GetNbSteps();
		nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
		itsNumericals.PushNbStepsTrace(nbSteps);


		/// Second oscillation integrations
		tmin = tmax;
		tmax += nextOscilSize;
		scalet=0.5*nextOscilSize;
		if(itsNumericals.IsStdFormula() && isNextSpeedLimit)
			nbPts = numericals[HWSVNumericals::NextLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::NextNbSteps];
		GaussLegendre_Coefficients GL2(nbPts);
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		PsiX0.resize(nbPts);
		PsiY0.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		t0=tmin+scalet;

		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL2.get_point(i) * scalet;

		if(itsNumericals.GetSolverType()==ARM_ODEFunc::RK4Constant)
		{
			nbStepsIncr=static_cast<int>(floor(itsNumericals.GetSolverParam(HWSVNumericals::RK4FactorTime)*(ImU[nbPts-1]-ImU[0])));
			itsNumericals.SetSolverParam(HWSVNumericals::RK4NbStepsPerYear,itsNumericals.GetSolverParam(HWSVNumericals::RK4NbStepsPerYear) + nbStepsIncr);
		}

		UpdateUDependentSystemData(ImU,itsNumericals.IsStdFormula(),isSOFormula);

		itsNumericals.SolveSystem(riccatiSystem,solverVars,PsiX1,PsiY1,PsiX0,PsiY0,var0);

		itsNumericals.IntegrateSystem(GL2,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,
			localIntegral_1,localIntegral_2);

		Integral_1 += localIntegral_1;
		Integral_2 += localIntegral_2;

	//		fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);

		itsNumericals.AddOscilTraceNbPoints(nbPts);
		nbSteps=riccatiSystem.GetNbSteps();
		nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
		itsNumericals.PushNbStepsTrace(nbSteps);


		/// Next oscillation integrations
		size_t iterIdx=2;
		tmin = tmax;
		scalet=0.5*lastOscilSize;
		if(itsNumericals.IsStdFormula() && isLastSpeedLimit)
			nbPts = numericals[HWSVNumericals::LastLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::LastNbSteps];
		GaussLegendre_Coefficients GL(nbPts);
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		PsiX0.resize(nbPts);
		PsiY0.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		while( iterIdx<OSCILLATION_NBMAX &&
			  (fabs(localIntegral_1) > numericals[HWSVNumericals::IntegrationPrecision] ||
			   fabs(localIntegral_2) > numericals[HWSVNumericals::IntegrationPrecision]) )
		{
			tmax += lastOscilSize;
			t0=tmin+scalet;

			for(i=0;i<nbPts;++i)
				ImU[i] = t0 + GL.get_point(i) * scalet;

			if(itsNumericals.GetSolverType()==ARM_ODEFunc::RK4Constant)
			{
				nbStepsIncr=static_cast<int>(floor(itsNumericals.GetSolverParam(HWSVNumericals::RK4FactorTime)*(ImU[nbPts-1]-ImU[0])));
				itsNumericals.SetSolverParam(HWSVNumericals::RK4NbStepsPerYear,itsNumericals.GetSolverParam(HWSVNumericals::RK4NbStepsPerYear) + nbStepsIncr);
			}

			UpdateUDependentSystemData(ImU,itsNumericals.IsStdFormula(),isSOFormula);

			itsNumericals.SolveSystem(riccatiSystem,solverVars,PsiX1,PsiY1,PsiX0,PsiY0,var0);

			itsNumericals.IntegrateSystem(GL,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,
				localIntegral_1,localIntegral_2);

			Integral_1 += localIntegral_1;
			Integral_2 += localIntegral_2;

	//		fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);

			itsNumericals.AddOscilTraceNbPoints(nbPts);
			nbSteps=riccatiSystem.GetNbSteps();
			nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
			itsNumericals.PushNbStepsTrace(nbSteps);

			tmin = tmax;
			++iterIdx;
		}

		itsNumericals.SetOscilTrace(HWSVNumericals::FirstPeriod,firstOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::NextPeriod,nextOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::LastPeriod,lastOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::NbPeriods,iterIdx);

	//	fclose(f);

		if(iterIdx >= OSCILLATION_NBMAX)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
				" : max iteration reached in option pricing with oscillating integral");
		}


		/// Compute option price
		if(itsNumericals.IsStdFormula())
		{
			pi_1	= 0.5 + Integral_1 / ARM_NumericConstants::ARM_PI;
			pi_2	= 0.5 + Integral_2 / ARM_NumericConstants::ARM_PI;
			expect	= pi_1 - K * pi_2 ;
			
			if (RecPay == K_RCV) expect = expect - 1 + K;

			(*values)[stateIdx] = (zcStart0/K)*expect*payNotional;

			/// Reset turns number for further use
			for(i=0;i<itsFunctionDatas.size();++i)
			{
				itsPrevValues1[i]=0.0;
				itsPrevValues2[i]=0.0;
				itsAngularShifts1[i]=0.0;
				itsAngularShifts2[i]=0.0;
			}
		}
		else
		{
			pi		= Integral_1 / (ARM_NumericConstants::ARM_PI * sqrt(K));
			expect	= 1 - pi ;
			
			if (RecPay == K_PAY) expect = expect + 1/K - 1;

			(*values)[stateIdx] = zcStart0*expect*payNotional;
		}
	}
	return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeAnalyticalOptionPrice
///	Returns: ARM_VectorPtr
///	Action : Compute option price by Gauss-Legendre
///			 numerical integration and analytical
///			 stepwise constant coefficients Riccatis
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeAnalyticalOptionPrice(	double evalTime,
														const ARM_VectorPtr& F0,
														double T0,
														double Te,
														const ARM_GP_Vector& newStrikes,
														int	RecPay,
														double payNotional,
														const ARM_PricingStatesPtr& states) const
{
	size_t modelIdx	= GetModelNb();

	///// Constant Param
	double yfT0 = T0/K_YEAR_LEN;
	double mrs=GetMrs();
	double expT0 = exp (- mrs * yfT0);

	/// Set non SO pricing in numerical stuff
	bool isSOFormula = false;
	itsNumericals.SetIsSOFormula(isSOFormula);

	const ARM_GP_Vector& numericals = itsNumericals.GetFormulaParams();

	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	size_t stateIdx,nbStates=states->size();
	ARM_VectorPtr values(new ARM_GP_Vector(nbStates,0.0));

	ARM_GP_Vector ImU,PsiX1,PsiY1,PsiX0,PsiY0;
	PsiX1.reserve(maxNbPts);
	PsiY1.reserve(maxNbPts);
	PsiX0.reserve(maxNbPts);
	PsiY0.reserve(maxNbPts);

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,T0,states);

	ARM_GP_Vector schedule;
	ComputeRiccatiSchedule(Te,schedule,Te,isSOFormula,evalTime);
	ARM_GP_Vector times(schedule.size()-1);
	size_t i;
	for(i=0;i+1<schedule.size();++i)
		times[i] = schedule[i+1];
	itsNumericals.SetScheduleTrace(times);

	double F00,zcStart0,var0;
	double Integral_1,Integral_2,localIntegral_1,localIntegral_2=0.0;
	double K,lnK,abslnK;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,pi,pi_1,pi_2,expect;
	ARM_IntVector nbSteps;
	size_t nbPts,iterIdx;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		/// V(evalTime)
		var0		= states->GetModelState(stateIdx,modelIdx+V_VARIABLE);
		F00			= (*F0)[stateIdx];
		zcStart0	= (*zcStart)[stateIdx];

		InitAnalyticalData(evalTime,F00,Te,expT0,isSOFormula,schedule);

		Integral_1 = 0.0, Integral_2 = 0.0;

		K		= newStrikes[stateIdx];
		lnK		= log(K);
		abslnK	= fabs(lnK);


		/// 1st oscillation or 1st interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isFirstSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		firstOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// 2nd oscillation or 2nd interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::NextLimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isNextSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		nextOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// Next oscillations or next intervals of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LastLimitStep];
		oscilSpeed			= (abslnK<oscilSpeedLimit ? oscilSpeedLimit : abslnK);
		isLastSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		lastOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;


		itsNumericals.SetOscilTrace(HWSVNumericals::NbPoints,0);
		itsNumericals.ResetNbStepsTrace();

		/// First oscillation integration
		if(itsNumericals.IsStdFormula() && isFirstSpeedLimit)
			nbPts = numericals[HWSVNumericals::FirstLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::FirstNbSteps];
		GaussLegendre_Coefficients GL1(nbPts);
		tmin = 0.0;
		tmax = firstOscilSize;
		scalet=0.5*firstOscilSize;

		t0=tmin+scalet;
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL1.get_point(i) * scalet;

		if(itsNumericals.IsStdFormula())
		{
			PsiX0.resize(nbPts);
			PsiY0.resize(nbPts);
			itsNumericals.SolveAnalyticalHestonSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,PsiX0,PsiY0,var0);
			itsNumericals.IntegrateHestonSystem(GL1,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,localIntegral_1,localIntegral_2);
			Integral_2 += localIntegral_2;
		}
		else
		{
			itsNumericals.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,var0);
			itsNumericals.IntegrateLewisSystem(GL1,ImU,scalet,lnK,PsiX1,PsiY1,localIntegral_1);
		}


		Integral_1 += localIntegral_1;

		itsNumericals.AddOscilTraceNbPoints(nbPts);


		/// Second oscillation integrations
		tmin = tmax;
		tmax += nextOscilSize;
		scalet=0.5*nextOscilSize;
		if(itsNumericals.IsStdFormula() && isNextSpeedLimit)
			nbPts = numericals[HWSVNumericals::NextLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::NextNbSteps];
		GaussLegendre_Coefficients GL2(nbPts);
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		t0=tmin+scalet;

		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL2.get_point(i) * scalet;

		if(itsNumericals.IsStdFormula())
		{
			PsiX0.resize(nbPts);
			PsiY0.resize(nbPts);
			itsNumericals.SolveAnalyticalHestonSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,PsiX0,PsiY0,var0);
			itsNumericals.IntegrateHestonSystem(GL2,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,localIntegral_1,localIntegral_2);
			Integral_2 += localIntegral_2;
		}
		else
		{
			itsNumericals.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,var0);
			itsNumericals.IntegrateLewisSystem(GL2,ImU,scalet,lnK,PsiX1,PsiY1,localIntegral_1);
		}


		Integral_1 += localIntegral_1;

		itsNumericals.AddOscilTraceNbPoints(nbPts);


		/// Next oscillation integrations
		iterIdx=2;
		tmin = tmax;
		scalet=0.5*lastOscilSize;
		if(itsNumericals.IsStdFormula() && isLastSpeedLimit)
			nbPts = numericals[HWSVNumericals::LastLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::LastNbSteps];
		GaussLegendre_Coefficients GL(nbPts);
		ImU.resize(nbPts);
		PsiX1.resize(nbPts);
		PsiY1.resize(nbPts);
		PsiX0.resize(nbPts);
		PsiY0.resize(nbPts);
		while( iterIdx<OSCILLATION_NBMAX &&
			  (fabs(localIntegral_1) > numericals[HWSVNumericals::IntegrationPrecision] ||
			   fabs(localIntegral_2) > numericals[HWSVNumericals::IntegrationPrecision]) )
		{
			tmax += lastOscilSize;
			t0=tmin+scalet;

			for(i=0;i<nbPts;++i)
				ImU[i] = t0 + GL.get_point(i) * scalet;

			if(itsNumericals.IsStdFormula())
			{
				itsNumericals.SolveAnalyticalHestonSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,PsiX0,PsiY0,var0);
				itsNumericals.IntegrateHestonSystem(GL,ImU,scalet,lnK,PsiX1,PsiY1,PsiX0,PsiY0,localIntegral_1,localIntegral_2);
				Integral_2 += localIntegral_2;
			}
			else
			{
				itsNumericals.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX1,PsiY1,var0);
				itsNumericals.IntegrateLewisSystem(GL,ImU,scalet,lnK,PsiX1,PsiY1,localIntegral_1);
			}


			Integral_1 += localIntegral_1;

			itsNumericals.AddOscilTraceNbPoints(nbPts);

			tmin = tmax;
			++iterIdx;
		}

		itsNumericals.SetOscilTrace(HWSVNumericals::FirstPeriod,firstOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::NextPeriod,nextOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::LastPeriod,lastOscilSize);
		itsNumericals.SetOscilTrace(HWSVNumericals::NbPeriods,iterIdx);

		if(iterIdx >= OSCILLATION_NBMAX)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
				" : max iteration reached in option pricing with oscillating integral");
		}


		/// Compute option price
		if(itsNumericals.IsStdFormula())
		{
			pi_1	= 0.5 + Integral_1 / ARM_NumericConstants::ARM_PI;
			pi_2	= 0.5 + Integral_2 / ARM_NumericConstants::ARM_PI;
			expect	= pi_1 - K * pi_2 ;
			
			if(RecPay == K_RCV) expect = expect - 1 + K;

			(*values)[stateIdx] = zcStart0/K * expect * payNotional;
		}
		else
		{
			pi		= Integral_1 / (ARM_NumericConstants::ARM_PI * sqrt(K));
			expect	= 1 - pi ;
			
			if(RecPay == K_PAY) expect = expect + 1/K - 1;

			(*values)[stateIdx] = zcStart0 * expect * payNotional;
		}

	} // for nbStates

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: MaxRate
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::MaxRate(		
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes,
        const ARM_GP_Vector& fwdEndTimes,
        const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& floatPayTimes,
        const ARM_GP_Vector& floatPayPeriods,
        const ARM_GP_Vector& margin,
        double firstReset,
		double firstStart,
		const ARM_VectorPtr& firstRate,
		int MaxOrMin,
		int ResetFreq,
		const ARM_VectorPtr& strikes,
		int CapOrFloor,
		double RhoMinMax,
		bool IsAccrued,
		double MinAccrued,
		double MaxAccrued,
        const ARM_PricingStatesPtr& states) const
{

	if(fabs(floatStartTime - fwdStartTimes[0]) >= 7
	|| fabs(floatEndTime - fwdEndTimes[fwdEndTimes.size()-1]) >= 7
	|| firstRate->size() != states->size())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : input datas inconsistency" );
	}
	size_t i;
	for(i=0;i<margin.size();++i)
	{
		if(fabs(margin[i]) > K_DOUBLE_TOL)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : MaxRate() only available for null margin" );
		}
	}

	/// Convert reset frequency to average lag between to resets
	double maxLag;
	switch(ResetFreq)
	{
		case K_DAILY		: maxLag=365.25/250.0;	break;
		case K_WEEKLY		: maxLag=7.0;			break;
		case K_MONTHLY		: maxLag=365.25/12.0;	break;
		case K_BIMONTHLY	: maxLag=365.25/6.0;	break;
		case K_QUARTERLY	: maxLag=365.25/4.0;	break;
		case K_SEMIANNUAL	: maxLag=365.25/2.0;	break;
		case K_ANNUAL		: maxLag=365.25;		break;
		default :
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : not supported reset frequency in MaxRate()" );
	}

	/// Check if it is necessary to clone !!
	ARM_VectorPtr St = ARM_VectorPtr(new ARM_GP_Vector(*firstRate));

	double tStart=firstReset,tEnd=evalTime;
	double tStartEnd = tEnd-tStart;
	double nbMaxSteps = ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;
	if(maxLag > 0)
		nbMaxSteps = tStartEnd/maxLag;

/****
if(MaxOrMin>0)
{
FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
fprintf(f,"\n\n ---> t=%8.1lf\tT=%8.1lf\n",firstReset,evalTime);
for(i=0;i<nbStates;++i)
	fprintf(f,"St=\t%10.7lf\tST=\t%10.7lf\n",(*St)[i],(*ST)[i]);
fclose(f);
}
****/

	size_t nbStates = states->size();

	/// For each path a conditional SLN approximation is computed
	/// Each forward rate is SLN : shift=1/IT and vol=ZcVolStart-ZcVolEnd
	/// and we get the equivalent SLN on swap rate
	/// Floating leg is assumed to be at fixed leg frequency for efficiency

	double mrs = GetMrs();

	/// Realized variance may be computed using the saved ones
	/// from a previous event time to the current one.
	/// Check that that MaxRate uses it correctly
	const ARM_GP_Vector& modelSchedule = GetModelSchedule();
	size_t nbEvents = modelSchedule.size();
	for(int varIdx=itsVarIdx;varIdx>=0;--varIdx)
	{
		if(modelSchedule[varIdx]<firstReset+K_NEW_DOUBLE_TOL)
			if(modelSchedule[varIdx]>firstReset-K_NEW_DOUBLE_TOL)
				break;
	}
	if(varIdx<=0 || evalTime != modelSchedule[itsVarIdx])
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : Mismatch in start/end for realized variance" );
	size_t startVarIdx = varIdx,endVarIdx = itsVarIdx;

	size_t k,nbFixFlows = fixPayTimes.size();

	/// Compute SLN parameter related to St at previous event time (=firstReset named tStart) and
	/// to ST at current event time (=evalTime named tEnd)
	/// For St, schedules are arbitrary shifted by tEnd-tStart because keyword doesn't give them
	ARM_PricingStatesPtr statesStart = itsRealizedStates[startVarIdx];
	ARM_VectorPtr B0Start	= GetDiscountFunctor()->DiscountFactor(curveName,tStart,floatStartTime-tStartEnd,statesStart);
	ARM_VectorPtr BNStart	= GetDiscountFunctor()->DiscountFactor(curveName,tStart,floatEndTime-tStartEnd,statesStart);
	double nextExpStart,prevExpStart = exp(-mrs*(floatStartTime-tStartEnd)/K_YEAR_LEN);

	ARM_PricingStatesPtr statesEnd = states;
	ARM_VectorPtr B0End = GetDiscountFunctor()->DiscountFactor(curveName,tEnd,floatStartTime,statesEnd);
	ARM_VectorPtr BNEnd = GetDiscountFunctor()->DiscountFactor(curveName,tEnd,floatEndTime,statesEnd);
	double nextExpEnd,prevExpEnd = exp(-mrs*floatStartTime/K_YEAR_LEN);

	ARM_GP_Vector volFactorStart(nbFixFlows),volFactorEnd(nbFixFlows);
	vector<ARM_VectorPtr> BStart(nbFixFlows),BEnd(nbFixFlows);

	for(k=0;k<nbFixFlows;++k)
	{
		BStart[k]			= GetDiscountFunctor()->DiscountFactor(curveName,tStart,fixPayTimes[k]-tStartEnd,statesStart);
		nextExpStart		= exp(-mrs*(fixPayTimes[k]-tStartEnd)/K_YEAR_LEN);
		volFactorStart[k]	= (prevExpStart - nextExpStart)/mrs;
		prevExpStart		= nextExpStart;

		BEnd[k]				= GetDiscountFunctor()->DiscountFactor(curveName,tEnd,fixPayTimes[k],statesEnd);
		nextExpEnd			= exp(-mrs*fixPayTimes[k]/K_YEAR_LEN);
		volFactorEnd[k]		= (prevExpEnd - nextExpEnd)/mrs;
		prevExpEnd			= nextExpEnd;
	}

	ARM_VectorPtr ST = ARM_VectorPtr(new ARM_GP_Vector(nbStates));

	ARM_GP_Vector condShiftStart(nbStates);
	ARM_GP_Vector condVarStart(nbStates);
	ARM_GP_Vector condShiftEnd(nbStates);
	ARM_GP_Vector condVarEnd(nbStates);

	double O1Start,O1End,sigmaBetaStart,sigmaBetaEnd,sigmaStart,sigmaEnd;
	double betaStart,betaEnd,Bkm1Start,Bkm1End,BkStart,BkEnd;
	double localRealizedVar;

/****
FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
fprintf(f,"\n\n ---> t=%8.1lf\tT=%8.1lf\n",firstReset,evalTime);
****/

	for(i=0;i<nbStates;++i)
	{
		O1Start			= 0.0;
		sigmaStart		= 0.0;
		sigmaBetaStart	= 0.0;
		Bkm1Start		= (*B0Start)[i];

		O1End			= 0.0;
		sigmaEnd		= 0.0;
		sigmaBetaEnd	= 0.0;
		Bkm1End			= (*B0End)[i];
		for(k=0;k<nbFixFlows;++k)
		{
			BkStart			= (*(BStart[k]))[i];
			O1Start			+= fixPayPeriods[k]*BkStart;
			sigmaStart		+= Bkm1Start*volFactorStart[k];
			sigmaBetaStart	+= volFactorStart[k];
			Bkm1Start		= BkStart;

			BkEnd			= (*(BEnd[k]))[i];
			O1End			+= fixPayPeriods[k]*BkEnd;
			sigmaEnd		+= Bkm1End*volFactorEnd[k];
			sigmaBetaEnd	+= volFactorEnd[k];
			Bkm1End			= BkEnd;
		}

		localRealizedVar	= itsRealizedVar(endVarIdx,i)-itsRealizedVar(startVarIdx,i);

		sigmaBetaStart		/= nbFixFlows;
		condVarStart[i]		= sigmaBetaStart*sigmaBetaStart*localRealizedVar;
		betaStart			= sigmaBetaStart/sigmaStart*((*B0Start)[i]-(*BNStart)[i]);
		condShiftStart[i]	= (*St)[i]*(1/betaStart-1);

		(*ST)[i]			= ((*B0End)[i]-(*BNEnd)[i])/O1End;
		sigmaBetaEnd		/= nbFixFlows;
		condVarEnd[i]		= sigmaBetaEnd*sigmaBetaEnd*localRealizedVar;
		betaEnd				= sigmaBetaEnd/sigmaEnd*((*B0End)[i]-(*BNEnd)[i]);
		condShiftEnd[i]		= (*ST)[i]*(1/betaEnd-1);

/****
ARM_GP_VectorPtr STRef = SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods,
	fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, true, states);
sigmaStart = sqrt(condVarStart[i]/tStartEnd*K_YEAR_LEN);
sigmaEnd = sqrt(condVarEnd[i]/tStartEnd*K_YEAR_LEN);
fprintf(f,"shiftStart=\t%10.7lf\tbetaStart=\t%10.7lf\tvolStart=\t%10.7lf\tvolLNStart=\t%10.7lf\tshiftEnd=\t%10.7lf\tbetaEnd=\t%10.7lf\tvolEnd=\t%10.7lf\tvolLNEnd=\t%10.7lf\t",
		condShiftStart[i],betaStart,sigmaStart,sigmaStart/betaStart,
		condShiftEnd[i],betaEnd,sigmaEnd,sigmaEnd/betaEnd);
fprintf(f,"ErrSt(bp)=\t%10.7lf\tErrST(bp)=\t%10.7lf\n",
		10000*(((*B0Start)[i]-(*BNStart)[i])/O1Start-(*St)[i]),10000*((*ST)[i]-(*STRef)[i]));
****/

	}

//fclose(f);


	ARM_VectorPtr result(new ARM_GP_Vector(nbStates));

	/// Compute bucket independant uniform draws (but not full "pathorder" compliant
	/// because the transposer mixes paths)
	ARM_RandUniform_NRRan2 randGen(-156);
	ARM_VectorPtr buckets = GetNumMethod()->GetBuckets();
	size_t nbBuckets = buckets->size();
	size_t offset=0,idx=GetNumMethod()->GetBucketIndex();
	if((*buckets)[idx] != nbStates)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : Inconsistency of nbre of paths per bucket" );

	for(i=0;i<idx;++i)
		offset += (*buckets)[i];
	size_t nbTotalPaths = offset;
	for(;i<nbBuckets;++i)
		nbTotalPaths += (*buckets)[i];
	ARM_GP_Vector totalUnif(2*nbTotalPaths);
	randGen.draw(totalUnif);

	ARM_GP_Vector unif(2*nbStates);
	for(i=0;i<2*nbStates;++i)
		unif[i]=totalUnif[offset+i];

	/// Draw MinMax rate then correct it w.r.t. observation frequency
	if(MaxOrMin == 1 || MaxOrMin == -1)
	{
		/// The MaxRate or MinRate is generated
		double u,b,c,d,xt,xT;
		double factorStart,factorEnd,xStart,xEnd,x,xMinMax;
		double blendX = 0.0;
		for(i=0;i<nbStates;++i) 
		{
			xMinMax = MaxOrMin*((*ST)[i] - (*St)[i]) > 0 ? (*ST)[i] : (*St)[i];

			if(IsAccrued)
			{
				if(MaxOrMin>0)
					xMinMax = MaxAccrued < xMinMax ? xMinMax : MaxAccrued;
				else
					xMinMax = MinAccrued > xMinMax ? xMinMax : MinAccrued;
			}

			if(nbMaxSteps < 1.1)
			{
				/// Only 2 resets : St and ST
				(*result)[i] = xMinMax;
			}
			else
			{
				factorStart = 1.0;
				factorEnd = 1.0;
				if(nbMaxSteps < ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER)
				{
					// For Max
					factorStart = 1-0.5826*sqrt(condVarStart[i]/nbMaxSteps);
					factorEnd = 1-0.5826*sqrt(condVarEnd[i]/nbMaxSteps);
				}
				if(MaxOrMin<0)
				{
					// For Min
					factorStart = 2-factorStart;
					factorEnd = 2-factorEnd;
				}

				u = MaxOrMin > 0 ? 1-unif[i] : unif[nbStates+i];

				/// Shift values then draw corrected MinMax and go back to unshifted world
				xt = log((*St)[i] + condShiftStart[i]);
				xT = log((*ST)[i] + condShiftStart[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarStart[i]*log(u);
				d	= b*b -4*c;
				xStart = factorStart*exp(0.5*(-b+MaxOrMin*sqrt(d))) - condShiftStart[i];

				xt = log((*St)[i] + condShiftEnd[i]);
				xT = log((*ST)[i] + condShiftEnd[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarEnd[i]*log(u);
				d	= b*b -4*c;
				xEnd = factorEnd*exp(0.5*(-b+MaxOrMin*sqrt(d))) - condShiftEnd[i];

				x = blendX*xStart + (1-blendX)*xEnd;

				/// MinRate < MIN(St,ST) & MaxRate > MAX(St,ST)
				(*result)[i] = MaxOrMin*(x - xMinMax) < 0 ? xMinMax : x;
			}
		}
	}
	else
	{
		double RhoMinMaxAnti = sqrt(1.0-RhoMinMax*RhoMinMax);

		double u,v,x,y;
		double b,c,d,xt,xT;
		double factorStartMax,factorEndMax,xStartMax,xEndMax,xMax;
		double factorStartMin,factorEndMin,xStartMin,xEndMin,xMin;;
		double xMinMin,xMaxMax;
		double blendX = 0.0;
		for(i=0;i<nbStates;++i) 
		{
			xMaxMax = (*ST)[i] > (*St)[i] ? (*ST)[i] : (*St)[i];
			xMinMin = (*ST)[i] < (*St)[i] ? (*ST)[i] : (*St)[i];

			if(IsAccrued)
			{
				xMaxMax = MaxAccrued < xMaxMax ? xMaxMax : MaxAccrued;
				xMinMin = MinAccrued > xMinMin ? xMinMin : MinAccrued;
			}

			if(nbMaxSteps < 1.1)
			{
				/// Only 2 resets : St and ST
				xMin = xMinMin;
				xMax = xMaxMax;
			}
			else
			{
				/// Generated correlated uniforms
				x	= ARM_GaussianAnalytics::cdfNormal_Inv(unif[i]);
				y	= ARM_GaussianAnalytics::cdfNormal_Inv(unif[nbStates+i]);
				u	= ARM_GaussianAnalytics::cdfNormal(x);
				v	= ARM_GaussianAnalytics::cdfNormal(x*RhoMinMax + y*RhoMinMaxAnti);

				/// Computed correction factor w.r.t. reset frequency
				factorStartMax = 1.0;
				factorEndMax = 1.0;
				factorStartMin = 1.0;
				factorEndMin = 1.0;
				if(nbMaxSteps < ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER)
				{
					// For Max
					factorStartMax = 1-0.5826*sqrt(condVarStart[i]/nbMaxSteps);
					factorEndMax = 1-0.5826*sqrt(condVarEnd[i]/nbMaxSteps);

					// For Min
					factorStartMin = 2-factorStartMax;
					factorEndMin = 2-factorEndMax;
				}


				/// Shift values then draw corrected MAX and go back to unshifted world
				xt = log((*St)[i] + condShiftStart[i]);
				xT = log((*ST)[i] + condShiftStart[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarStart[i]*log(1-u);
				d	= b*b -4*c;
				xStartMax = factorStartMax*exp(0.5*(-b+sqrt(d))) - condShiftStart[i];

				xt = log((*St)[i] + condShiftEnd[i]);
				xT = log((*ST)[i] + condShiftEnd[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarEnd[i]*log(1-u);
				d	= b*b -4*c;
				xEndMax = factorEndMax*exp(0.5*(-b+sqrt(d))) - condShiftEnd[i];

				xMax = blendX*xStartMax + (1-blendX)*xEndMax;

				/// MaxRate > MAX(St,ST)
				xMax = xMax < xMaxMax ? xMaxMax : xMax;



				/// Shift values then draw corrected MIN and go back to unshifted world
				xt = log((*St)[i] + condShiftStart[i]);
				xT = log((*ST)[i] + condShiftStart[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarStart[i]*log(v);
				d	= b*b -4*c;
				xStartMin = factorStartMin*exp(0.5*(-b-sqrt(d))) - condShiftStart[i];

				xt = log((*St)[i] + condShiftEnd[i]);
				xT = log((*ST)[i] + condShiftEnd[i]);
				b = -(xt + xT);
				c	= xt * xT +  0.5*condVarEnd[i]*log(v);
				d	= b*b -4*c;
				xEndMin = factorEndMin*exp(0.5*(-b-sqrt(d))) - condShiftEnd[i];

				xMin = blendX*xStartMin + (1-blendX)*xEndMin;

				/// MinRate < MIN(St,ST)
				xMin = xMin > xMinMin ? xMinMin : xMin;
			}

			x = CapOrFloor*(xMax-xMin-(*strikes)[i]);
			(*result)[i] = x > 0 ? x : 0;
		}
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeSwapRateF0
///	Returns: ARM_VectorPtr
///	Action : compute swap rate dynamics coefficients :
///			 * drift factor due to payment date
///			 * vol factor
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeSwapRateF0(	double evalTime,
												double startTime, 
												double payTime, 
												const ARM_GP_Vector& fixPayTimes,
												const ARM_GP_Vector& fixPayPeriods,
												double mrs,
												const ARM_PricingStatesPtr& states,
												ARM_VectorPtr& swapRate) const
{
	
	size_t fixIdx,nbPeriods = fixPayTimes.size();
	double endTime = fixPayTimes[nbPeriods - 1];

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,startTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,endTime,states);

	double yfs = startTime/K_YEAR_LEN;
	double expMrss = exp(-mrs*yfs);

	double yfe = endTime/K_YEAR_LEN;
	double expMrse = exp(-mrs*yfe);

	double yfp = payTime/K_YEAR_LEN;
	double expMrsp = exp(-mrs*yfp);

	double yff;
	vector<ARM_VectorPtr> zcFixPay(nbPeriods);
	ARM_GP_Vector expMrs(nbPeriods);
	for(fixIdx = 0; fixIdx<nbPeriods; ++fixIdx)
	{
		zcFixPay[fixIdx] = GetDiscountFunctor()->DiscountFactor("",evalTime,fixPayTimes[fixIdx],states);
		yff = fixPayTimes[fixIdx]/K_YEAR_LEN;
		expMrs[fixIdx] = exp (-mrs*yff);
	}

	size_t i,nbStates = zcStart->size();

	ARM_VectorPtr F0(new ARM_GP_Vector(2*nbStates,0.0));
	swapRate->resize(nbStates);

	size_t offset=0;
	double flow,O1,zcs,zce;
	for(i=0;i<nbStates;++i,offset+=2)
	{
		zcs = (*zcStart)[i];
		zce = (*zcEnd)[i];

		O1 = 0.0;
		for(fixIdx = 0; fixIdx<nbPeriods; ++fixIdx)
		{
			flow	= fixPayPeriods[fixIdx] * (*(zcFixPay[fixIdx]))[i];
			O1		+= flow;
			(*F0)[offset]	+= flow * expMrs[fixIdx];
		}
		(*swapRate)[i] = (zcs - zce)/O1;

		(*F0)[offset] /= O1;
		/// Drift factors due to proba change QO1 -> Qpay
		(*F0)[offset+1]	= (expMrsp - (*F0)[offset])/mrs;

		/// Swap rate vol factors
		(*F0)[offset] = ( (zcs*expMrss - zce*expMrse)/(zcs-zce) - (*F0)[offset] )/mrs;
	}

	return F0;
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV1F
///	Routines: VanillaSpreadOption
///	Returns : ARM_VectorPtr
///	Action  : Price an option on a.LongCMS - b.ShortCMS where
///			  LongCMS and ShortCMS are two swap rates paid at time Te
///			  Useful if degenerated in CMS caplet with b=0
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_HWSV1F::VanillaSpreadOptionLet(
				const string& curveName,
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
				const ARM_GP_Vector& strikes,
				double swapLongFloatStartTime,
				double swapLongFloatEndTime,
				const ARM_GP_Vector& swapLongFixPayTimes,
				const ARM_GP_Vector& swapLongFixPayPeriods,
				double swapShortFloatStartTime,
				double swapShortFloatEndTime,
				const ARM_GP_Vector& swapShortFixPayTimes,
				const ARM_GP_Vector& swapShortFixPayPeriods,
				const ARM_PricingStatesPtr& states) const
{
	/// Handle the case of dummy states
	if(states == ARM_PricingStatesPtr(NULL) && evalTime > K_NEW_DOUBLE_TOL)
// FIXMEFRED: mig.vc8 (30/05/2007 16:16:19):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/// If spot evaluation reset the pricing state
	ARM_PricingStatesPtr newStates(states);
	if(evalTime < K_NEW_DOUBLE_TOL)
		newStates = FirstPricingStates(1);
	size_t stateIdx,nbStates = newStates->size();

	/// Compute schedules if necessary (in case of a call through
	/// ARM_VanillaSpreadOptionArg which may be partially built !!)
	double longFloatStartTime = swapLongFloatStartTime;
	double longFloatEndTime = swapLongFloatEndTime;
	double shortFloatStartTime = swapShortFloatStartTime;
	double shortFloatEndTime = swapShortFloatEndTime;
	ARM_GP_Vector *longFixPayTimes=NULL;
	ARM_GP_Vector *longFixPayPeriods=NULL;
	ARM_GP_Vector *shortFixPayTimes=NULL;
	ARM_GP_Vector *shortFixPayPeriods=NULL;

	bool isLong		= coeffLong!=0.0;
	bool isShort	= coeffShort!=0.0;

	/// compute swap schedule since they have not been computed yet !
	ARM_Currency* ccy = GetCurrency(GetModelName());
	double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
	char fixCalendar[100];
	ccy->CalcFixPayCal(fixCalendar);
	int  fixFreq	 = ccy->GetFixedPayFreq();
	int  fixDayCount = ccy->GetFixedDayCount();

	if(isLong)
	{
		if(&swapLongFixPayTimes != NULL)
		{
			longFixPayTimes		= const_cast<ARM_GP_Vector*>(&swapLongFixPayTimes);
			longFixPayPeriods	= const_cast<ARM_GP_Vector*>(&swapLongFixPayPeriods);
		}
		else
		{
			ARM_SwapRatePtr longSwapRate(ARM_SwapRate::CreateSwapRate(asOf,
														asOf+swapLongFloatStartTime, 
														asOf+swapLongFloatEndTime, 
														fixDayCount, 
														fixFreq, 
														fixCalendar));

			longFloatStartTime	= longSwapRate->floatStartTime;
			longFloatEndTime	= longSwapRate->floatEndTime;
			longFixPayTimes		= & longSwapRate->fixPayTimes;
			longFixPayPeriods	= & longSwapRate->fixPayPeriods;
		}
	}
	
	if(isShort)
	{
		if(&swapShortFixPayTimes != NULL)
		{
			shortFixPayTimes	= const_cast<ARM_GP_Vector*>(&swapShortFixPayTimes);
			shortFixPayPeriods	= const_cast<ARM_GP_Vector*>(&swapShortFixPayPeriods);
		}
		else
		{
			ARM_SwapRatePtr shortSwapRate(ARM_SwapRate::CreateSwapRate(asOf,
														asOf+swapShortFloatStartTime, 
														asOf+swapShortFloatEndTime, 
														fixDayCount, 
														fixFreq, 
														fixCalendar));

			shortFloatStartTime	= shortSwapRate->floatStartTime;
			shortFloatEndTime	= shortSwapRate->floatEndTime;
			shortFixPayTimes	= & shortSwapRate->fixPayTimes;
			shortFixPayPeriods	= & shortSwapRate->fixPayPeriods;
		}
	}

	/// Float & fixed leg ends are equal
	if( (isLong && fabs(longFloatEndTime - (*longFixPayTimes)[longFixPayTimes->size()-1])>0.001) ||
	    (isShort && fabs(shortFloatEndTime - (*shortFixPayTimes)[shortFixPayTimes->size()-1])>0.001) )
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end dates are not matching in HWSV1F spread option pricing" );

	ARM_VectorPtr values;
	ARM_VectorPtr refPrices;

#ifdef PRICING_NB_TIMES
ARM_Timer timer;
timer.ClockStartTime();

size_t nbPrices = PRICING_NB_TIMES/10;
for(size_t tIdx=0;tIdx<nbPrices;++tIdx)
{
#endif

	/// Compute vol factors for long and short swap
	double mrs=GetMrs();

	ARM_VectorPtr longSwapRateValue(new ARM_GP_Vector(0)),shortSwapRateValue(new ARM_GP_Vector(0));

	ARM_VectorPtr longF0,shortF0;
	if(isLong)
		longF0 = ComputeSwapRateF0(evalTime,longFloatStartTime,payTime,
					*longFixPayTimes,*longFixPayPeriods,mrs,newStates,longSwapRateValue);

	if(isShort)
		shortF0 = ComputeSwapRateF0(evalTime,shortFloatStartTime,payTime,
					*shortFixPayTimes,*shortFixPayPeriods,mrs,newStates,shortSwapRateValue);


	/// Refresh linked standard H&W model for deterministic volatility reference price
	/// For deterministic correction, V(t=0) is assumed to be equal to Long Term Vol at t=0
	double LongTermVol0 = 1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
	{
		double LongTermVar0 = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
		LongTermVol0 = sqrt(LongTermVar0);
	}
	UpdateStdHWModel(LongTermVol0);
	ARM_HullWhite1F* refModel = dynamic_cast<ARM_HullWhite1F*>(&*(GetAnalyticalModel()));

	double yfTe = resetTime/K_YEAR_LEN;
	double expTe = exp(mrs * yfTe);
	double exp2Te = expTe*expTe;

    double var = static_cast<const ARM_ModelParamsHW1F* const>(GetAnalyticalModel()->GetModelParams())->StateLocalVariance(evalTime,resetTime,resetTime);

	ARM_GP_Vector F0(nbStates);
	ARM_GP_Vector Mu0(nbStates);
	ARM_GP_Vector stdDev(nbStates),newStrikes(nbStates),cvxFwdSpread(nbStates);
	ARM_IntVector statusITM(nbStates);

	size_t offset=0;
	double S1,S2,fwdSpread,volS1,volS2,F00,Mu00;
	double stdFwdFlow,maxDeviation;
	size_t nbDeepITM=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=2)
	{
		/// Spread vol factors
		S1 = isLong ? coeffLong * (*longSwapRateValue)[stateIdx] : 0.0;
		S2 = isShort ? coeffShort * (*shortSwapRateValue)[stateIdx] : 0.0;
		fwdSpread = S1 - S2;
		newStrikes[stateIdx] = strikes[stateIdx] - fwdSpread;

		/// Spread vol factors
		volS1 = isLong ? S1 * (*longF0)[offset] : 0.0;
		volS2 = isShort ? S2 * (*shortF0)[offset] : 0.0;

		F0[stateIdx] = volS1 - volS2;
		F00 = F0[stateIdx] * expTe;

		/// Drift factors under Qpay
		Mu0[stateIdx] = (isLong ? volS1 * (*longF0)[offset+1] : 0.0)
						- (isShort ? volS2 * (*shortF0)[offset+1] : 0.0);


		/// Check if computation is really necessary by computing H&W2F derministic vol datas
		/// then compare spread between convexified fwd and strike vs standard deviation
		stdDev[stateIdx] =	F00 * F00 * var;
		stdDev[stateIdx] = (stdDev[stateIdx]>0.0 ? sqrt(stdDev[stateIdx]) : 0.0);

		/// Compute convexified forward leverage spread
		Mu00 = Mu0[stateIdx] * exp2Te;
		cvxFwdSpread[stateIdx] = fwdSpread + Mu00 * var;

		stdFwdFlow = callPut*(cvxFwdSpread[stateIdx]-strikes[stateIdx]);
		maxDeviation = 16*stdDev[stateIdx];


		if(stdFwdFlow > maxDeviation)
		{
			statusITM[stateIdx] = DEEP_ITM;
			++nbDeepITM;
		}
		else if(stdFwdFlow < - maxDeviation)
			statusITM[stateIdx] = DEEP_OTM;
		else
			statusITM[stateIdx] = STD_ATM;

	} // for nbStates

	/// AsOf deep OTM evaluation : option is worth nothing
	if(nbStates==1 && statusITM[0] == DEEP_OTM)
		return ARM_VectorPtr(new ARM_GP_Vector(1,0.0));


	values = ComputeSpreadOptionPrice(evalTime,F0,Mu0,mrs,payTime,resetTime,newStrikes,callPut,notional,payPeriod,newStates,statusITM);


	/// Compute the reference price for correction using standard Hull & White 1F model
	/// Branching to a 1D numerical integration

	/// Discount factors at eval date will be computed through H&W SV 1F then
	/// change discount/forecast functors
    ARM_ZeroCurveFunctor* refDisFctor = refModel->GetDiscountFunctor();
    ARM_ZeroCurveFunctor* refFixFctor = refModel->GetFixingFunctor();
	refModel->SetDiscountFunctor(GetDiscountFunctor());
	refModel->SetFixingFunctor(GetFixingFunctor());

	const ARM_BoolVector& oldSOFormulaFlags = refModel->GetSOFormulaFlags();
	refModel->SetIsApproxSOFormula(false);
	refModel->SetDeepITMSOFormula(nbDeepITM > 0.90*nbStates); /// speed up VI computation if 90% of states are deep ITM
	refPrices = refModel->VanillaSpreadOptionLet(
			curveName,
			evalTime,
			callPut,
			startTime,
			endTime,
			resetTime,
			payTime,
			payPeriod,
			notional,
			coeffLong,
			coeffShort,
			strikes,
			longFloatStartTime,
			longFloatEndTime,
			*longFixPayTimes,
			*longFixPayPeriods,
			shortFloatStartTime,
			shortFloatEndTime,
			*shortFixPayTimes,
			*shortFixPayPeriods,
			newStates);

	/// Restore original functors & flags
	refModel->SetDiscountFunctor(refDisFctor);
	refModel->SetFixingFunctor(refFixFctor);
	refModel->SetSOFormulaFlags(oldSOFormulaFlags);

	/// Compute standard H&W1F price using same approximation used in stochastic volatility case
	ARM_VectorPtr zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,newStates);
	double approxPrice,ratio = notional * payPeriod;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		approxPrice = ratio * (*zcPay)[stateIdx] * VanillaOption_N(cvxFwdSpread[stateIdx], stdDev[stateIdx], strikes[stateIdx], 1.0, callPut);
		(*values)[stateIdx] = (*values)[stateIdx] + (*refPrices)[stateIdx] - approxPrice;
	}

#ifdef PRICING_NB_TIMES
}
timer.ClockEndTime();
FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
fprintf(f,"Duration = %10.5lf ms\n",timer.GetDuration()*1000.0/nbPrices);
fclose(f);
#endif

	return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV1F
///	Routine: ComputeSpreadOptionPrice
///	Returns: ARM_VectorPtr
///	Action : Compute spread option price by Gauss-Legendre
///			 numerical integration and Runge-Kutta or analytical
///			 Riccati equation solving
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV1F::ComputeSpreadOptionPrice(	double evalTime,
													const ARM_GP_Vector& F0,
													const ARM_GP_Vector& Mu0,
													double mrs,
													double Tp, double Te,
													const ARM_GP_Vector& newStrikes,
													int	callPut,
													double payNotional,
													double payPeriod,
													const ARM_PricingStatesPtr& states,
													const ARM_IntVector& statusITM) const
{

	/// Standard formula (Heston framework) needs to solve a coupled Riccati ODE systems
	/// and Runge-Kutta method is necessary.
	/// In case of (non standard) Lewis framework, if exponential sampling is not set (MaxDecay=0),
	/// Runge Kutta method will be used to solve a simple Riccati ODE system (not implemented yet !).
	/// In case of exponential sampling, this system is solved analitically
	bool isAnalytical = itsNumericalsSO.GetMaxDecay()>0 && !itsNumericalsSO.IsStdFormula();

	/// At the moment Numerical standard (Heston) formula or Analytical enhanced (Lewis) formula
	/// are supported
	if(!(itsNumericalsSO.IsStdFormula() || isAnalytical))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, " : Heston or Lewis + exp sampling only supported" );
	}

	size_t modelIdx	= GetModelNb();

	///// Constant Param
	double yfTp = Tp/K_YEAR_LEN;
	double expTp = exp(-mrs*yfTp);

	/// Set SO pricing in numerical stuff
	/// Only Heston like formula available at the moment
	bool isSOFormula	= true;
	itsNumericalsSO.SetIsSOFormula(isSOFormula);
	ARM_GP_Vector numericals = itsNumericalsSO.GetFormulaParams();

	/// Fourier inversion at u = phi + i.factor. Only used for Lewis formula
	double factor = numericals[HWSVNumericals::InvFourierPoint];
	if(callPut==K_PUT)
		factor = -factor; /// if max(K-S,0), inversion is allowed only for a negative imaginary point

	size_t fctMultiplier = 4 * (itsNumericalsSO.IsStdFormula() ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);


	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	double stateF0,stateMu0;
	size_t stateIdx,nbStates=states->size();
	ARM_VectorPtr values(new ARM_GP_Vector(nbStates,0.0));

	vector<ARM_GP_Vector> solverVars(RK5_NBTMP_VAR);
	size_t nbFcts = maxNbPts * fctMultiplier;
    size_t i;
	for(i=0;i<solverVars.size();++i)
		solverVars[i].reserve(nbFcts);

	ARM_GP_Vector ImU,PsiX,PsiY,PsiTildeX,PsiTildeY;
	PsiX.reserve(maxNbPts);
	PsiY.reserve(maxNbPts);
	PsiTildeX.reserve(maxNbPts);
	PsiTildeY.reserve(maxNbPts);

	ARM_VectorPtr zcPay	= GetDiscountFunctor()->DiscountFactor("",evalTime,Tp,states);

	/// Initialise a RK solver and its internal variable set (for optimisation purpose)
	ARM_RiccatiHWSV1F_SO riccatiSystem(itsSystemDatas,itsSystemDatasSO,itsNumericalsSO.IsStdFormula());

	ARM_GP_Vector schedule;
	if(isAnalytical)
	{
		ComputeRiccatiSchedule(Te,schedule,Te,isSOFormula,evalTime);
		ARM_GP_Vector times(schedule.size()-1);
		for(i=0;i+1<schedule.size();++i)
			times[i] = schedule[i+1];
		itsNumericalsSO.SetScheduleTrace(times);
	}


	double var0,ratio = payPeriod * payNotional;
	double Integral_1,Integral_2,localIntegral_1,localIntegral_2=0.0;
	double K,absK;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,pi_1,pi_2,expect;
	ARM_IntVector nbSteps;
	size_t nbPts,iterIdx;
	double psiTildeZero=0.0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		if(statusITM[stateIdx]==DEEP_OTM)
		{
			/// Deep OTM => do nothing because is worth 0
			continue;
		}

		/// V(evalTime)
		var0 = states->GetModelState(stateIdx,modelIdx+V_VARIABLE);

		/// The option is worth something
		stateF0		= F0[stateIdx];
		stateMu0	= Mu0[stateIdx];

		if(isAnalytical)
			InitAnalyticalData(evalTime,stateF0,Te,expTp,isSOFormula,schedule,stateMu0);
		else
			InitSystemData(evalTime,stateF0,mrs,Te,expTp,maxNbPts,itsNumericalsSO.IsStdFormula(),isSOFormula,stateMu0);

		Integral_1 = 0.0, Integral_2 = 0.0;

		K=newStrikes[stateIdx];
		absK=fabs(K);

		/// 1st oscillation or 1st interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LimitStep];
		oscilSpeed			= (absK<oscilSpeedLimit ? oscilSpeedLimit : absK);
		isFirstSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		firstOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// 2nd oscillation or 2nd interval of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::NextLimitStep];
		oscilSpeed			= (absK<oscilSpeedLimit ? oscilSpeedLimit : absK);
		isNextSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		nextOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

		/// Next oscillations or next intervals of integration if large oscillation period
		oscilSpeedLimit		= ARM_NumericConstants::ARM_2_PI/numericals[HWSVNumericals::LastLimitStep];
		oscilSpeed			= (absK<oscilSpeedLimit ? oscilSpeedLimit : absK);
		isLastSpeedLimit	= (oscilSpeed == oscilSpeedLimit);
		lastOscilSize		= ARM_NumericConstants::ARM_2_PI/oscilSpeed;

//FILE* f=fopen("c:\\temp\\dumpHW1FSV.txt","a");
//fprintf(f,"Theo oscil size=\t%15.10lf\n\n",ARM_NumericConstants::ARM_2_PI/absK);

/****
double u=0.0;
double uStep = numericals[HWSVNumericals::LimitStep]/numericals[HWSVNumericals::FirstNbSteps];
ImU.resize(1);
PsiX.resize(1);
PsiY.resize(1);
double f2 = factor*factor;
for(size_t ii=0;ii<numericals[HWSVNumericals::FirstNbSteps];++ii)
{
	ImU[0]=u;
	itsNumericalsSO.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX,PsiY,var0,factor);
	double costK=cos(u*K);
	double sintK=sin(u*K);
	double den = u*u+f2;
	den *= den;
	double re1 = (u*u-f2)/den;
	double re2 = PsiX[0]*costK-PsiY[0]*sintK;
	double im1 = 2*factor*u/den;
	double im2 = PsiY[0]*costK+PsiX[0]*sintK;
	double x=re1*re2+im1*im2;
	fprintf(f,"u=\t%15.10lf\tRe1=\t%15.10lf\tRe2=\t%15.10lf\tIm1\t%15.10lf\tIm2\t%15.10lf\tf(u)=\t%15.10lf\n",u,re1,re2,im1,im2,x);
	u+=uStep;
}
fprintf(f,"\n\n");
fclose(f);
return values;
****/

		/// Reset computation trace
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NbPoints,0);
		itsNumericalsSO.ResetNbStepsTrace();

		/// Point u=(0,0)
		nbPts=1;
		ImU.resize(nbPts);
		PsiX.resize(nbPts);
		PsiY.resize(nbPts);
		PsiTildeX.resize(nbPts);
		PsiTildeY.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		ImU[0]=0.0;

		if(itsNumericalsSO.IsStdFormula())
		{
			UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);
			itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

			psiTildeZero = PsiTildeX[0];

			if(statusITM[stateIdx]==DEEP_ITM)
			{
				/// Deep ITM => compute the instrinsic value and back to next loop
				(*values)[stateIdx] = (*zcPay)[stateIdx] * ratio * callPut * (psiTildeZero-K);
				continue;
			}
		}


		/// First oscillation integration
		if(itsNumericalsSO.IsStdFormula() && isFirstSpeedLimit)
			nbPts = numericals[HWSVNumericals::FirstLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::FirstNbSteps];
		GaussLegendre_Coefficients GL1(nbPts);
		tmin = 0.0;
		tmax = firstOscilSize;
		scalet=0.5*firstOscilSize;

		t0=tmin+scalet;
		ImU.resize(nbPts);
		PsiX.resize(nbPts);
		PsiY.resize(nbPts);
		PsiTildeX.resize(nbPts);
		PsiTildeY.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL1.get_point(i) * scalet;

		if(isAnalytical)
		{
			itsNumericalsSO.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX,PsiY,var0,factor);
			itsNumericalsSO.IntegrateLewisSystem(GL1,ImU,scalet,K,PsiX,PsiY,localIntegral_1,factor);
		}
		else
		{
			UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

			itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

			itsNumericalsSO.IntegrateSystem(GL1,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
				localIntegral_2,localIntegral_1);

			Integral_2 += localIntegral_2;
		}

		Integral_1 += localIntegral_1;


//fprintf(f,"First Oscil min=\t%15.10lf\tmax=\t%15.10lf\n",tmin,tmax);
//fprintf(f,"localI=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//fprintf(f,"I=\t%15.10lf\t%15.10lf\n\n",Integral_1,Integral_2);


		itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
		if(!isAnalytical)
		{
			nbSteps=riccatiSystem.GetNbSteps();
			nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
			itsNumericalsSO.PushNbStepsTrace(nbSteps);
		}


		/// Second oscillation integrations
		tmin = tmax;
		tmax += nextOscilSize;
		scalet=0.5*nextOscilSize;
		if(itsNumericalsSO.IsStdFormula() && isNextSpeedLimit)
			nbPts = numericals[HWSVNumericals::NextLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::NextNbSteps];
		GaussLegendre_Coefficients GL2(nbPts);
		ImU.resize(nbPts);
		PsiX.resize(nbPts);
		PsiY.resize(nbPts);
		PsiTildeX.resize(nbPts);
		PsiTildeY.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		t0=tmin+scalet;

		for(i=0;i<nbPts;++i)
			ImU[i] = t0 + GL2.get_point(i) * scalet;

		if(isAnalytical)
		{
			itsNumericalsSO.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX,PsiY,var0,factor);
			itsNumericalsSO.IntegrateLewisSystem(GL2,ImU,scalet,K,PsiX,PsiY,localIntegral_1,factor);
		}
		else
		{
			UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

			itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

			itsNumericalsSO.IntegrateSystem(GL2,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
				localIntegral_2,localIntegral_1);

			Integral_2 += localIntegral_2;
		}

		Integral_1 += localIntegral_1;

//fprintf(f,"Second Oscil min=\t%15.10lf\tmax=\t%15.10lf\n",tmin,tmax);
//fprintf(f,"localI=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//fprintf(f,"I=\t%15.10lf\t%15.10lf\n\n",Integral_1,Integral_2);

		itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
		if(!isAnalytical)
		{
			nbSteps=riccatiSystem.GetNbSteps();
			nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
			itsNumericalsSO.PushNbStepsTrace(nbSteps);
		}


		/// Next oscillation integrations
		iterIdx=2;
		tmin = tmax;
		scalet=0.5*lastOscilSize;
		if(itsNumericalsSO.IsStdFormula() && isLastSpeedLimit)
			nbPts = numericals[HWSVNumericals::LastLimitNbSteps];
		else
			nbPts = numericals[HWSVNumericals::LastNbSteps];
		GaussLegendre_Coefficients GL(nbPts);
		ImU.resize(nbPts);
		PsiX.resize(nbPts);
		PsiY.resize(nbPts);
		PsiTildeX.resize(nbPts);
		PsiTildeY.resize(nbPts);
		nbFcts = nbPts * fctMultiplier;
		for(i=0;i<solverVars.size();++i)
			solverVars[i].resize(nbFcts);
		while( iterIdx<OSCILLATION_NBMAX &&
			  (fabs(localIntegral_1) > numericals[HWSVNumericals::IntegrationPrecision] ||
			   fabs(localIntegral_2) > numericals[HWSVNumericals::IntegrationPrecision]) )
		{
			tmax += lastOscilSize;
			t0=tmin+scalet;

			for(i=0;i<nbPts;++i)
				ImU[i] = t0 + GL.get_point(i) * scalet;

			if(isAnalytical)
			{
				itsNumericalsSO.SolveAnalyticalLewisSystem(itsAnalyticalDatas,ImU,PsiX,PsiY,var0,factor);
				itsNumericalsSO.IntegrateLewisSystem(GL,ImU,scalet,K,PsiX,PsiY,localIntegral_1,factor);
			}
			else
			{
				UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

				itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

				itsNumericalsSO.IntegrateSystem(GL,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
					localIntegral_2,localIntegral_1);

				Integral_2 += localIntegral_2;
			}

			Integral_1 += localIntegral_1;

//fprintf(f,"#%3d Oscil min=\t%15.10lf\tmax=\t%15.10lf\n",iterIdx+1,tmin,tmax);
//fprintf(f,"localI=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

			itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
			if(!isAnalytical)
			{
				nbSteps=riccatiSystem.GetNbSteps();
				nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
				itsNumericalsSO.PushNbStepsTrace(nbSteps);
			}

			tmin = tmax;
			++iterIdx;
		}

		itsNumericalsSO.SetOscilTrace(HWSVNumericals::FirstPeriod,firstOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NextPeriod,nextOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::LastPeriod,lastOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NbPeriods,iterIdx);

//fclose(f);

		if(iterIdx >= OSCILLATION_NBMAX)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
				" : max iteration reached in option pricing with oscillating integral");
		}

		/// Compute option price
		if(itsNumericalsSO.IsStdFormula())
		{
			pi_1	= 0.5*psiTildeZero + Integral_1 / ARM_NumericConstants::ARM_PI;
			pi_2	= 0.5 + Integral_2 / ARM_NumericConstants::ARM_PI;
			expect	= pi_1 - K * pi_2 ;
			
			if(callPut == K_PUT) expect += K - psiTildeZero;

			(*values)[stateIdx] = (*zcPay)[stateIdx] * ratio * expect;
		}
		else
		{
			expect = Integral_1 * exp(-factor*K) / ARM_NumericConstants::ARM_PI;

			if(callPut == K_CALL) expect = -expect;

			(*values)[stateIdx] = (*zcPay)[stateIdx] * ratio * expect;
		}

	} // for nbStates

    return values;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

