/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV2F.cpp
 *  \brief Hull & White Stochastic Volatility 2 factors model
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date September 2006
 */

/*--------------------------------------------------------------
  --------------------------------------------------------------

			2 factors Hull & White with stochastic vol

dX1		= [phi11(t)+phi12(t) - mrs1.X1]dt + sig1(t).sqrt(V(t)).dW1
dX2		= [phi12(t)+phi22(t) - mrs2.X2]dt + sig2(t).sqrt(V(t)).dW2
dV		= kappa.[theta(t)-V(t)]dt + nu(t).sqrt(V(t)).dW3

with 

dphiij	= [rhoij(t).sigi(t).sigj(t) - (mrsi+mrsj)].dt

<dWi,dWj> = rhoij(t).dt
--------------------------------------------------------------
--------------------------------------------------------------*/

/// this header comes first as it include some preprocessor constants
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/hwsv2f.h"
#include "gpmodels/hwsv.h"
#include "gpmodels/modelparamshwsv.h"
#include "gpmodels/hw.h"
#include "gpmodels/hw2f.h"



/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/timer.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/irrate.h"

/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

/// gpnumlib
#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/numfunction.h"
#include "gpnumlib/rungekutta.h"


/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )

/// Array aliases
const int X1_VARIABLE		= 0;
const int X2_VARIABLE		= 1;
const int V_VARIABLE		= 2;
const int PHI11_VARIABLE	= 3;
const int PHI12_VARIABLE	= 4;
const int PHI22_VARIABLE	= 5;

const int X1_X2_CORREL		= 0;
const int X1_V_CORREL		= 1;
const int X2_V_CORREL		= 2;

const int X1_VALUE			= 0;
const int X2_VALUE			= 1;
const int NB_X_VALUES		= 2;

const int REAL_X1_X1_VALUE	= 0;
const int REAL_X1_X2_VALUE	= 1;
const int REAL_X2_X2_VALUE	= 2;
const int IMAG_X1_X1_VALUE	= 3;
const int IMAG_X1_X2_VALUE	= 4;
const int IMAG_X2_X2_VALUE	= 5;
const int NB_CPLX_X_VALUES	= 6;

const int DEEP_OTM			= -1;
const int DEEP_ITM			= 1;
const int STD_ATM			= 0;

//#define PRICING_NB_TIMES	1000

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV2F
///	Routine: Constructors & affectation operator
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_RiccatiHWSV2F::ARM_RiccatiHWSV2F(const vector< ARM_SystemData >& systemDatas, bool isStdFormula)
: ARM_RiccatiHWSV(ARM_ODEFunc::RK5Adaptative,isStdFormula), itsSystemDatas(systemDatas)

{
	size_t sysSize = systemDatas.size() > 0 ? systemDatas.size()-1 : 0;
	SetSystemIdx(sysSize);
}

ARM_RiccatiHWSV2F::ARM_RiccatiHWSV2F(const ARM_RiccatiHWSV2F& rhs)
: ARM_RiccatiHWSV(rhs),itsSystemDatas(rhs.itsSystemDatas)
{}

ARM_RiccatiHWSV2F& ARM_RiccatiHWSV2F::operator=(const ARM_RiccatiHWSV2F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_RiccatiHWSV2F();
		new (this) ARM_RiccatiHWSV2F (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV2F
///	Routine: derivs
///	Returns: 
///	Action : Compute derivatives of the Riccati system
////////////////////////////////////////////////////
void ARM_RiccatiHWSV2F::derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const
{
	++itsNbDerivCalls;

	/// For each postive number phi, we compute derivatives of :
	/// for U=(1.0,phi) or U=(phi,0.5) :
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


	/// Search system index (starting with the current one)
	LocateSystemIdx(t);

	size_t fctMultiplier = 4 * (itsIsStdFormula ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);
	size_t i,j,nbPhis = (yt->size()) / fctMultiplier;
	size_t fctIncr = fctMultiplier;

//FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
//fprintf(f,"#26 : t=%8.3lf\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\n",t,
//		(*yt)[fctIncr*26],(*yt)[fctIncr*26+1],(*yt)[fctIncr*26+2],(*yt)[fctIncr*26+3],
//		(*yt)[fctIncr*26+4],(*yt)[fctIncr*26+5],(*yt)[fctIncr*26+6],(*yt)[fctIncr*26+7]);
//fprintf(f,"#27 : t=%8.3lf\tBeta2=(%lf,%lf)\tGamma(%lf,%lf)\tBeta2=(%lf,%lf)\tGamma(%lf,%lf)\n\n",t,
//		(*yt)[fctIncr*27],(*yt)[fctIncr*27+1],(*yt)[fctIncr*27+2],(*yt)[fctIncr*27+3],
//		(*yt)[fctIncr*27+4],(*yt)[fctIncr*27+5],(*yt)[fctIncr*27+6],(*yt)[fctIncr*27+7]);
//fprintf(f,"t,Idx (%3d)=\t%8.3lf\t%3d\n",itsSystemDatas.size()-1,t,itsSystemIdx);
//fclose(f);

	double exp1t		= exp(itsSystemDatas[itsSystemIdx].itsMrs1*t);
	double expSq1t		= exp1t*exp1t;
	double exp2t		= exp(itsSystemDatas[itsSystemIdx].itsMrs2*t);
	double expSq2t		= exp2t*exp2t;
	double exp12t		= exp1t*exp2t;

	double At			= itsSystemDatas[itsSystemIdx].itsModA;

	double kappaTheta	= itsSystemDatas[itsSystemIdx].itsKappaTheta;

	double BtUx0,BtUx1,BtUx;
	double BtUy,Ctx,Cty,alpha2x0,alpha2y0,alpha2x1,alpha2y1,alpha2x,alpha2y;
	if(itsIsStdFormula)
	{
		/// Standard computation with two ODE system to solve at
		/// point U=(1.0,phi) and U=(0.0,phi)
		BtUx0	= itsSystemDatas[itsSystemIdx].itsModB1 +
				  exp1t * itsSystemDatas[itsSystemIdx].itsModB21 +
				  exp2t * itsSystemDatas[itsSystemIdx].itsModB22;

		BtUx1	= BtUx0 +
				  itsSystemDatas[itsSystemIdx].itsModB21x * exp1t +
				  itsSystemDatas[itsSystemIdx].itsModB22x * exp2t;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y(X1_VALUE,i) * exp1t +
					  itsSystemDatas[itsSystemIdx].itsModB2y(X2_VALUE,i) * exp2t;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X2_X2_VALUE,i) * expSq2t;

			Cty		= itsSystemDatas[itsSystemIdx].itsModC(IMAG_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(IMAG_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(IMAG_X2_X2_VALUE,i) * expSq2t;


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
		BtUx	= itsSystemDatas[itsSystemIdx].itsModB1 +
				  ( itsSystemDatas[itsSystemIdx].itsModB21 +
				    itsSystemDatas[itsSystemIdx].itsModB21x ) * exp1t +
				  ( itsSystemDatas[itsSystemIdx].itsModB22 +
				    itsSystemDatas[itsSystemIdx].itsModB22x ) * exp2t;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y(X1_VALUE,i) * exp1t +
					  itsSystemDatas[itsSystemIdx].itsModB2y(X2_VALUE,i) * exp2t;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X2_X2_VALUE,i) * expSq2t;


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
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_RiccatiHWSV2F_SO
///	Routine: derivs
///	Returns: 
///	Action : Compute derivatives of the Riccati system
////////////////////////////////////////////////////
void ARM_RiccatiHWSV2F_SO::derivs(double t, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const
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

//FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
//fprintf(f,"#26 : t=%8.3lf\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\tAlpha2=(%lf,%lf)\tGamma(%lf,%lf)\n",t,
//		(*yt)[fctIncr*26],(*yt)[fctIncr*26+1],(*yt)[fctIncr*26+2],(*yt)[fctIncr*26+3],
//		(*yt)[fctIncr*26+4],(*yt)[fctIncr*26+5],(*yt)[fctIncr*26+6],(*yt)[fctIncr*26+7]);
//fprintf(f,"#27 : t=%8.3lf\tBeta2=(%lf,%lf)\tGamma(%lf,%lf)\tBeta2=(%lf,%lf)\tGamma(%lf,%lf)\n\n",t,
//		(*yt)[fctIncr*27],(*yt)[fctIncr*27+1],(*yt)[fctIncr*27+2],(*yt)[fctIncr*27+3],
//		(*yt)[fctIncr*27+4],(*yt)[fctIncr*27+5],(*yt)[fctIncr*27+6],(*yt)[fctIncr*27+7]);
//fprintf(f,"t,Idx (%3d)=\t%8.3lf\t%3d\n",itsSystemDatas.size()-1,t,itsSystemIdx);
//fclose(f);

	double exp1t		= exp(itsSystemDatas[itsSystemIdx].itsMrs1*t);
	double expSq1t		= exp1t*exp1t;
	double exp2t		= exp(itsSystemDatas[itsSystemIdx].itsMrs2*t);
	double expSq2t		= exp2t*exp2t;
	double exp12t		= exp1t*exp2t;

	double At			= itsSystemDatas[itsSystemIdx].itsModA;

	double Bt			= itsSystemDatasSO[itsSystemIdx].itsModB;

	double kappaTheta	= itsSystemDatas[itsSystemIdx].itsKappaTheta;

	double BtUx,BtTildeUx,BtUy,Ctx,Cty,Dtx,Dty,alpha2x,alpha2y,alphaTilde2x,alphaTilde2y;
	if(itsIsStdFormula)
	{
		/// Standard computation with two ODE system to solve at
		/// point U=(0.0,phi)
		BtUx		= itsSystemDatas[itsSystemIdx].itsModB1 +
					  exp1t * itsSystemDatas[itsSystemIdx].itsModB21 +
				      exp2t * itsSystemDatas[itsSystemIdx].itsModB22;
		BtTildeUx	= exp1t * itsSystemDatas[itsSystemIdx].itsModB21x +
				      exp2t * itsSystemDatas[itsSystemIdx].itsModB22x;

		Dtx		= itsSystemDatasSO[itsSystemIdx].itsModD11 * expSq1t +
				  itsSystemDatasSO[itsSystemIdx].itsModD12 * exp12t +
				  itsSystemDatasSO[itsSystemIdx].itsModD22 * expSq2t;

//FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
//fprintf(f,"idx=%2d\tt=%8.3lf\t[%8.3lf,%8.3lf]\tA=%15.10lf\tB=%15.10lf\tBU=%15.10lf\tBTilU=%15.10lf\tC=%15.10lf\tD=%15.10lf\n",
//		itsSystemIdx,t*365.0,
//		itsSystemDatas[itsSystemIdx==0?0:itsSystemIdx-1].itsYf*365.0,itsSystemDatas[itsSystemIdx].itsYf*365.0,
//		At,Bt,BtUx,BtTildeUx,
//		itsSystemDatas[itsSystemIdx].itsModC11+itsSystemDatas[itsSystemIdx].itsModC12+
//		itsSystemDatas[itsSystemIdx].itsModC22,Dtx);
//fclose(f);




		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y(X1_VALUE,i) * exp1t +
					  itsSystemDatas[itsSystemIdx].itsModB2y(X2_VALUE,i) * exp2t;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X2_X2_VALUE,i) * expSq2t;

			Cty		= itsSystemDatas[itsSystemIdx].itsModC(IMAG_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(IMAG_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(IMAG_X2_X2_VALUE,i) * expSq2t;

			Dty		= itsSystemDatasSO[itsSystemIdx].itsModD(REAL_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatasSO[itsSystemIdx].itsModD(REAL_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatasSO[itsSystemIdx].itsModD(REAL_X2_X2_VALUE,i) * expSq2t;


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
				  ( itsSystemDatas[itsSystemIdx].itsModB21 +
				    itsSystemDatas[itsSystemIdx].itsModB21x ) * exp1t +
				  ( itsSystemDatas[itsSystemIdx].itsModB22 +
				    itsSystemDatas[itsSystemIdx].itsModB22x ) * exp2t;

		for(i=0,j=0;i<nbPhis;++i,j+=fctIncr)
		{
			BtUy	= itsSystemDatas[itsSystemIdx].itsModB2y(X1_VALUE,i) * exp1t +
					  itsSystemDatas[itsSystemIdx].itsModB2y(X2_VALUE,i) * exp2t;

			Ctx		= itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X1_VALUE,i) * expSq1t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X1_X2_VALUE,i) * exp12t +
					  itsSystemDatas[itsSystemIdx].itsModC(REAL_X2_X2_VALUE,i) * expSq2t;


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
///	Class  : ARM_HWSV2F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV2F::ARM_HWSV2F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params,const ARM_GP_Vector& solverParams,
					   int formulaType,const ARM_GP_Vector& formulaParams,
					   int formulaTypeSO,const ARM_GP_Vector& formulaParamsSO,
					   double maxDecay,double maxDecaySO)
:ARM_HWSV(zc,params,ARM_ODEFunc::RK5Adaptative,solverParams,formulaType,formulaParams,formulaTypeSO,formulaParamsSO,maxDecay,maxDecaySO)
{
	ARM_PricingModelPtr hwModel(ARM_PricingModelPtr (static_cast< ARM_PricingModel* >(new ARM_HullWhite2F(  zc ) ) ) );
	SetAnalyticalModel(hwModel);
	UpdateStdHWModel();
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HWSV2F::ARM_HWSV2F(const ARM_HWSV2F& rhs)
: ARM_HWSV(rhs) 
{
	itsSystemDatas	= rhs.itsSystemDatas;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HWSV2F& ARM_HWSV2F::operator=(const ARM_HWSV2F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_HWSV2F();
		new (this) ARM_HWSV2F (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV2F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_HWSV2F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "2F HW Stochastic Volatility Model\n";
    os << indent << "---------------------------------\n\n";

    os << indent << "Numericals for caplet/swaption\n";
    os << indent << "------------------------------\n";
    os << itsNumericals.toString(indent,nextIndent);
    os << indent << "\n\n";

    os << indent << "Numericals for CMS spread option\n";
    os << indent << "--------------------------------\n";
    os << itsNumericalsSO.toString(indent,nextIndent);
    os << indent << "\n\n";

    os << ARM_PricingModel::toString(indent);

    return os.str();
}
////////////////////////////////////////////////////
///	Class   : ARM_HWSV2F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_HWSV2F::Clone() const
{
	return new ARM_HWSV2F(*this);
}

///////////////////////////////////////////////////
///	Class   : ARM_HWSV2F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HWSV2F::FirstPricingStates( size_t bucketSize ) const
{
	double LongTermVar0 = 1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		LongTermVar0 = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(0.0);

	size_t nbModelStates=ModelStatesSize();
	size_t nbNumStates=FactorCount();
	ARM_PricingStatesPtr states(new ARM_PricingStates(bucketSize,nbModelStates,0,nbNumStates));
	int statesNb = states->size();
	for( size_t i=0;i<statesNb; ++i )
	{
		states->SetModelState(i,V_VARIABLE,LongTermVar0); // V(0)=LongTermVar(0)
	}
	return states;
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_HWSV2F::ValidateModelParams(const ARM_ModelParams& params) const
{
    if(!params.DoesModelParamExist(ARM_ModelParamType::MeanReversion) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::Volatility) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolatilityRatio) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::MeanReversionSpread) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::Correlation) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolMeanReversion))
    {
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "At least 1 Model Param is not of a good type!");
    }
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_HWSV2F
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_HWSV2F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    ARM_GP_Vector  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
    case ARM_ModelParamType::Volatility:
    case ARM_ModelParamType::VolatilityRatio:
	case ARM_ModelParamType::Correlation:
	case ARM_ModelParamType::VolOfVol:
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
    case ARM_ModelParamType::MeanReversionSpread:
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
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
            "Unknown type... Model Param Not Supported by HWSV2F" );
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of ND matrix
///	Action : Compute local variances of the state variable 
/// between each time step
////////////////////////////////////////////////////
void ARM_HWSV2F::NumMethodStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelRank= GetModelRank();
    double nextStep,step=timeSteps[0];

	size_t offsetIndex	= (nbSteps-1)*modelRank;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : localDrifts.size() != offsetIndex" );
#endif

	localVariances.resize((nbSteps-1)*(modelRank+1));

	size_t i,j,k,l,nbFactors = FactorCount();

	const ARM_MultiCurve* correlCurves = static_cast<const ARM_CorrelMatParam&>(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation)).GetMultiCurve();
	ARM_GP_Matrix correlMatrix(nbFactors,nbFactors,1.0);
	ARM_GP_Vector correlValues;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];

		/// Interpolate Correlation
		correlValues = correlCurves->Interpolate(step);

		/// [i] => local variance from ti->ti+1
		for (j=0,k=0;k<nbFactors;++k)
		{
			correlMatrix(k,k)=1.0;
			for(l=k+1;l<nbFactors;++l,++j)
			{
				correlMatrix(k,l)=correlValues[j];
				correlMatrix(l,k)=correlMatrix(k,l);
			}
		}
		localVariances[offsetIndex+i] = static_cast<ARM_GP_Matrix*>(correlMatrix.Clone());
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: GetMrs
///	Returns: void
///	Action : fill MRSs with 0 limit value
////////////////////////////////////////////////////
void ARM_HWSV2F::GetMrs(double& mrs1, double& mrs2) const
{
	mrs1 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
    if(fabs(mrs1)<=K_NEW_DOUBLE_TOL)
        mrs1=(mrs1>=0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

	mrs2 = mrs1 + GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversionSpread).GetValueAtPoint(0);
    if(fabs(mrs2)<=K_NEW_DOUBLE_TOL)
        mrs2=(mrs2>=0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

#define LOG_MOMENT_MATCHING

void ARM_HWSV2F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : time index bigger than max size!");
#endif

	double time = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt = (nextTime - time)/K_YEAR_LEN;
	size_t nbFactors= FactorCount();
	size_t nbStates = states->size();
	size_t modelIdx	= GetModelNb();

	if(time<K_NEW_DOUBLE_TOL)
			time+=INTERPOL_EPS;

	/// Mean Reversions are explicitly constant
	double mrs1,mrs2;
	GetMrs(mrs1,mrs2);

	double mrs11 = 2*mrs1;
	double mrs12 = mrs1+mrs2;
	double mrs22 = 2*mrs2;

	/// Volatilities
	double vol1 = GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->Interpolate(time);
	double vol2 = vol1 * GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio).ToCurveModelParam().GetCurve()->Interpolate(time);

	/// Correlations
	ARM_GP_Vector correls(static_cast<const ARM_CorrelMatParam&>(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation)).GetMultiCurve()->Interpolate(time));
	double rho12 = correls[0]; // X1 vs X2

	/// Volatility of variance
	double VoV = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(time);

	double theta=1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		theta = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(time);

	//// Mean reversion of variance
	double kappa = GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()->Interpolate(time);
	double kappadt = kappa*dt;

#ifdef LOG_MOMENT_MATCHING
	double c1 = exp(-kappadt);
	double c2 = 0.5*VoV*VoV/kappa*(1-c1*c1);
#endif


//FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
//fprintf(f,"t=\t%8.3lf\n",time);

	/// Euler diffusion scheme with local lognormal correction
	/// on variance to ensure positive values
	double x1State,x2State,vState,phi11State,phi12State,phi22State,vol1State,vol2State;
	double x1,x2,var,stdDev,dx1,dx2,dv,phi11,phi12,phi22,meanVi,varVi;
	double sdt = sqrt(dt);
	for( size_t i=0;i<nbStates; ++i )
	{
		x1		= states->GetModelState(i,modelIdx+X1_VARIABLE);
		x2		= states->GetModelState(i,modelIdx+X2_VARIABLE);
		var		= states->GetModelState(i,modelIdx+V_VARIABLE);
		phi11	= states->GetModelState(i,modelIdx+PHI11_VARIABLE);
		phi12	= states->GetModelState(i,modelIdx+PHI12_VARIABLE);
		phi22	= states->GetModelState(i,modelIdx+PHI22_VARIABLE);

		stdDev	= sqrt(var);
		vol1State = vol1 * stdDev;
		vol2State = vol2 * stdDev;

		/// Standard correlated gaussian vector value
		dx1 = states->GetNumMethodState(i,modelIdx+X1_VARIABLE);
		dx2 = states->GetNumMethodState(i,modelIdx+X2_VARIABLE);
		dv	= states->GetNumMethodState(i,modelIdx+V_VARIABLE);

		/// Phi(i,j)
		phi11State = phi11 + (vol1State*vol1State - mrs11*phi11)*dt;
		phi12State = phi12 + (vol1State*vol2State*rho12 - mrs12*phi12)*dt;
		phi22State = phi22 + (vol2State*vol2State - mrs22*phi22)*dt;

		states->SetModelState(i,modelIdx+PHI11_VARIABLE,phi11State);
		states->SetModelState(i,modelIdx+PHI12_VARIABLE,phi12State);
		states->SetModelState(i,modelIdx+PHI22_VARIABLE,phi22State);

		phi11State = 0.5*(phi11State+phi11);
		phi12State = 0.5*(phi12State+phi12);
		phi22State = 0.5*(phi22State+phi22);

		/// Rate factors
		x1State = x1 + (phi11State+phi12State - mrs1*x1)*dt + sdt*vol1State*dx1;
		x2State = x2 + (phi12State+phi22State - mrs2*x2)*dt + sdt*vol2State*dx2;

		states->SetModelState(i,modelIdx+X1_VARIABLE,x1State);
		states->SetModelState(i,modelIdx+X2_VARIABLE,x2State);


		/// Variance : lognormal + moment matching or ln(Var) diffusion
#ifdef LOG_MOMENT_MATCHING
		//// Andreasen
		meanVi = theta + (var-theta)*c1;
		if(meanVi < 0.0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : variance process is negative in HWSV2F model!");
		varVi = 1 + c2*var/(meanVi*meanVi);
		varVi = log(varVi);

		vState = meanVi*exp(-0.5*varVi+ sqrt(varVi)*dv);
#else
		if(var<K_NEW_DOUBLE_TOL)
			var = K_NEW_DOUBLE_TOL;
		vState = var * exp((kappa*(theta-var) - 0.5*VoV*VoV)/var*dt + sdt*VoV/stdDev*dv);
#endif

		states->SetModelState(i,modelIdx+V_VARIABLE,vState);

/***
		fprintf(f,"#%5d\tx1(t)=\t15.10lf\tx2(t)=\t15.10lf\tV(t)=\t15.10lf\tphi11(t)=\t15.10lf\tphi12(t)=\t15.10lf\tphi22(t)=\t15.10lf\n",
					i,
					states->GetModelState(i,modelIdx+X1_VARIABLE),
					states->GetModelState(i,modelIdx+X2_VARIABLE),
					states->GetModelState(i,modelIdx+V_VARIABLE),
					states->GetModelState(i,modelIdx+PHI11_VARIABLE),
					states->GetModelState(i,modelIdx+PHI12_VARIABLE),
					states->GetModelState(i,modelIdx+PHI22_VARIABLE));

fclose(f);
***/
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ComputeRiccatiSchedule
///	Returns: 
///	Action : Compute a schedule based on model schedule
///			 that starts at Tstart and ends at T
////////////////////////////////////////////////////
void ARM_HWSV2F::ComputeRiccatiSchedule(double T, ARM_GP_Vector& schedule, double Tref, bool isSOFormula, double Tstart) const
{
	schedule.empty();
	ARM_GP_Vector modelSched = static_cast<const ARM_ModelParamsHWSV* const>(GetModelParams())->GetSchedule();
	if(modelSched[0] != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : HWSV2F model schedule must start at spot time (t=0)" );

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


	double maxDecay=(isSOFormula ? CC_Min(itsNumericalsSO.GetMaxDecay(),itsNumericals.GetMaxDecay())
								 : itsNumericals.GetMaxDecay());
	if(maxDecay>0)
	{
		/// Add points to keep |exp(-(MRSi+MRSj).(Tref-t))-exp(-(MRSi+MRSj).(Tref-ti))| < MaxDecay for all t in ]ti,ti+1]
		ARM_GP_Vector additionalPoints;
		double mrs1,mrs2;
		GetMrs(mrs1,mrs2);
		double factor11=-2*mrs1/K_YEAR_LEN;
		double factor12=-(mrs1+mrs2)/K_YEAR_LEN;
		double factor22=-2*mrs2/K_YEAR_LEN;
		double factor = fabs(factor11) > fabs(factor12) ? factor11 : factor12;
		factor = fabs(factor22) > fabs(factor) ? factor22 : factor;
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
///	Class  : ARM_HWSV2F
///	Routine: InitAnalyticalData
///	Returns: 
///	Action : Initialise model dependent constant set
///			 for Riccati analytical solving
///			 Model parameters are assumed to be stepwise
///			 right constant. Csts saved at [i] are then
///			 related to ]ti, ti+1] interval
///			 exp1T0 = exp(-lambda1*T0/365.0) where
///			 exp2T0 = exp(-lambda2*T0/365.0) where
///			 T0 = swap or libor rate start time (lag from asOfDate)
///			 Te = swap or libor rate reset time (lag from asOfDate)
////////////////////////////////////////////////////
void ARM_HWSV2F::InitAnalyticalData(double evalTime,const ARM_VectorPtr& F0, double Te, double exp1T0, double exp2T0,
									bool isSOFormula, ARM_GP_Vector& schedule) const
{
	double Tref = Te;
	if(schedule.size()==0)
		ComputeRiccatiSchedule(Te,schedule,Tref,isSOFormula,evalTime);

    size_t nbTimes = schedule.size();
    itsAnalyticalDatas.resize(nbTimes-1);

	double lambda1,lambda2;
	GetMrs(lambda1,lambda2);
	double lambda1Yf = lambda1/K_YEAR_LEN;
	double lambda2Yf = lambda2/K_YEAR_LEN;

	double F10	= (*F0)[X1_VARIABLE];
	double F20	= (*F0)[X2_VARIABLE];

	double exp1Tref		= exp(lambda1Yf*Tref);
	double exp1Tref2	= exp1Tref*exp1Tref;
	double F10exp1Tref	= F10*exp1Tref;
	exp1T0				*= exp1Tref;

	double exp2Tref		= exp(lambda2Yf*Tref);
	double exp2Tref2	= exp2Tref*exp2Tref;
	double F20exp2Tref	= F20*exp2Tref;
	exp2T0				*= exp2Tref;

	double exp12T0		= exp1T0*exp2T0;

	double cst11C		= 0.5*F10exp1Tref*F10exp1Tref;
	double cst12C		= F10exp1Tref*F20exp2Tref;
	double cst22C		= 0.5*F20exp2Tref*F20exp2Tref;

	const ARM_Curve* nuCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve();
	const ARM_Curve* volCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve();
	const ARM_Curve* kappaCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve();
	ARM_Curve* thetaCurve=NULL;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		thetaCurve = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve();

	const ARM_CorrelMatParam& correls = static_cast<const ARM_CorrelMatParam&>(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation));
	ARM_GP_Vector rhos;

	double t=schedule[0],nextt,nu,nu2,vol1,vol2,rhoNuVol1,rhoNuVol2,kappa,theta;
	double exp1t=exp(-lambda1Yf*(Tref-t)),exp1Nextt,exp1Proxy,exp1ProxyVol;
	double exp2t=exp(-lambda2Yf*(Tref-t)),exp2Nextt,exp2Proxy,exp2ProxyVol;
    for(size_t i=0;i+1<nbTimes;++i)
    {
		nextt	= schedule[i+1];
		nu		= nuCurve->Interpolate(nextt);
		nu2		= nu*nu;
		vol1	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve()->Interpolate(nextt);
		vol2	= vol1*GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio).ToCurveModelParam().GetCurve()->Interpolate(nextt);
		rhos	= correls.GetMultiCurve()->Interpolate(nextt);

		rhoNuVol1= nu*vol1*rhos[X1_V_CORREL];
		rhoNuVol2= nu*vol2*rhos[X2_V_CORREL];

		kappa	= kappaCurve->Interpolate(nextt);
		theta	= (thetaCurve!=NULL ? thetaCurve->Interpolate(nextt) : 1.0);

		exp1Nextt	= exp(-lambda1Yf*(Tref-nextt));
		exp1Proxy	= 0.5*(exp1t+exp1Nextt);

		exp2Nextt	= exp(-lambda2Yf*(Tref-nextt));
		exp2Proxy	= 0.5*(exp2t+exp2Nextt);


		itsAnalyticalDatas[i].itsdt	= (nextt-t)/K_YEAR_LEN;
		itsAnalyticalDatas[i].itsNu2= nu2;
		itsAnalyticalDatas[i].itsA	= -0.5*nu2;
		itsAnalyticalDatas[i].itsB1	= kappa + rhoNuVol1*(1-exp1T0*exp1Proxy)/lambda1
											+ rhoNuVol2*(1-exp2T0*exp2Proxy)/lambda2;
		itsAnalyticalDatas[i].itsB2	= rhoNuVol1*F10exp1Tref*exp1Proxy + rhoNuVol2*F20exp2Tref*exp2Proxy;

		exp1ProxyVol	= exp1Proxy*vol1;
		exp2ProxyVol	= exp2Proxy*vol2;
		itsAnalyticalDatas[i].itsC2	= cst11C*exp1ProxyVol*exp1ProxyVol +
									  cst12C*exp1ProxyVol*exp2ProxyVol*rhos[X1_X2_CORREL] +
									  cst22C*exp2ProxyVol*exp2ProxyVol;
		/// To be updated !!
		if(isSOFormula)
			itsAnalyticalDatas[i].itsC1	= 0.0;
		else
			itsAnalyticalDatas[i].itsC1	= 0.0;

		itsAnalyticalDatas[i].itsD	= -kappa*theta;

		t=nextt;
		exp1t=exp1Nextt;
		exp2t=exp2Nextt;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: InitSystemData
///	Returns: 
///	Action : Initialise model dependent constant set
///			 for Riccati system solving
///			 Model parameters are assumed to be stepwise
///			 right constant. Csts saved at [i] are then
///			 related to ]ti, ti+1] interval
///			 exp1Tn = exp(-mrs1*Tn/365.0) & exp2Tn = exp(-mrs2*Tn/365.0) where
///			 Tn = numeraire time (lag from asOfDate)
///			 Te = swap or libor rate reset time (lag from asOfDate)
////////////////////////////////////////////////////
void ARM_HWSV2F::InitSystemData(double evalTime,const ARM_VectorPtr& F0, double mrs1, double mrs2, double Te, double exp1Tn, double exp2Tn, size_t maxSize,
								bool isStdFormula, bool isSOFormula, const ARM_GP_MatrixPtr& Mu0) const
{
	ARM_GP_Vector schedule;
	ComputeRiccatiSchedule(Te,schedule,Te,isSOFormula,evalTime);

    size_t nbTimes = schedule.size();
    itsSystemDatas.resize(nbTimes);
    //itsSystemDatas.resize(nbTimes-1);

	double cst1B = -exp1Tn/mrs1;
	double cst2B = -exp2Tn/mrs2;
	double cst11C = 0.5 * (*F0)[X1_VARIABLE] * (*F0)[X1_VARIABLE];
	double cst12C = (*F0)[X1_VARIABLE] * (*F0)[X2_VARIABLE];
	double cst22C = 0.5 * (*F0)[X2_VARIABLE] * (*F0)[X2_VARIABLE];

	double cst11D,cst12D,cst22D;
	if(isSOFormula)
	{
		cst11D = - (*Mu0)(X1_VARIABLE,X1_VARIABLE);
		cst12D = - ((*Mu0)(X1_VARIABLE,X2_VARIABLE)+(*Mu0)(X2_VARIABLE,X1_VARIABLE));
		cst22D = - (*Mu0)(X2_VARIABLE,X2_VARIABLE);
		itsSystemDatasSO.resize(nbTimes);
		//itsSystemDatasSO.resize(nbTimes-1);
	}

	const ARM_Curve* nuCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve();
	const ARM_Curve* volCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam().GetCurve();
	const ARM_Curve* volRatioCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio).ToCurveModelParam().GetCurve();
	const ARM_Curve* kappaCurve	= GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve();
	ARM_Curve* thetaCurve=NULL;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
		thetaCurve = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve();

	const ARM_CorrelMatParam& correls = static_cast<const ARM_CorrelMatParam&>(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation));
	ARM_GP_Vector rhos;
	double t,teps,nu,vol1,vol2,rhoNuVol1,rhoNuVol2,kappa,theta,v11,v12,v22,nu2;
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
		vol1	= volCurve->Interpolate(teps);
		vol2	= vol1 * volRatioCurve->Interpolate(teps);
		rhos	= correls.GetMultiCurve()->Interpolate(teps);

		rhoNuVol1= nu*vol1*rhos[X1_V_CORREL];
		rhoNuVol2= nu*vol2*rhos[X2_V_CORREL];

		kappa	= kappaCurve->Interpolate(teps);
		theta	= (thetaCurve!=NULL ? thetaCurve->Interpolate(teps) : 1.0);


		itsSystemDatas[i].itsTime		= t;
		itsSystemDatas[i].itsYf			= t/K_YEAR_LEN;
		itsSystemDatas[i].itsMrs1		= mrs1;
		itsSystemDatas[i].itsMrs2		= mrs2;
		itsSystemDatas[i].itsKappaTheta	= kappa*theta;
		nu2 = nu*nu;
		itsSystemDatas[i].itsModA		= -0.5*nu2;
		itsSystemDatas[i].itsModB1		= kappa + rhoNuVol1/mrs1 + rhoNuVol2/mrs2;
		itsSystemDatas[i].itsModB21		= rhoNuVol1*cst1B;
		itsSystemDatas[i].itsModB22		= rhoNuVol2*cst2B;
		itsSystemDatas[i].itsModB21x	= (isStdFormula ? -rhoNuVol1 * (*F0)[X1_VARIABLE] : -0.5*rhoNuVol1 * (*F0)[X1_VARIABLE]);
		itsSystemDatas[i].itsModB22x	= (isStdFormula ? -rhoNuVol2 * (*F0)[X2_VARIABLE] : -0.5*rhoNuVol2 * (*F0)[X2_VARIABLE]);
		v11 = vol1*vol1;
		v12 = vol1*vol2*rhos[X1_X2_CORREL];
		v22 = vol2*vol2;
		itsSystemDatas[i].itsModC11	= cst11C*v11;
		itsSystemDatas[i].itsModC12	= cst12C*v12;
		itsSystemDatas[i].itsModC22	= cst22C*v22;

		/// Reserve memory space
		itsSystemDatas[i].itsModB2y.reserve(NB_X_VALUES,maxSize);
		itsSystemDatas[i].itsModC.reserve(NB_CPLX_X_VALUES,maxSize);

		if(isSOFormula)
		{
			itsSystemDatasSO[i].itsModB		= -nu2;
			itsSystemDatasSO[i].itsModD11	= cst11D*v11;
			itsSystemDatasSO[i].itsModD12	= cst12D*v12;
			itsSystemDatasSO[i].itsModD22	= cst22D*v22;
			itsSystemDatasSO[i].itsModD.reserve(NB_CPLX_X_VALUES,maxSize);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: UpdateUDependentSystemData
///	Returns: 
///	Action : Update Riccati constants which depend
///			 on the complex number. Imaginary parts are
///			 input because real parts are already known (0 or 1)
////////////////////////////////////////////////////
void ARM_HWSV2F::UpdateUDependentSystemData(const ARM_GP_Vector& ImU, bool isStdFormula, bool isSOFormula) const
{
    size_t i,nbTimes = itsSystemDatas.size();
	size_t j,nbU = ImU.size();
	double b21x,b22x,b21y,b22y,c11,c12,c22,y,y2,c11y,c12y,c22y,d11,d12,d22;
	if(isStdFormula)
	{
		for(i=0;i<nbTimes;++i)
		{
			itsSystemDatas[i].itsModB2y.resize(NB_X_VALUES,nbU);
			itsSystemDatas[i].itsModC.resize(NB_CPLX_X_VALUES,nbU);
			if(isSOFormula)
			{
				itsSystemDatasSO[i].itsModD.resize(NB_CPLX_X_VALUES,nbU);
				d11 = itsSystemDatasSO[i].itsModD11;
				d12 = itsSystemDatasSO[i].itsModD12;
				d22 = itsSystemDatasSO[i].itsModD22;
			}

			b21x = itsSystemDatas[i].itsModB21x;
			b22x = itsSystemDatas[i].itsModB22x;
			c11 = itsSystemDatas[i].itsModC11;
			c12 = itsSystemDatas[i].itsModC12;
			c22 = itsSystemDatas[i].itsModC22;
			for(j=0;j<nbU;++j)
			{
				y = ImU[j];
				itsSystemDatas[i].itsModB2y(X1_VALUE,j) = b21x*y;
				itsSystemDatas[i].itsModB2y(X2_VALUE,j) = b22x*y;
				c11y = c11*y;
				c12y = c12*y;
				c22y = c22*y;
				itsSystemDatas[i].itsModC(REAL_X1_X1_VALUE,j)	= c11y*y;
				itsSystemDatas[i].itsModC(REAL_X1_X2_VALUE,j)	= c12y*y;
				itsSystemDatas[i].itsModC(REAL_X2_X2_VALUE,j)	= c22y*y;
				itsSystemDatas[i].itsModC(IMAG_X1_X1_VALUE,j)	= c11y;
				itsSystemDatas[i].itsModC(IMAG_X1_X2_VALUE,j)	= c12y;
				itsSystemDatas[i].itsModC(IMAG_X2_X2_VALUE,j)	= c22y;

				if(isSOFormula)
				{
					itsSystemDatasSO[i].itsModD(REAL_X1_X1_VALUE,j) = d11*y;
					itsSystemDatasSO[i].itsModD(REAL_X1_X2_VALUE,j) = d12*y;
					itsSystemDatasSO[i].itsModD(REAL_X2_X2_VALUE,j) = d22*y;
				}
			}
		}
	}
	else
	{
		for(i=0;i<nbTimes;++i)
		{
			itsSystemDatas[i].itsModB2y.resize(NB_X_VALUES,nbU);
			itsSystemDatas[i].itsModC.resize(NB_CPLX_X_VALUES,nbU);
			b21y = 2 * itsSystemDatas[i].itsModB21x;
			b22y = 2 * itsSystemDatas[i].itsModB22x;
			c11 = itsSystemDatas[i].itsModC11;
			c12 = itsSystemDatas[i].itsModC12;
			c22 = itsSystemDatas[i].itsModC22;
			if(isSOFormula)
			{
				itsSystemDatasSO[i].itsModD.resize(NB_CPLX_X_VALUES,nbU);
				d11 = itsSystemDatasSO[i].itsModD11;
				d12 = itsSystemDatasSO[i].itsModD12;
				d22 = itsSystemDatasSO[i].itsModD22;
			}
			for(j=0;j<nbU;++j)
			{
				y = ImU[j];
				y2 = y*y+0.25;
				itsSystemDatas[i].itsModB2y(X1_VALUE,j) = b21y*y;
				itsSystemDatas[i].itsModB2y(X2_VALUE,j) = b22y*y;
				itsSystemDatas[i].itsModC(REAL_X1_X1_VALUE,j)	= c11*y2;
				itsSystemDatas[i].itsModC(REAL_X1_X2_VALUE,j)	= c12*y2;
				itsSystemDatas[i].itsModC(REAL_X2_X2_VALUE,j)	= c22*y2;

				if(isSOFormula)
				{
					itsSystemDatasSO[i].itsModD(REAL_X1_X1_VALUE,j) = d11*y2;
					itsSystemDatasSO[i].itsModD(REAL_X1_X2_VALUE,j) = d12*y2;
					itsSystemDatasSO[i].itsModD(REAL_X2_X2_VALUE,j) = d22*y2;
				}
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: UpdateStdHWModel
///	Returns: void
///	Action : Update the deterministic volatility
///			 reference H&W model
////////////////////////////////////////////////////
void ARM_HWSV2F::UpdateStdHWModel(double volFactor) const
{
	ARM_ModelParamVector paramVector(5);
	ARM_CurveModelParam  paramVol(GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility).ToCurveModelParam());
	for(size_t i=0;i<paramVol.size();++i)
		paramVol.SetValueAtPoint(i,paramVol.GetValueAtPoint(i)*volFactor);
	ARM_CurveModelParam  paramMRS(GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam());
	ARM_CurveModelParam  paramVolRatio(GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio).ToCurveModelParam());
	ARM_CurveModelParam  paramMRSSpread(GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversionSpread).ToCurveModelParam());
	paramVector[0] = &paramVol;
	paramVector[1] = &paramMRS;
	paramVector[2] = &paramVolRatio;
	paramVector[3] = &paramMRSSpread;

	/// Restore time depedent correlations and extract X1 vs X2 correlation curve
	ARM_CurveModelParam	 paramCorrel(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation).ToCurveModelParam());
	const ARM_MultiCurve* correls = static_cast<const ARM_CorrelMatParam&>(GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation)).GetMultiCurve();
	ARM_GP_Vector correlTimes(correls->GetAbscisses());
	size_t nbTimes = correlTimes.size();
	ARM_GP_Vector correlValues(nbTimes);
	for(i=0;i<nbTimes;++i)
		correlValues[i]=(correls->GetOrdinate(i))[X1_X2_CORREL];
	paramCorrel.SetValuesAndTimes(&correlTimes,&correlValues);
	paramVector[4] = &paramCorrel;

    if(paramVolRatio.GetCurve()->GetAbscisses().size() <= 1 && nbTimes <= 1 )
        GetAnalyticalModel()->SetModelParams(ARM_ModelParamsHW2FStd(paramVector));
    else
        GetAnalyticalModel()->SetModelParams(ARM_ModelParamsHW2FExt(paramVector));
}


////////////////////////////////////////////////////
///	Class   : ARM_HWSV2F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::LocalDiscounts(
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

	size_t nbStates		= states->size();
	ARM_GP_Vector* result	= new ARM_GP_Vector(nbStates,0.0);
	size_t modelIdx			= GetModelNb();

    /// Simple mapping version  : r(t) = f(0,t) + X1(t) + X2(t)
    double X1,X2;
    dt /= K_YEAR_LEN;
	for( size_t i=0; i<nbStates; ++i )
	{
		X1 = states->GetModelState(i,modelIdx+X1_VARIABLE);
		X2 = states->GetModelState(i,modelIdx+X2_VARIABLE);
		(*result)[i] = discount * exp(-(X1+X2)*dt);
    }

	return ARM_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
	// Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
    double zct;
	ARM_GP_Vector* betas=NULL;

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,zcT) );
    }

    // State variable factors
    int i,nbStates=states->size();
    if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
    {
        if(evalTime < maturityTime)
        {
            // Could be optimized using phis(t)...?
            betas = ARM_ModelParamsHW2F::BetatT(GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion),
				GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversionSpread),evalTime,maturityTime);
            zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
        }
        else
            return ARM_VectorPtr(new ARM_GP_Vector(nbStates,1.0));
    }
	else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : only cash numeraire allowed");


	size_t modelNb = GetModelNb();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
	double temp = 0.;
	double zctT = zcT/zct;
	double beta1 = (*betas)[X1_VARIABLE];
	double beta2 = (*betas)[X2_VARIABLE];
	double beta11 = 0.5*beta1*beta1;
	double beta12 = beta1*beta2;
	double beta22 = 0.5*beta2*beta2;
    for(i=0;i<nbStates;++i)
    {
		(*values)[i] = zctT *
			exp( -beta1	* states->GetModelState(i,modelNb+X1_VARIABLE)
				 -beta2	* states->GetModelState(i,modelNb+X2_VARIABLE)
                 -beta11 * states->GetModelState(i,modelNb+PHI11_VARIABLE)
                 -beta12 * states->GetModelState(i,modelNb+PHI12_VARIABLE)
                 -beta22 * states->GetModelState(i,modelNb+PHI22_VARIABLE) );
	}
    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ComputeCapletF0
///	Returns: ARM_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeCapletF0(double evalTime,
								double fwdStartTime, 
								double fwdEndTime,
								double fwdPeriod,
								const ARM_GP_Vector& strikes,
								double mrs1,double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_GP_Vector& newStrikes) const
{	
	size_t stateIdx,nbStates = states->size();
	ARM_VectorPtr F0(new ARM_GP_Vector(NB_X_VALUES*nbStates,0.0));

	double yfs = fwdStartTime/K_YEAR_LEN;
	double yfe = fwdEndTime/K_YEAR_LEN;
	double F01 = ( exp(-mrs1*yfe) - exp(-mrs1*yfs) )/mrs1;
	double F02 = ( exp(-mrs2*yfe) - exp(-mrs2*yfs) )/mrs2;

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdStartTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,fwdEndTime,states);

	size_t offset=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=NB_X_VALUES)
	{
		(*F0)[offset+X1_VARIABLE] = F01;
		(*F0)[offset+X2_VARIABLE] = F02;
		newStrikes[stateIdx] = (1.0 + fwdPeriod * strikes[stateIdx]) * (*zcEnd)[stateIdx] / (*zcStart)[stateIdx];
		newStrikes[stateIdx] = 1.0 / newStrikes[stateIdx];
	}

	return F0;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::VanillaCaplet(
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
// FIXMEFRED: mig.vc8 (30/05/2007 16:16:00):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/// If spot evaluation reset the pricing state
	ARM_PricingStatesPtr newStates(states);
	if(evalTime < K_NEW_DOUBLE_TOL)
		newStates = FirstPricingStates(1);
	size_t nbStates = newStates->size();

	if(strikesPerState.size() != nbStates)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+ " : strike is missing in HWSV2F caplet pricing" );

	/// No payment lag adjustement computation
	if(fabs(payTime-fwdEndTime) > 5)
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+ " : convexity adjustment not available for HWSV2F model");
	}

	ARM_VectorPtr values;

#ifdef PRICING_NB_TIMES
ARM_Timer timer;
timer.ClockStartTime();

for(size_t tIdx=0;tIdx<PRICING_NB_TIMES;++tIdx)
{
#endif

	double mrs1,mrs2;
	GetMrs(mrs1,mrs2);

	ARM_GP_Vector newStrikes(nbStates);
	ARM_VectorPtr F0(ComputeCapletF0(evalTime,fwdStartTime,fwdEndTime,fwdPeriod,strikesPerState,mrs1,mrs2,newStates,newStrikes));

	if(nbStates==1 && newStrikes[0] <= 0.0)
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
		if(itsNumericals.GetMaxDecay()>0)
			/// Exponential terms are sampled then Riccati systems have analytical solutions
			values = ComputeAnalyticalOptionPrice(evalTime,F0,mrs1,mrs2,fwdStartTime,fwdResetTime,newStrikes,capFloor,payNotional*period/fwdPeriod,newStates);
		else
			/// Riccatti system are solved using the numerical Runge-Kutta method
			values = ComputeRungeKuttaOptionPrice(evalTime,F0,mrs1,mrs2,fwdStartTime,fwdResetTime,newStrikes,capFloor,payNotional*period/fwdPeriod,newStates);
	}


#ifdef PRICING_NB_TIMES
}
timer.ClockEndTime();
FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
fprintf(f,"Duration = %10.5lf ms\n",timer.GetDuration()*1000.0/PRICING_NB_TIMES);
fclose(f);
#endif

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ComputeSwaptionF0
///	Returns: ARM_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeSwaptionF0( double evalTime,
								double startTime, 
								double endTime,
								const ARM_GP_Vector& fixPayTimes,
								const ARM_GP_Vector& fixPayPeriods,
								int callPut,
								const ARM_GP_Matrix& strikes,
								bool isConstantNotional,
								const ARM_GP_Vector& fixNotional,
								const ARM_GP_Vector& floatNotional,
								int refNotionalIdx,
								double mrs1, double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_GP_Vector& newStrikes) const
{
/****/
	size_t stateIdx,nbStates = states->size();

	double yfs = startTime/K_YEAR_LEN;
	double expMrs1s = exp(-mrs1*yfs);
	double expMrs2s = exp(-mrs2*yfs);

	double yfe = endTime/K_YEAR_LEN;
	double expMrs1e = exp(-mrs1*yfe);
	double expMrs2e = exp(-mrs2*yfe);

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,startTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,endTime,states);

	size_t i,nbPeriods = fixPayTimes.size();

	if(endTime != fixPayTimes[nbPeriods - 1])
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : end time mismatch in bond strip vol computation");

	double yff;
	vector<ARM_VectorPtr> zcFixPay(nbPeriods);
	ARM_GP_Vector expMrs1(nbPeriods),expMrs2(nbPeriods);
	for(i = 0; i<nbPeriods; ++i)
	{
		zcFixPay[i] = GetDiscountFunctor()->DiscountFactor("",evalTime,fixPayTimes[i],states);
		yff = fixPayTimes[i]/K_YEAR_LEN;
		expMrs1[i] = exp (-mrs1*yff);
		expMrs2[i] = exp (-mrs2*yff);
	}

	ARM_VectorPtr F0(new ARM_GP_Vector(NB_X_VALUES*nbStates,0.0));
	size_t offset=0;

	double flow,newStrike,F01,F02;
	if(isConstantNotional)
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=NB_X_VALUES)
		{
			newStrike=0.0,F01=0.0,F02=0.0;
			for(i=0;i+1<nbPeriods;++i)
			{
				flow		= strikes(stateIdx,i) * fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F01	+= flow * (expMrs1[i] - expMrs1s);
				F02	+= flow * (expMrs2[i] - expMrs2s);
			}
			flow		= (1.0 + strikes(stateIdx,nbPeriods-1) * fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			newStrike	= (newStrike + flow) / (*zcStart)[stateIdx];
			F01	= (F01 + flow * (expMrs1e - expMrs1s)) / (*zcStart)[stateIdx];
			F02	= (F02 + flow * (expMrs2e - expMrs2s)) / (*zcStart)[stateIdx];

			newStrike	= 1.0 / newStrike;
			F01	*= (newStrike/mrs1);
			F02	*= (newStrike/mrs2);

			(*F0)[offset+X1_VARIABLE]	= F01;
			(*F0)[offset+X2_VARIABLE]	= F02;
			newStrikes[stateIdx] = newStrike;
		}
	}
	else
	{
		/// Case where first fixed flows have non null notional but floating leg ones are null
		for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=NB_X_VALUES)
		{
			newStrike=0.0,F01=0.0,F02=0.0;
			for(i=0;i<refNotionalIdx;++i)
			{
				if(fabs(fixNotional[i]) > K_NEW_DOUBLE_TOL)
				{
					flow		= strikes(stateIdx,i) * fixNotional[i] *
								  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
					newStrike	+= flow;
					F01	+= flow * (expMrs1[i] - expMrs1s);
					F02	+= flow * (expMrs2[i] - expMrs2s);
				}
			}
			for(i=refNotionalIdx;i+1<nbPeriods;++i)
			{
				flow		= floatNotional[i] - floatNotional[i+1] + strikes(stateIdx,i) * fixNotional[i] *
							  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F01	+= flow * (expMrs1[i] - expMrs1s);
				F02	+= flow * (expMrs2[i] - expMrs2s);
			}
			flow		= (floatNotional[nbPeriods-1] + strikes(stateIdx,nbPeriods-1) * fixNotional[nbPeriods-1] *
						  fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			double coefStart = (*zcStart)[stateIdx] * floatNotional[refNotionalIdx];
			newStrike	= (newStrike + flow)/coefStart;
			F01	= (F01 + flow * (expMrs1e - expMrs1s))/coefStart;
			F02	= (F02 + flow * (expMrs2e - expMrs2s))/coefStart;

			newStrike	= 1.0 / newStrike;
			F01	*= (newStrike/mrs1);
			F02	*= (newStrike/mrs2);

			(*F0)[offset+X1_VARIABLE]	= F01;
			(*F0)[offset+X2_VARIABLE]	= F02;
			newStrikes[stateIdx] = newStrike;
		}
	}

	return F0;
/****/
/****
	size_t stateIdx,nbStates = states->size();

	double yf = startTime/K_YEAR_LEN;
	double expMrs1s = exp(-mrs1*yf);
	double expMrs2s = exp(-mrs2*yf);

	yf = endTime/K_YEAR_LEN;
	double expMrs1e = exp(-mrs1*yf);
	double expMrs2e = exp(-mrs2*yf);

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,startTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,endTime,states);

	size_t i,nbPeriods = fixPayTimes.size();
	vector<ARM_VectorPtr> zcFixPay(nbPeriods);
	ARM_GP_Vector expMrs1(nbPeriods),expMrs2(nbPeriods);
	for(i = 0; i<nbPeriods; ++i)
	{
		zcFixPay[i] = GetDiscountFunctor()->DiscountFactor("",evalTime,fixPayTimes[i],states);
		expMrs1[i] = exp(-mrs1*(fixPayTimes[i]-evalTime)/K_YEAR_LEN);
		expMrs2[i] = exp(-mrs2*(fixPayTimes[i]-evalTime)/K_YEAR_LEN);
	}

	if(endTime != fixPayTimes[nbPeriods - 1])
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : end time mismatch in bond strip vol computation");

	ARM_VectorPtr F0(new ARM_GP_Vector(NB_X_VALUES*nbStates,0.0));

	double flow,newStrike,F1,F2;
	size_t offset=0;
	if(isConstantNotional)
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=NB_X_VALUES)
		{
			newStrike=0.0,F1=0.0,F2=0.0;
			for(i = 0; i<nbPeriods-1; i++)
			{
				flow		= strikes(stateIdx,i) * fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F1	+= flow * (expMrs1[i] - expMrs1s);
				F2	+= flow * (expMrs2[i] - expMrs2s);
			}
			flow		= (1.0 + strikes(stateIdx,nbPeriods-1) * fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			newStrike	= (newStrike + flow)/(*zcStart)[stateIdx];
			F1	= (F1 + flow * (expMrs1e - expMrs1s))/(*zcStart)[stateIdx];
			F2	= (F2 + flow * (expMrs2e - expMrs2s))/(*zcStart)[stateIdx];

			newStrike	= 1.0 / newStrike;
			F1	*= (newStrike/mrs1);
			F2	*= (newStrike/mrs2);

			(*F0)[offset+X1_VARIABLE]	= F1;
			(*F0)[offset+X2_VARIABLE]	= F2;
			newStrikes[stateIdx] = newStrike;
		}
	}
	else
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset+=NB_X_VALUES)
		{
			/// Case where first fixed flows have non null notional but floating leg ones are null
			newStrike=0.0,F1=0.0,F2=0.0;
			for(i = 0; i<refNotionalIdx; i++)
			{
				if(fabs(fixNotional[i]) > K_NEW_DOUBLE_TOL)
				{
					flow		= strikes(stateIdx,i) * fixNotional[i] *
								  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
					newStrike	+= flow;
					F1	+= flow * (expMrs1[i] - expMrs1s);
					F2	+= flow * (expMrs2[i] - expMrs2s);
				}
			}
			for(i = refNotionalIdx; i<nbPeriods-1; i++)
			{
				flow		= floatNotional[i] - floatNotional[i+1] + strikes(stateIdx,i) * fixNotional[i] *
							  fixPayPeriods[i] * (*(zcFixPay[i]))[stateIdx];
				newStrike	+= flow;
				F1	+= flow * (expMrs1[i] - expMrs1s);
				F2	+= flow * (expMrs2[i] - expMrs2s);
			}
			flow		= (floatNotional[nbPeriods-1] + strikes(stateIdx,nbPeriods-1) * fixNotional[nbPeriods-1] *
						  fixPayPeriods[nbPeriods - 1]) * (*zcEnd)[stateIdx];
			double coefStart = (*zcStart)[stateIdx] * floatNotional[refNotionalIdx];
			newStrike	= (newStrike + flow)/coefStart;
			F1	= (F1 + flow * (expMrs1e - expMrs1s))/coefStart;
			F2	= (F2 + flow * (expMrs2e - expMrs2s))/coefStart;

			newStrike	= 1.0 / newStrike;
			F1	*= (newStrike/mrs1);
			F2	*= (newStrike/mrs2);

			(*F0)[offset+X1_VARIABLE]	= F1;
			(*F0)[offset+X2_VARIABLE]	= F2;
			newStrikes[stateIdx] = newStrike;
		}
	}


	return F0;
****/
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
///			 If the swaption is standard, this method
///          calls ARM_HullWhite::VanillaSwaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::VanillaSwaption(
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
// FIXMEFRED: mig.vc8 (30/05/2007 16:16:12):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/// If spot evaluation reset the pricing state
	ARM_PricingStatesPtr newStates(states);
	if(evalTime < K_NEW_DOUBLE_TOL)
		newStates = FirstPricingStates(1);
	size_t i,nbStates = newStates->size();

	/// First row is always used
	if(strikesPerState.cols() < 1 || strikesPerState.rows() != nbStates)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : strike profile is missing in HWSV2F swaption pricing" );

	/// Float & fixed leg ends are equal
	if(fabs(floatEndTime - fixPayTimes[fixPayTimes.size()-1])>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end dates are not matching in HWSV2F swaption pricing" );

	int refNotionalIdx=0;
	if(isConstantNotional)
	{
		/// Same notional on float & fixed legs is constant
		if(floatNotional[0] != fixNotional[0] || fabs(floatNotional[0])<K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, " : same non null notional is required on float & fixed legs if bullet in HWSV2F swaption pricing" );
	}
	else
	{
		/// Same notional profile size on float & fixed legs
		if(floatNotional.size() != fixNotional.size())
			ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fix legs are required to be of same frequency in HWSV2F swaption pricing" );

		/// Locate first non null floating leg notional
		for(refNotionalIdx=0;refNotionalIdx<floatNotional.size();++refNotionalIdx)
		{
			if(fabs(floatNotional[refNotionalIdx]) >= K_NEW_DOUBLE_TOL)
				break;
		}
		if(refNotionalIdx >= floatNotional.size()) 
			ARM_THROW( ERR_INVALID_ARGUMENT, " : can't locate a non null notional on floating leg in HWSV2F swaption pricing" );
	}
	double refNotional = floatNotional[refNotionalIdx];

	ARM_VectorPtr values;

#ifdef PRICING_NB_TIMES
ARM_Timer timer;
timer.ClockStartTime();

for(size_t tIdx=0;tIdx<PRICING_NB_TIMES;++tIdx)
{
#endif

	double mrs1,mrs2;
	GetMrs(mrs1,mrs2);

	/// Start time = start of the 1st period where floating leg notional is not null
	double startTime = (refNotionalIdx==0 ? floatStartTime : fixPayTimes[refNotionalIdx-1]);
	ARM_GP_Vector newStrikes(nbStates);
	ARM_VectorPtr F0(ComputeSwaptionF0(evalTime,startTime,floatEndTime,fixPayTimes,fixPayPeriods,callPut,
		strikesPerState,isConstantNotional,fixNotional,floatNotional,refNotionalIdx,mrs1,mrs2,newStates,newStrikes));

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
		if(itsNumericals.GetMaxDecay()>0)
			/// Exponential terms are sampled then Riccati systems have analytical solutions
			values = ComputeAnalyticalOptionPrice(evalTime,F0,mrs1,mrs2,startTime,swapResetTime,newStrikes,callPut,refNotional,newStates);
		else
			/// Riccatti system are solved using the numerical Runge-Kutta method
			values = ComputeRungeKuttaOptionPrice(evalTime,F0,mrs1,mrs2,startTime,swapResetTime,newStrikes,callPut,refNotional,newStates);

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


		/// Compute standard H&W1F price using same approximation used in stochastic volatility case
		double yfTe = swapResetTime/K_YEAR_LEN;
		double expMrs1e = exp(mrs1 * yfTe);
		double expMrs2e = exp(mrs2 * yfTe);

		ARM_GP_TriangularMatrix* covars = static_cast<const ARM_ModelParamsHW2FExt* const>(GetAnalyticalModel()->GetModelParams())->StateLocalVariance(evalTime,swapResetTime);
		double v11 = (*covars)(X1_VARIABLE,X1_VARIABLE);
		double v12 = (*covars)(X1_VARIABLE,X2_VARIABLE);
		double v22 = (*covars)(X2_VARIABLE,X2_VARIABLE);

		/// Free memory
		delete covars;

		ARM_VectorPtr zcStart = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,startTime,newStates);

		/// Payer swaption => Put on bond, Receiver swaption => Call on bond
		int optType = callPut==K_CALL ? K_PUT : K_CALL;
		double F10,F20,stdDev,approxPrice;
		size_t offset=0;
		for(i=0;i<nbStates;++i,offset+=NB_X_VALUES)
		{
			F10			= (*F0)[offset+X1_VARIABLE] * expMrs1e;
			F20			= (*F0)[offset+X2_VARIABLE] * expMrs2e;
			stdDev		= F10*(F10 * v11 + 2*F20 * v12) + F20*F20 * v22;
			stdDev		= (stdDev>0.0 ? sqrt(stdDev) : 0.0);
			approxPrice	= refNotional/newStrikes[i] * BlackSholes_Formula(1.0,stdDev,(*zcStart)[i],newStrikes[i],optType);

			/// Correct closed form formula with H&W2F std prices
			(*values)[i] += (*refPrices)[i] - approxPrice;
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
/***** Formula correction *****/
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
///	Class  : ARM_HWSV2F
///	Routine: ComputeRungeKuttaOptionPrice
///	Returns: ARM_VectorPtr
///	Action : Compute option price by Gauss-Legendre
///			 numerical integration and Runge-Kutta
///			 Riccati equation solving
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeRungeKuttaOptionPrice( double evalTime,
														const ARM_VectorPtr& F0,
														double mrs1, double mrs2,
														double T0, double Te,
														const ARM_GP_Vector& newStrikes,
														int	RecPay,
														double payNotional,
														const ARM_PricingStatesPtr& states) const
{
	size_t modelIdx	= GetModelNb();

	///// Constant Param
	double yfT0 = T0/K_YEAR_LEN;
	double exp1T0 = exp (- mrs1 * yfT0);
	double exp2T0 = exp (- mrs2 * yfT0);

	/// Set non SO pricing in numerical stuff
	bool isSOFormula = false;
	itsNumericals.SetIsSOFormula(isSOFormula);

	const ARM_GP_Vector& numericals = itsNumericals.GetFormulaParams();

	size_t fctMultiplier = 4 * (itsNumericals.IsStdFormula() ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);

	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	ARM_VectorPtr stateF0(new ARM_GP_Vector(NB_X_VALUES));

	size_t stateIdx,nbStates=states->size();
	ARM_VectorPtr values(new ARM_GP_Vector(nbStates,0.0));

	vector<ARM_GP_Vector> solverVars(RK5_NBTMP_VAR);
	size_t nbFcts = maxNbPts * fctMultiplier;
	size_t i;
	for(i=0;i<solverVars.size();++i)
		solverVars[i].reserve(nbFcts);

	ARM_GP_Vector ImU,PsiX1,PsiY1,PsiX0,PsiY0;
	PsiX1.reserve(maxNbPts);
	PsiY1.reserve(maxNbPts);
	PsiX0.reserve(maxNbPts);
	PsiY0.reserve(maxNbPts);

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,T0,states);

	/// Initialise a RK solver and its internal variable set (for optimisation purpose)
	ARM_RiccatiHWSV2F riccatiSystem(itsSystemDatas,itsNumericals.IsStdFormula());

	double var0;
	double Integral_1,Integral_2,localIntegral_1,localIntegral_2=0.0;
	double K,lnK,abslnK;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,pi,pi_1,pi_2,expect;
	ARM_IntVector nbSteps;
	size_t nbPts,iterIdx,offset=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset += NB_X_VALUES)
	{
		/// V(evalTime)
		var0 = states->GetModelState(stateIdx,modelIdx+V_VARIABLE);

		(*stateF0)[X1_VARIABLE] = (*F0)[offset+X1_VARIABLE];
		(*stateF0)[X2_VARIABLE] = (*F0)[offset+X2_VARIABLE];

		InitSystemData(evalTime,stateF0,mrs1,mrs2,Te,exp1T0,exp2T0,maxNbPts,itsNumericals.IsStdFormula(),isSOFormula);

		Integral_1 = 0.0, Integral_2 = 0.0;

		K = newStrikes[stateIdx];
		lnK = log(K);
		abslnK=fabs(lnK);

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


		/// Reset computation trace
		itsNumericals.SetOscilTrace(HWSVNumericals::NbPoints,0);
		itsNumericals.ResetNbStepsTrace();


//FILE *f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
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
			pi_1   = 0.5 + Integral_1 / ARM_NumericConstants::ARM_PI;
			pi_2	  = 0.5 + Integral_2 / ARM_NumericConstants::ARM_PI;
			expect = pi_1 - K * pi_2 ;
			
			if(RecPay == K_RCV) expect = expect - 1 + K;

			(*values)[stateIdx] = (*zcStart)[stateIdx]/K * expect * payNotional;
		}
		else
		{
			pi   = Integral_1 / (ARM_NumericConstants::ARM_PI * sqrt(K));
			expect = 1 - pi ;
			
			if(RecPay == K_PAY) expect = expect + 1/K - 1;

			(*values)[stateIdx] = (*zcStart)[stateIdx] * expect * payNotional;
		}

	} // for nbStates

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ComputeAnalyticalOptionPrice
///	Returns: ARM_VectorPtr
///	Action : Compute option price by Gauss-Legendre
///			 numerical integration and analytical
///			 stepwise constant coefficients Riccatis
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeAnalyticalOptionPrice(	double evalTime,
														const ARM_VectorPtr& F0,
														double mrs1, double mrs2,
														double T0, double Te,
														const ARM_GP_Vector& newStrikes,
														int	RecPay,
														double payNotional,
														const ARM_PricingStatesPtr& states) const
{
	size_t modelIdx	= GetModelNb();

	///// Constant Param
	double yfT0 = T0/K_YEAR_LEN;
	double exp1T0 = exp (- mrs1 * yfT0);
	double exp2T0 = exp (- mrs2 * yfT0);

	/// Set non SO pricing in numerical stuff
	bool isSOFormula = false;
	itsNumericals.SetIsSOFormula(isSOFormula);

	const ARM_GP_Vector& numericals = itsNumericals.GetFormulaParams();

	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	ARM_VectorPtr stateF0(new ARM_GP_Vector(NB_X_VALUES));

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

	double var0;
	double Integral_1,Integral_2,localIntegral_1,localIntegral_2=0.0;
	double K,lnK,abslnK;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,pi,pi_1,pi_2,expect;
	ARM_IntVector nbSteps;
	size_t nbPts,iterIdx,offset=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offset += NB_X_VALUES)
	{
		/// V(evalTime)
		var0 = states->GetModelState(stateIdx,modelIdx+V_VARIABLE);

		(*stateF0)[X1_VARIABLE] = (*F0)[offset+X1_VARIABLE];
		(*stateF0)[X2_VARIABLE] = (*F0)[offset+X2_VARIABLE];

		InitAnalyticalData(evalTime,stateF0,Te,exp1T0,exp2T0,isSOFormula,schedule);

		Integral_1 = 0.0, Integral_2 = 0.0;

		K = newStrikes[stateIdx];
		lnK = log(K);
		abslnK=fabs(lnK);

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


		/// Reset computation trace
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

			(*values)[stateIdx] = (*zcStart)[stateIdx]/K * expect * payNotional;
		}
		else
		{
			pi		= Integral_1 / (ARM_NumericConstants::ARM_PI * sqrt(K));
			expect	= 1 - pi ;
			
			if(RecPay == K_PAY) expect = expect + 1/K - 1;

			(*values)[stateIdx] = (*zcStart)[stateIdx] * expect * payNotional;
		}

	} // nbStates

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HWSV2F
///	Routine: ComputeSwapRateF0
///	Returns: ARM_VectorPtr
///	Action : compute swap rate dynamics coefficients :
///			 * drift factor due to payment date
///			 * vol factor
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeSwapRateF0(double evalTime,
								double startTime, 
								double payTime, 
								const ARM_GP_Vector& fixPayTimes,
								const ARM_GP_Vector& fixPayPeriods,
								double mrs1, double mrs2,
								const ARM_PricingStatesPtr& states,
								ARM_VectorPtr& swapRate) const
{
	
	size_t fixIdx,nbPeriods = fixPayTimes.size();
	double endTime = fixPayTimes[nbPeriods - 1];

	ARM_VectorPtr zcStart	= GetDiscountFunctor()->DiscountFactor("",evalTime,startTime,states);
	ARM_VectorPtr zcEnd		= GetDiscountFunctor()->DiscountFactor("",evalTime,endTime,states);

	double yfs = startTime/K_YEAR_LEN;
	double expMrs1s = exp(-mrs1*yfs);
	double expMrs2s = exp(-mrs2*yfs);

	double yfe = endTime/K_YEAR_LEN;
	double expMrs1e = exp(-mrs1*yfe);
	double expMrs2e = exp(-mrs2*yfe);

	double yfp = payTime/K_YEAR_LEN;
	double expMrs1p = exp(-mrs1*yfp);
	double expMrs2p = exp(-mrs2*yfp);

	double yff;
	vector<ARM_VectorPtr> zcFixPay(nbPeriods);
	ARM_GP_Vector expMrs1(nbPeriods),expMrs2(nbPeriods);
	for(fixIdx = 0; fixIdx<nbPeriods; ++fixIdx)
	{
		zcFixPay[fixIdx] = GetDiscountFunctor()->DiscountFactor("",evalTime,fixPayTimes[fixIdx],states);
		yff = fixPayTimes[fixIdx]/K_YEAR_LEN;
		expMrs1[fixIdx] = exp (-mrs1*yff);
		expMrs2[fixIdx] = exp (-mrs2*yff);
	}

	size_t i,nbStates = zcStart->size();

	ARM_VectorPtr F0(new ARM_GP_Vector(2*NB_X_VALUES*nbStates,0.0));
	swapRate->resize(nbStates);

	size_t offset=0,offsetX=2*NB_X_VALUES;
	double flow,O1,zcs,zce;
	for(i=0;i<nbStates;++i,offset+=offsetX)
	{
		zcs = (*zcStart)[i];
		zce = (*zcEnd)[i];

		O1 = 0.0;
		for(fixIdx = 0; fixIdx<nbPeriods; ++fixIdx)
		{
			flow	= fixPayPeriods[fixIdx] * (*(zcFixPay[fixIdx]))[i];
			O1		+= flow;
			(*F0)[offset+X1_VARIABLE]	+= flow * expMrs1[fixIdx];
			(*F0)[offset+X2_VARIABLE]	+= flow * expMrs2[fixIdx];
		}
		(*swapRate)[i] = (zcs - zce)/O1;

		(*F0)[offset+X1_VARIABLE] /= O1;
		(*F0)[offset+X2_VARIABLE] /= O1;

		/// Drift factors due to proba change QO1 -> Qpay
		(*F0)[offset+NB_X_VALUES+X1_VARIABLE]	= (expMrs1p - (*F0)[offset+X1_VARIABLE])/mrs1;
		(*F0)[offset+NB_X_VALUES+X2_VARIABLE]	= (expMrs2p - (*F0)[offset+X2_VARIABLE])/mrs2;

		/// Swap rate vol factors
		(*F0)[offset+X1_VARIABLE] = ( (zcs*expMrs1s - zce*expMrs1e)/(zcs-zce) - (*F0)[offset+X1_VARIABLE] )/mrs1;
		(*F0)[offset+X2_VARIABLE] = ( (zcs*expMrs2s - zce*expMrs2e)/(zcs-zce) - (*F0)[offset+X2_VARIABLE] )/mrs2;
	}

	return F0;
}



////////////////////////////////////////////////////
///	Class   : ARM_HWSV
///	Routines: VanillaSpreadOption
///	Returns : ARM_VectorPtr
///	Action  : Price an option on a.LongCMS - b.ShortCMS where
///			  LongCMS and ShortCMS are two swap rates paid at time Te
///			  Useful if degenerated in CMS caplet with b=0
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_HWSV2F::VanillaSpreadOptionLet(const string& curveName,
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
		ARM_THROW( ERR_INVALID_ARGUMENT, " : float & fixed leg end dates are not matching in HWSV2F spread option pricing" );

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
	double mrs1,mrs2;
	GetMrs(mrs1,mrs2);

	ARM_VectorPtr longSwapRateValue(new ARM_GP_Vector(0)),shortSwapRateValue(new ARM_GP_Vector(0));

	ARM_VectorPtr longF0,shortF0;
	if(isLong)
		longF0 = ComputeSwapRateF0(evalTime,longFloatStartTime,payTime,
					*longFixPayTimes,*longFixPayPeriods,mrs1,mrs2,newStates,longSwapRateValue);

	if(isShort)
		shortF0 = ComputeSwapRateF0(evalTime,shortFloatStartTime,payTime,
			*shortFixPayTimes,*shortFixPayPeriods,mrs1,mrs2,newStates,shortSwapRateValue);

	/// Refresh linked standard H&W model for deterministic volatility reference price
	/// For deterministic correction, V(t=0) is assumed to be equal to Long Term Vol at t=0
	double LongTermVol0 = 1.0;
	if(GetModelParams()->DoesModelParamExist( ARM_ModelParamType::LongTermVol))
	{
		double LongTermVar0 = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
		LongTermVol0 = sqrt(LongTermVar0);
	}
	UpdateStdHWModel(LongTermVol0);
	ARM_HullWhite2F* refModel = dynamic_cast<ARM_HullWhite2F*>(&*(GetAnalyticalModel()));

	double yfTe = resetTime/K_YEAR_LEN;
	double exp1Te = exp(mrs1 * yfTe);
	double exp2Te = exp(mrs2 * yfTe);
	double exp11Te = exp1Te*exp1Te;
	double exp12Te = exp1Te*exp2Te;
	double exp22Te = exp2Te*exp2Te;

    ARM_GP_TriangularMatrix* covars = static_cast<const ARM_ModelParamsHW2FExt* const>(GetAnalyticalModel()->GetModelParams())->StateLocalVariance(evalTime,resetTime);

	ARM_GP_Vector F0(NB_X_VALUES*nbStates);
	ARM_GP_Vector Mu0(2*NB_X_VALUES*nbStates);
	ARM_GP_Vector stdDev(nbStates),newStrikes(nbStates),cvxFwdSpread(nbStates);
	ARM_IntVector statusITM(nbStates);

	size_t i1,i2,i01,i02,i11,i12,i21,i22,offsetF=0,offsetSL=0;
	double S1,S2,fwdSpread,volS11,volS21,volS12,volS22,F10,F20,Mu110,Mu120,Mu220;
	double stdFwdFlow,maxDeviation;
	size_t nbDeepITM=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offsetF+=NB_X_VALUES,offsetSL+=2*NB_X_VALUES)
	{
		S1=coeffLong * (*longSwapRateValue)[stateIdx];
		S2=coeffShort * (*shortSwapRateValue)[stateIdx];
		fwdSpread = S1 - S2;
		newStrikes[stateIdx] = strikes[stateIdx] - fwdSpread;

		/// Spread vol factors
		i1 = offsetSL+X1_VARIABLE;
		i2 = offsetSL+X2_VARIABLE;
		volS11 = S1 * (*longF0)[i1];
		volS21 = S2 * (*shortF0)[i1];
		volS12 = S1 * (*longF0)[i2];
		volS22 = S2 * (*shortF0)[i2];

		i01 = offsetF+X1_VARIABLE;
		i02 = offsetF+X2_VARIABLE;
		F0[i01] = volS11 - volS21;
		F0[i02] = volS12 - volS22;
		F10 = F0[i01] * exp1Te;
		F20 = F0[i02] * exp2Te;

		/// Drift factors under Qpay
		i11 = i1, i12 = i2;
		i1 = offsetSL+NB_X_VALUES+X1_VARIABLE;
		i2 = offsetSL+NB_X_VALUES+X2_VARIABLE;
		i21 = i1, i22 = i2;
		Mu0[i11] = volS11 * (*longF0)[i1] - volS21 * (*shortF0)[i1];
		Mu0[i12] = volS11 * (*longF0)[i2] - volS21 * (*shortF0)[i2];
		Mu0[i21] = volS12 * (*longF0)[i1] - volS22 * (*shortF0)[i1];
		Mu0[i22] = volS12 * (*longF0)[i2] - volS22 * (*shortF0)[i2];


		/// Check if computation is really necessary by computing H&W2F derministic vol datas
		/// then compare spread between convexified fwd and strike vs standard deviation
		stdDev[stateIdx] =	F10 * (	F10 * (*covars)(X1_VARIABLE,X1_VARIABLE) +
									2*F20 * (*covars)(X1_VARIABLE,X2_VARIABLE) ) +
							F20 * F20 * (*covars)(X2_VARIABLE,X2_VARIABLE);
		stdDev[stateIdx] = (stdDev[stateIdx]>0.0 ? sqrt(stdDev[stateIdx]) : 0.0);

		/// Compute convexified forward leverage spread
		Mu110 = Mu0[i11] * exp11Te;
		Mu120 = (Mu0[i12] + Mu0[i21]) * exp12Te;
		Mu220 = Mu0[i22] * exp22Te;
		cvxFwdSpread[stateIdx] = fwdSpread + Mu110 * (*covars)(X1_VARIABLE,X1_VARIABLE) +
											 Mu120 * (*covars)(X1_VARIABLE,X2_VARIABLE) +
											 Mu220 * (*covars)(X2_VARIABLE,X2_VARIABLE);

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

	/// Free memory
	delete covars;

	/// AsOf deep ITM evaluation : option is worth nothing
	if(nbStates==1 && statusITM[0] == DEEP_OTM)
		return ARM_VectorPtr(new ARM_GP_Vector(1,0.0));

	/// Compute numerical price
	values = ComputeRungeKuttaSpreadOptionPrice(evalTime,F0,Mu0,mrs1,mrs2,payTime,resetTime,newStrikes,callPut,notional,payPeriod,newStates,statusITM);

	/// Compute the reference price for correction using standard Hull & White 2F model
	/// Branching to a 2D numerical integration

	/// Discount factors at eval date will be computed through H&W SV 2F then
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

	/// Compute standard H&W2F price using the same approximation as
	/// used in stochastic volatility case then correct SV price
	double approxPrice;
	ARM_VectorPtr zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,newStates);
	double ratio = notional * payPeriod;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx)
	{
		approxPrice = ratio * (*zcPay)[stateIdx] * VanillaOption_N(cvxFwdSpread[stateIdx], stdDev[stateIdx], strikes[stateIdx], 1.0, callPut);
		(*values)[stateIdx] += (*refPrices)[stateIdx] - approxPrice;
	}

	/// Check that std H&W2F approx formula gives same price
//	refModel->SetIsApproxSOFormula(true);
//	ARM_VectorPtr checkPrices = refModel->VanillaSpreadOptionLet(
//			curveName,
//			evalTime,
//			callPut,
//			startTime,
//			endTime,
//			resetTime,
//			payTime,
//			payPeriod,
//			notional,
//			coeffLong,
//			coeffShort,
//			strikes,
//			longFloatStartTime,
//			longFloatEndTime,
//			*longFixPayTimes,
//			*longFixPayPeriods,
//			shortFloatStartTime,
//			shortFloatEndTime,
//			*shortFixPayTimes,
//			*shortFixPayPeriods,
//			newStates);
//	refModel->SetSOFormulaFlags(oldSOFormulaFlags);


//FILE* f=fopen("c:\\temp\\dumpHW2FSV.txt","a");
//fprintf(f,"HW2FApproxCF=%15.10lf\tHW2ExactCF=%15.10lf\tHW2FSVApproxCF=%15.10lf\n",
//		approxPrice,(*refPrices)[0],(*values)[0]);
//fclose(f);


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
///	Class  : ARM_HWSV2F
///	Routine: ComputeRungeKuttaSpreadOptionPrice
///	Returns: ARM_VectorPtr
///	Action : Compute spread option price by Gauss-Legendre
///			 numerical integration and Runge-Kutta
///			 Riccati equation solving
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HWSV2F::ComputeRungeKuttaSpreadOptionPrice(	double evalTime,
																const ARM_GP_Vector& F0,
																const ARM_GP_Vector& Mu0,
																double mrs1, double mrs2,
																double Tp, double Te,
																const ARM_GP_Vector& newStrikes,
																int	callPut,
																double payNotional,
																double payPeriod,
																const ARM_PricingStatesPtr& states,
																const ARM_IntVector& statusITM) const
{
	size_t modelIdx	= GetModelNb();

	///// Constant Param
	double yfTp = Tp/K_YEAR_LEN;
	double exp1Tp = exp(-mrs1*yfTp);
	double exp2Tp = exp(-mrs2*yfTp);

	/// Set SO pricing in numerical stuff
	/// Only Heston like formula available at the moment
	bool isSOFormula	= true;
	itsNumericalsSO.SetIsSOFormula(isSOFormula);
	ARM_GP_Vector numericals = itsNumericalsSO.GetFormulaParams();

	size_t fctMultiplier = 4 * (itsNumericalsSO.IsStdFormula() ? HWSVNumericals::HestonFctMultiplier : HWSVNumericals::LewisFctMultiplier);


	size_t maxNbPts = numericals[HWSVNumericals::FirstNbSteps] < numericals[HWSVNumericals::NextNbSteps]
					? numericals[HWSVNumericals::NextNbSteps] : numericals[HWSVNumericals::FirstNbSteps];
	maxNbPts = maxNbPts < numericals[HWSVNumericals::LastNbSteps] ? numericals[HWSVNumericals::LastNbSteps]: maxNbPts;

	ARM_VectorPtr stateF0(new ARM_GP_Vector(NB_X_VALUES));
	ARM_GP_MatrixPtr stateMu0(new ARM_GP_Matrix(NB_X_VALUES,NB_X_VALUES));

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
	ARM_RiccatiHWSV2F_SO riccatiSystem(itsSystemDatas,itsSystemDatasSO,itsNumericalsSO.IsStdFormula());

	double var0,ratio = payPeriod * payNotional;
	double Integral_1,Integral_2,localIntegral_1,localIntegral_2=0.0;
	double K,absK;
	double oscilSpeedLimit,oscilSpeed,firstOscilSize,nextOscilSize,lastOscilSize;
	bool isFirstSpeedLimit,isNextSpeedLimit,isLastSpeedLimit;
	double tmin,tmax,scalet,t0,pi_1,pi_2,expect;
	ARM_IntVector nbSteps;
	size_t nbPts,iterIdx,offsetF=0,offsetMu=0;
	for(stateIdx=0;stateIdx<nbStates;++stateIdx,offsetF += NB_X_VALUES,offsetMu += 2*NB_X_VALUES)
	{
		if(statusITM[stateIdx]==DEEP_OTM)
		{
			/// Deep OTM => do nothing because is worth 0
			continue;
		}

		/// V(evalTime)
		var0 = states->GetModelState(stateIdx,modelIdx+V_VARIABLE);

		/// The option is worth something
		(*stateF0)[X1_VARIABLE] = F0[offsetF+X1_VARIABLE];
		(*stateF0)[X2_VARIABLE] = F0[offsetF+X2_VARIABLE];
		(*stateMu0)(X1_VARIABLE,X1_VARIABLE) = Mu0[offsetMu+X1_VARIABLE];
		(*stateMu0)(X1_VARIABLE,X2_VARIABLE) = Mu0[offsetMu+X2_VARIABLE];
		(*stateMu0)(X2_VARIABLE,X1_VARIABLE) = Mu0[offsetMu+NB_X_VALUES+X1_VARIABLE];
		(*stateMu0)(X2_VARIABLE,X2_VARIABLE) = Mu0[offsetMu+NB_X_VALUES+X2_VARIABLE];

		InitSystemData(evalTime,stateF0,mrs1,mrs2,Te,exp1Tp,exp2Tp,maxNbPts,itsNumericalsSO.IsStdFormula(),isSOFormula,stateMu0);

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


		/// Reset computation trace
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NbPoints,0);
		itsNumericalsSO.ResetNbStepsTrace();

//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

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

		UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

		itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

		/// To check PsiX[0]=1, PsiY[0]=0, PsiTildeY[0]=0 !!
		double psiTildeZero = PsiTildeX[0];

		if(statusITM[stateIdx]==DEEP_ITM)
		{
			/// Deep ITM => compute the instrinsic value and back to next loop
			(*values)[stateIdx] = (*zcPay)[stateIdx] * ratio * callPut * (psiTildeZero-K);
			continue;
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

		UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

		itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

		itsNumericalsSO.IntegrateSystem(GL1,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
			localIntegral_2,localIntegral_1);

		Integral_1 += localIntegral_1;
		Integral_2 += localIntegral_2;

//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

		itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
		nbSteps=riccatiSystem.GetNbSteps();
		nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
		itsNumericalsSO.PushNbStepsTrace(nbSteps);


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

		UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

		itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

		itsNumericalsSO.IntegrateSystem(GL2,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
			localIntegral_2,localIntegral_1);

		Integral_1 += localIntegral_1;
		Integral_2 += localIntegral_2;

//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

		itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
		nbSteps=riccatiSystem.GetNbSteps();
		nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
		itsNumericalsSO.PushNbStepsTrace(nbSteps);


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

			UpdateUDependentSystemData(ImU,itsNumericalsSO.IsStdFormula(),isSOFormula);

			itsNumericalsSO.SolveSystem(riccatiSystem,solverVars,PsiX,PsiY,PsiTildeX,PsiTildeY,var0);

			itsNumericalsSO.IntegrateSystem(GL,ImU,scalet,K,PsiX,PsiY,PsiTildeX,PsiTildeY,
				localIntegral_2,localIntegral_1);

			Integral_1 += localIntegral_1;
			Integral_2 += localIntegral_2;

//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",localIntegral_1,localIntegral_2);
//	fprintf(f,"I=\t%15.10lf\t%15.10lf\n",Integral_1,Integral_2);

			itsNumericalsSO.AddOscilTraceNbPoints(nbPts);
			nbSteps=riccatiSystem.GetNbSteps();
			nbSteps.push_back(riccatiSystem.GetNbDerivCalls());
			itsNumericalsSO.PushNbStepsTrace(nbSteps);

			tmin = tmax;
			++iterIdx;
		}

		itsNumericalsSO.SetOscilTrace(HWSVNumericals::FirstPeriod,firstOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NextPeriod,nextOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::LastPeriod,lastOscilSize);
		itsNumericalsSO.SetOscilTrace(HWSVNumericals::NbPeriods,iterIdx);

//	fclose(f);

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
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : SO + Lewis not available in H&WSV2F");
		}

	} // for nbStates

    return values;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

