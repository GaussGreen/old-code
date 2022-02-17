/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QGM1F.cpp
 *
 *  \brief Quadratic Gaussian Model 1 factor
 *
 *	\author  JM Prie, Amine Triki
 *	\version 1.0
 *	\date July 2004
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/QGM1F.h"
#include "gpmodels/ModelParamsQGM1F.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"

/// gpcalib
#include "gpcalib/modelfitter.h"
#include "gpcalib/calibmethod.h"

/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <inst/portfolio.h>
#include <inst/swaption.h>

/// nag
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nags.h"

#include <iomanip> /// for setprecision()
#include <algorithm>
CC_USING_NS(std,pair)
#include <set>
CC_USING_NS(std,set)

CC_BEGIN_NAMESPACE( ARM )

const short BCOEF_STORED    = 1;
const short B_STORED        = 2;
const short CACOEF_STORED   = 4;
const short VAR_STORED      = 8;

const double VOL_LIMIT      = 0.000001;
const double MRS_LIMIT      = 1.0e-6;
const double SKEW_LIMIT      = 1.0e-6;
const double DELTA_LIMIT    = 1.0e-10;

const double STDDEV_RATIO           = 5.0;
const double MIN_SLOPE_STEP_ZERO	= 0.001;
const double MIN_PRICE_ZERO	        = 0.0000001;
const double RATIO_NR_ROOT		    = 0.000001;
const double RATIO_NR_MAX_MOVE		= 0.01;
const double STEP_NR_1D		        = 0.000001;
const double MAX_NR_ITER		    = 20;

const int NBPOINT_YEAR_GL   = 2;
const int NBPOINT_MIN_GL    = 2;


//#define COMPUTATION_TIME_TEST

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_FunctType
///	Routine: GetBCoef, GetB
///	Returns: 
///	Action : Compute/get BCoef, B, CACoef or B2
////////////////////////////////////////////////////
double ARM_QGM1F::ARM_FunctType::GetBCoef(double tinf, double tsup, bool isStorable)
{
    if(!isStorable)
        return ComputeBCoef(tinf,tsup);
    else if(!(itsFlag&BCOEF_STORED))
    {
        itsBCoef = ComputeBCoef(tinf,tsup);
        itsFlag |= BCOEF_STORED;
    }

    return itsBCoef;
}

double ARM_QGM1F::ARM_FunctType::GetB(double tinf, double tsup, bool isStorable)
{
    if(!isStorable)
        return ComputeB(tinf,tsup);
    else if(!(itsFlag&B_STORED))
    {
        itsB = ComputeB(tinf,tsup);
        itsFlag |= B_STORED;
    }
    return itsB;
}

double ARM_QGM1F::ARM_FunctType::GetCACoef(double tinf, double tsup, bool isStorable)
{
    if(!isStorable)
        return ComputeCACoef(tinf,tsup);
    else if(!(itsFlag&CACOEF_STORED))
    {
        itsCACoef = ComputeCACoef(tinf,tsup);
        itsFlag |= CACOEF_STORED;
    }
    return itsCACoef;
}

double ARM_QGM1F::ARM_FunctType::GetVar(double tinf, double tsup, bool isStorable)
{
    if(!isStorable)
        return ComputeVar(tinf,tsup);
    else if(!(itsFlag&VAR_STORED))
    {
        itsVar = ComputeVar(tinf,tsup);
        itsFlag |= VAR_STORED;
    }
    return itsVar;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct0
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : No root case : 
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct0::InitCst(double t,double a)
{
    t /= K_YEAR_LEN;
    double invDelta = 1.0/itsData->itsDelta;
    itsCst = t-invDelta*atan((itsData->its2Vol2*a+itsData->itsMrs)*invDelta);
}

double ARM_QGM1F::ARM_Funct0::ComputeA(double t) const
{
    t /= K_YEAR_LEN;
    return (-itsData->itsMrs + itsData->itsDelta*tan(itsData->itsDelta*(t-itsCst)))/itsData->its2Vol2;
}

double ARM_QGM1F::ARM_Funct0::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return cos(itsData->itsDelta*(tsup-itsCst))/cos(itsData->itsDelta*(tinf-itsCst));
}

double ARM_QGM1F::ARM_Funct0::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return (sin(itsData->itsDelta*(tsup-itsCst))-sin(itsData->itsDelta*(tinf-itsCst)))/
           (cos(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}

double ARM_QGM1F::ARM_Funct0::ComputeCACoef(double tinf, double tsup) const
{
    return -log(ComputeBCoef(tinf,tsup));
}

double ARM_QGM1F::ARM_Funct0::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return cos(itsData->itsDelta*(tsup-itsCst))*sin(itsData->itsDelta*(tsup-tinf))*
           0.5*itsData->its2Vol2/(cos(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct1
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 1 root case :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct1::InitCst(double t,double a)
{
    t /= K_YEAR_LEN;
	itsCst=t+1.0/(itsData->its2Vol2*a+itsData->itsMrs);
}

double ARM_QGM1F::ARM_Funct1::ComputeA(double t) const
{
    t /= K_YEAR_LEN;
    return -(itsData->itsMrs + 1.0/(t-itsCst))/itsData->its2Vol2;
}

double ARM_QGM1F::ARM_Funct1::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return (tsup - itsCst)/(tinf - itsCst);
}

double ARM_QGM1F::ARM_Funct1::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double xinf = tinf-itsCst;
    double xsup = tsup-itsCst;
    return 0.5*(xsup*xsup/xinf-xinf);
}

double ARM_QGM1F::ARM_Funct1::ComputeCACoef(double tinf, double tsup) const
{
    return -log(ComputeBCoef(tinf,tsup));
}

double ARM_QGM1F::ARM_Funct1::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return 0.5*itsData->its2Vol2*(tsup-itsCst)*(tsup-tinf)/(tinf-itsCst);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct2Sup
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the greater root :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct2Sup::InitCst(double t,double a)
{
    if(fabs(a-itsData->itsRootSup) > 2*DELTA_LIMIT*itsData->its2Vol2)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Only constant solution (=greater root) is allowed for A(t,T)");

    itsCst = itsData->itsRootSup;
}

double ARM_QGM1F::ARM_Funct2Sup::ComputeA(double t) const
{
    return itsData->itsRootSup;
}

double ARM_QGM1F::ARM_Funct2Sup::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=itsData->its2Vol2*itsData->itsRootSup+itsData->itsMrs;
    return exp(-coef*(tsup-tinf));
}

double ARM_QGM1F::ARM_Funct2Sup::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=itsData->its2Vol2*itsData->itsRootSup+itsData->itsMrs;
    return (1-exp(-coef*(tsup-tinf)))/coef;
}

double ARM_QGM1F::ARM_Funct2Sup::ComputeCACoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return itsData->its2Vol2*itsData->itsRootSup*(tsup-tinf);
}

double ARM_QGM1F::ARM_Funct2Sup::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=2*itsData->itsDelta;
    return 0.5*itsData->its2Vol2*(1-exp(-coef*(tsup-tinf)))/coef;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct2Inf
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a constant solution
///          equals to the lower root :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct2Inf::InitCst(double t,double a)
{
    if(fabs(a-itsData->itsRootInf) > 2*DELTA_LIMIT*itsData->its2Vol2)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Only constant solution (=lower root) is allowed for A(t,T)");

    itsCst = itsData->itsRootInf;
}

double ARM_QGM1F::ARM_Funct2Inf::ComputeA(double t) const
{
    return itsData->itsRootInf;
}

double ARM_QGM1F::ARM_Funct2Inf::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=itsData->its2Vol2*itsData->itsRootInf+itsData->itsMrs;
    return exp(-coef*(tsup-tinf));
}

double ARM_QGM1F::ARM_Funct2Inf::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=itsData->its2Vol2*itsData->itsRootInf+itsData->itsMrs;
    return (1-exp(-coef*(tsup-tinf)))/coef;
}

double ARM_QGM1F::ARM_Funct2Inf::ComputeCACoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return itsData->its2Vol2*itsData->itsRootInf*(tsup-tinf);
}

double ARM_QGM1F::ARM_Funct2Inf::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    double coef=2*itsData->itsDelta;
    return 0.5*itsData->its2Vol2*(exp(coef*(tsup-tinf))-1)/coef;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct2Out
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct2Out::InitCst(double t,double a)
{
    t /= K_YEAR_LEN;
    itsCst = t-atanh(-itsData->itsDelta/(itsData->its2Vol2*a+itsData->itsMrs))/itsData->itsDelta;
}

double ARM_QGM1F::ARM_Funct2Out::ComputeA(double t) const
{
    t /= K_YEAR_LEN;
    return -(itsData->itsMrs + itsData->itsDelta/tanh((itsData->itsDelta*(t-itsCst))))/itsData->its2Vol2;
}

double ARM_QGM1F::ARM_Funct2Out::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return sinh(itsData->itsDelta*(tsup-itsCst))/sinh(itsData->itsDelta*(tinf-itsCst));
}

double ARM_QGM1F::ARM_Funct2Out::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return (cosh(itsData->itsDelta*(tsup-itsCst))-cosh(itsData->itsDelta*(tinf-itsCst)))/
           (sinh(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}

double ARM_QGM1F::ARM_Funct2Out::ComputeCACoef(double tinf, double tsup) const
{
    return -log(ComputeBCoef(tinf,tsup));
}

double ARM_QGM1F::ARM_Funct2Out::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return sinh(itsData->itsDelta*(tsup-itsCst))*sinh(itsData->itsDelta*(tsup-tinf))*
           0.5*itsData->its2Vol2/(sinh(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F::ARM_Funct2In
///	Routine: ComputeA, ComputeCst
///	Returns: 
///	Action : 2 roots case and a solution outside roots :
///          solve Cst so that A(t,Cst)=a,
///          compute A(t,Cst),
///          compute BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
///          compute B(t,u)=integ(s=t->u,BCoef(t,s)ds)
////////////////////////////////////////////////////
void ARM_QGM1F::ARM_Funct2In::InitCst(double t,double a)
{
    t /= K_YEAR_LEN;
    double invDelta = 1.0/itsData->itsDelta;
    itsCst = t-invDelta*atanh(-(itsData->its2Vol2*a+itsData->itsMrs)*invDelta);
}

double ARM_QGM1F::ARM_Funct2In::ComputeA(double t) const
{
    t /= K_YEAR_LEN;
    return -(itsData->itsMrs + itsData->itsDelta*tanh((itsData->itsDelta*(t-itsCst))))/itsData->its2Vol2;
}

double ARM_QGM1F::ARM_Funct2In::ComputeBCoef(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return cosh(itsData->itsDelta*(tsup-itsCst))/cosh(itsData->itsDelta*(tinf-itsCst));
}

double ARM_QGM1F::ARM_Funct2In::ComputeB(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return (sinh(itsData->itsDelta*(tsup-itsCst))-sinh(itsData->itsDelta*(tinf-itsCst)))/
           (cosh(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}

double ARM_QGM1F::ARM_Funct2In::ComputeCACoef(double tinf, double tsup) const
{
    return -log(ComputeBCoef(tinf,tsup));
}

double ARM_QGM1F::ARM_Funct2In::ComputeVar(double tinf, double tsup) const
{
    tinf /= K_YEAR_LEN;
    tsup /= K_YEAR_LEN;
    return cosh(itsData->itsDelta*(tsup-itsCst))*sinh(itsData->itsDelta*(tsup-tinf))*
           0.5*itsData->its2Vol2/(cosh(itsData->itsDelta*(tinf-itsCst))*itsData->itsDelta);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_QGM1F::CopyNoCleanUp(const ARM_QGM1F& rhs)
{
    itsFunctionDatas    = rhs.itsFunctionDatas;
    itsFunctions        = rhs.itsFunctions;

    itsTime             = rhs.itsTime;
    itsMaturity         = rhs.itsMaturity;
    itsA                = rhs.itsA;
    itsB                = rhs.itsB;
    itsC                = rhs.itsC;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QGM1F::ARM_QGM1F(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQGM1F& params) :
ARM_PricingModelIR(zc,&params)
{
    InitFunctionDatas();
    CC_ARM_SETNAME(ARM_QGM1F_MODEL);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QGM1F::ARM_QGM1F(const ARM_QGM1F& rhs)
: ARM_PricingModelIR(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_QGM1F::~ARM_QGM1F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QGM1F& ARM_QGM1F::operator=(const ARM_QGM1F& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_QGM1F::Clone() const
{
	return new ARM_QGM1F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::InitFunctionDatas()
{
    /// Compute a merged schedule of model parameters
    const ARM_GP_Vector& schedVol  = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedMrs  = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetAbscisses();
    const ARM_GP_Vector& schedSkew = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Skew)).GetCurve()->GetAbscisses();
    ARM_GP_Vector tmpSched(schedVol.size()+schedMrs.size()),sched(tmpSched.size()+schedSkew.size());
    CC_NS(std,merge)(schedVol.begin(),schedVol.end(),schedMrs.begin(),schedMrs.end(),tmpSched.begin());
    CC_NS(std,merge)(tmpSched.begin(),tmpSched.end(),schedSkew.begin(),schedSkew.end(),sched.begin());
    ARM_GP_Vector::iterator last=CC_NS(std,unique)(sched.begin(),sched.end());
	sched.resize( last-sched.begin() );
    if(sched[0] < K_NEW_DOUBLE_TOL && sched.size() > 1)
        sched.erase(sched.begin());

    /// Store parameters for each interval ]ti, ti+1]
    /// They are stepwise right constant
    size_t schedSize = sched.size();
    itsFunctionDatas.resize(schedSize);
    double t,mrs,vol,d2Vol2,skew,delta;
    for(size_t i=0;i<schedSize;++i)
    {
        t=sched[i];
        itsFunctionDatas[i].itsTime=t;
        mrs = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->Interpolate(t);
        if(fabs(mrs) < MRS_LIMIT)
		mrs = (mrs > 0.0 ? MRS_LIMIT : -MRS_LIMIT);
		itsFunctionDatas[i].itsMrs=mrs;
        vol = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->Interpolate(t);
        if(fabs(vol) < VOL_LIMIT)
            vol = (vol > 0.0 ? VOL_LIMIT : -VOL_LIMIT);
		d2Vol2=2*vol*vol;
	    itsFunctionDatas[i].its2Vol2=d2Vol2;
        skew = ((ARM_CurveModelParam&)GetModelParams()->GetModelParam(ARM_ModelParamType::Skew)).GetCurve()->Interpolate(t);
        if(fabs(skew) < SKEW_LIMIT)
        skew = (skew > 0.0 ? SKEW_LIMIT : -SKEW_LIMIT);
		delta = mrs*mrs + skew*d2Vol2;
        if(delta > DELTA_LIMIT)
        {
            delta=sqrt(delta);
            itsFunctionDatas[i].itsNbRoots=2;
            itsFunctionDatas[i].itsRootInf=-(mrs+delta)/d2Vol2;
            itsFunctionDatas[i].itsRootSup=(-mrs+delta)/d2Vol2;
        }
        else if(delta < -DELTA_LIMIT)
        {
            delta=sqrt(-delta);
            itsFunctionDatas[i].itsNbRoots=0;
            itsFunctionDatas[i].itsRootInf=0.0;
            itsFunctionDatas[i].itsRootSup=0.0;
        }
        else
        {
            delta=0.0;
            itsFunctionDatas[i].itsNbRoots=1;
            itsFunctionDatas[i].itsRootInf=-mrs/d2Vol2;
            itsFunctionDatas[i].itsRootSup=itsFunctionDatas[i].itsRootInf;
        }
        itsFunctionDatas[i].itsDelta=delta;   
    }

    /// Erase saved functions if any
    ARM_FunctionsMapIter firstMapElem = itsFunctions.begin();
    ARM_FunctionsMapIter lastMapElemn = itsFunctions.end();
    itsFunctions.erase(firstMapElem,lastMapElemn);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_QGM1F::ARM_FunctionsMapIter ARM_QGM1F::GenerateFunction(double maturity) const
{
    /// Localise time just before maturity
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < maturity;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;

    double initCond;
    vector< ARM_FunctTypePtr > function(iFirst+1);
    double T0=maturity;
    for(int i=iFirst;i>=0;--i)
    {
        /// Compute the initial condition
        if(i==iFirst)
            initCond = 0.0;
        else
        {
            T0 = itsFunctionDatas[i].itsTime;
            initCond = function[i+1]->ComputeA(T0);
        }

        if(itsFunctionDatas[i].itsNbRoots == 2)
        {
            /// Affect the right function depending of the initial condition
            if(itsFunctionDatas[i].itsRootInf < initCond && initCond < itsFunctionDatas[i].itsRootSup)
                function[i] = ARM_FunctTypePtr(new ARM_Funct2In(&(itsFunctionDatas[i])));
            else if(initCond < itsFunctionDatas[i].itsRootInf || itsFunctionDatas[i].itsRootSup < initCond)
                function[i] = ARM_FunctTypePtr(new ARM_Funct2Out(&(itsFunctionDatas[i])));
            else if(itsFunctionDatas[i].itsRootInf == initCond)
                function[i] = ARM_FunctTypePtr(new ARM_Funct2Inf(&(itsFunctionDatas[i])));
            else
                function[i] = ARM_FunctTypePtr(new ARM_Funct2Sup(&(itsFunctionDatas[i])));
        }
        else if(itsFunctionDatas[i].itsNbRoots == 1)
            function[i] = ARM_FunctTypePtr(new ARM_Funct1(&(itsFunctionDatas[i])));
        else
            function[i] = ARM_FunctTypePtr(new ARM_Funct0(&(itsFunctionDatas[i])));

        /// Compute integration constant
        function[i]->InitCst(T0,initCond);
    }

    pair< double,vector< ARM_FunctTypePtr > > value(maturity,function);
    pair< ARM_FunctionsMapIter,bool > result = itsFunctions.insert(value);
    if(!result.second)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": can't insert a new Riccati solution in the map");

    return result.first;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: ComputeA
///	Returns: double
///	Action : Compute A(t,T) using the input function
///          A(.,T).
////////////////////////////////////////////////////
double ARM_QGM1F::ComputeA(double time,const vector< ARM_FunctTypePtr >& function) const
{
    /// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < time;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;

    /// Call the right function
    ARM_FunctTypePtr funct=function[iFirst];
    if(funct == ARM_FunctTypePtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no function for this interval");

    return funct->ComputeA(time);
}

double ARM_QGM1F::A(double t, double T) const
{
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

    double AtT=ComputeA(t,found->second);

#ifdef COMPUTATION_TIME_TEST
    itsFunctions.erase(found);
#endif

    return AtT;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: ComputeB
///	Returns: double
///	Action : Compute B(t,T)
////////////////////////////////////////////////////
double ARM_QGM1F::ComputeB(double time,double maturity,vector< ARM_FunctTypePtr >& function) const
{
    if(time > maturity + K_NEW_DOUBLE_TOL)
    {
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": try to compute B(t,T) function with t > T");
    }
    else if(maturity - K_NEW_DOUBLE_TOL <= time)
        return 0.0;

    /// Localise indexes w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst,iLast;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < time;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;
    for(iLast=iFirst;iLast<nbFunc && itsFunctionDatas[iLast].itsTime + K_NEW_DOUBLE_TOL < maturity;++iLast);
    if(iLast>=nbFunc) iLast=nbFunc-1;

    /// Loop over intervals to compute B(t,T)
    double tinf=time,tsup;
    double BCoef=1.0,B=0.0;
    bool isStorable=false;
    for(int i=iFirst;i<iLast;++i)
    { 
        tsup=itsFunctionDatas[i].itsTime;
        B += BCoef * function[i]->GetB(tinf,tsup,isStorable);
        BCoef *= function[i]->GetBCoef(tinf,tsup,isStorable);
        tinf=tsup;
        isStorable = true;
    }
    if(maturity > tinf + K_NEW_DOUBLE_TOL)
        B += BCoef * function[i]->GetB(tinf,maturity,false);

    return B;
}

double ARM_QGM1F::B(double t, double T) const
{
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

    double BtT=ComputeB(t,T,found->second);

#ifdef COMPUTATION_TIME_TEST
    itsFunctions.erase(found);
#endif

    return BtT;

}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_QGM1F::IntegrateB2(double tinf, double tsup,
                             double time,vector< ARM_FunctTypePtr >& timeFunction,
                             double maturity,vector< ARM_FunctTypePtr >& maturityFunction) const
{
    int nbQuadPoints = static_cast<int>(floor((tsup-tinf)/K_YEAR_LEN * NBPOINT_YEAR_GL));
    if(nbQuadPoints < NBPOINT_MIN_GL)
        nbQuadPoints = NBPOINT_MIN_GL;

    GaussLegendre_Coefficients Quadrature(nbQuadPoints,tinf,tsup);

    double t,w,Bt,BT;
    double IntegB2=0.0;
    for(int i=0;i<nbQuadPoints;++i)
    {
        t   = Quadrature.get_point(i);
        w   = Quadrature.get_weight(i);

        Bt = ComputeB(t,time,timeFunction);
        BT = ComputeB(t,maturity,maturityFunction);
        IntegB2 += w * (BT*BT - Bt*Bt);

    }

    return IntegB2/K_YEAR_LEN;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_QGM1F::ComputeC(double time,vector< ARM_FunctTypePtr >& timeFunction,
                           double maturity,vector< ARM_FunctTypePtr >& maturityFunction) const
{
    /// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < time;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;

    /// Loop over intervals to compute CACoef
    double tinf=0.0,tsup;
    double CACoef=0.0,CBCoef=0.0;
    bool isStorable=true;
    for(int i=0;i<iFirst;++i)
    { 
        if( timeFunction[i] == ARM_FunctTypePtr(NULL) ||
            maturityFunction[i] == ARM_FunctTypePtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": no function for this interval");

        tsup=itsFunctionDatas[i].itsTime;

        /// Compute analytical A dependent part
        CACoef += timeFunction[i]->GetCACoef(tinf,tsup,isStorable) -
                  maturityFunction[i]->GetCACoef(tinf,tsup,isStorable);

        /// Compute B dependent part by a Gauss-Legendre integration
        if(tinf < tsup - K_NEW_DOUBLE_TOL)
            CBCoef += itsFunctionDatas[i].its2Vol2 *
                      IntegrateB2(tinf,tsup,time,timeFunction,maturity,maturityFunction);

        tinf    = tsup;
    }
    isStorable = (iFirst < nbFunc-1 && itsFunctionDatas[iFirst].itsTime - K_NEW_DOUBLE_TOL <= time) ? true : false;
    CACoef += timeFunction[iFirst]->GetCACoef(tinf,time,isStorable) -
              maturityFunction[iFirst]->GetCACoef(tinf,time,isStorable);

    if(tinf < time - K_NEW_DOUBLE_TOL)
        CBCoef += itsFunctionDatas[iFirst].its2Vol2*
                  IntegrateB2(tinf,time,time,timeFunction,maturity,maturityFunction);

    return 0.5*CACoef + 0.25*CBCoef;
}


double ARM_QGM1F::C(double t, double T) const
{
	ARM_FunctionsMap::iterator timeFound = itsFunctions.find(t);
    if(timeFound == itsFunctions.end())
        /// Insert the new function in the map
        timeFound = GenerateFunction(t);

    
	ARM_FunctionsMap::iterator maturityFound = itsFunctions.find(T);
    if(maturityFound == itsFunctions.end())
        /// Insert the new function in the map
        maturityFound = GenerateFunction(T);

    double CtT=ComputeC(t,timeFound->second,T,maturityFound->second);

#ifdef COMPUTATION_TIME_TEST
    itsFunctions.erase(timeFound);
    itsFunctions.erase(maturityFound);
#endif

    return CtT;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_QGM1F::IntegrateBBCoef(double tinf, double tsup,int tIdx,
            double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction) const
{
    int nbQuadPoints = static_cast<int>(floor((tsup-tinf)/K_YEAR_LEN * NBPOINT_YEAR_GL));
    if(nbQuadPoints < NBPOINT_MIN_GL)
        nbQuadPoints = NBPOINT_MIN_GL;

    GaussLegendre_Coefficients Quadrature(nbQuadPoints,tinf,tsup);

    double t,w,B,BCoef;
    double IntegBBCoef=0.0;
    for(int i=0;i<nbQuadPoints;++i)
    {
        t   = Quadrature.get_point(i);
        w   = Quadrature.get_weight(i);

        B = ComputeB(t,numeraireTime,numeraireTimeFunction);
        BCoef = numeraireTimeFunction[tIdx]->GetBCoef(t,tsup,false);
        IntegBBCoef += w * B * BCoef;
    }

    return IntegBBCoef/K_YEAR_LEN;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::ConditionalMean(double time,int timeIdx,double maturity,int maturityIdx,
    double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction,
	const ARM_PricingStatesPtr& states,double absoluteDrift) const
{
    /// Compute deterministic part using a Gauss-Legendre integration
    double tinf,tsup=maturity;
    bool isStorable=false;
    double BCoef=1.0,BBCoef=0.0;
    for(int idx=maturityIdx; idx>timeIdx;--idx)
    {
        tinf=itsFunctionDatas[idx-1].itsTime;
        BBCoef += BCoef * 0.5 * itsFunctionDatas[idx].its2Vol2 *
                  IntegrateBBCoef(tinf,tsup,idx,numeraireTime,numeraireTimeFunction);
        BCoef *= numeraireTimeFunction[idx]->GetBCoef(tinf,tsup,isStorable);
        tsup=tinf;
        isStorable=true;
    }
    if(time < tsup - K_NEW_DOUBLE_TOL)
    {
        BBCoef += BCoef * 0.5 * itsFunctionDatas[timeIdx].its2Vol2 *
                  IntegrateBBCoef(time,tsup,timeIdx,numeraireTime,numeraireTimeFunction);
        BCoef *= numeraireTimeFunction[timeIdx]->GetBCoef(time,tsup,false);
    }

    /// Compute stochastic part
    int nbStates=states->size();

    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
	bool isFuturState = (time > K_NEW_DOUBLE_TOL  && (nbStates > 1 || states->GetModelState(0,0) != 0.0));
    double x=0;
    for(int i=0;i<nbStates;++i)
    {
		if(isFuturState)
        x=states->GetModelState(i,0);/// 1 factor model
		x-=absoluteDrift;
        (*values)[i] = x*BCoef - BBCoef;
	}

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_QGM1F::Variance(double time,int timeIdx,double maturity,int maturityIdx,
    double numeraireTime,vector< ARM_FunctTypePtr >& numeraireTimeFunction) const
{
    double tinf,tsup=maturity;
    bool isStorable=false;
    double b,B2Coef=1.0,var=0.0;
    for(int i=maturityIdx; i>timeIdx;--i)
    {
        tinf=itsFunctionDatas[i-1].itsTime;
        var += B2Coef * numeraireTimeFunction[i]->GetVar(tinf,tsup,isStorable);
        b = numeraireTimeFunction[i]->GetBCoef(tinf,tsup,isStorable);
        B2Coef *= b*b;
        tsup=tinf;
        isStorable=true;
    }
    if(time < tsup - K_NEW_DOUBLE_TOL)
        var += B2Coef * numeraireTimeFunction[timeIdx]->GetVar(time,tsup,false);

    return var;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::XLaw(double time, double maturity, double numeraireTime,
                              const ARM_PricingStatesPtr& states,double absoluteDrift) const
{
	ARM_FunctionsMap::iterator numeraireTimeFound = itsFunctions.find(numeraireTime);
    if(numeraireTimeFound == itsFunctions.end())
        /// Insert the new function in the map
        numeraireTimeFound = GenerateFunction(numeraireTime);

    /// Locate indexes
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst,iLast;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < time;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;
    for(iLast=iFirst;iLast<nbFunc && itsFunctionDatas[iLast].itsTime + K_NEW_DOUBLE_TOL < maturity;++iLast);
    if(iLast>=nbFunc) iLast=nbFunc-1;

    ARM_VectorPtr means( ConditionalMean(time,iFirst,maturity,iLast,numeraireTime,numeraireTimeFound->second,states,absoluteDrift) );
    ARM_GP_Vector* lawParams = new ARM_GP_Vector(means->size()+1);
    for(int i=0;i<means->size();++i) (*lawParams)[i] = (*means)[i];
    (*lawParams)[means->size()] = Variance(time,iFirst,maturity,iLast,numeraireTime,numeraireTimeFound->second);

    return ARM_VectorPtr(lawParams);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::InitABC(const ARM_GP_Vector& times, const vector< ARM_GP_Vector >& maturities)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QGM1F::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in QGM model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the QGM model!");

        if(numeraire->GetType() != ARM_Numeraire::TerminalZc )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            ": only TerminalZc numeraire supported by QGM model at the moment!");
		
		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

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
        int nbDir = GetModelParams()->FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Pre-Processing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::PreProcessing(ARM_ModelFitter& modelFitter)
{
///After validate modelfitter, we call this function to manage correctly optimisation.
    GetModelParams()->PreProcessing(modelFitter,modelFitter.GetFactorNb());
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: Post Processing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
/// Just reinitialise function datas
    InitFunctionDatas();
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: AdviseCurrentCalib
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
    /// Just reinitialise function datas
    InitFunctionDatas();
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: AdviseCurrentCalibSecIndex
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM1F::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: AdviseBreakPointTimes
///	Returns: void
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_QGM1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio,
				  ARM_ModelParam* inputModelParam,
				  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int pfSize = portfolio->GetSize();
    ARM_GP_Vector sched;
    double prevTime = 0.0,endTime,payTime,curTime;
    ARM_Vector *endDates,*payDates,*resetDates;

    ARM_Security *sec;
    ARM_CapFloor *cap;
    ARM_Swaption *swaption;

    switch( modelParam->GetType() )
    {
	case ARM_ModelParamType::Volatility:
        {
            for(int i=0; i<pfSize; ++i)
            {
                /// Get the last expiry of the product
                sec = portfolio->GetAsset(i);
                if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
                {
					resetDates = cap->GetResetDates();
                    curTime = (*resetDates)[resetDates->size()-1] - asOfDate;
                }
				
                else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
                    curTime = swaption->GetExpiryDate().GetJulian() - asOfDate;
				
                else
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                    ": QGM1F may be calibrated only on cap/floor or swaption");
				
                if(curTime > prevTime)
                {
                    sched.push_back(curTime);
                    prevTime = curTime;
                }
            }
			modelParam->UpdateValues(&sched);
        }
        break;
		
	case ARM_ModelParamType::Skew:
		{
            for(int i=0; i<pfSize; ++i)
            {
                /// Get the very last Zc maturity of the product
                sec = portfolio->GetAsset(i);
                if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
                {
					endDates = cap->GetSwapLeg()->GetFwdRateEndDates();
                    endTime  = (*endDates)[endDates->size()-1] - asOfDate;
					payDates = cap->GetPaymentDates();
					payTime  = (*payDates)[payDates->size()-1] - asOfDate;
					curTime     = (endTime < payTime ? payTime : endTime);
                }
                else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
                {
					endTime  = swaption->GetFloatLeg()->GetEndDateNA().GetJulian() - asOfDate;
					payDates = swaption->GetFixedLeg()->GetFlowEndDates();
					//get the next payment date to avoid taking into account days adjustments in swaption
					if(((*payDates)[0]-swaption->GetExpiryDate().GetJulian())<15) 
						payTime  = (*payDates)[1] - asOfDate;
					else payTime  = (*payDates)[0] - asOfDate;
					
					curTime     = payTime;
                }
                else
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                    ": QGM1F may be calibrated only on cap/floor or swaption");
				
                if(curTime > prevTime)
                {
                    sched.push_back(curTime);
                    prevTime = curTime;
                }
            }
			modelParam->UpdateValues(&sched);
			
        }
        break;
		
	case ARM_ModelParamType::MeanReversion:
        {
			for(int i=0; i<pfSize; ++i)
			{
				/// Get the very last Zc maturity of the product
				sec = portfolio->GetAsset(i);
				if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
				{
					endDates = cap->GetSwapLeg()->GetFwdRateEndDates();
					endTime  = (*endDates)[endDates->size()-1] - asOfDate;
					payDates = cap->GetPaymentDates();
				}
				else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
				{
					endTime  = swaption->GetFloatLeg()->GetEndDateNA().GetJulian() - asOfDate;
					payDates = swaption->GetFixedLeg()->GetFlowEndDates();
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					": QGM1F may be calibrated only on cap/floor or swaption");
				
				payTime  = (*payDates)[payDates->size()-1] - asOfDate;
				
				curTime     = (endTime < payTime ? payTime : endTime);
				if(curTime > prevTime)
				{
					sched.push_back(curTime);
					prevTime = curTime;
				}
				modelParam->UpdateValues(&sched);
			}
        }
        break;
        
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            "Unknown type : QGM1F model only supports mean reversion, volatility and skew" );
    }

}



////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_QGM1F::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	double T	= GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);

    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0], nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(nbSteps-1,1,0.0) );
	absoluteDrifts	= ARM_GP_MatrixPtr( NULL );
	int index=0;

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];

		// QGM1F => a single state variable
		(*relativeDrifts)(i,0) = ComputeLocalDrift(step,nextStep,index,found->second);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F  //to test
///	Routine: ComputeLocalDrift
///	Returns: the Local Drift between step and next step
///         
///	Action : Compute local drift using ComputeBCoef
///          BCoef(t,u)=exp{-integ(s=t->u, 2vol2*A(s,Cst)+mrs)ds}
////////////////////////////////////////////////////

double ARM_QGM1F::ComputeLocalDrift(double t, double u, int lastIndex, const vector< ARM_FunctTypePtr >& function) const
{
	/// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst,iLast;
    for(iFirst=lastIndex;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < t;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;
	for(iLast=iFirst;iLast<nbFunc && itsFunctionDatas[iLast].itsTime + K_NEW_DOUBLE_TOL < u;++iLast);
    if(iLast>=nbFunc) iLast=nbFunc-1;
	double tinf=t,tsup;
	double result=1;
    
	/// Call the right function
    for(int i=iFirst; i<iLast; ++i)
    { 
        tsup=itsFunctionDatas[i].itsTime;
        result*=function[i]->ComputeBCoef(tinf,tsup);
		tinf=tsup;
    }
	if(u > tinf + K_NEW_DOUBLE_TOL)
        result *=function[i]->ComputeBCoef(tinf,u);
	lastIndex=iLast;
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: ComputeLocalVariance
///	Returns: 
///	Action : Compute local and Variance using ComputeVar
///          
////////////////////////////////////////////////////
double ARM_QGM1F::ComputeLocalVariance(double t, double u, int lastIndex, const vector< ARM_FunctTypePtr >& function) const
{
	/// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst,iLast;
    for(iFirst=lastIndex;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < t;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;
	for(iLast=iFirst;iLast<nbFunc && itsFunctionDatas[iLast].itsTime + K_NEW_DOUBLE_TOL < u;++iLast);
    if(iLast>=nbFunc) iLast=nbFunc-1;
	double tinf,tsup=u;
	double b,B2Coef=1.0,var=0.0;

    for(int i=iLast; i>iFirst;--i)
    {
        tinf=itsFunctionDatas[i-1].itsTime;
        var += B2Coef * function[i]->ComputeVar(tinf,tsup);
        b = function[i]->ComputeBCoef(tinf,tsup);
        B2Coef *= b*b;
        tsup=tinf;

    }
    if(t < tsup - K_NEW_DOUBLE_TOL)
    var += B2Coef * function[iFirst]->ComputeVar(t,tsup);

	lastIndex=iLast;
    return var;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F  
///	Routine: ComputeDeterministicDrift
///	Returns: the Deterministic Drift between 0 and maturity
///         
///	Action : Compute Deterministic Drift using Gauss Legendre Integration
////////////////////////////////////////////////////

double ARM_QGM1F::ComputeDeterministicDrift(double evalTime, double maturity, double T,vector< ARM_FunctTypePtr >& function) const
{
	/// Localise index w.r.t. schedule
    size_t nbFunc=itsFunctionDatas.size();
    int iFirst,iLast;
    for(iFirst=0;iFirst<nbFunc && itsFunctionDatas[iFirst].itsTime + K_NEW_DOUBLE_TOL < evalTime;++iFirst);
    if(iFirst>=nbFunc) iFirst=nbFunc-1;
	for(iLast=0;iLast<nbFunc && itsFunctionDatas[iLast].itsTime + K_NEW_DOUBLE_TOL < maturity;++iLast);
    if(iLast>=nbFunc) iLast=nbFunc-1;

	/// Compute deterministic part using a Gauss-Legendre integration
    double tinf,tsup=maturity;
    bool isStorable=false;
    double BCoef=1.0,BBCoef=0.0;
    for(int idx=iLast; idx>iFirst;--idx)
    {
        tinf=itsFunctionDatas[idx-1].itsTime;
        BBCoef += BCoef * 0.5 * itsFunctionDatas[idx].its2Vol2 *
                  IntegrateBBCoef(tinf,tsup,idx,T,function);
        BCoef *= function[idx]->GetBCoef(tinf,tsup,isStorable);
        tsup=tinf;
        isStorable=true;
    }
	if(evalTime < tsup - K_NEW_DOUBLE_TOL)
    {
        BBCoef += BCoef * 0.5 * itsFunctionDatas[iFirst].its2Vol2 *
                  IntegrateBBCoef(evalTime,tsup,iFirst,T,function);
        BCoef *= function[iFirst]->GetBCoef(evalTime,tsup,false);
    }
	return BBCoef;
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM1F::ModelStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;
#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelNb+1));
	int Index=0;

	// There is no variance in the model

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,1.0);

		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F //Tested
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local  variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM1F::NumMethodStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;
#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelNb+1));
	int Index=0;

	// All the variance is in the numerical method

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,ComputeLocalVariance(step,nextStep,Index,found->second));

		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_QGM1F::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM1F::NumMethodStateGlobalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& variances) const
{
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
    found = GenerateFunction(T);

	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex1	= (nbSteps-1)*modelNb;
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	variances.resize(nbSteps*(modelNb+1));

    variances[offsetIndex2+0]=new ARM_GP_TriangularMatrix(1,0.0);

	int Index=0,maturityIdx;
	size_t nbFunc=itsFunctionDatas.size();

    for(size_t i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
		for(maturityIdx=Index;maturityIdx<nbFunc && itsFunctionDatas[maturityIdx].itsTime + K_NEW_DOUBLE_TOL < nextStep;++maturityIdx)
			;
		
		if(maturityIdx>=nbFunc) 
			maturityIdx=nbFunc-1;
		
		/// [i+1] => variance from 0 -> ti+1
        /// we can't sum up local variance !
        variances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(1,Variance(0.0,0,nextStep,maturityIdx,T,found->second));

        step=nextStep;
		Index=maturityIdx;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

    // CurveName is not used (only for multi-currencies model)
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);

	if( evalTime <= K_NEW_DOUBLE_TOL ||
        states == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,zcT) );
    }

    int i,nbStates=states->size();
    if(    GetNumeraire() == ARM_NumerairePtr(NULL) 
		|| GetNumeraire()->GetType() == ARM_Numeraire::Cash)
    {
        if(evalTime < maturityTime)
        {
            double zctT = zcT/ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);

            /// Compute zero-coupon functions and save them
            double AtT=A(evalTime,maturityTime);
            double BtT=B(evalTime,maturityTime);
            double CtT=C(evalTime,maturityTime);

            ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
            double x;
            for(i=0;i<nbStates;++i)
            {
                x=states->GetModelState(i,0); /// 1 factor model
                (*values)[i]=zctT*exp(-x*(AtT*x + BtT) - CtT);
            }
        }
        else
            return ARM_VectorPtr(new ARM_GP_Vector(nbStates,1.0));
    }
    else
    {
        /// TerminalZc cases
		if(evalTime < maturityTime)
        {
            double zctT = zcT/ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);

            /// Compute zero-coupon functions and save them
            double AtT=A(evalTime,maturityTime);
            double BtT=B(evalTime,maturityTime);
            double CtT=C(evalTime,maturityTime);

            ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
            double x;
			double drift=ComputeDeterministicDrift(0,evalTime,T,found->second);
            for(i=0;i<nbStates;++i)
            {
                x=states->GetModelState(i,0); /// 1 factor model
				x-=drift;
                (*values)[i]=zctT*exp(-x*(AtT*x + BtT) - CtT);
            }
			return values;
        }
        else
            return ARM_VectorPtr(new ARM_GP_Vector(nbStates,1.0));
    }

    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
    for(i=0;i<nbStates;++i)
        (*values)[i]=0.0;

    return values;
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
	
	/// handle libor with no convexity correction
	if (	(	evalTime <= K_NEW_DOUBLE_TOL 
	 	||	states   == ARM_PricingStatesPtr(NULL) ) 
		&&  (-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5) )
    {
		ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		double libor = ((*ZcStart)[0]/(*ZcEnd)[0]-1.0)/period;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,libor) );
    }

    int i,nbStates=states->size();
    ARM_VectorPtr values( new ARM_GP_Vector(nbStates) );
    if(-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5)
    {
		ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
        /// No convexity
        for(i=0;i<nbStates;++i)
            (*values)[i]=((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
    }
    else
    {
		ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
        /// Specialised version to compute payment convexity    
		ARM_FunctionsMap::iterator found = itsFunctions.find(payTime);
	    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(payTime);

		double T=GetNumeraire()->GetMaturity();
		ARM_FunctionsMap::iterator found2 = itsFunctions.find(T);
	    if(found2 == itsFunctions.end())
        /// Insert the new function in the map
        found2 = GenerateFunction(T);


		double drift= ComputeDeterministicDrift(0,evalTime,T,found2->second);
		/// Compute X(Tstart) conditional mean (to X(t)) & variance under QTp probability
		ARM_VectorPtr lawParams( XLaw(evalTime,fwdResetTime,payTime,states,drift) );
		
	    // CurveName is not used (only for multi-currencies model)
	    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	    double zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
	    double zcTs=ZcCurve->DiscountPrice(fwdStartTime/K_YEAR_LEN);
	    double zcTe=ZcCurve->DiscountPrice(fwdEndTime/K_YEAR_LEN);
	    double zc0TsTe = zcTs/zcTe;
	    double zcTp = zcTe;
	    if(fwdEndTime != payTime)
	    zcTp = ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);
	    double zc0tTp = zcTp/zct;

		/// Compute functions for Zc maturing at forward start date
		double ATrTe=A(fwdResetTime,fwdEndTime);
		double BTrTe=B(fwdResetTime,fwdEndTime);
		double CTrTe=C(fwdResetTime,fwdEndTime);
		
		/// Compute functions for Zc maturing at forward start date
		double ATrTs=A(fwdResetTime,fwdStartTime);
		double BTrTs=B(fwdResetTime,fwdStartTime);
		double CTrTs=C(fwdResetTime,fwdStartTime);

		
		double var      = (*lawParams)[nbStates];
		/// To be consistent with Caplet Pricing
		double A = ATrTs - ATrTe;
		double B = BTrTs - BTrTe;
		double C = CTrTs - CTrTe;

		double stdDev       = sqrt(var);
		double stdDevMax    = STDDEV_RATIO * stdDev;
		double invStdDev    = 1.0/stdDev;
		double invVar       = 1.0/var;
		double inv2Var      = 0.5*invVar;
		
		double AA       = 2*(A + inv2Var);
		bool isANotNull = (A < -K_NEW_DOUBLE_TOL ||  K_NEW_DOUBLE_TOL < A);
		if(AA<0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Can't price the Libor with the current parameters");
		if(isANotNull)
		{
			double sAA		= sqrt(AA);
			double invSsAA	= invStdDev/sAA;
			
			for(i=0;i<nbStates;++i)
			{
				double mean = (*lawParams)[i];
				double x=states->GetModelState(i,0); /// 1 factor model
				double BB = mean*invVar - B;
				double shift = BB/AA;
				double payConvex = invSsAA*exp(0.5*shift*BB - C - mean*mean*inv2Var);
				double fwd = zc0TsTe*payConvex;
				(*values)[i] = ( fwd - 1 )/(period);
			}
		}
		else
		{
			if(-K_NEW_DOUBLE_TOL <= B && B <= K_NEW_DOUBLE_TOL)
            {
				for(i=0;i<nbStates;++i)
				/// Always exercised
				(*values)[i] = ( zc0TsTe*exp(-C) - 1 )/period;
            }
            else
            {
				for(i=0;i<nbStates;++i)
				{
					/// Exercise boundary with one bound
					double mean = (*lawParams)[i];
					double BB = mean*invVar - B;
					double shift = BB*var;
					double fwd = zc0TsTe*exp(0.5*shift*BB - C - mean*mean*inv2Var);
					
					(*values)[i] = ( fwd- 1)/period;
				}
			} // if B != 0			
		}
    }

    return values;
}


///////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QGM1F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
// FIXMEFRED: mig.vc8 (25/05/2007 15:27:39):cast
	return static_cast<ARM_PricingStatesPtr>(new ARM_PricingStates(bucketSize,1,0,1));
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_QGM1F::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new ARM_GP_Vector(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////

void ARM_QGM1F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
	const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb= FactorCount();
	size_t statesNb = states->size();
	double currentState,stdDev;
	size_t modelNb	= GetModelNb();
	
	for( size_t i=0;i<statesNb; ++i )
	{
		for( size_t j=0;  j<factorsNb; ++j )
		{
			currentState = 0.0;
			for( size_t k =0; k<=j; ++k )
			{
				stdDev			= GetLocalMatrixElemWithModelNb(localStdDev,timeIndex,modelNb,j,k);
				double gaussian = states->GetNumMethodState(i,modelNb+k);
				currentState  += stdDev*gaussian;
			}
			states->SetModelState(i,j+modelNb,currentState);
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::VanillaCaplet(
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
	/// Handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:27:49):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

    int nbStates=states->size();

#if defined(__GP_STRICT_VALIDATION)
    if(nbStates != strikesPerState.size())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": inconsistency between model & strikes states");
#endif
	
    double amount=payNotional;
    if(period != fwdPeriod)
        amount *= period/fwdPeriod;
	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	ARM_PricingStatesPtr NullStates = ARM_PricingStatesPtr(NULL);
    // CurveName is not used (only for multi-currencies model)
	ARM_VectorPtr zc0PayVector =GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,NullStates);
	ARM_VectorPtr zc0EvalVector =GetDiscountFunctor()->DiscountFactor(curveName,evalTime,evalTime,NullStates);
	ARM_VectorPtr zc0TsVector =GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,NullStates);
	ARM_VectorPtr zc0TeVector =GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,NullStates);
	double zc0TsTe = (*zc0TsVector)[0]/(*zc0TeVector)[0];
	double zc0tTp = amount*(*zc0PayVector)[0]/(*zc0EvalVector)[0];

    /// Compute functions for Zc maturing at forward start date
    double ATrTs=A(fwdResetTime,fwdStartTime);
    double BTrTs=B(fwdResetTime,fwdStartTime);
    double CTrTs=C(fwdResetTime,fwdStartTime);

    /// Compute functions for Zc maturing at forward end date
    double ATrTe=A(fwdResetTime,fwdEndTime);
    double BTrTe=B(fwdResetTime,fwdEndTime);
    double CTrTe=C(fwdResetTime,fwdEndTime);

    double AtTp=0.0,BtTp=0.0,CtTp=0.0;
    if(evalTime > K_NEW_DOUBLE_TOL)
    {
        /// Compute functions for Zc maturing at pay date
        AtTp=A(evalTime,payTime);
        BtTp=B(evalTime,payTime);
        CtTp=C(evalTime,payTime);
    }
    else
    {
    }

    double A = ATrTs-ATrTe;
    bool isANotNull = (A < -K_NEW_DOUBLE_TOL ||  K_NEW_DOUBLE_TOL < A);
    double A4=0.0,inv2A=0.0;
    if(isANotNull)
    {
        A4      = 4*A;
        inv2A   = 0.5/A;
    }
    double B    = BTrTs-BTrTe;
    double B2   = B*B;
    double C    = CTrTs-CTrTe;

	//Price Caplet under different numeraire
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(T);

	bool isFuturState = (evalTime > K_NEW_DOUBLE_TOL  && (nbStates > 1 || states->GetModelState(0,0) != 0.0));
	double drift=ComputeDeterministicDrift(0,evalTime,T,found->second);


	/// Compute X(Tr) conditional mean (to X(t)) & variance under QTp probability
    ARM_VectorPtr lawParams( XLaw(evalTime,fwdResetTime,payTime,states,drift) );

#if defined(__GP_STRICT_VALIDATION)
    if(nbStates+1 != lawParams->size())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": inconsistency between model & OU law parameters");
#endif

	double x=0;
    double BB,D,delta,K,x1,x2,mean,shift,fwd,zctTp;
    double rmin,rmax,rminmin,rmaxmax;

    double var      = (*lawParams)[nbStates];
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));


	///No Variance or Reset Time =Eval Time
	if(var < K_NEW_DOUBLE_TOL)
	{
		for(size_t i=0;i<nbStates;++i)
		{
			if(isFuturState)
			{
				x=states->GetModelState(i,0);// 1 factor model
				x-=drift;
			}
			else x=0;


			/// Zc(t,Tp)
		    if(x != 0.0)
		        zctTp = zc0tTp*exp(-x*(AtTp*x+BtTp)-CtTp);
		    else
		        zctTp = zc0tTp;

			fwd = zc0TsTe*exp(-x*(A*x+B)-C);
			K = 1 + fwdPeriod*strikesPerState[i];
			double result=zctTp*(fwd-K);
			if(capFloor == K_CAP)
				(*values)[i] = (result>0)?result:0;
			else
				(*values)[i] = (-result>0)?-result	:0;
		}
		return values;
	}
    double stdDev       = sqrt(var);
    double stdDevMax    = STDDEV_RATIO * stdDev;
    double invStdDev    = 1.0/stdDev;
    double invVar       = 1.0/var;
    double inv2Var      = 0.5*invVar;

    double AA       = 2*(A + inv2Var);
    double sAA		= (AA > 0 && isANotNull ? sqrt(AA) : 0.0);

    /// In case of sAA = 0.0, specific code is done
    double invSsAA  = 1.0;
    if(sAA > 0.0)
	    invSsAA	= invStdDev/sAA;


    /// Loop over X states at evalTime
    for(size_t i=0;i<nbStates;++i)
    {
        K = 1 + fwdPeriod*strikesPerState[i];

        mean = (*lawParams)[i];
		
		/// handle the case of negative strike
		if(K<=K_NEW_DOUBLE_TOL)
		{
			
			if(isANotNull)
			{
				if(AA <= 0)
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					": can't price the caplet with the model parameters");
				BB = mean*invVar - B;
				shift = BB/AA;
				fwd = zc0TsTe*invSsAA*exp(0.5*shift*BB - C - mean*mean*inv2Var);
				(*values)[i] = (capFloor == K_CAP)? zc0tTp*(fwd-K): 0.0;
			}	
			else
				(*values)[i] = (capFloor == K_CAP) ? (zctTp * ( zc0TsTe*exp(-C) - K )) : 0.0;
		}
		
		else
		{
			D = C + log(K/zc0TsTe);
				
			rminmin = mean - stdDevMax;
			rmaxmax = mean + stdDevMax;
			
			if(isFuturState)
				x=states->GetModelState(i,0); /// 1 factor model
			/// add the absolute drift
			x-=drift;
			
			/// Zc(t,Tp)
			if(x != 0.0)
				zctTp = zc0tTp*exp(-x*(AtTp*x+BtTp)-CtTp);
			else
				zctTp = zc0tTp;
			
			if(isANotNull)
			{
				if(AA <= 0)
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					": can't price the caplet with the model parameters");
				BB = mean*invVar - B;
				shift = BB/AA;
				fwd = zc0TsTe*invSsAA*exp(0.5*shift*BB - C - mean*mean*inv2Var);
				
				/// Exercise boundary with two bounds
				delta = B2 - A4*D;
				
				if(delta <= 0.0)
				{
					if(D < 0 && capFloor == K_CAP)
						/// Always exercised (inside the roots rminmin rmaxmax
						(*values)[i] = zctTp*(fwd-K);
					else if(D > 0 && capFloor == K_FLOOR)
						/// Always exercised
						(*values)[i] = zctTp*(K-fwd);
					else
						/// Never exercised
						(*values)[i] = 0.0;
				}
				else
				{
					
					delta = sqrt(delta);
					x1 = -(B+delta)*inv2A;
					x2 = (-B+delta)*inv2A;
					if(x1 < x2)
					{
						rmin=x1;
						rmax=x2;
					}
					else
					{
						rmin=x2;
						rmax=x1;
					}
					if((A > 0 && capFloor==K_CAP) || (A < 0 && capFloor==K_FLOOR))
					{
						/// Exercise occurs inside roots  [rmin,rmax]
						if(rminmin <= rmin)
						{
							if(rmax <= rmaxmax)
								(*values)[i] = zctTp*( fwd*( cdfNormal(sAA*(rmax-shift)) - cdfNormal(sAA*(rmin-shift)) )
								- K*(cdfNormal(invStdDev*(rmax-mean)) - cdfNormal(invStdDev*(rmin-mean)) ) );
							else
								(*values)[i] = zctTp*( fwd*cdfNormal(-sAA*(rmin-shift))
								- K*cdfNormal(-invStdDev*(rmin-mean)) );
						}
						else
							(*values)[i] = zctTp*( fwd*cdfNormal(sAA*(rmax-shift))
							- K*cdfNormal(invStdDev*(rmax-mean)) );
					}
					else
					{
						/// Exercise occurs outside roots ]-inf,rmin] U [rmax,+inf[
						if(rminmin <= rmin)
						{
							if(rmax <= rmaxmax)
								(*values)[i] = zctTp*( fwd*(cdfNormal(sAA*(rmin-shift)) + cdfNormal(-sAA*(rmax-shift)) )
								- K*(cdfNormal(invStdDev*(rmin-mean)) + cdfNormal(-invStdDev*(rmax-mean)) ) );
							else
								(*values)[i] = zctTp*( fwd*cdfNormal(sAA*(rmin-shift))
								- K*cdfNormal(invStdDev*(rmin-mean)) );
						}
						else
							(*values)[i] = zctTp*( fwd*cdfNormal(-sAA*(rmax-shift))
							- K*cdfNormal(-invStdDev*(rmax-mean)) );
					}
					if(capFloor == K_FLOOR)
						(*values)[i] = -(*values)[i];
					
				}
			} // if A != 0
			else
			{
				if(-K_NEW_DOUBLE_TOL <= B && B <= K_NEW_DOUBLE_TOL)
				{
					if(D < 0 && capFloor == K_CAP)
						/// Always exercised
						(*values)[i] = zctTp * ( zc0TsTe*exp(-C) - K );
					else if(D > 0 && capFloor == K_FLOOR)
						/// Always exercised
						(*values)[i] = zctTp * ( K - zc0TsTe*exp(-C) );
					else
						/// Never exercised
						(*values)[i] = 0.0;
				}
				else
				{
					/// Exercise boundary with one bound
					BB = mean*invVar - B;
					shift = BB*var;
					fwd = zc0TsTe*exp(0.5*shift*BB - C - mean*mean*inv2Var);
					rmin = - D/B;
					if(B>0)
					{
						(*values)[i] = capFloor * zctTp*( fwd*cdfNormal(capFloor*invStdDev*(rmin-shift))
							- K*cdfNormal(capFloor*invStdDev*(rmin-mean)) );
					}
					else if (B<0)
					{
						(*values)[i] = capFloor * zctTp*( fwd*cdfNormal(-capFloor*invStdDev*(rmin-shift))
							- K*cdfNormal(-capFloor*invStdDev*(rmin-mean)) );
					}
				} // if B != 0
			} // if A = 0
		}
    }

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: CapletImpliedVol
///	Returns:
///	Action : For test only
////////////////////////////////////////////////////
double ARM_QGM1F::CapletImpliedVol(
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		double strike,
        int capFloor,
		double price) const
{
    double amount=period*payNotional;

    // CurveName is not used (only for multi-currencies model)
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcTs=ZcCurve->DiscountPrice(fwdStartTime/K_YEAR_LEN);
    double zcTe=ZcCurve->DiscountPrice(fwdEndTime/K_YEAR_LEN);
    double fwd = (zcTs/zcTe-1.0)/fwdPeriod;
    double zcTp = zcTe;
    if(fwdEndTime != payTime)
        zcTp = ZcCurve->DiscountPrice(payTime/K_YEAR_LEN);

    price /= (amount*zcTp);

    double t=fwdResetTime/K_YEAR_LEN;
    double x=0.1*sqrt(t);
    double dx=1.0,dxmax,fx,dfx,h=0.00001;
    int maxIter=20;
    for(int i=0;i<maxIter;++i)
    {
        fx = BlackSholes_Formula(fwd,x,1.0,strike,capFloor)-price;
        dxmax = fabs(x);
        dxmax = (dxmax > 0.0001 ? 0.25*dxmax : 0.000025);
        if(fabs(fx) <= K_NEW_DOUBLE_TOL || fabs(dx) <= K_NEW_DOUBLE_TOL)
            break;
        dfx = (fx-BlackSholes_Formula(fwd,x-h,1.0,strike,capFloor)+price)/h;
        dx = -fx/dfx;
        if(dx > dxmax) x += dxmax;
        else if(dx < -dxmax) x -= dxmax;
        else x += dx;
    }

    return i==maxIter ? 100.0 : x/sqrt(t);
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::VanillaDigital(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:27:59):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

    double amount=payNotional*period;

    int nbStates=states->size();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));


    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: SwaptionPayOff
///	Returns:
///	Action :
////////////////////////////////////////////////////
double ARM_QGM1F::SwaptionRecPayoff(double x,
                                    const ARM_GP_Vector& ATTpT0,
                                    const ARM_GP_Vector& BTTpT0,
                                    const ARM_GP_Vector& CTTpT0,
                                    const ARM_GP_Vector& zc0TpT0,
                                    const ARM_GP_Vector& zcCoef,
                                    bool isSynchroEnd) const
{
    int i,nbFlows = zcCoef.size();
    double zctTpT0,price = 0.0;
    /// Fixed flows
	for(i=0;i<nbFlows;++i)
    {
        zctTpT0 = zc0TpT0[i] * exp(-x*(ATTpT0[i]*x+BTTpT0[i])-CTTpT0[i]);
        price += zcCoef[i]*zctTpT0;
    }

    /// Nominal first & last flows
    if(isSynchroEnd)
        price += zctTpT0 - 1.0;
    else
        price += zc0TpT0[nbFlows] * exp(-x*(ATTpT0[nbFlows]*x+BTTpT0[nbFlows])-CTTpT0[nbFlows]) - 1.0;

    return price;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: LocateExerciseSignedZeros
///	Returns:
///	Action :
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_QGM1F::LocateExerciseSignedZeros(double firstRoot,double endRoot,
                                                 const ARM_GP_Vector& ATTpT0,
                                                 const ARM_GP_Vector& BTTpT0,
                                                 const ARM_GP_Vector& CTTpT0,
                                                 const ARM_GP_Vector& zc0TpT0,
                                                 const ARM_GP_Vector& zcCoef,
                                                 bool isSynchroEnd, int payRec,
                                                 double Amax,double Bmax,
                                                 double minStep,double maxMove) const
{
    ARM_GP_Vector* zeros = new ARM_GP_Vector;

    double root0 = firstRoot;

    zeros->push_back(firstRoot);

    /// If isInf is true search by increasing values else by decreasing values
    bool isInf = (firstRoot < endRoot);

	double step,stepMax;
    double root,price,price1D;
	int i;
    while( isInf ==  (root0 < endRoot) )
    {
        price = SwaptionRecPayoff(root0,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);

		/// Locate zero with max slope
        if(price > 0.0)
            step = price/(((root0 < 0 ? -root0 : root0)*Amax + Bmax)*(price+1.0));
        else
            step = -price/((root0 < 0 ? -root0 : root0)*Amax + Bmax);

        root0 = root0 + (isInf ? step : -step);

        if(-MIN_SLOPE_STEP_ZERO <= step && step <= MIN_SLOPE_STEP_ZERO)
        {
            /// Switch to a NR search to quickly get the zero if any
            root = root0 + (isInf ? minStep : -minStep);
            step = 1.0;
			for(i=0;i<MAX_NR_ITER && (isInf == (root >= root0));++i)
			{
				price = SwaptionRecPayoff(root,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);

                stepMax = fabs(root);
                if(stepMax < maxMove) stepMax = maxMove;

				if(-MIN_PRICE_ZERO <= price && price <= MIN_PRICE_ZERO ||
                   -minStep <= step && step <= minStep )
				{
                    if( (isInf && (*zeros)[zeros->size()-1] < root - minStep) ||
                        (!isInf && root + minStep < (*zeros)[zeros->size()-1]) )
                    {
					    /// A new root is found
					    zeros->push_back(root);
                        root0 = root;
                    }
					break;
				}
				price1D = ( SwaptionRecPayoff(root+STEP_NR_1D,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd)
							- price)/STEP_NR_1D;
				step = -price/price1D;
				stepMax *= 0.5;
				if(step > stepMax)
					root += stepMax;
				else
				{
					if(step < -stepMax)
						root -= stepMax;
					else
						root += step;
				}
			}
            if( i==MAX_NR_ITER && (isInf == (root < firstRoot) || isInf == (root > endRoot)) )
            {
                /// No more root
                break;
            }
            root0 = root0 + (isInf ? MIN_SLOPE_STEP_ZERO : -MIN_SLOPE_STEP_ZERO);

        } // NR search
    } // root < rmaxmax

    return ARM_GP_VectorPtr(zeros);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: LocateExerciseZeros
///	Returns:
///	Action :
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_QGM1F::LocateExerciseZeros(double rminmin,double rmaxmax,
                                                 const ARM_GP_Vector& ATTpT0,
                                                 const ARM_GP_Vector& BTTpT0,
                                                 const ARM_GP_Vector& CTTpT0,
                                                 const ARM_GP_Vector& zc0TpT0,
                                                 const ARM_GP_Vector& zcCoef,
                                                 bool isSynchroEnd, int payRec,
                                                 double Amax,double Bmax,
                                                 bool& isFirstExer) const
{
    /// Save if the first region is exercised or not
    double price = SwaptionRecPayoff(rminmin,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);
    /// cheks that rminmin is not a root
	while(fabs(price)<0.0001) 
	{
		rminmin-=0.01;
		price = SwaptionRecPayoff(rminmin,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);
	}

	isFirstExer = (payRec==K_RCV && price > 0.0) || (payRec==K_PAY && price < 0.0);     
    isFirstExer = (payRec==K_RCV && price > 0.0) || (payRec==K_PAY && price < 0.0);

    double minStep = RATIO_NR_ROOT*(rmaxmax-rminmin);
	if(minStep < K_NEW_DOUBLE_TOL) minStep=K_NEW_DOUBLE_TOL;
    double maxMove = RATIO_NR_MAX_MOVE*(rmaxmax-rminmin);
	if(maxMove < K_NEW_DOUBLE_TOL) maxMove=K_NEW_DOUBLE_TOL;

    /// Locate zeros on negative side
    ARM_GP_VectorPtr zerosNeg;
    if(rminmin < 0.0)
        zerosNeg = LocateExerciseSignedZeros(rminmin,0.0,ATTpT0,BTTpT0,CTTpT0,
            zc0TpT0,zcCoef,isSynchroEnd,payRec,Amax,Bmax,minStep,maxMove);
    else
    {
        zerosNeg=ARM_GP_VectorPtr(new ARM_GP_Vector(1));
        (*zerosNeg)[0]=rminmin;
    }

    /// Locate zeros on positive side
    ARM_GP_VectorPtr zerosPos;
    if(rmaxmax > 0.0)
        zerosPos = LocateExerciseSignedZeros(rmaxmax,0.0,ATTpT0,BTTpT0,CTTpT0,
            zc0TpT0,zcCoef,isSynchroEnd,payRec,Amax,Bmax,minStep,maxMove);
    else
    {
        zerosPos=ARM_GP_VectorPtr(new ARM_GP_Vector(1));
        (*zerosPos)[0]=rmaxmax;
    }

    /// Merge both sides
    double lastZeroNeg = (*zerosNeg)[zerosNeg->size()-1];
    double maxDiffRoot=3*minStep;
    for(int i=zerosPos->size()-1;i>=0;--i)
    {
        if(lastZeroNeg + maxDiffRoot < (*zerosPos)[i])
            zerosNeg->push_back((*zerosPos)[i]);
    }

    return zerosNeg;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::VanillaSwaption(
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

	if (!isConstantNotional)
		return VariableNotionalSwaption (curveName, 
										 evalTime, 
										 swapResetTime, 
										 fixNotional, 
										 floatNotional, 
										 floatStartTime, 
										 floatEndTime, 
										 fixPayTimes, 
										 fixPayPeriods, 
										 strikesPerState, 
										 callPut, 
										 states, 
										 isConstantNotional, 
										 isConstantSpread, 
										 isConstantStrike);

	/// TO BE UPDATED
	/// Check that the notional is constant
	double swapNotional = fixNotional[0];
	if (!(isConstantNotional&&isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );

	
    if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + 
            "Swaption pricing not implemented for differents discount & fixing curves" );

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:28:14):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

    int i,j,nbStates=states->size();

#if defined(__GP_STRICT_VALIDATION)
    if(nbStates != strikesPerState.GetRowsNb())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": inconsistency between model & strikes states");
#endif

    int payRec=K_PAY;
    if(callPut == K_PUT)
        payRec=K_RCV;

	//Price the Swaption under a different Numeraire
	//Take into account the deterministic drift not included in Pricing States
	double T=GetNumeraire()->GetMaturity();
	ARM_FunctionsMap::iterator found = itsFunctions.find(T);
	if(found == itsFunctions.end())
    /// Insert the new function in the map
    found = GenerateFunction(T);

	bool isFuturState = (evalTime > K_NEW_DOUBLE_TOL  && (nbStates > 1 || states->GetModelState(0,0) != 0.0));	
	double drift=ComputeDeterministicDrift(0,evalTime,T,found->second);

	/// Compute X(T) conditional mean (to X(t)) & variance under QStart probability
    ARM_VectorPtr lawParams( XLaw(evalTime,swapResetTime,floatStartTime,states,drift) );

#if defined(__GP_STRICT_VALIDATION)
    if(nbStates+1 != lawParams->size())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": inconsistency between model & OU law parameters");
#endif

    double var          = (*lawParams)[nbStates];
    double stdDev       = sqrt(var);
    double stdDevMax    = STDDEV_RATIO * stdDev;
    double invStdDev    = 1.0/stdDev;
    double invVar       = 1.0/var;
    double inv2Var      = 0.5*invVar;

    // CurveName is not used (only for multi-currencies model)
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zct  = ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
    double zcT0 = ZcCurve->DiscountPrice(floatStartTime/K_YEAR_LEN);
    double invZcT0 = 1.0/zcT0;
    double zc0tT0  = swapNotional * zcT0/zct;

    double AtT0=0.0,BtT0=0.0,CtT0=0.0;
    if(evalTime > 0.0)
    {
        /// Compute functions for Zc functions at swap start date
        AtT0=A(evalTime,floatStartTime);
        BtT0=B(evalTime,floatStartTime);
        CtT0=C(evalTime,floatStartTime);
    }
    else
    {
#if defined(__GP_STRICT_VALIDATION)
        if(nbStates != 1 && states->GetModelState(0,0) != 0.0)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                ": a single state equal to 0 is required for an as of date pricing");
#endif
    }

	/// Compute Zc functions at swap start dates
    double ATT0=A(swapResetTime,floatStartTime);
    double BTT0=B(swapResetTime,floatStartTime);
    double CTT0=C(swapResetTime,floatStartTime);

	/// Compute Zc functions at payment dates
	size_t nbFlows=fixPayTimes.size();
	ARM_GP_Vector ATTpT0(nbFlows+1),BTTpT0(nbFlows+1),CTTpT0(nbFlows+1);
    ARM_GP_Vector zc0TpT0(nbFlows+1);
    double Aabs,Amax=-1.0e+15,Babs,Bmax=-1.0e+15;
    bool isSynchroEnd = (-K_NEW_DOUBLE_TOL <= floatEndTime - fixPayTimes[nbFlows-1] &&
                         floatEndTime - fixPayTimes[nbFlows-1] <= K_NEW_DOUBLE_TOL);
    ARM_GP_Vector AA(nbFlows+1),sAA(nbFlows+1),invSsAA(nbFlows+1),BB(nbFlows+1);
    double payTime;
	for(i=0;i<nbFlows+1;++i)
	{
        if(i==nbFlows)
        {
            if(isSynchroEnd)
            {
                /// Same end for floating & fixed legs
                ATTpT0[nbFlows] = ATTpT0[nbFlows-1];
                AA[nbFlows] = AA[nbFlows-1];
                sAA[nbFlows] = sAA[nbFlows-1];
                invSsAA[nbFlows] = invSsAA[nbFlows-1];
                BTTpT0[nbFlows] = BTTpT0[nbFlows-1];
                CTTpT0[nbFlows] = CTTpT0[nbFlows-1];
                zc0TpT0[nbFlows] = zc0TpT0[nbFlows-1];
                break;
            }
            else
                payTime = floatEndTime;
        }
        else
            payTime = fixPayTimes[i];

		ATTpT0[i]=A(swapResetTime,payTime)-ATT0;
        Aabs = fabs(ATTpT0[i]);
        if(Amax < Aabs) Amax = Aabs;
        AA[i] = 2*(ATTpT0[i] + inv2Var);

#if defined(__GP_STRICT_VALIDATION)
        if( (ATTpT0[i] < -K_NEW_DOUBLE_TOL ||  K_NEW_DOUBLE_TOL < ATTpT0[i]) && AA[i] <= 0)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                ": can't price swaption because of A(t,Ti)-A(t,T0) and 0.5/var");
#endif

        if(AA[i] > 0)
        {
            sAA[i] = sqrt(AA[i]);
            invSsAA[i] = invStdDev/sAA[i];
        }
        else
        {
            /// => ATTpT0[i] is null
            sAA[i] = invStdDev;
	        invSsAA[i]	= 1.0;
        }

		BTTpT0[i]=B(swapResetTime,payTime)-BTT0;
        Babs = fabs(BTTpT0[i]);
        if(Bmax < Babs) Bmax = Babs;

		CTTpT0[i]=C(swapResetTime,payTime)-CTT0;

        zc0TpT0[i] = ZcCurve->DiscountPrice(payTime/K_YEAR_LEN)*invZcT0;

	}
    Amax = 2.0*Amax;

    /// Loop over X states at evalTime
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
    double zctT0,price;
    ARM_GP_Vector zcCoef(nbFlows),fwd(nbFlows+1),shift(nbFlows+1);
    double mean,rminmin,rmaxmax,rmin,rmax;
    double x,y,z,v,zctTpT0;

	///No Variance or Reset Time =Eval Time
	if(var < K_NEW_DOUBLE_TOL)
	{
		for(size_t i=0;i<nbStates;++i)
		{
			price =0.;
			if(isFuturState)
			{
				x=states->GetModelState(i,0);// 1 factor model
				x-=drift;
			}
			else x=0;

			for(j=0;j<nbFlows;++j) zcCoef[j] = fixPayPeriods[j]*strikesPerState(i,j);
		    /// Zc(t,T0)
			if(x != 0.0)
				zctT0 = zc0tT0*exp(-x*(AtT0*x+BtT0)-CtT0);
			else
			    zctT0 = zc0tT0;

			price = SwaptionRecPayoff(x,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);
			if(payRec == K_PAY && price<0) price = -price;
			else if (payRec == K_RCV && price>0);
			else price=0;
			(*values)[i] = zctT0*price;
		}		
		return values;
	}




    for(i=0;i<nbStates;++i)
    {
        mean = (*lawParams)[i];
        rminmin = mean - stdDevMax;
        rmaxmax = mean + stdDevMax;

        x=states->GetModelState(i,0); /// 1 factor model
		//Add the absolute drift
		x-=drift;

        /// Zc(t,T0)
        if(x != 0.0)
            zctT0 = zc0tT0*exp(-x*(AtT0*x+BtT0)-CtT0);
        else
            zctT0 = zc0tT0;


        /// Locate zeros of the exercice function
	    for(j=0;j<nbFlows;++j) zcCoef[j] = fixPayPeriods[j]*strikesPerState(i,j);

        bool isFirstExer;
        ARM_GP_VectorPtr zeros( LocateExerciseZeros(rminmin,rmaxmax,ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd,payRec,Amax,Bmax,isFirstExer) );
        int nbRegions = zeros->size()-1;
		double localPrice;
		
        if(nbRegions >= 2)
        {
            /// Compute intermediate variables for integration inside exercise regions
            y = mean*invVar;
            z = mean*mean*inv2Var;
			for(j=0;j<nbFlows;++j)
            {
                v = y - BTTpT0[j];
                shift[j] = (AA[j] > 0 ? v/AA[j] : v*var);
                zctTpT0 = zc0TpT0[j] * invSsAA[j] * exp(0.5*shift[j]*v - CTTpT0[j] - z);
                fwd[j] = zcCoef[j] *  zctTpT0;
            }
            if(isSynchroEnd)
                fwd[nbFlows-1] += zctTpT0;
            else
            {
                v = y - BTTpT0[nbFlows];
                shift[nbFlows] = (AA[nbFlows] > 0 ? v/AA[nbFlows] : v*var);
                fwd[nbFlows] = zc0TpT0[nbFlows] * invSsAA[nbFlows] * exp(0.5*shift[nbFlows]*v - CTTpT0[nbFlows] - z);
            }
			
            price = 0.0;
            for(int h=0;h<nbRegions;++h)
            {
                rmin=(*zeros)[h];
                rmax=(*zeros)[h+1];
				localPrice = SwaptionRecPayoff(0.5*(rmin+rmax),ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);
				if((payRec==K_RCV && localPrice > 0.0) || (payRec==K_PAY && localPrice < 0.0))
				{
					if(rmin <= rminmin + K_NEW_DOUBLE_TOL)
					{
						for(j=0;j<nbFlows;++j)
							price += fwd[j] * cdfNormal(sAA[j]*(rmax-shift[j]));
						if(!isSynchroEnd)
							price += fwd[nbFlows] * cdfNormal(sAA[nbFlows]*(rmax-shift[nbFlows]));
						
						price -= cdfNormal(invStdDev*(rmax-mean));
					}
					else if(rmax < rmaxmax - K_NEW_DOUBLE_TOL)
					{
						for(j=0;j<nbFlows;++j)
							price += fwd[j] * ( cdfNormal(sAA[j]*(rmax-shift[j]))
							- cdfNormal(sAA[j]*(rmin-shift[j])) );
						if(!isSynchroEnd)
							price += fwd[nbFlows] * ( cdfNormal(sAA[nbFlows]*(rmax-shift[nbFlows]))
							- cdfNormal(sAA[nbFlows]*(rmin-shift[nbFlows])) );
						
						price -= ( cdfNormal(invStdDev*(rmax-mean)) - cdfNormal(invStdDev*(rmin-mean)) );
					}
					else
					{
						for(j=0;j<nbFlows;++j)
							price += fwd[j] * cdfNormal(-sAA[j]*(rmin-shift[j]));
						if(!isSynchroEnd)
							price += fwd[nbFlows] * cdfNormal(-sAA[nbFlows]*(rmin-shift[nbFlows]));
						
						price -= cdfNormal(-invStdDev*(rmin-mean));
					}
				}
			}
			
            if(payRec == K_PAY) price = -price;
			
        } //nbRegions >= 2

        else 
		{
			rmin=(*zeros)[0];
            rmax=(*zeros)[1];
			localPrice = SwaptionRecPayoff(0.5*(rmin+rmax),ATTpT0,BTTpT0,CTTpT0,zc0TpT0,zcCoef,isSynchroEnd);
			if((payRec==K_RCV && localPrice > 0.0) || (payRec==K_PAY && localPrice < 0.0))
			{
				price = 0.0;
				/// The only one region is exercised
				z = mean*mean*inv2Var;
				for(j=0;j<nbFlows;++j)
				{
					BB[j] = mean*invVar - BTTpT0[j];
					shift[j] = BB[j]/AA[j];
					zctTpT0= zc0TpT0[j]*invSsAA[j]*exp(0.5*shift[j]*BB[j] - CTTpT0[j]- z);
					fwd[j]=zcCoef[j]*zctTpT0;
					price +=fwd[j];
				}
				if(isSynchroEnd)
				{
					fwd[nbFlows-1] += zctTpT0;
					price+=zctTpT0;;
				}
				else
				{
					BB[nbFlows] = mean*invVar - BTTpT0[nbFlows];
					shift[nbFlows] = (AA[nbFlows] > 0 ? BB[nbFlows]/AA[nbFlows] : BB[nbFlows]*var);
					fwd[nbFlows] = zc0TpT0[nbFlows] * invSsAA[nbFlows] * exp(0.5*shift[nbFlows]*BB[nbFlows] - CTTpT0[nbFlows] - z);
					price += fwd[nbFlows] * cdfNormal(sAA[nbFlows]*(rmax-shift[nbFlows]));
				}           
				price = (price - 1);
				if(payRec == K_PAY) price = -price;
			}
			else
				/// The only one region is not exercised
				price = 0.0;

		}
        /// NPV at time t
        (*values)[i] = zctT0*price;
    }

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: SwaptionImpliedVol
///	Returns: 
///	Action : for test only
////////////////////////////////////////////////////
double ARM_QGM1F::SwaptionImpliedVol(
		double swapResetTime,
        double swapNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double strike,
        int callPut,
		double price) const
{
    // CurveName is not used (only for multi-currencies model)
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcTs=ZcCurve->DiscountPrice(floatStartTime/K_YEAR_LEN);
    double zcTe=ZcCurve->DiscountPrice(floatEndTime/K_YEAR_LEN);
    double O1=0.0;
    for(int i=0;i<fixPayTimes.size();++i)
        O1 += ZcCurve->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN)*fixPayPeriods[i];
    double fwd = (zcTs - zcTe)/O1;

    price /= (O1*swapNotional);

    double t=swapResetTime/K_YEAR_LEN;
    double x=0.1*sqrt(t);
    double dx=1.0,dxmax,fx,dfx,h=0.00001;
    int maxIter=20;
    for(i=0;i<maxIter;++i)
    {
        fx = BlackSholes_Formula(fwd,x,1.0,strike,callPut)-price;
        dxmax = fabs(x);
        dxmax = (dxmax > 0.0001 ? 0.25*dxmax : 0.000025);
        if(fabs(fx) <= K_NEW_DOUBLE_TOL || fabs(dx) <= K_NEW_DOUBLE_TOL)
            break;
        dfx = (fx-BlackSholes_Formula(fwd,x-h,1.0,strike,callPut)+price)/h;
        dx = -fx/dfx;
        if(dx > dxmax) x += dxmax;
        else if(dx < -dxmax) x -= dxmax;
        else x += dx;
    }

    return i==maxIter ? 100.0 : x/sqrt(t);
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM1F
///	Routine: VariableNotionalSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM1F::VariableNotionalSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
	/// some validations...
	if (!isConstantSpread)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F: variable spread not supported)" );

	if (!isConstantStrike)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F: variable strike not supported)" );

	if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Swaption pricing not implemented for differents discount & fixing curves" );

	if (states->size() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  multi state not supported" );

	if (strikesPerState.rows() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  strike per state not supported" );

	if (strikesPerState.cols() != fixPayTimes.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  strikesPerState: bad size" );
	
	if ( fabs(floatEndTime - fixPayTimes[fixPayTimes.size()-1])>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption: float and fixed end dates are not matching" );

	if (floatNotional.size() != fixNotional.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption: float an fix legs are required to be of same frequency." );
	
	
	/// compute coupons
	size_t i, j, size_cpn   = fixPayTimes.size();
	vector<double> coupons (size_cpn);

	for (i=0; i<size_cpn-1; i++)
		coupons[i] = floatNotional[i] + strikesPerState(0,i) * fixPayPeriods[i] * fixNotional[i] - floatNotional[i+1] ;
	
	coupons[size_cpn-1] = floatNotional[size_cpn-1] + fixPayPeriods[size_cpn-1] * strikesPerState(0,size_cpn-1) * fixNotional[size_cpn-1];

	/// bond option strike
	double strike = floatNotional[0];
			
	/// precomputations for numerical integral
	vector<double> dfs (size_cpn);
	vector<double> A_(size_cpn), B_ (size_cpn), C_ (size_cpn);
	double Astart, Bstart, Cstart;
	double startdf = GetZeroCurve()->DiscountPrice(floatStartTime/K_YEAR_LEN);

	Astart = A (swapResetTime, floatStartTime);
	Bstart = B (swapResetTime, floatStartTime);
	Cstart = C (swapResetTime, floatStartTime);

	for (i=0; i<size_cpn; i++)
	{
		dfs[i] = GetZeroCurve()->DiscountPrice (fixPayTimes[i]/K_YEAR_LEN) / startdf;
		A_[i]   = A (swapResetTime, fixPayTimes[i]) - Astart;
		B_[i]   = B (swapResetTime, fixPayTimes[i]) - Bstart;
		C_[i]   = C (swapResetTime, fixPayTimes[i]) - Cstart;
	}

	ARM_FunctionsMap::iterator found = itsFunctions.find(floatStartTime);
    if(found == itsFunctions.end())
        /// Insert the new function in the map
        found = GenerateFunction(floatStartTime);
	
	ARM_VectorPtr lawParams( XLaw(evalTime, swapResetTime, floatStartTime, states) );

	//-------------------------------
	// numerical integral
	//-------------------------------
	double EX   = (*lawParams)[0];
	double varX = (*lawParams)[1];
	double invVarX = 1.0 / varX;
	double factor  =  1.0 / sqrt(2.0 * ARM_NumericConstants::ARM_PI * varX) ;
		
	double	nstdev	= 6.0 ;
	int		nx		= 151;
	double  dX      = 2.0 * nstdev * sqrt(varX) / (nx - 1.0);
	double  Xmin    = EX -  nstdev * sqrt(varX) ;
		
	double expect (0.0);
	double X (Xmin), X2;
	double argexp, coeff, payoff, df;
	double sgn = -(double)callPut; // payer swaption = put bond option

	for (j=0; j<nx; j++)
	{
		X += dX;
		coeff = (j==0||j==nx-1) ? 0.5 * dX : dX ;
		argexp = -0.5 * (X-EX) * (X-EX) * invVarX ;
		X2 = X * X;
	
		payoff = 0.0;
		for (i=0; i<size_cpn; i++)
		{
			df = dfs[i] * exp( - A_[i] * X2 - B_[i] * X - C_[i] );
			payoff += coupons[i] * df; 
		}
		
		payoff -= strike;
		payoff *= sgn;
				
		if (payoff<0.0) payoff = 0.0;
		
		expect += coeff * payoff * factor * exp (argexp);
	}

	//-------------------------------
	// end numerical integral
	//-------------------------------

	double price = startdf * expect;
	return ARM_VectorPtr (new ARM_GP_Vector(1, price));
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routine : toString
///	Returns : string
///	Action  : object dump into a string
////////////////////////////////////////////////////
string ARM_QGM1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << indent << "QGM1F Model\n";
    os << indent << "-----------\n";

    os << ARM_PricingModel::toString(indent);

    os << indent << "Function Datas\n";
    os << indent << setw(5)     << "time";
    os << indent << setw(10)    << "mrs";
    os << indent << setw(10)    << "2Vol2";
    os << indent << setw(10)    << "Delta";
    os << indent << setw(10)    << "NbR";
    os << indent << setw(10)    << "RInf";
    os << indent << setw(10)    << "Rsup";
    os << endl;
    for(int i=0; i<itsFunctionDatas.size();++i)
    {
        os << indent << fixed << setw(5) << setprecision(0) << itsFunctionDatas[i].itsTime;
        os << indent << fixed << setw(10) << setprecision(4) << itsFunctionDatas[i].itsMrs;
        os << indent << fixed << setw(10) << setprecision(6) << itsFunctionDatas[i].its2Vol2;
        os << indent << fixed << setw(10) << setprecision(6) << itsFunctionDatas[i].itsDelta;
        os << indent << dec << setw(10) << itsFunctionDatas[i].itsNbRoots;
        os << indent << fixed << setw(10) << setprecision(3) << itsFunctionDatas[i].itsRootInf;
        os << indent << fixed << setw(10) << setprecision(3) << itsFunctionDatas[i].itsRootSup;
        os << endl;
    }

    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_QGM1F::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsQGM1F* modelParamsQGM1F = dynamic_cast<const ARM_ModelParamsQGM1F*>(&params);
	if( !modelParamsQGM1F )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsQGM" );
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routine : ValidateCalibMethod
///	Returns : void
///	Action  : validate calibMethod
////////////////////////////////////////////////////

void ARM_QGM1F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	/// checks that if the CalibParam that we are trying 
	/// to calibrate exist into Model
	ARM_ModelParamVector CalibParams = calibMethod.GetCalibParams();
	size_t sizeCalibParams = CalibParams.size();
	for(size_t i=0; i< sizeCalibParams; ++i)
	{
		if(!CalibParams[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "You are trying to validate calibMethod with modelParam NULL, please advice!" );

		if( !GetModelParams()->DoesModelParamExist(ARM_ModelParamType::ParamNb(CalibParams[i]->GetType())) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": a Calib Method should contain a model Param of pricing model!");
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM1F
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_QGM1F::VanillaSpreadOptionLet(const string& curveName,
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
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSpreadOption : unimplemented function for ARM_QGM1F Model!");
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

