/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsHW1F.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date September 2003
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/modelparamshw1f.h"

/// gpbase headers
#include "gpbase/gpmatrix.h"
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/comparisonfunctor.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparamutil.h"
#include "gpinfra/curvemodelparam.h"

/// gpmodel 
#include "gpmodels/TargetFuncHW.h"

/// gpcalib
#include "gpcalib/modelfitter.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"

/// kernel
#include <inst/portfolio.h>


#include <algorithm>
CC_USING_NS( std, sort )

#include <functional>
CC_USING_NS( std, ptr_fun )

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW1F
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1F::ARM_ModelParamsHW1F( const ARM_ModelParamsHW1F& rhs )
: ARM_ModelParamsHW(rhs),  itsVolatilityType(rhs.itsVolatilityType )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW1F::ARM_ModelParamsHW1F( const ARM_ModelParamVector& params, int volatilityType )
: ARM_ModelParamsHW(params), itsVolatilityType(volatilityType )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1F::~ARM_ModelParamsHW1F()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW1F& ARM_ModelParamsHW1F::operator=(const ARM_ModelParamsHW1F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW::operator=(rhs);
		itsVolatilityType = rhs.itsVolatilityType;
		/// Copy class attributes if any
	}
	return *this;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW1FStd
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1FStd::ARM_ModelParamsHW1FStd( const ARM_ModelParamsHW1FStd& rhs )
: ARM_ModelParamsHW1F(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW1FStd::ARM_ModelParamsHW1FStd( const ARM_ModelParamVector& params, int volatilityType )
: ARM_ModelParamsHW1F(params,volatilityType )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1FStd::~ARM_ModelParamsHW1FStd()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW1FStd& ARM_ModelParamsHW1FStd::operator=(const ARM_ModelParamsHW1FStd& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW1F::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1FStd
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsHW1FStd::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsHW1FStd\n";
    os << "----------------------\n";
    os << ARM_ModelParams::toString();
    return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1FStd
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHW1FStd::Clone() const
{
	return new ARM_ModelParamsHW1FStd(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: Variance
///	Returns: value of the variance
///	Action : Integrate sigma(t)^2*exp(scale*t)
///          between [a,b] for a stepwise right
///          constant sigma curve
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::Variance(double a,double b,double scale) const
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

    ARM_GP_Vector sigmaValues( ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
    int i,nbU = sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(a <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbU;++i)
        if(b <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double U,lastU,newU,sig;
    double timeScale=scale/K_YEAR_LEN;
    bool notTiny=(fabs(scale)>0.1*K_NEW_DOUBLE_TOL);

    double value=0.0;
    lastU=notTiny? exp(timeScale*a):a;

    /// Between a and Un(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        U=(i<nb ? sigmaTimes[i] : b);
        sig=sigmaValues[i<nbU ? i : nbU-1];
		if(fabs(sig) < ARM_ModelParamsHW::VOL_LIMIT)
        sig = (sig > 0.0 ? ARM_ModelParamsHW::VOL_LIMIT : -ARM_ModelParamsHW::VOL_LIMIT);
        newU=notTiny?exp(timeScale*U):U;
        value += sig*sig*(newU-lastU);
        lastU=newU;
    }

    return value/(notTiny ? scale : K_YEAR_LEN);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: StateLocalDrift
///	Returns: value of the drift
///	Action : Relative drift of the state variable
///          from a to b>=a
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::StateLocalDrift(double a,double b) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

    if(fabs(MRSValue)>K_NEW_DOUBLE_TOL)
        return exp(-MRSValue*(b-a)/K_YEAR_LEN);
    else
        return 1.0;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
///          (alias phi(a,b)) that is :
///          Integ(a->b) { sigma(s)^2.exp(-2.MRS.(c-s)).ds }
///          Note that Variance(a,b,scale) computes
///          the integral a -> b of sigma(s)^2.exp(scale.s).ds
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::StateLocalVariance(double a,double b,double c) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	double localVariance;

    if(fabs(MRSValue)<=K_NEW_DOUBLE_TOL)
	{
		double scale = K_NEW_DOUBLE_TOL*0.01;
		return Variance(a,b,scale);
	}
	else
	{
		double scale = 2.0*MRSValue;
		localVariance=Variance(a,b,scale)*exp(-scale*c/K_YEAR_LEN);
		return localVariance;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: FwdZcLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
///          <=> FwdZcLocalCovariance(a,b,T1,T2,T1,T2)
///           but specialised to save a function call
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

    /// May cause regression in calibration when MRS is set to 0
//    if(fabs(MRSValue)<=ARM_ModelParamsHW::MrsMinValue)
//        MRSValue=(MRSValue>0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    if(fabs(MRSValue)<=K_NEW_DOUBLE_TOL)
        MRSValue=(MRSValue>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double localVariance=Variance(a,b,2.0*MRSValue);

    double MRS=MRSValue/K_YEAR_LEN;
    double scale=(exp(-MRS*T1)-exp(-MRS*T2))/MRSValue;

    return scale*scale*localVariance;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: FwdZcLocalcovariance
///	Returns: value of the covariance
///	Action : Covariance in [a,b] of Zc(.,T1)/Zc(.,U1)
///          and Zc(.,T2)/Zc(.,U2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,ARM_GP_Vector& vars) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

    /// May cause regression in calibration when MRS is set to 0
//    if(fabs(MRSValue)<=ARM_ModelParamsHW::MrsMinValue)
//        MRSValue=(MRSValue>0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    if(fabs(MRSValue)<=K_NEW_DOUBLE_TOL)
        MRSValue=(MRSValue>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double localVariance=Variance(a,b,2.0*MRSValue);

    double MRS=MRSValue/K_YEAR_LEN;
    double invMRSValue=1/MRSValue;
    double betaT1U1=(exp(-MRS*T1)-exp(-MRS*U1))*invMRSValue;
    double betaT2U2=(exp(-MRS*T2)-exp(-MRS*U2))*invMRSValue;

    if(vars.size()>=2)
    {
        /// Capitalise variance and betas computation
        vars[0]=betaT1U1*betaT1U1*localVariance;
        vars[1]=betaT2U2*betaT2U2*localVariance;
    }

    return betaT1U1*betaT2U2*localVariance;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: ZcVarianceSpread
///	Returns: value of the varaince spread
///	Action : Integ{ s=0->t1 VolZc(s,T1)^2 }
///         -Integ{ s=0->t2 VolZc(s,T2)^2 }
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd:: ZcVarianceSpread(double t1,double t2,double T1,double T2) const
{
    double MRSValue = ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue=(MRSValue>0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    double t = (t1<t2 ? t1 : t2);

    double totalVar1=Variance(0.0,t,MRSValue);
    double totalVar2=Variance(0.0,t,2.0*MRSValue);

    double MRS      = MRSValue/K_YEAR_LEN;
    double expT1    = exp(-MRS*T1);
    double expT2    = exp(-MRS*T2);
    double expT12   = expT1*expT1;
    double expT22   = expT2*expT2;

    // Better than (totalVar2*(expT12-expT22) - 2.0*totalVar1*(expT1-expT2))/(MRSValue*MRSValue);
    double betaT2T1 = BetatT(T2,T1);
    double value = betaT2T1*expT2*(2.0*totalVar1-(expT1+expT2)*totalVar2)/MRSValue;

    if(fabs(t1-t2) > K_NEW_DOUBLE_TOL)
    {
        /// To avoid ill behaviour of exp() limit MRS to 10-4
        if(fabs(MRSValue)<=0.0001)
            MRSValue=(MRSValue>0 ? 0.0001 : -0.0001);
        if(t1 < t2)
        {
            expT2 = exp(-MRSValue*T2/K_YEAR_LEN);
            value -= (Variance(t1,t2,0.0) - 2.0*expT2*Variance(t1,t2,MRSValue) + expT2*expT2*Variance(t1,t2,2.0*MRSValue))/(MRSValue*MRSValue);
        }
        else
        {
            expT1 = exp(-MRSValue*T1/K_YEAR_LEN);
            value += (Variance(t2,t1,0.0) - 2.0*expT1*Variance(t2,t1,MRSValue) + expT1*expT1*Variance(t2,t1,2.0*MRSValue))/(MRSValue*MRSValue);
        }
    }

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: StateZcCovariance
///	Returns: value of the covariance
///	Action : Covariance in [0,t] between the
///          state variable and Zc(.,T)
///
///          !!!! Not tested !!!!
///
///          Could be tested if Monte-Carlo method
///          and Zc numeraire are used
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FStd::StateZcCovariance(double t,double T) const
{
    double MRSValue= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

    /// May cause regression in calibration when MRS is set to 0
//    if(fabs(MRSValue)<=ARM_ModelParamsHW::MrsMinValue)
//        MRSValue=(MRSValue>0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    if(fabs(MRSValue)<=K_NEW_DOUBLE_TOL)
        MRSValue=(MRSValue>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double totalVar1=Variance(0.0,t,MRSValue);
    double totalVar2=Variance(0.0,t,2.0*MRSValue);
    double invMRS=1.0/MRSValue;
    double MRS=MRSValue/K_YEAR_LEN;

    return exp(-MRS*t)*invMRS*(exp(-MRS*T)*totalVar2-totalVar1);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW1FExt
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1FExt::ARM_ModelParamsHW1FExt( const ARM_ModelParamsHW1FExt& rhs )
: ARM_ModelParamsHW1F(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW1FExt::ARM_ModelParamsHW1FExt( const ARM_ModelParamVector& params, int volatilityType )
: ARM_ModelParamsHW1F(params,volatilityType )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW1FExt::~ARM_ModelParamsHW1FExt()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW1FExt& ARM_ModelParamsHW1FExt::operator=(const ARM_ModelParamsHW1FExt& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW1F::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1FExt
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsHW1FExt::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsHW1FExt\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();
    return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1FExt
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHW1FExt::Clone() const
{
	return new ARM_ModelParamsHW1FExt(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: Lambda
///	Returns: value of Lambda(t)
///	Action : Compute Lambda(t)=Integ{a->b>=a,MRS(s)ds}
///          for a stepwise right constant MRS curve
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::Lambda(double a,double b) const
{
    if(b<=a+K_NEW_DOUBLE_TOL)
        return 0.0;

    ARM_GP_Vector MRSTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetAbscisses()); // Vi
    ARM_GP_Vector MRSValues(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()); // MRS(Vi)

    /// Find Vn(a)-1 < a <= Vn(a) and Vn(b)-1 < b <= Vn(b)
    int i,nbV=MRSTimes.size();
    for(i=0;i<nbV;++i)
        if(a <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbV;++i)
        if(b <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double V,lastV,MRS;

    double value=0.0;
    lastV=a;

    /// Between a and Vn(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        V=(i<nb ? MRSTimes[i] : b);
        MRS=MRSValues[i<nbV ? i : nbV-1];
        if(fabs(MRS)<=K_NEW_DOUBLE_TOL)
            MRS=(MRS>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
        value += MRS*(V-lastV);
        lastV=V;
    }

    return value/K_YEAR_LEN;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: ScaledBetatT
///	Returns: value of Beta(a,b,scale)
///	Action : Compute Beta(a,b,scale)=Integ{a->b>=a,exp[scale*(Lambda(v)-Lambda(a))]dv}
///          with Lambda(v)=Integ{0->v,MRS(s)ds}
///          exp{scale*[Lambda(b)-Lambda(a)]}
///          is also return by reference
///          MRS is a stepwise right constant curve
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt:: ScaledBetatT(double a,double b,double scale,double& expLambda) const
{
    ARM_GP_Vector MRSTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetAbscisses()); // Vi
    ARM_GP_Vector MRSValues(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()); // MRS(Vi)

    /// Find Vn(a)-1 < a <= Vn(a) and Vn(b)-1< b <= Vn(b)
    int i,nbV=MRSTimes.size();
    for(i=0;i<nbV;++i)
        if(a <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbV;++i)
        if(b <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double V,lastV,MRS;
    bool notTiny;

    const double tinyLevel=0.1*K_NEW_DOUBLE_TOL;
    const double invYear=1.0/K_YEAR_LEN;

    double expx;

    expLambda=1.0;
    lastV=a;

    double value=0.0;

    /// Between a and Vn(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        V=(i<nb ? MRSTimes[i] : b);
        MRS=(scale==-1 ? -MRSValues[i<nbV ? i : nbV-1] : scale*MRSValues[i<nbV ? i : nbV-1]);
        notTiny=(fabs(MRS)>tinyLevel);
        expx=exp(MRS*(V-lastV)*invYear);
        value += expLambda*(notTiny ? (expx-1.0)/MRS : (V-lastV)*invYear);
        expLambda *= expx;
        lastV=V;
    }

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: Variance
///	Returns: value of the integration
///	Action : Integrate sigma(t)^2*exp(2.0*[Lambda(t)-Lambda(a)])
///          between [a,b] for stepwise right
///          constant sigma curve
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::Variance(double a,double b) const
{
    ARM_GP_Vector sigmaTimes( ( (ARM_CurveModelParam&) GetModelParam( GetVolatilityType() )).GetCurve()->GetAbscisses());  // Ui
    ARM_GP_Vector sigmaValues(( (ARM_CurveModelParam&) GetModelParam( GetVolatilityType() )).GetCurve()->GetOrdinates()); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
    int i,nbU=sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(a <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbU;++i)
        if(b <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double U,lastU,sig;
    double expx,expLambda;

    double value=0.0;

    expLambda=1.0;
    lastU=a;

    /// Between a and Un(b)-1 then b
    for(i=na;i<=nb;i++)
    {
        U=(i<nb ? sigmaTimes[i] : b);
        sig=sigmaValues[i<nbU ? i : nbU-1];
		if(fabs(sig) < ARM_ModelParamsHW::VOL_LIMIT)
        sig = (sig > 0.0 ? ARM_ModelParamsHW::VOL_LIMIT : -ARM_ModelParamsHW::VOL_LIMIT);
        value += sig*sig*expLambda*ScaledBetatT(lastU,U,2.0,expx);
        expLambda *= expx;
        lastU=U;
    }

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: VecBetatT
///	Returns: a vector of Beta(Vj,T)
///	Action : Let t<=T, compute a vector of Beta(Vj,T)
///          for Vj=Vn(t) if Vn(t)<=T else t and down to V0=0
///          using a stepwise right constant MRS curve
///
///          !!!! Not tested !!!!
///
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW1FExt:: VecBetatT(double t,double T) const
{
    ARM_GP_Vector MRSTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetAbscisses()); // Vi
    ARM_GP_Vector MRSValues(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()); // MRS(Vi)

    /// Find Vn(t)-1 < t <= Vn(t)
    int i,nbV=MRSTimes.size();
    for(i=0;i<nbV;++i)
        if(t <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nt=i;

    double V,nextV,MRS;

    const double tinyLevel=0.1*K_NEW_DOUBLE_TOL;
    const double invYear=1.0/K_YEAR_LEN;

    double expx,expLambda;

    ARM_GP_Vector* values=new ARM_GP_Vector(nt+1);

    if(nt<nbV-1)
    {
        nextV=MRSTimes[nt];
        nextV=(nextV<=T+K_NEW_DOUBLE_TOL ? nextV : t);
    }
    else
        nextV=t;

    (*values)[nt]=ScaledBetatT(nextV,T,-1.0,expLambda);

    for(i=nt-1;i>=0;--i)
    {
        V=MRSTimes[i];
        MRS=MRSValues[i+1<nbV ? i+1 : nbV-1];
        if(fabs(MRS)>tinyLevel)
        {
            expx=exp(-MRS*(nextV-V)*invYear);
            (*values)[i] = (*values)[i+1]*expx + (1.0-expx)/MRS;
        }
        else
            (*values)[i] = (*values)[i+1] + (nextV-V)*invYear;

        nextV=V;
    }

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: IntegLambdaBeta
///	Returns: value of the integration
///	Action : Integrate exp(Lambda(s))*beta(s,T)
///          between [a,b] for stepwise right
///          constant MRS curve
///          Needs a vector of beta(Vj,T)
///
///          !!!! Not tested !!!!
///
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::IntegLambdaBeta(double a,double b,double T,ARM_GP_Vector& betatTs) const
{
    if(b<=a+K_NEW_DOUBLE_TOL)
        return 0.0;

    ARM_GP_Vector MRSTimes(( (ARM_CurveModelParam&)  GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetAbscisses()); // Vi
    ARM_GP_Vector MRSValues(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()); // MRS(Vi)

    /// Find Vn(a)-1 < a <= Vn(a) and Vn(b)-1 < b <= Vn(b)
    int i,nbV=MRSTimes.size();
    for(i=0;i<nbV;++i)
        if(a <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbV;++i)
        if(b <= MRSTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    const double tinyLevel=0.1*K_NEW_DOUBLE_TOL;
    const double invYear=1.0/K_YEAR_LEN;

    double V,lastV,MRS;
    double x,expLambda,invMRS;

    double value=0.0;

    lastV=a;
    expLambda=exp(Lambda(0.0,a));

    /// Between 0 and Vn(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        V=(i<nb ? MRSTimes[i] : b);
        MRS=MRSValues[i<nbV ? i : nbV-1];
        x=(V-lastV)*invYear;
        if(fabs(MRS)>tinyLevel)
        {
            x=exp(MRS*x);
            invMRS=1.0/MRS;
            value += expLambda*invMRS*( (x-1.0)*invMRS + (betatTs[i]-invMRS)*0.5*(x-1/x) );
            expLambda *= x;
        }
        else
            value += expLambda*x*( 0.5*x + betatTs[i] );

        lastV=V;
    }

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: BetatT
///	Returns: value of beta(t,T)
///	Action : Compute Beta(t,T)=Integ{t->T>=t,exp[-(Lambda(v)-Lambda(t))]dv}
///          with Lambda(v)=Integ{0->v,MRS(s)ds}
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::BetatT(double a,double b) const
{
    double unusedExpLambda;

    return ScaledBetatT(a,b,-1.0,unusedExpLambda);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: StateLocalDrift
///	Returns: value of the drift
///	Action : Relative drift of the state variable
///          from a to b>=a
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::StateLocalDrift(double a,double b) const
{
    return exp(-Lambda(a,b));
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
///          that is : Integ(a->b) { sigma(s)^2.exp(-2.Lambda(s,c)).ds }
///          with Lambda(s,c) = Integ(s->c) { MRS(u).du }
///          Note that Variance(a,b) computes the integral a -> b of
///          sigma(s)^2.exp(2.Lambda(a,s)).ds
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::StateLocalVariance(double a,double b,double c) const
{
    return exp(-2.0*Lambda(a,c))*Variance(a,b);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: FwdZcLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
///          <=> FwdZcLocalCovariance(a,b,T1,T2,T1,T2)
///           but specialised to save a function call
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    double explambda = exp(-2.0*Lambda(a,T1));
    double var = Variance(a,b);
    var *= explambda; 
    double betaT1T2=BetatT(T1,T2);

    double result = var*betaT1T2*betaT1T2;   

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: FwdZcLocalcovariance
///	Returns: value of the covariance
///	Action : Covariance in [a,b] of Zc(.,T1)/Zc(.,U1)
///          and Zc(.,T2)/Zc(.,U2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,ARM_GP_Vector& vars) const
{
    double betaT1U1=BetatT(T1,U1);
    double betaT2U2=BetatT(T2,U2);

    if(vars.size()>=2)
    {
        /// Capitalise variance and betas computation
        vars[0]=exp(-2*Lambda(a,T1))*Variance(a,b)*betaT1U1*betaT1U1;
        vars[1]=exp(-2*Lambda(a,T2))*Variance(a,b)*betaT2U2*betaT2U2;
    }

    return exp(-(Lambda(a,T1)+Lambda(a,T2)))*Variance(a,b)*betaT1U1*betaT2U2;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: ZcVarianceSpread
///	Returns: value of the variance spread
///	Action : Integ{ s=0->t VolZc(s,T1)^2 - VolZc(s,T2)^2 }
///
///          !!!! Not tested !!!!
///
///          Could be tested if Cash numeraire is used
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt:: ZcVarianceSpread(double t1,double t2,double T1,double T2) const
{

    if(fabs(t1-t2) > K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " t1!=t2 : not implemented case !" );

    double t = (t1<t2 ? t1 : t2);

    ARM_GP_Vector sigmaTimes(( (ARM_CurveModelParam&) GetModelParam(GetVolatilityType()  )).GetCurve()->GetAbscisses());  // Ui
    ARM_GP_Vector sigmaValues(( (ARM_CurveModelParam&) GetModelParam(GetVolatilityType() )).GetCurve()->GetOrdinates());  // Sigma(Ui)

    /// Find Un(t)-1 < t <= Un(t)
    int i,nbU=sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(t <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nt=i;

    double U,lastU,sig;

    /// Compute a set of Beta(Vj,T) for V0=0 to Vn(t) if Vn(t)<=T else t
    ARM_GP_Vector* betatT1s=VecBetatT(t,T1);
    ARM_GP_Vector* betatT2s=VecBetatT(t,T2);

    lastU=0.0;

    double value=0.0;

    /// Between 0 and Un(t)-1 then t
    for(i=0;i<=nt;i++)
    {
        U=(i<nt ? sigmaTimes[i] : t);
        sig=sigmaValues[i<nbU ? i : nbU-1];
		if(fabs(sig) < ARM_ModelParamsHW::VOL_LIMIT)
        sig = (sig > 0.0 ? ARM_ModelParamsHW::VOL_LIMIT : -ARM_ModelParamsHW::VOL_LIMIT);
        value += sig*sig*( IntegLambdaBeta(lastU,U,T1,*betatT1s)+IntegLambdaBeta(lastU,U,T1,*betatT2s) );
        lastU=U;
    }

    delete betatT1s;
    delete betatT2s;

    return -exp(-Lambda(0.0,T1))*BetatT(T1,T2)*value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FExt
///	Routine: StateZcCovariance
///	Returns: value of the covariance
///	Action : Covariance in [0,t] between the state variable and Zc(.,T)
///
///          !!!! Not tested !!!!
///
///          Could be tested if Monte-Carlo method
///          and Zc numeraire are used
////////////////////////////////////////////////////
double ARM_ModelParamsHW1FExt::StateZcCovariance(double t,double T) const
{
    ARM_GP_Vector sigmaTimes(( (ARM_CurveModelParam&) GetModelParam(GetVolatilityType() )).GetCurve()->GetAbscisses());  // Ui
    ARM_GP_Vector sigmaValues(( (ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetOrdinates());  // Sigma(Ui)

    /// Find Un(t)-1 < t <= Un(t)
    int i,nbU=sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(t <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nt=i;

    double U,lastU,sig;

    /// Compute set of Beta(Vj,T) for V0=0 to Vn(t) if Vn(t)<=T else t
    ARM_GP_Vector* betatTs=VecBetatT(t,T);

    lastU=0.0;

    double value=0.0;

    /// Between 0 and Un(t)-1 then t
    for(i=0;i<=nt;i++)
    {
        U=(i<nt ? sigmaTimes[i] : t);
        sig=sigmaValues[i<nbU ? i : nbU-1];
		if(fabs(sig) < ARM_ModelParamsHW::VOL_LIMIT)
        sig = (sig > 0.0 ? ARM_ModelParamsHW::VOL_LIMIT : -ARM_ModelParamsHW::VOL_LIMIT);
        value += sig*sig*IntegLambdaBeta(lastU,U,T,*betatTs);
        lastU=U;
    }

    delete betatTs;

    return -exp(-Lambda(0.0,t))*value;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_CalibParamsHW1FExt
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


////////////////////////////////////////////////////
///	Class  : ARM_CalibParamsHW1FExt
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_CalibParamsHW1FExt::ARM_CalibParamsHW1FExt( const ARM_CalibParamsHW1FExt& rhs )
: ARM_ModelParamsHW1FExt(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CalibParamsHW1FExt
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_CalibParamsHW1FExt::ARM_CalibParamsHW1FExt( const ARM_ModelParamVector& params, int volatilityType )
: ARM_ModelParamsHW1FExt(params,volatilityType)
{}


////////////////////////////////////////////////////
///	Class  : ARM_CalibParamsHW1FExt
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_CalibParamsHW1FExt::~ARM_CalibParamsHW1FExt()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1FStd
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_CalibParamsHW1FExt& ARM_CalibParamsHW1FExt::operator=(const ARM_CalibParamsHW1FExt& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW1FExt::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_CalibParamsHW1FExt
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_CalibParamsHW1FExt::Clone() const
{
	return new ARM_CalibParamsHW1FExt(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_CalibParamsHW1FExt
///	Routines: FwdZcLocalVariance
///	Returns :
///	Action  : Calculate a local volatility of fwd ZC
////////////////////////////////////////////////////
double ARM_CalibParamsHW1FExt::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

    double integb,value,
        intega = 0.0;
    if(DoesModelParamExist( GetVolatilityType() ))
    {
        const ARM_CurveModelParam& foundVol = (const ARM_CurveModelParam&) GetModelParam( GetVolatilityType() );
        integb = foundVol.GetCurve()->Interpolate(b);
        if(fabs(a) > K_NEW_DOUBLE_TOL)
            intega = foundVol.GetCurve()->Interpolate(a);

        value = integb -  intega;
    }
    else
    {
        CC_Ostringstream msg;
			msg << " No modelParam for volatility type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }

    double integlambda = 0.0;

    if (DoesModelParamExist( ARM_ModelParamType::MeanReversion ))
    {
        const ARM_CurveModelParam& foundMRS = (const ARM_CurveModelParam&) GetModelParam( ARM_ModelParamType::MeanReversion );
        integlambda = ExpLambdaIntegral(T1,T2,&foundMRS.GetCurve()->GetAbscisses(), &foundMRS.GetCurve()->GetOrdinates());
    }
    else
    {
        CC_Ostringstream msg;
			msg << " No modelParam for mean reversion type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }

    double result = value*integlambda*integlambda;

    return result;    
}

////////////////////////////////////////////////////
///	Class  : ARM_CalibParamsHW1FExt
///	Routine: ScaledBetatT
///	Returns: value of Beta(a,b,scale)
///	Action : Compute Beta(a,b,scale)=Integ{a->b>=a,exp[scale*(Lambda(v)-Lambda(a))]dv}
///          with Lambda(v)=Integ{0->v,MRS(s)ds}
///          exp{scale*[Lambda(b)-Lambda(a)]}
///          is also return by reference
///          MRS is a stepwise right constant curve
////////////////////////////////////////////////////
double ARM_CalibParamsHW1FExt::ExpLambdaIntegral(double a,double b,ARM_GP_Vector* times, 
                                                               ARM_GP_Vector* values) const
{
    /// Find Vn(a)-1 < a <= Vn(a) and Vn(b)-1< b <= Vn(b)
    int i,nbV=times->size();
    for(i=0;i<nbV;++i)
        if(a <= (*times)[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbV;++i)
        if(b <= (*times)[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double V,lastV,MRS;
    bool notTiny;

    const double tinyLevel=0.1*K_NEW_DOUBLE_TOL;
    const double invYear=1.0/K_YEAR_LEN;

    double expx;

    double expLambda=1.0;
    lastV=a;

    double value=0.0;

    /// Between a and Vn(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        V=(i<nb ? (*times)[i] : b);
        MRS= -(*values)[i<nbV ? i : nbV-1];
        notTiny=(fabs(MRS)>tinyLevel);
        expx=exp(MRS*(V-lastV)*invYear);
        value += expLambda*(notTiny ? (expx-1.0)/MRS : (V-lastV)*invYear);
        expLambda *= expx;
        lastV=V;
    }

    return value;
}



////////////////////////////////////////////////////
///	Class  : ARM_CalibParamsHW1FExt
///	Routine: SquaredIntegral
///	Returns: value of the integration
///	Action : Integrate sigma(t)^2*exp(2.0*[Lambda(t)-Lambda(a)])
///          between [a,b] for stepwise right
///          constant sigma curve
////////////////////////////////////////////////////

double ARM_CalibParamsHW1FExt::SquaredIntegral(double a,double b,
                                               ARM_GP_Vector* times, 
                                               ARM_GP_Vector* values) const
{
    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
    int i,nbU=times->size();
    for(i=0;i<nbU;++i)
        if(a <= (*times)[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbU;++i)
        if(b <= (*times)[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double U,lastU,sig;
    double expx,expLambda;

    double value=0.0;

    expLambda=1.0;
    lastU=a;

    /// Between a and Un(b)-1 then b
    for(i=na;i<=nb;i++)
    {
        U=(i<nb ? (*times)[i] : b);
        sig=(*values)[i<nbU ? i : nbU-1];
		if(fabs(sig) < ARM_ModelParamsHW::VOL_LIMIT)
        sig = (sig > 0.0 ? ARM_ModelParamsHW::VOL_LIMIT : -ARM_ModelParamsHW::VOL_LIMIT);
        value += sig*sig*expLambda*GetScaledBetatT(lastU,U,2.0,expx);
        expLambda *= expx;
        lastU=U;
    }

    return value;
}

////////////////////////////////////////////////////
///	Class   : ARM_CalibParamsHW1FExt
///	Routines: PreProcessing
///	Returns :
///	Action  : Calculate the integrals lieing to 
/// volatility and mean reversion
////////////////////////////////////////////////////

void ARM_CalibParamsHW1FExt::PreProcessing(ARM_ModelFitter& modelFitter, int factorNb)
{
    ARM_ModelParamVector CalibParamVector = modelFitter.GetCalibParams();

    ARM_GP_Vector* timesVolatility = NULL;
    ARM_ModelParamVector::iterator foundVol = CC_NS( std, find_if) ( CalibParamVector.begin(), CalibParamVector.end(), 
                        FindModelParamWEnumUnaryVersion( (ARM_ModelParamType::ParamNb) GetVolatilityType() ) );
	ARM_CurveModelParam* foundVolCalibParam = NULL;

    if(foundVol!= CalibParamVector.end())
    {
		foundVolCalibParam = dynamic_cast<ARM_CurveModelParam*>(*foundVol);
		if( foundVolCalibParam  )
			timesVolatility = &foundVolCalibParam->GetCurve()->GetAbscisses();
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " : should be a curve Calib Param!" );
	}
    else
    {
        CC_Ostringstream msg;
		msg << " No modelParam for volatility type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }

    ARM_GP_Vector* timesMeanReversion = NULL;
    ARM_ModelParamVector::iterator foundMRS = CC_NS( std, find_if) ( CalibParamVector.begin(), CalibParamVector.end(), 
        FindModelParamWEnumUnaryVersion( ARM_ModelParamType::MeanReversion ) );
	ARM_CurveModelParam* foundMRSCalibParam = NULL;

    if(foundMRS!= CalibParamVector.end())
	{
		foundMRSCalibParam = dynamic_cast<ARM_CurveModelParam*>(*foundMRS);
		if( foundMRSCalibParam )
			timesMeanReversion = &foundMRSCalibParam->GetCurve()->GetAbscisses();
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " : should be a curve Calib Param!" );
	}
    else
    {
        CC_Ostringstream msg;
			msg << " No modelParam for mean reversion type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }
    
    size_t size = timesVolatility->size();
    ARM_GP_Vector varvector(size);
    ARM_GP_Vector lowerbound(size);
    ARM_GP_Vector upperbound(size);
    size_t i;

    for( i= 0; i<size; ++i)  
    {
        double explambda = exp(-2.0*GetLambda(0,(*timesMeanReversion)[i]));
        double var = SquaredIntegral(0,(*timesVolatility)[i],timesVolatility,&foundVolCalibParam->GetCurve()->GetOrdinates());
        varvector[i] = explambda*var;
        lowerbound[i] = SquaredIntegral(0,(*timesVolatility)[i],timesVolatility,&*(foundVolCalibParam->GetLowerBound()));
        upperbound[i] = SquaredIntegral(0,(*timesVolatility)[i],timesVolatility,&*(foundVolCalibParam->GetUpperBound()));

    }     
    foundVolCalibParam->GetCurve()->SetOrdinates(varvector);
    foundVolCalibParam->SetLowerBound(&lowerbound);
    foundVolCalibParam->SetUpperBound(&upperbound);


    ARM_ModelParamVector::iterator foundVol1 = CC_NS( std, find_if) (begin(),end(), 
		FindModelParamWEnumUnaryVersion( (ARM_ModelParamType::ParamNb) GetVolatilityType() ) );

    if(DoesModelParamExist( GetVolatilityType() ))
    {
        const ARM_CurveModelParam& foundVol1 = (const ARM_CurveModelParam&) GetModelParam( GetVolatilityType() );
        foundVol1.GetCurve()->SetOrdinates(varvector);
    }
    else
    {
        CC_Ostringstream msg;
			msg << " No modelParam for volatility type,  please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_CalibParamsHW1FExt
///	Routines: VarianceToSigma
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

void ARM_CalibParamsHW1FExt::VarianceToSigma(ARM_ModelParamVector& CalibParamVector,ARM_PricingModel* model)
{  

    ARM_ModelParamVector::iterator foundVol = CC_NS( std, find_if) ( CalibParamVector.begin(), CalibParamVector.end(), 
		FindModelParamWEnumUnaryVersion( (ARM_ModelParamType::ParamNb) GetVolatilityType() ) );
	ARM_CurveModelParam* foundVolCalibParam = NULL;

	/// test that it exists and is the correct type!
    if(foundVol!= CalibParamVector.end() )
	{
		if( !(foundVolCalibParam = dynamic_cast<ARM_CurveModelParam*>(*foundVol)) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " : should be a curve Calib Param!" );
	}
    else
    {
        CC_Ostringstream msg;
			msg << " No calibParam for volatility type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }

    VarDiffFunc targetfunc(this,model, (ARM_ModelParamType::ParamNb) GetVolatilityType());    
    const double target   = 0.0;
    const double ftol      = 1e-6; 
	const double xtol      = 1e-6; 
    T_SmoothNewtonRaphsonSolver<VarDiffFunc> solver(targetfunc,target,ftol,xtol,ARM_MAX_ITER);  
       
    size_t size2 = (*foundVol)->size();
    targetfunc.SetParamType( (ARM_ModelParamType::ParamNb) GetVolatilityType());
    int i;  
    for( i= 0; i < size2; ++i)
    {
        double lag = foundVolCalibParam->GetCurve()->GetAbscisses()[i];
        int spotDays = model->GetZeroCurve()->GetCurrencyUnit()->GetSpotDays();
        char Calendar;
        model->GetZeroCurve()->GetCurrencyUnit()->CalcFloatPayCal(&Calendar);
        ARM_Date date = ARM_Date(lag + model->GetAsOfDate().GetJulian()).NextBusinessDay(spotDays,&Calendar);
        double lagmean = date - model->GetAsOfDate();
        double lambda = GetLambda(0,lagmean);
        lambda = exp(-2.0*lambda);

        double target =foundVolCalibParam->GetCurve()->GetOrdinates()[i];
        target /= lambda;
        targetfunc.SetParams(i,target, lag);
        double guess_initial = foundVolCalibParam->GetCurve()->GetOrdinates()[i];
        solver.setInitialGuess(guess_initial);
        double root = fabs(solver.Solve());
        foundVolCalibParam->SetValueAtPoint(i, root); 
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_CalibParamsHW1FExt
///	Routines: PostProcessing
///	Returns :
///	Action  :  Call Newton raphson to update sigma and lambda curves
////////////////////////////////////////////////////

void ARM_CalibParamsHW1FExt::PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb)
{
    ARM_ModelParamVector CalibParamVector = modelFitter.GetCalibParams();
    ARM_ModelParamVector::iterator foundVol = CC_NS( std, find_if) ( CalibParamVector.begin(), CalibParamVector.end(), 
		FindModelParamWEnumUnaryVersion( (ARM_ModelParamType::ParamNb) GetVolatilityType() ) );
	
	ARM_CurveModelParam* foundVolCalibParam = NULL;
	
	/// test that it exists and is the correct type!
    if(foundVol!= CalibParamVector.end() )
	{
		if( !(foundVolCalibParam = dynamic_cast<ARM_CurveModelParam*>(*foundVol)) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " : should be a curve Calib Param!" );
	}
    else
    {
        CC_Ostringstream msg;
		msg << " No calibParam for volatility type please, advice!!" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, msg.str() );
    }
    
    size_t size2 = (*foundVol)->size();
    ARM_GP_Vector value(size2);
    int i;  
    double lastU	= 0.0;
    double sum		= 0.0;
    double expLambda= 1.0;
    double expx;

    for(i=0; i<size2; i++)
    {
        double U = foundVolCalibParam->GetCurve()->GetAbscisses()[i];
        int spotDays = model->GetZeroCurve()->GetCurrencyUnit()->GetSpotDays();
        char Calendar;
        model->GetZeroCurve()->GetCurrencyUnit()->CalcFloatPayCal(&Calendar);
        ARM_Date date = ARM_Date(U + model->GetAsOfDate().GetJulian()).NextBusinessDay(spotDays,&Calendar);
        double lagmean = date - model->GetAsOfDate();
        double lambda = GetLambda(0,lagmean);
        lambda = exp(-2.0*lambda);
        double target =foundVolCalibParam->GetCurve()->GetOrdinates()[i]/lambda;
		
        double value = expLambda*GetScaledBetatT(lastU,U,2.0,expx);
        expLambda *= expx;
		
        double dist = target - sum;
        double sigma;
        if(dist > K_FRM_TOL)
			sigma = sqrt((target - sum )/value);
        else
            sigma = 1.e-5;
		
        foundVolCalibParam->SetValueAtPoint(i,sigma);
        sum += sigma*sigma*value;
        lastU=U;
    }
}



////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FStateCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between states variables of 2 H&W1F models.
///           So integrate sigma1(t)*sigma2(t)*exp(-MRS1*(T-t))*exp(-MRS2*(T-t))
///           between [a,b] for a stepwise right constant sigma curves
///           (only valid for cst MRSs !)
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FStateCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(lhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");

    ARM_GP_Vector sigmaValues1( ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes1(  ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)->
    ARM_GP_Vector sigmaValues2( ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes2(  ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

	/// GetMeanReversion 1 and 2
    double MRSValue1 = ((ARM_CurveModelParam&) lhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    double MRSValue2 = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

	return HW1FStateCovarianceWithVec(sigmaTimes1,sigmaValues1,sigmaTimes2,sigmaValues2,MRSValue1,MRSValue2,a,b,T);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FStateCovarianceWithVec
///	Returns : double
///	Action  : Compute the covariance between states variables of 2 H&W1F models.
///           So integrate sigma1(t)*sigma2(t)*exp(-MRS1*(T-t))*exp(-MRS2*(T-t))
///           between [a,b] for a stepwise right constant sigma curves
///           (only valid for cst MRSs !)
///			  It uses ARM_GP_Vector
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FStateCovarianceWithVec(
		const ARM_GP_Vector& sigmaTimes1,
		const ARM_GP_Vector& sigmaValues1,
		const ARM_GP_Vector& sigmaTimes2,
		const ARM_GP_Vector& sigmaValues2,
		double MRSValue1,
		double MRSValue2,
		double a,
		double b,
		double T)
{
	/// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(sigmaTimes1,a);
	int nb1 = lower_boundPosWithPrecision(sigmaTimes1,b);
	int na2 = lower_boundPosWithPrecision(sigmaTimes2,a);
	int nb2 = lower_boundPosWithPrecision(sigmaTimes2,b);
    int nbTimes1 = sigmaTimes1.size();
    int nbTimes2 = sigmaTimes2.size();

	if(fabs(MRSValue1)<=K_NEW_DOUBLE_TOL)
        MRSValue1=(MRSValue1>=0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

	if(fabs(MRSValue2)<=K_NEW_DOUBLE_TOL)
        MRSValue2=(MRSValue2>=0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

	/// mean reversion
	double timeScale=(MRSValue1+MRSValue2)/K_YEAR_LEN;
    bool notTiny=(fabs(timeScale)>0.1*K_NEW_DOUBLE_TOL);

    /// Between a and Un(b)-1 then b
	size_t i1,i2;
	double value,lastU,UMiddle,U,lastExpScale,newExpScale,sig1,sig2;

	value = 0.0;
	lastU = a;
	lastExpScale=notTiny ? exp(timeScale*a) : 1.0;

    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
        sig1=sigmaValues1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;
        
		U=(i1<nb1 ? sigmaTimes1[i1] : b);

	    sig2=sigmaValues2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		
		UMiddle = lastU;
		while( i2<nb2 && sigmaTimes2[i2]<U )
		{
			UMiddle = sigmaTimes2[i2];
			if( UMiddle > lastU )
			{
				newExpScale=notTiny ? exp(timeScale*UMiddle) : 1.0;
				value += sig1*sig2*(notTiny ? (newExpScale-lastExpScale): UMiddle-lastU);
				lastExpScale=newExpScale;
			}
			lastU = UMiddle;
			
			++i2;

            sig2=sigmaValues2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		}
		
		/// last bit
		if( UMiddle < U )
		{
	        newExpScale=exp(timeScale*U);
			value += sig1*sig2*(notTiny ? newExpScale-lastExpScale: U-UMiddle);
		}

        lastExpScale=newExpScale;
        lastU=U;
    }

    if(notTiny)
        return value * exp(-timeScale*T) / (MRSValue1+MRSValue2);
    else
        return value / K_YEAR_LEN;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FZcCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between Zc1(.,T) & Zc2(.,T) : these Zc are
///           related to two different markets modelised both by H&W 1F models.
///           So : integrate sigma1(t)*sigma2(t)*beta1(t,T)*beta2(t,T)
///           between [a,b] for a stepwise right constant sigma curves
///           with beta(t,T)=(1-exp(-MRS*(T-t))/MRS (only valid for cst MRSs !)
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FZcCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T1, double T2 )
{
	if (T2 == -1)
		T2 = T1;

    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(lhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");

    ARM_GP_Vector sigmaValues1( ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes1(  ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)->
    ARM_GP_Vector sigmaValues2( ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes2(  ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(sigmaTimes1,a);
	int nb1 = lower_boundPosWithPrecision(sigmaTimes1,b);
	int na2 = lower_boundPosWithPrecision(sigmaTimes2,a);
	int nb2 = lower_boundPosWithPrecision(sigmaTimes2,b);
    int nbTimes1 = sigmaTimes1.size();
    int nbTimes2 = sigmaTimes2.size();

	/// GetMeanReversion 1 and 2
    double MRSValue1 = ((ARM_CurveModelParam&) lhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue1)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue1=(MRSValue1>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    double MRSValue2 = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue2)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue2=(MRSValue2>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    double MRSValueSum = MRSValue1 + MRSValue2;
    if(fabs(MRSValueSum)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValueSum=(MRSValueSum>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    double invMRS1 = 1.0/MRSValue1;
    double invMRS2 = 1.0/MRSValue2;
    double invMRSSum = 1.0/MRSValueSum;

    /// No optimisation w.r.t. very low value for MRSs because too many cases !!
	double timeScale1=MRSValue1/K_YEAR_LEN;
	double timeScale2=MRSValue2/K_YEAR_LEN;
	double timeScale12=timeScale1+timeScale2;

    /// Between a and Un(b)-1 then b
	size_t i1,i2;
	double value,lastU,UMiddle,U,sig1,sig2;
    double newExpScale1,newExpScale2,newExpScale12;

	value = 0.0;
	lastU = a;
	double lastExpScale1    = exp(timeScale1*(a-T1));
	double lastExpScale2    = exp(timeScale2*(a-T2));
    double lastExpScale12   = lastExpScale1*lastExpScale2;

    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
        sig1=sigmaValues1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;
        
		U=(i1<nb1 ? sigmaTimes1[i1] : b);

	    sig2=sigmaValues2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		
		UMiddle = lastU;
		while( i2<nb2 && sigmaTimes2[i2]<U )
		{
			UMiddle = sigmaTimes2[i2];
			if( UMiddle > lastU )
			{
				newExpScale1    = exp(timeScale1*(UMiddle-T1));
				newExpScale2    = exp(timeScale2*(UMiddle-T2));
                newExpScale12   = newExpScale1*newExpScale2;

                value += sig1*sig2*( (UMiddle-lastU)/K_YEAR_LEN - (newExpScale1-lastExpScale1)*invMRS1
                                     - (newExpScale2-lastExpScale2)*invMRS2
                                     + (newExpScale12-lastExpScale12)*invMRSSum );

				lastExpScale1   = newExpScale1;
				lastExpScale2   = newExpScale2;
				lastExpScale12  = newExpScale12;
			}
			lastU = UMiddle;
			
			++i2;

            sig2=sigmaValues2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		}
		
		/// last bit
		if( UMiddle < U )
		{
	        newExpScale1    = exp(timeScale1*(U-T1));
	        newExpScale2    = exp(timeScale2*(U-T2));
	        newExpScale12   = newExpScale1*newExpScale2;

            value += sig1*sig2*( (U-UMiddle)/K_YEAR_LEN - (newExpScale1-lastExpScale1)*invMRS1
                                 - (newExpScale2-lastExpScale2)*invMRS2
                                 + (newExpScale12-lastExpScale12)*invMRSSum );
		}

        lastExpScale1   = newExpScale1;
        lastExpScale2   = newExpScale2;
        lastExpScale12  = newExpScale12;
        lastU=U;
    }

    return value * invMRS1 * invMRS2;
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FStateZcCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between the state variable and a Zc(.,T) of 2 H&W1F models.
///           So : integrate -sigma1(t)*exp(-MRS1*(T1-t))*sigma2(t)*beta2(t,T2)
///           between [a,b] for a stepwise right constant sigma curves
///           with beta(t,T2)=(1-exp(-MRS*(T2-t))/MRS (only valid for cst MRSs !)
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FStateZcCovariance( const ARM_ModelParamsHW1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T1, double T2 )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(lhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");

    ARM_GP_Vector sigmaValues1( ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes1(  ((ARM_CurveModelParam&) lhs->GetModelParam(lhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)->
    ARM_GP_Vector sigmaValues2( ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes2(  ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(sigmaTimes1,a);
	int nb1 = lower_boundPosWithPrecision(sigmaTimes1,b);
	int na2 = lower_boundPosWithPrecision(sigmaTimes2,a);
	int nb2 = lower_boundPosWithPrecision(sigmaTimes2,b);
    int nbTimes1 = sigmaTimes1.size();
    int nbTimes2 = sigmaTimes2.size();

	/// GetMeanReversion 1 and 2
    double MRSValue1 = ((ARM_CurveModelParam&) lhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue1)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue1=(MRSValue1>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    double MRSValue2 = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue2)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue2=(MRSValue2>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    double MRSValueSum = MRSValue1 + MRSValue2;
    if(fabs(MRSValueSum)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValueSum=(MRSValueSum>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    double invMRS1 = 1.0/MRSValue1;
    double invMRS2 = 1.0/MRSValue2;
    double invMRSSum = 1.0/(MRSValue1+MRSValue2);

    /// No optimisation w.r.t. very low value for MRSs because too many cases !!
	double timeScale1=MRSValue1/K_YEAR_LEN;
	double timeScale2=MRSValue2/K_YEAR_LEN;

    /// Between a and Un(b)-1 then b
	size_t i1,i2;
	double value,lastU,UMiddle,U,sig1,sig2;
    double newExpScale1,newExpScale2;

	value = 0.0;
	lastU = a;
	double lastExpScale1    = exp(timeScale1*(a-T1));
	double lastExpScale2    = exp(timeScale2*(a-T2));

    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
        sig1=sigmaValues1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;
        
		U=(i1<nb1 ? sigmaTimes1[i1] : b);

	    sig2=sigmaValues2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		
		UMiddle = lastU;
		while( i2<nb2 && sigmaTimes2[i2]<U )
		{
			UMiddle = sigmaTimes2[i2];
			if( UMiddle > lastU )
			{
				newExpScale1    = exp(timeScale1*(UMiddle-T1));
				newExpScale2    = exp(timeScale2*(UMiddle-T2));

                value -= sig1*sig2*( (newExpScale1-lastExpScale1)*invMRS1
                                     - (newExpScale1*newExpScale2-lastExpScale1*lastExpScale2)*invMRSSum );

				lastExpScale1   = newExpScale1;
				lastExpScale2   = newExpScale2;
			}
			lastU = UMiddle;
			
			++i2;

            sig2=sigmaValues2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		}
		
		/// last bit
		if( UMiddle < U )
		{
	        newExpScale1    = exp(timeScale1*(U-T1));
	        newExpScale2    = exp(timeScale2*(U-T2));

            value -= sig1*sig2*( (newExpScale1-lastExpScale1)*invMRS1
                                 - (newExpScale1*newExpScale2-lastExpScale1*lastExpScale2)*invMRSSum );
		}

        lastExpScale1   = newExpScale1;
        lastExpScale2   = newExpScale2;
        lastU=U;
    }

    return value * invMRS2;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FEqFxZcCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between a pure lognormal spot equity or forex and
///           a Zc(.,T) modelised by a H&W 1F model (volZc(t,T)=-sigma2(t)*beta2(t,T))
///           So : integrate sigma1(t)*(-sigma2(t)*beta2(t,T))
///           between [a,b] for a stepwise right constant sigma curve
///           with beta(t,T)=(1-exp(-MRS*(T-t))/MRS (only valid for cst MRS !)
///           Specialized function of HW1FStateZcCovariance() for a 1st H&W model with
///           a MRS set to 0
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FEqFxZcCovariance( const ARM_CurveModelParam& eqFxVol, const ARM_ModelParamsHW1F* rhs, double a, double b, double T )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");

    ARM_GP_Vector sigmaValues1( eqFxVol.GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes1(  eqFxVol.GetCurve()->GetAbscisses() ); // Sigma(Ui)->
    ARM_GP_Vector sigmaValues2( ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes2(  ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(sigmaTimes1,a);
	int nb1 = lower_boundPosWithPrecision(sigmaTimes1,b);
	int na2 = lower_boundPosWithPrecision(sigmaTimes2,a);
	int nb2 = lower_boundPosWithPrecision(sigmaTimes2,b);
    int nbTimes1 = sigmaTimes1.size();
    int nbTimes2 = sigmaTimes2.size();

	/// GetMeanReversion
    double MRSValue = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue=(MRSValue>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

	/// mean reversion
    double invMRS = 1.0/MRSValue;
	double timeScale=MRSValue/K_YEAR_LEN;

    /// Between a and Un(b)-1 then b
	size_t i1,i2;
	double value,lastU,UMiddle,U,lastExpScale,newExpScale,sig1,sig2;

	value = 0.0;
	lastU = a;
	lastExpScale=exp(timeScale*(a-T));

    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
        sig1=sigmaValues1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;
        
		U=(i1<nb1 ? sigmaTimes1[i1] : b);

	    sig2=sigmaValues2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		
		UMiddle = lastU;
		while( i2<nb2 && sigmaTimes2[i2]<U )
		{
			UMiddle = sigmaTimes2[i2];
			if( UMiddle > lastU )
			{
	            newExpScale = exp(timeScale*(UMiddle-T));

                value -= sig1*sig2*( (UMiddle-lastU)/K_YEAR_LEN - (newExpScale-lastExpScale)*invMRS );

				lastExpScale = newExpScale;
			}
			lastU = UMiddle;
			
			++i2;

            sig2=sigmaValues2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		}
		
		/// last bit
		if( UMiddle < U )
		{
	        newExpScale = exp(timeScale*(U-T));

            value -= sig1*sig2*( (U-UMiddle)/K_YEAR_LEN - (newExpScale-lastExpScale)*invMRS );
		}

        lastExpScale = newExpScale;
        lastU=U;
    }

    return value * invMRS;
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: HW1FEqFxStateCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between a pure lognormal spot equity or forex and
///           the state variable of a H&W 1F model
///           So : integrate sigma1(t)*sigma2(t)*exp(-MRS*(T-t))
///           between [a,b] for a stepwise right constant sigma curve
///           Specialized function of HW1FStateZcCovariance() for a 1st H&W model with
///           a MRS set to 0
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::HW1FEqFxStateCovariance( const ARM_CurveModelParam& eqFxVol, const ARM_ModelParamsHW1F* rhs, double a, double b, double T )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");

    ARM_GP_Vector sigmaValues1( eqFxVol.GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes1(  eqFxVol.GetCurve()->GetAbscisses() ); // Sigma(Ui)->
    ARM_GP_Vector sigmaValues2( ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates() ); // Ui
    ARM_GP_Vector sigmaTimes2(  ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)

	/// GetMeanReversion
    double MRSValue = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

	return HW1FStateCovarianceWithVec(sigmaTimes1,sigmaValues1,sigmaTimes2,sigmaValues2,0,MRSValue,a,b,T);
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsHW1F::ModelParamsTimeSteps() const
{
    ARM_GP_Vector sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve()->GetAbscisses() ); // Sigma(Ui)
	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW1F
///	Routines: IsLn
///	Returns : bool
///	Action  : Check if the model is lognormal
////////////////////////////////////////////////////

bool ARM_ModelParamsHW1F::IsLn() const
{
	ARM_Curve* fxSkewParam;

	if (DoesModelParamExist( ARM_ModelParamType::Beta ))
	{
		fxSkewParam = GetModelParam(ARM_ModelParamType::Beta).ToCurveModelParam().GetCurve();
	}
	else if (DoesModelParamExist( ARM_ModelParamType::QParameter ))
	{
		fxSkewParam = GetModelParam(ARM_ModelParamType::QParameter).ToCurveModelParam().GetCurve();
	}
	else
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, "There is no skew parameter in this model param");
	}

	return *fxSkewParam == 1.0;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: BetatT (STATIC VERSION)
///	Returns: value fo beta(t,T)
///	Action : beta(t,T)=(1-exp(-MRS*(T-t))/MRS
///                   = Integ{t->T,exp(-MRS*(u-t))du}
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::BetatT(const  ARM_ModelParam& mrsParam, double t,double T)
{
	double MRSValue = static_cast<const ARM_CurveModelParam&>(mrsParam).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue)>K_NEW_DOUBLE_TOL)
        return (1.0-exp(-MRSValue*(T-t)/K_YEAR_LEN))/MRSValue;
    else
        return (T-t)/K_YEAR_LEN;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW1F
///	Routine: DerivBetatT (STATIC VERSION)
///	Returns: value fo d(beta(t,T))/dt
///	Action : = exp(-MRS*(T-t))
////////////////////////////////////////////////////
double ARM_ModelParamsHW1F::DerivBetatT(const  ARM_ModelParam& mrsParam, double t,double T)
{
	double MRSValue = static_cast<const ARM_CurveModelParam&>(mrsParam).GetCurve()->GetOrdinates()[0];
    return (exp(-MRSValue*(T-t)/K_YEAR_LEN));
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

