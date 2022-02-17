/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsHW2F.cpp
 *
 *  \brief file for the model parameters of the HW2F model
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/ModelParamsHW2F.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/cloneutilityfunc.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW2F
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2F::ARM_ModelParamsHW2F( const ARM_ModelParamsHW2F& rhs )
: ARM_ModelParamsHW(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW2F::ARM_ModelParamsHW2F( const ARM_ModelParamVector& params )
: ARM_ModelParamsHW(params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2F::~ARM_ModelParamsHW2F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW2F& ARM_ModelParamsHW2F::operator=(const ARM_ModelParamsHW2F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW2F
///	Routines: ModelParamsTimeSteps
///	Returns : ARM_GP_VectorPtr
///	Action  : volatility discretization Time Steps
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_ModelParamsHW2F::ModelParamsTimeSteps() const
{
    ARM_GP_Vector sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses() ); // Sigma(Ui)
	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW2FStd
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2FStd::ARM_ModelParamsHW2FStd( const ARM_ModelParamsHW2FStd& rhs )
: ARM_ModelParamsHW2F(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW2FStd::ARM_ModelParamsHW2FStd( const ARM_ModelParamVector& params )
: ARM_ModelParamsHW2F(params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2FStd::~ARM_ModelParamsHW2FStd()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW2FStd& ARM_ModelParamsHW2FStd::operator=(const ARM_ModelParamsHW2FStd& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW2FStd
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsHW2FStd::toString() const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsHW2FStd\n";
    os << "-------------------\n";
    os << ARM_ModelParams::toString();
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW2FStd
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHW2FStd::Clone() const
{
	return new ARM_ModelParamsHW2FStd(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: Variance
///	Returns: Vector of variances
///	Action : For a list of scaling factors,
///          integrate sigma(t)^2*exp(scale*t)
///          between [a,b] for a stepwise right
///          constant sigma curve 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW2FStd::Variance(double a,double b,ARM_GP_Vector& scale) const
{
    int nbScales=scale.size();

    if(b - K_NEW_DOUBLE_TOL <= a)
        return new ARM_GP_Vector(nbScales,0.0);

    ARM_GP_Vector sigmaValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates()); // Ui
    ARM_GP_Vector sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses()); // Sigma(Ui)

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
    int i,j,nbU=sigmaTimes.size();
    for(i=0;i<nbU;++i)
        if(a <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int na=i;
    for(;i<nbU;++i)
        if(b <= sigmaTimes[i] + K_NEW_DOUBLE_TOL) break;
    int nb=i;

    double U,lastU,sig;
    ARM_GP_Vector timeScale(nbScales);
    ARM_GP_Vector newExpScale(nbScales);
    ARM_GP_Vector lastExpScale(nbScales);
    vector< bool > notTiny(nbScales);
    for(i=0;i<nbScales;++i)
    {
        timeScale[i]=scale[i]/K_YEAR_LEN;
        lastExpScale[i]=exp(timeScale[i]*a);
        notTiny[i]=(fabs(scale[i])>0.1*K_NEW_DOUBLE_TOL);
    }

    ARM_GP_Vector* variances = new ARM_GP_Vector(nbScales,0.0);

    lastU=a;

    /// Between a and Un(b)-1 then b
    for(i=na;i<=nb;++i)
    {
        U=(i<nb ? sigmaTimes[i] : b);
        sig=sigmaValues[i<nbU ? i : nbU-1];
        sig *= sig;
        for(j=0;j<nbScales;++j)
        {
            newExpScale[j]=exp(timeScale[j]*U);
            (*variances)[j] += sig*(notTiny[j] ? newExpScale[j]-lastExpScale[j] : U-lastU);
            lastExpScale[j]=newExpScale[j];
        }
        lastU=U;
    }

    for(j=0;j<nbScales;++j)
        (*variances)[j] /= (notTiny[j] ? scale[j] : K_YEAR_LEN);

    return variances;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: StateLocalDrift
///	Returns: vector of drifts
///	Action : Relative drift of state variables
///          from a to b>=a 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW2FStd::StateLocalDrift(double a,double b) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];

    ARM_GP_Vector* localDrift = new ARM_GP_Vector(2,1.0);

    double timeba=(b-a)/K_YEAR_LEN;
    if(fabs(MRS1Value)>K_NEW_DOUBLE_TOL)
        (*localDrift)[0] = exp(-MRS1Value*timeba);

    if(fabs(MRS2Value)>K_NEW_DOUBLE_TOL)
        (*localDrift)[1] = exp(-MRS2Value*timeba);

    return localDrift;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: StateLocalDrift
///	Returns: vector of drifts
///	Action : Relative drift of state variables
///          from a to b>=a 
////////////////////////////////////////////////////

ARM_GP_MatrixPtr ARM_ModelParamsHW2FStd::StateLocalDriftVec( const ARM_GP_Vector& timeSteps ) const
{
	size_t nbSteps	= timeSteps.size();
	double step		= timeSteps[0],nextStep;
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
	ARM_GP_Matrix* result = new ARM_GP_Matrix(nbSteps-1,2, 1.0 );

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
	    double dt=(nextStep-step)/K_YEAR_LEN;

		/// [i] => local variance from ti->ti+1
	    if(fabs(MRS1Value)>K_NEW_DOUBLE_TOL)
		    (*result)(i,0) = exp(-MRS1Value*dt);
		if(fabs(MRS2Value)>K_NEW_DOUBLE_TOL)
			(*result)(i,1)= exp(-MRS2Value*dt);

		step=nextStep;
	}

// FIXMEFRED: mig.vc8 (25/05/2007 15:27:18):cast
	return static_cast<ARM_GP_MatrixPtr>(result);
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: StateLocalVariance
///	Returns: A lower triangular matrix of
///          variances/covariance
///	Action : Compute the variances/covariances matrix
///          of state variables between [a,b]
////////////////////////////////////////////////////
ARM_GP_TriangularMatrix* ARM_ModelParamsHW2FStd::StateLocalVariance(double a,double b) const
{
	double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
		MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
	
	double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
	if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
		MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
	
	ARM_GP_Vector scale(3);
	scale[0] = 2.0*MRS1Value;
	scale[1] = 2.0*MRS2Value;
	scale[2] = MRS1Value + MRS2Value;
	
	// Compute weighted volatility integrals
	ARM_GP_Vector* localVariance=Variance(a,b,scale);
	
	double volRatio=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates()[0];
	double correl=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates()[0];
	
	ARM_GP_TriangularMatrix* varCovar = new ARM_GP_TriangularMatrix(2);
	
	double timeb=b/K_YEAR_LEN;
	double expMRS1 = exp(-MRS1Value*timeb);
	double expMRS2 = volRatio * exp(-MRS2Value*timeb);
	
	// Compute variances/covariance lower triangular matrix
	(*varCovar)(0,0) = (*localVariance)[0] * expMRS1 * expMRS1;
	(*varCovar)(1,1) = (*localVariance)[1] * expMRS2 * expMRS2;
	(*varCovar)(1,0) = (*localVariance)[2] * correl  * expMRS1 * expMRS2;
	(*varCovar)(0,1) = (*localVariance)[2] * correl  * expMRS1 * expMRS2;
	delete localVariance;
	
	return varCovar;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: FwdZcLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
///          <=> FwdZcLocalCovariance(a,b,T1,T2,T1,T2)
///           but specialised to save a function call
///           and redundant exp()
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FStd::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);


    ARM_GP_Vector scale(3);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* localVariance=Variance(a,b,scale);

    double volRatio =((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates()[0];
    double correl   =((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates()[0];

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;
    double betaSpread1=(exp(-MRS1*T1)-exp(-MRS1*T2))/MRS1Value;
    double betaSpread2=volRatio*(exp(-MRS2*T1)-exp(-MRS2*T2))/MRS2Value;

    double value = betaSpread1*( (*localVariance)[0]*betaSpread1 + 2.0*correl*(*localVariance)[2]*betaSpread2 ) +
        (*localVariance)[1]*betaSpread2*betaSpread2;

    delete localVariance;

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: FwdZcLocalCovariance
///	Returns: value of the covariance
///	Action : Covariance in [a,b] of Zc(.,T1)/Zc(.,U1)
///          and Zc(.,T2)/Zc(.,U2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FStd::FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,ARM_GP_Vector& vars) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    ARM_GP_Vector scale(3);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* localVariance=Variance(a,b,scale);

    double volRatio=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates()[0];
    double correl=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates()[0];

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;
    double coef1=1/MRS1Value;
    double coef2=volRatio/MRS2Value;
    double beta1T1U1=(exp(-MRS1*T1)-exp(-MRS1*U1))*coef1;
    double beta1T2U2=(exp(-MRS1*T2)-exp(-MRS1*U2))*coef1;
    double beta2T1U1=(exp(-MRS2*T1)-exp(-MRS2*U1))*coef2;
    double beta2T2U2=(exp(-MRS2*T2)-exp(-MRS2*U2))*coef2;

    if(vars.size()>=2)
    {
        /// Capitalise variance and betas computation
        vars[0]=beta1T1U1*( (*localVariance)[0]*beta1T1U1 + 2.0*correl*(*localVariance)[2]*beta2T1U1 ) +
            (*localVariance)[1]*beta2T1U1*beta2T1U1;
        vars[1]=beta1T2U2*( (*localVariance)[0]*beta1T2U2 + 2.0*correl*(*localVariance)[2]*beta2T2U2 ) +
            (*localVariance)[1]*beta2T2U2*beta2T2U2;
    }

    double value = (*localVariance)[0]*beta1T1U1*beta1T2U2 + (*localVariance)[1]*beta2T1U1*beta2T2U2 +
        correl*((*localVariance)[2])*(beta1T1U1*beta2T2U2 + beta1T2U2*beta2T1U1);

    delete localVariance;

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: ZcVarianceSpread
///	Returns: value of the varaince spread
///	Action : Variance spread in [0,t] between
///          Zc(.,T1) and Zc(.,T2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FStd::ZcVarianceSpread(double t1, double t2, double T1,double T2) const
{
    if(fabs(t1-t2) > K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " t1!=t2 : not implemented case !" );

    double t = (t1<t2 ? t1 : t2);

    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);


    ARM_GP_Vector scale(5);
    scale[0] = MRS1Value;
    scale[1] = 2.0*MRS1Value;
    scale[2] = MRS2Value;
    scale[3] = 2.0*MRS2Value;
    scale[4] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* totalVar=Variance(0.0,t,scale);

    double volRatio=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates()[0];
    double correl=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates()[0];


    double invMRS1=1.0/MRS1Value;
    double inv2MRS1=invMRS1*invMRS1;

    double invMRS2=volRatio/MRS2Value;
    double inv2MRS2=invMRS2*invMRS2;

    double invMRS12=correl*invMRS1*invMRS2;

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;

    double exp1T1=exp(-MRS1*T1);
    double exp1T2=exp(-MRS1*T2);
    double exp2T1=exp(-MRS2*T1);
    double exp2T2=exp(-MRS2*T2);

    double value1 = (*totalVar)[0]*(exp1T1-exp1T2)*(inv2MRS1+invMRS12);
    double value2 = (*totalVar)[2]*(exp2T1-exp2T2)*(inv2MRS2+invMRS12);
    double value12 = (*totalVar)[4]*(exp1T1*exp2T1-exp1T2*exp2T2)*invMRS12;
    exp1T1 *= exp1T1;
    exp1T2 *= exp1T2;
    exp2T1 *= exp2T1;
    exp2T2 *= exp2T2;

    double value = (*totalVar)[1]*(exp1T1-exp1T2)*inv2MRS1 
		+ (*totalVar)[3]*(exp2T1-exp2T2)*inv2MRS2 + 
        2.0*(value12 - value1 - value2);

    delete totalVar;
    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FStd
///	Routine: StateZcCovariance
///	Returns: value of the covariance
///	Action : Covariance in [0,t] between the
///          state variable and Zc(.,T)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FStd::StateZcCovariance(double t,double T) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);


    ARM_GP_Vector scale(4);
    scale[0] = MRS1Value;
    scale[1] = 2.0*MRS1Value;
    scale[2] = MRS2Value;
    scale[3] = 2.0*MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* totalVar=Variance(0.0,t,scale);

    double volRatio=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates()[0];
    double correl=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates()[0];

    double invMRS1=1.0/MRS1Value;
    double invMRS2=volRatio/MRS2Value;

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;

    double exp1T=exp(-MRS1*T);
    double exp2T=exp(-MRS2*T);

    double exp1t=invMRS1*exp(-MRS1*t);
    double exp2t=invMRS2*exp(-MRS2*t);

    double value=exp1t*(exp1T*(*totalVar)[1]-(*totalVar)[0]) + exp2t*(exp2T*(*totalVar)[3]-(*totalVar)[2]);

    delete totalVar;

    return value;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// ARM_ModelParamsHW2FExt
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2FExt::ARM_ModelParamsHW2FExt( const ARM_ModelParamsHW2FExt& rhs )
: ARM_ModelParamsHW2F(rhs)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsHW2FExt::ARM_ModelParamsHW2FExt( const ARM_ModelParamVector& params )
: ARM_ModelParamsHW2F(params)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsHW2FExt::~ARM_ModelParamsHW2FExt()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsHW2FExt& ARM_ModelParamsHW2FExt::operator=(const ARM_ModelParamsHW2FExt& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParamsHW::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW2FExt
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsHW2FExt::toString() const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsHW2FExt\n";
    os << "-------------------\n";
    os << ARM_ModelParams::toString();
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsHW2FExt
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsHW2FExt::Clone() const
{
	return new ARM_ModelParamsHW2FExt(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: Variance
///	Returns: Vector of variances
///	Action : For a list of scaling factors,
///          integrate sigma(t)^2*exp(scale*t)
///          between [a,b] for a stepwise right
///          constant sigma curve 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW2FExt::Variance(double a,double b,ARM_GP_Vector& scale) const
{
    int nbScales=scale.size();
	if (nbScales%3)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Pb with scale size!" );
	int size=nbScales/3;

    if(b - K_NEW_DOUBLE_TOL <= a)
        return new ARM_GP_Vector(nbScales,0.0);

    ARM_GP_Vector sigmaValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates()); // Ui
    ARM_GP_Vector sigmaTimes(  ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses()); // Sigma(Ui)
	ARM_GP_Vector volRatioValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetOrdinates());
    ARM_GP_Vector volRatioTimes(  ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::VolatilityRatio)).GetCurve()->GetAbscisses());
	ARM_GP_Vector correlValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetOrdinates());
    ARM_GP_Vector correlTimes(  ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Correlation)).GetCurve()->GetAbscisses());
	
	int na1 = lower_boundPosWithPrecision(sigmaTimes,a);
	int nb1 = lower_boundPosWithPrecision(sigmaTimes,b);
	int na2 = lower_boundPosWithPrecision(volRatioTimes,a);
	int nb2 = lower_boundPosWithPrecision(volRatioTimes,b);
	int na3 = lower_boundPosWithPrecision(correlTimes,a);
	int nb3 = lower_boundPosWithPrecision(correlTimes,b);

    int nbTimes1 = sigmaTimes.size();
    int nbTimes2 = volRatioTimes.size();
	int nbTimes3 = correlTimes.size();

	ARM_GP_Vector timeScale(nbScales);
    ARM_GP_Vector newExpScale(nbScales);
    ARM_GP_Vector lastExpScale(nbScales);
    vector< bool > notTiny(nbScales);

	int i,j,k;
    for(i=0;i<nbScales;++i)
    {
        timeScale[i]	=scale[i]/K_YEAR_LEN;
		notTiny[i]		=(fabs(scale[i])>0.1*K_NEW_DOUBLE_TOL);
        lastExpScale[i]	=notTiny[i] ? exp(timeScale[i]*a) : 1.0;
    }

    /// Between a and Un(b)-1 then b
	size_t i1,i2,i3;
	double U,lastU;
	double U1,U2,U3,sig1,sig2,sig3;

	lastU = a;

	ARM_GP_Vector* variances = new ARM_GP_Vector(nbScales,0.0);

	bool done = false;

	i1 = na1;
	i2 = na2;
	i3 = na3;
    while (!done)
    {
        sig1=sigmaValues[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;
        
		U1=(i1<nb1 ? sigmaTimes[i1] : b);

	    sig2=volRatioValues[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2*sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT / sig1;

		U2=(i2<nb2 ? volRatioTimes[i2] : b);

		sig3=correlValues[i3<nbTimes3 ? i3 : nbTimes3-1];
		U3=(i3<nb3 ? correlTimes[i3] : b);

		U = CC_Min<double>(CC_Min<double>(CC_Min<double>(U1,U2),U3),b);
		
		for (j=0;j<size;j++)
		{
			for (k=0;k<3;k++)
				newExpScale[j*3+k]	= notTiny[j*3+k] ? exp(timeScale[j*3+k]*U) : 1.0;

			(*variances)[j*3+0]	+= sig1*sig1*			( notTiny[j*3+0] ? (newExpScale[j*3+0]-lastExpScale[j*3+0]): U-lastU);
			(*variances)[j*3+1]	+= sig1*sig1*sig2*sig2*	( notTiny[j*3+1] ? (newExpScale[j*3+1]-lastExpScale[j*3+1]): U-lastU);
			(*variances)[j*3+2]	+= sig1*sig1*sig2*sig3*	( notTiny[j*3+2] ? (newExpScale[j*3+2]-lastExpScale[j*3+2]): U-lastU);

			for (k=0;k<3;k++)
				lastExpScale[j*3+k]	= newExpScale[j*3+k];
		}

		if (U1<U+K_NEW_DOUBLE_TOL) i1++;
		if (U2<U+K_NEW_DOUBLE_TOL) i2++;
		if (U3<U+K_NEW_DOUBLE_TOL) i3++;

		lastU = U;
		done = (U>b-K_NEW_DOUBLE_TOL);
	}

    for(j=0;j<nbScales;++j)
        (*variances)[j] /= (notTiny[j] ? scale[j] : K_YEAR_LEN);

    return variances;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: StateLocalDrift
///	Returns: vector of drifts
///	Action : Relative drift of state variables
///          from a to b>=a 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW2FExt::StateLocalDrift(double a,double b) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];

    ARM_GP_Vector* localDrift = new ARM_GP_Vector(2,1.0);

    double timeba=(b-a)/K_YEAR_LEN;
    if(fabs(MRS1Value)>K_NEW_DOUBLE_TOL)
        (*localDrift)[0] = exp(-MRS1Value*timeba);

    if(fabs(MRS2Value)>K_NEW_DOUBLE_TOL)
        (*localDrift)[1] = exp(-MRS2Value*timeba);

    return localDrift;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: StateLocalDrift
///	Returns: vector of drifts
///	Action : Relative drift of state variables
///          from a to b>=a 
////////////////////////////////////////////////////

ARM_GP_MatrixPtr ARM_ModelParamsHW2FExt::StateLocalDriftVec( const ARM_GP_Vector& timeSteps ) const
{
	size_t nbSteps	= timeSteps.size();
	double step		= timeSteps[0],nextStep;
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
	ARM_GP_Matrix* result = new ARM_GP_Matrix(nbSteps-1,2, 1.0 );

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
	    double dt=(nextStep-step)/K_YEAR_LEN;

		/// [i] => local variance from ti->ti+1
	    if(fabs(MRS1Value)>K_NEW_DOUBLE_TOL)
		    (*result)(i,0) = exp(-MRS1Value*dt);
		if(fabs(MRS2Value)>K_NEW_DOUBLE_TOL)
			(*result)(i,1)= exp(-MRS2Value*dt);

		step=nextStep;
	}

// FIXMEFRED: mig.vc8 (25/05/2007 15:27:27):cast
	return static_cast<ARM_GP_MatrixPtr>(result);
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: StateLocalVariance
///	Returns: A lower triangular matrix of
///          variances/covariance
///	Action : Compute the variances/covariances matrix
///          of state variables between [a,b]
////////////////////////////////////////////////////
ARM_GP_TriangularMatrix* ARM_ModelParamsHW2FExt::StateLocalVariance(double a,double b) const
{
	double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
		MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
	
	double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
	if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
		MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);
	
	ARM_GP_Vector scale(3);
	scale[0] = 2.0*MRS1Value;
	scale[1] = 2.0*MRS2Value;
	scale[2] = MRS1Value + MRS2Value;
	
	// Compute weighted volatility integrals
	ARM_GP_Vector* localVariance=Variance(a,b,scale);
	
	ARM_GP_TriangularMatrix* varCovar = new ARM_GP_TriangularMatrix(2);
	
	double timeb=b/K_YEAR_LEN;
	double expMRS1 = exp(-MRS1Value*timeb);
	double expMRS2 = exp(-MRS2Value*timeb);
	
	// Compute variances/covariance lower triangular matrix
	(*varCovar)(0,0) = (*localVariance)[0] * expMRS1 * expMRS1;
	(*varCovar)(1,1) = (*localVariance)[1] * expMRS2 * expMRS2;
	(*varCovar)(1,0) = (*localVariance)[2] * expMRS1 * expMRS2;
	(*varCovar)(0,1) = (*varCovar)(1,0);
	delete localVariance;
	
	return varCovar;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: FwdZcLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
///          <=> FwdZcLocalCovariance(a,b,T1,T2,T1,T2)
///           but specialised to save a function call
///           and redundant exp()
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FExt::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);


    ARM_GP_Vector scale(3);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* localVariance=Variance(a,b,scale);

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;
    double betaSpread1=(exp(-MRS1*T1)-exp(-MRS1*T2))/MRS1Value;
    double betaSpread2=(exp(-MRS2*T1)-exp(-MRS2*T2))/MRS2Value;

    double value = betaSpread1*( (*localVariance)[0]*betaSpread1 + 2.0*(*localVariance)[2]*betaSpread2 ) +
        (*localVariance)[1]*betaSpread2*betaSpread2;

    delete localVariance;

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: FwdZcLocalCovariance
///	Returns: value of the covariance
///	Action : Covariance in [a,b] of Zc(.,T1)/Zc(.,U1)
///          and Zc(.,T2)/Zc(.,U2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FExt::FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,ARM_GP_Vector& vars) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    ARM_GP_Vector scale(3);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* localVariance=Variance(a,b,scale);

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;
    double coef1=1/MRS1Value;
    double coef2=1/MRS2Value;
    double beta1T1U1=(exp(-MRS1*T1)-exp(-MRS1*U1))*coef1;
    double beta1T2U2=(exp(-MRS1*T2)-exp(-MRS1*U2))*coef1;
    double beta2T1U1=(exp(-MRS2*T1)-exp(-MRS2*U1))*coef2;
    double beta2T2U2=(exp(-MRS2*T2)-exp(-MRS2*U2))*coef2;

    if(vars.size()>=2)
    {
        /// Capitalise variance and betas computation
        vars[0]=beta1T1U1*( (*localVariance)[0]*beta1T1U1 + 2.0*(*localVariance)[2]*beta2T1U1 ) +
            (*localVariance)[1]*beta2T1U1*beta2T1U1;
        vars[1]=beta1T2U2*( (*localVariance)[0]*beta1T2U2 + 2.0*(*localVariance)[2]*beta2T2U2 ) +
            (*localVariance)[1]*beta2T2U2*beta2T2U2;
    }

    double value = (*localVariance)[0]*beta1T1U1*beta1T2U2 + (*localVariance)[1]*beta2T1U1*beta2T2U2 +
        ((*localVariance)[2])*(beta1T1U1*beta2T2U2 + beta1T2U2*beta2T1U1);

    delete localVariance;

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: ZcVarianceSpread
///	Returns: value of the varaince spread
///	Action : Variance spread in [0,t] between
///          Zc(.,T1) and Zc(.,T2)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FExt::ZcVarianceSpread(double t1, double t2, double T1,double T2) const
{
    if(fabs(t1-t2) > K_NEW_DOUBLE_TOL)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " t1!=t2 : not implemented case !" );

    double t = (t1<t2 ? t1 : t2);

    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    ARM_GP_Vector scale(6);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;
    scale[3] = MRS1Value;
    scale[4] = MRS2Value;
	scale[5] = MRS1Value + MRS2Value;

    // Compute weighted volatility integrals
    ARM_GP_Vector* totalVar=Variance(0.0,t,scale);

    double invMRS1=1.0/MRS1Value;
    double inv2MRS1=invMRS1*invMRS1;

    double invMRS2=1.0/MRS2Value;
    double inv2MRS2=invMRS2*invMRS2;

    double invMRS12=invMRS1*invMRS2;

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;

    double exp1T1=exp(-MRS1*T1);
    double exp1T2=exp(-MRS1*T2);
    double exp2T1=exp(-MRS2*T1);
    double exp2T2=exp(-MRS2*T2);

    double value1 = (*totalVar)[3]*(exp1T1-exp1T2)*(inv2MRS1+invMRS12);
    double value2 = (*totalVar)[4]*(exp2T1-exp2T2)*(inv2MRS2+invMRS12);
    double value12 = (*totalVar)[2]*(exp1T1*exp2T1-exp1T2*exp2T2)*invMRS12;
    exp1T1 *= exp1T1;
    exp1T2 *= exp1T2;
    exp2T1 *= exp2T1;
    exp2T2 *= exp2T2;

    double value = (*totalVar)[0]*(exp1T1-exp1T2)*inv2MRS1 
		+ (*totalVar)[1]*(exp2T1-exp2T2)*inv2MRS2 + 
        2.0*(value12 - value1 - value2);

    delete totalVar;
    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2FExt
///	Routine: StateZcCovariance
///	Returns: value of the covariance
///	Action : Covariance in [0,t] between the
///          state variable and Zc(.,T)
////////////////////////////////////////////////////
double ARM_ModelParamsHW2FExt::StateZcCovariance(double t,double T) const
{
    double MRS1Value=((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS1Value)<=K_NEW_DOUBLE_TOL)
        MRS1Value=(MRS1Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

    double MRS2Value=MRS1Value + ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRS2Value)<=K_NEW_DOUBLE_TOL)
        MRS2Value=(MRS2Value>0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);


    ARM_GP_Vector scale(6);
    scale[0] = 2.0*MRS1Value;
    scale[1] = 2.0*MRS2Value;
    scale[2] = MRS1Value + MRS2Value;
    scale[3] = MRS1Value;
    scale[4] = MRS2Value;
	scale[5] = MRS1Value + MRS2Value;


    // Compute weighted volatility integrals
    ARM_GP_Vector* totalVar=Variance(0.0,t,scale);

    double invMRS1=1.0/MRS1Value;
    double invMRS2=1.0/MRS2Value;

    double MRS1=MRS1Value/K_YEAR_LEN;
    double MRS2=MRS2Value/K_YEAR_LEN;

    double exp1T=exp(-MRS1*T);
    double exp2T=exp(-MRS2*T);

    double exp1t=invMRS1*exp(-MRS1*t);
    double exp2t=invMRS2*exp(-MRS2*t);

    double value=exp1t*(exp1T*(*totalVar)[0]-(*totalVar)[3]) + exp2t*(exp2T*(*totalVar)[1]-(*totalVar)[4]);

    delete totalVar;

    return value;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsHW2F
///	Routine: BetatT (STATIC VERSION)
///	Returns: vector of betai(t,T)
///	Action : betai(t,T)=(1-exp(-MRSi*(T-t))/MRS
///                   = Integ{t->T,exp(-MRSi*(u-t))du}
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsHW2F::BetatT(const ARM_ModelParam& mrsParam,const ARM_ModelParam& mrsSpreadParam,double t,double T)
{
    double MRS1Value= static_cast<const ARM_CurveModelParam&>(mrsParam).GetCurve()->GetOrdinates()[0];
    double MRS2Value=MRS1Value +  static_cast<const ARM_CurveModelParam&>(mrsSpreadParam).GetCurve()->GetOrdinates()[0];

    ARM_GP_Vector* beta=new ARM_GP_Vector(2);

    if(fabs(MRS1Value)>K_NEW_DOUBLE_TOL)
        (*beta)[0] = (1.0-exp(-MRS1Value*(T-t)/K_YEAR_LEN))/MRS1Value;
    else
        (*beta)[0] = (T-t)/K_YEAR_LEN;

    if(fabs(MRS2Value)>K_NEW_DOUBLE_TOL)
        (*beta)[1] = (1.0-exp(-MRS2Value*(T-t)/K_YEAR_LEN))/MRS2Value;
    else
        (*beta)[1] = (T-t)/K_YEAR_LEN;

    return beta;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

