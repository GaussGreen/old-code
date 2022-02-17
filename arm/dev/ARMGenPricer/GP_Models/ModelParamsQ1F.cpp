/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQ1F.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


/// this headers has to come first
/// as it contains pragma for stupid warning
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/ModelParamsQ1F.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/ostringstream.h"
#include "gpbase/vectormanip.h"
#include "gpbase/comparisonfunctor.h"

/// gpinfra
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"

/// ARM Kernel
#include "glob/expt.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsQ1F::ARM_ModelParamsQ1F( const ARM_ModelParamVector& params )
: ARM_ModelParamsHW1FStd(params, ARM_ModelParamType::QVol)
{
	if( params.size() != 3 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected 3 model parameters: Q Vol, Q param, Mean Reversion!" );

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_ModelParamsQ1F::ValidateModelParams() const
{	
    /// checks that the q model is of size 1 since the current model is with cst q!
	if( !DoesModelParamExist(ARM_ModelParamType::QParameter) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: has to have a q parameter!");

	/// check the vol type
	if( !DoesModelParamExist(ARM_ModelParamType::QVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: requires a Q vol!");

	/// check that there is a mean reversion
	if( !DoesModelParamExist(ARM_ModelParamType::MeanReversion) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: requires mean reversion!");

	if( GetModelParam(ARM_ModelParamType::MeanReversion).size() > 1)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: constant mean reversion only supported!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsQ1F::ARM_ModelParamsQ1F( const ARM_ModelParamsQ1F& rhs )
: ARM_ModelParamsHW1FStd(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsQ1F::~ARM_ModelParamsQ1F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsQ1F& ARM_ModelParamsQ1F::operator=(const ARM_ModelParamsQ1F& rhs)
{
	if(this != &rhs)
		ARM_ModelParamsHW1FStd::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsQ1F::Clone() const
{
	return new ARM_ModelParamsQ1F(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ModelParamsQ1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "ARM_ModelParamsQ1F\n";
    os << "----------------------\n";
    os << ARM_ModelParams::toString();
    return os.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: GetModelCurves
///	Returns : void
///	Action  : Get sigma & q curves merging schedules if necessary
////////////////////////////////////////////////////
size_t ARM_ModelParamsQ1F::GetModelCurves(ARM_GP_Vector& times, ARM_GP_Vector& sigmas, ARM_GP_Vector& qs, ARM_GP_Vector& dqs) const
{
    /// Get model curves
    const ARM_Curve& sigmaCurve = * ((ARM_CurveModelParam&) GetModelParam(GetVolatilityType())).GetCurve();
    const ARM_GP_Vector& sigmaValues    = sigmaCurve.GetOrdinates();
    const ARM_GP_Vector& sigmaTimes     = sigmaCurve.GetAbscisses();

    const ARM_Curve& qCurve = * ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter)).GetCurve();
    const ARM_GP_Vector& qValues        = qCurve.GetOrdinates();
    const ARM_GP_Vector& qTimes         = qCurve.GetAbscisses();

    const ARM_CurveInterpolator* interpolator    = qCurve.GetInterpolator();
    if( !dynamic_cast<const ARM_LinInterpCstExtrapolDble*>(interpolator) &&
        !dynamic_cast<const ARM_StepUpRightOpenCstExtrapolDble*>(interpolator) )
		ARM_THROW(ERR_INVALID_ARGUMENT, ": Q model requires linear or constant right interpolated Q parameter");

    /// Merge sigma and q schedules
    ARM_GP_VectorPtr mergedTimes(MergeSortedVectorNoDuplicates(sigmaTimes,qTimes));

    size_t i,nbTimes=mergedTimes->size();
    times.resize(nbTimes);
    sigmas.resize(nbTimes);
    qs.resize(nbTimes);
    dqs.resize(nbTimes+1);

    double lastt=0.0,t,dt;
    if((*mergedTimes)[0] < lastt - K_NEW_DOUBLE_TOL)
		ARM_THROW(ERR_INVALID_ARGUMENT, ": first time lag in curve schedule must be >= 0.0");

    double q1,q2;

    for(i=0;i<nbTimes;++i)
    {
        t = (*mergedTimes)[i];
        times[i]    = t;
        sigmas[i]   = sigmaCurve.Interpolate(t);
        qs[i]       = qCurve.Interpolate(t);

        /// Compute qs and qds such that q(t) = qs + dqs.(t-tp+1) for t in ]tp,tp+1[
        /// Use t+dt & t+2dt to avoid interpolation rules (necessarily right cst or linear)
        dt      = (t-lastt)/3.0;
        q1      = qCurve.Interpolate(lastt+dt);
        q2      = qCurve.Interpolate(lastt+2*dt);
        dqs[i]  = (q2-q1)*K_YEAR_LEN/dt;
        lastt   = t;
    }
    dt      = t;
    q1      = qCurve.Interpolate(t+dt);
    q2      = qCurve.Interpolate(t+2*dt);
    dqs[i]  = (q2-q1)*K_YEAR_LEN/dt;

    return nbTimes;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: Q1FStateCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between states variables of 2 Q1F models.
///           So integrate q1(t)*sigma1(t)*q2(t)*sigma2(t)*exp(-MRS1*(T-t))*exp(-MRS2*(T-t))
///           between [a,b] for a stepwise right constant sigma curve and a cst or linear q curve
///           If isrhsQ is false then 2nd process is assumed to be a simple H&W1F
///           (only valid for cst MRSs !)
////////////////////////////////////////////////////
double ARM_ModelParamsQ1F::Q1FStateCovariance( const ARM_ModelParamsQ1F* lhs, const ARM_ModelParamsQ1F* rhs, double a, double b, double T, bool isrhsQ )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

    /// 1st model curves
    ARM_GP_Vector sigmas1,qs1,dqs1,times1;
    size_t nbTimes1 = lhs->GetModelCurves(times1,sigmas1,qs1,dqs1);

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(times1,a);
	int nb1 = lower_boundPosWithPrecision(times1,b);

	/// MRS1
    double MRSValue1 = ((ARM_CurveModelParam&) lhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue1)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue1=(MRSValue1>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    /// 2nd model curves
    ARM_GP_Vector sigmas2,qs2,dqs2,times2;
    int na2,nb2;
    size_t nbTimes2;
    double MRSValue2;
    if(rhs != lhs)
    {
        nbTimes2 = rhs->GetModelCurves(times2,sigmas2,qs2,dqs2);

        /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	    na2 = lower_boundPosWithPrecision(times2,a);
	    nb2 = lower_boundPosWithPrecision(times2,b);

	    /// MRS2
        MRSValue2 = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
        if(fabs(MRSValue2)<=ARM_ModelParamsHW::MrsMinValue)
            MRSValue2=(MRSValue2>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    }
    else
    {
        sigmas2     = sigmas1;
        qs2         = qs1;
        dqs2        = dqs1;
        times2      = times1;
        na2         = na1;
        nb2         = nb1;
        nbTimes2    = nbTimes1;
        MRSValue2   = MRSValue1;
    }

    double MRSValueSum = MRSValue1 + MRSValue2;
    if(fabs(MRSValueSum)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValueSum=(MRSValueSum>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    double MRSValueSum2 = MRSValueSum*MRSValueSum;

	size_t i1,i2,ii1,ii2;
	double value,lastU,UMiddle,U,lastExpScale,newExpScale,yfLastU,yfU;
    double sig1,sig2,q1,q2,dq1,dq2,qU1,qU2;
    double A,B,C;

	value = 0.0;
	lastU = a;
    yfLastU = a/K_YEAR_LEN;
	lastExpScale=exp(MRSValueSum*yfLastU);

    /// Between a and Un(b)-1 then b
    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
		U=(i1<nb1 ? times1[i1] : b);

        sig1    = sigmas1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;

        ii1     = i1<nbTimes1 ? i1 :nbTimes1-1;
        qU1     = times1[ii1]/K_YEAR_LEN;
        q1      = qs1[ii1];
        dq1     = dqs1[i1<nbTimes1 ? i1 :nbTimes1];
        

	    sig2    = sigmas2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;

        ii2     = i2<nbTimes2 ? i2 :nbTimes2-1;
        qU2     = times2[ii2]/K_YEAR_LEN;
        if(isrhsQ)
        {
            q2  = qs2[ii2];
            dq2 = dqs2[i2<nbTimes2 ? i2 :nbTimes2];
        }
        else
        {
            q2  = 1.0;
            dq2 = 0.0;
        }

        A   = dq1*dq2;
        B   = dq1*(-dq2*qU2+q2) + dq2*(-dq1*qU1+q1);
        C   = (-dq1*qU1+q1)*(-dq2*qU2+q2);
		
		UMiddle = lastU;
		while( i2<nb2 && times2[i2]<U )
		{
			UMiddle = times2[i2];
			if( UMiddle > lastU )
			{
                yfU = UMiddle/K_YEAR_LEN;
				newExpScale=exp(MRSValueSum*yfU);
				value += sig1 * sig2 * (
                    newExpScale * ( yfU*(A*yfU+B)+C - (2*A*yfU+B)/MRSValueSum + 2*A/MRSValueSum2 ) - 
                    lastExpScale* ( yfLastU*(A*yfLastU+B)+C - (2*A*yfLastU+B)/MRSValueSum + 2*A/MRSValueSum2 ) );
				lastExpScale=newExpScale;
			}
			lastU   = UMiddle;
            yfLastU = yfU;
			
			++i2;

            sig2    = sigmas2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;

            ii2     = i2<nbTimes2 ? i2 :nbTimes2-1;
            qU2     = times2[ii2]/K_YEAR_LEN;
            if(isrhsQ)
            {
                q2  = qs2[ii2];
                dq2 = dqs2[i2<nbTimes2 ? i2 :nbTimes2];
            }
            else
            {
                q2  = 1.0;
                dq2 = 0.0;
            }

            A   = dq1*dq2;
            B   = dq1*(-dq2*qU2+q2) + dq2*(-dq1*qU1+q1);
            C   = (-dq1*qU1+q1)*(-dq2*qU2+q2);
		}
		
		/// last bit
		if( UMiddle < U )
		{
            yfU = U/K_YEAR_LEN;
	        newExpScale=exp(MRSValueSum*yfU);
			value += sig1 * sig2 * (
                newExpScale * ( yfU*(A*yfU+B)+C - (2*A*yfU+B)/MRSValueSum + 2*A/MRSValueSum2 ) - 
                lastExpScale* ( yfLastU*(A*yfLastU+B)+C - (2*A*yfLastU+B)/MRSValueSum + 2*A/MRSValueSum2 ) );
		}

        lastExpScale=newExpScale;
        lastU   = U;
        yfLastU = yfU;
    }

    return value * exp(-MRSValueSum*T/K_YEAR_LEN) / MRSValueSum;
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F
///	Routines: Q1FStateZcCovariance (static routine)
///	Returns : double
///	Action  : Compute the covariance between the state variable and a Zc(.,T) of a Q1F & H&W1F models.
///           So : integrate q1(t)*sigma1(t)*exp(-MRS1*(T1-t))*(-sigma2(t)*beta2(t,T2))
///           between [a,b] for a stepwise right constant sigma curves
///           with beta(t,T2)=(1-exp(-MRS*(T2-t))/MRS (only valid for cst MRSs !)
////////////////////////////////////////////////////
double ARM_ModelParamsQ1F::Q1FStateZcCovariance( const ARM_ModelParamsQ1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T1, double T2 )
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	/// validation
	if( !dynamic_cast<const ARM_ModelParamsHW1FStd*>(rhs) )
		ARM_THROW(ERR_INVALID_ARGUMENT, "Currently only supported is ARM_ModelParamsHW1FStd");


    /// 1st model curves
    ARM_GP_Vector sigmas1,qs1,dqs1,times1;
    size_t nbTimes1 = lhs->GetModelCurves(times1,sigmas1,qs1,dqs1);

    /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	int na1 = lower_boundPosWithPrecision(times1,a);
	int nb1 = lower_boundPosWithPrecision(times1,b);

	/// MRS1
    double MRSValue1 = ((ARM_CurveModelParam&) lhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    if(fabs(MRSValue1)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValue1=(MRSValue1>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    /// 2nd model curves
    ARM_GP_Vector sigmas2,times2;
    int na2,nb2;
    size_t nbTimes2;
    double MRSValue2;
    if(rhs != lhs)
    {
        sigmas2 = ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetOrdinates();
        times2  = ((ARM_CurveModelParam&) rhs->GetModelParam(rhs->GetVolatilityType())).GetCurve()->GetAbscisses();
        nbTimes2 = times2.size();

        /// Find Un(a)-1 < a <= Un(a) and Un(b)-1< b <= Un(b)
	    na2 = lower_boundPosWithPrecision(times2,a);
	    nb2 = lower_boundPosWithPrecision(times2,b);

	    /// MRS2
        MRSValue2 = ((ARM_CurveModelParam&) rhs->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
        if(fabs(MRSValue2)<=ARM_ModelParamsHW::MrsMinValue)
            MRSValue2=(MRSValue2>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);
    }
    else
    {
        sigmas2     = sigmas1;
        times2      = times1;
        na2         = na1;
        nb2         = nb1;
        nbTimes2    = nbTimes1;
        MRSValue2   = MRSValue1;
    }

    
    double MRSValueSum = MRSValue1 + MRSValue2;
    if(fabs(MRSValueSum)<=ARM_ModelParamsHW::MrsMinValue)
        MRSValueSum=(MRSValueSum>=0 ? ARM_ModelParamsHW::MrsMinValue : -ARM_ModelParamsHW::MrsMinValue);

    double invMRS1 = 1.0/MRSValue1;
    double invMRS2 = 1.0/MRSValue2;
    double invMRSSum = 1.0/(MRSValue1+MRSValue2);

    /// Between a and Un(b)-1 then b
	size_t i1,i2,ii1;
	double value,lastU,UMiddle,U,yfLastU,yfU,yfT1,yfT2;
    double newExpScale1,newExpScale2;

    double sig1,sig2,q1,dq1,qU1;
    double A,B;

	value = 0.0;
	lastU   = a;
    yfLastU = a/K_YEAR_LEN;
    yfT1    = T1/K_YEAR_LEN;
    yfT2    = T2/K_YEAR_LEN;
	double lastExpScale1    = exp(MRSValue1*(yfLastU-yfT1));
	double lastExpScale2    = exp(MRSValue2*(yfLastU-yfT2));

    for( i1=na1,i2=na2;i1<=nb1; ++i1 )
    {
		U=(i1<nb1 ? times1[i1] : b);

        sig1    = sigmas1[i1<nbTimes1 ? i1 : nbTimes1-1];
		if(fabs(sig1)<ARM_ModelParamsHW::VOL_LIMIT)
			sig1 = ARM_ModelParamsHW::VOL_LIMIT;

        ii1     = i1<nbTimes1 ? i1 :nbTimes1-1;
        qU1     = times1[ii1]/K_YEAR_LEN;
        q1      = qs1[ii1];
        dq1     = dqs1[i1<nbTimes1 ? i1 :nbTimes1];
        

	    sig2    = sigmas2[i2<nbTimes2 ? i2 : nbTimes2-1];
		if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
			sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		
        A   = dq1;
        B   = -dq1*qU1+q1;

		UMiddle = lastU;
		while( i2<nb2 && times2[i2]<U )
		{
			UMiddle = times2[i2];
			if( UMiddle > lastU )
			{
                yfU = UMiddle/K_YEAR_LEN;
				newExpScale1    = exp(MRSValue1*(yfU-yfT1));
				newExpScale2    = exp(MRSValue2*(yfU-yfT2));

                value -= sig1 * sig2 * (
                    invMRS1*( newExpScale1*(A*(yfU-invMRS1)+B) - lastExpScale1*(A*(yfLastU-invMRS1)+B) )
                    -invMRSSum*( newExpScale1*newExpScale2*(A*(yfU-invMRSSum)+B) - lastExpScale1*lastExpScale2*(A*(yfLastU-invMRSSum)+B) ) );

				lastExpScale1   = newExpScale1;
				lastExpScale2   = newExpScale2;
			}
			lastU   = UMiddle;
            yfLastU = yfU;
			
			++i2;

            sig2    = sigmas2[i2<nbTimes2 ? i2 :nbTimes2-1];
			if(fabs(sig2)<ARM_ModelParamsHW::VOL_LIMIT)
				sig2 = ARM_ModelParamsHW::VOL_LIMIT;
		}
		
		/// last bit
		if( UMiddle < U )
		{
            yfU = U/K_YEAR_LEN;
	        newExpScale1    = exp(MRSValue1*(yfU-yfT1));
	        newExpScale2    = exp(MRSValue2*(yfU-yfT2));

            value -= sig1 * sig2 * (
                invMRS1*( newExpScale1*(A*(yfU-invMRS1)+B) - lastExpScale1*(A*(yfLastU-invMRS1)+B) )
                -invMRSSum*( newExpScale1*newExpScale2*(A*(yfU-invMRSSum)+B) - lastExpScale1*lastExpScale2*(A*(yfLastU-invMRSSum)+B) ) );
		}

        lastExpScale1   = newExpScale1;
        lastExpScale2   = newExpScale2;
        lastU   = U;
        yfLastU = yfU;
    }

    return value * invMRS2;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

