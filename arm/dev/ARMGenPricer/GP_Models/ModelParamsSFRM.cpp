/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file ModelParamsSFRM.cpp
 *  \brief
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 *
 */



/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/ModelParamsSFRM.h"
#include "gpmodels/SFRM.h"
#include "gpmodels/VanillaSwaptionArgSFRM.h"
#include "gpmodels/VanillaDigitalArgSFRM.h"

/// gpbase std
#include "gpbase/ostringstream.h"
#include "gpbase/checkarg.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/correlmatparam.h"

/// gpinfra
#include "gpcalib/vanillaarg.h"
#include "gpcalib/VanillaCap.h"
#include "gpcalib/modelfitter.h"

/// kernel
#include <inst/irindex.h>
#include "gpbase/datestrip.h"
#include <inst/portfolio.h>
#include <inst/swaption.h>

/// standard libraries
#include <cmath>
#include <memory>


CC_USING_NS(std,auto_ptr)

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegrateSquaredExpMRV,
///	Returns: double
///	Action :  compute the sum from time1 ti time2 for
///           exp(2.0*lambda*u) du
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegrateSquaredExpMRV( double time1, double time2 ) const
{
    double meanReversion = (( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates())[0];
    if( fabs(meanReversion) < K_NEW_DOUBLE_TOL )
        return time2 - time1;
    else
        return (exp(2.0*meanReversion*time2)- exp(2.0*meanReversion*time1))/(2.0*meanReversion);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ExpFunction,
///	Returns:double 
///	Action : compute at time 
///          exp(lambda*time) 
//////////// ////////////////////////////////////////
double ARM_ModelParamsSFRM::ExpFunction( double time ) const
{
    double meanReversion = (( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates())[0];
    if( fabs(meanReversion) < K_NEW_DOUBLE_TOL )
        return 1.0;
    else
        return exp(meanReversion*time);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSFRM::ARM_ModelParamsSFRM( const ARM_ModelParamVector& params, 
	ARM_IRIndex* index, 
	size_t factorsNb)
:
	ARM_ModelParams(params), 
	itsFactorsNb(factorsNb),
	itsIRIndex(static_cast< ARM_IRIndex * >(index->Clone())),
    itsMuCoeff(vector<ARM_VectorPtr>(0)),
    itsOneMu(ARM_VectorPtr(NULL)),
    itsFwdValues(NULL),
	/// by default the correlation type is set to correlation
	itsCorrelType(ARM_ModelParamType::Correlation) 
{
	ValidateModelParams();
	UpdateCorrelation( ARM_ModelParamType::Correlation );
	UpdateCorrelation( ARM_ModelParamType::BrownianCorrelation );
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsSFRM
///	Routines: UpdateCorrelation
///	Returns :
///	Action  : update if necessary the correlation according to its type
////////////////////////////////////////////////////
	void ARM_ModelParamsSFRM::UpdateCorrelation( ARM_ModelParamType::ParamNb correlType )
{
    /// ACP for the model factors!
    if( DoesModelParamExist(correlType) )
	{
		if( itsFactorsNb > 1 )
		{
			((ARM_CorrelMatParam&)GetModelParam( correlType) ).SetFactorCount(itsFactorsNb);
			((ARM_CorrelMatParam&)GetModelParam( correlType) ).UpdateRealizedCorrel();
			itsCorrelType = correlType;
		}
		else
		{
			/// remove the correlation if existing!
			DeleteModelParam( correlType);
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsSFRM
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the SFRM model
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::ValidateModelParams() const
{	
	if(!DoesModelParamExist(ARM_ModelParamType::Volatility))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model should contain volatility!");
	
	if(!DoesModelParamExist(ARM_ModelParamType::MeanReversion))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model should contain mean reversion!");

    /// checks that the mean reversion is of size 1 since the current model is with cst mean reversion
	if(DoesModelParamExist(ARM_ModelParamType::MeanReversion) && GetModelParam( ARM_ModelParamType::MeanReversion).size() != 1)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": currently, allow only SFRM model with constant mean reversion!");

    // check it there is a shift curve or a beta curve but not both
	if(DoesModelParamExist(ARM_ModelParamType::Shift) && DoesModelParamExist(ARM_ModelParamType::Beta)) 
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model should contain a shift or a beta (but not both)!" );
    
	if(!DoesModelParamExist(ARM_ModelParamType::Shift) && !DoesModelParamExist(ARM_ModelParamType::Beta))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model should contain a shift or a beta!");

	if(itsFactorsNb>1 )
	{
		if( !DoesModelParamExist( ARM_ModelParamType::Correlation) && !DoesModelParamExist(ARM_ModelParamType::BrownianCorrelation) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model with more than 1 factor should contain a correlation or brownian correlation matrix!");
		if( DoesModelParamExist(ARM_ModelParamType::Correlation) && DoesModelParamExist(ARM_ModelParamType::BrownianCorrelation) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an SFRM Model cannot have a brownian and a standard correlation!");
	}

}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: Copy constructor
///	Returns: 
///	Action : copy the object
////////////////////////////////////////////////////
ARM_ModelParamsSFRM::ARM_ModelParamsSFRM( const ARM_ModelParamsSFRM& rhs )
:	ARM_ModelParams(rhs), 
	itsFactorsNb( rhs.itsFactorsNb ),
    itsOneMu(rhs.itsOneMu),
    itsMuCoeff(rhs.itsMuCoeff),
	itsFwdValues(NULL),
	itsCorrelType(rhs.itsCorrelType)
{
    itsFwdValues = rhs.itsFwdValues ? (ARM_GP_Vector*)rhs.itsFwdValues->Clone() : NULL;
	itsIRIndex = rhs.itsIRIndex?(ARM_IRIndex*)rhs.itsIRIndex->Clone() : NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: oeprator=
///	Returns: 
///	Action : assignment operator 
////////////////////////////////////////////////////
ARM_ModelParamsSFRM& ARM_ModelParamsSFRM::operator=( const ARM_ModelParamsSFRM& rhs )
{
	if( this != &rhs )
	{
		ARM_ModelParams::operator=(rhs);
		itsFactorsNb	= rhs.itsFactorsNb;
		itsIRIndex		= rhs.itsIRIndex?(ARM_IRIndex*)rhs.itsIRIndex->Clone() : NULL;
        itsMuCoeff      = rhs.itsMuCoeff;
        itsOneMu        = rhs.itsOneMu;
	    itsFwdValues	= rhs.itsFwdValues ? (ARM_GP_Vector*)rhs.itsFwdValues->Clone() : NULL;
		itsCorrelType	= rhs.itsCorrelType;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_ModelParamsSFRM::~ARM_ModelParamsSFRM()
{
    delete itsIRIndex;
    itsIRIndex = NULL;

    delete itsFwdValues;
    itsFwdValues = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedVol
///	Returns: double
///	Action : computes the integrated vol
/// IntegratedVol = Sqrt(Integral Sigma(u,T)^2)*du )
///		-the part on du and use IntegratedLocalVariance
///		-the part on T = MaturityTerm(T)^2
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegratedVol(double s, double t, double T) const
{
	double integratedLocalVariance = IntegratedLocalVariance(s,t);

	if( integratedLocalVariance <= K_NEW_DOUBLE_TOL )
		integratedLocalVariance = K_NEW_DOUBLE_TOL;

	return sqrt(integratedLocalVariance) * MaturityTerm(T);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedVariance
///	Returns: double
///	Action : computes the integrated vol
/// IntegratedVol = Sum from s to t of Sigma(u,T)^2)*du
///		-the part on du and use IntegratedLocalVariance
///		-the part on T = MaturityTermSquared
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegratedVariance(double s, double t, double T) const
{
	return IntegratedLocalVariance(s,t) * MaturityTermSquared(T);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
/// IntegratedCovariance = Integral from s to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
///		-the part on du and use IntegratedLocalVariance
///		-the part on T = MaturityTermSquared
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegratedCovariance(double s, double t, double T1, double T2) const
{
    double integratedLocalVariance = IntegratedLocalVariance(s,t);
    double maturityTermT12 = MaturityTerm(T1) * MaturityTerm(T2);
   
    CC_NS(std,auto_ptr)<ARM_GP_Vector> correlT1(InterpolateCorrelation(T1));
    CC_NS(std,auto_ptr)<ARM_GP_Vector> correlT2(InterpolateCorrelation(T2));

    double sum = 0.0;
    for(size_t i=0;i<itsFactorsNb;++i)
        sum += (*correlT1)[i]*(*correlT2)[i];

    return integratedLocalVariance* maturityTermT12* sum;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
///          between L(.,T1,T2) and L(.,T1,T3)
/// IntegratedCovariance = Integral from s to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
///		-the part on du and use IntegratedLocalVariance
///		-the part on T = MaturityTermSquared
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegratedCovariance(double s, double t, double T1, double T2, double T3) const
{
    double integratedLocalVariance = IntegratedLocalVariance(s,t);

    /// T2 is always greater than T3
    if(T2 < T3)
    {
        double x=T3;
        T3=T2;
        T2=x;
    }

    /// Get SFRM tenor in year fraction
    double stdTerm=itsIRIndex->GetYearTerm()*K_YEAR_LEN;
    double term1=T2-T1, term2=T3-T1;

    double coef1,coef2,u1,u2,cor;
    size_t i;
    CC_NS(std,auto_ptr)<ARM_GP_Vector> correl1;
    CC_NS(std,auto_ptr)<ARM_GP_Vector> correl2;


    if( T2 <= T1 + K_NEW_DOUBLE_TOL )
    {
        coef1 = MaturityTerm(T1);
        return integratedLocalVariance * coef1 * coef1;
    }

    bool shortLibor;
    if( (shortLibor = (T3 <= T1 + K_NEW_DOUBLE_TOL)) )
    {
        coef2 = MaturityTerm(T1);
        correl2 = CC_NS(std,auto_ptr)<ARM_GP_Vector>(InterpolateCorrelation(T1));
        term2 = K_YEAR_LEN;
    }

    double t1=0.0,t2=0.0;
    double maturityCoef=0.0;
    while(t1 < term1)
    {
        u1 = T1 + t1;
        coef1 = MaturityTerm(u1) * (t1 + stdTerm <= term1 ? stdTerm : term1-t1);
        correl1 = CC_NS(std,auto_ptr)<ARM_GP_Vector>(InterpolateCorrelation(u1));
        if(shortLibor)
        {
            cor = 0.0;
            for(i=0;i<itsFactorsNb;++i)
                cor += (*correl1)[i]*(*correl2)[i];

            maturityCoef += coef1*coef2*cor;
        }
        else
        {
            while(t2 < term2)
            {
                u2 = T2 + t2;
                coef2 = MaturityTerm(u2) * (t2 + stdTerm <= term2 ? stdTerm : term2-t2);
                correl2 = CC_NS(std,auto_ptr)<ARM_GP_Vector>(InterpolateCorrelation(u2));

                cor = 0.0;
                for(i=0;i<itsFactorsNb;++i)
                    cor += (*correl1)[i]*(*correl2)[i];

                maturityCoef += coef1*coef2*cor;

                t2 += stdTerm;
            }
        }

        t1 += stdTerm;
    }


    return integratedLocalVariance * maturityCoef / (term1*term2);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: InterpolateCorrelation
///	Returns: double
///	Action : Compute local covariance from fromTime to toTime
///          between Fwd(resetTime1) and Fwd(resetTime2)
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsSFRM::InterpolateCorrelation(double T) const
{
    if(itsFactorsNb<=1)
        return new ARM_GP_Vector(1,1.0);
    ARM_GP_Vector correl			  = ((ARM_CorrelMatParam&)GetModelParam( itsCorrelType ) ).GetMultiCurve()->Interpolate(T);
    ARM_GP_Vector* renormalizedcorrel = new ARM_GP_Vector(itsFactorsNb);
	double sum = 0.;
    size_t i;
	for(i=0;i<itsFactorsNb;++i)
		sum	  += correl[i]*correl[i];

	/// renormalize for the correlation
	if(sum>K_FRM_TOL)
        for(i=0;i<itsFactorsNb;++i)
			(*renormalizedcorrel)[i] = correl[i]/sqrt(sum);

    return renormalizedcorrel;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedCorrelation
///	Returns: double
///	Action : Compute local covariance from fromTime to toTime
///          between Fwd(resetTime1) and Fwd(resetTime2)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::IntegratedCorrelation( double s, double t, double T1, double T2) const
{
    double correl = 0.0;
    if(fabs(t-s)<K_NEW_DOUBLE_TOL)
    {
        double normevol1 = VolatilityFunction(0.0,T1);
        double normevol2 = VolatilityFunction(0.0,T2);
        if(normevol1*normevol2 < K_NEW_DOUBLE_TOL)
            return correl;

        CC_NS(std,auto_ptr)<ARM_GP_Vector> instvol1(InstantaneousVolatility(0.0,T1));
        CC_NS(std,auto_ptr)<ARM_GP_Vector> instvol2(InstantaneousVolatility(0.0,T2));
        double sum = 0.0;
        for(size_t i=0;i<itsFactorsNb;++i)
            sum += (*instvol1)[i]*(*instvol2)[i];
        correl = sum/(normevol1*normevol2);
    }
    else
    {
        double sqrtVarianceT1 = sqrt(IntegratedVariance(s,t,T1));
        double sqrtVarianceT2 = sqrt(IntegratedVariance(s,t,T2));
        if(sqrt(sqrtVarianceT1*sqrtVarianceT1) < K_NEW_DOUBLE_TOL)
            return correl;
        double covariance  = IntegratedCovariance(s,t,T1,T2);
        correl = covariance/ (sqrtVarianceT1*sqrtVarianceT2);      
    }
    //////////////////   
    return correl;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: InstantaneousCorrelation
///	Returns: double
///	Action : Compute local covariance from fromTime to toTime
///          between Fwd(resetTime1) and Fwd(resetTime2)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::InstantaneousCorrelation(double T1, double T2) const
{
    double correl = 0.0;
    ARM_GP_Vector correlT1 = ((ARM_CorrelMatParam&)GetModelParam( itsCorrelType )).GetMultiCurve()->Interpolate(T1);
    ARM_GP_Vector correlT2 = ((ARM_CorrelMatParam&)GetModelParam( itsCorrelType )).GetMultiCurve()->Interpolate(T2);

    double sum = 0.0;
    for(size_t i=0;i<itsFactorsNb;++i)
        sum += correlT1[i]*correlT2[i];
    correl = sum;
    //////////////////   
    return correl;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: SwaptionVol
///	Returns: ARM_VectorPtr
///	Action : computes the swaption vol per factor
///		using coefficients from the relation vol swap vol fra
///		and the integrated vol per factor
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::IntegratedSwaptionVolVec(double s, double t, 
                                const ARM_VanillaSwaptionArg& arg,
                                const ARM_VectorPtr Mu)
{
    ARM_VectorPtr VolResult(new ARM_GP_Vector(itsFactorsNb));
	size_t datesSize=arg.GetFloatResetTimes()->size();
	double Ti;
	ARM_VectorPtr factor;

	for(size_t i=0;i<datesSize;++i)
	{
		Ti = (*arg.GetFloatResetTimes())[i];
		factor = IntegratedVolPerFactor(s,t,Ti);
		for(size_t j=0;j<itsFactorsNb;++j)
            (*VolResult)[j] += (*factor)[j]*(*Mu)[i];        
	}
	return VolResult;
}


////////////////////////////////////////////////////
///	Class  : DumpSwaptionVolsAndTimes
///	Routine: SwaptionVol
///	Returns: ARM_VectorPtr
///	Action : Dump times and volatilites informations
/// from a swaption
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::DumpSwaptionVolsAndTimes(
                                const ARM_VanillaSwaptionArg& arg,
								ARM_VectorPtr times,
                                ARM_VectorPtr volatilities)
{
	size_t datesSize=arg.GetFloatResetTimes()->size();
	double Ti;	
	ARM_VectorPtr factor;

	times->resize(datesSize);
	volatilities->resize(datesSize);

	double expiry = arg.GetExpiry();

	for(size_t i=0;i<datesSize;++i)
	{
		Ti = (*arg.GetFloatResetTimes())[i];
		(*times)[i] = Ti;
		factor = IntegratedVolPerFactor(0,expiry,Ti);
		for(size_t j=0;j<itsFactorsNb;++j)
            (*volatilities)[i] += (*factor)[j]*(*factor)[j];        
		(*volatilities)[i] = sqrt((*volatilities)[i]);

	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedCapVolVec
///	Returns: ARM_VectorPtr
///	Action : computes the integrated CapVolVec
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::IntegratedCapVolVec(double s, double t,
	const ARM_VanillaCapDigitalArg& arg)
{
	ARM_VectorPtr result(new ARM_GP_Vector(itsFactorsNb));

	for(size_t i=0;i< arg.GetResetTimes()->size(); ++i)
		*result += *IntegratedVolPerFactor(s,t,(*arg.GetResetTimes())[i]);
	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: LocalVolatity
///	Returns: double
///	Action : computes the vol observed at time s
///		for a vanilla arg
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::LocalVolatity(double s, double t, const ARM_VanillaArg& arg, const ARM_VectorPtr Mu)
{
	if(t-s<-K_NEW_DOUBLE_TOL)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"trying to compute a local volatility with t<=s!");

	if (t > s)
	{
		ARM_VanillaSwaptionArg* swaptionArg = dynamic_cast<ARM_VanillaSwaptionArg*>(const_cast<ARM_VanillaArg*>(&arg));
		ARM_VectorPtr integratedVol;
		if(swaptionArg)
		{
			integratedVol=IntegratedSwaptionVolVec(s,t,*swaptionArg,Mu);
		}
		else
		{
			ARM_VanillaCapDigitalArg* capArg = dynamic_cast<ARM_VanillaCapDigitalArg*>(const_cast<ARM_VanillaArg*>(&arg));
			if(capArg)
				integratedVol=IntegratedCapVolVec(s,t,*capArg);
			else
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"unknown vanilla arg!");
		}

		double volSquared=0;
		for( size_t i=0; i<integratedVol->size(); ++i )
			volSquared += (*integratedVol)[i]*(*integratedVol)[i];

		double vol=sqrt(volSquared/(t-s)*K_YEAR_LEN);
		if( vol < K_NEW_DOUBLE_TOL )
			return K_NEW_DOUBLE_TOL;
		else
			return vol;
	}
	else
	{
		return 0.0;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: InstantaneousVolatility
///	Returns: double
///	Action : computes the instantaneous Volatility
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_ModelParamsSFRM::InstantaneousVolatility(double t, double T) const
{
    double normeVol = VolatilityFunction(t,T);
    ARM_GP_Vector* result = new ARM_GP_Vector(itsFactorsNb);
    if(itsFactorsNb==1)
		(*result)[0]=normeVol;
    else
    {
        CC_NS(std,auto_ptr)<ARM_GP_Vector> correl(InterpolateCorrelation(T));
		for(size_t i=0;i<itsFactorsNb;++i)
			(*result)[i]=normeVol*(*correl)[i];
    }
    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedVolPerFactor
///	Returns: ARM_VectorPtr
///	Action : computes the integrated variance per factor
/// in one factor, returns the same as IntergratedVariance
///	which is the norm L1 of IntegratedVariancePerFactor!
///	otherwise in multi-factor, allocate the IntergratedVariance
///	accross factors according to their correlation!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::IntegratedVolPerFactor(double s, double t, double T)
{
	ARM_VectorPtr result(new ARM_GP_Vector(itsFactorsNb));
	double integratedVol = IntegratedVol(s,t,T);
	if(itsFactorsNb==1)
		(*result)[0]=integratedVol;
	else
	{
        CC_NS(std,auto_ptr)<ARM_GP_Vector> correl(InterpolateCorrelation(T));
        for(size_t i=0;i<itsFactorsNb;++i)
	        (*result)[i]=integratedVol*(*correl)[i];
	}
	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedVariancePerFactor
///	Returns: ARM_VectorPtr
///	Action : computes the integrated variance per factor
/// in one factor, returns the same as IntergratedVariance
///	which is the norm L1 of IntegratedVariancePerFactor!
///	otherwise in multi-factor, allocate the IntergratedVariance
///	accross factors according to their correlation!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::IntegratedVariancePerFactor(double s, double t, double T)
{
	ARM_VectorPtr result(new ARM_GP_Vector(itsFactorsNb));
	double integratedVariance = IntegratedVariance(s,t,T);
	if(itsFactorsNb==1)
		(*result)[0]=integratedVariance;
	else
	{
        CC_NS(std,auto_ptr)<ARM_GP_Vector> correl(InterpolateCorrelation(T));
        for(size_t i=0;i<itsFactorsNb;++i)
	        (*result)[i]=integratedVariance*(*correl)[i]*(*correl)[i];
	}
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: FwdsVarCorrelCoeffs
///	Returns: ARM_VectorPtr
///	Action : Calculate in the correlation between two forwards  
///          the part depending on fwd reset dates    
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::FwdsVarCorrelCoeffs(double T_1, double T_2)
{
	ARM_VectorPtr result(new ARM_GP_Vector(itsFactorsNb));
	double coeff_1 = MaturityTerm(T_1);
	double coeff_2 = MaturityTerm(T_2);
	double coeff = coeff_1*coeff_2;

	if(itsFactorsNb==1)
		(*result)[0]= coeff;
	else
	{
        CC_NS(std,auto_ptr)<ARM_GP_Vector> correl_1(InterpolateCorrelation(T_1));
		CC_NS(std,auto_ptr)<ARM_GP_Vector> correl_2(InterpolateCorrelation(T_2));
        for(size_t i=0;i<itsFactorsNb;++i)
	        (*result)[i]  =coeff*(*correl_1)[i]*(*correl_2)[i];
	}
	return result;
}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: IntegratedTreeStatesVariancePerFactor
///	Returns: ARM_VectorPtr
///	Action : computes the integrated variance per factor for the tree states
/// in one factor, returns the same as IntergratedVariance
///	which is the norm L1 of IntegratedVariancePerFactor!
///	otherwise in multi-factor, allocate the IntergratedVariance
///	accross factors according to their correlation!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ModelParamsSFRM::IntegratedTreeStatesVariancePerFactor(double s, double t)
{
	ARM_VectorPtr result(new ARM_GP_Vector(itsFactorsNb));
	double integratedVariance = IntegratedLocalVariance(s,t);
	if(itsFactorsNb==1)
		(*result)[0]=integratedVariance;
	else
	{
        for(size_t i=0;i<itsFactorsNb;++i)
	        (*result)[i]  =integratedVariance;
	}
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ShiftValue
///	Returns: 
///	Action : Get the shift at a given date!
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::ShiftValue(double t)
{
#ifdef __GP_STRICT_VALIDATION
    if(!DoesModelParamExist(ARM_ModelParamType::Shift))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                        "Trying to compute a shift value while the shift has not been set!");
#endif
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift)).GetCurve()->Interpolate(t);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: BetaValue
///	Returns: 
///	Action : Get the beta at a given date!
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::BetaValue(double t)
{
#ifdef __GP_STRICT_VALIDATION
    if(!DoesModelParamExist(ARM_ModelParamType::Beta))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                        "Trying to compute a shift value while the shift has not been set!");
#endif
	return ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->Interpolate(t);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ComputeVolSwapvolFRA
///	Returns: void
///	Action : computes the vol swap vol fra relationship
/// \latexonly
///
///		VolForward VolSwap Relationship
///
///		The connecting relation between swap and forward rate in the SFRM framework is given by:
///		\begin{equation*}
///		\vec{\sigma}_{S_{i,j}}^{BS}(t)=\sum\limits_{p=i}^{p=j}\mu _{p}^{i,j}(t)\vec{%
///		\sigma}_{F_{p}}^{BS}(t) 
///		\end{equation*}
///		
///		where
///		
///		\begin{eqnarray*}
///		\mu _{p}^{i,j}(t) &=&\frac{S_{i,j}(t)}{S_{i,j}(t)+m_{S_{,ji}}}\frac{%
///		a_{p}(F_{p}(t)+m_{p})}{(B(t,T_{0})-B(t,T_{j+1}))(1+a_{p}F_{p}(t))} \\
///		&&\left[ B(t,T_{j})+\frac{S_{i,j}(t)}{S_{i+k(i),j}(t)}\left(
///		B(t,T_{i+k(i)-1})-B(t,T_{j+1})\right) \right]
///		\end{eqnarray*}
///		
///		where $k(i)$ is the integer part of $\frac{p-i}{f}$.
///		
///		where $f$ is the number of \ variable flows in fix flow, $m_{p}$ the shift
///		parameter of the forward rate $F_{p}(.)$ and $m_{S_{i,j}}$the shift
///		parameter of the swap rate $S_{i,j}(.)$
///		
///		In the pratice, we approximate $\mu _{p}^{i,j}(t)$ \ by $\mu_{p}^{i,j}(0)$.
///
///	\endlatexonly
///////////////////////////////////////////////////////

ARM_VectorPtr ARM_ModelParamsSFRM::ComputeVolSwapvolFRA( 
	            const ARM_VanillaSwaptionArg& arg, 
	            const ARM_PricingModel& model ,
				bool isConstantNotional,
				bool useFixFrequency)
{
	int mode = 0;
	
	size_t nbVarFlow = arg.GetFloatResetTimes()->size();
	size_t nbFixFlow = arg.GetFixPayTimes()->size();

	ARM_GP_Vector* Mu = new ARM_GP_Vector(nbVarFlow);
	double averageShift= AverageShift(*arg.GetFloatResetTimes());

	ARM_GP_Vector DFFloatStarti(nbVarFlow+1);
	int i;

	//Determination of fix leg
	double DFix=0.0;
	double SwapNotional = 1;
	ARM_GP_Vector FixFlow(nbFixFlow);
	ARM_GP_Vector FixStarti(nbFixFlow+1);
	ARM_GP_Vector DFFixStarti(nbFixFlow+1);


	FixStarti[0]	= (*arg.GetFloatStartTimes())[0];
	DFFixStarti[0]	= model.GetZeroCurve()->DiscountPrice(FixStarti[0]/K_YEAR_LEN);
	for(i=0;i<nbFixFlow;++i)
	{
		if (!isConstantNotional)
			SwapNotional	= (*arg.GetFixNotional())[i]; 
		FixStarti[i+1]  = (*arg.GetFixPayTimes())[i];
		DFFixStarti[i+1]= model.GetZeroCurve()->DiscountPrice(FixStarti[i+1]/K_YEAR_LEN);
		FixFlow[i]		= SwapNotional * DFFixStarti[i+1]*(*arg.GetFixPayPeriods())[i];
		DFix		   += FixFlow[i];
	}

	double stubFix=0;

	/// see if there is a stub
	double asOfDate  = model.GetAsOfDate().GetJulian();
	ARM_Date startDate( asOfDate+(*arg.GetFloatStartTimes())[0]);
	ARM_Date firstEndDate(asOfDate + (*arg.GetFixPayTimes())[0]);
	ARM_Date stdStartDate(firstEndDate);
	ARM_Currency* pCcy = const_cast<ARM_PricingModel&>(model).GetZeroCurve()->GetCurrencyUnit();
	int fixedFreq = arg.GetFixFrequency();
	stdStartDate.AddPeriod(-fixedFreq,pCcy->GetCcyName());

	double alphaFix=0;
	bool stub = false;
	if(fabs(stdStartDate-startDate)>FRMVOL_LAG_THRESHOLD && nbFixFlow > 1)
	{
		if (!isConstantNotional)
			SwapNotional	= (*arg.GetFixNotional())[0]; 
		int fixDayCount = pCcy->GetFixedDayCount();
		double delta0	= CountYears(fixDayCount,stdStartDate,startDate);
		alphaFix		= delta0*DFFixStarti[1]*SwapNotional;
		stubFix			= delta0*(DFFixStarti[1]-DFFixStarti[0])*SwapNotional;
	}

	/////////////////////////mode0 Forward At 0//////////////////////////
	if(mode ==0)
	{
		double DFloat = 0.0;
		if (isConstantNotional)
			DFloat= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN)
				-model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);
		else
		{
			int NbVarFlow = arg.GetFloatEndTimes()->size();
			ARM_GP_Vector VarFlow(NbVarFlow);
			ARM_GP_Vector VarStarti(NbVarFlow+1);
			ARM_GP_Vector VarEndi(NbVarFlow+1);
			
			for(i=0;i<NbVarFlow;++i)
			{
				SwapNotional = (*arg.GetFloatNotional())[i];
				VarEndi[i+1]  = (*arg.GetFloatEndTimes())[i];
				VarStarti[i+1]  = (*arg.GetFloatStartTimes())[i];
				double DFl        = model.GetZeroCurve()->DiscountPrice(VarStarti[i+1]/K_YEAR_LEN);
				double DFlplus1   = model.GetZeroCurve()->DiscountPrice(VarEndi[i+1]/K_YEAR_LEN);
				double s = (DFl-DFlplus1)*SwapNotional; 
				VarFlow[i]=s;
				DFloat+=s;
			}
		}

		/// get all the zero coupons		
		for(i=0; i<nbVarFlow; ++i)
			DFFloatStarti[i] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[i]/K_YEAR_LEN);
		DFFloatStarti[i]= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);
		
		double Sij = DFloat/(DFix+stubFix);
// FIXMEFRED: mig.vc8 (25/05/2007 15:46:15):cast
		((ARM_VanillaSwaptionArgSFRM&)arg).SetSwapFwd(static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,Sij)));
		double Spj,DFpj,DFixTmp=DFix,Fwdi;
		double Cst0= Sij/(Sij+averageShift);
		double Cst = Cst0/DFloat;
		double deltai, shifti;

		size_t j=0;
		Spj	= Sij;
		DFpj= DFFloatStarti[0];
		double DFkj = 0.0;

		SwapNotional = 1.0;
		for(i=0;i<nbVarFlow;++i)
		{
			if (!isConstantNotional)
				SwapNotional	= (*arg.GetFloatNotional())[i]; 
			
			if( j<arg.GetFixPayTimes()->size() && 
				fabs( (*arg.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				Spj		= DFloat/DFix;
				DFpj    = DFFixStarti[j];
				if (!isConstantNotional)
					DFkj    = DFpj * (*arg.GetFixNotional())[j];
				else
					DFkj	= DFpj;
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			deltai	= (*arg.GetFloatIntTerms())[i];
			Fwdi	= (DFFloatStarti[i]/DFFloatStarti[i+1]-1.0)/deltai;
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);
			(*Mu)[i]	= Cst*deltai*(Fwdi+shifti)
				/(1.0+deltai*Fwdi)*(DFkj+(Sij-Spj)*DFixTmp+alphaFix*Sij*DFFixStarti[1]);
			
			DFloat -= (DFFloatStarti[i]-DFFloatStarti[i+1])*SwapNotional;	
		}
	}

	/////////////////////////mode 1 Forwards expected under terminal probability//////////////////////////
	else if (mode == 1)
	{	
		double asof = model.GetAsOfDate().GetJulian();
		double fromTime = arg.GetEvalTime();
		double toTime = arg.GetExpiry(); 
		ARM_GP_Vector Fwdvectori(nbVarFlow);
		ARM_GP_Vector DFFloatStartAtZeroi(nbVarFlow+1);
		double Fwdi, deltai,shifti;

		for(i=0; i<nbVarFlow; ++i)
			DFFloatStartAtZeroi[i] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[i]/K_YEAR_LEN);
		DFFloatStartAtZeroi[i]= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

		DFFloatStarti[0] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN);


		DFFloatStarti[nbVarFlow] = DFFloatStartAtZeroi[nbVarFlow];
		
		for(i=(nbVarFlow-1);i>=0;--i)
		{
			deltai	= (*arg.GetFloatIntTerms())[i];
			Fwdi	= (DFFloatStartAtZeroi[i]/DFFloatStartAtZeroi[i+1]-1.0)/deltai;
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);			

			double resetTimei = (*arg.GetFloatResetTimes())[i];
			double drift =0.0;
			for(int k=i+1;k<nbVarFlow;++k)
			{
				double resetTimek = (*arg.GetFloatResetTimes())[k];
				double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
				double deltak	= (*arg.GetFloatIntTerms())[k];
				double Fwdk	= Fwdvectori[k];//(DFFloatStartAtZeroi[k]/DFFloatStartAtZeroi[k+1]-1.0)/deltak;
				double shiftk	= ShiftValue(resetTimek);
				double coeffk	=deltak*(Fwdk+shiftk)/(1.0+deltak*Fwdk);
				drift += coeffk*FWDCovar_i_k;
			}			

			double DFi = DFFloatStarti[i]; 
			double expectShiftFwdi= (Fwdi+shifti)*exp(-drift);
			Fwdvectori[i]= expectShiftFwdi-shifti;
			DFFloatStarti[i] = DFFloatStarti[i+1]*(1.0+deltai*Fwdvectori[i]);
		}

		//////////////////////////////////////////
		FixStarti[0]	= (*arg.GetFloatStartTimes())[0];
		DFFixStarti[0]	= model.GetZeroCurve()->DiscountPrice(FixStarti[0]/K_YEAR_LEN);
		int periods = nbVarFlow/nbFixFlow;
		double DFix1 = 0.0;
		for(i=0;i<nbFixFlow;++i)
		{
			double coeffi =1.0;
			for(int k=i*periods;k<periods*(i+1);++k)
			{
				double deltak	= (*arg.GetFloatIntTerms())[k];
				double Fwdk	= Fwdvectori[k];
				coeffi*= (1/(1+deltak*Fwdk));
			}
			DFFixStarti[i+1]= DFFixStarti[i]*coeffi;
			FixFlow[i]		= DFFixStarti[i+1]*(*arg.GetFixPayPeriods())[i];
			DFix1		   += FixFlow[i];
		}
		/////////////////////////////////////////////////////

		DFix = DFix1;

		double DFloat= DFFloatStarti[0]-DFFloatStarti[nbVarFlow];

		double Sij = DFloat/(DFix+stubFix);

		double DFloatAtZero= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN)
			-model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

		((ARM_VanillaSwaptionArgSFRM&)arg).SetSwapFwd(ARM_VectorPtr(new ARM_GP_Vector(1,DFloatAtZero/(DFix+stubFix))));
		double Spj,DFpj,DFixTmp=DFix;
		double Cst0= Sij/(Sij+averageShift);
		double Cst = Cst0/DFloat;
		
		size_t j=0;
		Spj	= Sij;
		DFpj= DFFloatStarti[0];
					
		for(i=0;i<nbVarFlow;++i)
		{
			if( j<arg.GetFixPayTimes()->size() && 
				fabs( (*arg.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				Spj		= DFloat/DFix;
				DFpj    = DFFixStarti[j];
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			deltai	= (*arg.GetFloatIntTerms())[i];
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);
						
			DFloat -= (DFFloatStarti[i]-DFFloatStarti[i+1]);
		
			double DFi = DFFloatStarti[i]; 
			double epsilon_i = ((Sij-Spj)*DFixTmp+alphaFix*Sij*DFFixStarti[1])/DFi+(DFpj-DFi)/DFi;
			double expectFwdi= Fwdvectori[i];
			double expectShiftFwdi= expectFwdi+shifti;

			//approx1
			double poderation_1_i = Cst*deltai*expectShiftFwdi/(1.0+deltai*expectFwdi)*DFi;
			double approx1 = poderation_1_i*(1+epsilon_i);

			//approx2
			double DFiplus1 = DFFloatStarti[i+1]; 
			double poderation_2_i = Cst*deltai*expectShiftFwdi*DFiplus1;
			double approx2 = poderation_2_i*(1+epsilon_i);
			
			(*Mu)[i] = approx1;
		}
	}
	/////////////////////////Forwards expected under the swap probability//////////////////////////
	else 
	{
		double TotalDFix = DFix;
		double DFixTmp = DFix;
		double asof = model.GetAsOfDate().GetJulian();
		double fromTime = arg.GetEvalTime();
		double toTime = arg.GetExpiry();

		ARM_GP_Vector Fwdvectori(nbVarFlow);
		ARM_GP_Vector DFFloatStartAtZeroi(nbVarFlow+1);
		double Fwdi, deltai,shifti;

		for(i=0; i<nbVarFlow; ++i)
			DFFloatStartAtZeroi[i] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[i]/K_YEAR_LEN);
		DFFloatStartAtZeroi[i]= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

		DFFloatStarti[0] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN);

		size_t j=0;	

		//// calcul des drift des fra 
		for(i=0;i<nbVarFlow;++i)
		{
			deltai	= (*arg.GetFloatIntTerms())[i];
			Fwdi	= (DFFloatStartAtZeroi[i]/DFFloatStartAtZeroi[i+1]-1.0)/deltai;
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);			
			
			if( j<arg.GetFixPayTimes()->size() && 
				fabs( (*arg.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			double resetTimei = (*arg.GetFloatResetTimes())[i];
			double drift =0.0;
			for(int k=0;k<i+1;++k)
			{
				double resetTimek = (*arg.GetFloatResetTimes())[k];
				double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
				double deltak	= (*arg.GetFloatIntTerms())[k];
				double Fwdk	= (DFFloatStartAtZeroi[k]/DFFloatStartAtZeroi[k+1]-1.0)/deltak;
				double shiftk	= ShiftValue(resetTimek);
				double coeffk	=deltak*(Fwdk+shiftk)/(1.0+deltak*Fwdk)*(1-DFixTmp/TotalDFix);
				drift += coeffk*FWDCovar_i_k;
			}
			for(k=i+1;k<nbVarFlow;++k)
			{
				double resetTimek = (*arg.GetFloatResetTimes())[k];
				double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
				double deltak	= (*arg.GetFloatIntTerms())[k];
				double Fwdk	= (DFFloatStartAtZeroi[k]/DFFloatStartAtZeroi[k+1]-1.0)/deltak;
				double shiftk	= ShiftValue(resetTimek);
				double coeffk	=deltak*(Fwdk+shiftk)/(1.0+deltak*Fwdk)*(-DFixTmp/TotalDFix);
				drift += coeffk*FWDCovar_i_k;
			}	

			double DFi = DFFloatStarti[i]; 
			double expectShiftFwdi= (Fwdi+shifti)*exp(drift);
			Fwdvectori[i]= expectShiftFwdi-shifti;
		}

		// calcul des coeff vol swap vol fra 
		DFix = TotalDFix;
		double DFloat= DFFloatStartAtZeroi[0]-DFFloatStartAtZeroi[nbVarFlow];

		double Sij = DFloat/(DFix+stubFix);
		//double DFloatAtZero= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN)
			//-model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

		((ARM_VanillaSwaptionArgSFRM&)arg).SetSwapFwd(ARM_VectorPtr(new ARM_GP_Vector(1,Sij)));
		double Spj,DFpj;
		double Cst0= Sij/(Sij+averageShift);
		double Cst = Cst0/DFloat;		
		
		Spj	= Sij;
		DFpj= DFFloatStarti[0];

		for(i=0;i<nbVarFlow;++i)
		{
			if( j<arg.GetFixPayTimes()->size() && 
				fabs( (*arg.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			deltai	= (*arg.GetFloatIntTerms())[i];
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);
			double Fwdi= Fwdvectori[i];
			double expectShiftFwdi= Fwdi+shifti;
			double fwdDFi = DFFloatStartAtZeroi[i+1]/DFFloatStartAtZeroi[i];

			(*Mu)[i] = deltai*expectShiftFwdi*fwdDFi*((DFFloatStartAtZeroi[0]/DFloat)-(1-DFixTmp/TotalDFix));
		}
	}

// FIXMEFRED: mig.vc8 (30/05/2007 16:38:25):cast
	return static_cast<ARM_VectorPtr>(Mu);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ComputeVolSensiVolFra
///	Returns: ARM_VectorPtr
///	Action : computes mui of the volatility of 
///          swaption1 (arg) with expected forwards 
///          and swap1 
///          under the swap2 (argProba) probability 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ModelParamsSFRM::ComputeVolSensiVolFra( 
	            const ARM_VanillaSwaptionArg& arg, 
	            const ARM_PricingModel& model )
{
	int mode = 0;

	size_t nbVarFlow = arg.GetFloatStartTimes()->size();
	size_t nbFixFlow = arg.GetFixPayTimes()->size();

	ARM_GP_Vector* Mu = new ARM_GP_Vector(nbVarFlow);
	double averageShift= AverageShift(*arg.GetFloatResetTimes());

	ARM_GP_Vector DFFloatStarti(nbVarFlow+1);
	int i;

	//Determination of fix leg
	double DFix=0.0;
	ARM_GP_Vector FixFlow(nbFixFlow);
	ARM_GP_Vector FixStarti(nbFixFlow+1);
	ARM_GP_Vector DFFixStarti(nbFixFlow+1);

	FixStarti[0]	= (*arg.GetFloatStartTimes())[0];
	DFFixStarti[0]	= model.GetZeroCurve()->DiscountPrice(FixStarti[0]/K_YEAR_LEN);
	for(i=0;i<nbFixFlow;++i)
	{
		FixStarti[i+1]  = (*arg.GetFixPayTimes())[i];
		DFFixStarti[i+1]= model.GetZeroCurve()->DiscountPrice(FixStarti[i+1]/K_YEAR_LEN);
		FixFlow[i]		= DFFixStarti[i+1]*(*arg.GetFixPayPeriods())[i];
		DFix		   += FixFlow[i];
	}

	double stubFix=0;

	/// see if there is a stub
	double asOfDate  = model.GetAsOfDate().GetJulian();
	ARM_Date startDate( asOfDate+(*arg.GetFloatStartTimes())[0]);
	ARM_Date firstEndDate(asOfDate + (*arg.GetFixPayTimes())[0]);
	ARM_Date stdStartDate(firstEndDate);
	ARM_Currency* pCcy = const_cast<ARM_PricingModel&>(model).GetZeroCurve()->GetCurrencyUnit();
	int fixedFreq = arg.GetFixFrequency();
	stdStartDate.AddPeriod(-fixedFreq,pCcy->GetCcyName());

	double alphaFix=0;
	if(fabs(stdStartDate-startDate)>FRMVOL_LAG_THRESHOLD && nbFixFlow > 1)
	{
		int fixDayCount = pCcy->GetFixedDayCount();
		double delta0	= CountYears(fixDayCount,stdStartDate,startDate);
		alphaFix		= delta0*DFFixStarti[1];
		stubFix			= delta0*(DFFixStarti[1]-DFFixStarti[0]);
	}

	/////////////////////////mode0 Forward At 0//////////////////////////
	if(mode ==0)
	{
		double DFloat= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[0]/K_YEAR_LEN)
			-model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

		/// get all the zero coupons		
		for(i=0; i<nbVarFlow; ++i)
			DFFloatStarti[i] = model.GetZeroCurve()->DiscountPrice((*arg.GetFloatStartTimes())[i]/K_YEAR_LEN);
		DFFloatStarti[i]= model.GetZeroCurve()->DiscountPrice((*arg.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);
		
		double totalFix = DFix+stubFix;
		double Sij = DFloat/(DFix+stubFix);
		((ARM_VanillaSwaptionArgSFRM&)arg).SetSwapFwd(ARM_VectorPtr(new ARM_GP_Vector(1,Sij)));
		double Spj,DFpj,DFixTmp=DFix,Fwdi;
		double Cst0= Sij/(Sij+averageShift);
		double Cst = Cst0/DFloat;
		double deltai, shifti;

		size_t j=0;
		Spj	= Sij;
		DFpj= DFFloatStarti[0];

		for(i=0;i<nbVarFlow;++i)
		{
			if( j<arg.GetFixPayTimes()->size() && 
				fabs( (*arg.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			deltai	= (*arg.GetFloatIntTerms())[i];
			Fwdi	= (DFFloatStarti[i]/DFFloatStarti[i+1]-1.0)/deltai;
			shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);
			(*Mu)[i]	= -deltai*(Fwdi+shifti)/(1.0+deltai*Fwdi)*(DFixTmp/totalFix);						
		}
	}

// FIXMEFRED: mig.vc8 (30/05/2007 16:38:33):cast
	return static_cast<ARM_VectorPtr>(Mu);
}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ComputeVolSwapvolFRAUnderSwapProba
///	Returns: ARM_VectorPtr
///	Action : computes mui of the volatility of 
///          swaption1 (arg) with expected forwards 
///          and swap1 
///          under the swap2 (argProba) probability 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ModelParamsSFRM::ComputeVolSwapvolFRAUnderSwapProba(const ARM_VanillaSwaptionArg& arg,
												 const ARM_PricingModel& model, 
												 const ARM_VanillaSwaptionArg& argProba )
{
	double asOfDate  = model.GetAsOfDate().GetJulian();

	//-------------------------------------determination des forwards sous arg proba-------------------------

	size_t nbVarFlow = argProba.GetFloatStartTimes()->size();
	size_t nbFixFlow = argProba.GetFixPayTimes()->size();	
	double averageShift= AverageShift(*argProba.GetFloatResetTimes());
	ARM_GP_Vector DFFloatStarti(nbVarFlow+1);
	size_t i;

	//Determination of fix leg
	double DFix=0.0;
	ARM_GP_Vector FixFlow(nbFixFlow);
	ARM_GP_Vector FixStarti(nbFixFlow+1);
	ARM_GP_Vector DFFixStarti(nbFixFlow+1);

	FixStarti[0]	= (*argProba.GetFloatStartTimes())[0];
	DFFixStarti[0]	= model.GetZeroCurve()->DiscountPrice(FixStarti[0]/K_YEAR_LEN);
	for(i=0;i<nbFixFlow;++i)
	{
		FixStarti[i+1]  = (*argProba.GetFixPayTimes())[i];
		DFFixStarti[i+1]= model.GetZeroCurve()->DiscountPrice(FixStarti[i+1]/K_YEAR_LEN);
		FixFlow[i]		= DFFixStarti[i+1]*(*argProba.GetFixPayPeriods())[i];
		DFix		   += FixFlow[i];
	}

	double stubFix=0;
	/// see if there is a stub	
	ARM_Date startDate( asOfDate+(*argProba.GetFloatStartTimes())[0]);
	ARM_Date firstEndDate(asOfDate + (*argProba.GetFixPayTimes())[0]);
	ARM_Date stdStartDate(firstEndDate);
	ARM_Currency* pCcy = const_cast<ARM_PricingModel&>(model).GetZeroCurve()->GetCurrencyUnit();
	int fixedFreq = argProba.GetFixFrequency();
	stdStartDate.AddPeriod(-fixedFreq,pCcy->GetCcyName());

	double alphaFix=0;
	if(fabs(stdStartDate-startDate)>FRMVOL_LAG_THRESHOLD && nbFixFlow > 1)
	{
		int fixDayCount = pCcy->GetFixedDayCount();
		double delta0	= CountYears(fixDayCount,stdStartDate,startDate);
		alphaFix		= delta0*DFFixStarti[1];
		stubFix			= delta0*(DFFixStarti[1]-DFFixStarti[0]);
	}

	double TotalDFix = DFix;
	double DFixTmp = DFix;
	double fromTime = argProba.GetEvalTime();
	double toTime = argProba.GetExpiry();

	ARM_GP_Vector Fwdvectori(nbVarFlow);
	ARM_GP_Vector DFFloatStartAtZeroi(nbVarFlow+1);
	double Fwdi, deltai,shifti;

	for(i=0; i<nbVarFlow; ++i)
		DFFloatStartAtZeroi[i] = model.GetZeroCurve()->DiscountPrice((*argProba.GetFloatStartTimes())[i]/K_YEAR_LEN);
	DFFloatStartAtZeroi[i]= model.GetZeroCurve()->DiscountPrice((*argProba.GetFloatEndTimes())[nbVarFlow-1]/K_YEAR_LEN);

	DFFloatStarti[0] = model.GetZeroCurve()->DiscountPrice((*argProba.GetFloatStartTimes())[0]/K_YEAR_LEN);
	size_t j=0;	

	//// calcul des drift des fra 
	for(i=0;i<nbVarFlow;++i)
	{
		deltai	= (*argProba.GetFloatIntTerms())[i];
		Fwdi	= (DFFloatStartAtZeroi[i]/DFFloatStartAtZeroi[i+1]-1.0)/deltai;
		shifti	= ShiftValue((*argProba.GetFloatResetTimes())[i]);			
		
		if( j<argProba.GetFixPayTimes()->size() && 
			fabs( (*argProba.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
		{
			DFixTmp	= DFix;
			DFix   -= FixFlow[j];
			++j;
		}

		double resetTimei = (*argProba.GetFloatResetTimes())[i];
		double drift =0.0;
		for(int k=0;k<i+1;++k)
		{
			double resetTimek = (*argProba.GetFloatResetTimes())[k];
			double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
			double deltak	= (*argProba.GetFloatIntTerms())[k];
			double Fwdk	= (DFFloatStartAtZeroi[k]/DFFloatStartAtZeroi[k+1]-1.0)/deltak;
			double shiftk	= ShiftValue(resetTimek);
			double coeffk	=deltak*(Fwdk+shiftk)/(1.0+deltak*Fwdk)*(1-DFixTmp/TotalDFix);
			drift += coeffk*FWDCovar_i_k;
		}
		for(k=i+1;k<nbVarFlow;++k)
		{
			double resetTimek = (*argProba.GetFloatResetTimes())[k];
			double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
			double deltak	= (*argProba.GetFloatIntTerms())[k];
			double Fwdk	= (DFFloatStartAtZeroi[k]/DFFloatStartAtZeroi[k+1]-1.0)/deltak;
			double shiftk	= ShiftValue(resetTimek);
			double coeffk	=deltak*(Fwdk+shiftk)/(1.0+deltak*Fwdk)*(-DFixTmp/TotalDFix);
			drift += coeffk*FWDCovar_i_k;
		}	

		double expectShiftFwdi= (Fwdi+shifti)*exp(drift);
		Fwdvectori[i]= expectShiftFwdi-shifti;
		DFFloatStarti[i+1] = DFFloatStarti[i]/(1.0+deltai*Fwdvectori[i]);
	}

	//-------------------------------------------calcul du drift du swap2 under swapProba--------------------

	//calcul des flux fix de swap2
	size_t nbVarFlow2 = arg.GetFloatStartTimes()->size();
	size_t nbFixFlow2 = arg.GetFixPayTimes()->size();	
	double averageShift2= AverageShift(*arg.GetFloatResetTimes());

	double DFix2=0.0;
	ARM_GP_Vector FixFlow2(nbFixFlow2);
	ARM_GP_Vector FixStarti2(nbFixFlow2+1);
	ARM_GP_Vector DFFixStarti2(nbFixFlow2+1);

	FixStarti2[0]	= (*arg.GetFloatStartTimes())[0];
	DFFixStarti2[0]	= model.GetZeroCurve()->DiscountPrice(FixStarti2[0]/K_YEAR_LEN);
	for(i=0;i<nbFixFlow2;++i)
	{
		FixStarti2[i+1]  = (*arg.GetFixPayTimes())[i];
		DFFixStarti2[i+1]= model.GetZeroCurve()->DiscountPrice(FixStarti2[i+1]/K_YEAR_LEN);
		FixFlow2[i]		= DFFixStarti2[i+1]*(*arg.GetFixPayPeriods())[i];
		DFix2		   += FixFlow2[i];
	}

	double stubFix2=0;
	/// see if there is a stub	
	ARM_Date startDate2( asOfDate+(*arg.GetFloatStartTimes())[0]);
	ARM_Date firstEndDate2(asOfDate + (*arg.GetFixPayTimes())[0]);
	ARM_Date stdStartDate2(firstEndDate2);
	stdStartDate2.AddPeriod(-fixedFreq,pCcy->GetCcyName());

	double alphaFix2=0;
	if(fabs(stdStartDate-startDate)>FRMVOL_LAG_THRESHOLD && nbFixFlow2 > 1)
	{
		int fixDayCount = pCcy->GetFixedDayCount();
		double delta0	= CountYears(fixDayCount,stdStartDate2,startDate2);
		alphaFix2		= delta0*DFFixStarti2[1];
		stubFix2		= delta0*(DFFixStarti2[1]-DFFixStarti2[0]);
	}

	double TotalDFix2 = DFix2;
	double DFixTmp2 = DFix2;
	j=0;
	
	//la suite du calcul du drift
	ARM_VectorPtr mu2 = ((ARM_VanillaSwaptionArgSFRM&)arg).GetMu();

	double driftSwap2 =0.0;
	for(i=0;i<nbVarFlow;++i)
	{		
		deltai	= (*argProba.GetFloatIntTerms())[i];
		Fwdi	= Fwdvectori[i];
		shifti	= ShiftValue((*argProba.GetFloatResetTimes())[i]);
		double resetTimei = (*argProba.GetFloatResetTimes())[i];
		//calcul de alpha_i
		double alpha_i;
		if(i<nbVarFlow2)
		{
			if( j<argProba.GetFixPayTimes()->size() && 
					fabs( (*argProba.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				DFixTmp2	= DFix2;
				DFix2   -= FixFlow[2];
				++j;
			}

			double deltai	= (*argProba.GetFloatIntTerms())[i];
			alpha_i	= deltai*(Fwdi+shifti)/(1.0+deltai*Fwdi)*(DFixTmp2/TotalDFix2-DFixTmp/TotalDFix);
		}
		else
		{
			if( j<argProba.GetFixPayTimes()->size() && 
					fabs( (*argProba.GetFloatStartTimes())[i]-FixStarti[j] ) < FRMVOL_LAG_THRESHOLD )
			{
				DFixTmp	= DFix;
				DFix   -= FixFlow[j];
				++j;
			}

			double deltai	= (*argProba.GetFloatIntTerms())[i];
			alpha_i	= deltai*(Fwdi+shifti)/(1.0+deltai*Fwdi)*(-DFixTmp/TotalDFix);
		}

		for(int k=0;k<nbVarFlow2;++k)
		{
			double mu_k = (*mu2)[k];
			double resetTimek = (*arg.GetFloatResetTimes())[k];			
			double FWDCovar_i_k = IntegratedCovariance(fromTime,toTime,resetTimei,resetTimek);
			driftSwap2 += alpha_i*mu_k*FWDCovar_i_k;
		}
	}
	
	//-------------------------------------------calcul des Mui de swap2-------------------------------------

	ARM_GP_Vector* Mu = new ARM_GP_Vector(nbVarFlow2);
	// calcul des coeff vol swap vol fra 
	DFix2 = TotalDFix2;
	double DFloat2= DFFloatStarti[0]-DFFloatStarti[nbVarFlow2];

	double Sij = DFloat2/(DFix2+stubFix2);
	((ARM_VanillaSwaptionArgSFRM&)arg).SetSwapFwd(ARM_VectorPtr(new ARM_GP_Vector(1,Sij)));
	double Spj,DFpj;
	double Cst0= Sij*exp(driftSwap2)/(Sij*exp(driftSwap2)+averageShift);
	double Cst = Cst0/DFloat2;		
	
	Spj	= Sij;
	DFpj= DFFloatStarti[0];

	for(i=0;i<nbVarFlow2;++i)
	{
		if( j<arg.GetFixPayTimes()->size() && 
			fabs( (*arg.GetFloatStartTimes())[i]-FixStarti2[j] ) < FRMVOL_LAG_THRESHOLD )
		{
			Spj		= DFloat2/DFix2;
			DFpj    = DFFixStarti2[j];
			DFixTmp2	= DFix2;
			DFix2   -= FixFlow2[j];
			++j;
		}

		deltai	= (*arg.GetFloatIntTerms())[i];
		shifti	= ShiftValue((*arg.GetFloatResetTimes())[i]);
			
		DFloat2 -= (DFFloatStarti[i]-DFFloatStarti[i+1]);
	
		double DFi = DFFloatStarti[i]; 
		double epsilon_i = ((Sij-Spj)*exp(driftSwap2)*DFixTmp+alphaFix*Sij*DFFixStarti[1])/DFi+(DFpj-DFi)/DFi;
		double expectFwdi= Fwdvectori[i];
		double expectShiftFwdi= expectFwdi+shifti;

		//approx1
		double poderation_1_i = Cst*deltai*expectShiftFwdi/(1.0+deltai*expectFwdi)*DFi;
		double approx1 = poderation_1_i*(1+epsilon_i);

		//approx2
		double DFiplus1 = DFFloatStarti[i+1]; 
		double poderation_2_i = Cst*deltai*expectShiftFwdi*DFiplus1;
		double approx2 = poderation_2_i*(1+epsilon_i);
		
		(*Mu)[i] = approx1;
	}
// FIXMEFRED: mig.vc8 (30/05/2007 16:38:38):cast
	return static_cast<ARM_VectorPtr>(Mu);

}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: AverageShift
///	Returns: double
///	Action : computes the average shift over the times
////////////////////////////////////////////////////
double ARM_ModelParamsSFRM::AverageShift(const ARM_GP_Vector& times )
{
	/// computes the averageShift
	double averageShift=0;
	size_t i,Size=times.size();
	for(i=0;i<Size;++i)
		averageShift+= ShiftValue(times[i]);
	averageShift /=Size;
	return averageShift;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: PreComputeFwds
///	Returns: double
///	Action :  computes and stores fwds values 
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::PreComputeFwds(const ARM_SFRM& model,const ARM_Portfolio& shiftConvPort)
{
    ARM_GP_Vector breakPointTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
    size_t fwdsSize = breakPointTimes.size();

    if(!itsFwdValues || fwdsSize != itsFwdValues->size())
    {
        /// computes the corresponding forwards
        delete itsFwdValues;
        itsFwdValues = new ARM_GP_Vector(fwdsSize);

		ARM_IRIndex* pIndex	= GetIRIndex();
			
		//// AAAAAAAAAAAAARRRRRRRRRGGGGGGGGGGG ugly const_cast
        double asOfDate		= model.GetAsOfDate().GetJulian();
		ARM_Currency* pCcy	= const_cast<ARM_SFRM&>(model).GetZeroCurve()->GetCurrencyUnit();
		int floatResetFreq	= pIndex->GetResetFrequency();
		int resetGap		= pIndex->GetResetGap();
        int fwdRule			= pCcy->GetFwdRule();

		char* resetCalendar = pCcy->GetResetCalName( GetDefaultIndexFromCurrency(pCcy->GetCcyName()) ); 
		
		for(size_t i=0; i<fwdsSize; ++i )
		{
			if((&shiftConvPort)
			   && (shiftConvPort.size() == fwdsSize)
			   && (i < fwdsSize) 
			   && dynamic_cast<ARM_CapFloor*>(shiftConvPort.GetAsset(i)))        
				(*itsFwdValues)[i] = ((ARM_CapFloor *) shiftConvPort.GetAsset(i))->GetSwapLeg()->GetFwdRates()->Elt(0)/CC_NS(ARM_Constants,rateBase);
			else if((&shiftConvPort)
			   && (shiftConvPort.size() == fwdsSize)
			   && (i < fwdsSize) 
			   && dynamic_cast<ARM_Swaption *>(shiftConvPort.GetAsset(i)))   
				(*itsFwdValues)[i] = ((ARM_Swaption *) shiftConvPort.GetAsset(i))->GetFloatLeg()->GetFwdRates()->Elt(0)/CC_NS(ARM_Constants,rateBase);
			else if(&model)
			{
				/// If there is no calibration portfolio we use the float leg date 
				/// strip of the model

				ARM_Date fwdResetDate(asOfDate+breakPointTimes[i]);
				ARM_Date fwdStartDate(fwdResetDate);
				fwdStartDate.GapBusinessDay(-resetGap,resetCalendar);
				ARM_Date fwdEndDate(fwdStartDate);
				fwdEndDate.AddPeriod(floatResetFreq);
				fwdEndDate.GoodBusinessDay(fwdRule, resetCalendar);
			
				(*itsFwdValues)[i] = const_cast<ARM_SFRM&>(model).GetZeroCurve()->ForwardYield(fwdStartDate,fwdEndDate,K_MM_RATE,fwdRule)/100.0;
			}
			else
			{
				CC_Ostringstream os;
					os  << " To Pre-Compute, make sur you use a portfolio or pricing model";
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
        
			}
		}
		/// Free memory !!
		delete resetCalendar;        
    }
	
}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ConvertToBetaParam
///	Returns: 
///	Action : convert the beta to a shift param
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::ConvertToBetaParam(const ARM_SFRM& model,const ARM_Portfolio& shiftConvPort)
{
    if( DoesModelParamExist(ARM_ModelParamType::Shift) )
	{
		/// Pre_Compute the fwds values
		PreComputeFwds(model,shiftConvPort);
		
		ARM_GP_Vector breakPointTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
		size_t fwdsSize = breakPointTimes.size();
		ARM_GP_Vector values(fwdsSize);
		
        for(size_t i=0; i<fwdsSize; ++i )
        {
            double shift = ShiftValue(breakPointTimes[i]);	
            double fwd = (*itsFwdValues)[i];
			if( shift < -fwd - K_NEW_DOUBLE_TOL )
				shift = -fwd + K_NEW_DOUBLE_TOL;
			
			values[i]= fwd/(fwd + shift);
		}
        string interpolatorName = ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Shift)).GetCurve()->GetInterpolator()->toString();
        ARM_CurveModelParam betaModelParam(ARM_ModelParamType::Beta, &values,&breakPointTimes,
                                        "Beta", interpolatorName);
		
        /// this is the choice taken in the SFRM model in the kernel
        /// so to match it, we dit it as well!         
        SetModelParam(&betaModelParam);
    }

}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: ConvertToShiftParam
///	Returns: 
///	Action : convert the beta to a shift param
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::ConvertToShiftParam(const ARM_SFRM& model, const ARM_Portfolio& shiftConvPort)
{
    if( DoesModelParamExist(ARM_ModelParamType::Beta) )
	{
		/// Pre_Compute the fwds values
		PreComputeFwds(model,shiftConvPort);
		
		ARM_GP_Vector breakPointTimes(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
		size_t fwdsSize = breakPointTimes.size();
		ARM_GP_Vector values(fwdsSize);
		
        for(size_t i=0; i<fwdsSize; ++i )
        {
            double beta = BetaValue(breakPointTimes[i]);			    
			/// test for nul values!
			if( fabs(beta) < K_NEW_DOUBLE_TOL )
				beta = (beta<0.0? -1:1) * K_NEW_DOUBLE_TOL;
			
			values[i]= ((*itsFwdValues)[i])*(1.0-beta)/beta;
		}
        string interpolatorName = ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Beta)).GetCurve()->GetInterpolator()->toString();
        ARM_CurveModelParam shiftModelParam(ARM_ModelParamType::Shift, &values,&breakPointTimes,
                                        "SHIFT", interpolatorName);
		
        /// this is the choice taken in the SFRM model in the kernel
        /// so to match it, we dit it as well!         
        SetModelParam(&shiftModelParam);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: PostProcessing
///	Returns: 
///	Action : to delete localy variable
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::PostProcessing()
{
    delete itsFwdValues;
    itsFwdValues = NULL;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRM
///	Routine: PostProcessing in calibration
///	Returns: 
///	Action : update the realized correlation
////////////////////////////////////////////////////
void ARM_ModelParamsSFRM::PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model, int factorNb)
{
	ARM_ModelParamVector::const_iterator foundCorrel = modelFitter.FindCalibParamWType( itsCorrelType );
	if( foundCorrel != modelFitter.UnknownCalibParamIterator() ) 
	{
		((ARM_CorrelMatParam&)GetModelParam( itsCorrelType )).UpdateRealizedCorrel();
	}
}


/// Discretization Dates of the model params
ARM_GP_VectorPtr ARM_ModelParamsSFRM::ModelParamsTimeSteps() const
{
	ARM_GP_Vector sigmaTimes(GetVolCurve()->GetAbscisses());
	return ARM_GP_VectorPtr( static_cast<ARM_GP_Vector*> (sigmaTimes.Clone()) );
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
