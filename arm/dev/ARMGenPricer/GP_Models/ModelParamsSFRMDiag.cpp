/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsSFRMDiag.cpp
 *	\brief
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 *
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/ModelParamsSFRMDiag.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/gpvector.h"
#include "gpbase/gpmatrix.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/pricingmodel.h"

/// gpcalib
#include "gpcalib/bootstrap1d.h"

/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSFRMDiag::ARM_ModelParamsSFRMDiag( const ARM_ModelParamVector& params, ARM_IRIndex* index, size_t factorsNb)
:	ARM_ModelParamsSFRM( params, index, factorsNb )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: PostProcessing
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_ModelParamsSFRMDiag::PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model , int factorNb)
{
	ARM_ModelParamsSFRM::PostProcessing(modelFitter,model,factorNb);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: PreProcessing
///	Returns: void
///	Action : initializes the squared integrated vol (size)
////////////////////////////////////////////////////
void ARM_ModelParamsSFRMDiag::PreProcessing(ARM_ModelFitter& modelFitter,int factorNb)
{
    /// we put a typeid to implement in a sense a double dispatcher...
    /// we did not use a dynamic to avoid throwing exception of type std::bad_cast
    if( typeid(modelFitter) == typeid(ARM_Bootstrap1D) && !modelFitter.GetPortfolio().IsNull() && modelFitter.GetPortfolio()->GetSize() &&
        modelFitter.GetPortfolio()->GetAsset(0)->GetName() == ARM_SWAPTION)
    {
        modelFitter.SetCalibDirection(CalibDirection_Backward);
    }

}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: InstantaneousVolatility
///	Returns: double
///	Action : computes the instantaneous Volatility
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::VolatilityFunction(double t, double T) const
{
    double exp = (fabs((T-t)/K_YEAR_LEN))< K_NEW_DOUBLE_TOL ? 1.0: ExpFunction(-(T-t)/K_YEAR_LEN);
    return exp * GetVolCurve()->Interpolate(T);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: IntegratedLocalVariance
///	Returns: double
///	Action : computes the integrated vol
/// basically Sum(s,t) exp(2*MeanRev*u) du
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::IntegratedLocalVariance(double s, double t) const
{
	//// get time in increasing order!
    if(s>=t)
		CC_NS(std,swap)(s,t);
	return IntegrateSquaredExpMRV(s/K_YEAR_LEN,t/K_YEAR_LEN);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: MaturityTerm
///	Returns: double
///	Action : exp(-MeanRev*T)*Sigma(T)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::MaturityTerm(double T) const
{
    double exp = (T/K_YEAR_LEN)< K_NEW_DOUBLE_TOL ? 1.0: ExpFunction(-T/K_YEAR_LEN);
    return exp * GetVolCurve()->Interpolate(T);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: MaturityTermSquared
///	Returns: double
///	Action : exp(-2*MeanRev*T)*Sigma^2(T)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::MaturityTermSquared(double T) const
{
    double exp2	= (T/K_YEAR_LEN)< K_NEW_DOUBLE_TOL ? 1.0: ExpFunction(-2 * T/K_YEAR_LEN);
	double sigma=GetVolCurve()->Interpolate(T);
    return exp2 * sigma*sigma;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: VolatilitySpotSquared
///	Returns: double
///	Action : exp(2*MeanRev*s)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::VolatilitySpotSquared(double s) const
{
    return s<K_NEW_DOUBLE_TOL? 1.0: ExpFunction(2.*s/K_YEAR_LEN);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: VarianceToTime
///	Returns: double
///	Action : return t such that var(t)=var
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMDiag::VarianceToTime(double var,double minTime,double maxTime) const
{
    /// Convert the input variance to a variance per factor (identical for each factor)
    var /= FactorCount();

    double mrs = (( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates())[0];
    if( fabs(mrs) < K_NEW_DOUBLE_TOL )
        return var*K_YEAR_LEN;
    else
    {
        double scale = 2*mrs;
        return log(scale*var+1)/scale*K_YEAR_LEN;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: Copy constructor, assignment operator and destructor
///	Returns: 
///	Action : copy the object
////////////////////////////////////////////////////
ARM_ModelParamsSFRMDiag::ARM_ModelParamsSFRMDiag( const ARM_ModelParamsSFRMDiag& rhs )
:	ARM_ModelParamsSFRM(rhs)
{}


ARM_ModelParamsSFRMDiag& ARM_ModelParamsSFRMDiag::operator=( const ARM_ModelParamsSFRMDiag& rhs )
{
	if( this != &rhs )
		ARM_ModelParamsSFRM::operator =(rhs);
	return *this;
}

ARM_ModelParamsSFRMDiag::~ARM_ModelParamsSFRMDiag()
{}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsSFRMDiag
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsSFRMDiag::Clone() const
{
	return new ARM_ModelParamsSFRMDiag(*this);
}

///////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMDiag
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ModelParamsSFRMDiag::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n =======> Model Params SFRM Diag <====== \n";
    os << "---------------------\n\n";
    os << ARM_ModelParams::toString();
	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
