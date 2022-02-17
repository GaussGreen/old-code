/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsSFRMRow.cpp
 *  \brief From ROW of model params for SFRM model!
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/ModelParamsSFRMRow.h"

#include "gpbase/ostringstream.h"
#include "gpbase/checkarg.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/gpvector.h"

#include "gpcalib/modelfitter.h"

#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/pricingmodel.h"

#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSFRMRow::ARM_ModelParamsSFRMRow( const ARM_ModelParamVector& params, 
       ARM_IRIndex* index,
       size_t factorsNb)
:	ARM_ModelParamsSFRM( params, index, factorsNb ), itsSquaredIntegratedVol()
{
    InitSquaredIntegratedVol();
    PreComputeSquaredIntegratedVol();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: InitSquaredIntegratedVol
///	Returns: void
///	Action : initializes the squared integrated vol (size)
////////////////////////////////////////////////////
void ARM_ModelParamsSFRMRow::InitSquaredIntegratedVol()
{
	/// Initializes the squaredIntegratedVol
	std::vector<double> times(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
    std::vector<double> values(times.size());
	itsSquaredIntegratedVol = ARM_CurvePtr(new ARM_Curve(times,values, new ARM_LinInterpCstExtrapolDble ));
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: PreComputeSquaredIntegratedVol
///	Returns: void
///	Action : computes the squared integrated vol at break point times
////////////////////////////////////////////////////
void ARM_ModelParamsSFRMRow::PreComputeSquaredIntegratedVol()
{
	/// Computes the squared integrated Vol
	std::vector<double> timesValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
	std::vector<double> volValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates());

#if defined(__GP_STRICT_VALIDATION)
	CC_NS(ARM_Check,CheckSameArgSize)(volValues,timesValues,"volValues","timesValues");
#endif


    /// computes Sum(0,i) exp(2*meanRev*u)*sigma^2(u)*du
	size_t volSize=volValues.size();
	std::vector<double> squaredVol(volSize);
	double value;

	double time = 0.0;
	for(size_t i=0;i<volSize;++i)
	{
        double nextTime = timesValues[i]/K_YEAR_LEN;
		double integrate= IntegrateSquaredExpMRV(time,nextTime);
		value = volValues[i]*volValues[i]*integrate;
		time = nextTime;
#if defined(__GP_STRICT_VALIDATION)
		if(value<0.0)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "negative vol squared!");
#endif
		if(i>0)
			value+=squaredVol[i-1];
		squaredVol[i]=value;
	}

	itsSquaredIntegratedVol->SetOrdinates(squaredVol);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: PostProcessing
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_ModelParamsSFRMRow::PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model , int FactorNb)
{
	ARM_ModelParamVector::const_iterator foundMeanRev = modelFitter.FindCalibParamWType( ARM_ModelParamType::MeanReversion );
	ARM_ModelParamVector::const_iterator foundVol = modelFitter.FindCalibParamWType( ARM_ModelParamType::Volatility );
	
	if(    (foundMeanRev!= modelFitter.UnknownCalibParamIterator() ) 
		||   (foundVol!= modelFitter.UnknownCalibParamIterator()) )
		PreComputeSquaredIntegratedVol();

	ARM_ModelParamsSFRM::PostProcessing(modelFitter,model,FactorNb);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_ModelParams::iterator ARM_ModelParamsSFRMRow::SetModelParamValue( int paramType, 
	size_t i,
	double value, 
	double time,
	double tenor )
{
    ARM_ModelParams::iterator found = ARM_ModelParamsSFRM::SetModelParamValue(paramType,i,value, time );

    /// should update the vol values
    if( paramType == ARM_ModelParamType::Volatility )
        PreComputeSquaredIntegratedVol();
    return found;
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: Copy constructor, assignment operator, destructor
///	Returns: 
///	Action : copy the object
////////////////////////////////////////////////////
ARM_ModelParamsSFRMRow::ARM_ModelParamsSFRMRow( const ARM_ModelParamsSFRMRow& rhs )
:	ARM_ModelParamsSFRM(rhs), 
	itsSquaredIntegratedVol( (ARM_Curve*) rhs.itsSquaredIntegratedVol->Clone() )
{}


ARM_ModelParamsSFRMRow& ARM_ModelParamsSFRMRow::operator=( const ARM_ModelParamsSFRMRow& rhs )
{
	if( this != &rhs )
	{
		ARM_ModelParamsSFRM::operator =(rhs);
		itsSquaredIntegratedVol = ARM_CurvePtr( (ARM_Curve*) rhs.itsSquaredIntegratedVol->Clone() );
	}
	return *this;
}

ARM_ModelParamsSFRMRow::~ARM_ModelParamsSFRMRow()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: InstantaneousVolatility
///	Returns: double
///	Action : computes the instantaneous Volatility
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::VolatilityFunction(double t, double T) const
{
    double exp = (T/K_YEAR_LEN)< K_NEW_DOUBLE_TOL ? 1.0: ExpFunction(-T/K_YEAR_LEN);
    return exp * GetVolCurve()->Interpolate(t);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: IntegratedLocalVariance
///	Returns: double
///	Action : computes the integrated vol
/// basically Sum(s,t) sigma(u)^2*exp(2*MeanRev*u) du
///	and with itsSquaredIntegratedVol containing the Function Sum(0,ti) for ti
/// the break point time of the volatility
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::IntegratedLocalVariance(double s, double t) const
{
	//// get time in increasing order!
    if(s>=t)
		CC_NS(std,swap)(s,t);

	/// integrates sum from s to t of sigma(u)^2*exp(2*MeanRev*u)du
	/// Chales relation: sum(s,t) = Sum(0,t)-Sum(0,s)
	double sqr = IntegratedLocalVarianceFromZero(t)-IntegratedLocalVarianceFromZero(s);

	if(fabs(sqr) < K_NEW_DOUBLE_TOL )
		sqr=0;

#if defined(__GP_STRICT_VALIDATION)
	if ( sqr < 0.0 )	    
	   throw Exception(__LINE__, __FILE__,ERR_INVALID_ARGUMENT,"The integral of Sigma curve must be increasing ");
#endif
	
	return(sqr);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: IntegratedLocalVariance
///	Returns: double
///	Action : integrates sum(0,t) of sigma(u)^2*exp(2*MeanRev*u)du
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::IntegratedLocalVarianceFromZero(double t) const
{
    if (t<K_NEW_DOUBLE_TOL)
            return  0.0;

	double vol,Ti,sqrtt, sqrtTi;
	size_t Size = GetVolCurve()->size();

    std::vector<double> vect(GetVolCurve()->GetAbscisses());
	
	/// to get exactly the index lower than the index, we need to substract one!
	int index = lower_boundPosWithPrecision(vect,t,K_FRM_TOL)-1;
	
	/// t<t0?
	/// if so Sum(0,t) = 
	///		-if( MeanRev==0 ) Sigma(0)^2*t
	///		-else 1.0/(2.0*MeanRev) * Sigma(0)^2 * (exp(2.0*MeanRev*t)-1.0)
    if(index == -1)
	{
	    vol = (GetVolCurve()->GetOrdinates())[0];
	    sqrtt = vol*vol*IntegrateSquaredExpMRV(0.0, t/K_YEAR_LEN);
	}
	/// t>t0 : Use SquaredIntegratedVol from 0 to Ti and adds the sum from Ti to t
	else
	{
        double sqrTi = (itsSquaredIntegratedVol->GetOrdinates())[index];
		Ti			 = (GetVolCurve()->GetAbscisses())[index];

		/// t==Ti up to K_NEW_DOUBLE_TOL precision?
		if( fabs(t-Ti)<K_NEW_DOUBLE_TOL )
			sqrtt = sqrTi;	    

		/// otherwise
		/// Get the next one...
		/// and compute the remaining portion with same formula as above!
		/// sqrtt = sqrTi + sqrTitot (sqrtTi)
		else
		{
            vol		= GetVolCurve()->Interpolate(t);
			sqrtTi	= vol*vol*(IntegrateSquaredExpMRV(Ti/K_YEAR_LEN, t/K_YEAR_LEN));
			sqrtt	= sqrTi + sqrtTi;
		}
	}
	return sqrtt;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: MaturityTerm
///	Returns: 
///	Action : exp(-mean reversion * T/K_YEAR_LEN)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::MaturityTerm(double T) const
{
	return T<K_NEW_DOUBLE_TOL? 1.0: ExpFunction(-T/K_YEAR_LEN);
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: MaturityTermSquared
///	Returns: 
///	Action : exp(-2.*mean reversion * T/K_YEAR_LEN)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::MaturityTermSquared(double T) const
{
	return T<K_NEW_DOUBLE_TOL? 1.0: ExpFunction(-2.*T/K_YEAR_LEN);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: VolatilitySpotSquared
///	Returns: double
///	Action : exp(2*MeanRev*s)*Sigma^2(s)
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::VolatilitySpotSquared(double s) const
{
    double exp2	= (s/K_YEAR_LEN)< K_NEW_DOUBLE_TOL ? 1.0: ExpFunction(2 * s/K_YEAR_LEN);
	double sigma=GetVolCurve()->Interpolate(s);
    return exp2 * sigma*sigma;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: VarianceToTime
///	Returns: double
///	Action : return t such that var(t)=var
////////////////////////////////////////////////////
double ARM_ModelParamsSFRMRow::VarianceToTime(double var,double minTime,double maxTime) const
{
    /// Convert the input variance to a variance per factor (identical for each factor)
    var /= FactorCount();

    double mrs = (( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates())[0];
    double vol2,scale = 2*mrs;

    const std::vector<double>& times = GetVolCurve()->GetAbscisses();
    const std::vector<double>& vars = itsSquaredIntegratedVol->GetOrdinates();
	int index = lower_boundPosWithPrecision(vars,var,K_NEW_DOUBLE_TOL)-1;
    if(index == -1)
    {
        /// t < times[0]
        vol2 = (GetVolCurve()->GetOrdinates())[0];
        vol2 *= vol2;
        if( fabs(mrs) < K_NEW_DOUBLE_TOL )
            return var/vol2*K_YEAR_LEN;
        else
            return log(scale*var/vol2+1)/scale*K_YEAR_LEN;
    }
    else
    {
		double time = times[index];
        double varErr = var-vars[index];
		if( varErr < - K_NEW_DOUBLE_TOL)
        {
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+ " : Can't inverse variance to get time" );
        }
        else if( varErr <= K_NEW_DOUBLE_TOL )
            /// t = times[index]
			return 	time;
        else
		{
            /// t > times[index]
            vol2 = GetVolCurve()->Interpolate(time + 10*K_NEW_DOUBLE_TOL);
            vol2 *= vol2;
            if( fabs(mrs) < K_NEW_DOUBLE_TOL )
                return time + varErr/vol2*K_YEAR_LEN;
            else
                return log(scale*varErr/vol2+exp(scale*time/K_YEAR_LEN))/scale*K_YEAR_LEN;
		}
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: Clone
///	Returns: 
///	Action : Standard ARM_Object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsSFRMRow::Clone() const
{
	return new ARM_ModelParamsSFRMRow(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ModelParamsSFRMRow::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n =======> Model Params SFRM Row <====== \n";
    os << "---------------------\n\n";
    os << ARM_ModelParams::toString();

	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
