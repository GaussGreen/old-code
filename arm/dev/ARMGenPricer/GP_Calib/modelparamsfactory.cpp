/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file modelparamsfactory.cpp
 *  \brief
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 *
 */
#include "gpbase/curve.h"
#include "gpbase/checkarg.h"

#include "gpcalib/modelparamsfactory.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpbase/singleton.h"
#include "gpbase/ostringstream.h"
#include "gpbase\numericconstant.h"
#include "gpbase/gpmatrixlinalg.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamFactoryImp
///	Routine: CreateModelParam
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_ModelParam* ARM_ModelParamFactoryImp::CreateModelParam(
	int type,
	ARM_GP_Vector* values,
	ARM_GP_Vector* breakPointTimes,
	const string& name,
    const string& interpolatorName,
	ARM_GP_Vector* lowerbound,
	ARM_GP_Vector* upperbound,
	bool adviseBreakPointTimes,
	const string & currency)
{
	if( values && breakPointTimes)
	{
		if(type == ARM_ModelParamType::Correlation && values->size() > breakPointTimes->size())
		{
			if ( values->size() == breakPointTimes->size() * breakPointTimes->size() )
				/// Time independent correlation matrix between S1(reset=ti) and S2(reset=tj) where Sk=CMS or Libor
				// (ti) is given by breakPointTimes vector
				return new ARM_CorrelMatNoCalibParam( breakPointTimes, values, interpolatorName,lowerbound, upperbound, adviseBreakPointTimes );
			else
			{
				/// Time dependent correlations between N processes Xk
				/// rho(k,l)(ti)=instantaneous correlation at time ti between Xk and Xl
				size_t i,nbTimes = breakPointTimes->size();
				size_t j,nbCurves = (size_t)(floor((double)(values->size())/nbTimes));
// FIXMEFRED: mig.vc8 (22/05/2007 18:03:32): fabs(int) doesnt exist.
				if(abs(nbCurves*breakPointTimes->size() != values->size()))
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : a correlation curve is truncated" );

// FIXMEFRED: mig.vc8 (22/05/2007 18:05:12): sqrt(1.+8.
				double x = 0.5*(1+sqrt(1.+8.*nbCurves));
				size_t k,l,nbFactors = (size_t)(floor(x));

				if(fabs(x-nbFactors)>K_NEW_DOUBLE_TOL)
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : curves are missing to make a correlation matrix" );

				ARM_GP_Vector firstValues(nbTimes),valueAtpoint(nbCurves);
				ARM_GP_T_Vector<ARM_GP_Vector> curvesValues(nbTimes);
				ARM_GP_Matrix correlMatrix(nbFactors,nbFactors);
				ARM_GP_Matrix* eigenVectors;
				ARM_GP_Vector eigenValues(nbFactors);
				for(i=0;i<nbTimes;++i)
				{
					/// The single curve is set to the first correlation curve
					firstValues[i] = (*values)[i*nbCurves];

					for(j=0;j<nbCurves;++j)
						valueAtpoint[j] = (*values)[i*nbCurves+j];

					/// Check that input correlations make a
					/// symetric definite postive matrix
					for(j=0,k=0;k<nbFactors;++k)
					{
						correlMatrix(k,k)=1.0;
						for(l=k+1;l<nbFactors;++l,++j)
						{
							correlMatrix(k,l)=valueAtpoint[j];
							correlMatrix(l,k)=correlMatrix(k,l);
						}
					}
					eigenVectors = ACPTransformation(&correlMatrix,eigenValues,nbFactors);
					for(k=0;k<nbFactors;++k)
					{
						if(eigenValues[k]<0.0)
						{
							delete eigenVectors;
							CC_Ostringstream msge;
							msge << ARM_USERNAME << " : correl inconsistency at time=" << (*breakPointTimes)[i];
							ARM_THROW( ERR_INVALID_ARGUMENT, msge.str());
						}
					}
					delete eigenVectors;

					/// Set multi curves value
					curvesValues[i] = valueAtpoint;
				}
				ARM_CorrelMatParam* correlParam = new ARM_CorrelMatParam( ARM_ModelParamType::Correlation, breakPointTimes, &firstValues, interpolatorName,lowerbound, upperbound, adviseBreakPointTimes );
				ARM_MultiCurve *curves = new ARM_MultiCurve(*breakPointTimes,curvesValues,static_cast<ARM_MultiCurveInterpolator*>(correlParam->GetMultiCurve()->GetInterpolator()->Clone()));
				correlParam->SetMultiCurve(curves); // not cloned !
				return correlParam;
			}

		}
		else if( type == ARM_ModelParamType::BrownianCorrelation )
		{
			return new ARM_BrownianCorrelationParam( breakPointTimes, values, interpolatorName,lowerbound, upperbound, adviseBreakPointTimes );
		}
		else
		{
			CC_NS(ARM_Check,CheckSameArgSize)(*values,*breakPointTimes,"values","breakPointTimes");
				return new ARM_CurveModelParam( (ARM_ModelParamType::ParamNb) type, values, breakPointTimes, name, interpolatorName,lowerbound, upperbound, adviseBreakPointTimes, currency );
		}
	}
	else if(type != ARM_ModelParamType::Correlation && type != ARM_ModelParamType::BrownianCorrelation)
    {
        /// A model param could be built without values nor times (but risky way to create it !)
	    return new ARM_CurveModelParam( (ARM_ModelParamType::ParamNb) type, NULL, NULL, name, interpolatorName,NULL, NULL, adviseBreakPointTimes, currency );
    }
    else
	{
		CC_Ostringstream os;
		os << "break points and values are missing" ;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );

	}
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamFactoryImp
///	Routine: CreateTrigoCorrelMatParam
///	Returns: 
///	Action : With datestrip or reset dates
////////////////////////////////////////////////////

ARM_ModelParam* ARM_ModelParamFactoryImp::CreateTrigoCorrelMatParam( 
	double theta,
	const ARM_DateStrip& datestrip,
	double asOf,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes )
{
	return new ARM_CorrelTrigoMatParam( theta, datestrip, asOf, interpolatorName, lowerBound, upperBound, adviseBreakPointTimes );
}

ARM_ModelParam* ARM_ModelParamFactoryImp::CreateTrigoCorrelMatParam( 
	double theta,
	const ARM_GP_Vector& resetDates,
	double asOf,
	const string& interpolatorName,
	ARM_GP_Vector* lowerBound,
	ARM_GP_Vector* upperBound,
	bool adviseBreakPointTimes )
{
	return new ARM_CorrelTrigoMatParam( theta, resetDates, asOf, interpolatorName, lowerBound, upperBound, adviseBreakPointTimes );
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamFactoryImp
///	Routine: CreateSurfaceModelParam
///	Returns: 
///	Action : Constructor for a surface model param or surface calib model param
////////////////////////////////////////////////////

ARM_SurfaceModelParam* ARM_ModelParamFactoryImp::CreateSurfaceModelParam(
	int modelParamType,
	ARM_Surface* surface,
	const string& modelParamName,
	double lowerBound,
	double upperBound,
	bool adviseBreakPointTimes )
{
	return new ARM_SurfaceModelParam( 
			(ARM_ModelParamType::ParamNb) modelParamType, surface, modelParamName,lowerBound, upperBound, adviseBreakPointTimes  );
}


ARM_SingletonHolder<ARM_ModelParamFactoryImp> ARM_ModelParamFactory;


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
