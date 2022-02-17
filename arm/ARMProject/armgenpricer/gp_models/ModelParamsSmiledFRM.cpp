/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// gpmodels
#include "gpmodels/ModelParamsSmiledFRM.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/vectormanip.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/numericconstant.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"




/// standard libraries
#include <cmath>




CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSmiled::ARM_ModelParamsSmiled(	const ARM_ModelParamVector& params, 
													size_t dim, 
													CorrelType correlType, 
													bool allowInterpol,
													double recorrel)
:	ARM_ModelParams(params),
	itsVol(0),
	itsVolCalib(0),
	itsResetDates(0),
	itsGlobalCovariance(0),
	itsGlobalCovarianceCalib(0),
	itsGlobalCovTimes(0),
	itsNbFwds(0),
	itsNbFactors(dim),
	itsIndexFrom(0),
	itsIndexTo(0),
	itsCorrelType(correlType),
	itsAllowInterpol(allowInterpol),
	itsUseModifiedCovForCalib(true)
{
	if (!UsingRecorrel())
	{
		std::vector<double> abs(1,0.0);
		std::vector<double> ord(1,recorrel);

		ARM_CurveModelParam reCorrelModeParam(ARM_ModelParamType::ReCorrelation, &ord, &abs);

		SetModelParam(&reCorrelModeParam,ARM_ModelParamType::ReCorrelation);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsSmiled::ARM_ModelParamsSmiled( const ARM_ModelParamsSmiled& rhs )
:	ARM_ModelParams(rhs),
	itsResetDates(rhs.itsResetDates),
	itsGlobalCovTimes(CreateClone(rhs.itsGlobalCovTimes)),
	itsNbFwds(rhs.itsNbFwds),
	itsNbFactors(rhs.itsNbFactors),
	itsIndexFrom(rhs.itsIndexFrom),
	itsIndexTo(rhs.itsIndexTo),
	itsCorrelType(rhs.itsCorrelType),
	itsAllowInterpol(rhs.itsAllowInterpol),
	itsUseModifiedCovForCalib(rhs.itsUseModifiedCovForCalib),
	itsCorrelMatrix(rhs.itsCorrelMatrix)
{
	DuplicateCloneablePointorVectorInPlace<ARM_Curve>( rhs.itsVol, itsVol );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsGlobalCovariance, itsGlobalCovariance );
	DuplicateCloneablePointorVectorInPlace<ARM_Curve>( rhs.itsVolCalib, itsVolCalib );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsGlobalCovarianceCalib, itsGlobalCovarianceCalib );
}



////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsSmiled::~ARM_ModelParamsSmiled()
{
	DeletePointorVector<ARM_Curve>( itsVol );
	DeletePointorVector<ARM_GP_Matrix>( itsGlobalCovariance );
	DeletePointorVector<ARM_Curve>( itsVolCalib );
	DeletePointorVector<ARM_GP_Matrix>( itsGlobalCovarianceCalib );

	delete itsGlobalCovTimes;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: setResetDates
///	Returns: 
///	Action : set the reset dates
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::setVolatilitiesAndCorrelMatrix(ARM_GP_Matrix& volatilities, ARM_MatrixVector& matrixVector )
{
	itsVolatilities = CreateClonedPtr(&volatilities);
	itsCorrelMatrix = matrixVector;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: setVolatilities
///	Returns: 
///	Action : set volatilities
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::setVolatilities(ARM_GP_Matrix& volatilities)
{
	itsVolatilities = CreateClonedPtr(&volatilities);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: setResetDates
///	Returns: 
///	Action : set the reset dates
////////////////////////////////////////////////////
void ARM_ModelParamsSmiled::setResetDates(const ARM_GP_VectorPtr& v)
{
	itsResetDates = v;

	if ((itsCorrelMatrix.size() != 0) && (itsCorrelMatrix.size() != itsResetDates->size()))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "correl matrix and reset dates should have the same size!");


	itsNbFwds = v->size();
	if (itsNbFwds<itsNbFactors)
		itsNbFactors = itsNbFwds;

	BuildVolatilityStructure();


	if ((itsCorrelType==CorrelMatrix)||(itsCorrelType==Fwd))
	{
		UpdateVolatilityFromVolatilites();
	}
	else
	{
		UpdateVolatilityFromHump();
	}

	PreComputeIntegratedCovariance();

}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetCorrel
///	Returns: 
///	Action : Returns instantaneous correlation 
///		between gaussian factors driving Libor i and Libor j
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::InstantCorrel(double t,size_t i,size_t j,size_t index) const
{
	double beta = GetCorrelCurve()->Interpolate(t);
	double dt	= ((*itsResetDates)[i]-(*itsResetDates)[j])/K_YEAR_LEN;
	double dT	= ((*itsResetDates)[itsNbFwds-1]-(*itsResetDates)[0])/K_YEAR_LEN;
	double correl = 1.;
	double recorrel = 0.0;
	if (UsingRecorrel())
		recorrel = GetReCorrelCurve()->Interpolate(t);

	if (itsCorrelType==Beta)
	{
		if (!UsingRecorrel())
			correl = exp(-beta*fabs(dt));
		else
		{
			
			double delta = ( (*itsResetDates)[i] + (*itsResetDates)[j] - 2. * (*itsResetDates)[index] ) / K_YEAR_LEN;
			if (delta<K_DOUBLE_TOL)
				correl = exp(-beta*fabs(dt));
			else
				correl = exp(-beta*fabs(dt)/pow(delta,recorrel));
		}
	}
	else if (itsCorrelType==Theta)
	{
		double ratio = itsNbFwds==1?1:fabs(dt)/dT;	
		correl = cos(beta*ratio*ARM_NumericConstants::ARM_PI);	
	}
	else if (itsCorrelType==CorrelMatrix)
	{
		correl = (*itsCorrelMatrix[index])(i,j);
	}
	return correl;
}

double ARM_ModelParamsSmiled::InstantCorrel(double t, int i, int j) const
{
	double beta = GetCorrelCurve()->Interpolate(t);
	double dt	= ((*itsResetDates)[i]-(*itsResetDates)[j])/K_YEAR_LEN;
	double dT	= ((*itsResetDates)[itsNbFwds-1]-(*itsResetDates)[0])/K_YEAR_LEN;
	double correl = 1.;
	double recorrel = 0.0;
	if (UsingRecorrel())
		recorrel = GetReCorrelCurve()->Interpolate(t);

	if (itsCorrelType==Beta)
	{
		if (!UsingRecorrel())
			correl = exp(-beta*fabs(dt));
		else
		{
			
			double delta = ( (*itsResetDates)[i] + (*itsResetDates)[j] - 2. * t ) / K_YEAR_LEN;
			if (delta<K_DOUBLE_TOL)
				correl = exp(-beta*fabs(dt));
			else
				correl = exp(-beta*fabs(dt)/pow(delta,recorrel));
		}
	}
	else
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+"ARM_ModelParamsSmiled::InstantCorrel(t,i,j) not defined if correl type != from Beta");
	}

	return correl;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: BuildVolatilityStructure
///	Returns: 
///	Action : build volatility structure (step up right at reset dates)
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::BuildVolatilityStructure()
{
	//in case of previous computing
	DeletePointorVector<ARM_Curve>( itsVol );
	itsVol.resize(itsNbFwds);

	size_t dec	= (*itsResetDates)[0]>K_DOUBLE_TOL?1:0;
	size_t n	= itsNbFwds + dec;

	std::vector<double> dates(n,0.);
	std::vector<double> values(n,1.);
	for (size_t k = dec ; k < n ; k++ )
		dates[k] = (*itsResetDates)[k-dec];

	for (k=0;k<itsNbFwds;k++)
	{
		ARM_StepUpRightOpenCstExtrapolDble* interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
		itsVol[k] = new ARM_Curve(dates ,values, interpolator);
	}

	//in case of previous computing
	DeletePointorVector<ARM_Curve>( itsVolCalib );
	itsVolCalib.resize(itsNbFwds);

	for (k=0;k<itsNbFwds;k++)
	{
		ARM_StepUpRightOpenCstExtrapolDble* interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
		itsVolCalib[k] = new ARM_Curve(dates ,values, interpolator);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: StationaryConstraint
///	Returns: 
///	Action : Integral of h² between a and b
////////////////////////////////////////////////////

double ARM_ModelParamsSmiled::StationaryConstraint(double a, double b)
{
	if (a<0.)
		return 0.;

	ARM_Curve* pcurve = ( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Hump) ).GetCurve();
	const std::vector<double>& x = pcurve->GetAbscisses();
	const std::vector<double>& y = pcurve->GetOrdinates();

//Localization of a and b
	size_t index_a = IndexOfFirstHigherInVector_DefaultLast(a,x);
	size_t index_b = IndexOfFirstHigherInVector_DefaultLast(b,x);
	int inter   = index_b - index_a;
	int k;
	
	std::vector<double> xx(inter+2,0.);
	std::vector<double> yy(inter+2,1.);

	xx[0] = a; yy[0] = pcurve->Interpolate(a);
	for (k=0;k<inter;k++)
	{
		xx[k+1] = x[index_a + k];
		yy[k+1] = y[index_a + k];
	}
	xx[inter+1] = b;yy[inter+1] = pcurve->Interpolate(b);

//We now have to integrate a quadratic function of xx -> simpson
	double sum=0.,yymid;
	for (k=0;k<inter+1;k++)
	{
		yymid =	(yy[k]+yy[k+1])*(yy[k]+yy[k+1]);
		sum	 +=	(yy[k]*yy[k]+yy[k+1]*yy[k+1]+yymid)*(xx[k+1]-xx[k])/6.;
	}
	
	#if defined(__GP_STRICT_VALIDATION)
	if ( sum < 0.0 )	    
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled::StationaryConstraint : The integral should be positive");
	#endif
	
	return sum/(b-a);		
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: UpdateVolatilityFromHump
///	Returns: 
///	Action : updates volatility structure
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::UpdateVolatilityFromHump()
{
	if( DoesModelParamExist(ARM_ModelParamType::Hump) )
	{
		SBGM_VolCurveVector::iterator iterVol		= itsVol.begin();
		SBGM_VolCurveVector::iterator iterVolCalib	= itsVolCalib.begin();
		std::vector<double>::iterator iterReset			= itsResetDates->begin();

		for (;iterVol!=itsVol.end();++iterVol,++iterVolCalib,++iterReset)
		{
			double sum1=0.,sum2=0;
			std::vector<double> value	= (*iterVol)->GetOrdinates();
			std::vector<double> date	= (*iterVol)->GetAbscisses();
			for (size_t k=0 ; k<value.size()-1 ; k++ )
			{
				double from	= ((*iterReset)-date[k+1])/K_YEAR_LEN;
				double to	= ((*iterReset)-date[k])/K_YEAR_LEN;
				value[k+1] = sqrt( StationaryConstraint(from,to) );
				if (from>=0.)
				{
					sum1+=value[k+1]*value[k+1]*(to-from);
					sum2+=to-from;
				}
			}
			if (sum2>0.)
			{
				value[0]=0.;
				for ( k=0 ; k<value.size()-1 ; k++ )
					value[k+1]*=sqrt(sum2/sum1);
			}
			else
			{
				value[0]=0.;
				for ( k=0 ; k<value.size()-1 ; k++ )
					value[k+1]*=0.;
			}

			if (itsAllowInterpol) 
				for ( k=0 ; k<value.size()-1 ; k++ ) 
					value[k+1]=(value[k+1]>0.?value[k+1]:1.);

			(*iterVol)->SetOrdinates(value);
			(*iterVolCalib)->SetOrdinates(value);
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: UpdateVolatilityFromHump
///	Returns: 
///	Action : updates volatility structure
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::UpdateVolatilityFromVolatilites()
{
	SBGM_VolCurveVector::iterator iterVol		= itsVol.begin();
	SBGM_VolCurveVector::iterator iterVolCalib	= itsVolCalib.begin();
	std::vector<double>::iterator iterReset			= itsResetDates->begin();

	double step = 0.0, nextStep;

	for (size_t i=0;iterVol!=itsVol.end();++iterVol,++iterVolCalib,++iterReset,++i)
	{
		nextStep = (*iterReset);
		std::vector<double> value	= (*iterVol)->GetOrdinates();

		for ( size_t k=0 ; k<value.size()-1 ; k++ ) 
			value[k+1]=(*itsVolatilities)(k,i);

		(*iterVol)->SetOrdinates(value);
		(*iterVolCalib)->SetOrdinates(value);

		step = nextStep;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_ModelParams::iterator ARM_ModelParamsSmiled::SetModelParamValue( int paramType, 
	size_t i,
	double value, 
	double time,
	double tenor,
	size_t factorNb)
{
    if( paramType == ARM_ModelParamType::BetaCorrelation )
	{
		ARM_ModelParams::iterator found = ARM_ModelParams::SetModelParamValue(paramType, i, value, time);
		if (!itsUseModifiedCovForCalib) UpdateIntegratedCovarianceModel();
		if (itsUseModifiedCovForCalib)	UpdateIntegratedCovarianceCalib();
	    return found;
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled::SetModelParamValue : should be called for beta correl only ");
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: AdviseCurrentCalib
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::AdviseCurrentCalib(double from, double to)
{
	itsIndexFrom = IndexOfLastLowerEqInVector(from,*itsGlobalCovTimes);
	itsIndexTo = IndexOfFirstHigherEqInVector_DefaultLast(to,*itsGlobalCovTimes);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: ForOptimize
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::ForOptimize()
{
	if (!itsUseModifiedCovForCalib) UpdateIntegratedCovarianceModel();
	if (itsUseModifiedCovForCalib)	UpdateIntegratedCovarianceCalib();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: PreComputeIntegratedCovariance
///	Returns: 
///	Action : computes the squared integrated covariances
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::PreComputeIntegratedCovariance()
{
	size_t nbProcs = itsNbFwds;
	
	if (itsGlobalCovTimes) delete itsGlobalCovTimes;
	itsGlobalCovTimes = MergeSortedVectorNoDuplicates( GetVolCurve(0)->GetAbscisses(), GetVolCurve(0)->GetAbscisses());

	size_t nbSteps = itsGlobalCovTimes->size();
	
	DeletePointorVector<ARM_GP_Matrix>( itsGlobalCovariance );
	itsGlobalCovariance.resize(nbSteps);
	
	for(size_t index=0;index<nbSteps ;++index)
	{
		itsGlobalCovariance[index] = new ARM_GP_Matrix(itsNbFwds,itsNbFwds);
	}

	DeletePointorVector<ARM_GP_Matrix>( itsGlobalCovarianceCalib );
	itsGlobalCovarianceCalib.resize(nbSteps);
	
	for(index=0;index<nbSteps ;++index)
		itsGlobalCovarianceCalib[index] = new ARM_GP_Matrix(itsNbFwds,itsNbFwds);

	AdviseCurrentCalib( (*itsGlobalCovTimes)[0], (*itsGlobalCovTimes)[nbSteps-1] );
	UpdateIntegratedCovarianceModel();
	UpdateIntegratedCovarianceCalib();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: UpdateIntegratedCovarianceModel
///	Returns: 
///	Action : computes the squared integrated covariances at break point times
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::UpdateIntegratedCovarianceModel()
{
	if (!itsGlobalCovTimes)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,UpdateIntegratedCovarianceModel : covariance structure not initialized!! ");

	double voli,volj;
	size_t i,j;
	
	for ( size_t index = itsIndexFrom ; index <= itsIndexTo ; index++ )
	{
		if (index==0){
			for (i=0;i<itsNbFwds;i++)
			{
				for (j=0;j<=i;j++)
				{
					(*itsGlobalCovariance[index])(i,j)=0.0;
					(*itsGlobalCovariance[index])(j,i)=0.0;
				}
			}
		}else{
			//ne marche que pour beta stepupright
			size_t volIndex = IndexOfFirstHigherEqInVector_DefaultLast((*itsGlobalCovTimes)[index],GetVolCurve(0)->GetAbscisses());
			for (i=0;i<itsNbFwds;i++){
				voli = GetVolCurve(i)->GetOrdinates()[volIndex];
				for (j=0;j<=i;j++){
					volj = GetVolCurve(j)->GetOrdinates()[volIndex];
					(*itsGlobalCovariance[index])(i,j) = (*itsGlobalCovariance[index-1])(i,j)
						+ InstantCorrel((*itsGlobalCovTimes)[index],i,j,index-1) * voli * volj * ((*itsGlobalCovTimes)[index]-(*itsGlobalCovTimes)[index-1])/K_YEAR_LEN;
					(*itsGlobalCovariance[index])(j,i) = (*itsGlobalCovariance[index])(i,j);
				}
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: UpdateIntegratedCovarianceCalib
///	Returns: 
///	Action : computes the squared integrated covariances at break point times
////////////////////////////////////////////////////

void ARM_ModelParamsSmiled::UpdateIntegratedCovarianceCalib()
{
	if (!itsGlobalCovTimes)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,UpdateIntegratedCovarianceCalib : covariance structure not initialized!! ");

	double voli,volj;
	size_t nbProcs	 = itsNbFwds;
	size_t i,j,size  = itsGlobalCovTimes->size();
	
	for ( size_t index = itsIndexFrom ; index <= itsIndexTo ; index++ )
	{
		if (index==0){
			for (i=0;i<itsNbFwds;i++){
				for (j=0;j<=i;j++){
					(*itsGlobalCovarianceCalib[index])(i,j)=(*itsGlobalCovarianceCalib[index])(j,i)=0.0;
				}
			}
		}else{
			//ne marche que pour beta stepupright
			size_t volIndex = IndexOfFirstHigherEqInVector_DefaultLast((*itsGlobalCovTimes)[index],GetVolCurve(0)->GetAbscisses());
			
			for (i=0;i<itsNbFwds;i++){
				voli = GetVolCurveCalib(i)->GetOrdinates()[volIndex];
				for (j=0;j<=i;j++){
					volj = GetVolCurveCalib(j)->GetOrdinates()[volIndex];
					(*itsGlobalCovarianceCalib[index])(i,j) = (*itsGlobalCovarianceCalib[index])(j,i) = 
						(*itsGlobalCovarianceCalib[index-1])(i,j)
						+ InstantCorrel((*itsGlobalCovTimes)[index],i,j,index-1) * voli * volj * 
						((*itsGlobalCovTimes)[index]-(*itsGlobalCovTimes)[index-1])/K_YEAR_LEN;
				}
			}
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
/// IntegratedCovariance = Integral from s to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::IntegratedCovariance(double s, double t, size_t i, size_t j) const
{
    if(s>=t)
		CC_NS(std,swap)(s,t);

	double sqr = IntegratedCovarianceFromZero(t,i,j)-IntegratedCovarianceFromZero(s,i,j);

	if(fabs(sqr) < K_NEW_DOUBLE_TOL )
		sqr=0;

	return(sqr);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
/// IntegratedCovariance = Integral from s to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
////////////////////////////////////////////////////
void ARM_ModelParamsSmiled::AdviseSwaptionApprox(CalibProxy proxy)
{
	switch (proxy)
	{
		case ATM:
			itsUseModifiedCovForCalib = false;
			break;
		case AtmBlack:
			itsUseModifiedCovForCalib = false;
			break;
		case MomentMatching:
			itsUseModifiedCovForCalib = false;
			break;
		case LocalVolatility:
			itsUseModifiedCovForCalib = true;
			break;
		case LocalVolatilityBlack:
			itsUseModifiedCovForCalib = true;
			break;
		case LocalVolatilityWithRescaling:
			itsUseModifiedCovForCalib = true;
			break;
		case LocalVolatilityBlackWithRescaling:
			itsUseModifiedCovForCalib = true;
			break;
		case EffectiveSkew:
			itsUseModifiedCovForCalib = true;
			break;
		case EffectiveSkewWithRescaling:
			itsUseModifiedCovForCalib = true;
			break;
		case GaussBasketAtm:
			itsUseModifiedCovForCalib = false;
			break;
		case GaussBasketMoneyness:
			itsUseModifiedCovForCalib = false;
			break;
		default:
			itsUseModifiedCovForCalib = false;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
/// IntegratedCovariance = Integral from s to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentVol(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight,const std::vector<double>& coeff, CalibProxy proxy) const
{
	switch (proxy)
	{
		case ATM:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case AtmBlack:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case MomentMatching:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case LocalVolatility:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case LocalVolatilityBlack:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case LocalVolatilityWithRescaling:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case LocalVolatilityBlackWithRescaling:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case EffectiveSkew:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case EffectiveSkewWithRescaling:
			return GetEquivalentVolCalib(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case GaussBasketAtm:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		case GaussBasketMoneyness:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
			break;
		default:
			return GetEquivalentVolModel(evalTime,swapResetTime,i,j,weight,coeff);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetEquivalentCov
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentCov(double evalTime, double swapResetTime, const std::vector<double>& weight1, const std::vector<double>& weight2, const std::vector<double>& coeff, CalibProxy proxy) const
{
	switch (proxy)
	{
		case ATM:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
			break;
		case AtmBlack:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
			break;
		case MomentMatching:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
			break;
		case LocalVolatility:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case LocalVolatilityBlack:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case LocalVolatilityWithRescaling:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case LocalVolatilityBlackWithRescaling:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case EffectiveSkew:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case EffectiveSkewWithRescaling:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovarianceCalib);
			break;
		case GaussBasketAtm:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
			break;
		case GaussBasketMoneyness:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
			break;
		default:
			return GetEquivalentCov(evalTime,swapResetTime,weight1,weight2,coeff,itsGlobalCovariance);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetEquivalentVol
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentCov(double evalTime, double swapResetTime, const std::vector<double>& weight1, const std::vector<double>& weight2, const std::vector<double>& coeff, const ARM_MatrixVector& mat) const
{
	if (!itsGlobalCovTimes)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,GetEquivalentVol : covariance structure not initialized!! ");

	size_t	size	= weight1.size();
	double	cov		= 0;

	size_t index	= IndexOfLastLowerEqInVector(swapResetTime,*itsGlobalCovTimes);
	
	if (swapResetTime == (*itsGlobalCovTimes)[index])
		for (size_t m=0;m<size;m++)
		{
			for (size_t n=0;n<size;n++)
				cov += weight1[m] * weight2[n] * coeff[m] * coeff[n] *(*mat[index])(m,n);
		}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,GetEquivalentVol : should not be there");
	
	return cov;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetEquivalentVol
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentVol(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight,const std::vector<double>& coeff, const ARM_MatrixVector& mat) const
{
	if (!itsGlobalCovTimes)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,GetEquivalentVol : covariance structure not initialized!! ");

	size_t	size	= weight.size();
	double	cov		= 0;

	size_t index	= IndexOfLastLowerEqInVector(swapResetTime,*itsGlobalCovTimes);
	
	if (swapResetTime == (*itsGlobalCovTimes)[index])
		for (size_t m=0;m<size;m++)
		{
			if (fabs(weight[m])> K_NEW_DOUBLE_TOL )
			{
				for (size_t n=0;n<m;n++)
					cov += 2 * weight[m] * weight[n] * coeff[m] * coeff[n] *(*mat[index])(m,n);
				cov += weight[m] * weight[n] * coeff[m] * coeff[n] *(*mat[index])(m,n);
			}
		}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,GetEquivalentVol : should not be there");
	
	#if defined(__GP_STRICT_VALIDATION)
	if ( cov < 0.0 )	    
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsSmiled,GetEquivalentVol : Covariance must be positive");
	#endif

	return sqrt( cov / (swapResetTime-evalTime) * K_YEAR_LEN);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetEquivalentVolCalib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentVolCalib(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight,const std::vector<double>& coeff) const
{
	return GetEquivalentVol(evalTime,swapResetTime,i,j,weight,coeff,itsGlobalCovarianceCalib);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: GetEquivalentVolModel
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::GetEquivalentVolModel(double evalTime, double swapResetTime, size_t i, size_t j, const std::vector<double>& weight,const std::vector<double>& coeff) const
{
	return GetEquivalentVol(evalTime,swapResetTime,i,j,weight,coeff,itsGlobalCovariance);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: IntegratedCovariance
///	Returns: double
///	Action : computes the integrated local covariance
/// IntegratedCovariance = Integral from 0 to t of 
///                         Scalar product(Sigma(u,T1).Sigma(u,T2))*du)
////////////////////////////////////////////////////
double ARM_ModelParamsSmiled::IntegratedCovarianceFromZero(double t, size_t i, size_t j) const
{
	if (t<K_NEW_DOUBLE_TOL)
	{
		return  0.0;
	}

	size_t	k = IndexOfLastLowerEqInVector(t,*itsGlobalCovTimes);

	double covToTk = (*itsGlobalCovariance[k])(i,j);
	double Tk	   = (*itsGlobalCovTimes)[k];
	double covTotal;

	if( fabs(t-Tk)<K_NEW_DOUBLE_TOL )
		covTotal = covToTk;	    
	else
	{
		double voli		= GetVolCurve(i)->Interpolate(t);
		double volj		= GetVolCurve(j)->Interpolate(t);
		double correlij	= InstantCorrel(t,i,j,k);
		double covRes	= voli*volj*correlij*(t-Tk)/K_YEAR_LEN;
		covTotal		= covToTk + covRes;
	}

	return covTotal;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSmiled
///	Routine: toString
///	Returns: string 
///	Action : computes integrated underlying process variance
////////////////////////////////////////////////////
string ARM_ModelParamsSmiled::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "ARM_ModelParamsSmiled\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();

	os << "Volatilities\n";
	os << CC_NS(std,fixed);

	size_t nbproc = itsVol.size();
	if (nbproc>0)
	{
		size_t nbdates = GetVolCurve(0)->GetAbscisses().size();
		os << "Dates" << "\t";
		for (size_t j=0;j<nbproc;j++)
		{
			os << CC_NS(std,setprecision)(0) << (*itsResetDates)[j] << "\t";
			os << "\t";
		}
		os << "\n";
		for (size_t i=0;i<nbdates;i++)
		{
			os << CC_NS(std,setprecision)(0) << GetVolCurve(0)->GetAbscisses()[i] << "\t";
			for (size_t j=0;j<nbproc;j++)
			{
				os << CC_NS(std,setprecision)(6) << GetVolCurve(j)->GetOrdinates()[i];
				os << "\t";
			}
			os << "\n";
		}
	}

	return os.str();
}


CC_END_NAMESPACE()