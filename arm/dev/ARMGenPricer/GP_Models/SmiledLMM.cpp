/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/SmiledLMM.h"
#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SmiledLMM::ARM_SmiledLMM( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, size_t timeStepsNb,size_t gridSize,double stdDevNb,bool skipPDE, bool allowInterpol, ARM_ModelParamsSmiled::CalibProxy calibProxy)
:	ARM_SmiledMM(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy),
	itsWeight(0),
	itsDelta(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SmiledLMM::ARM_SmiledLMM(const ARM_SmiledLMM& rhs)
:	ARM_SmiledMM(rhs),
	itsWeight(rhs.itsWeight),
	itsDelta(rhs.itsDelta)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SmiledLMM::~ARM_SmiledLMM()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SmiledLMM& ARM_SmiledLMM::operator=(const ARM_SmiledLMM& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_SmiledLMM();
		new (this) ARM_SmiledLMM (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledLMM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledLMM::DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if (itsAllowInterpol)
		return ForwardDiscountFactor(curveName,evalTime,evalTime,maturityTime,states);
	else
		return DiscountFactorNoInterpol(curveName,evalTime,maturityTime,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledLMM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledLMM::DiscountFactorNoInterpol( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();
	
	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
// FIXMEFRED: mig.vc8 (30/05/2007 16:34:29):cast
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( maturityTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: MaturityTime after schedule's last date");


	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );
	
	if( DoesResetDateExist( evalTime ) )
	{
		ARM_GP_Vector::iterator iter;
		ARM_GP_Vector::iterator iter1;
		ARM_GP_Vector::iterator iter2;

		size_t ResetIdx = IdxFromValueWithTol( itsResetDates, evalTime, 1. );
		size_t StartIdx = ResetIdx;

		if ( maturityTime >= itsStartDates[itsStartDates.size()-1] )
		{
			StartIdx = itsStartDates.size()-1;
			while (!IsOnSamePath( ResetIdx, StartIdx )) {StartIdx--;}
			
		}
		else
		{
			if ( itsStartDates[StartIdx] > maturityTime )
			{
			}
			else {
				while (itsStartDates[StartIdx] <= maturityTime) { StartIdx++;};
				StartIdx -= 1;
			}
		}

		if ( IsOnSamePath( ResetIdx, StartIdx ) && FastCompute() )
		{
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, StartIdx, modelNb, states);

			for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
				(*iter) *= (*iter1);
		}
		else
		{
			ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, StartIdx, itsResetDates.size(), modelNb, states);
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, itsResetDates.size(), modelNb, states);

			for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
				(*iter) *= (*iter2)/(*iter1);
		}

		if (maturityTime != itsStartDates[StartIdx])
		{
			double theta = (maturityTime - itsStartDates[StartIdx])/360.;
			double delta = itsDelta[StartIdx];
			double mu =		GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
						/	GetZeroCurve()->DiscountPrice(itsStartDates[StartIdx]/K_YEAR_LEN)
						*	( 1 + theta * itsProcess[StartIdx]->GetFwdRate());
			
			ARM_VectorPtr rate = itsProcess[StartIdx]->Rate(evalTime,states,modelNb,StartIdx);

			for (iter = result->begin(),iter2 = rate->begin() ; iter != result->end() ; ++iter,++iter2)
				(*iter) *= mu / ( 1. + theta * (*iter2) );
		}
		
		return static_cast<ARM_VectorPtr>(result);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime not in schedule!!");		
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledLMM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledLMM::ForwardDiscountFactor( const string& curveName, double evalTime, double startTime, double endTime, const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN) ));
	}

	if( endTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: endTime after schedule's last date");

	if( evalTime < itsResetDates[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime before first reset");

	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );

	ARM_GP_Vector::iterator iter;
	ARM_GP_Vector::iterator iter1;
	ARM_GP_Vector::iterator iter2;

	size_t FirstIdx = IndexOfFirstHigherEqInVector_DefaultLast(startTime, itsStartDates); //cas où start dans dernière période!
	if (startTime > itsStartDates[FirstIdx])
		FirstIdx++;

	size_t LastIdx  = IndexOfLastLowerEqInVector(endTime,itsStartDates);
	
	if ( (LastIdx<FirstIdx) || ( IsOnSamePath( FirstIdx, LastIdx ) ) && FastCompute())
	{
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, LastIdx, modelNb, states);

		for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
			(*iter) *= (*iter1);
	}
	else
	{
		ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, LastIdx, itsResetDates.size(), modelNb, states);
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, itsResetDates.size(), modelNb, states);

		for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
			(*iter) *= (*iter2)/(*iter1);
	}

	if ( (LastIdx < FirstIdx) || (startTime < itsStartDates[FirstIdx]))
	{
		double mat				= LastIdx < FirstIdx?endTime:itsStartDates[FirstIdx];
		double theta			= (mat - startTime)/360.;
		size_t idx				= FirstIdx>0?FirstIdx-1:FirstIdx;
		double delta1			= itsDelta[idx];
		double fwd1				= itsProcess[idx]->GetFwdRate();
		double fwd				= (GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)/GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)-1.)/theta;//( 1. - weight2 ) * fwd1 + weight2 * fwd2;

		ARM_VectorPtr rate1	= itsProcess[idx]->Rate(evalTime,states,modelNb,idx,true);
		
		double mu				=		GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)
									/	GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)
									*	( 1 + theta * fwd);

		for (iter = result->begin(),iter1 = rate1->begin(); iter != result->end() ; ++iter,++iter1)
				(*iter) *= mu / ( 1. + theta * fwd/fwd1 * (*iter1) );
	}

	if ( (endTime > itsStartDates[LastIdx]) && (LastIdx >= FirstIdx))
	{
			double theta = (endTime - itsStartDates[LastIdx])/360.;
			double delta = itsDelta[LastIdx];
				
			double mu =		GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)
						/	GetZeroCurve()->DiscountPrice(itsStartDates[LastIdx]/K_YEAR_LEN)
						*	( 1 + theta * itsProcess[LastIdx]->GetFwdRate());
					
			ARM_VectorPtr rate = itsProcess[LastIdx]->Rate(evalTime,states,modelNb,LastIdx,true);

			for (iter = result->begin(),iter2 = rate->begin() ; iter != result->end() ; ++iter,++iter2)
				(*iter) *= mu / ( 1. + theta * (*iter2) );
	}

	return static_cast<ARM_VectorPtr>(result);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledLMM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledLMM::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t IdxFrom, size_t IdxTo , size_t modelNb, const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr rate;
	ARM_GP_Vector* result	= new ARM_GP_Vector( states->size(),1.0 );
		
	for (size_t k = IdxFrom; k < IdxTo ; k++ )
	{
		if ( itsWeight(IdxFrom,k) == 1. )
		{
			rate = itsProcess[k]->Rate(evalTime,states,modelNb,k);

			ARM_GP_Vector::iterator iter2=rate->begin();

			for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter, ++iter2)
				(*iter) /= ( 1. + itsDelta[k] * (*iter2));
		}
	}
	return static_cast<ARM_VectorPtr>(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledLMM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double currTime  = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime  = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb = FactorCount();
	size_t statesNb  = states->size();
	size_t modelNb	 = GetModelNb();

	size_t j,k;
	size_t fwdsNb  = itsResetDates.size();
	size_t iLast   = fwdsNb-1;

	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetDates[iFirst]<nextTime)
	{
		++iFirst;
	}
	if (itsAllowInterpol && iFirst>0) iFirst--;

	//parait débile mais indispensable pour les échéanciers non triviaux!!!
	size_t iFirstMartingale = iFirst;
	while( itsEndDates[iFirstMartingale] < itsEndDates[iLast]) iFirstMartingale++;
	
	const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();

	ARM_VectorPtrVector dxRate(iLast-iFirst+1);
	ARM_VectorPtrVector rate(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
	{
		rate[k-iFirst] = itsProcess[k]->Rate(currTime,states,modelNb,k,true);
		dxRate[k-iFirst] = itsProcess[k]->DxRate(currTime,states,modelNb,k,true);
	}

	int i;
	double inc,modStd;
	double lib,dlib;

	ARM_GP_Vector* eigen = itsEigenValues[timeIndex];
	ARM_GP_Matrix E(fwdsNb,factorsNb,0.);

	vector<ARM_GP_Vector::iterator> rateIter(iLast-iFirst+1);
	vector<ARM_GP_Vector::iterator> dxRateIter(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
	{
		rateIter[k-iFirst] = rate[k-iFirst]->begin();
		dxRateIter[k-iFirst] = dxRate[k-iFirst]->begin();
	}
	
	for(k = 0 ; k < statesNb ; ++k )
	{
		for ( i=iFirstMartingale ; i<fwdsNb ; i++ )
			for (j = 0 ; j < factorsNb ; j++ )
				E(i,j) = 0.;
		
		for (i = iLast; i >= iFirst; i--)
		{
			inc			= 0;
			size_t prec	= i+1;

			if (itsEndDates[i] == itsEndDates[iLast])
				prec = i;
			else
				while (!IsOnSamePath(i,prec)) prec++;

			lib  = itsDelta[i] * (*(rateIter[i-iFirst]++));
			dlib = itsDelta[i] * (*(dxRateIter[i-iFirst]++));

			j = 0;
			ARM_GP_Vector::iterator iter = eigen->begin();
			for (; iter != eigen->end() ; ++iter,++j )
			{
				modStd	=	(*modelLocalVar[timeIndex])(i,j);
				inc		+=	modStd * ( states->GetNumMethodState(k,j+modelNb) + (*iter) * E(prec,j));
				E(i,j)	=	E(prec,j) - modStd * dlib / ( 1 + lib);
			}
			states->SetModelState(k,i+modelNb,states->GetModelState(k,i+modelNb)+inc);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: ComputeDZCDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

double ARM_SmiledLMM::DZCDRate( double maturity, size_t i) const
{
	if ( maturity < itsStartDates[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM,DZCDLib : maturity < first startDate" );

	size_t StartIdx = 0;

	if ( maturity >= itsStartDates[itsStartDates.size()-1] )
		StartIdx = itsStartDates.size()-1;
	else
	{
		if ( itsStartDates[StartIdx] > maturity)
		{
		}
		else 
		{
			while (itsStartDates[StartIdx] <= maturity) { StartIdx++;};
			StartIdx -= 1;
		}
	}

	double res = 0.;

	if ( IsOnSamePath(0,StartIdx) && IsOnSamePath(0,i) )
	{
		if ( i < StartIdx )
			res =	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsDelta[i] / (1 + itsDelta[i] * itsProcess[i]->GetFwdRate() ) ;
		else
			if ( StartIdx == i )
			{
				double theta = (maturity - itsStartDates[i])/360.;
				res =	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
						* theta / (1 + theta * itsProcess[i]->GetFwdRate() );
			}
	}
	else
	{
		if ( IsOnSamePath(0,i))
			res = - GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsDelta[i] / (1 + itsDelta[i] * itsProcess[i]->GetFwdRate() ) ;
		else 
			if ( IsOnSamePath(StartIdx,i))
				res = GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsDelta[i] / (1 + itsDelta[i] * itsProcess[i]->GetFwdRate() ) ;

		if ( StartIdx == i )
			{
				double theta = (maturity - itsStartDates[i])/360.;
				res +=	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
						* theta / (1 + theta * itsProcess[i]->GetFwdRate() );
			}
	}

	return res;
}



////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter + reset/start/end dates to DF map
////////////////////////////////////////////////////
	
void ARM_SmiledLMM::setNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{

	if (numericalModelFitter)
	{
		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		size_t size = densityVector.size();

		ARM_GP_Vector resetDates(size); 
		ARM_GP_Vector startDates(size); 
		ARM_GP_Vector endDates(size); 
		ARM_GP_Vector delta(size); 

		for (size_t k = 0 ; k < size ; k++)
		{
			resetDates[k] = densityVector[k]->getResetDate();
			startDates[k] = densityVector[k]->getStartDate();
			endDates[k] = densityVector[k]->getEndDate();
			if (densityVector[k]->getInterestTerms().size()!=1)
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::setNumericalModelFitter : only LIBOR allowed!" );
			delta[k] = densityVector[k]->getInterestTerms()[0];
		}
		
		double asof = GetAsOfDate().GetJulian();

		setResetTimes( resetDates - asof);
		setStartTimes( startDates - asof);
		setEndTimes  ( endDates   - asof);
		setDelta	 ( delta);

		computeWeights();
		itsCalibrated = false;
		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);
	}
	itsNumericalModelFitter = numericalModelFitter;

}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledLMM
///	Routine: computeWeights
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
	
void ARM_SmiledLMM::computeWeights() 
{
	size_t size		= itsResetDates.size();
	double terminal	= itsEndDates[size-1];
	double term;

	itsWeight.resize( size, size);

	for (size_t j=0;j<size;j++)
	{
		size_t k = 0;
		term = itsEndDates[j];
		while (term<terminal)
		{
			while (itsStartDates[k]<term)
			{
				itsWeight(j,k)=0.;
				k++;
			}
			itsWeight(j,k)=1.;
			term = itsEndDates[k];
			k++;
		}
		itsWeight(j,j)=1.;
	}


}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

