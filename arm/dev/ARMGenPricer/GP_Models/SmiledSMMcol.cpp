/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/datestrip.h"

#include "gpmodels/SmiledSMMcol.h"
#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SmiledSMMcol::ARM_SmiledSMMcol( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, size_t timeStepsNb,size_t gridSize,double stdDevNb,bool skipPDE, bool allowInterpol, ARM_ModelParamsSmiled::CalibProxy calibProxy,bool cache)
:	ARM_SmiledMM(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy,cache),
	itsTheta(0),
	itsPreviousRate(0),
	itsNextRate(0),
	itsLastRate(0),
	itsPreviousRate2(0),
	itsNextRate2(0),
	itsLastRate2(0),
	itsExtraLibor(0),
	itsWeight(0),
	itsWeightZC(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SmiledSMMcol::ARM_SmiledSMMcol(const ARM_SmiledSMMcol& rhs)
:	ARM_SmiledMM(rhs),
	itsTheta(rhs.itsTheta),
	itsPreviousRate(rhs.itsPreviousRate),
	itsNextRate(rhs.itsNextRate),
	itsLastRate(rhs.itsLastRate),
	itsPreviousRate2(rhs.itsPreviousRate2),
	itsNextRate2(rhs.itsNextRate2),
	itsLastRate2(rhs.itsLastRate2),
	itsExtraLibor(rhs.itsExtraLibor),
	itsWeight(rhs.itsWeight),
	itsWeightZC(rhs.itsWeightZC)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SmiledSMMcol::~ARM_SmiledSMMcol()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_SmiledSMMcol& ARM_SmiledSMMcol::operator=(const ARM_SmiledSMMcol& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_SmiledSMMcol();
		new (this) ARM_SmiledSMMcol (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledSMMcol::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	if (itsCache)
	{
		MCModelStatesFromToNextTime_std(states,timeIndex);
		MCModelStatesFromToNextTime_evalOnly(states,timeIndex+1);
	}
	else
	{
		MCModelStatesFromToNextTime_std(states,timeIndex);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: MCModelStatesFromToNextTime_fast
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledSMMcol::MCModelStatesFromToNextTime_std(ARM_PricingStatesPtr& states,int timeIndex) const
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

	size_t iFirstMartingale = iFirst;
	while( itsNextRate[iFirstMartingale] > -1) iFirstMartingale++;
	
	const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();

	ARM_VectorPtrVector dxRate(iLast-iFirst+1);
	ARM_VectorPtrVector rate(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
	{
		rate[k-iFirst] = itsProcess[k]->Rate(currTime,states,modelNb,k,true);
		dxRate[k-iFirst] = itsProcess[k]->DxRate(currTime,states,modelNb,k,true);
	}

	int i;
	double next,std,sensi0,coeff1,coeff2,coeff3,theta=1.,nextRate;

	ARM_GP_Vector* eigen = itsEigenValues[timeIndex];
	ARM_GP_Matrix* covar = modelLocalVar[timeIndex];
	ARM_GP_Matrix E(fwdsNb,factorsNb,0.);
	ARM_GP_Vector G(fwdsNb,0.);
	ARM_GP_Vector D(fwdsNb,0.);
	ARM_GP_Vector S(fwdsNb,0.);

	vector<ARM_GP_Vector::iterator> rateIter(iLast-iFirst+1);
	vector<ARM_GP_Vector::iterator> dxRateIter(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
	{
		rateIter[k-iFirst] = rate[k-iFirst]->begin();
		dxRateIter[k-iFirst] = dxRate[k-iFirst]->begin();
	}
	
	for(k = 0 ; k < statesNb ; ++k )
	{
		for ( i = iFirstMartingale ; i < fwdsNb ; i++ )
		{
			G(i) = itsTheta[i];
			for (j = 0 ; j < factorsNb ; j++ )
				E(i,j) = 0.;
		}

		for (i = iLast; i >= iFirst; i--)
		{
			next		= states->GetModelState(k,i+modelNb);
			j			= 0;
			S[i]		= (*(rateIter[i-iFirst]++));
			D[i]		= (*(dxRateIter[i-iFirst]++));
			nextRate	= i>itsNextRate[i]?i:itsNextRate[i];

			theta = (itsPreviousRate[i]>-1?itsTheta[itsPreviousRate[i]]:1);
		
			if (itsExtraLibor[i]>-1)
				coeff3 = theta * itsTheta[itsExtraLibor[i]] ;
	
			coeff1		= ( 1 + theta *  S[i]) * G(i);
			coeff2		= theta * G(i);
			
			if (itsLastRate[i]>-1)
			{
				sensi0		= coeff1;
				if ( isDiag(itsPreviousRate[i]) )
				{
					sensi0		= sensi0 + theta;
					if (itsExtraLibor[i]>-1)
						sensi0  = sensi0 + coeff3 * S[itsExtraLibor[i]];
				}
			}
			else				
				sensi0	= theta * (1. + itsTheta[i] * S[i] ) ;
				
			
			for (ARM_GP_Vector::iterator iter = eigen->begin(); iter != eigen->end() ; ++iter,++j )
			{
				std		=	(*covar)(i,j);
				next	+=	(*covar)(i,j) * ( states->GetNumMethodState(k,j+modelNb) - (*iter) * E(nextRate,j));

				if (itsLastRate[i]>-1)
				{
					E(i,j)	=	coeff1/sensi0 * E(nextRate,j) +	coeff2/sensi0 * std * D[i];
					if (itsExtraLibor[i]>-1 && isDiag(itsPreviousRate[i])) 
						E(i,j)	=	E(i,j) +  coeff3/sensi0 * (*covar)(itsExtraLibor[i],j) * D[itsExtraLibor[i]];
				}
				else
				{
					E(i,j) = theta * itsTheta[i]/ sensi0 * std * D[i] ;
				}
			}

			if (itsPreviousRate[i]>-1)
				G(itsPreviousRate[i]) = sensi0;
			
			states->SetModelState(k,i+modelNb,next);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: MCModelStatesFromToNextTime_std
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledSMMcol::MCModelStatesFromToNextTime_evalOnly(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double currTime  = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime  = currTime;
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

	size_t iFirstMartingale = iFirst;
	while( itsNextRate[iFirstMartingale] > -1) iFirstMartingale++;
	
	const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();

	ARM_VectorPtrVector rate(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
		rate[k-iFirst] = itsProcess[k]->Rate(currTime,states,modelNb,k,true);

	int i,nextRate;
	double next,sensi0,coeff1,coeff2,coeff3,theta=1.,zc;

	ARM_GP_Vector G(fwdsNb,0.);
	ARM_GP_Vector P(fwdsNb,0.);
	ARM_GP_Vector S(fwdsNb,0.);

	vector<ARM_GP_Vector::iterator> rateIter(iLast-iFirst+1);
	for (k = iFirst;k<=iLast;k++)
		rateIter[k-iFirst] = rate[k-iFirst]->begin();
	
	for(k = 0 ; k < statesNb ; ++k )
	{
		for ( i = iFirstMartingale ; i < fwdsNb ; i++ )
			G(i) = itsTheta[i];
		
		for (i = iLast; i >= iFirst; i--)
		{
			next		= states->GetModelState(k,i+modelNb);
			j			= 0;
			S[i]		= (*(rateIter[i-iFirst]++));
			nextRate	= i>itsNextRate[i]?i:itsNextRate[i];

			theta = (itsPreviousRate[i]>-1?itsTheta[itsPreviousRate[i]]:1);
		
			if (itsExtraLibor[i]>-1)
				coeff3 = theta * itsTheta[itsExtraLibor[i]] ;
	
			coeff1		= ( 1 + theta *  S[i]) * G(i);
			coeff2		= theta * G(i);
			
			if (itsLastRate[i]>-1)
			{
				sensi0		= coeff1;
				if ( isDiag(itsPreviousRate[i]) )
				{
					sensi0		= sensi0 + theta;
					if (itsExtraLibor[i]>-1)
						sensi0  = sensi0 + coeff3 * S[itsExtraLibor[i]];
				}
				else
				{
					if (itsNextRate2[i]>-1)
					{	
						double teta = itsTheta[itsPreviousRate[itsNextRate2[i]]];
						sensi0 = sensi0 + (theta - teta ) * P(itsNextRate2[i]);
					}
					else
					{
						double teta = itsTheta[itsNextRate2[itsPreviousRate[i]]];
						sensi0 = sensi0 + (theta - teta);
					}
				}
			}
			else
			{
				sensi0	= theta * (1. + itsTheta[i] * S[i] ) ;
			}

			if (itsLastRate[i]>-1)
			{
				if ( isDiag(i) )
					zc = itsExtraLibor[i]>-1?P(itsExtraLibor[i])+S(i)*G(i):1+S(i)*G(i);
				else
					zc = itsNextRate2[i]>-1?P(itsNextRate2[i])+S(i)*G(i):1+S(i)*G(i);
			}
			else
				zc = 1+S(i)*G(i);
			P(i) = zc;
				
			if (itsPreviousRate[i]>-1)
				G(itsPreviousRate[i]) = sensi0;
			
			states->SetModelState(k,i+modelNb+fwdsNb,zc);
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledSMMcol::DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
// FIXMEFRED: mig.vc8 (30/05/2007 16:32:52):cast
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( maturityTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: MaturityTime after schedule's last date");


	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );
	
	if( DoesResetDateExist( evalTime ) )
	{
		size_t StartIdx = IdxFromValueWithTol( itsResetDates, evalTime, 1. );
		size_t EndIdx   = IndexOfLastLowerEqInVector_DefaultFirst( maturityTime, itsStartDates);
		if (EndIdx<StartIdx) // between reset and start
			EndIdx++;

		ARM_VectorPtr fwdEnd = ForwardDiscountFactorFromIdx(curveName, evalTime, EndIdx, states);
		ARM_VectorPtr fwdStart  = ForwardDiscountFactorFromIdx(curveName, evalTime, StartIdx, states);

		ARM_GP_Vector::iterator iResult = result->begin();
		ARM_GP_Vector::iterator iFwdStart = fwdStart->begin();
		ARM_GP_Vector::iterator iFwdEnd = fwdEnd->begin();

		for (; iResult != result->end() ; ++iResult, ++iFwdEnd, ++iFwdStart)
			(*iResult) = (*iFwdEnd)/(*iFwdStart);

		if (maturityTime != itsStartDates[EndIdx])
		{
			double endDate = (EndIdx+1<itsResetDates.size()?itsStartDates[EndIdx+1]:itsEndDates[EndIdx]);
			double theta = (maturityTime - itsStartDates[EndIdx])/360.;
			double delta = (endDate - itsStartDates[EndIdx])/360.;

			double fwd = (	GetZeroCurve()->DiscountPrice(itsStartDates[EndIdx]/K_YEAR_LEN)/
							GetZeroCurve()->DiscountPrice(endDate/K_YEAR_LEN)	
							- 1 )/delta;
							
			double mu =		GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
						/	GetZeroCurve()->DiscountPrice(itsStartDates[EndIdx]/K_YEAR_LEN)
						*	( 1 + theta * fwd);

			ARM_VectorPtr fwdEndPlus = ForwardDiscountFactorFromIdx(curveName, evalTime, EndIdx+1, states);

			ARM_GP_Vector::iterator iResult = result->begin();
			ARM_GP_Vector::iterator iFwdEnd = fwdEnd->begin();
			ARM_GP_Vector::iterator iFwdEndPlus = fwdEndPlus->begin();

			for (; iResult != result->end() ; ++iResult, ++iFwdEnd, ++iFwdEndPlus)
				(*iResult) *= mu / ( 1. + theta/delta * ( (*iFwdEnd) / (*iFwdEndPlus) - 1 ));
		}
		return static_cast<ARM_VectorPtr>(result);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::DiscountFactor: evalTime not in schedule!!");		
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledSMMcol::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t Idx, const ARM_PricingStatesPtr& states) const
{
	if (itsCache)
		return ForwardDiscountFactorFromIdx_cache(curveName,evalTime,Idx,states);
	else
		return ForwardDiscountFactorFromIdx_std(curveName,evalTime,Idx,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledSMMcol::ForwardDiscountFactorFromIdx_cache( const string& curveName, double evalTime, size_t Idx, const ARM_PricingStatesPtr& states) const
{
	size_t fwdsNb  = itsResetDates.size();
	size_t modelNb = GetModelNb();
	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.);
	
	if (Idx<fwdsNb)
	{
		ARM_GP_Vector::iterator iResult;
		iResult = result->begin();
		size_t k = 0;
		for (; iResult != result->end() ; ++iResult, ++k)
			(*iResult) = states->GetModelState(k,Idx+modelNb+fwdsNb);
	}

	return ARM_GP_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledSMMcol::ForwardDiscountFactorFromIdx_std( const string& curveName, double evalTime, size_t Idx, const ARM_PricingStatesPtr& states) const
{
	size_t fwdsNb  = itsResetDates.size();
	size_t modelNb = GetModelNb();
	ARM_GP_VectorPtr init = ARM_GP_VectorPtr(new ARM_GP_Vector( states->size(),1.));
	ARM_GP_VectorPtr initZC = ARM_GP_VectorPtr(new ARM_GP_Vector( states->size(),1.));
	
	if (Idx==fwdsNb)
		return init;


	ARM_GP_Vector::iterator iRate;
	ARM_GP_Vector::iterator iResult;
	ARM_GP_Vector::iterator iInit;
	ARM_GP_Vector::iterator iInitZC;
		
	ARM_GP_VectorPtr rate;
	ARM_GP_Vector* result;
	
	if (itsLastRate[Idx]>-1)
	{
		result = new ARM_GP_Vector( states->size(),itsTheta[itsLastRate[Idx]]);

		if (itsExtraLibor[Idx]>-1)
		{
			rate = itsProcess[itsExtraLibor[Idx]]->Rate(evalTime,states,modelNb,itsExtraLibor[Idx]);
			iRate = rate->begin();
			iInit = init->begin();
			iResult = result->begin();
			iInitZC = initZC->begin();
			for (; iInit != init->end() ; ++iInit, ++iRate,++iResult,++iInitZC)
			{
				(*iInit)	=  1. + itsTheta[itsExtraLibor[Idx]] * (*iRate);
				(*iResult)	*= (*iInit);
				(*iInitZC)	=  (*iInit);
			}

		}

		size_t nextCMS = itsLastRate2[Idx];
		for (int k = itsLastRate[Idx]; k > Idx; k-- )
		{
			if (itsLastRate[k] == itsLastRate[Idx])
			{
				rate = itsProcess[k]->Rate(evalTime,states,modelNb,k);

				if (nextCMS==k)
				{
					iRate = rate->begin();
					iInitZC = initZC->begin();
					iResult = result->begin();
					for (; iResult != result->end() ; ++iResult, ++iRate,++iInitZC)
						(*iInitZC) = (*iInitZC) + (*iRate) * (*iResult);
					nextCMS = itsPreviousRate2[nextCMS];
				}

				if (itsNextRate2[itsPreviousRate[k]]==itsExtraLibor[itsPreviousRate[k]])
				{
					iRate = rate->begin();
					iInit = init->begin();
					iResult = result->begin();
					for (; iResult != result->end() ; ++iResult, ++iRate,++iInit)
						(*iResult) = ( 1. + itsTheta[itsPreviousRate[k]] * (*iRate)) * (*iResult) + itsTheta[itsPreviousRate[k]] * (*iInit);
				}
				else
				{
					iRate = rate->begin();
					iResult = result->begin();
					for (; iResult != result->end() ; ++iResult, ++iRate)
						(*iResult) = ( 1. + itsTheta[itsPreviousRate[k]] * (*iRate)) * (*iResult);
				}
			}
		}

		rate = itsProcess[Idx]->Rate(evalTime,states,modelNb,Idx);

		iRate = rate->begin();
		iInitZC = initZC->begin();
		iResult = result->begin();
		for (; iResult != result->end() ; ++iResult, ++iRate,++iInitZC)
			(*iResult) =  (*iInitZC) + (*iRate) * (*iResult);
	}
	else
	{
		rate = itsProcess[Idx]->Rate(evalTime,states,modelNb,Idx);
		result = new ARM_GP_Vector( states->size(),1.);

		iRate = rate->begin();
		iResult = result->begin();
		for (; iResult != result->end() ; ++iResult, ++iRate)
			(*iResult) +=  itsTheta[Idx] * (*iRate);
	}

	return ARM_GP_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter + reset/start/end dates to DF map
////////////////////////////////////////////////////
	
void ARM_SmiledSMMcol::setNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{
	if (numericalModelFitter)
	{
		checkNumericalModelFitter(numericalModelFitter);

		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		size_t size = densityVector.size();

		ARM_GP_Vector resetDates(size); 
		ARM_GP_Vector startDates(size); 
		ARM_GP_Vector endDates(size); 
		ARM_GP_Vector theta(size); 

		for (size_t k = 0 ; k < size ; k++)
		{
			resetDates[k]	= densityVector[k]->getResetDate();
			startDates[k]	= densityVector[k]->getStartDate();
			endDates[k]		= densityVector[k]->getEndDate();
			theta[k]		= densityVector[k]->getInterestTerms()[0];
		}
		
		double asof = GetAsOfDate().GetJulian();

		setResetTimes( resetDates - asof);

		setStartTimes( startDates - asof);
		setEndTimes  ( endDates   - asof);
		setTheta	 ( theta);

		itsCalibrated = false;
		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);
	}
	itsNumericalModelFitter = numericalModelFitter;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter + reset/start/end dates to DF map
////////////////////////////////////////////////////
	
void ARM_SmiledSMMcol::checkNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{
	if (numericalModelFitter)
	{
		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		size_t size = densityVector.size();

		if (size>0)
		{
			if (densityVector[0]->getInterestTerms().size()>1)
			{
				//deduce nb of family from nb of Libor at end
				size_t nbLibor = 0;
				while (densityVector[size-1-nbLibor]->getInterestTerms().size()==1) nbLibor++;
				size_t nbFamily = (nbLibor+1)/2;

				int rateFreq = densityVector[0]->getFrequency();
				int calibSchedFreq = CC_Min(CC_Round( 365. / (numericalModelFitter->getCalibSchedule()->GetFlowEndDates()->Elt(0) - numericalModelFitter->getCalibSchedule()->GetFlowStartDates()->Elt(0)) ),12) ;
				
				if (nbFamily != calibSchedFreq/rateFreq)
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: pb with frequency!");
		
				//get enddate for every diagonal structure
				ARM_GP_Vector EndDates(nbFamily);
				for (size_t k=0;k<nbFamily;k++)
				{
					EndDates[k] = densityVector[size-nbLibor+k]->getEndDate();
					if (k>0)
						if (EndDates[k]==EndDates[k-1])
							ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: Date structure not allowed!");
				}

				//check that each rate has coherent enddate
				for (k=0;k<nbFamily-1;k++)
					if (densityVector[size-nbFamily+k]->getEndDate()!=EndDates[nbFamily-1])
						ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: Libor at end should end at last date!");

				//look for first real CMS
				if (densityVector[0]->getEndDate()==EndDates[nbFamily-1-size%nbFamily])
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: first rate is diagonal!");

				k=nbFamily;
				while (densityVector[size-k]->getEndDate()==EndDates[nbFamily-1-k%nbFamily]) 
					k++;
				//size-k est l'indice du premier non terminal

				size_t nbCMS = size-k+1;
				size_t tenor = densityVector[size-k+1]->getInterestTerms().size();
				
				for (;k<size;k++)
					if (densityVector[size-k]->getInterestTerms().size()!=tenor)
						ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: Date structure not allowed!");

				//Dates structures for PreviousRate, NextRate and ExtraLibor
				size_t nbSwap = size-(nbFamily-1);
				size_t nbBeg = nbSwap%nbFamily;
				ARM_IntVector previous(size);
				ARM_IntVector next(size);
				ARM_IntVector last(size);
				ARM_IntVector previous2(size);
				ARM_IntVector next2(size);
				ARM_IntVector last2(size);
				ARM_IntVector extra(size);

				for (k=0;k<size;k++)
				{
					previous[k]	= (k<nbFamily?-1:k-nbFamily);
					next[k]		= (k<size-nbFamily?k+nbFamily:-1);
					size_t index = (k+nbFamily-nbBeg) % nbFamily;
					if (k<nbSwap)
					{
						last[k]	 = size-nbLibor+index;
						extra[k] = (densityVector[k]->getEndDate()== densityVector[size-1]->getEndDate()?-1:size-(nbFamily-1)+index);
						next2[k] = k+tenor*nbFamily>size-1?extra[k]:k+tenor*nbFamily;
						previous2[k]  = k<tenor*nbFamily?-1:k-tenor*nbFamily;
					}
					else
					{
						last[k]	 = -1;
						extra[k] = -1;
						next2[k] = -1;
						previous2[k]  = -1;
						last2[k] = -1;
					}
				}
				for (k=0;k<nbCMS;k++)
					extra[nbCMS-1-k]=extra[last[nbCMS-1-k]];

				for (k=0;k<nbSwap;k++)
				{
					last2[k]=k;
					while (next2[last2[k]]!=extra[k]) 
						last2[k]=next2[last2[k]];
					if (last2[k]==k)
						last2[k]=next2[k];
				}

				setLinks(previous,next,last,previous2,next2,last2,extra);
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: first rate is Libor!");
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: numericalModelFitter has size 0!");
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledSMMcol::checkNumericalModelFitter: numericalModelFitter is null");
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: computeDWeights
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
void ARM_SmiledSMMcol::computeDWeights()
{
	int i,nextRate;
	double sensi0,zc,theta=1.;

	size_t k;
	size_t fwdsNb  = itsResetDates.size();
	int iLast   = fwdsNb-1;
	int iFirst  = 0;

	size_t iFirstMartingale = iFirst;
	while( itsNextRate[iFirstMartingale] > -1) iFirstMartingale++;
	
	ARM_GP_Vector G(fwdsNb,0.);
	ARM_GP_Vector P(fwdsNb,0.);
	ARM_GP_Vector S(fwdsNb,0.);

	itsWeight.resize(fwdsNb+1,fwdsNb);
	itsWeightZC.resize(fwdsNb+1,fwdsNb);
	
	for ( i = iFirstMartingale ; i < fwdsNb ; i++ )
		G(i) = itsTheta[i];
	
	for (i = iLast; i >= iFirst; i--)
	{
		S[i]		= itsProcess[i]->GetFwdRate();
		nextRate	= i>itsNextRate[i]?i:itsNextRate[i];
		theta = (itsPreviousRate[i]>-1?itsTheta[itsPreviousRate[i]]:1);
			
		if (isExtraLibor(i))
		{
			sensi0	= theta * (1. + itsTheta[i] * S[i] );
			itsWeight(i,i) += theta * itsTheta[i];
		}
		else
		{
			sensi0		= ( 1 + theta *  S[i]) * G(i);

			itsWeight(i,i) += theta * G(i);
			for ( k = i+1; k < fwdsNb; k++)
				itsWeight(i,k) += ( 1 + theta *  S[i]) * itsWeight(nextRate,k);
			
			if ( isDiag(itsPreviousRate[i]) )
			{
				sensi0		+= theta;
				if (itsExtraLibor[i]>-1)
				{
					sensi0  += theta * itsTheta[itsExtraLibor[i]] * S[itsExtraLibor[i]];
					itsWeight(i,itsExtraLibor[i]) += theta * itsTheta[itsExtraLibor[i]];
				}
			}
			else
			{
				if (itsNextRate2[i]>-1)
				{
					double teta = itsTheta[itsPreviousRate[itsNextRate2[i]]];
					sensi0 = sensi0 + (theta - teta ) * P(itsNextRate2[i]);

					for ( k = i+1; k < fwdsNb; k++)
						itsWeight(i,k) += (theta - teta ) * itsWeightZC(itsNextRate2[i],k);
				}
				else
				{
					double teta = itsTheta[itsNextRate2[itsPreviousRate[i]]];
					sensi0 = sensi0 + (theta - teta );
				}
			}
		}

		nextRate = itsNextRate[i]>-1?itsNextRate[i]:fwdsNb;
		if ( isDiag(i) && (itsExtraLibor[i]>-1))
		{
			zc = P(itsExtraLibor[i])+S(i)*G(i);
			itsWeightZC(i,i) = G(i);
			for ( k = i+1; k < fwdsNb; k++)
				itsWeightZC(i,k) = itsWeightZC(itsExtraLibor[i],k)+S(i)*itsWeight(nextRate,k);
		}
		else if (itsNextRate2[i]>-1)
		{
			zc = P(itsNextRate2[i])+S(i)*G(i);
			itsWeightZC(i,i) = G(i);
			for ( k = i+1; k < fwdsNb; k++)
				itsWeightZC(i,k) = itsWeightZC(itsNextRate2[i],k)+S(i)*itsWeight(nextRate,k);
		}
		else
		{
			zc = 1+S(i)*G(i);
			itsWeightZC(i,i) = G(i);
			for ( k = i+1; k < fwdsNb; k++)
				itsWeightZC(i,k) = S(i)*itsWeight(nextRate,k);
		}
		P(i) = zc;
				
		if (itsPreviousRate[i]>-1)
			G(itsPreviousRate[i]) = sensi0;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledSMMcol
///	Routine: ComputeDZCDRate
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

double ARM_SmiledSMMcol::DZCDRate( double maturityTime, size_t i) const
{
	if( maturityTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: MaturityTime after schedule's last date");

	double result;

	size_t EndIdx   = IndexOfLastLowerEqInVector_DefaultFirst( maturityTime, itsStartDates);
	
	if (maturityTime != itsStartDates[EndIdx])
	{
		double endDate = (EndIdx+1<itsResetDates.size()?itsStartDates[EndIdx+1]:itsEndDates[EndIdx]);
		double theta = (maturityTime - itsStartDates[EndIdx])/360.;
		double delta = (endDate - itsStartDates[EndIdx])/360.;

		double end	   = itsWeightZC(EndIdx,i);
		double endplus = itsWeightZC(EndIdx+1,i);

		result = ( 1 - theta/delta) * end + theta * delta * endplus;
	}
	else
		result	= itsWeightZC(EndIdx,i);

	return result*GetZeroCurve()->DiscountPrice(getTerminalTime()/K_YEAR_LEN);
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledSMMcol
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SmiledSMMcol::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
  	os << ARM_SmiledMM::toString(indent);
	os << "\n\n";

	os << indent << "Rate Structure\n";
	os << "Rate" << "\t"<< "Prev" << "\t"<< "Next" << "\t"<< "Last" << "\t";
	os << "Prev2" << "\t"<< "Next2" << "\t"<< "Last2" << "\t"<< "Extra"<< "\t" <<"Theta" << "\n";
	for (size_t k = 0;k<itsPreviousRate.size();k++)
	{
		os << CC_NS(std,setprecision)(0) << k << "\t";
		os << itsPreviousRate[k] << "\t"<< itsNextRate[k]  << "\t"<< itsLastRate[k]  << "\t";
		os << itsPreviousRate2[k] << "\t"<< itsNextRate2[k]  << "\t"<< itsLastRate2[k]  << "\t"<< itsExtraLibor[k] << "\t"<< itsTheta[k] << "\n";
	}

	os << "\n\n";
	os << indent << "Level weights\n";
	os << "Level> " << "\t";
		
	for (size_t j=0;j<itsWeight.rows();j++)
		os << CC_NS(std,setprecision)(0) << j << "\t";
	os << "\n";
	for (size_t i=0;i<itsWeight.cols();i++)
	{
		os << CC_NS(std,setprecision)(0) << i<< "\t";
		for (size_t j=0;j<itsWeight.rows();j++)
			os << CC_NS(std,setprecision)(4) << itsWeight(j,i)<< "\t";
		os << "\n";
	}

	os << "\n\n";
	os << indent << "FwdZC weights\n";
	os << "FwdZC> " << "\t";
		
	for (j=0;j<itsWeightZC.rows();j++)
		os << CC_NS(std,setprecision)(0) << j << "\t";
	os << "\n";
	for (i=0;i<itsWeightZC.cols();i++)
	{
		os << CC_NS(std,setprecision)(0) << i<< "\t";
		for (j=0;j<itsWeightZC.rows();j++)
			os << CC_NS(std,setprecision)(4) << itsWeightZC(j,i)<< "\t";
		os << "\n";
	}

	return os.str();	
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
