/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file pricingmodelinflation.cpp
 *  \brief interest rate version of a pricing model!
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 */


#include "gpbase/datestrip.h"
#include "gpinfra/pricingmodelinflation.h"
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinflation/infidx.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelInflation::ARM_PricingModelInflation( const ARM_InfCurvPtr& infc, 
						const ARM_ModelParams* params )
:	ARM_PricingModel(ARM_ZeroCurvePtr(NULL),params), itsInfCurve(infc), itsIRModel(NULL), itsCorrelWithBonds(0)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_PricingModelInflation::ARM_PricingModelInflation(const ARM_PricingModelInflation& rhs)
: ARM_PricingModel(rhs), itsInfCurve(rhs.itsInfCurve), itsIRModel(rhs.itsIRModel), itsCorrelWithBonds(rhs.itsCorrelWithBonds)
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_PricingModelInflation::~ARM_PricingModelInflation()
{}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_PricingModelInflation& ARM_PricingModelInflation::operator=(const ARM_PricingModelInflation& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModel::operator=(rhs);
        // Copy class attributes if any
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: DiscountFactor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelInflation::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
        const ARM_PricingStatesPtr& states) const
{

#ifdef __GP_STRICT_VALIDATION

	if( evalTime > maturityTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a discountFactor in the past.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	if( !itsIRModel )
	{
		CC_Ostringstream os;
		os << "An inflation model cannot work without an IRModel\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#endif

	return itsIRModel->DiscountFactor( curveName, evalTime,maturityTime,states);

}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: DefaultAnnuity
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_PricingModelInflation::DefaultAnnuity(const string& curveName, double evalTime,const ARM_DateStripPtr& PayDateStrip, 
													const ARM_PricingStatesPtr& states) const
{
	ARM_GP_Vector * PayDates = PayDateStrip->GetFlowEndDates();
	ARM_GP_Vector * InterestTerms = PayDateStrip->GetInterestTerms();
	double AsOfDate = GetAsOfDate().GetJulian();

	if( !itsIRModel )
	{
		/// In this case, we used the default discount factor (which assumes
		/// that forward rates happen

		size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		size_t nbFlows = PayDates->size();
		double dayCountCorrection;  /// A365, etc.
		ARM_GP_VectorPtr DFs;
		size_t i;

		ARM_GP_Vector * result = new ARM_GP_Vector( nbStates );
		ARM_GP_Vector::iterator resultIter,DFsIter;

		/// Loop on the cash flows
		for( i = 0 ; i<nbFlows ; i++)
		{
			dayCountCorrection = InterestTerms->Elt(i);
			DFs = DiscountFactor( curveName, evalTime, PayDates->Elt(i) - AsOfDate, states);
			DFsIter = DFs->begin();

			/// Loop on the states vector
			for( resultIter = result->begin() ; resultIter != result->end() ; resultIter++ , DFsIter++ )
				(*resultIter) += dayCountCorrection * (*DFsIter);
		}

		return ARM_GP_VectorPtr ( result );

	}
	else
	{
		/// In this case, we have an IRModel. We ask it to compute the annuity

		/// Clone and change of newPayDates (this is a bit porky)
		ARM_GP_Vector * newPayDates = new ARM_GP_Vector( *PayDates );
		for( ARM_GP_Vector::iterator iter = newPayDates->begin() ; iter != newPayDates->end() ; iter++ )
			*iter -= AsOfDate;
		ARM_GP_VectorPtr temp( newPayDates ); 

		return itsIRModel->Annuity( curveName, evalTime,*newPayDates,*InterestTerms,states);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: DefaultYoYSwapLeg
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::DefaultYoYSwapLeg(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	ARM_GP_Vector * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	ARM_GP_Vector * numCPIDates = numDateStrip->GetResetDates();
	ARM_GP_Vector * denomCPIDates = denomDateStrip->GetResetDates();
	ARM_GP_Vector * paymentDates = numDateStrip->GetFlowEndDates();
	ARM_GP_Vector * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double spread = itsSpread;
	double AsOfDate = GetAsOfDate().GetJulian();

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new ARM_GP_Vector( nbStates ) );
	ARM_GP_VectorPtr numDFs, CPIRatios;
	ARM_GP_Vector::iterator iterResult, iterDFs, iterCPIRatios; 

	// Compute the floating leg
	for( i = 0 ; i< nbFlows; i++ )
	{
		numDFs = DiscountFactor( irCurveName, evalTime, (*paymentDates)[i] - AsOfDate, states );
		CPIRatios = ForwardCPIRatio( infCurveName, evalTime, (*numCPIDates)[i]-(*denomCPIDates)[i], (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate,states );
		iterDFs = numDFs->begin();
		iterCPIRatios = CPIRatios->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterDFs++, iterCPIRatios++ )
			(*iterResult) += ((*iterCPIRatios)-1.+spread)*(*iterDFs) * dayCountCorrection;
	}


	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: DefaultOATSwapLeg
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::DefaultOATSwapLeg(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	const ARM_DateStripPtr& numDateStrip,
	const ARM_DateStripPtr& denomDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;

	/// This vector contains the num CPIs fixingDates (and not really the FlowStartDates!)
	ARM_GP_Vector * numCPIDCFDates = numDateStrip->GetFlowStartDates();
	ARM_GP_Vector * numCPIDates = numDateStrip->GetResetDates();
	ARM_GP_Vector * denomCPIDates = denomDateStrip->GetResetDates();
	ARM_GP_Vector * paymentDates = numDateStrip->GetFlowEndDates();
	ARM_GP_Vector * NumInterestTerms = numDateStrip->GetInterestTerms();
	double dayCountCorrection;
	double AsOfDate = GetAsOfDate().GetJulian();
	double denomCPIDate = (*denomCPIDates)[0];

	size_t i, nbFlows = numCPIDates->size();
	ARM_GP_VectorPtr result( new ARM_GP_Vector( nbStates ) );
	ARM_GP_VectorPtr numDFs, CPIRatios;
	ARM_GP_Vector::iterator iterResult, iterDFs, iterCPIRatios; 

	// Compute the floating leg
	for( i = 0 ; i< nbFlows; i++ )
	{
		numDFs = DiscountFactor( irCurveName, evalTime, (*paymentDates)[i] - AsOfDate, states );
		CPIRatios = ForwardCPIRatio( infCurveName, evalTime, (*numCPIDates)[i]-denomCPIDate, (*numCPIDates)[i]-AsOfDate, (*numCPIDCFDates)[i]-AsOfDate,states );
		iterDFs = numDFs->begin();
		iterCPIRatios = CPIRatios->begin();
		dayCountCorrection = NumInterestTerms->Elt(i);

		for( iterResult = result->begin() ; iterResult != result->end() ; iterResult++, iterDFs++, iterCPIRatios++ )
			(*iterResult) += ((*iterCPIRatios)+itsSpread)*(*iterDFs) * dayCountCorrection;
	}


	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: YoYSwapRate
///	Returns: ARM_GP_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::YoYSwapRate(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	const ARM_DateStripPtr& FloationgNumDateStrip,
	const ARM_DateStripPtr& FloationgDenomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result = DefaultYoYSwapLeg(irCurveName, infCurveName, evalTime, FloationgNumDateStrip,FloationgDenomDateStrip, itsSpread,states);
	ARM_GP_VectorPtr annuity = DefaultAnnuity(irCurveName, evalTime,fixedDateStrip, states);

	ARM_GP_Vector::iterator iter,iterAnnuity;

	iterAnnuity = annuity->begin();
	for( iter = result->begin() ; iter != result->end() ; iter++,iterAnnuity++ )
		(*iter) /= (*iterAnnuity);

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: OATSwapRate
///	Returns: ARM_GP_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::OATSwapRate(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	const ARM_DateStripPtr& FloationgNumDateStrip,
	const ARM_DateStripPtr& FloationgDenomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result = DefaultOATSwapLeg(irCurveName, infCurveName, evalTime, FloationgNumDateStrip,FloationgDenomDateStrip, itsSpread, states);
	ARM_GP_VectorPtr annuity = DefaultAnnuity(irCurveName, evalTime,fixedDateStrip, states);

	ARM_GP_Vector::iterator iter,iterAnnuity;

	iterAnnuity = annuity->begin();
	for( iter = result->begin() ; iter != result->end() ; iter++,iterAnnuity++ )
		(*iter) /= (*iterAnnuity);

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: YoYSwap
///	Returns: ARM_GP_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::YoYSwap(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	const ARM_DateStripPtr& FloationgNumDateStrip,
	const ARM_DateStripPtr& FloationgDenomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsSpread,
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result = DefaultYoYSwapLeg(irCurveName, infCurveName, evalTime, FloationgNumDateStrip,FloationgDenomDateStrip, itsSpread,states);
	ARM_GP_VectorPtr annuity = DefaultAnnuity(irCurveName, evalTime,fixedDateStrip, states);

	ARM_GP_Vector::iterator iter,iterAnnuity;

	iterAnnuity = annuity->begin();
	for( iter = result->begin() ; iter != result->end() ; iter++,iterAnnuity++ )
		(*iter) -= Strike*(*iterAnnuity);

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingModelInflation
///	Routine: OATSwap
///	Returns: ARM_GP_VectorPtr
///	Action : 
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_PricingModelInflation::OATSwap(
	const string& irCurveName,
	const string& infCurveName, 
	double evalTime,
	double Strike,
	double FloatMargin, 
	const ARM_DateStripPtr& FloationgNumDateStrip,
	const ARM_DateStripPtr& FloationgDenomDateStrip,
	const ARM_DateStripPtr& fixedDateStrip,
	double itsCoupon,
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result = DefaultOATSwapLeg(irCurveName, infCurveName, evalTime, FloationgNumDateStrip,FloationgDenomDateStrip, FloatMargin, states);
	ARM_GP_VectorPtr annuity = DefaultAnnuity(irCurveName, evalTime,fixedDateStrip, states);

	ARM_GP_Vector::iterator iter,iterAnnuity;

	iterAnnuity = annuity->begin();
	for( iter = result->begin() ; iter != result->end() ; iter++,iterAnnuity++ )
		(*iter) -= Strike*(*iterAnnuity);

	return result;
}

CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
