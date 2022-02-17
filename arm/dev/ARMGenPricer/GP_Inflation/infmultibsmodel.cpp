/*!
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *	\file infmultibsmodel.cpp
 *
 *  \brief inflation multi Black Scholes model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2004
 */

#include <glob/firsttoinc.h>

/// gpinflation
#include "gpinflation/infmultibsmodel.h"
#include "gpinflation/infidx.h"
#include "gpinflation/infcurv.h"
#include "gpinflation/infbsmodel.h"
#include "gpinflation/infswaption.h"
#include "gpinflation/assetinfo.h"

/// gpclosedforms
#include "gpclosedforms/spreadoption_lognormal_interface.h"

/// gpbase
#include "gpbase/gpvector.h"

/// kernel
#include <crv/correlmanager.h>

CC_USING_NS(std,pair)

CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: Constructor
///	Returns: 
///	Action : Fills the model name map
////////////////////////////////////////////////////
ARM_InfMultiBSModel::ARM_InfMultiBSModel( const vector<ARM_InfBSModelPtr>& models )
:	ARM_InfBSModel( *models[0])
{
	for( size_t i=0; i<models.size(); ++i)
	{
		/// validate that all model have the same curve!
		if( strcmp( models[i]->GetZeroCurve()->GetCurrencyUnit()->GetCcyName(), models[0]->GetZeroCurve()->GetCurrencyUnit()->GetCcyName() ) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": only single interest rates inflation model allowed! found " 
				+ models[0]->GetZeroCurve()->GetCurrencyUnit()->GetCcyName() + " and "
				+ models[i]->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());

		if( !itsModelMap.insert( pair<const string, ARM_InfBSModelPtr >( models[i]->GetInfFwdCurv()->GetInfIdxName(), models[i] ) ).second )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
				ARM_USERNAME + ": duplicated model for key " +  models[i]->GetInfFwdCurv()->GetInfIdxName() );
	}
	SetName(ARM_INFMULTIBSMODEL);
};


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: Copy Constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfMultiBSModel::ARM_InfMultiBSModel(  const ARM_InfMultiBSModel& rhs )
:	ARM_InfBSModel( rhs ), itsModelMap( rhs.itsModelMap )
{}



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: operator=
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfMultiBSModel& ARM_InfMultiBSModel::operator=( const ARM_InfMultiBSModel& rhs )
{
	if( this != &rhs )
	{
		ARM_InfBSModel::operator =( rhs );
		itsModelMap = rhs.itsModelMap;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_InfMultiBSModel::~ARM_InfMultiBSModel()
{}



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: Clone
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Object* ARM_InfMultiBSModel::Clone()
{
	return new ARM_InfMultiBSModel(*this);
}



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: View
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_InfMultiBSModel::View(char* id, FILE* ficOut )
{
    FILE* fOut;
    char fOutName[200];
	
	/// first determine that the file is not already opened
    if ( NULL == ficOut )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

	/// loop over the models!
	fprintf( fOut," =========> Inflation Multi Black Scholes\n\n" );
	fprintf( fOut," Quick View\n\n" );
	
	SimpleInfBSModelMap::iterator iter;
	size_t i;
	for( iter= itsModelMap.begin(), i=0; iter != itsModelMap.end(); ++iter, ++i )
		fprintf( fOut," %d) Model for %s\n", i, (*iter).second->GetInfFwdCurv()->GetInfIdxName().c_str() );

	fprintf( fOut,"\n\n\n Detailled model\n" );
	for( iter= itsModelMap.begin(), i=0; iter != itsModelMap.end(); ++iter, ++i )
	{
		fprintf( fOut,"\n\n\n\n" );
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut,"				%d) Model for %s\n", i, (*iter).second->GetInfFwdCurv()->GetInfIdxName().c_str() );
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut," ========================================================================================================\n");
		fprintf( fOut,"\n\n" );

		(*iter).second->View( id, fOut );
	}

	/// to allow to have nested view
    if ( NULL == ficOut )
       fclose(fOut);
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: GetCorrespondingInflationModel
///	Returns: 
///	Action : Find the corresponding model
////////////////////////////////////////////////////
ARM_InfBSModelPtr ARM_InfMultiBSModel::GetCorrespondingInflationModel( ARM_InfIdx* infIdx ) const
{
	SimpleInfBSModelMap::const_iterator found = itsModelMap.find( infIdx->GetIndexName() );
	if( found == itsModelMap.end() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": could not find a model for " + infIdx->GetIndexName() );
	return (*found).second;
}



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: FwdCPIRatio
///	Returns: return a vector containing various information about the fwd cpi ratio pricing
///	Action : prices a fwd cpi ratio
////////////////////////////////////////////////////
ARM_GP_Vector ARM_InfMultiBSModel::FwdCPIRatio( 
	const ARM_Date& numDate,
	const ARM_Date& denomDate,
	const ARM_Date& paymentDate,
	double multiple,
	double spread,
	long dailyInterpType,
	double denomFixing,
	ARM_InfIdx* infIdx )
{
	return GetCorrespondingInflationModel( infIdx )->FwdCPIRatio( 
		numDate, denomDate, paymentDate, multiple,
		spread, dailyInterpType, denomFixing, infIdx );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: View
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::DiscountedCPI(
	const ARM_Date& resetDate, 
	const ARM_Date& paymentDate, 
	long dailyInterpType,
	const string& CPILag,
	const string& DCFLag,
	ARM_InfIdx* infIdx )
{
	return GetCorrespondingInflationModel( infIdx )->DiscountedCPI(
		resetDate, paymentDate, dailyInterpType, CPILag, DCFLag, infIdx );
}

////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: GetModelTimeWPublishLag
///	Returns: double
///	Action : Get the model time with a publishing lag
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx )
{
	return GetCorrespondingInflationModel( infIdx )->GetModelTimeWPublishLag( date, infIdx );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: SingleAssetOptionPrice
///	Returns: double (the price)
///	Action : computes a single asset option
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::SingleAssetOptionPrice(
		double CPIForward,
		double strike,
		int callput,
		double discounting, 
		ARM_InfIdx* infIdx,
		ARM_Object* optionContext,
		StoreInfoObj& storeInfo	)
{
	return GetCorrespondingInflationModel( infIdx )->SingleAssetOptionPrice(
		CPIForward, strike, callput, discounting, infIdx,
		optionContext, storeInfo );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: TwoAssetsOptionPrice
///	Returns: double (the price)
///	Action : computes a two asset options with correlation
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::TwoAssetsOptionPrice(
	double CPIForward,
	double secondAssetFwd, 
	double strike,
	int callput, 
	double discounting, 
	ARM_InfIdx* infIdx,
	ARM_IRIndex* secondIndex,
	ARM_Object* optionContext, 
	StoreInfoObj& storeInfo	)
{
	if( ARM_InfSwaptionContext* infSwaptionContext = dynamic_cast<ARM_InfSwaptionContext*>(optionContext) )
	{
		double tenor	= infSwaptionContext->GetOptionMaturity();
		double volAsset1= ARM_InfBSModel::ComputeVolForSwaption(strike, 
			infSwaptionContext->GetExpiryAsset1(),
			infSwaptionContext->GetTenorAsset1(),
			GetCorrespondingInflationModel( infIdx )->GetInfSwoptVolCurve(),
			"inflation" );
		
		ARM_VolCurve* secondAssetVolCurve = NULL;
		string secondAssetString, mktTag;
		double correl;

		if( dynamic_cast<ARM_InfIdx*>(secondIndex) )
		{
			secondAssetVolCurve = GetCorrespondingInflationModel( (ARM_InfIdx*) secondIndex )->GetInfSwoptVolCurve();
			secondAssetString	= "inflation";
			mktTag				= "INF/INF_SWOPT";
			correl = GetInfInfCorrel( infSwaptionContext->GetOptionMaturity(), 
				infSwaptionContext->GetTenorAsset1(), infIdx, (ARM_InfIdx*) secondIndex, mktTag );
		}
		else
		{
			secondAssetVolCurve = GetIRSwoptVolCurve();
			secondAssetString = "interest rates";
			mktTag				= "INF/IR_SWOPT";
			correl = GetInfIRCorrel( infSwaptionContext->GetOptionMaturity(), 
				infSwaptionContext->GetTenorAsset1(), infIdx, secondIndex, mktTag );
		}

		double volAsset2= ARM_InfBSModel::ComputeVolForSwaption( 
			strike,
			infSwaptionContext->GetExpiryAsset2(),
			infSwaptionContext->GetTenorAsset2(),
			secondAssetVolCurve,
			secondAssetString );

		double optionType= 1; 

		/// store the appropriate information
		double data[10];
		data[0]= CPIForward * CC_NS( ARM_Constants, rateBase );
		data[1]= volAsset1  * CC_NS( ARM_Constants, volBase );
		data[2]= infSwaptionContext->GetTenorAsset1();
		data[3]= secondAssetFwd * CC_NS( ARM_Constants, rateBase );
		data[4]= volAsset2 * CC_NS( ARM_Constants, volBase );
		data[5]= infSwaptionContext->GetTenorAsset2();
		data[6]= correl * CC_NS( ARM_Constants, correlBase );
		data[7]= discounting;
		data[8]= infSwaptionContext->GetOptionMaturity();
		data[9]= strike * CC_NS( ARM_Constants, rateBase );
		storeInfo.Store( data );

		/// Uses the closed formula library on the spread option valuation
		int n = 120;
		double value = Export_LogNormal_SpreadOption(CPIForward, secondAssetFwd, volAsset1, volAsset2, correl, strike, tenor, callput, optionType, n);
		return ( discounting * value * CC_NS( ARM_Constants, rateBase ));
	}
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			ARM_USERNAME + ": only swaption requires two asset pricing!");
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: GetInfIRCorrel
///	Returns: double
///	Action : Get the inflatin interest rate correlation
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::GetInfIRCorrel( double TjFromAsOf, double TiFromAsOf, 
	ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const
{
	return GetCorrespondingInflationModel( infIdx )->GetInfIRCorrel( TjFromAsOf, TiFromAsOf, 
		infIdx, otherIndex, mktTag  );
}


////////////////////////////////////////////////////
///	Class  : ARM_InfBSModel
///	Routine: GetInfIRCorrelMatrix
///	Returns: 
///	Action : gets the inflation interest rates correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfMultiBSModel::GetInfIRCorrelMatrix( ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const
{
	return GetCorrespondingInflationModel( infIdx )->GetInfIRCorrelMatrix( infIdx, otherIndex, mktTag  );
}



////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: GetInfInfCorrelMatrix
///	Returns: ARM_CorrelMatrix*
///	Action : Get the inflatin interest rate correlation
////////////////////////////////////////////////////
ARM_CorrelMatrix* ARM_InfMultiBSModel::GetInfInfCorrelMatrix( 
	ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag) const
{
	string indexName1	= infIdx1->GetIndexName();
	string indexName2	= infIdx2->GetIndexName();
	string intraMktTag	= indexName1 < indexName2 ? indexName1 + "_" + indexName2 :	indexName2 + "_"+ indexName1;

	/*if( GetCorrelManager() )
		return GetCorrelManager()->ComputeCorrelData( mktTag,  intraMktTag );
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": correlmanager not found! " );*/

	if( GetCorrespondingInflationModel( infIdx1 )->GetCorrelManager())
		if(GetCorrespondingInflationModel( infIdx1 )->GetCorrelManager()->GetCorrelData( mktTag,  intraMktTag ))
			return GetCorrespondingInflationModel( infIdx1 )->GetCorrelManager()->GetCorrelData( mktTag,  intraMktTag );
		else if(GetCorrespondingInflationModel( infIdx2 )->GetCorrelManager())
			if (GetCorrespondingInflationModel( infIdx2 )->GetCorrelManager()->GetCorrelData( mktTag,  intraMktTag ))
				return GetCorrespondingInflationModel( infIdx2 )->GetCorrelManager()->GetCorrelData( mktTag,  intraMktTag );
			else
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": correl matrix " + intraMktTag + "not found! " );
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": correlmanager not found for ! "+indexName2+"Model" );
	else
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			ARM_USERNAME + ": correlmanager not found! " +indexName1+"Model");
			
}


////////////////////////////////////////////////////
///	Class  : ARM_InfMultiBSModel
///	Routine: GetInfInfCorrel
///	Returns: double
///	Action : Get the inflatin interest rate correlation
////////////////////////////////////////////////////
double ARM_InfMultiBSModel::GetInfInfCorrel( double TjFromAsOf, double TiFromAsOf, 
	ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag ) const
{
	return GetInfInfCorrelMatrix(infIdx1, infIdx2, mktTag )->ComputeCorrelData( TjFromAsOf, TiFromAsOf ) 
		/ CC_NS(ARM_Constants,correlBase);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

