/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfHybModel.h
 *  \brief file for the inflation hybrid model
 *	\author  E Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPINFLATION_INFMULTIBSMODEL_H
#define _INGPINFLATION_INFMULTIBSMODEL_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "infbsmodel.h"
#include "gpbase/port.h"
#include "typedef.h"
#include <string>
#include <map>
CC_USING_NS( std, map )
CC_USING_NS( std, less )
CC_USING_NS( std, string )


CC_BEGIN_NAMESPACE( ARM )



/// forward declaration
class ARM_InfBSModel;
typedef map< string, ARM_InfBSModelPtr, less< string> > SimpleInfBSModelMap;

class ARM_InfMultiBSModel : public ARM_InfBSModel
{
	/// the equivalent of a model name map object but with old type model!
	SimpleInfBSModelMap itsModelMap;
	
	
public:
	ARM_InfMultiBSModel( const vector<ARM_InfBSModelPtr>& models );
	ARM_InfMultiBSModel( const ARM_InfMultiBSModel& rhs );
	ARM_InfMultiBSModel& operator=( const ARM_InfMultiBSModel& rhs );
	virtual ~ARM_InfMultiBSModel();

	/// standard ARM_Support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// pricing function
	/// standard pricing of forward CPI
	virtual ARM_GP_Vector FwdCPIRatio( 
		const ARM_Date& numDate,
		const ARM_Date& denomDate,
		const ARM_Date& paymentDate,
		double multiple,
		double spread,
		long dailyInterpType,
		double denomFixing		= GETDEFAULTVALUE,
		ARM_InfIdx* infIdx		= NULL );

	virtual double DiscountedCPI(
		const ARM_Date& resetDate, 
		const ARM_Date& paymentDate, 
		long dailyInterpType,
		const string& CPILag,
		const string& DCFLag,
		ARM_InfIdx* infIdx		= NULL );

	virtual ARM_INF_PRICING_INFO CanPriceInflation() const { return PRICE_FWDNOPTION; }
	virtual double GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx );
	virtual double SingleAssetOptionPrice(
		double CPIForward,
		double strike,
		int callPut,
		double discounting, 
		ARM_InfIdx* infIdx,
		ARM_Object* optionContext,
		StoreInfoObj& storeInfo	);

	virtual double TwoAssetsOptionPrice(
		double CPIForward,
		double secondAssetFwd, 
		double strike,
		int callPut, 
		double discounting, 
		ARM_InfIdx* infIdx,
		ARM_IRIndex* secondIndex,
		ARM_Object* optionContext, 
		StoreInfoObj& storeInfo	);

	virtual ARM_CorrelMatrix* GetInfIRCorrelMatrix( ARM_InfIdx* infIdx, 
		ARM_IRIndex* otherIndex, const string& mktTag  ) const;
	virtual double GetInfIRCorrel( double TjFromAsOf, double TiFromAsOf, 
		ARM_InfIdx* infIdx, ARM_IRIndex* otherIndex, const string& mktTag  ) const;

	virtual ARM_CorrelMatrix* GetInfInfCorrelMatrix( ARM_InfIdx* infIdx1, 
		ARM_InfIdx* infIdx2, const string& mktTag ) const;
	virtual double GetInfInfCorrel( double TjFromAsOf, double TiFromAsOf, 
		ARM_InfIdx* infIdx1, ARM_InfIdx* infIdx2, const string& mktTag ) const;

	ARM_InfBSModelPtr GetCorrespondingInflationModel( ARM_InfIdx* infIdx ) const;
	ARM_InfBSModelPtr GetCorrespondingIRModel( ARM_IRIndex* irIndex ) const;

};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
