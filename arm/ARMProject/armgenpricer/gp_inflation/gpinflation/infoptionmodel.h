/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file InfHybModel.h
 *  \brief file to factorize code for the inflation option model
 *	\author  E Benhamou
 *	\version 1.0
 *	\date July 2004
 */


#ifndef _INGPINFLATION_INFOPTIONMODEL_H
#define _INGPINFLATION_INFOPTIONMODEL_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "infcurvmodel.h"

class ARM_IRIndex;

CC_BEGIN_NAMESPACE( ARM )

class ARM_InfIdx;

struct StoreInfoObj
{
	///forces redefinition
	virtual void Store( double* datas )=0;
};


//////////////////////////////////////////////////////////
/// \class dummy class for pricing of simple inflation option
/// \author  Eric Benhamou
/// \version 1.0
/// \date August 2003
////////////////////////////////////////////////////////////
class InfOptionModel : public InfFwdModel
{
public:
	InfOptionModel( ARM_InfCurv* infFwdCurv = NULL ) : InfFwdModel( infFwdCurv ) {} ;
	InfOptionModel( const InfOptionModel& rhs ) : InfFwdModel( rhs ) {}
	InfOptionModel& operator=( const InfOptionModel& rhs )
	{
		if( this != &rhs )
			InfFwdModel::operator=( rhs );
		return *this;
	}
	virtual ~InfOptionModel(){};

	virtual double SingleAssetOptionPrice(double CPIForward, double strike,
		int callPut, double discounting, ARM_InfIdx* infIdx, ARM_Object* optionContext, StoreInfoObj& storeInfo	) = 0;
	
	virtual double TwoAssetsOptionPrice(double CPIForward, double secondAssetFwd, double strike,
		int callPut, double discounting, ARM_InfIdx* infIdx, ARM_IRIndex* secondIndex,
		ARM_Object* optionContext, StoreInfoObj& storeInfo	) = 0;

	virtual double GetModelTimeWPublishLag( const ARM_Date& date, ARM_InfIdx* infIdx) = 0;
};



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
