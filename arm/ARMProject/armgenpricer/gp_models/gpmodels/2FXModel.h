/*!
 *
 * Copyright (c) CDC IXIS CIB 2005 Paris
 *
 *	\file 2FXModel.h
 *
 *  \brief
 *
 *  \brief 2 fx multi assets model
 *
 *	\author  K.BELKHEIR
 *	\version 1.0
 *	\date Octobre 2006
 */



#ifndef _INGPMODELS_2FXMODEL_H
#define _INGPMODELS_2FXMODEL_H

#include "gpbase/port.h"
#include "gpbase/assignop.h"
#include "MultiAssets.h"
#include "gpmodels/EqFxBase.h"

#include "gpinfra/modelnamemap.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricingContext;
class ARM_GaussReplic2D;

////////////////////////////////////////////////////////////////////////////////
//
// This class defines a multi asset with 2 models:
// _ a 1 EqFXBase model associated to a FX1 rate
// _ a 1 EqFXBase model associated to a FX2 rate
// The two FX must have a common currency
//
//
////////////////////////////////////////////////////////////////////////////////

class ARM_2FXModel :	public ARM_MultiAssetsModel
{
public:
	//standard constructor
	ARM_2FXModel(const ARM_ModelNameMap&	modelNameMap, 
		const ARM_CurveMatrix& correlationMatrix  = ARM_CurveMatrix());

	//copy constructor
	ARM_2FXModel(const ARM_2FXModel& rhs);
	ASSIGN_OPERATOR(ARM_2FXModel)
	virtual ~ARM_2FXModel()
	{
	}

	//Initialize the model

    virtual ARM_PricingStatesPtr Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos);

	///Financial vectorial functions

	// Call function (vectorial strike version)
	virtual ARM_VectorPtr CallVectorial(
		const string& modelName,
        double evalTime,
	    double expiryTime,
	    double settlementTime,
		const ARM_GP_Vector& strikePerState,
		int callPut,
	    double payTime,
        const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// CallSpread function (vectorial strike version)
	virtual ARM_VectorPtr CallSpreadVectorial(
		const string& Model1Name,
		const string& Model2Name,
		double evalTime,
		double expiryTime,
		double settlementTime1,
		double settlementTime2,
		double payTime,
		const ARM_GP_Vector& strikePerState,
		double alpha,
		double beta,	
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	// RangeAccrual function (vectorial strike version)
	virtual ARM_VectorPtr ARM_2FXModel::RangeAccrualVectorial(
		const string& model1Name,
		double evalTime,
		double startTime,
		double endTime,
		const ARM_GP_Vector& fixingTimes,
		double payTime,
		const ARM_GP_Vector& downBarrierVect,
		const ARM_GP_Vector& upBarrierVect,
		const ARM_GP_Vector& notionalVect,
		const ARM_PricingStatesPtr& states,
        ARM_PricingContext* context=NULL) const;

	//Price of the call quanto
	double CallPrice( const ARM_EqFxBase& payoffModel, 
		const ARM_EqFxBase& quantoModel,
		ARM_GaussReplic2D* gaussReplic,
		double  fwd_payoff,
		double  fwd_quanto,
		double	expiryTime,
		const ARM_GP_Matrix& IntegParameters) const;

    /// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_2FXModel(*this);}
	virtual string ExportShortName() const { return "L2FXM";}
	virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
