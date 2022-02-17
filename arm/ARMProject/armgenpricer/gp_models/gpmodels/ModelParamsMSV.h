/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ModelParamsMSV.h,v $
 * Initial revision
 *
 *
 */



/*! \file ModelParamsMSV.h
 *
 *  \brief 
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_MODELPARAMSMSV_H
#define _INGPMODELS_MODELPARAMSMSV_H

#include "gpbase/port.h"
#include "gpinfra/ModelParams.h"
#include "gpbase/gpvector.h"


CC_BEGIN_NAMESPACE( ARM )


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMSV
// \brief Interface class for model parameters of MSV models
//-----------------------------------------------------------------------------
class ARM_ModelParamsMSV : public ARM_ModelParams 
{
protected :
	mutable std::vector<double> itsSchedule;
public:
	ARM_ModelParamsMSV( const ARM_ModelParamsMSV& rhs );
	ARM_ModelParamsMSV( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsMSV();
    ARM_ModelParamsMSV& operator = (const ARM_ModelParamsMSV& rhs);

	/// Discretization Dates of the model params
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const = 0;

	/// Accessors
	const std::vector<double>& GetSchedule() const { return itsSchedule;}

	virtual void GenerateSchedule();

	/// How many factors?
	virtual size_t FactorCount() const = 0;

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

