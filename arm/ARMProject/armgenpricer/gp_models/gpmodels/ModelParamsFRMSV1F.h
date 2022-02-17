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


#ifndef _INGPMODELS_MODELPARAMSFRMSV1F_H
#define _INGPMODELS_MODELPARAMSFRMSV1F_H

#include "gpbase/port.h"
#include "gpinfra/ModelParams.h"
#include "gpinfra/modelparamsvec.h"
#include "gpbase/gpvector.h"



CC_BEGIN_NAMESPACE( ARM )

class ARM_ModelParamsSFRM;
class ARM_VolParams;
//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMSV
// \brief Interface class for model parameters of MSV models
//-----------------------------------------------------------------------------
class ARM_ModelParamsFRMSV1F : public ARM_ModelParamsVec 
{
private:
	size_t itsSFRM_param_index;
	size_t itsSV_param_index;
	void ValidateModelParams();
public:
	ARM_ModelParamsFRMSV1F( const ARM_ModelParamsFRMSV1F& rhs );
	ARM_ModelParamsFRMSV1F( const vector<ARM_ModelParams*>& paramsVec );
	virtual ~ARM_ModelParamsFRMSV1F();
    ARM_ModelParamsFRMSV1F& operator = (const ARM_ModelParamsFRMSV1F& rhs);
	virtual ARM_Object* Clone() const;


	/// Discretization Dates of the model params
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 2; }
	
	ARM_ModelParamsSFRM* GetSFRMParams() const;
	ARM_VolParams* GetVolParams() const;

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

