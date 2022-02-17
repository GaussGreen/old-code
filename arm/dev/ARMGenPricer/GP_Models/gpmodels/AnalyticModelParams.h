/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: SurfaceModelParam.h,v $
 * Revision 1.1  2004/07/30 09:52:19  ebenhamou, ocroissant
 * Initial revision
 *
 *
 */

/*! \file SurfaceModelParam.h
 *
 *  \brief class for model params that are according to a surface
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_ANALYTICMODELPARAMS_H
#define _INGPMODELS_ANALYTICMODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_AnalyticModelParams : public ARM_ModelParams
{
public:
	ARM_AnalyticModelParams( const ARM_ModelParamVector& params = ARM_ModelParamVector() );
	ARM_AnalyticModelParams( const ARM_AnalyticModelParams& rhs );
	ARM_AnalyticModelParams& operator=( const ARM_AnalyticModelParams& rhs );
	virtual ~ARM_AnalyticModelParams() {};
	void ValidateSurfaceModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const ;
	void ValidateCurveModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const ;
	void ValidateSurfaceListModelParam( const string& modelParamsName, ARM_ModelParamType::ParamNb paramNb ) const ;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
