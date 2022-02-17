/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file typedef.h
 *  \brief global typedef of namespace ARM
 *	\author  E Benhamou
 *	\version 1.0
 *	\date June 2004
 */


#ifndef _INGPMODEL_TYPEDEF_H
#define _INGPMODEL_TYPEDEF_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/countedptr.h"
#include "gpbase/curvetypedef.h"

/// gpmodels
#include "enummodel.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_ModelParams;
class ARM_ModelParamsQ1F;
class ARM_CEV_ModelParams;
class ARM_BS_ModelParams;
class ARM_SABR_ModelParams;
class ARM_Heston_ModelParams;
class ARM_MarkovSV;
class ARM_SmiledModel_Fx;
class ARM_2IRFXModel;
class ARM_LocalFunctional;

/// template class
template <typename T> class ARM_ModelParams_EqFxBase_T;
template <typename T> class ARM_ModelParams_Fx_T;
template <typename T> class ARM_ModelParams_Eq_T;


///struct
struct ARM_PricingModelType;
typedef ARM_PricingModelType::ModelType ARM_ModelType;

struct ARM_MultiAssetsModelType;
typedef ARM_MultiAssetsModelType::MultiAssetsModelType ARM_MultiAssetsType;

/// THE TYPEDEF
typedef ARM_CountedPtr< ARM_Curve >					ARM_CurvePtr;

/// BS Model
typedef ARM_ModelParams_Eq_T<ARM_BS_ModelParams>			ARM_ModelParamsBS_Eq;
typedef ARM_ModelParams_Fx_T<ARM_BS_ModelParams>			ARM_ModelParamsBS_Fx;

/// Q1F Model
typedef ARM_ModelParams_Eq_T<ARM_ModelParamsQ1F>			ARM_ModelParamsQ1F_Eq;

/// SABR Model
//typedef ARM_ModelParams_EqFx_T<ARM_SABR_ModelParams>  ARM_ModelParamsSABR_EqFx;
typedef ARM_ModelParams_Eq_T<ARM_SABR_ModelParams>	    ARM_ModelParamsSABR_Eq;
typedef ARM_ModelParams_Fx_T<ARM_SABR_ModelParams>		ARM_ModelParamsSABR_Fx;

// Mixture
typedef ARM_ModelParams_Fx_T<ARM_ModelParams>			ARM_ModelParamsMixture_Fx;

/// MSV Model
typedef ARM_CountedPtr< ARM_MarkovSV >					ARM_MarkovSVPtr;

typedef ARM_CountedPtr<ARM_SmiledModel_Fx>				ARM_SmiledModel_FxPtr;
typedef ARM_CountedPtr<ARM_2IRFXModel>					ARM_2IRFXModelPtr;

typedef ARM_CountedPtr<ARM_LocalFunctional>				ARM_LocalFunctionalPtr;



CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
