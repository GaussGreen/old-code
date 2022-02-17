/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_local_gp_model.cpp
 *
 *  \brief file for the model part of the generic pricer local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include "firstToBeIncluded.h"

///////////////////////////////////////////////
/// WARNING include these headers FIRST
/// because it uses the new style headers
/// CCxl uses cstyle headers and leads to conflict
/// if defined first
///////////////////////////////////////////////

/// gpbase
#include <GP_Base\gpbase\autocleaner.h>
#include <GP_Base\gpbase\singleton.h>
#include <GP_Base\gpbase\curve.h>
#include <GP_Base\gpbase\curvematrix.h>
#include <GP_Base\gpbase\gpmatrix.h>
#include <GP_Base\gpbase\datestrip.h>
#include <GP_Base\gpbase\surface.h>
#include <GP_Base\gpbase\gplinalgtypedef.h>
#include <GP_Base\gpbase\gpvector.h>
#include <GP_Base\gpbase\surfacetypedef.h>
#include <GP_Base\gpbase\cloneutilityfunc.h>
#include <GP_Base\gpbase\stringmanip.h>
#include <GP_Base\gpbase\env.h>

/// gpclosedform
#include "gpclosedforms/gaussian_integrals.h"

/// gpinfra
#include <GP_Infra\gpinfra\typedef.h>
#include <GP_Infra\gpinfra\argconvdefault.h>
#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\modelparams.h>
#include <GP_Infra\gpinfra\pricingstates.h>
#include <GP_Infra\gpinfra\numeraire.h>
#include <GP_Infra\gpinfra\numerairefactory.h>
#include <GP_Infra\gpinfra\modelparamtype.h>
#include <GP_Infra\gpinfra\curvemodelparam.h>
#include <GP_Infra\gpinfra\surfacemodelparam.h>
#include <GP_Infra\gpinfra\surfacelistmodelparam.h>
#include <GP_Infra\gpinfra\correlmatparam.h>
#include <GP_Infra\gpinfra\modelnamemap.h>
#include <GP_Infra\gpinfra\mktdatamanagerrep.h>
#include <GP_Infra\gpinfra\modelparamsvec.h>

/// gpcalib
#include <GP_Calib\gpcalib\modelparamsfactory.h>
#include <GP_Calib\gpcalib\vanillaarg.h>
#include <GP_Calib\gpcalib\vanillaswaption.h>
#include <GP_Calib\gpcalib\kerneltogp.h>
#include <GP_Calib\gpcalib\densityfunctors.h>
#include <GP_Calib\gpcalib\vanillasecuritydensity.h>
#include <GP_Calib\gpcalib\vanilladensityfactory.h>

/// gpmodels
#include <GP_Models\gpmodels\IRFwdMod.h>
#include <GP_Models\gpmodels\InfFwdMod.h>
#include <GP_Models\gpmodels\InflationEquity.h>
#include <GP_Models\gpmodels\InflationEquityModelParams.h>
#include <GP_Models\gpmodels\HW1F.h>
#include <GP_Models\gpmodels\MSV1F.h>
#include <GP_Models\gpmodels\HWSV1F.h>
#include <GP_Models\gpmodels\HWSV2F.h>
#include <GP_Models\gpmodels\ModelParamsHWSV.h>
#include <GP_Models\gpmodels\FRMSV1F.h>
#include <GP_Models\gpmodels\ModelParamsHW1F.h>
#include <GP_Models\gpmodels\HW2F.h>
#include <GP_Models\gpmodels\ModelparamsHW2F.h>
#include <GP_Models\gpmodels\MarkovFunctional.h>
#include <GP_Models\gpmodels\ModelParamsMF.h>
#include <GP_Models\gpmodels\ModelParamsHW1F.h>
#include <GP_Models\gpmodels\ModelParamsMSV1F.h>
#include <GP_Models\gpmodels\ModelParamsFRMSV1F.h>
#include <GP_Models\gpmodels\VolDiffusionParams.h>
#include <GP_Models\gpmodels\SFRM.h>
#include <GP_Models\gpmodels\VanillaSwaptionArgSFRM.h>
#include <GP_Models\gpmodels\ModelParamsSFRMFactory.h>
#include <GP_Models\gpmodels\ModelParamsSmiledFRM.h>
#include <GP_Models\gpmodels\HybridBasisFwdIR.h>
#include <GP_Models\gpmodels\ForwardMarginIR.h>
#include <GP_Models\gpmodels\ForwardMarginBasis.h>
#include <GP_Models\gpmodels\ForwardForex.h>
#include <GP_Models\gpmodels\QGM1F.h>
#include <GP_Models\gpmodels\ModelParamsQGM1F.h>
#include <GP_Models\gpmodels\QGM2F.h>
#include <GP_Models\gpmodels\ModelParamsQGM2F.h>
#include <GP_Models\gpmodels\Q1FAna_Model.h>
#include <GP_Models\gpmodels\Q1FAna_ModelParams.h>
#include <GP_Models\gpmodels\Q1F.h>
#include <GP_Models\gpmodels\ModelParamsQ1F.h>
#include <GP_Models\gpmodels\Heston_ModelParams.h>
#include <GP_Models\gpmodels\Heston_Model.h>
#include <GP_Models\gpmodels\ShiftedHeston_ModelParams.h>
#include <GP_Models\gpmodels\ShiftedHeston_Model.h>
#include <GP_Models\gpmodels\BS_ModelParams.h>
#include <GP_Models\gpmodels\BS_Model.h>
#include <GP_Models\gpmodels\CEV_ModelParams.h>
#include <GP_Models\gpmodels\CEV_Model.h>
#include <GP_Models\gpmodels\Merton_ModelParams.h>
#include <GP_Models\gpmodels\Merton_Model.h>
#include <GP_Models\gpmodels\Normal_ModelParams.h>
#include <GP_Models\gpmodels\Normal_Model.h>
#include <GP_Models\gpmodels\SABR_ModelParams.h>
#include <GP_Models\gpmodels\SABR_Model.h>
#include <GP_Models\gpmodels\SLN_ModelParams.h>
#include <GP_Models\gpmodels\SLN_Model.h>
#include <GP_Models\gpmodels\QNF.h>
#include <GP_Models\gpmodels\ModelParamsQNF.h>
#include <GP_Models\gpmodels\Q1F_FX.h>
#include <GP_Models\gpmodels\EqFx_ModelFactory.h>
#include <GP_Models\gpmodels\MultiAssetsFactory.h>
#include <GP_Models\gpmodels\MultiAssets.h>
#include <GP_Models\gpmodels\SABR_Eq.h>
#include <GP_Models\gpmodels\2IRFXModel.h>
#include <GP_Models\gpmodels\1IRFXModel.h>
#include <GP_Models\gpmodels\HWHWQtoModel.h>
#include <GP_Models\gpmodels\HWHW2FQtoModel.h>
#include <GP_Models\gpmodels\ModelParams_EqFxBase.h>
#include <GP_Models\gpmodels\Local_Model.h>
#include <GP_Models\gpmodels\Local_Normal_ModelParams.h>
#include <GP_Models\gpmodels\Local_Normal_Model.h>
#include <GP_Models\gpmodels\Local_SLN_ModelParams.h>
#include <GP_Models\gpmodels\Local_SLN_Model.h>
#include <GP_Models\gpmodels\MarketIRModel.h>
#include <GP_Models\gpmodels\SmiledFRM.h>
#include <GP_Models\gpmodels\SmiledFRMfactory.h>
#include <GP_Models\gpmodels\SmiledMM.h>
#include <GP_Models\gpmodels\ModelParams_EqFxBase.h>
#include <GP_Models\gpmodels\argconvdefault.h>
#include <GP_Models\gpmodels\Smiled_FX.h>
#include <GP_Models\gpmodels\Mixture_FX.h>
#include <GP_Models\gpmodels\ModelParamsSVBGM.h>
#include <GP_Models\gpmodels\SVBGM.h>
#include <GP_Models\gpmodels\modelparamsbgmsv1f.h>
#include <GP_Models\gpmodels\BGMSV1F.h>
#include <GP_Models\gpmodels\EqFxBase.h>
#include <GP_Models\gpmodels\bisvmm.h>
#include <GP_Models\gpModels\EqModelParams.h>
#include <GP_Models\gpModels\EqModVolSto.h>
#include <GP_Models\gpmodels\svmmspread.h>
#include <GP_Models\gpmodels\np1irnfx.h>
#include <GP_Models\gpmodels\HWxSVMMSpread.h>
#include <GP_Models\gpmodels\HWSBGMQtoModel.h>
#include <GP_Models\gpmodels\LN_Fx.h>
#include <GP_Models\gpmodels\HWSVBGMQtoModel.h>
#include <GP_Models\gpmodels\BGMSV2F.h>
#include <GP_Models\gpmodels\2IRFXSV.h>

///	GP_Calculators
#include <GP_Calculators\gpcalculators\gencalculator.h>

/// GP_Inflation
#include <GP_Inflation\gpinflation\infcurv.h>



/// gphelp
#include <GP_Help\gphelp\crmcookies.h>

/// interface
#include "ARM_local_gp_model.h"
#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"
#include "ARM_local_glob.h"
#include "ARM_local_wrapper.h"
#include "zerointspreaded.h"

/// kernel
#include "irindex.h"
#include "expt.h"
#include "refvalue.h"
#include "portfolio.h"
#include "swaption.h"

/// STL
#include <vector>
#include <string>
#include <ctime>
#include <memory>

//// using the namespace directive to access ARM object!
using ARM::ARM_DateStrip;
using ARM::ARM_IrFwdMod;
using ARM::ARM_InfFwdMod;
using ARM::ARM_InflationEquityMod;
using ARM::ARM_ModelParamsInflationEquity;
using ARM::ARM_PricingModel;
using ARM::ARM_PricingModelIR;
using ARM::ARM_ModelParam;
using ARM::ARM_ModelParamFactory;
using ARM::ARM_VanillaDensityFactor;
using ARM::ARM_VanillaSecurityDensity;
using ARM::ARM_HullWhite1F;
using ARM::ARM_MarkovSV1F;
using ARM::ARM_FRMSV1F;
using ARM::ARM_HWSV1F;
using ARM::ARM_HWSV2F;
using ARM::ARM_ModelParamsHWSV;
using ARM::ARM_ModelParamsHW1F;
using ARM::ARM_ModelParamsMSV1F;
using ARM::ARM_ModelParamsFRMSV1F;
using ARM::ARM_VolParams;
using ARM::ARM_ModelParamsHW1FStd;
using ARM::ARM_ModelParamsHW1FExt;
using ARM::ARM_HullWhite2F;
using ARM::ARM_ModelParamsHW2F;
using ARM::ARM_ModelParamsHW2FStd;
using ARM::ARM_ModelParamsHW2FExt;
using ARM::ARM_VectorPtr;
using ARM::ARM_Numeraire;
using ARM::ARM_NumeraireFactory;
using ARM::ARM_SFRM;
using ARM::ARM_VanillaSwaptionArgSFRM;
using ARM::ARM_AutoCleaner;
using ARM::ARM_ModelParamsSFRM;
using ARM::ARM_ModelParamsSFRMFactory;
using ARM::ARM_ModelParamsSmiled;
using ARM::ARM_ModelParamsSmiled_Fx;
using ARM::ARM_SingletonHolder;
using ARM::ARM_CRMCookies;
using ARM::ARM_HybridBasisFwdIR;
using ARM::ARM_ForwardMarginIR;
using ARM::ARM_ForwardMarginBasis;
using ARM::ARM_ForwardForex;
using ARM::ARM_PricingModelPtr;
using ARM::ARM_PricingModelIRPtr;
using ARM::ARM_QGM1F;
using ARM::ARM_ModelParamsQGM1F;
using ARM::ARM_QGM1F;
using ARM::ARM_QGM2F;
using ARM::ARM_ModelParamsQGM2F;
using ARM::ARM_ModelParamsQGM2FOneDim;
using ARM::ARM_PricingStates;
using ARM::ARM_PricingStatesPtr;
using ARM::ARM_ModelParamsVec;
using ARM::ARM_ArgConv_ModelParam;
using ARM::ARM_ArgConv_BumpParam;
using ARM::ARM_BumpParamType;
using ARM::ARM_EqFxBase;
using ARM::GaussLegendre_Coefficients;
using ARM::ARM_ModelParamType;
using ARM::ARM_DataType;
using ARM::ARM_ParamType;
using ARM::ARM_CurveModelParam;
using ARM::ARM_Q1FAna_ModelParams;
using ARM::ARM_Q1FAna_Model;
using ARM::ARM_ModelParamsQ1F;
using ARM::ARM_QModel1F;
using ARM::ARM_GP_Matrix;
using ARM::ARM_ModelParams;
using ARM::ARM_SurfaceModelParam;
using ARM::ARM_Surface;
using ARM::ARM_Heston_Model;
using ARM::ARM_Heston_ModelParams;
using ARM::ARM_ShiftedHeston_Model;
using ARM::ARM_ShiftedHeston_ModelParams;
using ARM::ARM_BS_Model;
using ARM::ARM_BS_ModelParams;
using ARM::ARM_CEV_Model;
using ARM::ARM_CEV_ModelParams;
using ARM::ARM_MarkovFunctional;
using ARM::ARM_ModelParamsMF;
using ARM::ARM_Merton_Model;
using ARM::ARM_Merton_ModelParams;
using ARM::ARM_Normal_Model;
using ARM::ARM_Normal_ModelParams;
using ARM::ARM_SABR_Model;
using ARM::ARM_SABR_ModelParams;
using ARM::ARM_SLN_Model;
using ARM::ARM_SLN_ModelParams;
using ARM::ARM_BoolVector;
using ARM::std::vector<double>;
using ARM::ARM_GP_VectorPtr;
using ARM::ARM_ModelParamsQNF;
using ARM::ARM_QModelNF;
using ARM::ARM_CorrelMatParam;
using ARM::ARM_SurfacePtr;
using ARM::ARM_SurfacePtrVector;
using ARM::ARM_SurfaceListModelParam;
using ARM::ARM_Curve;
using ARM::ARM_CurveMatrix;
using ARM::ARM_InfCurv;
using ARM::ARM_GP_CurvePtr;
using ARM::ARM_ModelNameMap;
using ARM::ARM_StringVector;
using ARM::ARM_StringVectorVector;
using ARM::ARM_EqFx_ModelFactory;
using ARM::CreateClonedPtr;
using ARM::CreateClone;
using ARM::ARM_ZeroCurvePtr;
using ARM::ARM_MultiAssetsFactory;
using ARM::ARM_MultiAssetsModel;
using ARM::ARM_2IRFXModel;
using ARM::ARM_1IRFXModel;
using ARM::ARM_HWHWQtoModel;
using ARM::ARM_HWHW2FQtoModel;
using ARM::ARM_ModelParamsSABR_Eq;
using ARM::ARM_NumerairePtr;
using ARM::ARM_SABR_Eq;
using ARM::ARM_ModelParams_Eq;
using ARM::ARM_USERNAME;
using ARM::ARM_Local_Model;
using ARM::ARM_Local_Normal_Model;
using ARM::ARM_Local_Normal_ModelParams;
using ARM::ARM_Local_SLN_Model;
using ARM::ARM_Local_SLN_ModelParams;
using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_MarketIRModel;
using ARM::ARM_ConverterFromKernel;
using ARM::ARM_VanillaSwaptionArg;
using ARM::ARM_SmiledFRM;
using ARM::ARM_SmiledFRMfactory;
using ARM::ARM_SmiledFRMfactoryImp;
using ARM::ARM_SmiledMM;
using ARM::ARM_ArgConv_MMCalibProxy;
using ARM::ARM_ArgConv_MMCalibPattern;
using ARM::ARM_ArgConv_MMCorrelType;
using ARM::ARM_SmiledModel_Fx;
using ARM::ARM_MixtureModel_Fx;
using ARM::ARM_ParamsMixture_Fx;
using ARM::ARM_DensityFunctor;
using ARM::ARM_SVBGM;
using ARM::ARM_ModelParamsSVBGM;
using ARM::ARM_BGMSV1F;
using ARM::ARM_ModelParamsBGMSV1F;
using ARM::ARM_ParamsMixture_Fx;
using ARM::ARM_BiSVMM;
using ARM::ARM_EQHWSV_ModelParams;
using ARM::ARM_EQHWSV_NumMethods;
using ARM::ARM_EQHWSV;
using ARM::ARM_SVMMSpread;
using ARM::ARM_NP1IRNFXModel;
using ARM::ARM_HWxSVMMSpread;
using ARM::ARM_HWSBGMQtoModel;
using ARM::ARM_LN_Fx;
using ARM::ARM_HWSVBGMQtoModel;
using ARM::StrUpper;
using ARM::ARM_BGMSV2F;
using ARM::ARM_ModelParamsBGMSV2F;
using ARM::ARM_ArgConv_HestonMCScheme;
using ARM::ARM_2IRFXSV;
using ARM::ARM_GenCalculator;

/// typedef for modelType
typedef ARM::ARM_EqFx_ModelFactoryImp::ModelType ModelType;

/// STL
using std::string;
using std::vector;
using ARM::ARM_SurfaceListModelParam;

////////////////////////////////////////////
//// Function to create an ir forward model
////////////////////////////////////////////
extern long ARMLOCAL_IRFwd_Create(
	long zeroCurveId,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_IrFwdMod* mod= NULL;

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		};

		mod = new ARM_IrFwdMod( CreateClonedPtr( curve) );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an inf forward model
////////////////////////////////////////////
extern long ARMLOCAL_InfFwd_Create(
	long infCurveId,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_InfFwdMod* mod= NULL;

	try
	{
		ARM_InfCurv* InfCurve = NULL;

		if( !GetObjectFromId( &InfCurve, infCurveId, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: Inflation fwd curve is not of a good type");
			return ARM_KO;
		};

		mod = new ARM_InfFwdMod( CreateClonedPtr( InfCurve ) );
		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an inflation equity model
////////////////////////////////////////////
extern long ARMLOCAL_InflationEquityModel_Create(
	long infCurveId,
	double publicationLag,
	long param1Id,
	long param2Id,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_InflationEquityMod* mod= NULL;

	try
	{
		ARM_InfCurv* InfCurve = NULL;

		if( !GetObjectFromId( &InfCurve, infCurveId, ARM_INFCURV ) )
		{
			result.setMsg ("ARM_ERR: Inflation fwd curve is not of a good type");
			return ARM_KO;
		};

		//// Checks Model Parameters

        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);

        modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

		mod = new ARM_InflationEquityMod( CreateClonedPtr(InfCurve),ARM_PricingModelIRPtr(NULL), publicationLag);
		mod->SetModelParams(ARM_ModelParamsInflationEquity(modelParams));

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a Model Parameter
////////////////////////////////////////////
extern long ARMLOCAL_ModelParam_Create(
    long                    modelParamType,
    const VECTOR<double>&   paramTimes,
    const VECTOR<double>&   paramValues,
    const CCString&         modelParamName,
    const VECTOR<double>&   C_LowerBoundary,
    const VECTOR<double>&   C_UpperBoundary,
    const CCString&         C_InterpolMethodName,
	bool					C_AdviseBreakPointTimes,
	const CCString&			C_Currency,
    ARM_result&	            result, 
	long                    objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ModelParam* modelParam= NULL;
    std::vector<double>& times = NULL; 
    std::vector<double>& values = NULL; 

    std::vector<double>& lowerBound = NULL;
    std::vector<double>& upperBound = NULL;
    
	try
	{ 
        times   = paramTimes.empty() ?  NULL : new std::vector<double>(paramTimes);
		ARM_AutoCleaner<std::vector<double>> HoldTimes(times);
        values  = paramValues.empty() ?  NULL : new std::vector<double>(paramValues);
		ARM_AutoCleaner<std::vector<double>> HoldValues(values);


        lowerBound  = C_LowerBoundary.empty() ?  NULL : new std::vector<double>(C_LowerBoundary);
		ARM_AutoCleaner<std::vector<double>> HoldLowerBound(lowerBound);
        upperBound  = C_UpperBoundary.empty() ?  NULL : new std::vector<double>(C_UpperBoundary);
		ARM_AutoCleaner<std::vector<double>> HoldUpperBound(upperBound);

        modelParam = ARM_ModelParamFactory.Instance()->CreateModelParam(
			modelParamType,
            values,
            times,
            CCSTringToSTLString(modelParamName),
            CCSTringToSTLString(C_InterpolMethodName),
            lowerBound,
            upperBound,
			C_AdviseBreakPointTimes,
			CCSTringToSTLString(C_Currency) );

		/// assign object
		if( !assignObject( modelParam, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an SABR equity model
////////////////////////////////////////////
extern long ARMLOCAL_SABREquityModel_Create(
	const long& ZcCurveId,
	const double& Spot,
	const vector<long>& paramsIds,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_SABR_Eq* mod= NULL;

	try
	{
		ARM_ZeroCurve* ZcCurve = NULL;
		if( !GetObjectFromId( &ZcCurve, ZcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		//// Checks Model Parameters
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(6);
		if( paramsIds.size() != modelParams.size() )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": paramsIds.size() != modelParams.size() !");

		for( size_t i=0; i<modelParams.size(); ++i )
		{
		    if( !GetObjectFromId( &(modelParams[i]), paramsIds[i], ARM_MODELPARAM ) )
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << "ARM_ERR: model parameter " << i << " is not of a good type";
	            ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
			}
		}


		ARM_ModelParamsSABR_Eq modelparams_Eq(modelParams, CreateClonedPtr(ZcCurve), Spot );
		mod = new ARM_SABR_Eq( CreateClonedPtr(ZcCurve), &modelparams_Eq);

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a Cst Model Parameter
////////////////////////////////////////////
extern long ARMLOCAL_CstModelParam_Create(
    long                    modelParamType,
    double					paramValue,
	bool					adviseBreakPointTimes,
	ARM_result&             result, 
	long                    objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelParam* modelParam	= NULL;
    
	try
	{ 
        modelParam = new ARM_CurveModelParam(
			(ARM_ParamType)modelParamType,
            paramValue);

		/// assign object
		if( !assignObject( modelParam, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a Trigo Correl Parameter
////////////////////////////////////////////
extern long ARMLOCAL_TrigoCorrelParam_Create(
	double asOfDate,
    long dateStripId,
	double theta,
	const CCString& interpolatorName,
    const vector<double>& lowerBound,
    const vector<double>& upperBound,
	ARM_result& result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelParam* modelParam	= NULL;
    
	try
	{ 
		ARM_DateStrip* dateStrip= NULL;
		std::vector<double>& datesVector=NULL;
		std::vector<double>& julianDatesVector=NULL;

		if( !GetObjectFromId( &dateStrip, dateStripId, ARM_DATESTRIP ) )
		{
			dateStrip = NULL;
			if (!(datesVector = dynamic_cast<std::vector<double> *>(LOCAL_PERSISTENT_OBJECTS->GetObject(dateStripId))))
			{
				result.setMsg ("ARM_ERR: date strip is not of a good type");
				return ARM_KO;
			}
		};

		double julianAsOfDate = XLDateToJulian(asOfDate);

		if (datesVector)
		{
			julianDatesVector = new std::vector<double>(datesVector->size());
			for (size_t i = 0; i < datesVector->size(); ++i)
				(*julianDatesVector)[i] = XLDateToJulian((*datesVector)[i]);
		}

		std::vector<double>& lowerBoundVec = lowerBound.empty()? NULL : new std::vector<double>( lowerBound );
		std::vector<double>& upperBoundVec = upperBound.empty()? NULL : new std::vector<double>( upperBound );

		if (dateStrip && !datesVector)
		{
			modelParam = ARM_ModelParamFactory.Instance()->CreateTrigoCorrelMatParam(
				theta, 
				*dateStrip,
				julianAsOfDate,
				CCSTringToSTLString(interpolatorName),
				/// ugly const cast!
				lowerBoundVec,
				upperBoundVec );
		}
		else if (!dateStrip && datesVector)
		{
			modelParam = ARM_ModelParamFactory.Instance()->CreateTrigoCorrelMatParam(
				theta, 
				*julianDatesVector,
				julianAsOfDate,
				CCSTringToSTLString(interpolatorName),
				/// ugly const cast!
				lowerBoundVec,
				upperBoundVec );
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Either a datestrip or a dates vector.");

		/// assign object
		if( !assignObject( modelParam, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////
/// Function to return the model param from a model
////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelParam(
	long modelId,
    long modelParamType,
	long dataType,
	double index,
	long factorNb,
	VECTOR<double>& data,
	long& rows,
	long& cols,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_PricingModel* model= NULL;
		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};

		std::vector<double>& dataVec= NULL;
		ARM_ModelParam& modelParam =  model->GetModelParams()->GetModelParam( (ARM_ParamType) modelParamType, factorNb );
		
		if( ARM_CurveModelParam* curveModelParam = dynamic_cast<ARM_CurveModelParam*>(&modelParam) )
		{
			if( ARM_CorrelMatParam* correlMatParam = dynamic_cast<ARM_CorrelMatParam*>(&modelParam) )
			{
				dataVec = correlMatParam->GetData( (ARM_DataType) dataType , rows, cols);
				data.resize(dataVec->size() );
				std::copy( dataVec->begin(), dataVec->end(), data.begin() );
			}
			else
			{
				dataVec = curveModelParam->GetData( (ARM_DataType) dataType );
				data.resize(dataVec->size() );
				std::copy( dataVec->begin(), dataVec->end(), data.begin() );
				rows = dataVec->size();
				cols = 1;
			}
		}
		else if( ARM_SurfaceModelParam* surfaceModelParam = dynamic_cast<ARM_SurfaceModelParam*>(&modelParam) )
		{
			dataVec = surfaceModelParam->GetData( (ARM_DataType) dataType, rows, cols );
			data.resize(dataVec->size() );
			std::copy( dataVec->begin(), dataVec->end(), data.begin() );
		}
		else if( ARM_SurfaceListModelParam* surfaceListModelParam = dynamic_cast<ARM_SurfaceListModelParam*>(&modelParam) )
		{
			int taille = surfaceListModelParam->GetSurfaceList().size();
			/*
			if (index > taille-1)
			{
				result.setMsg ("ARM_ERR: index (starting at 0) is greater than the size of the surface list");
				return ARM_KO;
			};
			*/
			
			dataVec = surfaceListModelParam->GetData( (ARM_DataType) dataType, index, rows, cols );
			data.resize(dataVec->size());
			std::copy( dataVec->begin(), dataVec->end(), data.begin());

		}
		else
		{
			result.setMsg ("ARM_ERR: unknown type for model param: accepted curve and surface");
			return ARM_KO;
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////
/// Function to return the model param object  from a model
////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelParamId(
	const long& modelId,
    const long& modelParamType,
	const long& factorNb,
	ARM_result& result ,
	long                    objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_PricingModel* model= NULL;
		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};
		ARM_ModelParam* modelParam =  static_cast<ARM_ModelParam*> (model->GetModelParams()->GetModelParam( (ARM_ParamType) modelParamType, factorNb ).Clone());
		/// assign object
		if( !assignObject( modelParam, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////
/// Function to Get a modelNameMap from an Hybrid Model
////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelMap(
	const long& modelId,
	ARM_result& result ,
	long                    objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelNameMap* modelNameMap = NULL;

	try
	{
		ARM_PricingModel* model= NULL;

		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};

		if (ARM_MultiAssetsModel* modelMultiAssets = dynamic_cast<ARM_MultiAssetsModel*> (model))
		{
			modelNameMap = static_cast<ARM_ModelNameMap*>(modelMultiAssets->GetModelMap()->Clone());
		}
		else
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}
		/// assign object
		if( !assignObject( modelNameMap, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}
	catch(Exception& x)
	{
		delete modelNameMap;
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////
/// Function to Set a modelNameMap for an Hybrid Model
////////////////////////////////////
extern long ARMLOCAL_PricingModel_SetModelMap(
	long modelId,
	long modelMapId,
	ARM_result& result ,
	long                    objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_MultiAssetsModel* newModel;

	try
	{
		ARM_PricingModel* model= NULL;
		//Get Model
		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};
		
		// Get ModelMap
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelMapId);
		ARM_ModelNameMap*  modelMap = dynamic_cast<ARM_ModelNameMap*>(armObj);

		if( !modelMap)
		{
			result.setMsg ("ARM_ERR: model map is not of a good type");
			return ARM_KO;
		};

		newModel= dynamic_cast<ARM_MultiAssetsModel*>( model->Clone() ) ;
		if(	newModel )
		{
			newModel->SetModelMap(modelMap);
		}
		else
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an HW1F model
////////////////////////////////////////////
extern long ARMLOCAL_HW1FModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	vector<double>& flags,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_HullWhite1F* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "HW1F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
        modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
        
		ARM_BoolVector soFlags(flags.size());
		for(size_t i=0;i<flags.size();++i)
			soFlags[i] = flags[i] != 0.0;
		mod = new ARM_HullWhite1F( CreateClonedPtr( curve ), NULL, soFlags );
        
		ARM_CurveModelParam* meanReversionParam;
        if(modelParams[0]->GetType() == ARM_ModelParamType::MeanReversion)
            meanReversionParam=(ARM_CurveModelParam*) modelParams[0];
        else if(modelParams[1]->GetType() == ARM_ModelParamType::MeanReversion)
            meanReversionParam=(ARM_CurveModelParam*) modelParams[1];
        else
        {
            delete mod;
		    result.setMsg ("ARM_ERR: mean reversion parameter is missing");
		    return ARM_KO;
        }

        if(meanReversionParam->GetCurve()->GetAbscisses().size() <= 1)
            mod->SetModelParams(ARM_ModelParamsHW1FStd(modelParams));
        else if(meanReversionParam->GetCurve()->GetAbscisses().size() > 1)
            mod->SetModelParams(ARM_ModelParamsHW1FExt(modelParams));
        else
        {
            delete mod;
		    result.setMsg ("ARM_ERR: mean reversion parameter is not consistent");
		    return ARM_KO;
        }

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an MF1F model
////////////////////////////////////////////
extern long ARMLOCAL_MarkovFunctionalModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long		param2Id,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_MarkovFunctional* mod= NULL;

	try
	{

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
        modelParams[0] = NULL;
		modelParams[1] = NULL;

	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		

		/// param2Id provided
		if( param2Id != ARM_NULL_OBJECT )
		{		
			if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
		}
		/// param2Id not provided --> default = 0 mean reversion
		else
		{
			modelParams[1] = new ARM_CurveModelParam ( ARM_ModelParamType::MeanReversion, 0.0);
		}

		/// Create Markov Functional
		mod = new ARM_MarkovFunctional( CreateClonedPtr( curve ) );
		/// Sets model params
        mod->SetModelParams(ARM_ModelParamsMF(modelParams));

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an MSV1F model
////////////////////////////////////////////
extern long ARMLOCAL_MSV1FModel_Create(
	const long&        zeroCurveId,
	const vector<long>&		modelParamVec,
	const double&       forwardTerm,
	const long &        IRIndexId,
	const string&		fwdType,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(6);
	ARM_ModelParamsMSV1F* MSVModelParams = NULL;
	ARM_MarkovSV1F* mod= NULL;

	try
	{
		bool isSwapRate = true;
		if (fwdType == "N" || fwdType == "FALSE")
			isSwapRate = false;


		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_IRIndex* IRIndex= NULL;
		if( !GetObjectFromId( &IRIndex, IRIndexId, ARM_IRINDEX ) )
		{
			if (isSwapRate)
			{
				result.setMsg ("ARM_ERR: Need IRIndex for Swaptions");
				return ARM_KO;
			}
		}
        
		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: MSV model expected 6 model params!");
			return ARM_KO;
		}
		
		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

		MSVModelParams = new ARM_ModelParamsMSV1F(modelParams,IRIndex);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_ModelParamsMSV1F> Hold(MSVModelParams);
		mod = new ARM_MarkovSV1F( CreateClonedPtr( curve ), MSVModelParams ,forwardTerm,isSwapRate);

		/// Sets model params
        //mod->SetModelParams(ARM_ModelParamsMSV1F(modelParams));

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an FRMSV model
////////////////////////////////////////////
extern long ARMLOCAL_FRMSVModel_Create(
	const long&        zeroCurveId,
	const vector<long>&		modelParamVec,
	const vector<long>&		modelParamVec2,
	const long &        IRIndexId,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(3);
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams2(4);
	CC_STL_VECTOR( ARM_ModelParams* ) modelParams3(2);
	ARM_ModelParamsFRMSV1F* FRMSVModelParams = NULL;
	ARM_FRMSV1F* mod= NULL;

	try
	{
		int row_diag = 1;

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_IRIndex* IRIndex= NULL;
		if( !GetObjectFromId( &IRIndex, IRIndexId, ARM_IRINDEX ) )
		{
			result.setMsg ("ARM_ERR: Need IRIndex for Swaptions");
			return ARM_KO;
		}
        
		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: FRM SV model expected 3 model params!");
			return ARM_KO;
		}
		
		if( modelParamVec2.size() != modelParams2.size() )
		{
			result.setMsg ("ARM_ERR: FRM SV model expected 4 model params!");
			return ARM_KO;
		}
        size_t i;
		for( i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

		
		for(i=0; i<modelParams2.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec2[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams2[i] = modelParam;
		}


		ARM_ModelParamsSFRM* sfrmParams = ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(modelParams,&*IRIndex,1,row_diag);
		ARM_VolParams* volParams = new ARM_VolParams(modelParams2);

		modelParams3[0] = sfrmParams;
		modelParams3[1] = volParams;
		
		FRMSVModelParams = new ARM_ModelParamsFRMSV1F(modelParams3);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_ModelParamsFRMSV1F> Hold(FRMSVModelParams);
		mod = new ARM_FRMSV1F( CreateClonedPtr( curve ), FRMSVModelParams);

		/// Sets model params
        mod->SetModelParams(ARM_ModelParamsFRMSV1F(modelParams3));

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an HWSV1F model
////////////////////////////////////////////
extern long ARMLOCAL_HWSV1FModel_Create(
	const long&				zeroCurveId,
	const vector<long>&		modelParamVec,
	const long&				solverType,
	const vector<double>&	solverParams,
	const long&				formulaType,
	const vector<double>&	formulaParams,
	const long&				formulaTypeSO,
	const vector<double>&	formulaParamsSO,
	const double&			maxDecay,
	const double&			maxDecaySO,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	CC_STL_VECTOR( ARM_ModelParam* ) modelParams;
	ARM_ModelParamsHWSV* HWSVModelParams = NULL;
	ARM_HWSV1F* mod= NULL;

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
		if( modelParamVec.size() < 3 || modelParamVec.size() > 6 )
		{
			result.setMsg ("ARM_ERR: HWSV model expects 4 to 6 model params!");
			return ARM_KO;
		}
		
		for( size_t i=0; i<modelParamVec.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams.push_back(modelParam);
		}

		HWSVModelParams = new ARM_ModelParamsHWSV(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_ModelParamsHWSV> Hold(HWSVModelParams);
		std::vector<double> solverDatas(solverParams);
		std::vector<double> formulaDatas(formulaParams);
		std::vector<double> formulaDatasSO(formulaParamsSO);
		mod = new ARM_HWSV1F( CreateClonedPtr( curve ), HWSVModelParams, solverType, solverDatas, formulaType, formulaDatas, formulaTypeSO, formulaDatasSO, maxDecay, maxDecaySO);


		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an HWSV2F model
////////////////////////////////////////////
extern long ARMLOCAL_HWSV2FModel_Create(
	const long&				zeroCurveId,
	const vector<long>&		modelParamVec,
	const vector<double>&	solverParams,
	const long&				formulaType,
	const vector<double>&	formulaParams,
	const long&				formulaTypeSO,
	const vector<double>&	formulaParamsSO,
	const double&			maxDecay,
	const double&			maxDecaySO,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	CC_STL_VECTOR( ARM_ModelParam* ) modelParams;
	ARM_ModelParamsHWSV* HWSVModelParams = NULL;
	ARM_HWSV2F* mod= NULL;

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
		if( modelParamVec.size() < 7 || modelParamVec.size() > 8 )
		{
			result.setMsg ("ARM_ERR: HWSV model expects 7 or 8 model params!");
			return ARM_KO;
		}
		
		for( size_t i=0; i<modelParamVec.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams.push_back(modelParam);
		}

		HWSVModelParams = new ARM_ModelParamsHWSV(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_ModelParamsHWSV> Hold(HWSVModelParams);
		std::vector<double> solverDatas(solverParams);
		std::vector<double> formulaDatas(formulaParams);
		std::vector<double> formulaDatasSO(formulaParamsSO);
		mod = new ARM_HWSV2F( CreateClonedPtr( curve ), HWSVModelParams, solverDatas, formulaType, formulaDatas, formulaTypeSO, formulaDatasSO, maxDecay, maxDecaySO);

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an HW2F model
////////////////////////////////////////////
extern long ARMLOCAL_HW2FModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	long        param3Id,
	long        param4Id,
	long        param5Id,
	vector<double>& flags,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_HullWhite2F* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "HW2F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(5);

        modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[2] = NULL;
	    if( !GetObjectFromId( &(modelParams[2]), param3Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[3] = NULL;
	    if( !GetObjectFromId( &(modelParams[3]), param4Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[4] = NULL;
	    if( !GetObjectFromId( &(modelParams[4]), param5Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }


		ARM_BoolVector soFlags(flags.size());
		for(size_t i=0;i<flags.size();++i)
			soFlags[i] = flags[i] != 0.0;
		mod = new ARM_HullWhite2F( CreateClonedPtr( curve ), NULL, soFlags );

		size_t k;
		ARM_CurveModelParam* volratioParam = 0;
		ARM_CurveModelParam* correlationParam = 0;
		for (k=0;k<5;k++)
		{
			if(modelParams[k]->GetType() == ARM_ModelParamType::VolatilityRatio)
				volratioParam=(ARM_CurveModelParam*) modelParams[k];
			if(modelParams[k]->GetType() == ARM_ModelParamType::Correlation)
				correlationParam=(ARM_CurveModelParam*) modelParams[k];
		}

        if ( volratioParam == 0 ||  correlationParam ==0)
		{
			delete mod;
		    result.setMsg ("ARM_ERR: volatility ratio or correlation is missing");
		    return ARM_KO;
        }

        if(volratioParam->GetCurve()->GetAbscisses().size() <= 1 && correlationParam->GetCurve()->GetAbscisses().size() <= 1 )
            mod->SetModelParams(ARM_ModelParamsHW2FStd(modelParams));
        else if(volratioParam->GetCurve()->GetAbscisses().size() > 1 || correlationParam->GetCurve()->GetAbscisses().size() > 1)
            mod->SetModelParams(ARM_ModelParamsHW2FExt(modelParams));
        else
        {
            delete mod;
		    result.setMsg ("ARM_ERR: volatility ratio or correlation is not consistent");
		    return ARM_KO;
        }

		/// creates the default numeraire Terminal ZC for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create an EQHWSV_ParamModel
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_ModelParamsCreate(
	const vector<long>&		modelParamVec,
	ARM_result&				result, 
	long					objId )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	ARM_EQHWSV_ModelParams*	modParam=NULL;

	try	{
		
		CC_STL_VECTOR(ARM_ModelParam* )	vecModParam;
		for( size_t i=0; i<modelParamVec.size(); ++i ){
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);

			if( !modelParam ){
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			vecModParam.push_back(modelParam);
		}

		modParam	= new ARM_EQHWSV_ModelParams( vecModParam );

		if( !assignObject( modParam, result, objId ) )	return ARM_KO; 
		else											return ARM_OK; 
	}
	
	catch(Exception& x){
		delete modParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to create an ARM_EQHWSV_NumMethod
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_NumMethodsCreate(
	const long&				intStep,
	const double&			imAxis,
	const double&			maxDecay,
	ARM_result&				result, 
	long					objId )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	ARM_EQHWSV_NumMethods*	numParam=NULL;
	
	try	{

	ARM_EQHWSV_NumMethods*	numParam= new ARM_EQHWSV_NumMethods(maxDecay,imAxis,intStep);

	if( !assignObject( numParam, result, objId ) )	return ARM_KO; 
	else											return ARM_OK; 

	}
	catch(Exception& x)	{
		delete numParam;
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to create an ARM_EQHWSV
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_Create(
	const long&				zeroCurveId,
	const long&				modelParamsId,
	const long&				numMethodsId,
	const double&			dilatation,
	ARM_result&				result, 
	long					objId )
{
	if( !GlobalPersistanceOk( result ) )	return ARM_KO;
	CCString msg ("");
	ARM_EQHWSV*	model=NULL;
	
	try	{

		ARM_Object* obj1	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamsId );
		ARM_EQHWSV_ModelParams* modParam = dynamic_cast<ARM_EQHWSV_ModelParams*>(obj1);
		if( !modParam ){
			result.setMsg ("ARM_ERR: model parameter is not of a good type");
			return ARM_KO;
		}

		ARM_Object* obj2	=  LOCAL_PERSISTENT_OBJECTS->GetObject(numMethodsId );
		ARM_EQHWSV_NumMethods* numMethod = dynamic_cast<ARM_EQHWSV_NumMethods*>(obj2);
		if( !numMethod ){
			result.setMsg ("ARM_ERR: numeric method is not of a good type");
			return ARM_KO;
		}

		ARM_Object* obj3	=  LOCAL_PERSISTENT_OBJECTS->GetObject(zeroCurveId );
		ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(obj3);
		if( !zcCurve ){
			result.setMsg ("ARM_ERR: numeric method is not of a good type");
			return ARM_KO;
		}


	ARM_EQHWSV*	model= new ARM_EQHWSV( zcCurve, modParam, numMethod, dilatation) ;

	if( !assignObject( model, result, objId ) )	return ARM_KO; 
	else										return ARM_OK; 

	}
	catch(Exception& x)	{
		delete model;
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to create an SFRM model
////////////////////////////////////////////
extern long ARMLOCAL_SFRMModel_Create(
	const long&		zeroCurveId,	/// interest rate curve
	const VECTOR<long >& paramsIdVec,/// volatility,mean reversion,shift,correlation
	const long&		volType,		/// DIAG or Row
	const double&	factorsNb,		/// nb of factors
	const long&     IRIndexId,		/// irindex
	const long&		shiftConvPortId,/// portfolio id for beta conversion to shift
	const bool&		nonParamDrift,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "SFRM" );

		ARM_SFRM* mod= NULL;

        ARM_Portfolio* shiftConvPort = NULL;                            
	    if( !GetObjectFromIdwNull( &shiftConvPort, shiftConvPortId, ARM_PORTFOLIO) )
	    {
		    result.setMsg ("ARM_ERR: shift conversion portefolio is not of a good type");
		    return ARM_KO;
	    }

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
    
		const size_t paramSize=paramsIdVec.size();
		CC_STL_VECTOR( ARM_ModelParam* ) modelParams(paramSize);

		size_t i;
		for(i=0; i<paramSize; ++i )
		{
			modelParams[i] = NULL;
			if( !GetObjectFromId( &(modelParams[i]), paramsIdVec[i], ARM_MODELPARAM ) )
			{
				char msg[150];
				sprintf( msg, "ARM_ERR: model parameter %d is not of a good type", i );
				result.setMsg( msg );
				return ARM_KO;
			}
		}

		ARM_IRIndex* IRIndex= NULL;
		if( !GetObjectFromId( &IRIndex, IRIndexId, ARM_IRINDEX ) )
		{
			result.setMsg ("ARM_ERR: irindex is not of a good type");
			return ARM_KO;
		}

		/// use the factory class to get the model params!
		ARM_ModelParamsSFRM* pSFRMModelParams = 
			ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(modelParams,IRIndex,factorsNb,volType);

		/// use auto_ptr for exception safety!
		CC_NS(std,auto_ptr)<ARM_ModelParamsSFRM> SFRMmodelParams(pSFRMModelParams);
		mod = new ARM_SFRM(CreateClonedPtr( curve ),*SFRMmodelParams,shiftConvPort, nonParamDrift);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SFRM> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function To Dump the Vol Swap / Vol FRA 
//// infos of a swaption
////////////////////////////////////////////
extern long ARMLOCAL_SFRMModel_VolSwapVolFRADump(
	const long&		swaptionId,	/// Swaption
	const long&		sfrmId,	/// SFMR Model
	VECTOR<double>& OuputMatrix,
	long& rowsOutput,
	long& colsOutput,
	ARM_result&		result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "SFRM" );

		ARM_SFRM* mod= NULL;

        ARM_Swaption* swaption = NULL;                            
	    if( !(swaption=dynamic_cast<ARM_Swaption*>(LOCAL_PERSISTENT_OBJECTS->GetObject(swaptionId))))
	    {
		    result.setMsg ("ARM_ERR: swaption is not of a good type");
		    return ARM_KO;
	    }

		ARM_SFRM* sfrm = NULL;
		if( !(sfrm=dynamic_cast<ARM_SFRM*>(LOCAL_PERSISTENT_OBJECTS->GetObject(sfrmId))))
		{
			result.setMsg ("ARM_ERR:sfrm model is not of a good type");
			return ARM_KO;
		}

		double asOf = sfrm->GetZeroCurve()->GetAsOfDate().GetJulian();

		ARM_VanillaSwaptionArg* swoptArg = static_cast<ARM_VanillaSwaptionArg*>(ARM_ConverterFromKernel::ConvertSecuritytoArgObject(swaption, asOf));

		ARM_ModelParamsSFRM* sfrmModelParams = static_cast<ARM_ModelParamsSFRM*>(sfrm->GetModelParams());

		ARM_VanillaSwaptionArgSFRM sfrmSwoptArg( *swoptArg );

		ARM_VectorPtr mu = sfrmModelParams->ComputeVolSwapvolFRA(sfrmSwoptArg,*sfrm);

		ARM_VectorPtr times(new std::vector<double>);
		ARM_VectorPtr volatilies(new std::vector<double>);

		sfrmModelParams->DumpSwaptionVolsAndTimes(*swoptArg, times, volatilies);

		rowsOutput = mu->size();
		colsOutput = 3;

		OuputMatrix.resize(rowsOutput*colsOutput);

		for (size_t i = 0; i < rowsOutput; ++i)
		{
			OuputMatrix[3*i] = (*times)[i];
			OuputMatrix[3*i+1] = (*mu)[i];
			OuputMatrix[3*i+2] = (*volatilies)[i];
		}

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to set a fix scheduler on a 
//// SFRM model
////////////////////////////////////////////
extern long ARMLOCAL_SetSFRMFixScheduler(
	const long&	SFRMModelId,	/// model Id
	const double&	startDate,	/// Start Date
	const double&	endDate,	/// End Date
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_SFRM* oldMod = NULL, *mod = NULL;
		/// test that it is an object is an ARM_FRM
		ARM_Object* armObj=  LOCAL_PERSISTENT_OBJECTS->GetObject(SFRMModelId);
		oldMod = dynamic_cast<ARM_SFRM*>(armObj);
		mod = static_cast<ARM_SFRM*>(oldMod->Clone());

		/// corresponds to 1-Jan-90 in julian date format
		const double MININUMDATE = 2447893.0;

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SFRM> Hold(mod);

		double julStartDate = XLDateToJulian(startDate);
		double julEndDate = XLDateToJulian(endDate);

		// If the dates are to small we keep the automatic mode
		if (julStartDate > MININUMDATE)
			mod->SetFixStartDate(ARM_Date(julStartDate));

		if (julEndDate > MININUMDATE)
			mod->SetFixEndDate(ARM_Date(julEndDate));

		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to create an hybrid Basis Fwd IR (BFIR) model
////////////////////////////////////////////
extern long ARMLOCAL_HybridBasisFwdIRModel_Create(
            long refIRModelId,
            long zcCurveId,
            long basisZcCurveId,
            long forexId,
            const vector< string >& modelNames,
	        ARM_result&	result, 
	        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "BFIR" );

        /// Use auto-cleaners to avoid memory leaks if an error happens
		ARM_PricingModel* refModel      = NULL;
		if( !GetObjectFromId( &refModel, refIRModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model IR is not of a good type");
			return ARM_KO;
		}
    
        if( !dynamic_cast< ARM_PricingModelIR* >(refModel) )
        {
			result.setMsg ("ARM_ERR: pricing model IR is not of a good type");
			return ARM_KO;
        }

        /// No sharing possible then clone
        ARM_PricingModelPtr refIRModel( static_cast< ARM_PricingModel* >(const_cast< ARM_PricingModel* >(refModel)->Clone()) );
        string refCcy(refIRModel->GetZeroCurve()->GetCurrencyUnit()->GetCcyName());


        /// Built the IR Margin model (curve is shared)
		ARM_ZeroCurve* zcCurve = NULL;
		if( !GetObjectFromId( &zcCurve, zcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: zc curve is not of a good type");
			return ARM_KO;
		}
        ARM_PricingModelPtr irMarginModel( static_cast< ARM_PricingModel* >(new ARM_ForwardMarginIR(CreateClonedPtr( zcCurve) ) ) );
        string irMarginCcy(zcCurve->GetCurrencyUnit()->GetCcyName());


        /// Built the Basis Margin model (curve is shared)
		ARM_ZeroCurve* basisZcCurve = NULL;
		if( !GetObjectFromId( &basisZcCurve, basisZcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: basis zc curve is not of a good type");
			return ARM_KO;
		}
        ARM_PricingModelPtr basisMarginModel( static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr( basisZcCurve ) )) );
        string basisMarginCcy(basisZcCurve->GetCurrencyUnit()->GetCcyName());

        /// Get the reference model of the Basis Margin one
        ARM_HybridBasisFwdIR::modelsAlias refBasisModelIdx;
        if(modelNames.size() < ARM_HybridBasisFwdIR::NbModels)
        {
			result.setMsg ("ARM_ERR: Not enough model names");
			return ARM_KO;
        }
        if(basisMarginCcy == refCcy)
            refBasisModelIdx = ARM_HybridBasisFwdIR::RefModel;
        else if(basisMarginCcy == irMarginCcy)
            refBasisModelIdx = ARM_HybridBasisFwdIR::IrMarginModel;
        else
        {
            string err=string("ARM_ERR: basis curve currency must be ");
            err += refCcy;
            err += " or ";
            err += irMarginCcy;
			result.setMsg (err.c_str());
			return ARM_KO;
        }

        /// Get domestic and foreign models
		ARM_Forex* forex = NULL;
		if( !GetObjectFromId( &forex, forexId, ARM_FOREX ) )
		{
			result.setMsg ("ARM_ERR: forex is not of a good type");
			return ARM_KO;
		}

        string foreignCcy(forex->GetMainCurrency()->GetCcyName());
        string domesticCcy(forex->GetMoneyCurrency()->GetCcyName());
        ARM_HybridBasisFwdIR::modelsAlias domesticModelIdx,foreignModelIdx,modelIdx;
        if(domesticCcy == basisMarginCcy)
        {
            domesticModelIdx    = ARM_HybridBasisFwdIR::BasisMarginModel;
            foreignModelIdx     = (foreignCcy == irMarginCcy ? ARM_HybridBasisFwdIR::IrMarginModel
                                                             : ARM_HybridBasisFwdIR::RefModel);
            modelIdx            = foreignModelIdx;
        }
        else if(foreignCcy == basisMarginCcy)
        {
            domesticModelIdx    = (domesticCcy == irMarginCcy ? ARM_HybridBasisFwdIR::IrMarginModel
                                                              : ARM_HybridBasisFwdIR::RefModel);
            foreignModelIdx     = ARM_HybridBasisFwdIR::BasisMarginModel;
            modelIdx            = domesticModelIdx;
       }
        else
		{
			result.setMsg ("ARM_ERR: forex not related to input curves");
			return ARM_KO;
		}
 
        /// Built the Forward Forex model (curves are shared)
        ARM_PricingModelPtr forexModel( static_cast< ARM_PricingModel* >( new ARM_ForwardForex(*forex,
            modelIdx == ARM_HybridBasisFwdIR::RefModel ? refIRModel->GetZeroCurve() : CreateClonedPtr( zcCurve ), CreateClonedPtr( basisZcCurve ) ) ) );

        /// Create the model table
        vector< ARM_PricingModelPtr > modelTable(ARM_HybridBasisFwdIR::NbModels);
        modelTable[ARM_HybridBasisFwdIR::RefModel]          = refIRModel;
        modelTable[ARM_HybridBasisFwdIR::IrMarginModel]     = irMarginModel;
        modelTable[ARM_HybridBasisFwdIR::BasisMarginModel]  = basisMarginModel;
        modelTable[ARM_HybridBasisFwdIR::ForexModel]        = forexModel;


        /// Create the BFIR model
		ARM_HybridBasisFwdIR* mod = new ARM_HybridBasisFwdIR(modelNames,modelTable,refBasisModelIdx,domesticModelIdx,foreignModelIdx);


		/// Assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to Set a zero curve
////////////////////////////////////////////
extern long ARMLOCAL_Model_SetZCCurve(
	long modelId,
	long zcCurveId,
    ARM_result&	result, 
    long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	ARM_PricingModel* model = NULL;

	try
	{
		ARM_PricingModel* previousModel = NULL;
		if( !GetObjectFromId( &previousModel, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: model is not of a good type");
			return ARM_KO;
		};

		ARM_ZeroCurve* zcCurve = NULL;
		if( !GetObjectFromId( &zcCurve, zcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: zero curve is not of a good type");
			return ARM_KO;
		};

		/// clone for safety
		model = (ARM_PricingModel*) previousModel->Clone();
		model->SetZeroCurve( CreateClonedPtr( zcCurve ) );
		model = (ARM_PricingModel*) model;

		ARM_CLASS_NAME className = model->GetName();
		/// assign object
		if( !assignObject( model, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete model;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create an QGM1F model
////////////////////////////////////////////
extern long ARMLOCAL_QGM1FModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	long        param3Id,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_QGM1F* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "QGM1F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(3);
        modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
        
        modelParams[2] = NULL;
	    if( !GetObjectFromId( &(modelParams[2]), param3Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
        
		mod = new ARM_QGM1F( CreateClonedPtr(curve), ARM_ModelParamsQGM1F(modelParams) );

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_QGM1F_Test(
	long QGM1FId,
    double t,
    double T,
    double Tn,
    double Xt,
    double Ts,
    double Te,
    double K,
    int capPayOrFloorRec,
    const std::vector<double>& Tp,
    const std::vector<double>& YF,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


    ARM_QGM1F* mod=NULL;

	try
	{
		mod = dynamic_cast<ARM_QGM1F *>(LOCAL_PERSISTENT_OBJECTS->GetObject(QGM1FId));
		if (!mod)
		{
			result.setMsg ("ARM_ERR: QGM1F model is not of a good type");
			return ARM_KO;
		}

        clock_t ts,te;

        int i,NB_LOOP=1000;

        /// A(t,T, B(t,T), C(t,T)
        double AValue,BValue,CValue;

        ts=clock();
        for(i=0;i<NB_LOOP;++i) AValue = mod->A(t,T);
        te=clock();
        double ATime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;

        ts=clock();
        for(i=0;i<NB_LOOP;++i) BValue=mod->B(t,T);
        te=clock();
        double BTime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;

        ts=clock();
        for(i=0;i<NB_LOOP;++i) CValue=mod->C(t,T);
        te=clock();
        double CTime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;

        /// Mean & Variance of X(T) knowing X(t)=Xt under forward Tn probability
        ARM_PricingStatesPtr states( new ARM_PricingStates(1,1) );
        states->SetModelState(0,0,Xt);
        ARM_VectorPtr lawValues;
        ts=clock();
        for(i=0;i<NB_LOOP;++i) lawValues = mod->XLaw(t,T,Tn,states);
        te=clock();
        double meanValue = (*lawValues)[0];
        double varianceValue = (*lawValues)[1];
        double lawTime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;
		std::vector<double> strikes(1,K);

        /// Caplet price knowing X(t)=Xt paying at Tn C=max(Libor(T,Ts,Te) - K , 0).(Te-Ts)/365 in bp
        ARM_VectorPtr capletValues;
        double a = (Te-Ts)/365.0;
        ts=clock();
        for(i=0;i<NB_LOOP;++i)
            capletValues = mod->VanillaCaplet("",t,Tn,a,10000.0,T,Ts,Te,a,strikes,capPayOrFloorRec,states);
        te=clock();
        double capletPrice = (*capletValues)[0];
        double capletTime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;
        double capletVol = mod->CapletImpliedVol(Tn,a,10000.0,T,Ts,Te,a,K,capPayOrFloorRec,capletPrice);

        /// Swaption price knowing X(t)=Xt paying fixed flows at Tp[i] on period YF[i] in bp
        /// Expiry is T and floating starts at Ts and ends at Te
        ARM_GP_Matrix strikesswap(1,1,K);
		std::vector<double> Notional(1,1);
        ARM_VectorPtr swaptionValues;
        ts=clock();
        for(i=0;i<NB_LOOP;++i)
			/////// SCOTCH Tp for floatResetDates, floatStartDates , floatEndDates, floatIntTerms
            swaptionValues = mod->VanillaSwaption("",t,T,Notional,Notional,Ts,Te,Tp,Tp,Tp,Tp,Tp,YF,strikesswap,capPayOrFloorRec,states);
        te=clock();
        double swaptionPrice = (*swaptionValues)[0];
        double swaptionTime = ((double)(te-ts))/CLOCKS_PER_SEC/NB_LOOP;
        double swaptionVol = mod->SwaptionImpliedVol(T,10000.0,Ts,Te,Tp,YF,K,capPayOrFloorRec,swaptionPrice);

		result.setArray(AValue,0);
		result.setArray(ATime*1000,1);

		result.setArray(BValue,2);
		result.setArray(BTime*1000,3);

		result.setArray(CValue,4);
		result.setArray(CTime*1000,5);

		result.setArray(meanValue,6);
		result.setArray(varianceValue,7);
		result.setArray(lawTime*1000,8);

		result.setArray(capletPrice,9);
		result.setArray(capletTime*1000,10);
		result.setArray(capletVol,11);

		result.setArray(swaptionPrice,12);
		result.setArray(swaptionTime*1000,13);
		result.setArray(swaptionVol,14);

        result.setDouble(15);
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
////////////////////////////////////////////
//// Function to create an QGM2F model
////////////////////////////////////////////
extern long ARMLOCAL_QGM2FModel_Create(
	long        zeroCurveId,
	const VECTOR<long >& paramsFactor1IdVect,
	const VECTOR<long >& paramsFactor2IdVect,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_QGM2F* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "QGM2F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
		/////////
		const size_t param1Size=paramsFactor1IdVect.size();
		const size_t param2Size=paramsFactor2IdVect.size();
			
		CC_STL_VECTOR( ARM_ModelParam* ) modelParamsFactor1Vect(param1Size);
		CC_STL_VECTOR( ARM_ModelParam* ) modelParamsFactor2Vect(param2Size);

		size_t i;
		for(i=0; i<param1Size; ++i )
		{
			/*modelParamsFactor1Vect[i] = NULL;
			if( !GetObjectFromId( &(modelParamsFactor1Vect[i]), paramsFactor1IdVect[i], ARM_MODELPARAM ) )
			{
				char msg[150];
				sprintf( msg, "ARM_ERR: model parameter %d is not of a good type", i );
				result.setMsg( msg );
				return ARM_KO;
			}*/
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(paramsFactor1IdVect[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParamsFactor1Vect[i] = modelParam;
		}

		ARM_ModelParamsQGM2FOneDim* modelparams1 = new ARM_ModelParamsQGM2FOneDim(modelParamsFactor1Vect);

		for(i=0; i<param2Size; ++i )
		{
			/*modelParamsFactor2Vect[i] = NULL;
			if( !GetObjectFromId( &(modelParamsFactor2Vect[i]), paramsFactor2IdVect[i], ARM_MODELPARAM ) )
			{
				char msg[150];
				sprintf( msg, "ARM_ERR: model parameter %d is not of a good type", i );
				result.setMsg( msg );
				return ARM_KO;
			}*/
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(paramsFactor2IdVect[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParamsFactor2Vect[i] = modelParam;
		}

		ARM_ModelParamsQGM2FOneDim* modelparams2 = new ARM_ModelParamsQGM2FOneDim(modelParamsFactor2Vect);

		vector<ARM_ModelParams*> paramsVec(2);
		paramsVec[0] = modelparams1;
		paramsVec[1] = modelparams2;
		
		       
		mod = new ARM_QGM2F( CreateClonedPtr(curve), ARM_ModelParamsQGM2F(paramsVec) );

		/// creates the default numeraire for the time being
	    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a Q1FAna model
////////////////////////////////////////////
extern long ARMLOCAL_Q1FAnaModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_Q1FAna_Model* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Q1F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
        modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }

        modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), param2Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
        
		mod = new ARM_Q1FAna_Model( CreateClonedPtr( curve ),ARM_Q1FAna_ModelParams(modelParams));

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a Q1F model
////////////////////////////////////////////
extern long ARMLOCAL_Q1FModel_Create(
	const long&         zeroCurveId,
	const vector<long>& modelParamVec,
	const bool&			degenerateInHW,
	ARM_result&	result, 
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_QModel1F* mod= NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Q1F" );

		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}
        
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(3);
		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Q1F model expected 3 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}
       
		mod = new ARM_QModel1F( CreateClonedPtr( curve ),ARM_ModelParamsQ1F(modelParams),degenerateInHW);

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}
}



////////////////////////////////////////////
//// Function to create a Model Parameter (with lower anf higher bound)
////////////////////////////////////////////
extern long ARMLOCAL_SurfaceParam_Create(
    const long&		modelParamType,
    const long&		surfaceId,
	const double&	LowerBoundary,
	const double&	UpperBoundary,
    const CCString&	modelParamName,
	const bool&		adviseBreakPointTimes,
    ARM_result&	    result, 
	long            objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ModelParam* modelParam= NULL;
	ARM_Surface* surface = NULL;
    
	try
	{ 
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &surface, surfaceId ," surface", result ) ) return ARM_KO;
        modelParam = ARM_ModelParamFactory.Instance()->CreateSurfaceModelParam( 
			(ARM_ModelParamType::ParamNb) modelParamType,surface,
            CCSTringToSTLString(modelParamName), LowerBoundary, UpperBoundary, adviseBreakPointTimes);

		/// assign object
		return  ( !assignObject( modelParam, result, objId ) )  ?  ARM_KO : ARM_OK; 
	}
	
	catch(Exception& x)
	{
		delete modelParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a SurfaceList Model Parameter
////////////////////////////////////////////
long ARMLOCAL_SurfaceListParam_Create(
    const long&		modelParamType,
	const VECTOR<double>& index,
    const VECTOR<long>&		surfaceIdList,
    const CCString&	modelParamName,
    ARM_result&	    result, 

	long            objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ModelParam* modelParam= NULL;
	ARM_SurfacePtrVector surfaceList;
    
	try
	{ 
		if( index.size() != surfaceIdList.size() )
		{
			result.setMsg ("ARM_ERR: Index and Surface List should have the same size!");
			return ARM_KO;
		}

		for( size_t i=0; i<index.size(); ++i )
		{
			ARM_Object* armObj = LOCAL_PERSISTENT_OBJECTS->GetObject(surfaceIdList[i]);
			ARM_Surface* surface = dynamic_cast<ARM_Surface*>(armObj);
			if( !surface )
			{
				char msg[255];
				sprintf( msg, "object %d is not a surface", i+1 );
				result.setMsg (msg);
				return ARM_KO;
			}

			ARM_SurfacePtr surfacePtr =  ARM_SurfacePtr( static_cast<ARM_Surface*>(surface->Clone()) );
			surfaceList.push_back(surfacePtr);
		}

        modelParam = new ARM_SurfaceListModelParam( 
			(ARM_ModelParamType::ParamNb) modelParamType,
			std::vector<double>( index ),
			surfaceList,
            CCSTringToSTLString(modelParamName) );

		/// assign object
		if( !assignObject( modelParam, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelParam;
		x.DebugPrint();
		ARM_RESULT();
	}
}



//////////////////////////////////////////////////
//// Function to create a Heston Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Heston_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* HestonModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(8);
	ARM_Heston_ModelParams* HestonModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Heston model expected 8 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        HestonModelParams = new ARM_Heston_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_Heston_ModelParams> Hold(HestonModelParams);
		HestonModel = new ARM_Heston_Model( CreateClonedPtr( curve ), *HestonModelParams );

		/// assign object
		if( !assignObject( HestonModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete HestonModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a ShiftedHeston Model
//////////////////////////////////////////////////

extern long ARMLOCAL_ShiftedHeston_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
	const bool&				isMC,
	const long&				nbSteps,
	const long&				nbSimulations,
	const long&				nbIntegrationSteps,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* ShiftedHestonModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(6);
	ARM_ShiftedHeston_ModelParams* ShiftedHestonModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Shifted Heston model expected 6 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        ShiftedHestonModelParams = new ARM_ShiftedHeston_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_ShiftedHeston_ModelParams> Hold(ShiftedHestonModelParams);
		ShiftedHestonModel = new ARM_ShiftedHeston_Model( CreateClonedPtr( curve ), *ShiftedHestonModelParams,
			nbIntegrationSteps,
			isMC,
			nbSteps,
			nbSimulations);

		/// assign object
		if( !assignObject( ShiftedHestonModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete ShiftedHestonModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to create a CEV Model
//////////////////////////////////////////////////

extern long ARMLOCAL_CEV_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* CEVModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(3);
	ARM_CEV_ModelParams* CEVModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: CEV model expected 3 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        CEVModelParams = new ARM_CEV_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_CEV_ModelParams> Hold(CEVModelParams);
		CEVModel = new ARM_CEV_Model( CreateClonedPtr( curve ), *CEVModelParams );

		/// assign object
		if( !assignObject( CEVModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete CEVModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to create a BS Model
//////////////////////////////////////////////////

extern long ARMLOCAL_BS_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* BSModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
	ARM_BS_ModelParams* BSModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		modelParams.resize(modelParamVec.size());
		if( modelParamVec.size() > modelParams.size() )
		{
			result.setMsg ("ARM_ERR: BS model expected less than 2 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParamVec.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        BSModelParams = new ARM_BS_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_BS_ModelParams> Hold(BSModelParams);
		BSModel = new ARM_BS_Model( CreateClonedPtr( curve ), *BSModelParams );

		/// assign object
		if( !assignObject( BSModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete BSModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a Merton Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Merton_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* MertonModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(4);
	ARM_Merton_ModelParams* MertonModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Merton model expected 4 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        MertonModelParams = new ARM_Merton_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_Merton_ModelParams> Hold(MertonModelParams);
		MertonModel = new ARM_Merton_Model( CreateClonedPtr( curve ), *MertonModelParams );

		/// assign object
		if( !assignObject( MertonModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete MertonModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}



//////////////////////////////////////////////////
//// Function to create a Normal Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Normal_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* NormalModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(1);
	ARM_Normal_ModelParams* NormalModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Normal model expected 4 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        NormalModelParams = new ARM_Normal_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_Normal_ModelParams> Hold(NormalModelParams);
		NormalModel = new ARM_Normal_Model( CreateClonedPtr( curve ), *NormalModelParams );

		/// assign object
		if( !assignObject( NormalModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete NormalModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a SABR Model
//////////////////////////////////////////////////

extern long ARMLOCAL_SABR_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
	const long&             impliedVolType,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* SABRModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(4);
	ARM_SABR_ModelParams* SABRModelParams = NULL;
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "SABR Analytic" );

		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: SABR model expected 4 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        SABRModelParams = new ARM_SABR_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SABR_ModelParams> Hold(SABRModelParams);
		SABRModel = new ARM_SABR_Model( CreateClonedPtr( curve ), *SABRModelParams,impliedVolType );

		/// assign object
		if( !assignObject( SABRModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete SABRModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a SLN Model
//////////////////////////////////////////////////

extern long ARMLOCAL_SLN_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* SLNModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
	ARM_SLN_ModelParams* SLNModelParams = NULL;
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "SLN Analytic" );

		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: SLN model expected 2 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        SLNModelParams = new ARM_SLN_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SLN_ModelParams> Hold(SLNModelParams);
		SLNModel = new ARM_SLN_Model( CreateClonedPtr( curve ), *SLNModelParams );

		/// assign object
		if( !assignObject( SLNModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete SLNModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//// Function to create a SLN Model
//////////////////////////////////////////////////

extern long ARMLOCAL_GetVarianceSqueeze(
	const long&			modelId,
	 long&				nbRows,
     long&				nbCols,
	 vector<double>&    values,
	 bool&            status,
    ARM_result&			result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Variance Squeeze Info" );

		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelId);
		ARM_Local_Model* model	= dynamic_cast<ARM_Local_Model*>(armObj);
		if( !model)
		{
			result.setMsg ("ARM_ERR: local model is not of a good type");
			return ARM_KO;
		}

		nbRows = model->GetVarianceSqueezInfo().GetRowsNb();
        nbCols = model->GetVarianceSqueezInfo().GetColsNb();

		values =model->GetVarianceSqueezInfo().GetValues();

		status = model->GetVarianceSqueezeStatus();

		string txt("");
		txt += !status ? "No variance squeeze" : "Oups:variance squeeze";

		result.setString(txt.c_str());

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a HW1F Model Param
//////////////////////////////////////////////////

long ARMLOCAL_HW1FModelParam_Create(
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
	ARM_ModelParamsHW1F* modelParamsHW1F = NULL;
   
	try
	{ 
		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: SLN model expected 2 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        int volatilityType = ARM_ModelParamType::Volatility;
		for( i=0; i<modelParams.size(); ++i )
            if( modelParams[i]->GetType() == ARM_ModelParamType::QVol )
                volatilityType = ARM_ModelParamType::QVol;

        modelParamsHW1F = new ARM_ModelParamsHW1FStd(modelParams,volatilityType);

		/// assign object
		if( !assignObject( modelParamsHW1F, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelParamsHW1F ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a QNF Model Param
//////////////////////////////////////////////////

extern long ARMLOCAL_QNFModelParam_Create(
	const long&				qparamId,
    const vector<long>&		modelParamVec,
	const long&				correlMatId,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParamsHW1F* ) modelParams;
	ARM_ModelParamsQNF* QNFModelParams = NULL;
   
	try
	{
		ARM_Object* armObjQPAram	=  LOCAL_PERSISTENT_OBJECTS->GetObject(qparamId);
		ARM_CurveModelParam* QmodelParam = dynamic_cast<ARM_CurveModelParam*>(armObjQPAram);
		if( !QmodelParam )
		{
			result.setMsg ("ARM_ERR: Q model parameter is not of a good type");
			return ARM_KO;
		}

		ARM_Object* armObjCorrelMat	=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlMatId);
		ARM_GP_Matrix* CorrelMat	= dynamic_cast<ARM_GP_Matrix*>(armObjCorrelMat);
		if( !CorrelMat )
		{
			result.setMsg ("ARM_ERR: correl matrix is not a gp matrix");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParamVec.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParamsHW1F* modelParam = dynamic_cast<ARM_ModelParamsHW1F*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams.push_back(modelParam);
		}
		QNFModelParams = new ARM_ModelParamsQNF(*QmodelParam,modelParams,*CorrelMat);

		/// assign object
		if( !assignObject( QNFModelParams, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete QNFModelParams ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a QNF Model
//////////////////////////////////////////////////

extern long ARMLOCAL_QNFModel_Create(
	const long&				curveId,
    const long&				QNFmodelParamId,
	const bool&				degenerateInHW,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* QNFModel = NULL;
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "QNF" );

		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		ARM_Object* armObj2			=  LOCAL_PERSISTENT_OBJECTS->GetObject(QNFmodelParamId);
		ARM_ModelParamsQNF* params	= dynamic_cast<ARM_ModelParamsQNF*>(armObj2);
		if( !params)
		{
			result.setMsg ("ARM_ERR: param is not of a good type");
			return ARM_KO;
		}

        QNFModel = new ARM_QModelNF(CreateClonedPtr( curve ), *params, degenerateInHW );

		/// assign object
		if( !assignObject( QNFModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete QNFModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create an FX Model
//////////////////////////////////////////////////

extern long ARMLOCAL_FXModel_Create(
	const long&				curveId,
	const vector<long>&		modelParamVec,
	const double&			spot,
	const long&				forCurveId,
	const string&			mcScheme,
	const long&				modelType,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* Model_FX = NULL;
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXModel" );

		ARM_ZeroCurve* zccurve = NULL;       
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &zccurve, curveId ,"Doemstic Zc", result ) ) return ARM_KO;

		size_t size = modelParamVec.size();	
        CC_STL_VECTOR( ARM_ModelParam* ) modelParams(size);
		for( size_t i=0; i<modelParamVec.size(); ++i )
			if( GetObjectFromIdWithDynamicCastCheckandMsge( &modelParams[i], modelParamVec[i] ,"model parameter", result ) ) return ARM_KO;

		ARM_ZeroCurve* forCurve = NULL;       
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &forCurve, forCurveId ,"Foreign Zc", result ) ) return ARM_KO;
		
		long mcSchemeType = ARM_ArgConv_HestonMCScheme.GetNumber(mcScheme);

        Model_FX = ARM_EqFx_ModelFactory.Instance()->CreateModel( CreateClonedPtr( zccurve ), 
																	modelParams, 
																	spot,
																	CreateClonedPtr( forCurve ),
																	ARM_GP_Matrix(),
																	(ModelType) modelType,
																	mcSchemeType);

		/// assign object
		if( !assignObject( Model_FX, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete Model_FX ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create an Equity Model
//////////////////////////////////////////////////

extern long ARMLOCAL_EqModel_Create(
	const long&				curveId,
	const vector<long>&		modelParamVec,
	const double&			spot,
	const long&				modelType,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* Model_Eq = NULL;
   
	try
	{ 
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Q1FX" );

		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* zccurve = dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !zccurve )
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

        CC_STL_VECTOR( ARM_ModelParam* ) modelParams;
		modelParams.resize(0);
		for( size_t i=0; i<modelParamVec.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams.push_back(modelParam);
		}
       
        Model_Eq = ARM_EqFx_ModelFactory.Instance()->CreateModel( CreateClonedPtr( zccurve ), 
											modelParams, spot,
											ARM_ZeroCurvePtr(NULL),
											ARM_GP_Matrix(),
											(ModelType) modelType );

		/// assign object
		if( !assignObject( Model_Eq, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete Model_Eq ;
		x.DebugPrint();
		ARM_RESULT();
	}
}



//////////////////////////////////////////////////
//// Function to create a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_ModelNameMap_Create(
	const vector<string>& names,
    const vector<long>&	modelsId,
	const vector< vector< string > >& otherModelNames,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelNameMap* modelNameMap = NULL;
   
	try
	{
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "ModelNameMap" );

		if( names.size() != modelsId.size() )
		{
			result.setMsg ("names.size() != modelsId.size()");
			return ARM_KO;
		}

		CC_STL_VECTOR( ARM_PricingModelPtr ) models(modelsId.size());
		for(size_t i=0; i<modelsId.size(); ++i )
		{
			ARM_Object* armObj	= LOCAL_PERSISTENT_OBJECTS->GetObject(modelsId[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>(armObj);
		
			if( !pricingModel)
			{
				CC_Ostringstream os;
				os << "model " << i+1 << " with name " << names[i] << " is not of good type";
				result.setMsg ( os.str().c_str() );
				return ARM_KO;
			}
			models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>(pricingModel->Clone() ) );
		}

        modelNameMap = new ARM_ModelNameMap( names, models, otherModelNames );

		/// assign object
		if( !assignObject( modelNameMap, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelNameMap ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a multi assets model
//////////////////////////////////////////////////
extern long ARMLOCAL_MultiAssetsModel_Create(
	const long&				modelNameMapId,
	const long&				correlMatId,
	const string&           MultiAssetsName,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* multiAssetsModel = NULL;
   
	try
	{
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "MultiAssetsModel" );

		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelNameMapId);
		ARM_ModelNameMap* modelNameMap = dynamic_cast<ARM_ModelNameMap*>(armObj);
		if( !modelNameMap )
		{
			result.setMsg ("ARM_ERR: model name map is not of a good type");
			return ARM_KO;
		}

		ARM_GP_Matrix CorrelMat;
		if (correlMatId != ARM_NULL_OBJECT)
		{
			ARM_GP_Matrix* mat = dynamic_cast<ARM_GP_Matrix*>(LOCAL_PERSISTENT_OBJECTS->GetObject(correlMatId));
			if( !mat )
			{
				result.setMsg ("ARM_ERR: correl matrix is not a gp matrix");
				return ARM_KO;
			}
			CorrelMat = *mat;
		}
		else
		{
			/// Default
			CorrelMat = ARM_GP_Matrix(modelNameMap->size(),modelNameMap->size(),0.0);
			for(size_t i=0;i<CorrelMat.rows();++i)
				CorrelMat(i,i)=1.0;
		}

		multiAssetsModel = ARM_MultiAssetsFactory.Instance()->CreateMultiAssetsModel( *modelNameMap, CorrelMat, MultiAssetsName );

		/// assign object
		if( !assignObject( multiAssetsModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete multiAssetsModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to set a model to a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_SetModelToModelMap(
    const long&	modelMapId,
	const string& name,
	const long& modelId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelNameMap* modelNameMap = NULL;
   
	try
	{
		ARM_PricingModel* model= NULL;
		if( !GetObjectFromId( &model, modelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
		}	
		
		// Get ModelMap
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelMapId);
		ARM_ModelNameMap*  modelMap = dynamic_cast<ARM_ModelNameMap*>(armObj->Clone());

		if( !modelMap)
		{
			result.setMsg ("ARM_ERR: model map is not of a good type");
			return ARM_KO;
		};
		
		ARM_PricingModelPtr newModel( dynamic_cast<ARM_PricingModel*> (model->Clone()) );
		newModel->SetModelName(name);
		(*modelMap)[name]->Model() = newModel;
		(*modelMap)[name]->ResetOtherModels();  /// do something better when OtherModelNames will be input !
		
		/// assign object
		if( !assignObject( modelMap, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete modelNameMap;
		x.DebugPrint();
		ARM_RESULT();
	}
}



//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_GetModelFromModelMap(
    const long&	modelMapId,
	const string& name,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_ModelNameMap* modelNameMap = NULL;
	ARM_PricingModelPtr model(NULL);
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		// Get ModelMap
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelMapId);
		ARM_ModelNameMap*  modelMap = dynamic_cast<ARM_ModelNameMap*>(armObj);

		if( !modelMap)
		{
			result.setMsg ("ARM_ERR: model map is not of a good type");
			return ARM_KO;
		};
		
		newModel = static_cast<ARM_PricingModel*> ((*modelMap)[name]->Model()->Clone());
	
		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}



//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_Create2IRFXModel(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		ARM_CurveMatrix* correlCurveMatrix = dynamic_cast<ARM_CurveMatrix*>(armObj2);
		if( !correlationMatrix && !correlCurveMatrix)
		{
			result.setMsg ("ARM_ERR: correlation matrix should be a matrix or a matrix curve!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = names.size();
		ARM_StringVectorVector depends(nbModels);
		depends[ARM_2IRFXModel::FxModel]       = ARM_StringVector(2);
		depends[ARM_2IRFXModel::FxModel][0]	   = names[ARM_2IRFXModel::DomBasisModel];
		depends[ARM_2IRFXModel::FxModel][1]	   = names[ARM_2IRFXModel::ForBasisModel];
		if (nbModels > ARM_2IRFXModel::DomBasisModel)
			depends[ARM_2IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::DomModel]);
		if (nbModels > ARM_2IRFXModel::ForBasisModel)
			depends[ARM_2IRFXModel::ForBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::ForModel]);
		if(nbModels + ARM_2IRFXModel::NbLocalModels >= ARM_2IRFXModel::NbModels +1)
			depends[ARM_2IRFXModel::FlooredFxLocalModel]	= ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
		if(nbModels +ARM_2IRFXModel::NbLocalModels  >= ARM_2IRFXModel::NbModels +2)  
			depends[ARM_2IRFXModel::CappedFxLocalModel]		= ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
		if(nbModels +ARM_2IRFXModel::NbLocalModels  >= ARM_2IRFXModel::NbModels +3)  
			depends[ARM_2IRFXModel::RedemptionFxLocalModel] = ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
		ARM_ModelNameMap modelMap( names, models, depends );
		if (correlationMatrix)
			newModel = new ARM_2IRFXModel( modelMap,*correlationMatrix );
		else
			newModel = new ARM_2IRFXModel( modelMap,*correlCurveMatrix );

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_Create1IRFXModel(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const long&		Model2IRFXId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = names.size();
		ARM_StringVectorVector depends(nbModels);
		depends[ARM_1IRFXModel::FxModel]       = ARM_StringVector(1);
		depends[ARM_1IRFXModel::FxModel][0]	   = names[ARM_1IRFXModel::DomBasisModel];
		if (nbModels > ARM_1IRFXModel::DomBasisModel)
			depends[ARM_1IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_1IRFXModel::DomModel]);
		ARM_ModelNameMap modelMap( names, models, depends );
		
		ARM_2IRFXModel* Model2IRFX = NULL;
		if (Model2IRFXId != ARM_NULL_OBJECT)
		{
			if( !(Model2IRFX = (ARM_2IRFXModel*) LOCAL_PERSISTENT_OBJECTS->GetObject(Model2IRFXId)) )
			{
				result.setMsg ("ARM_ERR: 2IRFX model is not of a good type");
				return ARM_KO;
			}
		}

		newModel = new ARM_1IRFXModel( modelMap,*correlationMatrix, CreateClonedPtr(Model2IRFX));

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to create HWHWQto
//////////////////////////////////////////////////
extern long ARMLOCAL_HWHWQtoModel_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const bool& fxFlag,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		ARM_CurveMatrix* correlCurveMatrix = dynamic_cast<ARM_CurveMatrix*>(armObj2);
		if( !correlationMatrix && !correlCurveMatrix)
		{
			result.setMsg ("ARM_ERR: correlation matrix should be a matrix or a matrix curve!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = names.size();
		ARM_StringVectorVector depends(nbModels);
		depends[ARM_HWHWQtoModel::FxModel]       = ARM_StringVector(2);
		depends[ARM_HWHWQtoModel::FxModel][0]	   = names[ARM_HWHWQtoModel::DomBasisModel];
		depends[ARM_HWHWQtoModel::FxModel][1]	   = names[ARM_HWHWQtoModel::ForBasisModel];

		if (nbModels > ARM_HWHWQtoModel::DomBasisModel)
			depends[ARM_HWHWQtoModel::DomBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::DomModel]);
		if (nbModels > ARM_HWHWQtoModel::ForBasisModel)
			depends[ARM_HWHWQtoModel::ForBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);
		if(nbModels + ARM_HWHWQtoModel::NbLocalModels >= ARM_HWHWQtoModel::NbModels +1)
			depends[ARM_HWHWQtoModel::ForLocalModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);
		
		ARM_ModelNameMap modelMap( names, models, depends );
		if (correlationMatrix)
			newModel = new ARM_HWHWQtoModel( modelMap, *correlationMatrix, fxFlag );
		else
			newModel = new ARM_HWHWQtoModel( modelMap, *correlCurveMatrix, fxFlag );

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create HWHW2FQto
//////////////////////////////////////////////////
extern long ARMLOCAL_HWHW2FQtoModel_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const bool& fxFlag,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		ARM_CurveMatrix* correlCurveMatrix = dynamic_cast<ARM_CurveMatrix*>(armObj2);
		if( !correlationMatrix && !correlCurveMatrix)
		{
			result.setMsg ("ARM_ERR: correlation matrix should be a matrix or a matrix curve!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = names.size();
		ARM_StringVectorVector depends(nbModels);
		depends[ARM_HWHW2FQtoModel::FxModel]       = ARM_StringVector(2);
		depends[ARM_HWHW2FQtoModel::FxModel][0]	   = names[ARM_HWHW2FQtoModel::DomBasisModel];
		depends[ARM_HWHW2FQtoModel::FxModel][1]	   = names[ARM_HWHW2FQtoModel::ForBasisModel];

		if (nbModels > ARM_HWHW2FQtoModel::DomBasisModel)
			depends[ARM_HWHW2FQtoModel::DomBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::DomModel]);
		if (nbModels > ARM_HWHW2FQtoModel::ForBasisModel)
			depends[ARM_HWHW2FQtoModel::ForBasisModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);
		if(nbModels + ARM_HWHW2FQtoModel::NbLocalModels >= ARM_HWHWQtoModel::NbModels +1)
			depends[ARM_HWHW2FQtoModel::ForLocalModel] = ARM_StringVector(1,names[ARM_HWHWQtoModel::ForModel]);
		
		ARM_ModelNameMap modelMap( names, models, depends );
		if (correlationMatrix)
			newModel = new ARM_HWHW2FQtoModel( modelMap, *correlationMatrix, fxFlag );
		else
			newModel = new ARM_HWHW2FQtoModel( modelMap, *correlCurveMatrix, fxFlag );

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_CreateFwdMarginModel(
	const long& C_basisZcCurveId, 
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
        /// Built the IR Margin model (curve is shared)
		ARM_ZeroCurve* basisCurve = NULL;
		if( !GetObjectFromId( &basisCurve, C_basisZcCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: zc curve is not of a good type");
			return ARM_KO;
		}
        newModel = new ARM_ForwardMarginBasis(CreateClonedPtr( basisCurve) );

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to set the reference model name to a multi asset model
//////////////////////////////////////////////////

extern long ARMLOCAL_SetRefModelNameToMultiAsset(
	const string& C_Name,
	const long& C_MultiAssetModelId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_MultiAssetsModel *oldModel = NULL, *newModel = NULL ;
   
	try
	{	
		ARM_Object* armObj=  LOCAL_PERSISTENT_OBJECTS->GetObject(C_MultiAssetModelId);
		oldModel = dynamic_cast<ARM_MultiAssetsModel*>(armObj);
		if( !oldModel )
		{
			result.setMsg ("ARM_ERR: multi asset model is not of good type!");
			return ARM_KO;
		}
		newModel = (ARM_MultiAssetsModel*) oldModel->Clone();
		ARM_ModelNameMap* modelMap = newModel->GetModelMap();
		newModel->SetRefModel( &*(*modelMap)[ C_Name ]->Model() );
		
		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to create a Local Normal Model
//////////////////////////////////////////////////

extern long ARMLOCAL_LocalNormal_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Local_Model* localNormalModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams (modelParamVec.size());
	ARM_Local_Normal_ModelParams* localNormalModelParams = NULL;
   
	try
	{ 
		ARM_Object* armObj		=  LOCAL_PERSISTENT_OBJECTS->GetObject(curveId);
		ARM_ZeroCurve* curve	= dynamic_cast<ARM_ZeroCurve*>(armObj);
		if( !curve)
		{
			result.setMsg ("ARM_ERR: curve is not of a good type");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        localNormalModelParams = new ARM_Local_Normal_ModelParams(modelParams);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_Local_Normal_ModelParams> Hold(localNormalModelParams);
		localNormalModel = new ARM_Local_Normal_Model( CreateClonedPtr( curve ), *localNormalModelParams );

		/// assign object
		if( !assignObject( localNormalModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (localNormalModel) delete localNormalModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to create a Local Shifted LogNormal Model
//////////////////////////////////////////////////

extern long ARMLOCAL_LocalSLN_Model_Create(
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Local_Model* localSLNModel = NULL;
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams (3);
	ARM_Local_SLN_ModelParams* localSLNModelParams = NULL;
   
	try
	{
		if( modelParamVec.size() != modelParams.size() )
		{
			result.setMsg ("ARM_ERR: Local SLN model expected 3 model params!");
			return ARM_KO;
		}

		for( size_t i=0; i<modelParams.size(); ++i )
		{
			ARM_Object* armObj	=  LOCAL_PERSISTENT_OBJECTS->GetObject(modelParamVec[i] );
			ARM_ModelParam* modelParam = dynamic_cast<ARM_ModelParam*>(armObj);
			if( !modelParam )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams[i] = modelParam;
		}

        localSLNModelParams = new ARM_Local_SLN_ModelParams(modelParams);

		/// to avoid memory leak put an autocleaner on model params
		ARM_AutoCleaner<ARM_Local_SLN_ModelParams> Hold(localSLNModelParams);
		localSLNModel = new ARM_Local_SLN_Model( *localSLNModelParams );

		/// assign object
		if( !assignObject( localSLNModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (localSLNModel) delete localSLNModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to calibrate a Local Normal Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Local_Model_Calibrate(
	const long&				multiAssetsModelId,	
	const string&			localModelName, // supposed to be a ARM_Local_Model
    const long&				portfolioId,
	const vector<double>&	evalDates,
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* object = NULL;
	ARM_MultiAssetsModel* newMultiAssets = NULL;
	   
	try
	{ 		
		// Get ARM_MultiAssetsModel
		ARM_PricingModel* oldModel= NULL;
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !GetObjectFromId( &oldModel, multiAssetsModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}
		
		ARM_MultiAssetsModel* oldMultiAssets = dynamic_cast< ARM_MultiAssetsModel* > (oldModel);

		if (!oldMultiAssets)
		{
			result.setMsg ("ARM_ERR: inputed model is required to be a MultiAssets model");
			return ARM_KO;
		}

		// Clone pricing model
		newMultiAssets = (ARM_MultiAssetsModel*)oldMultiAssets->Clone();

		// Get local model embedded in MultiAssets model
		// ... and test that it is of type ARM_Local_Model
		ARM_PricingModelPtr model = (*newMultiAssets->GetModelMap())[localModelName]->Model();
		ARM_Local_Model* localModel = dynamic_cast<ARM_Local_Model*> (&*model);

		if (!localModel)
		{
			result.setMsg ("ARM_ERR: model name provided does not correspond to a ARM_Local_Model");
			return ARM_KO;
		}

		// Get portfolio
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(portfolioId);
        ARM_StdPortfolio* portfolio = dynamic_cast< ARM_StdPortfolio* >(object);

        if(!portfolio)
		{
			//if it is not a portfolio, maybe it is a calculator
			ARM_GenCalculator* calculator = NULL;
			if( !GetObjectFromId( &calculator, portfolioId, ARM_GENCALCULATOR) )
	        {
		        result.setMsg ("ARM_ERR: invalid portfolio");
		        return ARM_KO;
	        }
			
			// convert dates into times
			std::vector<double> evalTimes (evalDates.size());

			for (size_t i(0); i<evalDates.size(); i++)
			{
				char date[11];
				Local_XLDATE2ARMDATE (evalDates[i], date);
				ARM_Date armdate (date);
				evalTimes[i] = model->GetTimeFromDate(armdate);
			}
		
			// local model calibration
			localModel->CalibrateLocalModel(*calculator, evalTimes);
		
			/// assign object
			if( !assignObject( newMultiAssets, result, objId ) ){
				return ARM_KO; }
			else{
				return ARM_OK; }
		}
		
		// convert dates into times
		std::vector<double> evalTimes (evalDates.size());

		for (size_t i(0); i<evalDates.size(); i++)
		{
			char date[11];
			Local_XLDATE2ARMDATE (evalDates[i], date);
			ARM_Date armdate (date);
			evalTimes[i] = model->GetTimeFromDate(armdate);
		}
		
		// local model calibration
		localModel->CalibrateLocalModel(*portfolio, evalTimes);
		
		/// assign object
		if( !assignObject( newMultiAssets, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (newMultiAssets) delete newMultiAssets ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_LocalModel_CalibrateFunctional(
	const long&				multiAssetsModelId,	
	const string&			localModelName, // supposed to be a ARM_Local_Model
    const VECTOR<long>&		securitiesId,
	const VECTOR<long>&		densitiesId,
	const bool&				rescaling,
	ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* object = NULL;
	ARM_MultiAssetsModel* newMultiAssets = NULL;
	   
	try
	{ 		
		// Get ARM_MultiAssetsModel
		ARM_PricingModel* oldModel= NULL;
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !GetObjectFromId( &oldModel, multiAssetsModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}
		
		ARM_MultiAssetsModel* oldMultiAssets = dynamic_cast< ARM_MultiAssetsModel* > (oldModel);

		if (!oldMultiAssets)
		{
			result.setMsg ("ARM_ERR: inputed model is required to be a MultiAssets model");
			return ARM_KO;
		}

		// Clone pricing model
		newMultiAssets = (ARM_MultiAssetsModel*)oldMultiAssets->Clone();

		// Get local model embedded in MultiAssets model
		// ... and test that it is of type ARM_Local_Model
		ARM_PricingModelPtr model = (*newMultiAssets->GetModelMap())[localModelName]->Model();
		ARM_Local_Model* localModel = dynamic_cast<ARM_Local_Model*> (&*model);

		if (!localModel)
		{
			result.setMsg ("ARM_ERR: model name provided does not correspond to a ARM_Local_Model");
			return ARM_KO;
		}

		size_t sizeSec = securitiesId.size();
		vector<ARM_Security*> securities(sizeSec, NULL);
        int i;
        for (i = 0; i < sizeSec; i++)
		{
            if( !GetObjectFromId( &securities[i], securitiesId[i], ARM_SECURITY) )
	        {
		        result.setMsg ("ARM_ERR: Securities Vector is not of a good type, please check all products");
		        return ARM_KO;
	        }
        }

		ARM_DensityFunctor* density;
		size_t sizeDens = densitiesId.size();
		vector<ARM_DensityFunctor*> densities(sizeDens, NULL);
		for (i = 0; i < sizeDens; i++)
		{
			if(densitiesId[i] != ARM_NULL_OBJECT)
			{
				density = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(densitiesId[i]));
				if( !density )
				{
					result.setMsg ("ARM_ERR: density is not of a good type");
					return ARM_KO;
				}
				densities[i] = density;
			}
		}

		// local model calibration
		localModel->CalibrateLocalModelFunctional(securities, densities, 501, 6, rescaling);
		
		/// assign object
		if( !assignObject( newMultiAssets, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (newMultiAssets) delete newMultiAssets ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_Local_Model_VarSqueeze(
	const long&				multiAssetsModelId,	
	const string&			localModelName, // supposed to be a ARM_Local_Model
    ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Object* object = NULL;
	   
	try
	{ 		
		// Get ARM_MultiAssetsModel
		ARM_PricingModel* oldModel= NULL;
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !GetObjectFromId( &oldModel, multiAssetsModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		}
		
		ARM_MultiAssetsModel* oldMultiAssets = dynamic_cast< ARM_MultiAssetsModel* > (oldModel);

		if (!oldMultiAssets)
		{
			result.setMsg ("ARM_ERR: inputed model is required to be a MultiAssets model");
			return ARM_KO;
		}

		
		// Get local model embedded in MultiAssets model
		// ... and test that it is of type ARM_Local_Model
		ARM_PricingModelPtr model = (*oldMultiAssets->GetModelMap())[localModelName]->Model();
		ARM_Local_Model* localModel = dynamic_cast<ARM_Local_Model*> (&*model);

		if (!localModel)
		{
			result.setMsg ("ARM_ERR: model name provided does not correspond to a ARM_Local_Model");
			return ARM_KO;
		}

		string txt;
		if(localModel->GetVarianceSqueezeStatus())
		{
			txt += "Warning : variance squeeze, see event viewer for details";
		}
		else
		{
			txt += "no variance squeeze";
		}
		result.setString(txt.c_str());
		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to create a ARM_MarketIRModel
//////////////////////////////////////////////////
extern long ARMLOCAL_MarketIRModel_Create(
	const long&				mktDataManagerId,
	const vector< string >& keys,
	const int&				vnsPricingMethod,
    ARM_result&				result, 
	long					objId)
{
		/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_MarketIRModel* mktIrModel = NULL;
	   
	try
	{
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		mktIrModel = new ARM_MarketIRModel (*mktDataManager, keys, (ARM_MarketIRModel::VnsPricingMethod)vnsPricingMethod);


		/// assign object
		if( !assignObject( mktIrModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	catch(Exception& x)
	{
		if (mktIrModel) delete mktIrModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a SmiledFRM model
////////////////////////////////////////////
extern long ARMLOCAL_SmiledFRMModel_Create(
	const long&		zeroCurveId,	/// interest rate curve
	const long&		correlParamId,	/// correlationcurve
	const long&		humpId,			/// hump (forward volatility specification)
	const int&	factorsNb,		/// nb of factorss
	const int&	timeStepsNb,	
	const int&	gridSize,		
	const double&	stdDevNb,		
	const bool&		skipPDE,			/// skips PDE when possible
	const string&	correlType,		/// switch to theta-like correl
	const bool&		allowInterpol,		/// allow interpol in reset for ZC
	const string&	swaptionApprox,
	const double&	recorrel,
	const bool&		rescalling,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_SmiledFRM* mod;
		mod = new ARM_SmiledFRM(CreateClonedPtr( curve ),NULL,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol
			,(ARM_ModelParamsSmiled::CalibProxy) ARM_ArgConv_MMCalibProxy.GetNumber(swaptionApprox),
			rescalling);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SmiledFRM> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );


		modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), correlParamId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), humpId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		mod->SetModelParams(ARM_ModelParamsSmiled(
			modelParams,
			factorsNb,
			(ARM_ModelParamsSmiled::CorrelType) ARM_ArgConv_MMCorrelType.GetNumber(correlType),
			allowInterpol,
			recorrel)
			);

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

////////////////////////////////////////////
//// Function to create a Smiled Market model
////////////////////////////////////////////
extern long ARMLOCAL_SmiledMarketModel_Create(
	const string&	calibPattern,	/// calib pattern
	const long&		zeroCurveId,	/// interest rate curve
	const long&		correlParamId,	/// correlationcurve
	const long&		humpId,			/// hump (forward volatility specification)
	const int&	factorsNb,		/// nb of factorss
	const int&	timeStepsNb,	
	const int&	gridSize,		
	const double&	stdDevNb,		
	const bool&		skipPDE,			/// skips PDE when possible
	const string&	correlType,		/// switch to theta-like correl
	const bool&		allowInterpol,		/// allow interpol in reset for ZC
	const string&	swaptionApprox,
	const double&	recorrel,		
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_SmiledMM* mod;
		mod = ARM_SmiledFRMfactory.Instance()->CreateSmiledMarketModel(
					(ARM_SmiledFRMfactoryImp::CalibPattern) ARM_ArgConv_MMCalibPattern.GetNumber(calibPattern),
					CreateClonedPtr( curve ),
					NULL,
					timeStepsNb,
					gridSize,
					stdDevNb,
					skipPDE,
					allowInterpol,
					(ARM_ModelParamsSmiled::CalibProxy) ARM_ArgConv_MMCalibProxy.GetNumber(swaptionApprox)
					);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SmiledMM> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );


		modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), correlParamId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), humpId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		mod->SetModelParams(ARM_ModelParamsSmiled(
			modelParams,
			factorsNb,
			(ARM_ModelParamsSmiled::CorrelType) ARM_ArgConv_MMCorrelType.GetNumber(correlType),
			allowInterpol,
			recorrel)
			);

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

//////////////////////////////////////////////////
//// Function to create a ARM_SmiledMarketModelDS
//////////////////////////////////////////////////
extern long ARMLOCAL_SmiledMarketModelDS_Create(
	const string&	calibPattern,
	const double& startDate,
	const double& endDate,
	const long& resetFreq,			
	const long& indexFreq,			
	const int& indexType,			
	const long& resetTiming,
	const long& dayCount,
	const CCString& resetCalendar,
	const long& fwdRule,
	const long& intRule,
	const long& stubRule,
	const long& resetGap,
	ARM_result& result,                                      
	long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_DateStrip* dateStrip = NULL;
	
	try
	{
		char myStartDate[20];
		char myEndDate[20];

		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

		dateStrip = ARM_SmiledFRMfactory.Instance()->CreateSmiledMarketModelDateStrip(
			(ARM_SmiledFRMfactoryImp::CalibPattern) ARM_ArgConv_MMCalibPattern.GetNumber(calibPattern),
			(ARM_Date) myStartDate,
			(ARM_Date) myEndDate,
			resetFreq,			
			indexFreq,			
			indexType,			
			resetTiming,
			dayCount,
			resetCalendar,
			fwdRule,
			intRule,
			stubRule,
			resetGap);

		/// assign object
		if( !assignObject( dateStrip, result, objId ) )
		{
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete dateStrip;

		x.DebugPrint();
		
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete dateStrip;

		result.setMsg ("ARM_ERR: unrecognized failure");
		return ARM_KO;
	}
}

////////////////////////////////////////////
//// Function to create a SmiledFRM model
////////////////////////////////////////////
extern long ARMLOCAL_SVBGMModel_Create(
	const long&		zeroCurveId,		/// interest rate curve
	const long&		shiftId,			// shift des taux
	const long&		alphaId,			// vol initial
	const long&		nuId,				// vol de vol
	const long&		rhoId,				// correl taux(i) vol(i)
	const long&		rrcorrelParamId,	// correl taux taux
	const long&		rvcorrelParamId,	// correl taux(i) vol(k)
	const long&		vvcorrelParamId,	// correl vol(i) vol(k)
	const double&	recorrel,
	const int&		factorsNb,
	const double&	minratio,
	const bool&		Proxy,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	int paramSize = 7;

	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(paramSize);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_SVBGM* mod;
		mod = new ARM_SVBGM(CreateClonedPtr( curve ),NULL,false,true,Proxy);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SVBGM> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );


		int k = 0;
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), shiftId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (shift) is not of a good type");
			return ARM_KO;
		}
		
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), alphaId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (alpha) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), nuId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (nu) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rhoId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rho) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rrcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rrcorrelParamId) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rvcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rvcorrelParamId) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), vvcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (vvcorrelParamId) is not of a good type");
			return ARM_KO;
		}

		mod->SetModelParams(ARM_ModelParamsSVBGM(
			modelParams,
			recorrel,
			factorsNb,
			minratio)
			);

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

extern long ARMLOCAL_BGMSV1FModel_Create(
	const long&		zeroCurveId,		/// interest rate curve
	const long&		shiftId,			// shift des taux
	const long&		levelId,	
	const long&		initVarId,	
	const long&		LongTermVarId,
	const long&		VarVolId,
	const long&		VarMeanRevId,
	const long&		RhoId,
	const long&		rrcorrelParamId,	// correl taux taux
	const bool&		LocalCalibration,
	const double&	recorrel,
	const int&		factorsNb,
	const double&	minratio,
	const bool&		localrhoCalib,
	const vector<double>& stddevCalib,
	const bool&		Proxy,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	int paramSize = 8;

	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(paramSize);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_BGMSV1F* mod;
		mod = new ARM_BGMSV1F(CreateClonedPtr( curve ),NULL,false,true,Proxy);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_BGMSV1F> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );


		int k = 0;
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), shiftId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (shift) is not of a good type");
			return ARM_KO;
		}
		
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), levelId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (level) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), initVarId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (long term var) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), LongTermVarId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (long term var) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), VarVolId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (var vol) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), VarMeanRevId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (var vol) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), RhoId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rho) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rrcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rrcorrelParamId) is not of a good type");
			return ARM_KO;
		}

		mod->SetModelParams(ARM_ModelParamsBGMSV1F(
			modelParams,
			std::vector<double>(stddevCalib),
			recorrel,
			factorsNb,
			minratio,
			localrhoCalib,
			LocalCalibration)
			);

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}

extern long ARMLOCAL_BGMSV2FModel_Create(
	const long&		zeroCurveId,		/// interest rate curve
	const long&		rrcorrelParamId,	// correl taux taux
	const double&	recorrel,
	const int&		factorsNb,
	const double&	minratio,
	const double&	v01,
	const double&	kappa1,
	const double&	rho1,
	const double&	v02,
	const double&	kappa2,
	const double&	rho2,
	const double&	shift,
	const bool&		LocalRho1Calib,
	const bool&		LocalRho2Calib,
	const vector<double>& stddevCalib,
	const bool&		Proxy,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	int paramSize = 1;

	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(paramSize);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_BGMSV2F* mod;
		mod = new ARM_BGMSV2F(CreateClonedPtr( curve ),NULL,false,true,Proxy);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_BGMSV2F> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );


		int k = 0;
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rrcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rrcorrelParamId) is not of a good type");
			return ARM_KO;
		}

		mod->SetModelParams(ARM_ModelParamsBGMSV2F(
			modelParams,
			v01, kappa1, rho1, 
			v02, kappa2, rho2,
			shift,
			std::vector<double>(stddevCalib),
			recorrel,
			factorsNb,
			minratio,
			LocalRho1Calib,
			LocalRho2Calib)
			);

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_SVMMSpreadModel_Create(
	const long&		zeroCurveId,		/// interest rate curve
	const long&		levelId,	
	const long&		initVarId,	
	const long&		LongTermVarId,
	const long&		VarVolId,
	const long&		VarMeanRevId,
	const long&		RhoId,
	const long&		rrcorrelParamId,	// correl taux taux
	const double&	recorrel,
	const int&		factorsNb,
	const double&	minratio,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	int paramSize = 7;

	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(paramSize);
		

	try
	{
		ARM_ZeroCurve* curve = NULL;
		if( !GetObjectFromId( &curve, zeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_SVMMSpread* mod;
		mod = new ARM_SVMMSpread(CreateClonedPtr( curve ),NULL);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SVMMSpread> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		int k = 0;
		
		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), levelId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (level) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), initVarId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (initial var) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), LongTermVarId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (long term var) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), VarVolId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (var vol) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), VarMeanRevId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (var vol) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), RhoId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rho) is not of a good type");
			return ARM_KO;
		}

		modelParams[k] = NULL;
		if( !GetObjectFromId( &(modelParams[k++]), rrcorrelParamId, ARM_MODELPARAM ) )
		{
			result.setMsg ("ARM_ERR: model parameter (rrcorrelParamId) is not of a good type");
			return ARM_KO;
		}
/*
		mod->SetModelParams(ARM_ModelParamsBGMSV1F(
			modelParams,
			std::vector<double>(0),
			recorrel,
			factorsNb,
			minratio)
			);*/

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_HWxSVMMSpread_Create(
	const vector<string>&	modelNames,
	const vector<long>&		modelIds,
	const long&				hw2fId,
	const vector<double>&	CorrIndexEndTimes,
	const double&			ConstantCrossCorrel,
	ARM_result&				result,
	long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( modelNames.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: modelNames.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		std::vector<double> corrIdxEndTimes(0);
		
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(hw2fId);
		ARM_HullWhite2F * hw2f = dynamic_cast<ARM_HullWhite2F*>(armObj);
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = modelNames.size();
		ARM_StringVectorVector depends(nbModels);
		depends[1]		= ARM_StringVector(1);
		depends[1][0]	= modelNames[0];

		ARM_ModelNameMap modelMap( modelNames, models, depends );
		

		if(hw2f != NULL && CorrIndexEndTimes.size() > 0)
		{
			double asOf = hw2f->GetAsOfDate().GetJulian();

			corrIdxEndTimes.resize(CorrIndexEndTimes.size());

			char myDate[20];

			for(int k = 0; k < (int)corrIdxEndTimes.size(); k++)
			{
				Local_XLDATE2ARMDATE(CorrIndexEndTimes[k], myDate);
				ARM_Date armDate(myDate);
				corrIdxEndTimes[k] = armDate.GetJulian() - asOf;
			}
		}

		newModel = new ARM_HWxSVMMSpread(modelMap, hw2f, corrIdxEndTimes, ConstantCrossCorrel);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a SmiledFX model
////////////////////////////////////////////
extern long ARMLOCAL_SmiledFXModel_Create(
	const long&		domZeroCurveId,	/// domestic interest rate curve
	const long&		forZeroCurveId,	/// foreign interest rate curve
	const double&	FXSpot,
	const long&		correlParamId,	/// correlationcurve
	const long&		humpId,			/// hump (forward volatility specification)
	const int&	factorsNb,			/// nb of factors for the simulation
	const int&	timeStepsNb,	
	const int&	gridSize,		
	const double&	stdDevNb,
	const bool&		skipPDE,			/// skips PDE when possible
	const string&	correlType,		/// switch to theta-like correl
	const double&	recorrel,
	const bool&		rescalling,
	const long&		Model2IRFXId,
	ARM_result&		result, 
	long			objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	CC_STL_VECTOR( ARM_ModelParam* ) modelParams(2);
		

	try
	{
		ARM_ZeroCurve* domCurve = NULL;
		if( !GetObjectFromId( &domCurve, domZeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: domestic interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		ARM_ZeroCurve* forCurve = NULL;
		if( !GetObjectFromId( &forCurve, forZeroCurveId, ARM_ZERO_CURVE ) )
		{
			result.setMsg ("ARM_ERR: foreign interest rate fwd curve is not of a good type");
			return ARM_KO;
		}

		modelParams[0] = NULL;
	    if( !GetObjectFromId( &(modelParams[0]), correlParamId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		modelParams[1] = NULL;
	    if( !GetObjectFromId( &(modelParams[1]), humpId, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		
		ARM_ModelParamsSmiled_Fx modelParamsObj(
			modelParams,
			ARM_ZeroCurvePtr(static_cast<ARM_ZeroCurve*>(domCurve->Clone())),
			ARM_ZeroCurvePtr(static_cast<ARM_ZeroCurve*>(forCurve->Clone())),
			FXSpot,
			factorsNb,
			(ARM_ModelParamsSmiled::CorrelType) ARM_ArgConv_MMCorrelType.GetNumber(correlType),
			recorrel);

		ARM_2IRFXModel* Model2IRFX = NULL;
		if (Model2IRFXId != ARM_NULL_OBJECT)
		{
			if( !(Model2IRFX = (ARM_2IRFXModel*) LOCAL_PERSISTENT_OBJECTS->GetObject(Model2IRFXId)) )
			{
				result.setMsg ("ARM_ERR: 2IRFX model is not of a good type");
				return ARM_KO;
			}
		}

		ARM_SmiledModel_Fx* mod = new ARM_SmiledModel_Fx(CreateClonedPtr( domCurve ),&modelParamsObj,timeStepsNb,gridSize,stdDevNb,skipPDE,rescalling,Model2IRFX);

		/// to avoid memory leak
		/// put an autocleaner on the model
		ARM_AutoCleaner<ARM_SmiledModel_Fx> Hold(mod);

		/// creates the default numeraire Terminal ZC for the time being
		ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
		mod->SetNumeraire( numeraire );

		/// ok we can safely release the pointor as it passed all the code!
		Hold.Release();

		/// assign object
		if( !assignObject( mod, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}


//////////////////////////////////////////////////
//// Function to calibrate a Mixture FX Model
//////////////////////////////////////////////////

extern long ARMLOCAL_MixtureFXModel_Calibrate(
	const double& fwd,
	const double& expiry,
	const int& callPut,
	const vector<double>& strikes,
	const vector<double>& vols,
	const vector<double>& decvol,
	const vector<double>& alpha,
	const vector<double>& lambda,
	vector<double>& outParams,
	ARM_result&		result)
{
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{

		ARM_MixtureModel_Fx::CalibMixture(
			fwd,
			expiry,
			callPut,
			strikes,
			vols,
			decvol,
			alpha,
			lambda,
			outParams);

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}


long ARM_ParamsMixtureFx_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_ParamsMixture_Fx*	paramsMixture_Fx = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator" );

		std::vector<double>& lags = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Lags").GetObjectId()));
		if (!lags)
		{
			result.setMsg ("ARM_ERR: lags should be a vector");
			return ARM_KO;
		}

		std::vector<double>& volATM = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("VolATM").GetObjectId()));
		if (!volATM)
		{
			result.setMsg ("ARM_ERR: Vol ATM should be a vector");
			return ARM_KO;
		}

		std::vector<double>& decVol = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("DecVol").GetObjectId()));
		if (!decVol)
		{
			result.setMsg ("ARM_ERR: Dec Vol should be a vector");
			return ARM_KO;
		}

		std::vector<double>& shift = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Shift").GetObjectId()));
		if (!shift)
		{
			result.setMsg ("ARM_ERR: Shift should be a vector");
			return ARM_KO;
		}

		std::vector<double>& lambda = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Lambda").GetObjectId()));
		if (!shift)
		{
			result.setMsg ("ARM_ERR: Lambda should be a vector");
			return ARM_KO;
		}

		string interpolName = genericParams->GetParamValue("Interpol").GetString();
        
		paramsMixture_Fx = 
			new ARM_ParamsMixture_Fx(
				*lags,
				*volATM,
				*decVol,
				*shift,
				*lambda,
				interpolName);

		// assign object
		if ( !assignObject( paramsMixture_Fx, result, objId ) )
		{
			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}
	}

	catch(Exception& x)
	{
		delete	paramsMixture_Fx;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARM_MixtureModelFx_CreateWithParamsFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_MixtureModel_Fx*	mixtureModel_Fx = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator" );

		ARM_ZeroCurve* domCrv = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("DomCrv").GetObjectId()));
		if (!domCrv)
		{
			result.setMsg ("ARM_ERR: DomCrv should be a zc curve");
			return ARM_KO;
		}

		ARM_ZeroCurve* forCrv = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("ForCrv").GetObjectId()));
		if (!forCrv)
		{
			result.setMsg ("ARM_ERR: ForCrv should be a zc curve");
			return ARM_KO;
		}

		double spot = genericParams->GetParamValue("Spot").GetDouble();
		
		ARM_ParamsMixture_Fx* params = dynamic_cast<ARM_ParamsMixture_Fx*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MixParams").GetObjectId()));
		if (!params)
		{
			result.setMsg ("ARM_ERR: Shift should be mixture parameter");
			return ARM_KO;
		}
        
		mixtureModel_Fx = 
			new ARM_MixtureModel_Fx(
				CreateClonedPtr(domCrv),
				CreateClonedPtr(forCrv),
				spot,
				params);

		// assign object
		if ( !assignObject( mixtureModel_Fx, result, objId ) )
		{
			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}
	}

	catch(Exception& x)
	{
		delete	mixtureModel_Fx;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARM_ModelBumpParams_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_PricingModel*	model = NULL;
	ARM_PricingModel* newModel = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " Model Params bumping" );

		//Model Pricing
		long modelId = genericParams->GetParamValue("ModelId").GetObjectId();
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, modelId,"Model pricing", result ) ) return ARM_KO;

		ARM_ParamType paramType = 	(ARM_ParamType)ARM_ArgConv_ModelParam.GetNumber(genericParams->GetParamValue("ParamType").GetString());
		long rowNb = genericParams->GetParamValue("RowNumber").GetDouble();
		long columnNb = genericParams->GetParamValue("ColumnNumber").GetDouble();
		double shift = genericParams->GetParamValue("Shift").GetDouble();
		ARM_BumpParamType isCumulative = (ARM_BumpParamType)ARM_ArgConv_BumpParam.GetNumber(genericParams->GetParamValue("IsCumulative").GetString());

		newModel = (ARM_PricingModel*)model->Clone();
		newModel->GetModelParams()->GetModelParam(paramType).BumpModelParam(rowNb,columnNb,shift,isCumulative);        

		// assign object
		return assignObject( newModel, result, objId ) ? ARM_OK : ARM_KO; 
	}

	catch(Exception& x)
	{
		delete	newModel;
		newModel = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_FXModelDensity_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_GP_Matrix* densityCurve = NULL;
	ARM_PricingModel*	model = NULL;
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " Model Density Create " );

		//Model Pricing
		long modelId = genericParams->GetParamValue("ModelId").GetObjectId();
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, modelId,"Model pricing", result ) ) return ARM_KO;

		string densityType = genericParams->GetParamValue("DensityType").GetString();
		double  expiry	   = genericParams->GetParamValue("Expiry").GetDouble();
		double  xmin	   = genericParams->GetParamValue("Xmin").GetDouble();
		double  xmax	   = genericParams->GetParamValue("Xmax").GetDouble();
		double  nbPoints   = genericParams->GetParamValue("NbPoints").GetDouble();

		ARM_EqFxBase* eqfxmodel = dynamic_cast<ARM_EqFxBase*>(model);
		if( !eqfxmodel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : model has to be EqFXBase Model");

		double fwd = eqfxmodel->ComputeFwdAtTime(expiry) ;
		
		GaussLegendre_Coefficients glc( nbPoints, xmin, xmax);
		ARM_DensityFunctor* densityFunctor = eqfxmodel->GetDensityFunctor();

		if(densityType=="Normale")
			densityFunctor->SetIsDirect(true);
		else if(densityType=="InvNormale")
			densityFunctor->SetIsDirect(false);
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : DensityType must be Normale or InvNormale");

		std::vector<double>& densityOrd = (std::vector<double>&) (eqfxmodel->NormalDistribution(glc,fwd,expiry))->Clone();
		ARM_GP_Matrix* densityCurve = new ARM_GP_Matrix(nbPoints,2,0.0);
		for (int i=0; i<nbPoints; i++)
		{
			(*densityCurve).Elt(i,0)=glc.get_point(i);
			(*densityCurve).Elt(i,1) = (*densityOrd).Elt(i);	
		}
		delete densityOrd;//useless now
		// assign object
		return assignObject( densityCurve, result, objId ) ? ARM_OK : ARM_KO; 
	}
	catch(Exception& x)
	{
		delete	densityCurve;
		densityCurve = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_VanillaDensity_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();
	ARM_PricingModel*	model = NULL;
	ARM_VanillaSecurityDensity* vanillaDensity = NULL;
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " Vanilla Density Create " );

		//Model Pricing
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, genericParams->GetParamValue("ModelId").GetObjectId(),"Model pricing", result ) ) return ARM_KO;

		string isDirectStr = genericParams->GetParamValue("IsDirect").GetString().substr(0,1);
		bool isDirect = isDirectStr == "y" || isDirectStr == "Y" ?  true : false;
		double  ExpiryDate	   = genericParams->GetParamValue("ExpiryDate").GetDouble();
		char myExpiryDate[20];
        // Convert date
		Local_XLDATE2ARMDATE(ExpiryDate, myExpiryDate);
		ARM_VanillaSecurityDensity* vanillaDensity = ARM_VanillaDensityFactor.Instance()->CreateVanillaDenstityFactory( *model,myExpiryDate,isDirect);
		
		// assign object
		return assignObject( vanillaDensity, result, objId ) ? ARM_OK : ARM_KO; 
	}
	catch(Exception& x)
	{
		delete	vanillaDensity;
		vanillaDensity = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_BiSVMM_Create(
	const vector<string>&	modelNames,
	const vector<long>&		modelIds,
	const long&				corrMatrixId,
	ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( modelNames.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: modelNames.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(corrMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = modelNames.size();
		ARM_StringVectorVector depends(nbModels);
		depends[1]		= ARM_StringVector(1);
		depends[1][0]	= modelNames[0];

		ARM_ModelNameMap modelMap( modelNames, models, depends );
		

		newModel = new ARM_BiSVMM( modelMap,*correlationMatrix);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
long ARMLOCAL_NP1IRNFXModel_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		int i;
		for( i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		int nbModels = names.size();

		if((names.size()+1)%3 != 0)
		{
			ARMTHROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : the model should count 3*N-1 model and not" << names.size() << ".");
		}

		int nbCcy = (names.size()+1)/3;
		int BasisModelIdx = 2*nbCcy-1;

		ARM_StringVectorVector depends(nbModels);

		if (nbModels > BasisModelIdx)
		{
		}

		depends[BasisModelIdx] = ARM_StringVector(1,names[0]);

		for (i = 0; i < nbCcy-1; ++i)
		{
			depends[nbCcy+i]       = ARM_StringVector(2);
			if (nbModels > BasisModelIdx)
			{
				depends[nbCcy+i][0]	   = names[BasisModelIdx];
				depends[nbCcy+i][1]	   = names[BasisModelIdx+1+i];
				depends[BasisModelIdx+1+i] = ARM_StringVector(1,names[1+i]);
			}
			else
			{
				depends[nbCcy+i][0]	   = names[0];
				depends[nbCcy+i][1]	   = names[BasisModelIdx+1+i];
			}
		}

		ARM_ModelNameMap modelMap( names, models, depends );

		newModel = new ARM_NP1IRNFXModel( modelMap,*correlationMatrix);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
long ARMLOCAL_2IRFXSV_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( names.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: names.size() != modelIds.size()");
			return ARM_KO;
		};

		int nbModels = modelIds.size();

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		int i;
		for( i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_GP_Matrix* correlationMatrix = dynamic_cast<ARM_GP_Matrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		ARM_StringVectorVector depends(nbModels);
		depends[ARM_2IRFXSV::FxModel]       = ARM_StringVector(2);
		depends[ARM_2IRFXSV::FxModel][0]	   = names[ARM_2IRFXSV::DomBasisModel];
		depends[ARM_2IRFXSV::FxModel][1]	   = names[ARM_2IRFXSV::ForBasisModel];
		if (nbModels > ARM_2IRFXSV::DomBasisModel)
			depends[ARM_2IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_2IRFXSV::DomModel]);
		if (nbModels > ARM_2IRFXSV::ForBasisModel)
			depends[ARM_2IRFXSV::ForBasisModel] = ARM_StringVector(1,names[ARM_2IRFXSV::ForModel]);

		ARM_ModelNameMap modelMap( names, models, depends );

		newModel = new ARM_2IRFXSV( modelMap,*correlationMatrix);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_NP1IRNFX_CalibrateFunctional(
	const long&				NP1IRNFXId,
	const VECTOR<double>&	ResetDates,
	const VECTOR<long>&		densitiesId,
	const long&				nbRows,
	const long&				nbCols,
	const double&			sizeGrid,
	const double&			nbStdDev,
	const int&				rescaling,
	ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_NP1IRNFXModel* newNP1IRNFXModel = NULL;
	   
	try
	{ 		
		int i;

		// Get NP1IRNFX
		ARM_NP1IRNFXModel* NP1IRNFXModel = dynamic_cast<ARM_NP1IRNFXModel*>(LOCAL_PERSISTENT_OBJECTS->GetObject(NP1IRNFXId));
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !NP1IRNFXModel )
		{
			result.setMsg ("ARM_ERR: NP1IRNFXId is not of a good type");
			return ARM_KO;
		}

		// Clone NP1IRNFX model
		newNP1IRNFXModel = CreateClone(NP1IRNFXModel);

		ARM_GP_VectorPtr ResetTimes(ResetDates.empty() ?  NULL : new std::vector<double>(ResetDates));

		for (i = 0; i<ResetTimes->size(); i++)
		{
			char date[11];
			Local_XLDATE2ARMDATE ((*ResetTimes)[i], date);
			ARM_Date armdate (date);
			(*ResetTimes)[i] = newNP1IRNFXModel->GetTimeFromDate(armdate);
		}

		ARM_DensityFunctor* density;
		vector<ARM_DensityFunctor*> densities(nbRows*nbCols, NULL);
		for (i = 0; i < nbRows*nbCols; i++)
		{
			if(densitiesId[i] != ARM_NULL_OBJECT)
			{
				density = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(densitiesId[i]));
				if( !density )
				{
					result.setMsg ("ARM_ERR: density is not of a good type");
					return ARM_KO;
				}
				densities[i] = density;
			}
		}

		// local model calibration
		newNP1IRNFXModel->CalibrateFunctional(
			*ResetTimes,
			densities,
			nbRows,
			nbCols,
			sizeGrid,
			nbStdDev,
			rescaling==1);
		
		/// assign object
		if( !assignObject( newNP1IRNFXModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (newNP1IRNFXModel) delete newNP1IRNFXModel ;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_2IRFXSV_CalibrateFunctional(
	const long&				Model2IRFXSVId,
	const VECTOR<double>&	ResetDates,
	const VECTOR<long>&		densitiesId,
	const double&			sizeGrid,
	const double&			nbStdDev,
	ARM_result&				result, 
	long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_2IRFXSV* newModel2IRFXSV = NULL;
	   
	try
	{ 		
		int i;

		// Get NP1IRNFX
		ARM_2IRFXSV* Model2IRFXSV = dynamic_cast<ARM_2IRFXSV*>(LOCAL_PERSISTENT_OBJECTS->GetObject(Model2IRFXSVId));
		/// test that it is an object derived from ARM_PRICINGMODEL 
		if( !Model2IRFXSV )
		{
			result.setMsg ("ARM_ERR:2IRFXSVId is not of a good type");
			return ARM_KO;
		}

		// Clone NP1IRNFX model
		newModel2IRFXSV = CreateClone(Model2IRFXSV);

		ARM_GP_VectorPtr ResetTimes(ResetDates.empty() ?  NULL : new std::vector<double>(ResetDates));

		for (i = 0; i<ResetTimes->size(); i++)
		{
			char date[11];
			Local_XLDATE2ARMDATE ((*ResetTimes)[i], date);
			ARM_Date armdate (date);
			(*ResetTimes)[i] = Model2IRFXSV->GetTimeFromDate(armdate);
		}

		ARM_DensityFunctor* density;
		vector<ARM_DensityFunctor*> densities(densitiesId.size(), NULL);
		for (i = 0; i < densitiesId.size(); i++)
		{
			if(densitiesId[i] != ARM_NULL_OBJECT)
			{
				density = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(densitiesId[i]));
				if( !density )
				{
					result.setMsg ("ARM_ERR: density is not of a good type");
					return ARM_KO;
				}
				densities[i] = density;
			}
		}

		// local model calibration
		newModel2IRFXSV->CalibrateFunctional(
			*ResetTimes,
			densities,
			sizeGrid,
			nbStdDev);
		
		/// assign object
		if( !assignObject( newModel2IRFXSV, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if (newModel2IRFXSV) delete newModel2IRFXSV ;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_HWSBGMQtoModel_Create(
	const vector<string>&	modelNames,
	const vector<long>&		modelIds,
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( modelNames.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: modelNames.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_CurveMatrix * correlationMatrix = dynamic_cast<ARM_CurveMatrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = modelNames.size();
		ARM_StringVectorVector depends(nbModels);
		depends[1]		= ARM_StringVector(1);
		depends[1][0]	= modelNames[0];

		ARM_ModelNameMap modelMap( modelNames, models, depends );
		

		newModel = new ARM_HWSBGMQtoModel( modelMap,*correlationMatrix);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}

}

extern long ARMLOCAL_HWSVBGMQtoModel_Create(
	const vector<string>&	modelNames,
	const vector<long>&		modelIds,
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_PricingModel* newModel = NULL;
   
	try
	{	
		/// load models objects!
		if( modelNames.size() != modelIds.size() )
		{
			result.setMsg ("ARM_ERR: modelNames.size() != modelIds.size()");
			return ARM_KO;
		};

		vector< ARM_PricingModelPtr > models( modelIds.size() );
		for( size_t i=0; i<modelIds.size(); ++i )
		{
			ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(modelIds[i]);
			ARM_PricingModel* pricingModel = dynamic_cast<ARM_PricingModel*>( armObj );
			if( !pricingModel )
			{
				result.setMsg ("ARM_ERR: model is not of the good type, expected gp pricing model object!");
				return ARM_KO;
			}
			else
				models[i] = ARM_PricingModelPtr( static_cast<ARM_PricingModel*>( pricingModel->Clone() ) );
		}

		ARM_Object* armObj2=  LOCAL_PERSISTENT_OBJECTS->GetObject(correlationMatrixId);
		ARM_CurveMatrix * correlationMatrix = dynamic_cast<ARM_CurveMatrix*>(armObj2);
		if( !correlationMatrix )
		{
			result.setMsg ("ARM_ERR: correlation matrix is not of the good type!");
			return ARM_KO;
		}
		
		/// FIXME, Ugly, We should be create ModelNamemap underExcel
		int nbModels = modelNames.size();
		ARM_StringVectorVector depends(nbModels);
		depends[1]		= ARM_StringVector(1);
		depends[1][0]	= modelNames[0];

		ARM_ModelNameMap modelMap( modelNames, models, depends );
		

		newModel = new ARM_HWSVBGMQtoModel( modelMap,*correlationMatrix);

		/// assign object
		if( !assignObject( newModel, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete newModel;
		x.DebugPrint();
		ARM_RESULT();
	}

}

long ARM_2IRFX_ComputeTimeLagFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " 2IRFX Model Compute Time Lag " );

		//Model Pricing
		long model2IRFXId = genericParams->GetParamValue("ModelId").GetObjectId();
		ARM_2IRFXModel* model;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, model2IRFXId,"Model pricing", result ) ) return ARM_KO;

		double  expiryTime = genericParams->GetParamValue("ExpiryTime").GetDouble();
		double  payTime	   = genericParams->GetParamValue("PayTime").GetDouble();
		string domForFlag  = genericParams->GetParamValue("DomForFlag").GetString();

		domForFlag = StrUpper(domForFlag);

		double drift = model->ComputeTimeLag(expiryTime,payTime,(domForFlag=="DOM"?true:false));

		ResizeRetValues(1,1);
		SetValue(0,0,drift);
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}


long ARM_2IRFX_ComputeFwdFxVolFunctor::operator()( ARM_result& result, long objId )
{
		/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	ARM_GenericParams* genericParams = GetGenericParams();
	
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " 2IRFX Model Compute Fwd Fx Vol" );

		//Model Pricing
		long model2IRFXId = genericParams->GetParamValue("ModelId").GetObjectId();
		ARM_2IRFXModel* model;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, model2IRFXId,"Model pricing", result ) ) return ARM_KO;

		std::vector<double>& resetTimes = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("ResetTimes").GetObjectId()));
		if (!resetTimes)
		{
			result.setMsg ("ARM_ERR: lags should be a vector");
			return ARM_KO;
		}

		std::vector<double>& settlTimes = dynamic_cast<std::vector<double>&>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("SettlTimes").GetObjectId()));
		if (!settlTimes)
		{
			result.setMsg ("ARM_ERR: lags should be a vector");
			return ARM_KO;
		}

		ARM_GP_Matrix volatilities;
		std::vector<double> totalVar;
		
		model->ComputeVolatilities( *resetTimes, *settlTimes, volatilities, totalVar);

		ResizeRetValues( totalVar.size(), 1);

		double value(0.);
		double step(0.), nextStep(0.);
		for (int i=0; i<totalVar.size(); i++)
		{	
			value = sqrt( totalVar(i)/( (*resetTimes)[i]/K_YEAR_LEN ) );
			SetValue( 0, i, value );
		}

	} catch (Exception& x) {
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}

long ARM_2IRFX_ComputeFwdFxModelParamFunctor::operator()( ARM_result& result, long objId )
{
		/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	ARM_GenericParams* genericParams = GetGenericParams();
	
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, " 2IRFX Model Compute Fwd Fx Model Parameters : Q and Vol" );

		//Model Pricing
		long model2IRFXId = genericParams->GetParamValue("ModelId").GetObjectId();
		ARM_2IRFXModel* model;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, model2IRFXId,"Model pricing", result ) ) return ARM_KO;

		double  evalTime = genericParams->GetParamValue("EvalTime").GetDouble();
		double  settlementTime	   = genericParams->GetParamValue("SettlementTime").GetDouble();
				
		ARM_VectorPtr data = model->ComputeFwdFXModelParam(evalTime, settlementTime);

		ResizeRetValues(1, 2);
		SetValue(0,0,(*data)[0]);
		SetValue(1,0,(*data)[1]);
		
	} catch (Exception& x) {
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}
