/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_local_gp_model.h
 *
 *  \brief file for the model part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARMLOCAL_GP_MODEL_H
#define ARMLOCAL_GP_MODEL_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include <GP_Base\gpbase\gplinalgtypedef.h>
#include "ARM_result.h"
using ARM::ARM_GP_Vector;

////////////////////////////////////////////
//// Function to create an ir forward model
////////////////////////////////////////////
extern long ARMLOCAL_IRFwd_Create(
	long zeroCurveId,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create an inf forward model
////////////////////////////////////////////
extern long ARMLOCAL_InfFwd_Create(
	long infCurveId,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create an inflation equity model
////////////////////////////////////////////
extern long ARMLOCAL_InflationEquityModel_Create(
	long infCurveId,
	double publicationLag,
	long param1Id,
	long param2Id,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create an SABR equity model
////////////////////////////////////////////
extern long ARMLOCAL_SABREquityModel_Create(
	const long& ZcCurveId,
	const double& Spot,
	const vector<long>& paramsIds,
	ARM_result&	result, 
	long objId );

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
    const CCString&         C_Currency,	
	ARM_result&             result, 
	long                    objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a Cst Model Parameter
////////////////////////////////////////////
extern long ARMLOCAL_CstModelParam_Create(
    long                    modelParamType,
    double					paramValue,
	bool					adviseBreakPointTimes,
	ARM_result&             result, 
	long                    objId = ARM_NULL_OBJECT_ID);


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
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to get a model param from a model
////////////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelParam(
	long modelId,
    long modelParamType,
	long dataType,
	double index,
	long factorNb,
	VECTOR<double>& data,
	long& rows,
	long& cols,
	ARM_result& result);


////////////////////////////////////////////
//// Function to get a model param object from a model
////////////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelParamId(
	const long& modelId,
    const long& modelParamType,
	const long& factorNb,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to get a model Map from a model
////////////////////////////////////////////
extern long ARMLOCAL_PricingModel_GetModelMap(
	const long& modelId,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to set a model Map to a model
////////////////////////////////////////////
extern long ARMLOCAL_PricingModel_SetModelMap(
	long modelId,
	long modelMapId,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a HW1F model
////////////////////////////////////////////
extern long ARMLOCAL_HW1FModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	vector<double>& flags,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a HW2F model
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
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a MF1F model
////////////////////////////////////////////
extern long ARMLOCAL_MarkovFunctionalModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a MSV1F model
////////////////////////////////////////////
extern long ARMLOCAL_MSV1FModel_Create(
	const long &       zeroCurveId,
	const vector<long>&		modelParamVec,
	const double&       forwardTerm,
	const long &        IRIndexId,
	const string&		fwdType,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a HWSV1F model
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
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a HWSV2F model
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
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create an EQHWSM_Create
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_Create(
	const long&				zeroCurveId,
	const long&				modelParamsId,
	const long&				numMethodsId,
	const double&			dilatation,
	ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// Function to create an EQHWSM_ModelParamsCreate
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_ModelParamsCreate(
	const vector<long>&		modelParamVec,
	ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// Function to create a EQHSVM_NumMethodCreate
////////////////////////////////////////////

extern long ARMLOCAL_EQHWSV_NumMethodsCreate(
	const long&				intStep,
	const double&			imAxis,
	const double&			maxDecay,
	ARM_result&				result, 
	long					objId= ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// Function to create a FRMSV model
////////////////////////////////////////////
extern long ARMLOCAL_FRMSVModel_Create(
	const long &       zeroCurveId,
	const vector<long>&		modelParamVec,
	const vector<long>&		modelParamVec2,
	const long &        IRIndexId,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create an SFRM Model
////////////////////////////////////////////
extern long ARMLOCAL_SFRMModel_Create(
	const long&		zeroCurveId,	/// interest rate curve
	const VECTOR<long >& paramsIdVec,/// volatility,mean reversion,shift,correlation
	const long&		volType,		/// DIAG or Row
	const double&	factorsNb,		/// nb of factors
	const long&		IRIndexId,		/// irindex
	const long&		shiftConvPortId,/// portfolio id for beta conversion to shift
	const bool&		nonParamDrift, /// drift diffusion flag
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
extern long ARMLOCAL_SFRMModel_VolSwapVolFRADump(
	const long&		swaptionId,	/// Swaption
	const long&		sfrmId,	/// SFRM Model
	VECTOR<double>& OuputMatrix,
	long& rowsOutput,
	long& colsOutput,
	ARM_result&		result);

////////////////////////////////////////////
//// Function to set a fix scheduler on a 
//// SFRM model
////////////////////////////////////////////
extern long ARMLOCAL_SetSFRMFixScheduler(
	const long&	SFRMModelId,	/// model Id
	const double&	StartDate,		/// startDate
	const double&	EndDate,		/// startDate
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


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
	        long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to set a ZC cUrve in An IR model
////////////////////////////////////////////
extern long ARMLOCAL_Model_SetZCCurve(
	long modelId,
	long zcCurveId,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a QGM1F model
////////////////////////////////////////////
extern long ARMLOCAL_QGM1FModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	long        param3Id,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);


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
    const ARM_GP_Vector& Tp,
    const ARM_GP_Vector& YF,
	ARM_result&	result);

////////////////////////////////////////////
//// Function to create a QGM2F model
////////////////////////////////////////////
extern long ARMLOCAL_QGM2FModel_Create(
	long        zeroCurveId,
	const VECTOR<long >& paramsFactor1IdVect,
	const VECTOR<long >& paramsFactor2IdVect,
	ARM_result& result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a Q1F Model
////////////////////////////////////////////
extern long ARMLOCAL_Q1FModel_Create(
	const long&         zeroCurveId,
	const vector<long>& modelParamVec,
	const bool&				degenerateInHW,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a Q1F Analytic Model
////////////////////////////////////////////
extern long ARMLOCAL_Q1FAnaModel_Create(
	long        zeroCurveId,
	long        param1Id,
	long        param2Id,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a surface model param
//////////////////////////////////////////////////
extern long ARMLOCAL_SurfaceParam_Create(
    const long&		modelParamType,
    const long&		surfaceId,
	const double&	LowerBoundary,
	const double&	UpperBoundary,
    const CCString&	modelParamName,
	const bool&		adviseBreakPointTimes,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a surface List model param
//////////////////////////////////////////////////
extern long ARMLOCAL_SurfaceListParam_Create(
    const long&		modelParamType,
	const vector<double>&		Index,
    const vector<long>&		surfaceIdList,
    const CCString&	modelParamName,
	ARM_result&	result, 
	long        objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a Heston Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Heston_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
	ARM_result&	result, 
	long        objId		= ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a Shifted Heston Model
//////////////////////////////////////////////////

extern long ARMLOCAL_ShiftedHeston_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
	const bool&				isMC,
	const long&				nbSteps,
	const long&				nbSimulations,
	const long&				nbIntegrationSteps,
	ARM_result&	result, 
	long        objId		= ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a CEV Model
//////////////////////////////////////////////////

extern long ARMLOCAL_CEV_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId= ARM_NULL_OBJECT_ID );

//////////////////////////////////////////////////
//// Function to create a BS Model
//////////////////////////////////////////////////

extern long ARMLOCAL_BS_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID );

//////////////////////////////////////////////////
//// Function to create a Merton Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Merton_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId  = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a Normal Model
//////////////////////////////////////////////////

extern long ARMLOCAL_Normal_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a SABR Model
//////////////////////////////////////////////////

extern long ARMLOCAL_SABR_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
	const long&             impliedVolType,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a SLN Model
//////////////////////////////////////////////////

extern long ARMLOCAL_SLN_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result,
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to get variance squeeze info froma local model
//////////////////////////////////////////////////
extern long ARMLOCAL_GetVarianceSqueeze(
	const long&			modelId,
	 long&				nbRows,
     long&				nbCols,
	 vector<double>&    values,
	 bool&					status,
    ARM_result&			result);

//////////////////////////////////////////////////
//// Function to create a HW1F Model Param
//////////////////////////////////////////////////

extern long ARMLOCAL_HW1FModelParam_Create(
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a QNF Model Param
//////////////////////////////////////////////////

extern long ARMLOCAL_QNFModelParam_Create(
	const long&				qparamId,
    const vector<long>&		modelParamVec,
	const long&				correlMatId,
    ARM_result&				result,
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a QNF Model
//////////////////////////////////////////////////

extern long ARMLOCAL_QNFModel_Create(
	const long&				curveId,
    const long&				QNFmodelParamId,
	const bool&				degenerateInHW,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);



//////////////////////////////////////////////////
//// Function to create a Q1F FX Model
//////////////////////////////////////////////////
extern long ARMLOCAL_FXModel_Create(
	const long&				curveId,
	const vector<long>&		modelParamVec,
	const double&			spot,
	const long&				forCurveId,
	const string&			mcScheme,
	const long&				modelType,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a Q1F Eq Model
//////////////////////////////////////////////////
extern long ARMLOCAL_EqModel_Create(
	const long&				curveId,
	const vector<long>&		modelParamVec,
	const double&			spot,
	const long&				modelType,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_ModelNameMap_Create(
	const vector<string>&				names,
    const vector<long>&					modelsId,
	const vector< vector< string > >&	otherModelNames,
    ARM_result&							result, 
	long								objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_MultiAssetsModel_Create(
	const long&				modelNameMapId,
	const long&				correlMatId,
	const string&           MultiAssetsName,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to set a model to a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_SetModelToModelMap(
    const long&	modelMapId,
	const string& name,
	const long& modelId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to Get a model from a model name map
//////////////////////////////////////////////////
extern long ARMLOCAL_GetModelFromModelMap(
    const long&	modelMapId,
	const string& name,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a 2 interest rate Fx Model
//////////////////////////////////////////////////
extern long ARMLOCAL_Create2IRFXModel(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a 2 interest rate Fx Model
//////////////////////////////////////////////////
extern long ARMLOCAL_Create1IRFXModel(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const long& ModelId2IRFXId,
    ARM_result&	result, 
	long objId );

//////////////////////////////////////////////////
//// Function to create a HWHWQto Model
//////////////////////////////////////////////////
extern long ARMLOCAL_HWHWQtoModel_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const bool& fxFlag,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a HWHW2FQto Model
//////////////////////////////////////////////////
extern long ARMLOCAL_HWHW2FQtoModel_Create(
	const vector< string >& names, 
	const vector< long >& modelIds, 
	const long& correlationMatrixId,
	const bool& fxFlag,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create an interest rate fwd margin model
//////////////////////////////////////////////////
extern long ARMLOCAL_CreateFwdMarginModel(
	const long& C_basisZcCurveId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to set the reference model name to a multi asset model
//////////////////////////////////////////////////
extern long ARMLOCAL_SetRefModelNameToMultiAsset(
	const string& C_Name,
	const long& C_MultiAssetModelId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a Local Normal Model
//////////////////////////////////////////////////
extern long ARMLOCAL_LocalNormal_Model_Create(
	const long&				curveId,
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

//////////////////////////////////////////////////
//// Function to create a Local Shifted LogNormal Model
//////////////////////////////////////////////////
extern long ARMLOCAL_LocalSLN_Model_Create(
    const vector<long>&		modelParamVec,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to calibrate a Local Model
//////////////////////////////////////////////////
extern long ARMLOCAL_Local_Model_Calibrate(
	const long&				multiAssetsModelId,
	const string&			localModelName, // supposed to be a ARM_Local_Model
    const long&				portfolioId,
	const vector<double>&	evalDates,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_Local_Model_VarSqueeze(
	const long&				multiAssetsModelId,
	const string&			localModelName,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_LocalModel_CalibrateFunctional(
	const long&				multiAssetsModelId,
	const string&			localModelName, // supposed to be a ARM_Local_Model
    const VECTOR<long>&		securitiesId,
	const VECTOR<long>&		densitiesId,
    const bool&				rescaling,
	ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a ARM_MarketIRModel
//////////////////////////////////////////////////
extern long ARMLOCAL_MarketIRModel_Create(
	const long&				mktDataManagerId,
	const vector< string >& keys,
	const int&				vnsPricingMethod,
    ARM_result&				result, 
	long					objId = ARM_NULL_OBJECT_ID);


//////////////////////////////////////////////////
//// Function to create a ARM_SmiledFRM
//////////////////////////////////////////////////
extern long ARMLOCAL_SmiledFRMModel_Create(
	const long&		zeroCurveId,	/// interest rate curve
	const long&		correlParamId,	/// correlationcurve
	const long&		humpId,			/// hump (forward volatility specification)
	const int&	factorsNb,		/// nb of factorss
	const int&	timeStepsNb,	
	const int&	gridSize,		
	const double&	stdDevNb,		
	const bool&		skipPDE,			/// skips PDE when possible
	const string&	switchToTheta,		/// switch to theta-like correl
	const bool&		allowInterpol,		/// allow interpol in reset for ZC
	const string&	swaptionApprox,
	const double&	recorrel,
	const bool&		rescalling,
	ARM_result&		result, 
	long			objId );

//////////////////////////////////////////////////
//// Function to create a ARM_SmiledMarketModel
//////////////////////////////////////////////////
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
	long			objId );

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
	long objId);

//////////////////////////////////////////////////
//// Function to create a ARM_SVBGM
//////////////////////////////////////////////////
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
	long			objId );

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
	const bool&		LocalRhoCalib,
	const vector<double>& stddevCalib,
	const bool&		Proxy,
	ARM_result&		result, 
	long			objId );

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
	long			objId );

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
	long			objId );


//////////////////////////////////////////////////
//// Function to create a ARMLOCAL_SmiledFXModel_Create
//////////////////////////////////////////////////

extern long ARMLOCAL_SmiledFXModel_Create(
	const long&		domZeroCurveId,	/// domestic interest rate curve
	const long&		forZeroCurveId,	/// foreign interest rate curve
	const double&	FXSpot,
	const long&		correlParamId,	/// correlationcurve
	const long&		humpId,			/// hump (forward volatility specification)
	const int&	factorsNb,			//  nb Factors for the simulation
	const int&	timeStepsNb,	
	const int&	gridSize,		
	const double&	stdDevNb,		
	const bool&		skipPDE,			/// skips PDE when possible
	const string&	correlType,		/// switch to theta-like correl
	const double&	recorrel,
	const bool&		rescalling,
	const long& Model2IRFXId,
	ARM_result&		result, 
	long			objId );

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
	ARM_result&		result);

class ARM_ParamsMixtureFx_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_ParamsMixtureFx_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


class ARM_MixtureModelFx_CreateWithParamsFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_MixtureModelFx_CreateWithParamsFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_ModelBumpParams_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_ModelBumpParams_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_FXModelDensity_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_FXModelDensity_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


class ARM_VanillaDensity_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_VanillaDensity_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

extern long ARMLOCAL_BiSVMM_Create(
	const vector<string>&	modelNames,
	const vector<long>&		Models,
	const long&				corrMatrix,
	ARM_result&				result, 
	long					objId );


extern long ARMLOCAL_NP1IRNFXModel_Create(
const vector< string >& names, 
const vector< long >& modelIds, 
const long& correlationMatrixId,
ARM_result&	result, 
long objId );

extern long ARMLOCAL_2IRFXSV_Create(
const vector< string >& names, 
const vector< long >& modelIds, 
const long& correlationMatrixId,
ARM_result&	result, 
long objId );


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
long					objId );


extern long ARMLOCAL_2IRFXSV_CalibrateFunctional(
const long&				Model2IRFXSVId,
const VECTOR<double>&	ResetDates,
const VECTOR<long>&		densitiesId,
const double&			sizeGrid,
const double&			nbStdDev,
ARM_result&				result, 
long					objId );

extern long ARMLOCAL_HWxSVMMSpread_Create(
	const vector<string>&	modelNames,
	const vector<long>&		Models,
	const long&				hw2fId,
	const vector<double>&	CorrIndexEndTimes,
	const double&			ConstantCrossCorrel,
	ARM_result&				result,
	long					objId);


extern long ARMLOCAL_HWSBGMQtoModel_Create(
	const vector<string>&	modelNames,
	const vector<long>&		ModelsId,
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);


extern long ARMLOCAL_HWSVBGMQtoModel_Create(
	const vector<string>&	modelNames,
	const vector<long>&		ModelsId,
	const long& correlationMatrixId,
    ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);

class ARM_2IRFX_ComputeTimeLagFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_2IRFX_ComputeTimeLagFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_2IRFX_ComputeFwdFxVolFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_2IRFX_ComputeFwdFxVolFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_2IRFX_ComputeFwdFxModelParamFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_2IRFX_ComputeFwdFxModelParamFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

#endif
