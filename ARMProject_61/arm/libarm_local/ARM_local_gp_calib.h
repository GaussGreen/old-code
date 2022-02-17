/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calib.h,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


/*! \file ARM_local_gp_calib.h
 *
 *  \brief file for addins of the generic calibration in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */

#ifndef ARMLOCAL_GP_CALIB_H
#define ARMLOCAL_GP_CALIB_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"


////////////////////////////////////////////
//// Function to calibrate
////////////////////////////////////////////
extern long ARMLOCAL_Calibrate(
   long         modelId,
   long         C_calibmethodId,
   ARM_result&  result, 
   long         objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a calibration method
////////////////////////////////////////////
extern long ARMLOCAL_CalibMethod_Create(const long& C_pfId,
    const VECTOR<long>&		C_calibparamsId,
	const CCString&         C_methodtype,
	const double&           C_max_iter,
	const CCString&         C_targetfunctype,
	const long&             C_linkedcalibmethod,
	const long&             C_previouscalibmethod,
	const double&			C_factorNb,
	const bool&				C_validate,
	ARM_result&	result, 
    long         objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////////////////////
//// Function to create a calibration method with descriptor
////////////////////////////////////////////////////////////
extern long ARMLOCAL_CalibMethodWithDescription_Create(const CCString&  C_methodtype,
	const long&			C_pfId,
	const VECTOR<long>&    C_calibparamsId,
	const long&            C_modelfitterdesId,
	const CCString&        C_targetfunctype,
	const long&            C_linkedcalibmethod, 
	const long&            C_previouscalibmethod,
	const double&          C_FactorNb,
	const double&          C_NbIteration,
	const bool&            C_Validate,
	ARM_result&	result,  
    long         objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////////////////////
//// Function to create a calibration method with descriptor
////////////////////////////////////////////////////////////
extern long ARMLOCAL_CalibMethod2D_Create(const long& C_pf1Id,
    const long& C_pf2Id,
	const VECTOR<long>&  C_calibparams1Id,
	const VECTOR<long>&  C_calibparams2Id,
    const long&          C_modelfitterdes1,
	const long&          C_modelfitterdes2,
    const long&          C_targetfunctype,
    const long&          C_linkedcalibmethod,
    const long&          C_previouscalibmethod,
	const long&          C_calib2DDirection,
    ARM_result&	         result, 
    long         objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to set detail flag on/off
////////////////////////////////////////////
extern long ARMLOCAL_SetDetailFlagToCalibMethod(
	long C_calibmethodId,
	bool detailFlag,
	ARM_result&	result );


////////////////////////////////////////////
//// Function to get duration from a calib method
////////////////////////////////////////////
extern long ARMLOCAL_GetDurationFromCalibMethod(
	long C_calibmethodId,
	ARM_result&	result );

////////////////////////////////////////////
//// Function to get Data from a calib method
////////////////////////////////////////////
extern long ARMLOCAL_DataFromCalibMethod(long calibMethodId, 
                                    const CCString& dataType,
                                    ARM_result& result,
                                    long        objId );

extern string DataCalibMethodGetClass(const CCString& dataType);
ARM_CLASS_NAME  GetClassNameForCalibMethodData(const string& typeToGet);

extern long ARMLOCAL_Optimizer_Create( 
	const CCString& algoType,
	const double& Max_iter, 
	const double& tol, 
	const double& stepMax, 
	const bool& localSearch,
	const bool& printLevel,
    ARM_result&	 result, 
    long         objId = ARM_NULL_OBJECT_ID );

extern long ARMLOCAL_Solver_Create( 
	const CCString& algoType,
	const double& Max_iter, 
	const double& xTol, 
	const double& fxTol, 
	const double& gradTol,
	const bool& printLevel,
    ARM_result&	 result, 
    long         objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// Function to create a numerical calibration method
////////////////////////////////////////////
extern long ARMLOCAL_NumericalCalibMethod_Create(
    const long&				C_CalibDateStripId,
    const VECTOR<long>&		C_VanillaSecDensitiesId,    
	const long&				C_PortfolioId,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_SLNDensityFunctor_Create(
	const double&			C_Volatility, 
	const double&			C_Shift, 
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MixtureDensityFunctor_Create(
	const double&			C_Volatility1, 
	const double&			C_Volatility2, 
	const double&			C_Alpha, 
	const double&			C_Lambda, 
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

long ARMLOCAL_MixtureDensityFunctor_CreateWithATMVol(
	const double&			C_Fwd,
	const double&			C_Maturity,
	const double&			C_ATMVol, 
	const double&			C_DecVol, 
	const double&			C_Alpha, 
	const double&			C_Lambda, 
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_SABRDensityFunctor_Create(
	const double&			C_Alpha, 
	const double&			C_Beta, 
	const double&			C_Rho, 
	const double&			C_Nu, 
	const long&				C_SabrType,
	const double&			C_GridSize,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_BiSABRDensityFunctor_Create(
	const double&			C_Alpha1, 
	const double&			C_Beta1, 
	const double&			C_Rho1, 
	const double&			C_Nu1, 
	const double&			C_Alpha2, 
	const double&			C_Beta2, 
	const double&			C_Rho2, 
	const double&			C_Nu2, 
	const double&			C_RhoS1S2,
	const double&			C_RhoS1V2,
	const double&			C_RhoS2V1,
	const double&			C_RhoV1V2,
	const long&				C_SabrType,
	const double&			C_GridSize,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

long ARMLOCAL_HestonDensityFunctor_Create(
	const double&			C_V0,
	const double&			C_Kappa,
	const double&			C_Theta,
	const double&			C_VVol,
	const double&			C_Rho,
	const double&			C_Shift,
	const double&			C_Level,
	const double&			C_Sigma,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_NormalHestonDensityFunctor_Create(
	const double&			C_Fwd,
	const double&			C_V0,
	const double&			C_Kappa,
	const double&			C_Theta,
	const double&			C_VVol,
	const double&			C_Rho,
	const double&			C_Level,
	ARM_result&				result,
	long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_SplineDensityFunctor_Create(
	const VECTOR<double>&	C_Money, 
	const VECTOR<double>&	C_Vol, 
	const string&				C_VolType, 
	const long&				C_SmileId,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_DensityFunctor_CallOption(
	const long&	C_DensityFunctorId, 
	const double& C_Forward, 
	const double& C_Strike, 
	const double& C_Maturity, 
	ARM_result&	result);

extern long ARMLOCAL_DensityFunctor_Quantile(
	const long&	C_DensityFunctorId, 
	const double& C_Forward, 
	const double& C_Proba, 
	const double& C_Maturity, 
	ARM_result&	result);

extern long ARMLOCAL_IrFwdDensityFunctor_Create(
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_VanillaSecurityDensity_Create(
	const double&			C_ResetDate, 
	const double&			C_StartDate, 
	const double&			C_EndDate,
	const long&				C_DensityFunctorId,
	const long&				C_Frequency,
	const long&				C_DayCount,
	const long&				C_StubRule,
	const double&			C_Weight,
	const double&			C_AdjFwdAdd,
	const double&			C_AdjFwdMult,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_VanillaSecurityDensitySpread_Create(
	const double&			C_ResetDate, 
	const double&			C_StartDate1, 
	const double&			C_EndDate1,
	const double&			C_StartDate2, 
	const double&			C_EndDate2,
	const long&				C_DensityFunctorId,
	const long&				C_Frequency1,
	const long&				C_DayCount1,
	const long&				C_Frequency2,
	const long&				C_DayCount2,
	const long&				C_StubRule,
	const double&			C_Weight,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

long ARMLOCAL_VanillaSecurityDensityFX_Create(
	const double&			C_ResetDate,
	const long&				C_DensityFunctorId,
	const long&				C_DomCurveId,
	const long&				C_ForCurveId,
	const double&			C_FXSpot,
	ARM_result&				result, 
    long					objId);

////////////////////////////////////////////
//// For HW2F only
////////////////////////////////////////////
extern long ARMLOCAL_CalibMethodHW2F_Create(
    const long&				C_PortfolioId1,
	const long&				param1Id,
	const long&				C_PortfolioId2,
	const long&				param2Id,
	const long&				C_PortfolioId3,
	const long&				param3Id,
	const int&				withOptim,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Basket decomp
////////////////////////////////////////////
extern long ARMLOCAL_BasketDecomp_Create(
    const VECTOR<long>&		C_securitiesId,
	const VECTOR<long>&		C_modelsId,
	const long&				C_datestripId,
	const long&				C_mkmoId,
	const VECTOR<double>&	C_weightsId,
	const double&			C_side,
	const string&			C_method,
	const string&			C_strike,
	const double&			C_notional,
	const long&				C_notionalId,
	const double&			C_fees,
	const long&				C_feesId,
	ARM_result&				result, 
    long					objId = ARM_NULL_OBJECT_ID);

extern string BsktGetTypeToClass(const string& typeToGet, long gcId);

extern long ARMLOCAL_Basket_Get(
	const long&		basketId,
	const string&	getType,
	ARM_result&		result, 
	long			objId );

////////////////////////////////////////////
//// Smile Viewer
////////////////////////////////////////////
extern long ARMLOCAL_SmileViewer_Create(
    const long&		C_securityId,
	const long&		C_modelId,
	const VECTOR<double>&	C_moneyness,
	const string&	C_moneyType,
	const VECTOR<double>&	C_strikes,
	ARM_result&		result, 
    long			objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_DensityFunctorGen_Create(
	const long&		C_securityId,
	const long&		C_modelId,
	const double&	C_decStrike,
	const bool&		C_isDirect,
	const double&	C_minProba,
	const double&	C_maxProba,
	ARM_result&		result, 
    long			objId = ARM_NULL_OBJECT_ID);

class ARM_QDensityFunctor_Create : public ARM_GenericAddinFunctor
{
public:
	ARM_QDensityFunctor_Create() {}

	virtual	long operator()( ARM_result& result, long objId );
};

#endif
