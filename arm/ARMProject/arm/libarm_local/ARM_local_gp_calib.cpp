/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calib.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */

 
/*! \file ARM_local_gp_calib.cpp,
 *
 *  \brief file for the generic calibration local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date June 2003
 */


#include "firstToBeIncluded.h"

///////////////////////////////////////////////
/// WARNING include these headers FIRST
/// because it uses the new style headers
/// CCxl uses cstyle headers and leads to conflict
/// if defined first
///////////////////////////////////////////////
#include <GP_Base\gpbase\curve.h>
#include <GP_Base\gpbase\stringmanip.h>
#include <GP_Base\gpbase\datestrip.h>

#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\pricingmodelir.h>
#include <GP_Infra\gpinfra\modelparam.h>
#include <GP_Infra\gpinfra\calibdirection.h>

#include <GP_Calib\gpcalib\calibmethod.h>
#include <GP_Calib\gpcalib\typedef.h>
#include <GP_Calib\gpcalib\argconvdefault.h>
#include <GP_Calib\gpcalib\modelfitterdes.h>
#include <GP_Calib\gpcalib\vanillasecuritydensity.h>
#include <GP_Calib\gpcalib\densityfunctors.h>
#include <GP_Calib\gpcalib\basket.h>
#include <GP_Calib\gpcalib\densityfunctors.h>

#include "ARM_local_gp_calib.h"
#include "ARM_local_class.h"
#include "ARM_local_persistent.h"
#include "ARM_result.h"
#include "ARM_local_glob.h"
#include "ARM_local_wrapper.h"

#include "glob\expt.h"
#include "inst\portfolio.h"
#include "inst\portfolio.h"
#include <mod/model.h>

//// using the namespace directive to access ARM object!
using ARM::ARM_ModelParam;
using ARM::ARM_ModelParamVector;
using ARM::ARM_CalibMethod;
using ARM::ARM_PortfolioVector;
using ARM::ARM_CalibMethodVector;
using ARM::ARM_PricingModel;
using ARM::ARM_PricingModelIR;
using ARM::ARM_CalibDirection;
using ARM::ARM_CalibrationTarget;
using ARM::ARM_CalibMethodType;
using ARM::ARM_MktTargetType;
using ARM::ARM_MethodType;
using ARM::ARM_SolverType;
using ARM::ARM_OptimizerType;
using ARM::ARM_BasketCalib;
using ARM::ARM_ArgConv_BasketCalibrationType;
using ARM::ARM_ArgConv_BasketCalibrationStrike;
using ARM::ARM_SmileViewer;
using ARM::ARM_GenDensityFunctor;
using ARM::ARM_ArgConv_MoneyType;

using ARM::ARM_VectorPtr;
using ARM::ARM_StdPortfolioPtr;
using ARM::ARM_CalibMethodPtr;
using ARM::ARM_CalibMethodPtrVector;
using ARM::stringGetUpper;
using ARM::ARM_ArgConv_CalibMethod;
using ARM::ARM_ArgConv_TargetFuncMethod;
using ARM::ARM_ArgConvReverse_CalibMethod;
using ARM::ARM_ArgConv_SolverTypeMethod;
using ARM::ARM_ArgConv_OptimizerTypeMethod;
using ARM::ARM_ModelFitterDes;
using ARM::ARM_ModelFitterDesPtr;
using ARM::ARM_DensityFunctor;
using ARM::ARM_VanillaSecurityDensity;
using ARM::ARM_VanillaSecurityDensitySpread;
using ARM::ARM_VanillaSecurityDensityFX;
using ARM::ARM_VanillaSecDensityPtr;
using ARM::ARM_VanillaSecDensityPtrVector;
using ARM::ARM_DensityFunctorPtr;
using ARM::std::vector<double>;
using ARM::ARM_GP_VectorPtr;
using ARM::ARM_ShiftedLNDensityFunctor;
using ARM::ARM_MixtureDensityFunctor;
using ARM::ARM_NoVolDensityFunctor;
using ARM::ARM_SABRDensityFunctor;
using ARM::ARM_BiSABRDensityFunctor;
using ARM::ARM_SplineDensityFunctor;
using ARM::ARM_NormalHestonDensityFunctor;
using ARM::ARM_DateStrip;
using ARM::ARM_DateStripPtr;
using ARM::ARM_ZeroCurvePtr;
using ARM::ARM_QDensityFunctor;
using ARM::ARM_HestonDensityFunctor;






////////////////////////////////////////////
//// Function to create a Calib method
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
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* method = NULL;
	try
	{        
        ARM_StdPortfolio* portfolio = NULL;   
	    if( !GetObjectFromId( &portfolio, C_pfId, ARM_PORTFOLIO) )
	    {
		    result.setMsg ("ARM_ERR: portefolio is not of a good type");
		    return ARM_KO;
	    }
           
        size_t size = C_calibparamsId.size();
        ARM_ModelParamVector calibParams;
        ARM_ModelParam* param = NULL;

        size_t i;
        for(i=0; i<size; ++i)
        {                
	        if( !GetObjectFromId( &param, C_calibparamsId[i], ARM_MODELPARAM) )
	        {
   		        result.setMsg ("ARM_ERR: calib param is not of a good type");
		        return ARM_KO;
	        }
            calibParams.push_back(param);
        }

        ARM_CalibMethod* linkedcalibmethod = NULL;
        if (C_linkedcalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &linkedcalibmethod, C_linkedcalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: linked calib method is not of a good type");
			    return ARM_KO;
            }
        }

        ARM_CalibMethod* previouscalibmethod = NULL;
        if (C_previouscalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &previouscalibmethod, C_previouscalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: previous calib method is not of a good type");
			    return ARM_KO;
            }
        }
        bool calibMethodIsShared = true;
		string MethodTypeStr = CCSTringToSTLString(C_methodtype);
		string MktTargetStr = CCSTringToSTLString(C_targetfunctype);
		ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);
		ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);

        method = new ARM_CalibMethod(ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio->Clone()),
                    calibParams,
                    methodType, 
                    (size_t)C_max_iter,
                    mktTargetType,
                    linkedcalibmethod,
                    previouscalibmethod,
                    calibMethodIsShared,
					(size_t)C_factorNb,
					1,
					C_validate);  

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete method;
		result.setMsg ("ARM_ERR: unrecognized failure in Creating A Calib Method");
		return ARM_KO;
	}
}

////////////////////////////////////////////////////////
//// Function to create a Calib method with description
///////////////////////////////////////////////////////
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
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* method = NULL;
	try
	{        
        ARM_StdPortfolio* portfolio = NULL;   
	    if( !GetObjectFromId( &portfolio, C_pfId, ARM_PORTFOLIO) )
	    {
		    result.setMsg ("ARM_ERR: portefolio is not of a good type");
		    return ARM_KO;
	    }
           
        size_t size = C_calibparamsId.size();
        ARM_ModelParamVector calibParams;
        ARM_ModelParam* param = NULL;

        size_t i;
        for(i=0; i<size; ++i)
        {                
	        if( !GetObjectFromId( &param, C_calibparamsId[i], ARM_MODELPARAM) )
	        {
   		        result.setMsg ("ARM_ERR: calib param is not of a good type");
		        return ARM_KO;
	        }
            calibParams.push_back(param);
        }

		
		ARM_ModelFitterDes* ModelFitterDes = NULL;   
		if( !GetObjectFromIdWithDynamicCastCheck( &ModelFitterDes, C_modelfitterdesId  ) )
		{
		    result.setMsg ("ARM_ERR: Model Fitter Descriptor is not of a good type");
			return ARM_KO;
		};


        ARM_CalibMethod* linkedcalibmethod = NULL;
        if (C_linkedcalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &linkedcalibmethod, C_linkedcalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: linked calib method is not of a good type");
			    return ARM_KO;
            }
        }

        ARM_CalibMethod* previouscalibmethod = NULL;
        if (C_previouscalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &previouscalibmethod, C_previouscalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: previous calib method is not of a good type");
			    return ARM_KO;
            }
        }
        bool calibMethodIsShared = true;
		string MethodTypeStr = CCSTringToSTLString(C_methodtype);
		string MktTargetStr = CCSTringToSTLString(C_targetfunctype);
        method = new ARM_CalibMethod(ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio->Clone()),
            calibParams,
			(ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr),
			ModelFitterDes, 
           (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr),
            linkedcalibmethod,
            previouscalibmethod,
            calibMethodIsShared,
			(size_t)C_FactorNb,
			(size_t)C_NbIteration,
			C_Validate);  

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete method;
		result.setMsg ("ARM_ERR: unrecognized failure in Creating A Calib Method");
		return ARM_KO;
	}
}


////////////////////////////////////////////////////////
//// Function to create a Calib method 2D with description
///////////////////////////////////////////////////////
extern long ARMLOCAL_CalibMethod2D_Create(const long& C_pf1Id,
		 const long& C_pf2Id,
         const VECTOR<long>&    C_calibparams1Id,
		 const VECTOR<long>&    C_calibparams2Id,
         const long&            C_modelfitterdes1Id,
		 const long&            C_modelfitterdes2Id,
         const long&            C_targetfunctype,
         const long&            C_linkedcalibmethod, 
         const long&            C_previouscalibmethod,
		 const long&            C_calib2DDirection,
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* method = NULL;
	try
	{        
        ARM_StdPortfolio* portfolio1 = NULL;   
	    if( !GetObjectFromId( &portfolio1, C_pf1Id, ARM_PORTFOLIO) )
	    {
		    result.setMsg ("ARM_ERR: first portefolio is not of a good type");
		    return ARM_KO;
	    }

		ARM_StdPortfolio* portfolio2 = NULL;   
	    if( !GetObjectFromId( &portfolio2, C_pf2Id, ARM_PORTFOLIO) )
	    {
		    result.setMsg ("ARM_ERR: second portefolio is not of a good type");
		    return ARM_KO;
	    }
           
        size_t size = C_calibparams1Id.size();
        ARM_ModelParamVector calibParams1;
        ARM_ModelParam* param1 = NULL;

        size_t i;
        for(i=0; i<size; ++i)
        {                
	        if( !GetObjectFromId( &param1, C_calibparams1Id[i], ARM_MODELPARAM) )
	        {
   		        result.setMsg ("ARM_ERR: calib param is not of a good type");
		        return ARM_KO;
	        }
            calibParams1.push_back(param1);
        }

		size = C_calibparams2Id.size();
        ARM_ModelParamVector calibParams2;
        ARM_ModelParam* param2 = NULL;

        for(i=0; i<size; ++i)
        {                
	        if( !GetObjectFromId( &param2, C_calibparams2Id[i], ARM_MODELPARAM) )
	        {
   		        result.setMsg ("ARM_ERR: calib param is not of a good type");
		        return ARM_KO;
	        }
            calibParams2.push_back(param2);
        }

		
		ARM_ModelFitterDes* ModelFitterDes1 = NULL;   
		if( !GetObjectFromIdWithDynamicCastCheck( &ModelFitterDes1, C_modelfitterdes1Id  ) )
		{
		    result.setMsg ("ARM_ERR: Model Fitter Descriptor is not of a good type");
			return ARM_KO;
		};

		ARM_ModelFitterDes* ModelFitterDes2 = NULL;   
		if( !GetObjectFromIdWithDynamicCastCheck( &ModelFitterDes2, C_modelfitterdes2Id  ) )
		{
		    result.setMsg ("ARM_ERR: Model Fitter Descriptor is not of a good type");
			return ARM_KO;
		};


        ARM_CalibMethod* linkedcalibmethod = NULL;
        if (C_linkedcalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &linkedcalibmethod, C_linkedcalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: linked calib method is not of a good type");
			    return ARM_KO;
            }
        }

        ARM_CalibMethod* previouscalibmethod = NULL;
        if (C_previouscalibmethod != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &previouscalibmethod, C_previouscalibmethod, ARM_CALIBMETHOD))
            {
			    result.setMsg ("ARM_ERR: previous calib method is not of a good type");
			    return ARM_KO;
            }
        }
        bool calibMethodIsShared = true;
       
		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
		delete method;
		result.setMsg ("ARM_ERR: unrecognized failure in Creating A Calib Method");
		return ARM_KO;
	}
}




////////////////////////////////////////////
//// Function to calibrate a pricing model
////////////////////////////////////////////
extern long ARMLOCAL_Calibrate(long modelId,
                                long C_calibmethodId,
                                 ARM_result& result, 
                                 long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_PricingModel* mod= NULL;
    ARM_CalibMethod* calibMethod = NULL;
    try
	{
		ARM_PricingModelIR* oldModel = NULL;
		if( !GetObjectFromId( &oldModel, modelId, ARM_PRICINGMODEL) )
		{
			result.setMsg ("ARM_ERR: pricing model of a good type");
			return ARM_KO;
		}
		
        if( !GetObjectFromId( &calibMethod, C_calibmethodId, ARM_CALIBMETHOD) )
		{
            result.setMsg ("ARM_ERR: calib method is not of a good type");
            return ARM_KO;
        }   
        bool calibMethodIsShared = true;
        calibMethod->SetIsCalibMethodShared (calibMethodIsShared);
        /// clone it to avoid side effect
		mod = (ARM_PricingModel*) oldModel->Clone();
		calibMethod->Calibrate(mod);

        /// get the real class Name
		ARM_CLASS_NAME className = mod->GetName();

		/// assign object 
		if( !assignObject( mod, result, objId ))        
			return ARM_KO;        
		else        
			return ARM_OK;        
    }
	
	catch(Exception& x)
	{
        delete mod;
		x.DebugPrint();
		ARM_RESULT();
	}

	/// catch the rest
	catch (...)
	{
        delete mod;
		result.setMsg ("ARM_ERR: unrecognized failure in Generic Calibration");
		return ARM_KO;
	}
}


////////////////////////////////////////////
//// Function to set detail flag on/off
////////////////////////////////////////////
extern long ARMLOCAL_SetDetailFlagToCalibMethod(
	long C_calibmethodId,
	bool detailFlag,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* calibmethod= NULL;

	try
	{   
	    if( !GetObjectFromId( &calibmethod, C_calibmethodId, ARM_CALIBMETHOD) )
	    {
		    result.setMsg ("ARM_ERR: calib method is not of a good type");
		    return ARM_KO;
	    }

		calibmethod->GetModelFitterDes()->SetGetDetails( detailFlag );

		string txt( "Detail Mode:" );
		txt += detailFlag? "On" : "Off";
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}




extern long ARMLOCAL_GetDurationFromCalibMethod(
	long C_calibmethodId,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* calibmethod= NULL;

	try
	{   
	    if( !GetObjectFromId( &calibmethod, C_calibmethodId, ARM_CALIBMETHOD) )
	    {
		    result.setMsg ("ARM_ERR: calib method is not of a good type");
		    return ARM_KO;
	    }

		double duration = calibmethod->GetDuration();
		result.setDouble(duration);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

/// General types to Get from a Calculator
#define CM_PORTFOLIO                    "PORTFOLIO"
#define CM_LINKEDCALIBMETHOD            "LINKEDCALIBMETHOD"
#define CM_PREVIOUSCALIBMATHOD          "PREVIOUSCALIBMATHOD"
#define CM_CALIBPARAMS                  "CALIBPARAMS"
#define CM_CALIBMETHODTYPE              "CALIBMETHODTYPE"


extern string DataCalibMethodGetClass(const CCString& dataType)
{
    
        /// Find the model name for interface
        string typestr = CCSTringToSTLString(dataType);
        string typeToGet(stringGetUpper(typestr));

        if(typeToGet == CM_PORTFOLIO)
        {
            return LOCAL_PF_CLASS;
        }
        else if (typeToGet == CM_LINKEDCALIBMETHOD  || typeToGet == CM_PREVIOUSCALIBMATHOD)
        {
            return LOCAL_CALIBMETHOD_CLASS;
        }
        else if (typeToGet == CM_CALIBPARAMS)
        {
            return LOCAL_MODELPARAM_CLASS;
        }
        else if (typeToGet == CM_CALIBMETHODTYPE)
        {
            return LOCAL_NOID_CLASS;
        }
        return LOCAL_ANY_CLASS;
 }

ARM_CLASS_NAME  GetClassNameForCalimMethodData(const string& typeToGet)
{
    if(typeToGet == CM_PORTFOLIO)
        {
            return ARM_PORTFOLIO;
        }
        else if (typeToGet == CM_LINKEDCALIBMETHOD  || typeToGet == CM_PREVIOUSCALIBMATHOD)
        {
            return ARM_CALIBMETHOD;
        }
        else if (typeToGet == CM_CALIBPARAMS)
        {
            return ARM_MODELPARAM;
        }
        else if (typeToGet == CM_CALIBMETHODTYPE)
        {
            return ARM_OBJECT;
        }
        return ARM_OBJECT;
}


extern long ARMLOCAL_DataFromCalibMethod(long calibMethodId, 
                                    const CCString& dataType,
                                    ARM_result& result,
                                    long        objId )
{
    /// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;
    /// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_CalibMethod* calibMethod = NULL;
    ARM_Object* object=NULL;
	ARM_CLASS_NAME className;

	try
	{
		if( !GetObjectFromId( &calibMethod, calibMethodId, ARM_CALIBMETHOD ) )
		{
			result.setMsg ("ARM_ERR: CalibMethod is not of a good type");
			return ARM_KO;
		};

        string typestr = CCSTringToSTLString(dataType);
        
        /// Extract the internal data

        string typeToGet(stringGetUpper(typestr));
        if(typeToGet == CM_PORTFOLIO)
        {
            ARM_StdPortfolio* pf = static_cast< ARM_StdPortfolio* >(calibMethod->GetPortfolio()->Clone());
            className = ARM_PORTFOLIO;
            object = (ARM_Object*)pf;
        }
        else if(typeToGet == CM_LINKEDCALIBMETHOD)
        {
            ARM_CalibMethod* linkedCalibMethod = static_cast< ARM_CalibMethod* >(calibMethod->GetlinkedMethod()->Clone());
            className = ARM_CALIBMETHOD;
            object = (ARM_Object*)linkedCalibMethod;
        }
        else if(typeToGet == CM_PREVIOUSCALIBMATHOD)
        {
            ARM_CalibMethod* previousCalibMethod = static_cast< ARM_CalibMethod* >(calibMethod->GetPreviousMethod()->Clone());
            className = ARM_CALIBMETHOD;
            object = (ARM_Object*)previousCalibMethod;
        }
        else if(typeToGet == CM_CALIBPARAMS)
        {
            ARM_ModelParam* calibParam = static_cast< ARM_ModelParam* >(calibMethod->GetCalibParam()->Clone());
            object = (ARM_Object*)calibParam;
        }
        else if(typeToGet == CM_CALIBMETHODTYPE)
        {
            CCString methodtype(ARM_ArgConvReverse_CalibMethod.GetString(calibMethod->GetMethodType()).c_str());

        }
        else
        {
            result.setMsg ("ARM_ERR: DataTy is not of a good type");
			return ARM_KO;
        }

        className = GetClassNameForCalimMethodData(typeToGet);

        if(object)
        {
		    /// Assign the object in ARM cache
		    if( !assignObject( object, result, objId ) ){
			    return ARM_KO; }
		    else{
			    return ARM_OK; }
        }
        else
        {
            return ARM_KO;
        }
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a Nag Optimizer
////////////////////////////////////////////
extern long ARMLOCAL_Optimizer_Create( 
	const CCString& algoType,
	const double& Max_iter, 
	const double& tol, 
	const double& stepMax, 
	const bool& localSearch,
	const bool& printLevel,
    ARM_result&	 result, 
    long         objId  )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ModelFitterDes* ModelFitterDes = NULL;
	try
	{
		string algoTypeStr = CCSTringToSTLString(algoType);
		ModelFitterDes = new ARM_ModelFitterDes((ARM_OptimizerType) ARM_ArgConv_OptimizerTypeMethod.GetNumber(algoTypeStr),
			(size_t)Max_iter,
			tol,
			stepMax,
			localSearch,
			printLevel);

		/// assign object
		if( !assignObject( ModelFitterDes, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete ModelFitterDes;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a Nag Optimizer
////////////////////////////////////////////
extern long ARMLOCAL_Solver_Create( 
	const CCString& algoType,
	const double& Max_iter, 
	const double& xTol, 
	const double& fxTol, 
	const double& gradTol,
	const bool& printLevel,
    ARM_result&	 result, 
    long         objId  )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ModelFitterDes* ModelFitterDes = NULL;
	try
	{
		string algoTypeStr = CCSTringToSTLString(algoType);
		ModelFitterDes = new ARM_ModelFitterDes((ARM_SolverType) ARM_ArgConv_SolverTypeMethod.GetNumber(algoTypeStr),
			(size_t)Max_iter,
			xTol,
			fxTol,
			gradTol,
			printLevel);

		/// assign object
		if( !assignObject( ModelFitterDes, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete ModelFitterDes;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// For Markov Functional
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////
//// Function to create a numerical calibration method
////////////////////////////////////////////
extern long ARMLOCAL_NumericalCalibMethod_Create(
    const long&				C_CalibDateStripId,
    const VECTOR<long>&		C_VanillaSecDensitiesId,
	const long&				C_PortfolioId,
    ARM_result&				result, 
    long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* method = NULL;
	try
	{
		size_t i=0;
	
		ARM_VanillaSecDensityPtrVector vanillaSecDensities( C_VanillaSecDensitiesId.size() );

        ARM_VanillaSecurityDensity* vanillaSecDensity	= NULL;
		ARM_DateStrip* calibDateStrip					= NULL;

		calibDateStrip = dynamic_cast<ARM_DateStrip*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_CalibDateStripId));

		if (!calibDateStrip)
		{
			result.setMsg ("ARM_ERR: DateStrip is not of good type!");
			return ARM_KO;
		}


		for( i=0; i<C_VanillaSecDensitiesId.size(); ++i )
		{
			if (C_VanillaSecDensitiesId[i] != ARM_NULL_OBJECT)
			{
				vanillaSecDensity = dynamic_cast<ARM_VanillaSecurityDensity*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_VanillaSecDensitiesId[i]));

			    if (!vanillaSecDensity)
				{
					result.setMsg ("ARM_ERR: vanilla security is not of good type!");
					return ARM_KO;
				}
				
				vanillaSecDensities[i] = ARM_VanillaSecDensityPtr( static_cast<ARM_VanillaSecurityDensity*> (vanillaSecDensity->Clone()) );
			}

		}

		/// Get portfolio if provided
		ARM_StdPortfolioPtr portfolioPtr(NULL);   

		if( C_PortfolioId != ARM_NULL_OBJECT )
		{
			ARM_StdPortfolio* portfolio = NULL;   
			if( !GetObjectFromId( &portfolio, C_PortfolioId, ARM_PORTFOLIO) )
			{
				result.setMsg ("ARM_ERR: portefolio is not of a good type");
				return ARM_KO;
			}
			portfolioPtr = ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio->Clone());
		}
		


		string MethodTypeStr = "Numerical";
		ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);

		string MktTargetStr = "UNKNOWN_TAR";
		ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);

		method = new ARM_CalibMethod (	portfolioPtr,
										ARM_ModelParamVector(0),
										methodType,
										100, /// Max Iter (Newton Raphson)
										mktTargetType,
										NULL, 
										NULL, 
										false, 
										0, /// C_factorNb
										1, /// NbIter
										true, 
										ARM_DateStripPtr((ARM_DateStrip*)calibDateStrip->Clone()), 
										vanillaSecDensities );

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Shifted Lognormal density fctor
////////////////////////////////////////////
long ARMLOCAL_SLNDensityFunctor_Create(
	const double&			C_Volatility, 
	const double&			C_Shift, 
	ARM_result&				result, 
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ShiftedLNDensityFunctor* method = NULL;
	try
	{

		method = new ARM_ShiftedLNDensityFunctor(C_Volatility, C_Shift);

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Mixture density fctor
////////////////////////////////////////////
long ARMLOCAL_MixtureDensityFunctor_Create(
	const double&			C_Volatility1, 
	const double&			C_Volatility2, 
	const double&			C_Alpha, 
	const double&			C_Lambda, 
	ARM_result&				result, 
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_MixtureDensityFunctor* densFunc = NULL;
	try
	{

		densFunc = new ARM_MixtureDensityFunctor(C_Volatility1, C_Volatility2, C_Alpha, C_Lambda);

		/// assign object
		if( !assignObject( densFunc, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete densFunc;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Mixture density fctor
////////////////////////////////////////////
long ARMLOCAL_MixtureDensityFunctor_CreateWithATMVol(
	const double&			C_Fwd,
	const double&			C_Maturity,
	const double&			C_ATMVol, 
	const double&			C_DecVol, 
	const double&			C_Alpha, 
	const double&			C_Lambda, 
	ARM_result&				result, 
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_MixtureDensityFunctor* densFunc = NULL;
	try
	{

		densFunc = new ARM_MixtureDensityFunctor(C_Fwd, C_Maturity, C_ATMVol, C_DecVol, C_Alpha, C_Lambda);

		/// assign object
		if( !assignObject( densFunc, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete densFunc;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// SABR density fctor
////////////////////////////////////////////
long ARMLOCAL_SABRDensityFunctor_Create(
	const double&			C_Alpha, 
	const double&			C_Beta, 
	const double&			C_Rho, 
	const double&			C_Nu, 
	const long&				C_SabrType,
	const double&			C_GridSize,
	ARM_result&				result, 
    long					objId )
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_SABRDensityFunctor* method = NULL;
	try
	{
		method = new ARM_SABRDensityFunctor(C_Alpha, C_Beta, C_Rho, C_Nu, C_SabrType, (size_t)C_GridSize);

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_BiSABRDensityFunctor_Create(
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
    long					objId )
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_BiSABRDensityFunctor* method = NULL;
	try
	{
		method = new ARM_BiSABRDensityFunctor(C_Alpha1, C_Beta1, C_Rho1, C_Nu1, C_Alpha2, C_Beta2, C_Rho2, C_Nu2, C_RhoS1S2, C_RhoS1V2, C_RhoS2V1, C_RhoV1V2, C_SabrType, (size_t)C_GridSize);

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Mixture density fctor
////////////////////////////////////////////
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
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_HestonDensityFunctor* densFunc = NULL;
	try
	{

		densFunc = new ARM_HestonDensityFunctor(C_V0, C_Kappa, C_Theta, C_Rho, C_VVol, C_Shift, C_Level, C_Sigma);

		/// assign object
		if( !assignObject( densFunc, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete densFunc;
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_NormalHestonDensityFunctor_Create(
	const double&			C_Fwd,
	const double&			C_V0,
	const double&			C_Kappa,
	const double&			C_Theta,
	const double&			C_VVol,
	const double&			C_Rho,
	const double&			C_Level,
	ARM_result&				result,
	long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_NormalHestonDensityFunctor * method = NULL;
	try
	{
		method = new ARM_NormalHestonDensityFunctor(C_Fwd, C_V0, C_Kappa, C_Theta, C_VVol, C_Rho, C_Level);

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// SPLINE density fctor
////////////////////////////////////////////
long ARMLOCAL_SplineDensityFunctor_Create(
	const VECTOR<double>&	C_Money, 
	const VECTOR<double>&	C_Vol, 
	const string&			C_VolType, 
	const long&				C_SmileId,
	ARM_result&				result, 
    long					objId )
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_SplineDensityFunctor* method = NULL;
	try
	{
		if (C_SmileId!=ARM_NULL_OBJECT)
		{
			ARM_SmileViewer* smile = NULL;
			smile = dynamic_cast<ARM_SmileViewer*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_SmileId));
			if (!smile)
			{
				result.setMsg ("ARM_ERR: Smile is not of good type!");
				return ARM_KO;
			}
			method = new ARM_SplineDensityFunctor(smile->GetMoney(), smile->GetVol(), smile->GetVolType());
		}
		else
		{
			std::vector<double> money(C_Money);
			std::vector<double> vol(C_Vol);
			method = new ARM_SplineDensityFunctor(money, vol, (ARM_SmileViewer::MoneyType) ARM_ArgConv_MoneyType.GetNumber(C_VolType));
		}

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_DensityFunctor_CallOption(
	const long&	C_DensityFunctorId, 
	const double& C_Forward, 
	const double& C_Strike, 
	const double& C_Maturity, 
	ARM_result&	result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId);
		ARM_DensityFunctor* fctor = NULL;
		if( !(fctor = dynamic_cast< ARM_DensityFunctor* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: fctor is not of a good type");
			return ARM_KO;
		};
		
		double price = fctor->Call_Option(C_Strike,C_Forward,C_Maturity);
		result.setDouble(price);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_DensityFunctor_Quantile(
	const long&	C_DensityFunctorId, 
	const double& C_Forward, 
	const double& C_Proba, 
	const double& C_Maturity, 
	ARM_result&	result)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId);
		ARM_DensityFunctor* fctor = NULL;
		if( !(fctor = dynamic_cast< ARM_DensityFunctor* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: fctor is not of a good type");
			return ARM_KO;
		};
		
		ARM_GP_VectorPtr pricePtr = fctor->Quantile(ARM_GP_VectorPtr(new std::vector<double>(1,C_Proba)),C_Forward,C_Maturity);
		result.setDouble((*pricePtr)[0]);
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}
////////////////////////////////////////////
//// No Vol density Fctor
////////////////////////////////////////////
long ARMLOCAL_IrFwdDensityFunctor_Create(
	ARM_result&				result, 
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_NoVolDensityFunctor* method = NULL;
	try
	{

		method = new ARM_NoVolDensityFunctor();

		/// assign object
		if( !assignObject( method, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Create a Vanilla Security with density
////////////////////////////////////////////
long ARMLOCAL_VanillaSecurityDensity_Create(
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
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_VanillaSecurityDensity* obj = NULL;
	ARM_DensityFunctor*	densityFunctor = NULL;
	
	try
	{	/// convert dates
		char _resetDate[20];
		char _startDate[20];
		char _endDate[20];


		Local_XLDATE2ARMDATE(C_ResetDate, _resetDate);
		Local_XLDATE2ARMDATE(C_StartDate, _startDate);
		Local_XLDATE2ARMDATE(C_EndDate,   _endDate);
		ARM_Date resetDate = (ARM_Date)_resetDate;
		ARM_Date startDate = (ARM_Date)_startDate;
		ARM_Date endDate   = (ARM_Date)_endDate;
		

		/// get density fctor
		densityFunctor = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId));

		if (!densityFunctor)
		{
			result.setMsg ("ARM_ERR: density functor is not of good type!");
			return ARM_KO;
		}

		/// create vanilla security with density		
		obj = new ARM_VanillaSecurityDensity( resetDate.GetJulian(), 
											  startDate.GetJulian(), 
											  endDate.GetJulian(), 
											  ARM_DensityFunctorPtr(static_cast<ARM_DensityFunctor*>(densityFunctor->Clone())),
											  C_Frequency,
											  C_DayCount,
											  C_StubRule,
											  C_Weight,
											  C_AdjFwdAdd,
											  C_AdjFwdMult);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete obj;
		delete densityFunctor;
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_VanillaSecurityDensitySpread_Create(
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
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_VanillaSecurityDensity* obj = NULL;
	ARM_DensityFunctor*	densityFunctor = NULL;
	
	try
	{	/// convert dates
		char _resetDate[20];
		char _startDate1[20];
		char _endDate1[20];
		char _startDate2[20];
		char _endDate2[20];


		Local_XLDATE2ARMDATE(C_ResetDate, _resetDate);
		Local_XLDATE2ARMDATE(C_StartDate1, _startDate1);
		Local_XLDATE2ARMDATE(C_EndDate1,   _endDate1);
		Local_XLDATE2ARMDATE(C_StartDate2, _startDate2);
		Local_XLDATE2ARMDATE(C_EndDate2,   _endDate2);
		ARM_Date resetDate = (ARM_Date)_resetDate;
		ARM_Date startDate1 = (ARM_Date)_startDate1;
		ARM_Date endDate1   = (ARM_Date)_endDate1;
		ARM_Date startDate2 = (ARM_Date)_startDate2;
		ARM_Date endDate2   = (ARM_Date)_endDate2;
		

		/// get density fctor
		densityFunctor = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId));

		if (!densityFunctor)
		{
			result.setMsg ("ARM_ERR: density functor is not of good type!");
			return ARM_KO;
		}

		/// create vanilla security with density		
		obj = new ARM_VanillaSecurityDensitySpread( resetDate.GetJulian(), 
											  startDate1.GetJulian(), 
											  endDate1.GetJulian(), 
											  startDate2.GetJulian(),
											  endDate2.GetJulian(),
											  ARM_DensityFunctorPtr(static_cast<ARM_DensityFunctor*>(densityFunctor->Clone())),
											  C_Frequency1,
											  C_DayCount1,
											  C_Frequency2,
											  C_DayCount2,
											  C_StubRule  );
		obj->setWeight(C_Weight);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete obj;
		delete densityFunctor;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Create a Vanilla Security FX with density
////////////////////////////////////////////
long ARMLOCAL_VanillaSecurityDensityFX_Create(
	const double&			C_ResetDate,
	const long&				C_DensityFunctorId,
	const long&				C_DomCurveId,
	const long&				C_ForCurveId,
	const double&			C_FXSpot,
	ARM_result&				result, 
    long					objId)
{
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_VanillaSecurityDensity* obj = NULL;
	ARM_DensityFunctor*	densityFunctor = NULL;
	ARM_ZeroCurve* domCurve = NULL;
	ARM_ZeroCurve* forCurve = NULL;
	
	try
	{	/// convert dates
		char _resetDate[20];


		Local_XLDATE2ARMDATE(C_ResetDate, _resetDate);
		ARM_Date resetDate = (ARM_Date)_resetDate;
		

		/// get density fctor
		densityFunctor = dynamic_cast<ARM_DensityFunctor*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_DensityFunctorId));

		if (!densityFunctor)
		{
			result.setMsg ("ARM_ERR: density functor is not of good type!");
			return ARM_KO;
		}

		domCurve = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_DomCurveId));

		if (!domCurve)
		{
			result.setMsg ("ARM_ERR: Domestic Curve is not of good type!");
			return ARM_KO;
		}

		forCurve = dynamic_cast<ARM_ZeroCurve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_ForCurveId));

		if (!forCurve)
		{
			result.setMsg ("ARM_ERR: Foreign Curve is not of good type!");
			return ARM_KO;
		}

		/// create vanilla security with density		
		obj = new ARM_VanillaSecurityDensityFX( resetDate.GetJulian(), 
											  ARM_DensityFunctorPtr(static_cast<ARM_DensityFunctor*>(densityFunctor->Clone())),
											  ARM_ZeroCurvePtr(static_cast<ARM_ZeroCurve*>(domCurve->Clone())),
											  ARM_ZeroCurvePtr(static_cast<ARM_ZeroCurve*>(forCurve->Clone())),
											  C_FXSpot);

		/// assign object
		if( !assignObject( obj, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete obj;
		delete densityFunctor;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
//// For HW2F
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

////////////////////////////////////////////
//// Function to create a numerical calibration method
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
    long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_CalibMethod* method1 = NULL;
	ARM_CalibMethod* method2 = NULL;
	ARM_CalibMethod* method3 = NULL;
	try
	{
		int N=0;
		ARM_StdPortfolioPtr portfolioPtr1(NULL);   
		ARM_StdPortfolio* portfolio1 = NULL;   
		if(!GetObjectFromId( &portfolio1, C_PortfolioId1, ARM_PORTFOLIO))
        {
		   result.setMsg ("ARM_ERR: portfolio is not of a good typ");
		   return ARM_KO;
        }
        portfolioPtr1 = ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio1->Clone());
		N=1;

		ARM_StdPortfolioPtr portfolioPtr2(NULL);   
		ARM_StdPortfolio* portfolio2 = NULL;   
		if (C_PortfolioId2 != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &portfolio2, C_PortfolioId2, ARM_PORTFOLIO))
            {
			    result.setMsg ("ARM_ERR: portfolio is not of a good typ");
			    return ARM_KO;
            }
			portfolioPtr2 = ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio2->Clone());
			N=2;
        }
		
		ARM_StdPortfolioPtr portfolioPtr3(NULL);   
		ARM_StdPortfolio* portfolio3 = NULL;   
		if (C_PortfolioId3 != ARM_NULL_OBJECT)
        {
            if(!GetObjectFromId( &portfolio3, C_PortfolioId3, ARM_PORTFOLIO))
            {
			    result.setMsg ("ARM_ERR: portfolio is not of a good typ");
			    return ARM_KO;
            }
			portfolioPtr3 = ARM_StdPortfolioPtr((ARM_StdPortfolio*)portfolio3->Clone());
			N=3;
        }
		
		ARM_ModelParamVector modelParams1;
       	ARM_ModelParam* modelParam1 = NULL;
	    if( !GetObjectFromId( &modelParam1, param1Id, ARM_MODELPARAM ) )
	    {
		    result.setMsg ("ARM_ERR: model parameter is not of a good type");
		    return ARM_KO;
	    }
		modelParams1.push_back(modelParam1);

		ARM_ModelParamVector modelParams2;
       	ARM_ModelParam* modelParam2 = NULL;
	    if (param2Id != ARM_NULL_OBJECT)
        {
			if( !GetObjectFromId( &modelParam2, param2Id, ARM_MODELPARAM ) )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams2.push_back(modelParam2);
		}
		else
		{
			if (N>1)
			{
				result.setMsg ("ARM_ERR: empty modelparam2!");
				return ARM_KO;
			}
		}

		ARM_ModelParamVector modelParams3;
       	ARM_ModelParam* modelParam3 = NULL;
	    if (param3Id != ARM_NULL_OBJECT)
        {
			if( !GetObjectFromId( &modelParam3, param3Id, ARM_MODELPARAM ) )
			{
				result.setMsg ("ARM_ERR: model parameter is not of a good type");
				return ARM_KO;
			}
			modelParams3.push_back(modelParam3);
		}
		else
		{
			if (N>2)
			{
				result.setMsg ("ARM_ERR: empty modelparam3!");
				return ARM_KO;
			}
		}

		string MethodTypeStr = "HW2FOnly";
		ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);

		string MktTargetStr = "UNKNOWN_TAR";
		ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);


		method1 = new ARM_CalibMethod (	portfolioPtr1,
										modelParams1,
										methodType,
										withOptim*10+N, /// moche!
										mktTargetType,
										NULL, 
										NULL, 
										false, 
										0, /// C_factorNb
										1, /// NbIter
										false);

		if (N>1)
			method2 = new ARM_CalibMethod (	portfolioPtr2,
										modelParams2,
										methodType,
										withOptim*10+N, /// moche!
										mktTargetType,
										method1, 
										NULL, 
										false, 
										0, /// C_factorNb
										1, /// NbIter
										false);
		if (N>2)
			method3 = new ARM_CalibMethod (	portfolioPtr3,
										modelParams3,
										methodType,
										withOptim*10+N, /// moche!
										mktTargetType,
										method2, 
										NULL, 
										false, 
										0, /// C_factorNb
										1, /// NbIter
										false);
		
		if (N>1)
			delete method1;
		if (N>2)
			delete method2;

		/// assign object
		if( !assignObject( N==1?method1:N==2?method2:method3, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete method1;
		delete method2;
		delete method3;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Basket decomp
////////////////////////////////////////////
extern long ARMLOCAL_BasketDecomp_Create(
    const VECTOR<long>&		C_securitiesId,
	const VECTOR<long>&		C_modelsId,
	const long&				C_datestripId,
	const long&				C_mkmoId,
	const VECTOR<double>&	C_weights,
	const double&			C_side,
	const string&			C_method,
	const string&			C_strike,
	const double&			C_notional,
	const long&				C_notionalId,
	const double&			C_fees,
	const long&				C_feesId,
	ARM_result&				result, 
    long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* feesProfile = NULL;
    bool isFees=false;

	try
	{

		int size = C_securitiesId.size();

		vector<ARM_Security*> securities(size, NULL);
        int i;
        for (i = 0; i < size; i++)
		{
            if( !GetObjectFromId( &securities[i], C_securitiesId[i], ARM_SECURITY) )
	        {
		        result.setMsg ("ARM_ERR: Products Vector is not of a good type, please check all products");
		        return ARM_KO;
	        }
        }

		vector<ARM_Model*> models(size, NULL);
        for (i = 0; i < size; i++)
		{
            if( !GetObjectFromId( &models[i], C_modelsId[i], ARM_MODEL) )
	        {
		        result.setMsg ("ARM_ERR: Models Vector is not of a good type, please check all models");
		        return ARM_KO;
	        }
        }

		ARM_DateStrip* ds = NULL;
		ds = dynamic_cast<ARM_DateStrip*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_datestripId));
		if (!ds)
		{
			result.setMsg ("ARM_ERR: DateStrip is not of good type!");
			return ARM_KO;
		}

		//Notional Curve
        if (C_notionalId != ARM_NULL_OBJECT)
        {
		    notionalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_notionalId));

		    if (!notionalProfile)
		    {
			    result.setMsg ("ARM_ERR: notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            notionalProfile = new ARM_ReferenceValue(C_notional);
            isNotional=true;
        }

		//Fees Curve
        if (C_feesId != ARM_NULL_OBJECT)
        {
		    feesProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_feesId));

		    if (!feesProfile)
		    {
			    result.setMsg ("ARM_ERR: Fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            feesProfile = new ARM_ReferenceValue(C_fees);
            isFees=true;
        }

		ARM_PricingModel* mkmo = NULL;
		mkmo = dynamic_cast<ARM_PricingModel*>(LOCAL_PERSISTENT_OBJECTS->GetObject(C_mkmoId));
		if (!ds)
		{
			result.setMsg ("ARM_ERR: Market IR model is not of good type!");
			return ARM_KO;
		}

		std::vector<double> weights(C_weights);

		ARM_BasketCalib* bskt = new ARM_BasketCalib(
			ds,
			*notionalProfile,
			*feesProfile,
			C_side,
			(ARM_BasketCalib::BasketType) ARM_ArgConv_BasketCalibrationType.GetNumber(C_method),
			(ARM_BasketCalib::BasketStrike) ARM_ArgConv_BasketCalibrationStrike.GetNumber(C_strike));

		bskt->Compute(securities,models,weights);
		bskt->Price(mkmo);

		 // Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

       	if (isFees)
			delete feesProfile;
		feesProfile = NULL;

		/// assign object
		if( !assignObject( bskt, result, objId ) ){
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

#define BSKT_ACCESS_VNS_PORT				"PORTFOLIO"
#define BSKT_ACCESS_BSKT_WEIGHTS			"BASKET"

extern string BsktGetTypeToClass(const string& typeToGet, long gcId)
{
	string typeToGetUpper = ARM::stringGetUpper(typeToGet);
	if (typeToGetUpper == BSKT_ACCESS_VNS_PORT)
    {
		return LOCAL_PF_CLASS;
	}
	if (typeToGetUpper == BSKT_ACCESS_BSKT_WEIGHTS)
    {
		return LOCAL_GP_MATRIX_CLASS;
	}
	else
	{
		return LOCAL_ANY_CLASS;
	}
}

///////////////////////////////////////////////
//// Function to get data from a Basket Decomp
///////////////////////////////////////////////
extern long ARMLOCAL_Basket_Get(
        const long& basketId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_BasketCalib* basket;
    ARM_Object* object=NULL;

	try
	{
		basket = dynamic_cast<ARM_BasketCalib*>(LOCAL_PERSISTENT_OBJECTS->GetObject(basketId));

		if (!basket)
		{
			result.setMsg ("ARM_ERR: Basket is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));

        if( typeToGet == BSKT_ACCESS_VNS_PORT)
        {
            object = basket->GetPortfolio()->Clone();

        }
		else if( typeToGet == BSKT_ACCESS_BSKT_WEIGHTS)
		{
			object = basket->GetBasketWeights().Clone();
		}
		else
		{
			return ARM_KO;
		}

		/// Assign the object in ARM cache
		if( !assignObject( object, result, objId ) )
        {
			delete object;
			return ARM_KO;
		}
		else
			return ARM_OK;

	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Smile Viewer
////////////////////////////////////////////
extern long ARMLOCAL_SmileViewer_Create(
    const long&		C_securityId,
	const long&		C_modelId,
	const VECTOR<double>&	C_moneyness,
	const string&	C_moneyType,
	const VECTOR<double>&	C_strikes,
	ARM_result&				result, 
    long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{

		ARM_Security* sec = NULL;
	    if( !GetObjectFromId( &sec, C_securityId, ARM_SECURITY ) )
	    {
		    result.setMsg ("ARM_ERR: security is not of a good type");
		    return ARM_KO;
	    }

		ARM_Model* mod = NULL;
	    if( !GetObjectFromId( &mod, C_modelId, ARM_MODEL ) )
	    {
		    result.setMsg ("ARM_ERR: model is not of a good type");
		    return ARM_KO;
	    }

		std::vector<double> moneyness(C_moneyness);
		std::vector<double> strikes(C_strikes);

		ARM_SmileViewer* smvw;

		if (moneyness.size()>0 || strikes.size()>0)
		{
			smvw = new ARM_SmileViewer(
				moneyness,
				(ARM_SmileViewer::MoneyType) ARM_ArgConv_MoneyType.GetNumber(C_moneyType),
				strikes);

			smvw->Compute(sec,mod);
		}
		else
		{
			result.setMsg ("ARM_ERR: moneyness OR strike should have size>0");
		    return ARM_KO;
		}
			
		
		/// assign object
		if( !assignObject( smvw, result, objId ) ){
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

extern long ARMLOCAL_DensityFunctorGen_Create(
    const long&		C_securityId,
	const long&		C_modelId,
	const double&	C_decStrike,
	const bool&		C_isDirect,
	const double&	C_minProba,
	const double&	C_maxProba,
	ARM_result&				result, 
    long					objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{

		ARM_Security* sec = NULL;
	    if( !GetObjectFromId( &sec, C_securityId, ARM_SECURITY ) )
	    {
		    result.setMsg ("ARM_ERR: security is not of a good type");
		    return ARM_KO;
	    }

		ARM_Model* mod = NULL;
	    if( !GetObjectFromId( &mod, C_modelId, ARM_MODEL ) )
	    {
		    result.setMsg ("ARM_ERR: model is not of a good type");
		    return ARM_KO;
	    }

		ARM_GenDensityFunctor* functor = new ARM_GenDensityFunctor((ARM_Security*)sec->Clone(),(ARM_Model*)mod->Clone(),C_decStrike, C_minProba,C_maxProba,C_isDirect);
		
		/// assign object
		if( !assignObject( functor, result, objId ) ){
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

//***********************************************************************
//**********  ARM_QDensityFunctorFunctor_Create  ************************
//***********************************************************************

long ARM_QDensityFunctor_Create::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_QDensityFunctor* densityFunctor = NULL;

	try
	{
		/// crm tracing
		densityFunctor = new ARM_QDensityFunctor(
			genericParams->GetParamValue("vol").GetDouble(), 
			genericParams->GetParamValue("q").GetDouble());

		// assign object
		return  ( !assignObject( densityFunctor, result, objId ) )  ?  ARM_KO : ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	densityFunctor;

		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}