/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_calculators_local.cpp,v $
 * Revision 1.1  2004/02/07 15:08:43  ebenhamou
 * Initial version
 *
 */


////////////////////////////////////////////////
/// This file is for global interpolation and
/// various conversion routine for the generic
/// pricer... Everything specific to the generic
/// pricer should be included here and not in ARM_local_interglob.cpp
/// to avoid exporting the generic pricer in the 
/// local dllarm project (activex project)
/// since the interglob.h file is included in the 
/// local dllarm project (activex project)
////////////////////////////////////////////////

#ifndef VBARM_LOCAL_EXPORTS

#include <ARM\libarm_local\firstToBeIncluded.h>
#include "ARM_gp_local_interglob.h"
#include <GP_Infra\gpinfra\argconvdefault.h>
#include <GP_ClosedForms\gpclosedforms\argconvdefault_cf.h>
#include <GP_Inflation\gpinflation\argconvdefault.h>
#include <GP_NumMethods\gpnummethods\argconvdefault.h>
#include <GP_NumLib\gpnumlib\argconvdefault.h>
#include <GP_Models\gpmodels\argconvdefault.h>

#include "expt.h"
#include "ARM_local_interglob.h"


long ARM_ConvGPNumeraire( const CCString& numeraireType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_Numeraire.GetNumber(CCSTringToSTLString(numeraireType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


long ARM_ConvGP_CFSpreadDigitalOption_Formula_DerivedStruct( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_SpreadDigitalOption_Formula_DerivedStruct.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGP_CFSABR_ImplicitVol_Formula_Extended_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_SABR_ImplicitVol_Formula_Extended_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFOptimization_ObjectiveFuntion_Algorithm_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_InOut_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_UpDown_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFBS_EuropeanBarriere_Formula_OptionType_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

     
long ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFBS_PartialTime_Barriere_End_Formula_OptionType_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag( const CCString& nType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConvGP_CFHeston_Vector_InterpolationMethod_Flag.GetNumber(CCSTringToSTLString(nType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPModelParam( const CCString& modelParamType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_ModelParam.GetNumber(CCSTringToSTLString(modelParamType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPModelParamDataType( const CCString& modelParamDataType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_ModelParamDataType.GetNumber(CCSTringToSTLString(modelParamDataType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


long ARM_ConvGPCalib2DDirection( const CCString& calibDirection, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_CalibDirection2DMethod.GetNumber(CCSTringToSTLString(calibDirection));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPCalibType( const CCString& calibType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_CalibMethod.GetNumber(CCSTringToSTLString(calibType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPTargetFuncType( const CCString& calibType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_TargetFuncMethod.GetNumber(CCSTringToSTLString(calibType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


long ARM_ConvGPInfSwoptComputationMethod( const CCString& type, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_InfSwoptComputationMethod.GetNumber(CCSTringToSTLString(type));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPSamplerType( const CCString& samplerType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_SamplerType.GetNumber(CCSTringToSTLString(samplerType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPTruncatorType( const CCString& truncatorType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_TruncatorType.GetNumber(CCSTringToSTLString(truncatorType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPSchedulerType( const CCString& schedulerType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_SchedulerType.GetNumber(CCSTringToSTLString(schedulerType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPReconnectorType( const CCString& reconnectorType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_ReconnectorType.GetNumber(CCSTringToSTLString(reconnectorType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPSmootherType( const CCString& smootherType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_SmootherType.GetNumber(CCSTringToSTLString(smootherType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


long ARM_ConvGPImpSamplerType( const CCString& impSamplerType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_ImpSamplerType.GetNumber(CCSTringToSTLString(impSamplerType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


long ARM_ConvGPPathSchemeType( const CCString& pathSchemeType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_PathSchemeType.GetNumber(CCSTringToSTLString(pathSchemeType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvGPODESolverType( const CCString& solverType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_ODESolverType.GetNumber(CCSTringToSTLString(solverType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}

long ARM_ConvHWSVFormulaType( const CCString& formulaType, ARM_result& result)
{
    try
    {
        return ARM::ARM_ArgConv_HWSVFormula.GetNumber(CCSTringToSTLString(formulaType));
    }

	catch(Exception& x)
	{
	    result.setRetCode (ARM_DEFAULT_ERR);
	    result.setMsgString(x.GetErrorString());
	    return ARM_DEFAULT_ERR;
	}
}


#endif